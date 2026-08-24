#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"

namespace {

constexpr int kMcpPdg = 1000222;
constexpr int kDyEmitterType = 2;
constexpr int kDyProductionMode = 2;

enum class GeometryId : int { ArgoNeuT = 0, TwoByTwo = 1 };
enum class SpectraMode : int { All = 0, Accepted = 1, None = 2 };

struct DyParameters {
    double mass = 0.1;
    double epsilon = 0.01;
};

class Sigma2qqbar2GammaMcpPair final : public Pythia8::Sigma2Process {
  public:
    explicit Sigma2qqbar2GammaMcpPair(const DyParameters& parameters)
        : parameters_(parameters) {}

    void sigmaKin() override {
        const double mass2 = parameters_.mass * parameters_.mass;
        physical_ = sH > 4.0 * mass2;
        if (!physical_) {
            sigma0_ = 0.0;
            return;
        }

        const double massRatio = mass2 / sH;
        beta_ = std::sqrt(1.0 - 4.0 * massRatio);
        cosTheta_ = std::clamp((tH - uH) / sH, -1.0, 1.0);
        longitudinal_ = 4.0 * massRatio;
        sigma0_ = M_PI * alpEM * alpEM / sH2
                * parameters_.epsilon * parameters_.epsilon * beta_;
    }

    double sigmaHat() override {
        if (!physical_ || id1 * id2 >= 0 || std::abs(id1) != std::abs(id2)) {
            return 0.0;
        }

        const int incomingId = std::abs(id1);
        if (incomingId < 1 || incomingId > 6) return 0.0;

        const double charge = coupSMPtr->ef(incomingId);
        const double cosTheta2 = cosTheta_ * cosTheta_;
        const double angular = 1.0 + cosTheta2
                             + longitudinal_ * (1.0 - cosTheta2);
        return charge * charge * sigma0_ * angular / 3.0;
    }

    void setIdColAcol() override {
        setId(id1, id2, kMcpPdg, -kMcpPdg);
        setColAcol(1, 0, 0, 1, 0, 0, 0, 0);
        if (id1 < 0) swapColAcol();
    }

    std::string name() const override {
        return "q qbar -> gamma* -> chi chibar";
    }

    int code() const override { return 9901; }
    std::string inFlux() const override { return "qqbarSame"; }
    bool isSChannel() const override { return true; }
    int id3Mass() const override { return kMcpPdg; }
    int id4Mass() const override { return kMcpPdg; }

  private:
    DyParameters parameters_;
    bool physical_ = false;
    double beta_ = 0.0;
    double cosTheta_ = 0.0;
    double sigma0_ = 0.0;
    double longitudinal_ = 0.0;
};

std::shared_ptr<Pythia8::SigmaProcess> makeDyProcess(
    const DyParameters& parameters) {
    return std::make_shared<Sigma2qqbar2GammaMcpPair>(parameters);
}

struct Options {
    int seed = 1;
    int nThreads = 1;
    int jobId = 0;
    long long nEvents = 10000;
    double mass = 0.1;
    double epsilon = 0.01;
    double mHatMin = -1.0;
    GeometryId geometry = GeometryId::TwoByTwo;
    SpectraMode spectraMode = SpectraMode::All;
    int spectraPrescale = 1;
    std::string output = "dy_mcp.root";
    std::string emitterName = "dy";
    std::string productionMode = "drell_yan";
    std::string beamConfig = "beam.config";
    std::string momentumConfig = "momentum.config";
    std::string productionConfig;
    std::string pdfSet = "5";
    bool hardOnly = false;
    bool quiet = false;
};

struct Projection {
    bool passed = false;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct SpectraRow {
    int event_index = 0;
    int seed = 0;
    int emitter_pdg = 0;
    int emitter_type = kDyEmitterType;
    int production_mode = kDyProductionMode;
    int geometry_id = 0;
    int mcp_pdg = 0;
    int mother_pdg = 0;
    int mother_index = -1;
    int accepted = 0;
    int passed_geometry = 0;
    int mother_status = 0;
    int mother_is_primary_like = 0;
    int mother_mother_pdg = 0;
    int from_feeddown_flag = 0;
    int id1 = 0;
    int id2 = 0;
    double mcp_mass_GeV = 0.0;
    double event_weight = 1.0;
    double x1 = 0.0;
    double x2 = 0.0;
    double q2_fac_GeV2 = 0.0;
    double mhat_GeV = 0.0;
    double px_GeV = 0.0;
    double py_GeV = 0.0;
    double pz_GeV = 0.0;
    double E_GeV = 0.0;
    double p_GeV = 0.0;
    double theta_rad = 0.0;
    double theta_x_rad = 0.0;
    double theta_y_rad = 0.0;
    double eta = 0.0;
    double phi_rad = 0.0;
    double x_at_detector_m = 0.0;
    double y_at_detector_m = 0.0;
    double z_detector_m = 0.0;
    double mother_px_GeV = 0.0;
    double mother_py_GeV = 0.0;
    double mother_pz_GeV = 0.0;
    double mother_E_GeV = 0.0;
    double mother_p_GeV = 0.0;
    double mother_theta_rad = 0.0;
    double pair_px_GeV = 0.0;
    double pair_py_GeV = 0.0;
    double pair_pz_GeV = 0.0;
    double pair_E_GeV = 0.0;
    double pair_mass_GeV = 0.0;
    double pair_pT_GeV = 0.0;
};

struct ThreadResult {
    int threadId = 0;
    int seed = 0;
    long long nEvents = 0;
    long long nFailures = 0;
    long long nMalformed = 0;
    long long nPairs = 0;
    long long nMcp = 0;
    long long nAccepted = 0;
    double sigmaMb = 0.0;
    double sigmaErrMb = 0.0;
    std::vector<SpectraRow> spectra;
};

std::string geometryName(GeometryId geometry) {
    return geometry == GeometryId::ArgoNeuT ? "argoneut" : "2x2";
}

GeometryId parseGeometry(const std::string& value) {
    if (value == "argoneut") return GeometryId::ArgoNeuT;
    if (value == "2x2" || value == "twobytwo") return GeometryId::TwoByTwo;
    throw std::runtime_error("unknown geometry: " + value);
}

Projection project(double px, double py, double pz, GeometryId geometry) {
    Projection result;
    if (pz <= 0.0) return result;

    if (geometry == GeometryId::ArgoNeuT) {
        result.z = 1000.0;
        result.x = px / pz * result.z;
        result.y = py / pz * result.z;
        result.passed = result.x >= -0.24 && result.x <= 0.24
                     && result.y >= -0.20 && result.y <= 0.20;
        return result;
    }

    result.z = 1040.0;
    result.x = px / pz * result.z;
    result.y = py / pz * result.z;
    const bool inX = (result.x >= -0.65 && result.x <= -0.05)
                  || (result.x >= 0.05 && result.x <= 0.65);
    result.passed = inX && result.y >= -0.70 && result.y <= 0.70;
    return result;
}

void printUsage(const char* program) {
    std::cout
        << "Usage: " << program
        << " <seed> <nThreads> <nEvents> <mcpMassGeV> <geometry>"
        << " <emitterName> <productionMode> <outputFile> [options]\n"
        << "  --epsilon VALUE     Generation charge (default: 0.01)\n"
        << "  --mhat-min GEV       Comparison cut; default is 2 * MCP mass\n"
        << "  --beam-config FILE   Default: beam.config\n"
        << "  --momentum-config F  Default: momentum.config\n"
        << "  --production-config FILE\n"
        << "  --pdf-set VALUE      Optional PYTHIA PDF:pSet override\n"
        << "  --write-spectra all|accepted|none\n"
        << "  --spectra-prescale N\n"
        << "  --hard-only          Disable ISR, beam remnants and hadronization\n"
        << "  --job-id N\n"
        << "  --quiet | --batch\n";
}

Options parseOptions(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        printUsage(argv[0]);
        std::exit(0);
    }
    if (argc < 9) {
        printUsage(argv[0]);
        throw std::runtime_error("not enough arguments");
    }

    Options options;
    options.seed = std::stoi(argv[1]);
    options.nThreads = std::stoi(argv[2]);
    options.nEvents = std::stoll(argv[3]);
    options.mass = std::stod(argv[4]);
    options.geometry = parseGeometry(argv[5]);
    options.emitterName = argv[6];
    options.productionMode = argv[7];
    options.output = argv[8];

    if (options.emitterName != "dy") {
        throw std::runtime_error("DY generator requires emitterName=dy");
    }
    if (options.productionMode != "drell_yan") {
        throw std::runtime_error("DY generator requires productionMode=drell_yan");
    }

    for (int i = 9; i < argc; ++i) {
        const std::string key = argv[i];
        auto value = [&]() {
            if (++i >= argc) throw std::runtime_error("missing value for " + key);
            return std::string(argv[i]);
        };

        if (key == "--epsilon") options.epsilon = std::stod(value());
        else if (key == "--mhat-min") options.mHatMin = std::stod(value());
        else if (key == "--beam-config") options.beamConfig = value();
        else if (key == "--momentum-config") options.momentumConfig = value();
        else if (key == "--production-config") options.productionConfig = value();
        else if (key == "--pdf-set") options.pdfSet = value();
        else if (key == "--mode") {
            const std::string mode = value();
            if (mode != "fixed-events") {
                throw std::runtime_error("DY production supports fixed-events mode only");
            }
        } else if (key == "--n-events") {
            options.nEvents = std::stoll(value());
        }
        else if (key == "--write-spectra") {
            const std::string mode = value();
            if (mode == "all") options.spectraMode = SpectraMode::All;
            else if (mode == "accepted") options.spectraMode = SpectraMode::Accepted;
            else if (mode == "none") options.spectraMode = SpectraMode::None;
            else throw std::runtime_error("unknown spectra mode: " + mode);
        } else if (key == "--spectra-prescale") {
            options.spectraPrescale = std::stoi(value());
        } else if (key == "--hard-only") {
            options.hardOnly = true;
        } else if (key == "--job-id") {
            options.jobId = std::stoi(value());
        } else if (key == "--quiet" || key == "--batch") {
            options.quiet = true;
        } else {
            throw std::runtime_error("unknown option: " + key);
        }
    }

    if (options.seed < 1 || options.seed > 900000000) {
        throw std::runtime_error("seed must be in [1, 900000000]");
    }
    if (options.nThreads <= 0 || options.nEvents <= 0
        || options.mass <= 0.0 || options.epsilon <= 0.0
        || options.spectraPrescale <= 0) {
        throw std::runtime_error("threads, events, mass, epsilon and prescale must be positive");
    }

    const double threshold = 2.0 * options.mass;
    if (options.mHatMin < 0.0) options.mHatMin = threshold;
    if (options.mHatMin < threshold) {
        throw std::runtime_error("mhat-min cannot be below 2 * MCP mass");
    }
    return options;
}

std::string setting(const std::string& key, const std::string& value) {
    return key + " = " + value;
}

void configurePythia(Pythia8::Pythia& pythia, const Options& options, int seed) {
    if (!pythia.readFile(options.beamConfig)) {
        throw std::runtime_error("cannot read " + options.beamConfig);
    }
    if (!pythia.readFile(options.momentumConfig)) {
        throw std::runtime_error("cannot read " + options.momentumConfig);
    }
    if (!options.productionConfig.empty()
        && !pythia.readFile(options.productionConfig)) {
        throw std::runtime_error("cannot read " + options.productionConfig);
    }

    pythia.readString("SoftQCD:all = off");
    pythia.readString("HardQCD:all = off");
    pythia.readString(setting("PhaseSpace:mHatMin", std::to_string(options.mHatMin)));
    pythia.settings.parm("PhaseSpace:pTHatMinDiverge", 0.0, true);

    pythia.readString(setting("PDF:pSet", options.pdfSet));

    std::ostringstream particle;
    particle << kMcpPdg << ":new = chi chibar 2 0 0 "
             << std::setprecision(12) << options.mass << " 0 0 0 0";
    pythia.readString(particle.str());
    pythia.readString(std::to_string(kMcpPdg) + ":mayDecay = off");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("ParticleDecays:tau0Max = 1e12");

    pythia.readString("Random:setSeed = on");
    pythia.readString(setting("Random:seed", std::to_string(seed)));
    pythia.readString("Next:numberCount = 0");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Init:showProcesses = off");
    pythia.readString("Init:showMultipartonInteractions = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showChangedParticleData = off");

    if (options.hardOnly) {
        pythia.readString("PartonLevel:all = off");
        pythia.readString("HadronLevel:all = off");
    }

    pythia.setSigmaPtr(makeDyProcess({options.mass, options.epsilon}));
}

bool keepSpectrum(const Options& options, long long index, bool accepted) {
    if (options.spectraMode == SpectraMode::None) return false;
    if (options.spectraMode == SpectraMode::Accepted && !accepted) return false;
    return options.spectraPrescale == 1 || index % options.spectraPrescale == 0;
}

int motherIndex(const Pythia8::Event& event, int index) {
    int mother = event[index].mother1();
    if (mother > 0 && mother < event.size()) return mother;
    mother = event[index].mother2();
    if (mother > 0 && mother < event.size()) return mother;
    return -1;
}

ThreadResult runThread(
    const Options& options,
    int threadId,
    long long targetEvents
) {
    ThreadResult result;
    result.threadId = threadId;
    result.seed = options.seed + threadId;

    Pythia8::Pythia pythia("", false);
    configurePythia(pythia, options, result.seed);
    if (!pythia.init()) throw std::runtime_error("PYTHIA initialization failed");

    while (result.nEvents < targetEvents) {
        if (!pythia.next()) {
            ++result.nFailures;
            continue;
        }

        const int eventIndex = static_cast<int>(
            threadId * 1000000000LL + result.nEvents
        );
        ++result.nEvents;

        std::vector<int> mcpIndices;
        for (int i = 0; i < pythia.event.size(); ++i) {
            if (pythia.event[i].isFinal()
                && std::abs(pythia.event[i].id()) == kMcpPdg) {
                mcpIndices.push_back(i);
            }
        }

        if (mcpIndices.size() != 2) {
            ++result.nMalformed;
            continue;
        }

        ++result.nPairs;
        const Pythia8::Vec4 pair = pythia.event[mcpIndices[0]].p()
                                + pythia.event[mcpIndices[1]].p();

        for (int index : mcpIndices) {
            const Pythia8::Particle& particle = pythia.event[index];
            ++result.nMcp;
            const Projection projection = project(
                particle.px(), particle.py(), particle.pz(), options.geometry
            );
            if (projection.passed) ++result.nAccepted;
            if (!keepSpectrum(options, result.nMcp, projection.passed)) continue;

            SpectraRow row;
            row.event_index = eventIndex;
            row.seed = result.seed;
            row.geometry_id = static_cast<int>(options.geometry);
            row.mcp_pdg = particle.id();
            row.mcp_mass_GeV = options.mass;
            row.accepted = projection.passed;
            row.passed_geometry = projection.passed;
            row.event_weight = pythia.info.weight();
            row.id1 = pythia.info.id1();
            row.id2 = pythia.info.id2();
            row.x1 = pythia.info.x1();
            row.x2 = pythia.info.x2();
            row.q2_fac_GeV2 = pythia.info.Q2Fac();
            row.mhat_GeV = pythia.info.mHat();
            row.px_GeV = particle.px();
            row.py_GeV = particle.py();
            row.pz_GeV = particle.pz();
            row.E_GeV = particle.e();
            row.p_GeV = particle.pAbs();
            row.theta_rad = particle.theta();
            row.theta_x_rad = particle.pz() != 0.0
                ? particle.px() / particle.pz() : 0.0;
            row.theta_y_rad = particle.pz() != 0.0
                ? particle.py() / particle.pz() : 0.0;
            row.eta = particle.eta();
            row.phi_rad = particle.phi();
            row.x_at_detector_m = projection.x;
            row.y_at_detector_m = projection.y;
            row.z_detector_m = projection.z;
            row.pair_px_GeV = pair.px();
            row.pair_py_GeV = pair.py();
            row.pair_pz_GeV = pair.pz();
            row.pair_E_GeV = pair.e();
            row.pair_mass_GeV = pair.mCalc();
            row.pair_pT_GeV = pair.pT();

            const int mother = motherIndex(pythia.event, index);
            row.mother_index = mother;
            if (mother >= 0) {
                const auto& parent = pythia.event[mother];
                row.mother_pdg = parent.id();
                row.mother_status = parent.status();
                row.mother_is_primary_like = parent.mother1() == 1
                                           || parent.mother1() == 2;
                row.mother_px_GeV = parent.px();
                row.mother_py_GeV = parent.py();
                row.mother_pz_GeV = parent.pz();
                row.mother_E_GeV = parent.e();
                row.mother_p_GeV = parent.pAbs();
                row.mother_theta_rad = parent.theta();
            }
            result.spectra.push_back(row);
        }
    }

    result.sigmaMb = pythia.info.sigmaGen();
    result.sigmaErrMb = pythia.info.sigmaErr();
    return result;
}

void branchSpectra(TTree& tree, SpectraRow& row) {
#define BRANCH_INT(name) tree.Branch(#name, &row.name, #name "/I")
#define BRANCH_DOUBLE(name) tree.Branch(#name, &row.name, #name "/D")
    BRANCH_INT(event_index); BRANCH_INT(seed); BRANCH_INT(emitter_pdg);
    BRANCH_INT(emitter_type); BRANCH_INT(production_mode);
    BRANCH_INT(geometry_id); BRANCH_INT(mcp_pdg); BRANCH_INT(mother_pdg);
    BRANCH_INT(mother_index); BRANCH_INT(accepted); BRANCH_INT(passed_geometry);
    BRANCH_INT(mother_status); BRANCH_INT(mother_is_primary_like);
    BRANCH_INT(mother_mother_pdg); BRANCH_INT(from_feeddown_flag);
    BRANCH_INT(id1); BRANCH_INT(id2);
    BRANCH_DOUBLE(mcp_mass_GeV); BRANCH_DOUBLE(event_weight);
    BRANCH_DOUBLE(x1); BRANCH_DOUBLE(x2); BRANCH_DOUBLE(q2_fac_GeV2);
    BRANCH_DOUBLE(mhat_GeV); BRANCH_DOUBLE(px_GeV); BRANCH_DOUBLE(py_GeV);
    BRANCH_DOUBLE(pz_GeV); BRANCH_DOUBLE(E_GeV); BRANCH_DOUBLE(p_GeV);
    BRANCH_DOUBLE(theta_rad); BRANCH_DOUBLE(theta_x_rad);
    BRANCH_DOUBLE(theta_y_rad); BRANCH_DOUBLE(eta); BRANCH_DOUBLE(phi_rad);
    BRANCH_DOUBLE(x_at_detector_m); BRANCH_DOUBLE(y_at_detector_m);
    BRANCH_DOUBLE(z_detector_m); BRANCH_DOUBLE(mother_px_GeV);
    BRANCH_DOUBLE(mother_py_GeV); BRANCH_DOUBLE(mother_pz_GeV);
    BRANCH_DOUBLE(mother_E_GeV); BRANCH_DOUBLE(mother_p_GeV);
    BRANCH_DOUBLE(mother_theta_rad); BRANCH_DOUBLE(pair_px_GeV);
    BRANCH_DOUBLE(pair_py_GeV); BRANCH_DOUBLE(pair_pz_GeV);
    BRANCH_DOUBLE(pair_E_GeV); BRANCH_DOUBLE(pair_mass_GeV);
    BRANCH_DOUBLE(pair_pT_GeV);
#undef BRANCH_INT
#undef BRANCH_DOUBLE
}

void writeOutput(const Options& options, const std::vector<ThreadResult>& results) {
    const std::filesystem::path outputPath(options.output);
    if (outputPath.has_parent_path()) {
        std::filesystem::create_directories(outputPath.parent_path());
    }

    TFile file(options.output.c_str(), "RECREATE");
    if (file.IsZombie()) {
        throw std::runtime_error("cannot create output: " + options.output);
    }

    TTree summary("mcp_summary", "MCP emitter summary");
    int run_id = options.jobId;
    int job_id = options.jobId;
    int thread_id = 0;
    int seed = 0;
    int emitter_pdg = 0;
    int parent_pdg = 0;
    int emitter_type = kDyEmitterType;
    int parent_type = kDyEmitterType;
    int production_mode = kDyProductionMode;
    int geometry_id = static_cast<int>(options.geometry);
    int is_kinematically_open = 1;
    int stop_reason = 0;
    int stopping_mode = 0;
    int spectra_prescale = options.spectraPrescale;
    int hard_only = options.hardOnly;
    char emitter_name[32] = "dy";
    char parent_name[32] = "dy";
    char production_mode_name[32] = "drell_yan";
    char geometry_name[32];
    char stopping_mode_name[32] = "fixed-events";
    char spectra_mode_name[32];
    std::snprintf(geometry_name, 32, "%s", geometryName(options.geometry).c_str());
    std::snprintf(
        spectra_mode_name,
        32,
        "%s",
        options.spectraMode == SpectraMode::All ? "all"
        : options.spectraMode == SpectraMode::Accepted ? "accepted" : "none"
    );

    double mcp_mass_GeV = options.mass;
    double mcp_mass = options.mass;
    double epsilon_gen = options.epsilon;
    double mhat_min_GeV = options.mHatMin;
    double emitter_per_event = 0.0;
    double parent_yield_per_event = 0.0;
    double acceptance_fraction = 0.0;
    double acceptance_uncertainty_binomial = 0.0;
    double accepted_mcp_per_event = 0.0;
    double sigma_gen_mb = 0.0;
    double sigma_err_mb = 0.0;
    double sigma_over_epsilon2_mb = 0.0;
    double accepted_mcp_sigma_mb = 0.0;
    double accepted_mcp_sigma_over_epsilon2_mb = 0.0;
    double weight_per_event_mb = 0.0;
    long long n_events_generated = 0;
    long long n_pythia_next_failures = 0;
    long long n_malformed_events = 0;
    long long n_emitter_record_entries = 0;
    long long n_emitter_decayed_to_mcp = 0;
    long long n_mcp_pairs = 0;
    long long n_mcp_pairing_anomalies = 0;
    long long n_emitter_total = 0;
    long long n_parent_total = 0;
    long long n_mcp_total = 0;
    long long n_mcp_accepted = 0;
    long long n_mcp_wrong_mother = 0;
    long long n_mcp_no_mother = 0;

#define BRANCH(name, address, leaf) summary.Branch(name, address, leaf)
    BRANCH("run_id", &run_id, "run_id/I");
    BRANCH("job_id", &job_id, "job_id/I");
    BRANCH("thread_id", &thread_id, "thread_id/I");
    BRANCH("seed", &seed, "seed/I");
    BRANCH("mcp_mass_GeV", &mcp_mass_GeV, "mcp_mass_GeV/D");
    BRANCH("mcp_mass", &mcp_mass, "mcp_mass/D");
    BRANCH("epsilon_gen", &epsilon_gen, "epsilon_gen/D");
    BRANCH("mhat_min_GeV", &mhat_min_GeV, "mhat_min_GeV/D");
    BRANCH("emitter_pdg", &emitter_pdg, "emitter_pdg/I");
    BRANCH("parent_pdg", &parent_pdg, "parent_pdg/I");
    BRANCH("emitter_name", emitter_name, "emitter_name/C");
    BRANCH("parent_name", parent_name, "parent_name/C");
    BRANCH("emitter_type", &emitter_type, "emitter_type/I");
    BRANCH("parent_type", &parent_type, "parent_type/I");
    BRANCH("production_mode", &production_mode, "production_mode/I");
    BRANCH("production_mode_name", production_mode_name, "production_mode_name/C");
    BRANCH("geometry_id", &geometry_id, "geometry_id/I");
    BRANCH("geometry_name", geometry_name, "geometry_name/C");
    BRANCH("is_kinematically_open", &is_kinematically_open, "is_kinematically_open/I");
    BRANCH("hard_only", &hard_only, "hard_only/I");
    BRANCH("n_events_generated", &n_events_generated, "n_events_generated/L");
    BRANCH("n_pythia_next_failures", &n_pythia_next_failures, "n_pythia_next_failures/L");
    BRANCH("n_malformed_events", &n_malformed_events, "n_malformed_events/L");
    BRANCH("n_emitter_record_entries", &n_emitter_record_entries, "n_emitter_record_entries/L");
    BRANCH("n_emitter_decayed_to_mcp", &n_emitter_decayed_to_mcp, "n_emitter_decayed_to_mcp/L");
    BRANCH("n_mcp_pairs", &n_mcp_pairs, "n_mcp_pairs/L");
    BRANCH("n_mcp_pairing_anomalies", &n_mcp_pairing_anomalies, "n_mcp_pairing_anomalies/L");
    BRANCH("n_emitter_total", &n_emitter_total, "n_emitter_total/L");
    BRANCH("n_parent_total", &n_parent_total, "n_parent_total/L");
    BRANCH("n_mcp_total", &n_mcp_total, "n_mcp_total/L");
    BRANCH("n_mcp_accepted", &n_mcp_accepted, "n_mcp_accepted/L");
    BRANCH("n_mcp_wrong_mother", &n_mcp_wrong_mother, "n_mcp_wrong_mother/L");
    BRANCH("n_mcp_no_mother", &n_mcp_no_mother, "n_mcp_no_mother/L");
    BRANCH("emitter_per_event", &emitter_per_event, "emitter_per_event/D");
    BRANCH("parent_yield_per_event", &parent_yield_per_event, "parent_yield_per_event/D");
    BRANCH("acceptance_fraction", &acceptance_fraction, "acceptance_fraction/D");
    BRANCH("acceptance_uncertainty_binomial", &acceptance_uncertainty_binomial, "acceptance_uncertainty_binomial/D");
    BRANCH("accepted_mcp_per_event", &accepted_mcp_per_event, "accepted_mcp_per_event/D");
    BRANCH("sigma_gen_mb", &sigma_gen_mb, "sigma_gen_mb/D");
    BRANCH("sigma_err_mb", &sigma_err_mb, "sigma_err_mb/D");
    BRANCH("sigma_over_epsilon2_mb", &sigma_over_epsilon2_mb, "sigma_over_epsilon2_mb/D");
    BRANCH("accepted_mcp_sigma_mb", &accepted_mcp_sigma_mb, "accepted_mcp_sigma_mb/D");
    BRANCH("accepted_mcp_sigma_over_epsilon2_mb", &accepted_mcp_sigma_over_epsilon2_mb, "accepted_mcp_sigma_over_epsilon2_mb/D");
    BRANCH("weight_per_event_mb", &weight_per_event_mb, "weight_per_event_mb/D");
    BRANCH("stop_reason", &stop_reason, "stop_reason/I");
    BRANCH("stopping_mode", &stopping_mode, "stopping_mode/I");
    BRANCH("stopping_mode_name", stopping_mode_name, "stopping_mode_name/C");
    BRANCH("spectra_mode_name", spectra_mode_name, "spectra_mode_name/C");
    BRANCH("spectra_prescale", &spectra_prescale, "spectra_prescale/I");
#undef BRANCH

    for (const auto& result : results) {
        thread_id = result.threadId;
        seed = result.seed;
        n_events_generated = result.nEvents;
        n_pythia_next_failures = result.nFailures;
        n_malformed_events = result.nMalformed;
        n_emitter_record_entries = result.nPairs;
        n_emitter_decayed_to_mcp = result.nPairs;
        n_mcp_pairs = result.nPairs;
        n_mcp_pairing_anomalies = result.nMalformed;
        n_emitter_total = result.nPairs;
        n_parent_total = result.nPairs;
        n_mcp_total = result.nMcp;
        n_mcp_accepted = result.nAccepted;
        emitter_per_event = n_events_generated > 0
            ? static_cast<double>(n_mcp_pairs) / n_events_generated : 0.0;
        parent_yield_per_event = emitter_per_event;
        acceptance_fraction = n_mcp_total > 0
            ? static_cast<double>(n_mcp_accepted) / n_mcp_total : 0.0;
        acceptance_uncertainty_binomial = n_mcp_total > 0
            ? std::sqrt(std::max(
                0.0,
                acceptance_fraction * (1.0 - acceptance_fraction) / n_mcp_total
              ))
            : 0.0;
        accepted_mcp_per_event = n_events_generated > 0
            ? static_cast<double>(n_mcp_accepted) / n_events_generated : 0.0;
        sigma_gen_mb = result.sigmaMb;
        sigma_err_mb = result.sigmaErrMb;
        sigma_over_epsilon2_mb = sigma_gen_mb
            / (options.epsilon * options.epsilon);
        accepted_mcp_sigma_mb = sigma_gen_mb * accepted_mcp_per_event;
        accepted_mcp_sigma_over_epsilon2_mb = sigma_over_epsilon2_mb
            * accepted_mcp_per_event;
        weight_per_event_mb = n_events_generated > 0
            ? sigma_gen_mb / n_events_generated : 0.0;
        summary.Fill();
    }

    TTree spectra("mcp_spectra", "MCP spectra");
    SpectraRow row;
    branchSpectra(spectra, row);
    for (const auto& result : results) {
        for (const auto& source : result.spectra) {
            row = source;
            spectra.Fill();
        }
    }

    TTree status("emitter_status_counts", "Emitter event-record status counts");
    int status_code = 0;
    long long status_count = 0;
    status.Branch("mcp_mass_GeV", &mcp_mass_GeV, "mcp_mass_GeV/D");
    status.Branch("emitter_pdg", &emitter_pdg, "emitter_pdg/I");
    status.Branch("production_mode", &production_mode, "production_mode/I");
    status.Branch("geometry_id", &geometry_id, "geometry_id/I");
    status.Branch("status_code", &status_code, "status_code/I");
    status.Branch("count", &status_count, "count/L");

    summary.Write();
    spectra.Write();
    status.Write();
    file.Close();
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const long long base = options.nEvents / options.nThreads;
        const long long remainder = options.nEvents % options.nThreads;
        std::vector<ThreadResult> results(options.nThreads);
        std::vector<std::thread> threads;

        if (options.nThreads == 1) {
            results[0] = runThread(options, 0, options.nEvents);
        } else {
            for (int threadId = 0; threadId < options.nThreads; ++threadId) {
                const long long target = base + (threadId < remainder ? 1 : 0);
                threads.emplace_back([&, threadId, target]() {
                    results[threadId] = runThread(options, threadId, target);
                });
            }
            for (auto& thread : threads) thread.join();
        }

        writeOutput(options, results);

        long long generated = 0;
        long long failures = 0;
        long long malformed = 0;
        long long pairs = 0;
        long long mcps = 0;
        long long accepted = 0;
        double sigmaSum = 0.0;
        double sigmaErr2 = 0.0;
        for (const auto& result : results) {
            generated += result.nEvents;
            failures += result.nFailures;
            malformed += result.nMalformed;
            pairs += result.nPairs;
            mcps += result.nMcp;
            accepted += result.nAccepted;
            sigmaSum += result.sigmaMb;
            sigmaErr2 += result.sigmaErrMb * result.sigmaErrMb;
        }

        const double sigmaMean = sigmaSum / options.nThreads;
        const double sigmaErrMean = std::sqrt(sigmaErr2) / options.nThreads;
        const double acceptedPerEvent = generated > 0
            ? static_cast<double>(accepted) / generated : 0.0;

        std::cout << std::setprecision(12)
                  << "model=dirac_fermion_photon_only\n"
                  << "production_mode=drell_yan\n"
                  << "events=" << generated << '\n'
                  << "pairs=" << pairs << '\n'
                  << "mcp=" << mcps << '\n'
                  << "accepted=" << accepted << '\n'
                  << "accepted_mcp_per_event=" << acceptedPerEvent << '\n'
                  << "pythia_next_failures=" << failures << '\n'
                  << "malformed_events=" << malformed << '\n'
                  << "sigma_gen_mb_mean=" << sigmaMean << '\n'
                  << "sigma_err_mb_mean=" << sigmaErrMean << '\n'
                  << "sigma_over_epsilon2_mb="
                  << sigmaMean / (options.epsilon * options.epsilon) << '\n'
                  << "accepted_mcp_sigma_over_epsilon2_mb="
                  << sigmaMean / (options.epsilon * options.epsilon)
                   * acceptedPerEvent << '\n'
                  << "geometry=" << geometryName(options.geometry) << '\n'
                  << "hard_only=" << options.hardOnly << '\n'
                  << "output=" << options.output << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 2;
    }
}

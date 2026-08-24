#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kAlpha = 1.0 / 137.0;
constexpr double kProtonMass = 0.9382720813;
constexpr double kLambdaQCD = 0.25;
constexpr int kMcpPdg = 1000222;
constexpr int kVirtualPhotonPdg = 22;
constexpr int kProductionModePB = 3;

enum class GeometryId { ArgoNeuT = 0, TwoByTwo = 1 };
enum class SpectraMode { All, Accepted, None };

struct Options {
    int seed = 12345;
    int n_threads = 1;
    long long n_trials = 100000;
    double mcp_mass = 0.1;
    GeometryId geometry = GeometryId::TwoByTwo;
    std::string emitter = "pb";
    std::string production_mode = "proton_bremsstrahlung";
    std::string output = "pb.root";
    double beam_energy = 120.0;
    double epsilon = 0.01;
    double sigma_inelastic_mb = 38.4538;
    double interaction_ratio = 1.0;
    double m_pair_max = 12.0;
    SpectraMode spectra_mode = SpectraMode::Accepted;
    int spectra_prescale = 1;
    int job_id = 0;
    bool quiet = false;
};

struct Resonance {
    double mass;
    double width;
    double coefficient;
};

constexpr std::array<Resonance, 6> kResonances{{
    {0.775, 0.1474, 0.616},
    {1.464, 0.4000, 0.223},
    {1.570, 0.1440, -0.339},
    {0.783, 0.0086, 1.011},
    {1.410, 0.2900, -0.881},
    {1.670, 0.3150, 0.369},
}};

struct Projection {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    int passed = 0;
};

struct SpectraRow {
    int event_index = 0;
    int seed = 0;
    int mcp_pdg = 0;
    int geometry_id = 0;
    int production_mode = kProductionModePB;
    int passed_geometry = 0;
    double mcp_mass_GeV = 0.0;
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
    double pair_mass_GeV = 0.0;
    double pair_E_GeV = 0.0;
    double pair_pt_GeV = 0.0;
    double z_fraction = 0.0;
    double qmin2_GeV2 = 0.0;
    double mc_weight_over_epsilon2 = 0.0;
    double probability_weight_over_epsilon2 = 0.0;
    double sigma_weight_over_epsilon2_mb = 0.0;
};

double uniformOpen(std::mt19937_64& rng) {
    static std::uniform_real_distribution<double> distribution(
        std::nextafter(0.0, 1.0), 1.0
    );
    return distribution(rng);
}

GeometryId parseGeometry(const std::string& value) {
    if (value == "2x2" || value == "twobytwo") return GeometryId::TwoByTwo;
    if (value == "argoneut") return GeometryId::ArgoNeuT;
    throw std::runtime_error("unknown geometry: " + value);
}

const char* geometryName(GeometryId geometry) {
    return geometry == GeometryId::TwoByTwo ? "2x2" : "argoneut";
}

SpectraMode parseSpectraMode(const std::string& value) {
    if (value == "all") return SpectraMode::All;
    if (value == "accepted") return SpectraMode::Accepted;
    if (value == "none") return SpectraMode::None;
    throw std::runtime_error("unknown spectra mode: " + value);
}

const char* spectraModeName(SpectraMode mode) {
    if (mode == SpectraMode::All) return "all";
    if (mode == SpectraMode::Accepted) return "accepted";
    return "none";
}

void printUsage(const char* executable) {
    std::cout
        << "Usage: " << executable
        << " <seed> <nThreads> <nTrials> <mcpMassGeV> <geometry>"
        << " <emitterName> <productionMode> <outputFile> [options]\n"
        << "Options:\n"
        << "  --beam-energy GEV             Beam total energy (default: 120)\n"
        << "  --epsilon VALUE               Reference MCP charge (default: 0.01)\n"
        << "  --sigma-inelastic-mb VALUE    Reference interaction cross section\n"
        << "  --interaction-ratio VALUE     sigma(s') / sigma(s) (default: 1)\n"
        << "  --m-pair-max GEV              Maximum virtual-photon mass\n"
        << "  --write-spectra MODE          all, accepted, or none\n"
        << "  --spectra-prescale N          Retain one of every N selected rows\n"
        << "  --job-id N                    Job identifier stored in ROOT\n"
        << "  --quiet, --batch              Compact terminal output\n"
        << "  --help                        Show this message\n";
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
    options.n_threads = std::stoi(argv[2]);
    options.n_trials = std::stoll(argv[3]);
    options.mcp_mass = std::stod(argv[4]);
    options.geometry = parseGeometry(argv[5]);
    options.emitter = argv[6];
    options.production_mode = argv[7];
    options.output = argv[8];

    for (int i = 9; i < argc; ++i) {
        const std::string key = argv[i];
        const auto value = [&]() {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + key);
            return std::string(argv[++i]);
        };

        if (key == "--beam-energy") options.beam_energy = std::stod(value());
        else if (key == "--epsilon") options.epsilon = std::stod(value());
        else if (key == "--sigma-inelastic-mb") {
            options.sigma_inelastic_mb = std::stod(value());
        } else if (key == "--interaction-ratio") {
            options.interaction_ratio = std::stod(value());
        } else if (key == "--m-pair-max") {
            options.m_pair_max = std::stod(value());
        } else if (key == "--write-spectra") {
            options.spectra_mode = parseSpectraMode(value());
        } else if (key == "--spectra-prescale") {
            options.spectra_prescale = std::stoi(value());
        } else if (key == "--job-id") {
            options.job_id = std::stoi(value());
        } else if (key == "--quiet" || key == "--batch") {
            options.quiet = true;
        } else if (key == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown option: " + key);
        }
    }

    if (options.n_threads != 1) {
        throw std::runtime_error("the standalone PB generator currently requires nThreads=1");
    }
    if (options.n_trials <= 0 || options.mcp_mass <= 0.0 ||
        options.beam_energy <= kProtonMass || options.epsilon <= 0.0 ||
        options.sigma_inelastic_mb <= 0.0 || options.interaction_ratio <= 0.0 ||
        options.spectra_prescale <= 0) {
        throw std::runtime_error("invalid non-positive option");
    }
    if (options.emitter != "pb" || options.production_mode != "proton_bremsstrahlung") {
        throw std::runtime_error(
            "this executable requires emitterName=pb and productionMode=proton_bremsstrahlung"
        );
    }
    if (options.m_pair_max <= 2.0 * options.mcp_mass) {
        throw std::runtime_error("m-pair-max must exceed 2*mcpMass");
    }
    return options;
}

Projection project(double px, double py, double pz, GeometryId geometry) {
    Projection result;
    if (pz <= 0.0) return result;

    if (geometry == GeometryId::ArgoNeuT) {
        result.z = 1000.0;
        result.x = px / pz * result.z;
        result.y = py / pz * result.z;
        result.passed = (
            result.x >= -0.24 && result.x <= 0.24 &&
            result.y >= -0.20 && result.y <= 0.20
        );
        return result;
    }

    result.z = 1040.0;
    result.x = px / pz * result.z;
    result.y = py / pz * result.z;
    const bool in_x = (
        (result.x >= -0.65 && result.x <= -0.05) ||
        (result.x >= 0.05 && result.x <= 0.65)
    );
    result.passed = in_x && result.y >= -0.70 && result.y <= 0.70;
    return result;
}

std::complex<double> protonFormFactor(double mass2) {
    std::complex<double> result(0.0, 0.0);
    for (const Resonance& resonance : kResonances) {
        const double m2 = resonance.mass * resonance.mass;
        const std::complex<double> denominator(
            m2 - mass2,
            -resonance.mass * resonance.width
        );
        result += resonance.coefficient * m2 / denominator;
    }
    return result;
}

double splittingKernel(double z, double pt2, double mass2) {
    const double zp = 1.0 - z;
    const double mp2 = kProtonMass * kProtonMass;
    const double h = pt2 + zp * mass2 + z * z * mp2;
    const double h2 = h * h;
    const double term1 = (1.0 + zp * zp) / z;
    const double term2 = -2.0 * z * zp * (
        (2.0 * mp2 + mass2) / h -
        2.0 * mp2 * mp2 * z * z / h2
    );
    const double term3 = 2.0 * z * zp * (z + zp * zp) * mp2 * mass2 / h2;
    const double term4 = 2.0 * z * zp * zp * mass2 * mass2 / h2;
    return kAlpha / (2.0 * kPi * h) * (term1 + term2 + term3 + term4);
}

double conversionKernel(double mcp_mass, double mass2) {
    const double ratio = mcp_mass * mcp_mass / mass2;
    const double beta2 = 1.0 - 4.0 * ratio;
    if (beta2 <= 0.0) return 0.0;
    return kAlpha / (3.0 * kPi * mass2) *
           std::sqrt(beta2) * (1.0 + 2.0 * ratio);
}

double truncatedCauchyPdf(
    double value,
    const Resonance& resonance,
    double minimum,
    double maximum
) {
    const double center = resonance.mass * resonance.mass;
    const double scale = resonance.mass * resonance.width;
    const double lower = std::atan((minimum - center) / scale);
    const double upper = std::atan((maximum - center) / scale);
    const double delta = value - center;
    return (scale / (delta * delta + scale * scale)) / (upper - lower);
}

double sampleTruncatedCauchy(
    std::mt19937_64& rng,
    const Resonance& resonance,
    double minimum,
    double maximum
) {
    const double center = resonance.mass * resonance.mass;
    const double scale = resonance.mass * resonance.width;
    const double lower = std::atan((minimum - center) / scale);
    const double upper = std::atan((maximum - center) / scale);
    return center + scale * std::tan(lower + uniformOpen(rng) * (upper - lower));
}

struct Sample {
    double value = 0.0;
    double density = 0.0;
};

Sample sampleZ(std::mt19937_64& rng) {
    constexpr double minimum = 1.0e-4;
    constexpr double maximum = 1.0 - minimum;
    const double log_range = std::log(maximum / minimum);
    double value = 0.0;
    if (uniformOpen(rng) < 0.5) {
        value = minimum * std::exp(uniformOpen(rng) * log_range);
    } else {
        const double complement = minimum * std::exp(uniformOpen(rng) * log_range);
        value = 1.0 - complement;
    }
    const double density = 0.5 / (value * log_range) +
                           0.5 / ((1.0 - value) * log_range);
    return {value, density};
}

Sample samplePt2(std::mt19937_64& rng) {
    constexpr double minimum = 1.0e-10;
    constexpr double maximum = 1.0;
    const double log_range = std::log(maximum / minimum);
    const double value = minimum * std::exp(uniformOpen(rng) * log_range);
    return {value, 1.0 / (value * log_range)};
}

Sample sampleMass2(
    std::mt19937_64& rng,
    double minimum,
    double maximum
) {
    constexpr double log_fraction = 0.5;
    const double log_range = std::log(maximum / minimum);
    double value = 0.0;

    if (uniformOpen(rng) < log_fraction) {
        value = minimum * std::exp(uniformOpen(rng) * log_range);
    } else {
        const std::size_t index = std::min<std::size_t>(
            static_cast<std::size_t>(uniformOpen(rng) * kResonances.size()),
            kResonances.size() - 1
        );
        value = sampleTruncatedCauchy(
            rng, kResonances[index], minimum, maximum
        );
    }

    const double log_density = 1.0 / (value * log_range);
    double resonance_density = 0.0;
    for (const Resonance& resonance : kResonances) {
        resonance_density += truncatedCauchyPdf(
            value, resonance, minimum, maximum
        );
    }
    resonance_density /= static_cast<double>(kResonances.size());
    return {
        value,
        log_fraction * log_density + (1.0 - log_fraction) * resonance_density
    };
}

double standardError(double sum, double sum2, long long count) {
    if (count < 2) return 0.0;
    const double n = static_cast<double>(count);
    const double numerator = std::max(0.0, sum2 - sum * sum / n);
    return std::sqrt(numerator / (n * (n - 1.0)));
}

double weightedAcceptanceError(
    double sum_pair,
    double sum_pair2,
    double sum_accepted,
    double sum_accepted2,
    double sum_pair_accepted,
    long long count
) {
    if (count < 2 || sum_pair <= 0.0) return 0.0;
    const double n = static_cast<double>(count);
    const double mean_x = sum_pair / n;
    const double mean_y = sum_accepted / n;
    const double var_mean_x = std::max(0.0, sum_pair2 - n * mean_x * mean_x) /
                              (n * (n - 1.0));
    const double var_mean_y = std::max(0.0, sum_accepted2 - n * mean_y * mean_y) /
                              (n * (n - 1.0));
    const double cov_mean = (sum_pair_accepted - n * mean_x * mean_y) /
                            (n * (n - 1.0));
    const double grad_x = -mean_y / (2.0 * mean_x * mean_x);
    const double grad_y = 1.0 / (2.0 * mean_x);
    const double variance = grad_x * grad_x * var_mean_x +
                            grad_y * grad_y * var_mean_y +
                            2.0 * grad_x * grad_y * cov_mean;
    return std::sqrt(std::max(0.0, variance));
}

double pseudorapidity(const TLorentzVector& particle) {
    const double momentum = particle.P();
    if (momentum <= std::abs(particle.Pz())) {
        return particle.Pz() >= 0.0
            ? std::numeric_limits<double>::infinity()
            : -std::numeric_limits<double>::infinity();
    }
    return 0.5 * std::log(
        (momentum + particle.Pz()) / (momentum - particle.Pz())
    );
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const double pair_mass2_min = 4.0 * options.mcp_mass * options.mcp_mass;
        const double pair_mass2_max = options.m_pair_max * options.m_pair_max;
        const double beam_momentum = std::sqrt(
            options.beam_energy * options.beam_energy -
            kProtonMass * kProtonMass
        );
        const double total_s = 2.0 * kProtonMass * kProtonMass +
                               2.0 * kProtonMass * options.beam_energy;

        const std::filesystem::path output_path(options.output);
        if (output_path.has_parent_path()) {
            std::filesystem::create_directories(output_path.parent_path());
        }

        TFile output(options.output.c_str(), "RECREATE");
        if (output.IsZombie()) {
            throw std::runtime_error("failed to create output ROOT file");
        }

        TTree spectra("mcp_spectra", "Weighted MCP spectra from proton bremsstrahlung");
        SpectraRow row;
        spectra.Branch("event_index", &row.event_index, "event_index/I");
        spectra.Branch("seed", &row.seed, "seed/I");
        spectra.Branch("mcp_pdg", &row.mcp_pdg, "mcp_pdg/I");
        spectra.Branch("geometry_id", &row.geometry_id, "geometry_id/I");
        spectra.Branch("production_mode", &row.production_mode, "production_mode/I");
        spectra.Branch("passed_geometry", &row.passed_geometry, "passed_geometry/I");
        spectra.Branch("mcp_mass_GeV", &row.mcp_mass_GeV, "mcp_mass_GeV/D");
        spectra.Branch("px_GeV", &row.px_GeV, "px_GeV/D");
        spectra.Branch("py_GeV", &row.py_GeV, "py_GeV/D");
        spectra.Branch("pz_GeV", &row.pz_GeV, "pz_GeV/D");
        spectra.Branch("E_GeV", &row.E_GeV, "E_GeV/D");
        spectra.Branch("p_GeV", &row.p_GeV, "p_GeV/D");
        spectra.Branch("theta_rad", &row.theta_rad, "theta_rad/D");
        spectra.Branch("theta_x_rad", &row.theta_x_rad, "theta_x_rad/D");
        spectra.Branch("theta_y_rad", &row.theta_y_rad, "theta_y_rad/D");
        spectra.Branch("eta", &row.eta, "eta/D");
        spectra.Branch("phi_rad", &row.phi_rad, "phi_rad/D");
        spectra.Branch("x_at_detector_m", &row.x_at_detector_m, "x_at_detector_m/D");
        spectra.Branch("y_at_detector_m", &row.y_at_detector_m, "y_at_detector_m/D");
        spectra.Branch("z_detector_m", &row.z_detector_m, "z_detector_m/D");
        spectra.Branch("pair_mass_GeV", &row.pair_mass_GeV, "pair_mass_GeV/D");
        spectra.Branch("pair_E_GeV", &row.pair_E_GeV, "pair_E_GeV/D");
        spectra.Branch("pair_pt_GeV", &row.pair_pt_GeV, "pair_pt_GeV/D");
        spectra.Branch("z_fraction", &row.z_fraction, "z_fraction/D");
        spectra.Branch("qmin2_GeV2", &row.qmin2_GeV2, "qmin2_GeV2/D");
        spectra.Branch(
            "mc_weight_over_epsilon2",
            &row.mc_weight_over_epsilon2,
            "mc_weight_over_epsilon2/D"
        );
        spectra.Branch(
            "probability_weight_over_epsilon2",
            &row.probability_weight_over_epsilon2,
            "probability_weight_over_epsilon2/D"
        );
        spectra.Branch(
            "sigma_weight_over_epsilon2_mb",
            &row.sigma_weight_over_epsilon2_mb,
            "sigma_weight_over_epsilon2_mb/D"
        );

        std::mt19937_64 rng(static_cast<std::uint64_t>(options.seed));
        long long n_valid = 0;
        long long n_mcp_accepted = 0;
        long long n_spectra_candidates = 0;
        long long n_spectra_written = 0;
        long long cut_pt_relative = 0;
        long long cut_pt_absolute = 0;
        long long cut_qmin = 0;
        long long cut_energy = 0;
        long long cut_sprime = 0;
        long long cut_kernel = 0;
        double sum_pair_weights = 0.0;
        double sum_pair_weights2 = 0.0;
        double sum_accepted_weights = 0.0;
        double sum_accepted_weights2 = 0.0;
        double sum_pair_accepted_weights = 0.0;

        for (long long trial = 0; trial < options.n_trials; ++trial) {
            const Sample z_sample = sampleZ(rng);
            const Sample pt2_sample = samplePt2(rng);
            const Sample mass2_sample = sampleMass2(
                rng, pair_mass2_min, pair_mass2_max
            );
            const double z = z_sample.value;
            const double zp = 1.0 - z;
            const double pt2 = pt2_sample.value;
            const double pt = std::sqrt(pt2);
            const double mass2 = mass2_sample.value;
            const double pair_mass = std::sqrt(mass2);
            const double kz = z * beam_momentum;
            const double pair_momentum2 = kz * kz + pt2;
            const double pair_energy = std::sqrt(pair_momentum2 + mass2);

            if (pt >= 0.1 * pair_energy) {
                ++cut_pt_relative;
                continue;
            }
            if (pt >= 1.0) {
                ++cut_pt_absolute;
                continue;
            }

            const double h = pt2 + zp * mass2 +
                             z * z * kProtonMass * kProtonMass;
            const double qmin2 = h * h /
                (4.0 * options.beam_energy * options.beam_energy *
                 z * z * zp * zp);
            if (!(qmin2 < kLambdaQCD * kLambdaQCD)) {
                ++cut_qmin;
                continue;
            }

            if (!(pair_energy > 5.0 * kProtonMass &&
                  options.beam_energy - pair_energy > 5.0 * kProtonMass &&
                  pair_energy > 5.0 * pair_mass &&
                  options.beam_energy - pair_energy > 5.0 * pair_mass)) {
                ++cut_energy;
                continue;
            }

            const double azimuth = 2.0 * kPi * uniformOpen(rng);
            const double kx = pt * std::cos(azimuth);
            const double ky = pt * std::sin(azimuth);
            const double residual_energy = options.beam_energy + kProtonMass - pair_energy;
            const double residual_px = -kx;
            const double residual_py = -ky;
            const double residual_pz = beam_momentum - kz;
            const double sprime = residual_energy * residual_energy -
                                  residual_px * residual_px -
                                  residual_py * residual_py -
                                  residual_pz * residual_pz;
            if (!(sprime > 4.0 * kProtonMass * kProtonMass && sprime < total_s)) {
                ++cut_sprime;
                continue;
            }

            const double fww = splittingKernel(z, pt2, mass2);
            const double conversion = conversionKernel(options.mcp_mass, mass2);
            const double form_factor2 = std::norm(protonFormFactor(mass2));
            const double proposal = z_sample.density *
                                    pt2_sample.density *
                                    mass2_sample.density;
            const double weight = options.interaction_ratio *
                                  fww * form_factor2 * conversion / proposal;
            if (!(std::isfinite(weight) && weight > 0.0)) {
                ++cut_kernel;
                continue;
            }

            ++n_valid;
            const double rest_energy = 0.5 * pair_mass;
            const double rest_momentum = std::sqrt(std::max(
                0.0,
                rest_energy * rest_energy - options.mcp_mass * options.mcp_mass
            ));
            const double cos_theta = 2.0 * uniformOpen(rng) - 1.0;
            const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
            const double decay_phi = 2.0 * kPi * uniformOpen(rng);
            TLorentzVector chi(
                rest_momentum * sin_theta * std::cos(decay_phi),
                rest_momentum * sin_theta * std::sin(decay_phi),
                rest_momentum * cos_theta,
                rest_energy
            );
            TLorentzVector chibar(-chi.Px(), -chi.Py(), -chi.Pz(), rest_energy);
            const TVector3 boost(kx / pair_energy, ky / pair_energy, kz / pair_energy);
            chi.Boost(boost);
            chibar.Boost(boost);

            const std::array<TLorentzVector, 2> particles{{chi, chibar}};
            const std::array<int, 2> pdgs{{kMcpPdg, -kMcpPdg}};
            int accepted_in_pair = 0;
            std::array<Projection, 2> projections;
            for (std::size_t index = 0; index < particles.size(); ++index) {
                projections[index] = project(
                    particles[index].Px(),
                    particles[index].Py(),
                    particles[index].Pz(),
                    options.geometry
                );
                accepted_in_pair += projections[index].passed;
            }
            n_mcp_accepted += accepted_in_pair;

            const double accepted_weight = weight * accepted_in_pair;
            sum_pair_weights += weight;
            sum_pair_weights2 += weight * weight;
            sum_accepted_weights += accepted_weight;
            sum_accepted_weights2 += accepted_weight * accepted_weight;
            sum_pair_accepted_weights += weight * accepted_weight;

            if (options.spectra_mode == SpectraMode::None) continue;
            for (std::size_t index = 0; index < particles.size(); ++index) {
                const bool selected = options.spectra_mode == SpectraMode::All ||
                                      projections[index].passed;
                if (!selected) continue;
                ++n_spectra_candidates;
                if ((n_spectra_candidates - 1) % options.spectra_prescale != 0) continue;

                const TLorentzVector& particle = particles[index];
                const Projection& projection = projections[index];
                row.event_index = static_cast<int>(trial);
                row.seed = options.seed;
                row.mcp_pdg = pdgs[index];
                row.geometry_id = static_cast<int>(options.geometry);
                row.production_mode = kProductionModePB;
                row.passed_geometry = projection.passed;
                row.mcp_mass_GeV = options.mcp_mass;
                row.px_GeV = particle.Px();
                row.py_GeV = particle.Py();
                row.pz_GeV = particle.Pz();
                row.E_GeV = particle.E();
                row.p_GeV = particle.P();
                row.theta_rad = particle.Theta();
                row.theta_x_rad = particle.Pz() != 0.0
                    ? particle.Px() / particle.Pz() : 0.0;
                row.theta_y_rad = particle.Pz() != 0.0
                    ? particle.Py() / particle.Pz() : 0.0;
                row.eta = pseudorapidity(particle);
                row.phi_rad = particle.Phi();
                row.x_at_detector_m = projection.x;
                row.y_at_detector_m = projection.y;
                row.z_detector_m = projection.z;
                row.pair_mass_GeV = pair_mass;
                row.pair_E_GeV = pair_energy;
                row.pair_pt_GeV = pt;
                row.z_fraction = z;
                row.qmin2_GeV2 = qmin2;
                row.mc_weight_over_epsilon2 = weight;
                row.probability_weight_over_epsilon2 =
                    weight * options.spectra_prescale /
                    static_cast<double>(options.n_trials);
                row.sigma_weight_over_epsilon2_mb =
                    row.probability_weight_over_epsilon2 * options.sigma_inelastic_mb;
                spectra.Fill();
                ++n_spectra_written;
            }
        }

        const double n = static_cast<double>(options.n_trials);
        double probability_pair_over_epsilon2 = sum_pair_weights / n;
        double probability_pair_error_over_epsilon2 = standardError(
            sum_pair_weights, sum_pair_weights2, options.n_trials
        );
        double accepted_mcp_probability_over_epsilon2 =
            sum_accepted_weights / n;
        double accepted_mcp_probability_error_over_epsilon2 = standardError(
            sum_accepted_weights, sum_accepted_weights2, options.n_trials
        );
        double weighted_acceptance = sum_pair_weights > 0.0
            ? sum_accepted_weights / (2.0 * sum_pair_weights)
            : 0.0;
        double weighted_acceptance_error = weightedAcceptanceError(
            sum_pair_weights,
            sum_pair_weights2,
            sum_accepted_weights,
            sum_accepted_weights2,
            sum_pair_accepted_weights,
            options.n_trials
        );
        double sigma_over_epsilon2_mb =
            options.sigma_inelastic_mb * probability_pair_over_epsilon2;
        double sigma_error_over_epsilon2_mb =
            options.sigma_inelastic_mb * probability_pair_error_over_epsilon2;
        double accepted_mcp_sigma_over_epsilon2_mb =
            options.sigma_inelastic_mb * accepted_mcp_probability_over_epsilon2;
        double accepted_mcp_sigma_error_over_epsilon2_mb =
            options.sigma_inelastic_mb * accepted_mcp_probability_error_over_epsilon2;
        double sigma_gen_mb = sigma_over_epsilon2_mb *
                                    options.epsilon * options.epsilon;
        double sigma_err_mb = sigma_error_over_epsilon2_mb *
                                    options.epsilon * options.epsilon;
        double effective_sample_size = sum_pair_weights2 > 0.0
            ? sum_pair_weights * sum_pair_weights / sum_pair_weights2
            : 0.0;

        TTree summary("mcp_summary", "FWW proton-bremsstrahlung MCP summary");
        int run_id = options.job_id;
        int job_id = options.job_id;
        int seed = options.seed;
        int emitter_pdg = kVirtualPhotonPdg;
        int production_mode = kProductionModePB;
        int geometry_id = static_cast<int>(options.geometry);
        int spectra_prescale = options.spectra_prescale;
        long long n_events_generated = options.n_trials;
        long long n_phase_space_valid = n_valid;
        long long n_mcp_pairs = n_valid;
        long long n_mcp_total = 2 * n_valid;
        char emitter_name[32] = "pb";
        char production_mode_name[32] = "proton_bremsstrahlung";
        char geometry_name[16];
        char model_name[64] = "conventional_fww_full_vmd";
        char spectra_mode_name[16];
        std::snprintf(
            geometry_name, sizeof(geometry_name), "%s", ::geometryName(options.geometry)
        );
        std::snprintf(
            spectra_mode_name,
            sizeof(spectra_mode_name),
            "%s",
            ::spectraModeName(options.spectra_mode)
        );
        double mcp_mass_GeV = options.mcp_mass;
        double epsilon_gen = options.epsilon;
        double beam_energy_GeV = options.beam_energy;
        double sqrt_s_GeV = std::sqrt(total_s);
        double sigma_inelastic_mb = options.sigma_inelastic_mb;
        double interaction_cross_section_ratio = options.interaction_ratio;
        double m_pair_max_GeV = options.m_pair_max;
        double acceptance_fraction = weighted_acceptance;
        double acceptance_uncertainty_mc = weighted_acceptance_error;

#define BRANCH(name, address, leaf) summary.Branch(name, address, leaf)
        BRANCH("run_id", &run_id, "run_id/I");
        BRANCH("job_id", &job_id, "job_id/I");
        BRANCH("seed", &seed, "seed/I");
        BRANCH("emitter_pdg", &emitter_pdg, "emitter_pdg/I");
        BRANCH("production_mode", &production_mode, "production_mode/I");
        BRANCH("geometry_id", &geometry_id, "geometry_id/I");
        BRANCH("spectra_prescale", &spectra_prescale, "spectra_prescale/I");
        BRANCH("emitter_name", emitter_name, "emitter_name/C");
        BRANCH("production_mode_name", production_mode_name, "production_mode_name/C");
        BRANCH("geometry_name", geometry_name, "geometry_name/C");
        BRANCH("model_name", model_name, "model_name/C");
        BRANCH("spectra_mode_name", spectra_mode_name, "spectra_mode_name/C");
        BRANCH("mcp_mass_GeV", &mcp_mass_GeV, "mcp_mass_GeV/D");
        BRANCH("epsilon_gen", &epsilon_gen, "epsilon_gen/D");
        BRANCH("beam_energy_GeV", &beam_energy_GeV, "beam_energy_GeV/D");
        BRANCH("sqrt_s_GeV", &sqrt_s_GeV, "sqrt_s_GeV/D");
        BRANCH("sigma_inelastic_mb", &sigma_inelastic_mb, "sigma_inelastic_mb/D");
        BRANCH(
            "interaction_cross_section_ratio",
            &interaction_cross_section_ratio,
            "interaction_cross_section_ratio/D"
        );
        BRANCH("m_pair_max_GeV", &m_pair_max_GeV, "m_pair_max_GeV/D");
        BRANCH("n_events_generated", &n_events_generated, "n_events_generated/L");
        BRANCH("n_phase_space_valid", &n_phase_space_valid, "n_phase_space_valid/L");
        BRANCH("n_mcp_pairs", &n_mcp_pairs, "n_mcp_pairs/L");
        BRANCH("n_mcp_total", &n_mcp_total, "n_mcp_total/L");
        BRANCH("n_mcp_accepted", &n_mcp_accepted, "n_mcp_accepted/L");
        BRANCH("n_spectra_written", &n_spectra_written, "n_spectra_written/L");
        BRANCH("cut_pt_relative", &cut_pt_relative, "cut_pt_relative/L");
        BRANCH("cut_pt_absolute", &cut_pt_absolute, "cut_pt_absolute/L");
        BRANCH("cut_qmin", &cut_qmin, "cut_qmin/L");
        BRANCH("cut_energy", &cut_energy, "cut_energy/L");
        BRANCH("cut_sprime", &cut_sprime, "cut_sprime/L");
        BRANCH("cut_kernel", &cut_kernel, "cut_kernel/L");
        BRANCH("sum_pair_weights", &sum_pair_weights, "sum_pair_weights/D");
        BRANCH("sum_pair_weights2", &sum_pair_weights2, "sum_pair_weights2/D");
        BRANCH("sum_accepted_weights", &sum_accepted_weights, "sum_accepted_weights/D");
        BRANCH("sum_accepted_weights2", &sum_accepted_weights2, "sum_accepted_weights2/D");
        BRANCH(
            "sum_pair_accepted_weights",
            &sum_pair_accepted_weights,
            "sum_pair_accepted_weights/D"
        );
        BRANCH(
            "probability_pair_over_epsilon2",
            &probability_pair_over_epsilon2,
            "probability_pair_over_epsilon2/D"
        );
        BRANCH(
            "probability_pair_error_over_epsilon2",
            &probability_pair_error_over_epsilon2,
            "probability_pair_error_over_epsilon2/D"
        );
        BRANCH(
            "accepted_mcp_probability_over_epsilon2",
            &accepted_mcp_probability_over_epsilon2,
            "accepted_mcp_probability_over_epsilon2/D"
        );
        BRANCH(
            "accepted_mcp_probability_error_over_epsilon2",
            &accepted_mcp_probability_error_over_epsilon2,
            "accepted_mcp_probability_error_over_epsilon2/D"
        );
        BRANCH("acceptance_fraction", &acceptance_fraction, "acceptance_fraction/D");
        BRANCH(
            "acceptance_uncertainty_mc",
            &acceptance_uncertainty_mc,
            "acceptance_uncertainty_mc/D"
        );
        BRANCH("sigma_gen_mb", &sigma_gen_mb, "sigma_gen_mb/D");
        BRANCH("sigma_err_mb", &sigma_err_mb, "sigma_err_mb/D");
        BRANCH(
            "sigma_over_epsilon2_mb",
            &sigma_over_epsilon2_mb,
            "sigma_over_epsilon2_mb/D"
        );
        BRANCH(
            "sigma_error_over_epsilon2_mb",
            &sigma_error_over_epsilon2_mb,
            "sigma_error_over_epsilon2_mb/D"
        );
        BRANCH(
            "accepted_mcp_sigma_over_epsilon2_mb",
            &accepted_mcp_sigma_over_epsilon2_mb,
            "accepted_mcp_sigma_over_epsilon2_mb/D"
        );
        BRANCH(
            "accepted_mcp_sigma_error_over_epsilon2_mb",
            &accepted_mcp_sigma_error_over_epsilon2_mb,
            "accepted_mcp_sigma_error_over_epsilon2_mb/D"
        );
        BRANCH(
            "effective_sample_size",
            &effective_sample_size,
            "effective_sample_size/D"
        );
#undef BRANCH

        summary.Fill();
        output.cd();
        summary.Write();
        spectra.Write();
        output.Close();

        std::cout << std::setprecision(12)
                  << "model=conventional_fww_full_vmd\n"
                  << "mcp_mass_GeV=" << options.mcp_mass << "\n"
                  << "trials=" << options.n_trials << "\n"
                  << "valid_phase_space=" << n_valid << "\n"
                  << "sigma_over_epsilon2_mb=" << sigma_over_epsilon2_mb << "\n"
                  << "sigma_error_over_epsilon2_mb=" << sigma_error_over_epsilon2_mb << "\n"
                  << "accepted_mcp_sigma_over_epsilon2_mb="
                  << accepted_mcp_sigma_over_epsilon2_mb << "\n"
                  << "weighted_geometric_acceptance=" << weighted_acceptance << "\n"
                  << "effective_sample_size=" << effective_sample_size << "\n"
                  << "output=" << options.output << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << "\n";
        return 1;
    }
}

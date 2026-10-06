#include <iostream>

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h>
#include <AxilrodTellerMutoFunctor2.h>
#include <AxilrodTellerMutoFunctor3.h>
#include "benchmark/benchmark.h"

#include <autopas/cells/FullParticleCell.h>
#include <random>
#include <CLI/CLI.hpp>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <optional>
#include <type_traits>

#ifdef ENABLE_ITT
#include <ittnotify.h>

struct ITTScope {
    inline ITTScope() { __itt_resume(); }
    inline ~ITTScope() { __itt_pause(); }
};
#else
struct ITTScope {
    inline ITTScope() = default;
    inline ~ITTScope() = default;
};
#endif

// type aliases for ease of use
using Particle = mdLib::MoleculeLJ;
using Cell = autopas::FullParticleCell<Particle>;

// some constants that define the benchmark
constexpr bool mixing{false};
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool globals{false};
using ATM = mdLib::AxilrodTellerMutoFunctor<Particle, mixing, functorN3Modes, globals>;
using ATMGlobals = mdLib::AxilrodTellerMutoFunctor<Particle, mixing, functorN3Modes, true>;
using ATM2 = mdLib::AxilrodTellerMutoFunctor2<Particle, mixing, functorN3Modes, globals>;
using ATM2Globals = mdLib::AxilrodTellerMutoFunctor2<Particle, mixing, functorN3Modes, true>;
using ATM3 = mdLib::AxilrodTellerMutoFunctor3<Particle, mixing, functorN3Modes, globals>;
using ATM3Globals = mdLib::AxilrodTellerMutoFunctor3<Particle, mixing, functorN3Modes, true>;

struct BenchmarkConfig {
    int64_t minParticles = 1;
    int64_t maxParticles = 512;
    std::vector<int64_t> particles = {};
    double cellSize = 3;
    double cutoff = 3;
    size_t cellPoolSize = 1000;
    std::vector<std::string> targetFunctors = {"all"};
    std::vector<std::string> targetKernels = {"all"};
    std::string newton3 = "on";
    uint32_t seed = 42;
    bool verify = false;
    bool verifyOnly = false;
    int64_t verifyParticles = 16;
    double verifyTolerance = 1e-10;
    std::string verifyBaseline;
    std::string verifyCandidate;
};

enum FunctorKernel
{
    AOS,
    SOASINGLE,
    SOAPAIR,
    SOATRIPLE
};

double distSquared(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    using autopas::utils::ArrayMath::sub;
    using autopas::utils::ArrayMath::dot;
    const auto c = sub(a, b); // 3 FLOPS
    return dot(c, c); // 3+2=5 FLOPs
}

template <typename FunctorType>
void generateParticles(FunctorType& functor, std::vector<Cell>& cells, const size_t numberOfParticles,
                       const double cellSize, const FunctorKernel kernel, const uint32_t seed)
{
    // generate randomly distributed particles with deterministic seed
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(0.0, cellSize);

    auto fillCellWithParticles = [&](Cell& cell, const double xShift, const double yShift, const double zShift, const size_t startId)
    {
        for (size_t i = 0; i < numberOfParticles; ++i)
        {
            Particle p{
                {dis(gen) + xShift, dis(gen) + yShift, dis(gen) + zShift},
                {0., 0., 0.},
                startId + i,
                0
            };
            cell.addParticle(p); // Add particle to the current cell
        }
    };
    switch (kernel)
    {
    case AOS:
        fillCellWithParticles(cells[0], 0., 0., 0., 0);
        break;
    case SOATRIPLE:
        fillCellWithParticles(cells[2], 0., cellSize, 0., 2 * numberOfParticles);
        functor.SoALoader(cells[2], cells[2]._particleSoABuffer, 0, false);
        [[fallthrough]];
    case SOAPAIR:
        fillCellWithParticles(cells[1], cellSize, 0., 0., numberOfParticles);
        functor.SoALoader(cells[1], cells[1]._particleSoABuffer, 0, false);
        [[fallthrough]];
    case SOASINGLE:
        fillCellWithParticles(cells[0], 0., 0., 0., 0);
        functor.SoALoader(cells[0], cells[0]._particleSoABuffer, 0, false);
        break;
    }
}

template <typename FunctorType>
void applyAoSFunctor(FunctorType& functor, Cell& cell, bool newton3)
{
    for (std::size_t i = 0; i < cell.size(); ++i)
    {
        for (std::size_t j = i + 1; j < cell.size(); ++j)
        {
            for (std::size_t k = j + 1; k < cell.size(); ++k)
            {
                functor.AoSFunctor(cell[i], cell[j], cell[k], newton3);
            }
        }
    }
}

template <typename FunctorType>
void applyFunctorOnParticles(FunctorType& functor, std::vector<Cell>& cells, const FunctorKernel kernel, bool newton3)
{
    switch (kernel)
    {
    case AOS:
        applyAoSFunctor(functor, cells[0], newton3);
        break;
    case SOASINGLE:
        functor.SoAFunctorSingle(cells[0]._particleSoABuffer, newton3);
        break;
    case SOAPAIR:
        functor.SoAFunctorPair(cells[0]._particleSoABuffer, cells[1]._particleSoABuffer, newton3);
        break;
    case SOATRIPLE:
        functor.SoAFunctorTriple(cells[0]._particleSoABuffer, cells[1]._particleSoABuffer, cells[2]._particleSoABuffer,
                                 newton3);
        break;
    }
}

std::tuple<size_t, size_t> countInteractions(std::vector<Cell>& cells, const double cutoff, const FunctorKernel kernel)
{
    size_t calcsDist{0};
    size_t calcsForce{0};
    const auto cutoffSquared{cutoff * cutoff};

    switch (kernel)
    {
    case AOS:
    case SOASINGLE:
        for (size_t i = 0; i < cells[0].size(); i++)
        {
            for (size_t j = i + 1; j < cells[0].size(); j++)
            {
                if (distSquared(cells[0][i].getR(), cells[0][j].getR()) > cutoffSquared) continue;
                for (size_t k = j + 1; k < cells[0].size(); k++)
                {
                    if (distSquared(cells[0][i].getR(), cells[0][k].getR()) > cutoffSquared) continue;
                    if (distSquared(cells[0][j].getR(), cells[0][k].getR()) > cutoffSquared) continue;
                    ++calcsForce;
                }
            }
        }
        calcsDist = cells[0].size() * (cells[0].size() - 1) * (cells[0].size() - 2) / 6;
        break;
    case SOAPAIR:
        for (size_t i = 0; i < cells[0].size(); i++)
        {
            for (size_t j = i + 1; j < cells[0].size(); j++)
            {
                if (distSquared(cells[0][i].getR(), cells[0][j].getR()) > cutoffSquared) continue;
                for (size_t k = 0; k < cells[1].size(); k++)
                {
                    if (distSquared(cells[0][i].getR(), cells[1][k].getR()) > cutoffSquared) continue;
                    if (distSquared(cells[0][j].getR(), cells[1][k].getR()) > cutoffSquared) continue;
                    ++calcsForce;
                }
            }
            for (size_t j = 0; j < cells[1].size(); j++)
            {
                if (distSquared(cells[0][i].getR(), cells[1][j].getR()) > cutoffSquared) continue;
                for (size_t k = j + 1; k < cells[1].size(); k++)
                {
                    if (distSquared(cells[0][i].getR(), cells[1][k].getR()) > cutoffSquared) continue;
                    if (distSquared(cells[1][j].getR(), cells[1][k].getR()) > cutoffSquared) continue;
                    ++calcsForce;
                }
            }
        }
        calcsDist = cells[0].size() * cells[1].size() * (cells[0].size() + cells[1].size() - 2) / 2;
        break;
    case SOATRIPLE:
        for (size_t i = 0; i < cells[0].size(); i++)
        {
            for (size_t j = 0; j < cells[1].size(); j++)
            {
                if (distSquared(cells[0][i].getR(), cells[1][j].getR()) > cutoffSquared) continue;
                for (size_t k = 0; k < cells[2].size(); k++)
                {
                    if (distSquared(cells[0][i].getR(), cells[2][k].getR()) > cutoffSquared) continue;
                    if (distSquared(cells[1][j].getR(), cells[2][k].getR()) > cutoffSquared) continue;
                    ++calcsForce;
                }
            }
        }
        calcsDist = cells[0].size() * cells[1].size() * cells[2].size();
        break;
    default:
        break;
    }
    return {calcsDist, calcsForce};
}


template <typename FunctorType, typename Factory>
static void BM_Functor(benchmark::State& state, Factory factory, FunctorKernel kernel, bool newton3, uint32_t seed, size_t poolSize)
{
    const auto numParticles = static_cast<size_t>(state.range(0));
    const auto cellSize = static_cast<double>(state.range(1));
    const auto cutoff = static_cast<double>(state.range(2));

    auto functor = factory(cutoff);

    std::size_t calcsDistTotal = 0;
    std::size_t calcsForceTotal = 0;
    std::vector<std::vector<Cell>> cellPool(poolSize, std::vector<Cell>{3});

    for (size_t poolIdx = 0; poolIdx < cellPool.size(); ++poolIdx) {
        generateParticles(functor, cellPool[poolIdx], numParticles, cellSize, kernel, seed + static_cast<uint32_t>(poolIdx));
    }

    size_t pool_index = 0;

    for (auto _ : state)
    {
        auto& currentCells = cellPool[pool_index];

        {
            ITTScope itt;
            applyFunctorOnParticles(functor, currentCells, kernel, newton3);
        }

        pool_index = (pool_index + 1) % poolSize;
    }

    state.SetComplexityN(numParticles);

    const auto iters = static_cast<double>(state.iterations());
    const auto avg = std::min({5.0, iters, static_cast<double>(poolSize)});
    // Count interactions for first 5 or fewer cells
    for (size_t i = 0; i < static_cast<size_t>(avg); ++i)
    {
        const auto [calcsDist, calcsForce] = countInteractions(cellPool[i], cutoff, kernel);
        calcsDistTotal += calcsDist;
        calcsForceTotal += calcsForce;
    }

    // Per-iteration averages and hit rate as user counters.
    const double avgDist = static_cast<double>(calcsDistTotal) / avg;
    const double avgForce = static_cast<double>(calcsForceTotal) / avg;
    const double hitRate = (avgDist > 0.0) ? (avgForce / avgDist * 100.0) : 0.0;

    using benchmark::Counter;
    auto roundToPrecision = [](const double x, const unsigned int precision)
    {
        return std::round(x * std::pow(10, precision)) / std::pow(10, precision);
    };
    state.counters["Hit%"] = roundToPrecision(hitRate, 2);
    if (avgDist > 0.0) {
        state.counters["Time/Triplet"] = Counter(avgDist, Counter::kIsIterationInvariantRate | Counter::kInvert, Counter::OneK::kIs1000);
        state.counters["Triplets/s"] = Counter(avgDist, Counter::kIsIterationInvariantRate, Counter::OneK::kIs1000);
    }
    if (avgForce > 0.0) {
        state.counters["Time/Interaction"] = Counter(avgForce, Counter::kIsIterationInvariantRate | Counter::kInvert, Counter::OneK::kIs1000);
    }
}

template <typename T>
struct is_calculate_globals : std::false_type {};

template <template <typename, bool, autopas::FunctorN3Modes, bool, bool> class Template,
          typename P, bool M, autopas::FunctorN3Modes N, bool G, bool F>
struct is_calculate_globals<Template<P, M, N, G, F>> : std::integral_constant<bool, G> {};

struct GlobalsOutput {
    double potentialEnergy = 0.0;
    double virial = 0.0;
};

struct FunctorInfo {
    std::string name;
    std::string description;
    bool calculatesGlobals = false;
    std::function<void(const std::string& kernelName, FunctorKernel kernel, const BenchmarkConfig& config)> registerBenchmark;
    std::function<std::optional<GlobalsOutput>(std::vector<Cell>& cells, FunctorKernel kernel, bool newton3, double cutoff)> runOnce;
};

class FunctorRegistry {
public:
    template <typename FunctorType, typename Factory>
    void registerFunctor(const std::string& name, const std::string& description, Factory functorFactory) {
        FunctorInfo info;
        info.name = name;
        info.description = description;
        constexpr bool hasGlobals = is_calculate_globals<FunctorType>::value;
        info.calculatesGlobals = hasGlobals;

        info.registerBenchmark = [name, functorFactory](const std::string& kernelName, FunctorKernel kernel, const BenchmarkConfig& config) {
            auto registerVariant = [&](bool n3, const std::string& suffix) {
                auto* b = benchmark::RegisterBenchmark(
                    name + "_" + kernelName + suffix,
                    [=](benchmark::State& state) {
                        BM_Functor<FunctorType>(state, functorFactory, kernel, n3, config.seed, config.cellPoolSize);
                    });

                if (!config.particles.empty()) {
                    for (int64_t p : config.particles) {
                        b->Args({p, static_cast<int64_t>(config.cellSize), static_cast<int64_t>(config.cutoff)});
                    }
                } else {
                    b->RangeMultiplier(2)
                     ->Ranges({{config.minParticles, config.maxParticles},
                               {static_cast<int64_t>(config.cellSize), static_cast<int64_t>(config.cellSize)},
                               {static_cast<int64_t>(config.cutoff), static_cast<int64_t>(config.cutoff)}});
                }
            };

            auto lowerN3 = config.newton3;
            for (char& c : lowerN3) c = std::tolower(static_cast<unsigned char>(c));

            if (lowerN3 == "both") {
                registerVariant(true, "_N3ON");
                registerVariant(false, "_N3OFF");
            } else if (lowerN3 == "off") {
                registerVariant(false, "_N3OFF");
            } else {
                registerVariant(true, "");
            }
        };

        info.runOnce = [functorFactory](std::vector<Cell>& cells, FunctorKernel kernel, bool newton3, double cutoff)
            -> std::optional<GlobalsOutput> {
            auto functor = functorFactory(cutoff);
            functor.initTraversal();
            if (kernel != AOS) {
                for (auto& cell : cells) {
                    if (!cell.isEmpty()) {
                        functor.SoALoader(cell, cell._particleSoABuffer, 0, false);
                    }
                }
            }
            applyFunctorOnParticles(functor, cells, kernel, newton3);
            if (kernel != AOS) {
                for (auto& cell : cells) {
                    if (!cell.isEmpty()) {
                        functor.SoAExtractor(cell, cell._particleSoABuffer, 0);
                    }
                }
            }
            functor.endTraversal(newton3);

            if constexpr (hasGlobals) {
                return GlobalsOutput{functor.getPotentialEnergy(), functor.getVirial()};
            } else {
                return std::nullopt;
            }
        };

        _functors[name] = std::move(info);
        _names.push_back(name);
    }

    const std::vector<std::string>& getNames() const { return _names; }

    bool has(const std::string& name) const {
        return _functors.contains(name);
    }

    const FunctorInfo& get(const std::string& name) const {
        return _functors.at(name);
    }

private:
    std::unordered_map<std::string, FunctorInfo> _functors;
    std::vector<std::string> _names;
};

void initRegistry(FunctorRegistry& reg) {
    reg.registerFunctor<ATM>("ATM", "Axilrod-Teller-Muto Baseline Functor", [](const double cutoff) {
        ATM f{cutoff};
        f.setParticleProperties(1.0);
        return f;
    });

    reg.registerFunctor<ATM2>("ATM2", "Axilrod-Teller-Muto Candidate Functor", [](const double cutoff) {
        ATM2 f{cutoff};
        f.setParticleProperties(1.0);
        return f;
    });

    reg.registerFunctor<ATM3>("ATM3", "Axilrod-Teller-Muto Hybrid Functor", [](const double cutoff) {
        ATM3 f{cutoff};
        f.setParticleProperties(1.0);
        return f;
    });

    reg.registerFunctor<ATMGlobals>("ATMGlobals", "ATM Functor with Globals calculation", [](const double cutoff) {
        ATMGlobals f{cutoff};
        f.setParticleProperties(1.0);
        return f;
    });

    reg.registerFunctor<ATM2Globals>("ATM2Globals", "ATM2 Functor with Globals calculation", [](const double cutoff) {
        ATM2Globals f{cutoff};
        f.setParticleProperties(1.0);
        return f;
    });

    reg.registerFunctor<ATM3Globals>("ATM3Globals", "ATM3 Functor with Globals calculation", [](const double cutoff) {
        ATM3Globals f{cutoff};
        f.setParticleProperties(1.0);
        return f;
    });
}

void registerFunctors(const BenchmarkConfig& config, const FunctorRegistry& registry)
{
    std::cout << "==========================================" << std::endl;
    std::cout << "AutoPas Functor Benchmark" << std::endl;
    std::cout << "AutoPas Branch: " << AUTOPAS_BRANCH << std::endl;
    std::cout << "AutoPas Commit: " << AUTOPAS_COMMIT << std::endl;
    std::cout << "==========================================" << std::endl;

    constexpr std::array kernels = {
        std::make_pair("AoS", AOS),
        std::make_pair("SoASingle", SOASINGLE),
        std::make_pair("SoAPair", SOAPAIR),
        std::make_pair("SoATriple", SOATRIPLE)
    };

    auto stringsAreEqual = [](const std::string& a, const std::string& b) {
        return std::ranges::equal(a, b,
                                  [](const char ca, const char cb) { return std::tolower(static_cast<unsigned char>(ca)) == std::tolower(static_cast<unsigned char>(cb)); });
    };

    // Resolve the requested functors
    std::vector<std::string> requestedFunctors;
    bool allFunctorsRequested = false;
    for (const auto& functor : config.targetFunctors) {
        if (stringsAreEqual(functor, "all")) {
            allFunctorsRequested = true;
            break;
        }
    }

    if (allFunctorsRequested) {
        for (const auto& functorName : registry.getNames()) {
            requestedFunctors.push_back(functorName);
        }
    } else {
        for (const auto& requestedFunctor : config.targetFunctors) {
            for (const auto& functorName : registry.getNames()) {
                if (stringsAreEqual(requestedFunctor, functorName)) {
                    if (std::ranges::find(requestedFunctors, functorName) == requestedFunctors.end()) {
                        requestedFunctors.push_back(functorName);
                    }
                }
            }
        }
    }

    // Resolve requested functor kernels
    std::vector<std::pair<std::string, FunctorKernel>> requestedKernels;
    bool allKernelsRequested = false;
    for (const auto& functorKernel : config.targetKernels) {
        if (stringsAreEqual(functorKernel, "all")) {
            allKernelsRequested = true;
            break;
        }
    }

    if (allKernelsRequested) {
        for (const auto& kernel : kernels) {
            requestedKernels.emplace_back(kernel);
        }
    } else {
        for (const auto& requestedKernel : config.targetKernels) {
            for (const auto& [kernelName, functorKernel] : kernels) {
                if (stringsAreEqual(requestedKernel, kernelName)) {
                    auto it = std::ranges::find_if(requestedKernels,
                                                   [&](const auto& pair) { return pair.second == functorKernel; });
                    if (it == requestedKernels.end()) {
                        requestedKernels.emplace_back(kernelName, functorKernel);
                    }
                }
            }
        }
    }

    // Register each (functor, kernel) pair cleanly
    for (const auto& functorName : requestedFunctors) {
        if (!registry.has(functorName)) continue;
        const auto& info = registry.get(functorName);
        for (const auto& [kernelName, kernel] : requestedKernels) {
            info.registerBenchmark(kernelName, kernel, config);
        }
    }
}

struct VerifyResult {
    bool passed = true;
    double maxAbsDiff = 0.0;
    double maxRelDiff = 0.0;
    double maxBaselineForce = 0.0;
    std::string firstMismatch;

    // Globals verification (energy & virial)
    bool verifiedGlobals = false;
    double baselineUpot = 0.0;
    double candidateUpot = 0.0;
    double diffUpot = 0.0;
    double relDiffUpot = 0.0;

    double baselineVirial = 0.0;
    double candidateVirial = 0.0;
    double diffVirial = 0.0;
    double relDiffVirial = 0.0;
};

VerifyResult verifyOneKernel(const FunctorInfo& baselineInfo, const FunctorInfo& candidateInfo,
                             const FunctorKernel kernel, const bool newton3, const size_t numParticles,
                             const double cellSize, const double cutoff, const uint32_t seed, const double tolerance)
{
    std::vector<Cell> cellsBaseline(3);
    std::vector<Cell> cellsCandidate(3);

    std::mt19937 genBaseline(seed);
    std::mt19937 genCandidate(seed);
    std::uniform_real_distribution<double> dis(0.0, cellSize);

    auto fillCells = [&](std::vector<Cell>& cells, std::mt19937& gen) {
        auto fillCell = [&](Cell& cell, const double xShift, const double yShift, const double zShift, const size_t startId) {
            for (size_t i = 0; i < numParticles; ++i) {
                Particle p{
                    {dis(gen) + xShift, dis(gen) + yShift, dis(gen) + zShift},
                    {0., 0., 0.},
                    startId + i,
                    0
                };
                cell.addParticle(p);
            }
        };

        switch (kernel) {
        case AOS:
        case SOASINGLE:
            fillCell(cells[0], 0., 0., 0., 0);
            break;
        case SOAPAIR:
            fillCell(cells[0], 0., 0., 0., 0);
            fillCell(cells[1], cellSize, 0., 0., numParticles);
            break;
        case SOATRIPLE:
            fillCell(cells[0], 0., 0., 0., 0);
            fillCell(cells[1], cellSize, 0., 0., numParticles);
            fillCell(cells[2], 0., cellSize, 0., 2 * numParticles);
            break;
        }
    };

    fillCells(cellsBaseline, genBaseline);
    fillCells(cellsCandidate, genCandidate);

    const auto baselineGlobals = baselineInfo.runOnce(cellsBaseline, kernel, newton3, cutoff);
    const auto candidateGlobals = candidateInfo.runOnce(cellsCandidate, kernel, newton3, cutoff);

    VerifyResult res;
    for (size_t c = 0; c < cellsBaseline.size(); ++c) {
        for (size_t p = 0; p < cellsBaseline[c].size(); ++p) {
            const auto fBaseline = cellsBaseline[c][p].getF();
            const auto fCandidate = cellsCandidate[c][p].getF();

            for (int d = 0; d < 3; ++d) {
                const double baselineVal = fBaseline[d];
                const double candidateVal = fCandidate[d];
                res.maxBaselineForce = std::max(res.maxBaselineForce, std::abs(baselineVal));

                double diff = std::abs(baselineVal - candidateVal);
                if (std::isnan(diff) || std::isinf(diff)) {
                    res.passed = false;
                    res.firstMismatch = "NaN or Inf in force values!";
                    return res;
                }

                res.maxAbsDiff = std::max(res.maxAbsDiff, diff);
                const double denom = std::max(std::abs(baselineVal), std::abs(candidateVal));
                if (denom > 1e-12) {
                    res.maxRelDiff = std::max(res.maxRelDiff, diff / denom);
                }

                if (diff > tolerance && res.firstMismatch.empty()) {
                    std::ostringstream ss;
                    ss << "Force mismatch at Cell " << c << ", Particle " << p
                       << ", Dim " << (d == 0 ? "X" : (d == 1 ? "Y" : "Z"))
                       << ": Baseline=" << baselineVal << ", Candidate=" << candidateVal
                       << " (Diff=" << diff << ")";
                    res.firstMismatch = ss.str();
                }
            }
        }
    }

    if (res.maxAbsDiff > tolerance) {
        res.passed = false;
    }

    // Verify globals if both baseline and candidate calculate globals
    if (baselineGlobals.has_value() && candidateGlobals.has_value()) {
        res.verifiedGlobals = true;
        res.baselineUpot = baselineGlobals->potentialEnergy;
        res.candidateUpot = candidateGlobals->potentialEnergy;
        res.diffUpot = std::abs(res.baselineUpot - res.candidateUpot);
        const double denomUpot = std::max(std::abs(res.baselineUpot), std::abs(res.candidateUpot));
        if (denomUpot > 1e-12) {
            res.relDiffUpot = res.diffUpot / denomUpot;
        }

        res.baselineVirial = baselineGlobals->virial;
        res.candidateVirial = candidateGlobals->virial;
        res.diffVirial = std::abs(res.baselineVirial - res.candidateVirial);
        const double denomVirial = std::max(std::abs(res.baselineVirial), std::abs(res.candidateVirial));
        if (denomVirial > 1e-12) {
            res.relDiffVirial = res.diffVirial / denomVirial;
        }

        if (std::isnan(res.diffUpot) || std::isinf(res.diffUpot)) {
            res.passed = false;
            if (res.firstMismatch.empty()) res.firstMismatch = "NaN or Inf in Potential Energy!";
        } else if (res.diffUpot > tolerance) {
            res.passed = false;
            if (res.firstMismatch.empty()) {
                std::ostringstream ss;
                ss << "Upot mismatch: Baseline=" << res.baselineUpot
                   << ", Candidate=" << res.candidateUpot << " (Diff=" << res.diffUpot << ")";
                res.firstMismatch = ss.str();
            }
        }

        if (std::isnan(res.diffVirial) || std::isinf(res.diffVirial)) {
            res.passed = false;
            if (res.firstMismatch.empty()) res.firstMismatch = "NaN or Inf in Virial!";
        } else if (res.diffVirial > tolerance) {
            res.passed = false;
            if (res.firstMismatch.empty()) {
                std::ostringstream ss;
                ss << "Virial mismatch: Baseline=" << res.baselineVirial
                   << ", Candidate=" << res.candidateVirial << " (Diff=" << res.diffVirial << ")";
                res.firstMismatch = ss.str();
            }
        }
    }

    return res;
}

bool runVerification(const BenchmarkConfig& config, const FunctorRegistry& registry) {
    std::string baselineName = config.verifyBaseline;
    std::string candidateName = config.verifyCandidate;

    auto stringsAreEqual = [](const std::string& a, const std::string& b) {
        return std::ranges::equal(a, b,
                                  [](const char ca, const char cb) {
                                      return std::tolower(static_cast<unsigned char>(ca)) == std::tolower(static_cast<unsigned char>(cb));
                                  });
    };

    std::vector<std::pair<std::string, std::string>> verifyPairs;

    if (!baselineName.empty() || !candidateName.empty()) {
        if (baselineName.empty()) {
            baselineName = registry.getNames()[0];
        }
        if (candidateName.empty()) {
            candidateName = (registry.getNames().size() >= 2) ? registry.getNames()[1] : baselineName;
        }
        verifyPairs.emplace_back(baselineName, candidateName);
    } else {
        // Collect requested targets
        std::vector<std::string> validTargets;
        for (const auto& targetName : config.targetFunctors) {
            for (const auto& registeredName : registry.getNames()) {
                if (stringsAreEqual(targetName, registeredName)) {
                    if (std::ranges::find(validTargets, registeredName) == validTargets.end()) {
                        validTargets.push_back(registeredName);
                    }
                }
            }
        }
        if (validTargets.empty()) {
            validTargets = registry.getNames();
        }

        // If target functors include both standard and globals variants, verify both pairs
        bool hasATM = std::ranges::find(validTargets, "ATM") != validTargets.end();
        bool hasATM2 = std::ranges::find(validTargets, "ATM2") != validTargets.end();
        bool hasATM3 = std::ranges::find(validTargets, "ATM3") != validTargets.end();
        bool hasATMGlobals = std::ranges::find(validTargets, "ATMGlobals") != validTargets.end();
        bool hasATM2Globals = std::ranges::find(validTargets, "ATM2Globals") != validTargets.end();
        bool hasATM3Globals = std::ranges::find(validTargets, "ATM3Globals") != validTargets.end();

        if (hasATM && hasATM2) {
            verifyPairs.emplace_back("ATM", "ATM2");
        }
        if (hasATM && hasATM3) {
            verifyPairs.emplace_back("ATM", "ATM3");
        }
        if (hasATMGlobals && hasATM2Globals) {
            verifyPairs.emplace_back("ATMGlobals", "ATM2Globals");
        }
        if (hasATMGlobals && hasATM3Globals) {
            verifyPairs.emplace_back("ATMGlobals", "ATM3Globals");
        }

        if (verifyPairs.empty()) {
            if (validTargets.size() >= 2) {
                verifyPairs.emplace_back(validTargets[0], validTargets[1]);
            } else if (registry.getNames().size() >= 2) {
                verifyPairs.emplace_back(registry.getNames()[0], registry.getNames()[1]);
            } else {
                verifyPairs.emplace_back(registry.getNames()[0], registry.getNames()[0]);
            }
        }
    }

    const std::vector<std::pair<std::string, FunctorKernel>> allPossibleKernels = {
        {"AoS", AOS},
        {"SoASingle", SOASINGLE},
        {"SoAPair", SOAPAIR},
        {"SoATriple", SOATRIPLE}
    };

    std::vector<std::pair<std::string, FunctorKernel>> testKernels;
    bool allKernelsReq = false;
    for (const auto& m : config.targetKernels) {
        if (stringsAreEqual(m, "all")) {
            allKernelsReq = true;
            break;
        }
    }
    if (allKernelsReq) {
        testKernels = allPossibleKernels;
    } else {
        for (const auto& req : config.targetKernels) {
            for (const auto& [kernelName, kernel] : allPossibleKernels) {
                if (stringsAreEqual(req, kernelName)) {
                    auto it = std::ranges::find_if(testKernels,
                                                   [&](const auto& p) { return p.second == kernel; });
                    if (it == testKernels.end()) {
                        testKernels.emplace_back(kernelName, kernel);
                    }
                }
            }
        }
    }

    std::vector<bool> n3Options;
    auto lowerN3 = config.newton3;
    for (char& c : lowerN3) c = std::tolower(static_cast<unsigned char>(c));

    if (lowerN3 == "off") {
        n3Options = {false};
    } else if (lowerN3 == "both" || config.verifyOnly) {
        n3Options = {true, false};
    } else {
        n3Options = {true};
    }

    bool allPassed = true;

    for (const auto& [baseName, candName] : verifyPairs) {
        if (!registry.has(baseName) || !registry.has(candName)) {
            std::cerr << "[VERIFY ERROR] Cannot run verification: Need both '" << baseName 
                      << "' and '" << candName << "' registered." << std::endl;
            return false;
        }

        const auto& baselineInfo = registry.get(baseName);
        const auto& candidateInfo = registry.get(candName);

        std::cout << "\n==========================================" << std::endl;
        std::cout << "Functor Correctness Verification" << std::endl;
        std::cout << "Baseline:  " << baseName << (baselineInfo.calculatesGlobals ? " (Globals: YES)" : " (Globals: NO)") << std::endl;
        std::cout << "Candidate: " << candName << (candidateInfo.calculatesGlobals ? " (Globals: YES)" : " (Globals: NO)") << std::endl;
        if (baselineInfo.calculatesGlobals && candidateInfo.calculatesGlobals) {
            std::cout << "Scope:     Forces, Potential Energy, Virial" << std::endl;
        } else if (baselineInfo.calculatesGlobals != candidateInfo.calculatesGlobals) {
            std::cout << "Scope:     Forces only (one functor does not calculate globals)" << std::endl;
        } else {
            std::cout << "Scope:     Forces only" << std::endl;
        }
        std::cout << "Particles/cell: " << config.verifyParticles 
                  << " | Cell size: " << config.cellSize 
                  << " | Cutoff: " << config.cutoff 
                  << " | Tol: " << std::scientific << std::setprecision(2) << config.verifyTolerance 
                  << std::defaultfloat << std::endl;
        std::cout << "==========================================" << std::endl;

        for (const bool n3 : n3Options) {
            for (const auto& [kernelName, kernel] : testKernels) {
                auto res = verifyOneKernel(baselineInfo, candidateInfo, kernel, n3,
                                           config.verifyParticles, config.cellSize,
                                           config.cutoff, config.seed, config.verifyTolerance);

                std::cout << "[VERIFY] " << std::left << std::setw(12) << kernelName
                          << " | N3: " << (n3 ? "ON " : "OFF") << " | ";

                if (res.passed) {
                    std::cout << "PASS (force diff: " << std::scientific << std::setprecision(2) << res.maxAbsDiff;
                    if (res.verifiedGlobals) {
                        std::cout << ", Upot diff: " << res.diffUpot
                                  << ", Virial diff: " << res.diffVirial;
                    }
                    std::cout << ")";
                } else {
                    allPassed = false;
                    std::cout << "FAIL! (force diff: " << std::scientific << std::setprecision(2) << res.maxAbsDiff
                              << ", tol: " << config.verifyTolerance;
                    if (res.verifiedGlobals) {
                        std::cout << ", Upot diff: " << res.diffUpot
                                  << ", Virial diff: " << res.diffVirial;
                    }
                    std::cout << ")";
                    if (!res.firstMismatch.empty()) {
                        std::cout << "\n         -> " << res.firstMismatch;
                    }
                }
                std::cout << std::defaultfloat << std::endl;
            }
        }
    }

    std::cout << "==========================================" << std::endl;
    if (allPassed) {
        std::cout << "Verification Result: ALL CHECKS PASSED" << std::endl;
    } else {
        std::cout << "Verification Result: FAILED" << std::endl;
    }
    std::cout << "==========================================\n" << std::endl;

    return allPassed;
}

void setupCLI(CLI::App& app, BenchmarkConfig& config, const FunctorRegistry& registry) {
    app.add_option("--min", config.minParticles, "Minimum number of particles")->default_val(1);
    app.add_option("--max", config.maxParticles, "Maximum number of particles")->default_val(512);
    app.add_option("-p,--particles", config.particles, "Comma-separated list of particle counts (overrides --min/--max)")
       ->delimiter(',')
       ->check(CLI::PositiveNumber);
    app.add_option("-c,--cell-size", config.cellSize, "Size of the simulation cell")->default_val(3);
    app.add_option("-r,--cutoff", config.cutoff, "Cutoff radius for interactions")->default_val(3);
    app.add_option("-s,--seed", config.seed, "Random seed for reproducible particle generation")->default_val(42);
    app.add_option("--pool-size", config.cellPoolSize, "Number of cell stencils pre-generated in pool")
       ->default_val(1000)
       ->check(CLI::PositiveNumber);

    app.add_flag("-v,--verify", config.verify, "Verify correctness between baseline and candidate functors before benchmarking");
    app.add_flag("--verify-only", config.verifyOnly, "Run correctness verification and exit immediately");
    app.add_option("--verify-baseline", config.verifyBaseline, "Baseline functor for verification (default: first from -f)");
    app.add_option("--verify-candidate", config.verifyCandidate, "Candidate functor for verification (default: second from -f)");
    app.add_option("--verify-particles", config.verifyParticles, "Number of particles per cell for verification")->default_val(16);
    app.add_option("--verify-tol", config.verifyTolerance, "Numerical tolerance for verification")->default_str("1e-10");

    auto validFunctors = registry.getNames();
    validFunctors.emplace_back("all");

    app.add_option("-f,--functor", config.targetFunctors, "Comma-separated list of functors to test")
           ->check(CLI::IsMember(validFunctors, CLI::ignore_case))
           ->delimiter(',');

    app.add_option("-k,--kernel", config.targetKernels, "Comma-separated list of kernels to test")
       ->check(CLI::IsMember({"AoS", "SoASingle", "SoAPair", "SoATriple", "all"}, CLI::ignore_case))
       ->delimiter(',');

    app.add_option("--n3", config.newton3, "Newton3 option: {on, off, both}")
       ->default_str("on")
       ->check(CLI::IsMember({"on", "off", "both"}, CLI::ignore_case));
    app.add_flag_callback("--no-n3", [&config]() { config.newton3 = "off"; }, "Disable Newton3 (alias for --n3 off)");
}

bool handleHelpFlag(const int argc, char** argv, const CLI::App& app) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            std::cout << app.help() << "\n";
            std::cout << "--- Google Benchmark Options ---\n";
            return true;
        }
    }
    return false;
}

int main(int argc, char** argv) {
#ifdef ENABLE_ITT
    // Start paused so setup, CLI parsing, container building, and verification are ignored by profilers
    __itt_pause();
#endif
    // Add the version info to the JSON metadata
    benchmark::AddCustomContext("AutoPas Branch", AUTOPAS_BRANCH);
    benchmark::AddCustomContext("AutoPas Commit", AUTOPAS_COMMIT);
    benchmark::MaybeReenterWithoutASLR(argc, argv);

    // Initialize Functor Registry
    FunctorRegistry registry;
    initRegistry(registry);

    // Read CLI arguments
    CLI::App app{"AutoPas 3-Body Functor Benchmark"};
    BenchmarkConfig config;
    setupCLI(app, config, registry);

    if (handleHelpFlag(argc, argv, app)) {
        benchmark::Initialize(&argc, argv);
        return 0;
    }

    benchmark::Initialize(&argc, argv);
    CLI11_PARSE(app, argc, argv);

    if (!config.particles.empty()) {
        std::ranges::sort(config.particles);
        auto [first, last] = std::ranges::unique(config.particles);
        config.particles.erase(first, last);
    }

    registerFunctors(config, registry);

    if (config.verify || config.verifyOnly) {
        const bool correctResults = runVerification(config, registry);
        if (not correctResults) {
            std::cerr << "[ERROR] Verification failed! Aborting." << std::endl;
            return 1;
        }
        if (config.verifyOnly) {
            return 0;
        }
    }

    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
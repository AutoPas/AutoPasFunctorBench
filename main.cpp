#include <iostream>

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h>
#include <AxilrodTellerMutoFunctor2.h>
#include "benchmark/benchmark.h"

#include <autopas/cells/FullParticleCell.h>
#include <random>
#include <CLI/CLI.hpp>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include <cctype>

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

struct BenchmarkConfig {
    int64_t minParticles = 1;
    int64_t maxParticles = 512;
    int64_t cellSize = 3;
    int64_t cutoff = 3;
    std::vector<std::string> targetFunctors = {"all"};
    std::vector<std::string> targetModes = {"all"};
    bool newton3 = true;
    uint32_t seed = 42;
};

enum FunctorMode
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
                       const double cellSize, const FunctorMode mode, const uint32_t seed)
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
    switch (mode)
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
void applyFunctorOnParticles(FunctorType& functor, std::vector<Cell>& cells, const FunctorMode mode, bool newton3)
{
    switch (mode)
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

std::tuple<size_t, size_t> countInteractions(std::vector<Cell>& cells, const double cutoff, const FunctorMode mode)
{
    size_t calcsDist{0};
    size_t calcsForce{0};
    const auto cutoffSquared{cutoff * cutoff};

    switch (mode)
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
static void BM_Functor(benchmark::State& state, Factory factory, FunctorMode functorMode, bool newton3, uint32_t seed)
{
    const auto numParticles = static_cast<size_t>(state.range(0));
    const auto cellSize = static_cast<double>(state.range(1));
    const auto cutoff = static_cast<double>(state.range(2));

    auto functor = factory(cutoff);

    std::size_t calcsDistTotal = 0;
    std::size_t calcsForceTotal = 0;
    constexpr size_t poolSize = 1000;
    std::vector<std::vector<Cell>> cellPool(poolSize, std::vector<Cell>{3});

    for (size_t poolIdx = 0; poolIdx < cellPool.size(); ++poolIdx) {
        generateParticles(functor, cellPool[poolIdx], numParticles, cellSize, functorMode, seed + static_cast<uint32_t>(poolIdx));
    }

    size_t pool_index = 0;

    for (auto _ : state)
    {
        auto& currentCells = cellPool[pool_index];

        applyFunctorOnParticles(functor, currentCells, functorMode, newton3);

        pool_index = (pool_index + 1) % poolSize;
    }

    state.SetComplexityN(numParticles);

    const auto iters = static_cast<double>(state.iterations());
    const auto avg = std::min(5.0, iters);
    // Count interactions for first 5 or fewer cells
    for (auto i = 0; i < avg; i++)
    {
        const auto [calcsDist, calcsForce] = countInteractions(cellPool[i], cutoff, functorMode);
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
    state.counters["HitRate [%]"] = roundToPrecision(hitRate, 2);
    state.counters["# of Triplets"] = avgDist;
    if (avgDist > 0.0) {
        state.counters["Time per Triplet"] = Counter(avgDist, Counter::kIsIterationInvariantRate | Counter::kInvert, Counter::OneK::kIs1000);
        state.counters["Triplets/s"] = Counter(avgDist, Counter::kIsIterationInvariantRate, Counter::OneK::kIs1000);
    }
    if (avgForce > 0.0) {
        state.counters["Time per Interaction"] = Counter(avgForce, Counter::kIsIterationInvariantRate | Counter::kInvert, Counter::OneK::kIs1000);
        state.counters["Interactions/s"] = Counter(avgForce, Counter::kIsIterationInvariantRate, Counter::OneK::kIs1000);
    }
    state.SetItemsProcessed(static_cast<int64_t>(avgDist * iters));
}

struct FunctorInfo {
    std::string name;
    std::string description;
    std::function<void(const std::string& modeName, FunctorMode mode, const BenchmarkConfig& config)> registerBenchmark;
    std::function<void(std::vector<Cell>& cells, FunctorMode mode, bool newton3, double cutoff)> runOnce;
};

class FunctorRegistry {
public:
    template <typename FunctorType, typename Factory>
    void registerFunctor(const std::string& name, const std::string& description, Factory functorFactory) {
        FunctorInfo info;
        info.name = name;
        info.description = description;

        info.registerBenchmark = [name, functorFactory](const std::string& modeName, FunctorMode mode, const BenchmarkConfig& config) {
            benchmark::RegisterBenchmark(
                "BM_" + name + "_" + modeName,
                [=](benchmark::State& state) {
                    BM_Functor<FunctorType>(state, functorFactory, mode, config.newton3, config.seed);
                })
                ->RangeMultiplier(2)
                ->Ranges({{config.minParticles, config.maxParticles},
                          {config.cellSize, config.cellSize},
                          {config.cutoff, config.cutoff}});
        };

        info.runOnce = [functorFactory](std::vector<Cell>& cells, FunctorMode mode, bool newton3, double cutoff) {
            auto functor = functorFactory(cutoff);
            for (auto& cell : cells) {
                functor.SoALoader(cell, cell._particleSoABuffer, 0, false);
            }
            applyFunctorOnParticles(functor, cells, mode, newton3);
            for (auto& cell : cells) {
                functor.SoAExtractor(cell, cell._particleSoABuffer, 0);
            }
        };

        _functors[name] = std::move(info);
        _names.push_back(name);
    }

    const std::vector<std::string>& getNames() const { return _names; }

    bool has(const std::string& name) const {
        return _functors.find(name) != _functors.end();
    }

    const FunctorInfo& get(const std::string& name) const {
        return _functors.at(name);
    }

private:
    std::unordered_map<std::string, FunctorInfo> _functors;
    std::vector<std::string> _names;
};

void initRegistry(FunctorRegistry& reg) {
    reg.registerFunctor<ATM>("ATM", "Axilrod-Teller-Muto Reference Functor", [](const double cutoff) {
        ATM f{cutoff};
        f.setParticleProperties(1.0);
        return f;
    });

    reg.registerFunctor<ATM2>("ATM2", "Axilrod-Teller-Muto Variation", [](const double cutoff) {
        ATM2 f{cutoff};
        f.setParticleProperties(1.0);
        return f;
    });

    reg.registerFunctor<ATMGlobals>("ATMGlobals", "ATM Functor with Globals calculation", [](const double cutoff) {
        ATMGlobals f{cutoff};
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

    constexpr std::array modes = {
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

    // Resolve requested functor modes
    std::vector<std::pair<std::string, FunctorMode>> requestedModes;
    bool allModesRequested = false;
    for (const auto& functorMode : config.targetModes) {
        if (stringsAreEqual(functorMode, "all")) {
            allModesRequested = true;
            break;
        }
    }

    if (allModesRequested) {
        for (const auto& mode : modes) {
            requestedModes.emplace_back(mode);
        }
    } else {
        for (const auto& requestedMode : config.targetModes) {
            for (const auto& [modeName, functorMode] : modes) {
                if (stringsAreEqual(requestedMode, modeName)) {
                    auto it = std::ranges::find_if(requestedModes,
                                                   [&](const auto& pair) { return pair.second == functorMode; });
                    if (it == requestedModes.end()) {
                        requestedModes.emplace_back(modeName, functorMode);
                    }
                }
            }
        }
    }

    // Register each (functor, mode) pair cleanly
    for (const auto& functorName : requestedFunctors) {
        if (!registry.has(functorName)) continue;
        const auto& info = registry.get(functorName);
        for (const auto& [modeName, mode] : requestedModes) {
            info.registerBenchmark(modeName, mode, config);
        }
    }
}

void setupCLI(CLI::App& app, BenchmarkConfig& config, const FunctorRegistry& registry) {
    app.add_option("--min", config.minParticles, "Minimum number of particles")->default_val(1);
    app.add_option("--max", config.maxParticles, "Maximum number of particles")->default_val(512);
    app.add_option("-c,--cell-size", config.cellSize, "Size of the simulation cell")->default_val(3);
    app.add_option("-r,--cutoff", config.cutoff, "Cutoff radius for interactions")->default_val(3);
    app.add_option("-s,--seed", config.seed, "Random seed for reproducible particle generation")->default_val(42);

    auto validFunctors = registry.getNames();
    validFunctors.emplace_back("all");

    app.add_option("-f,--functor", config.targetFunctors, "Comma-separated list of functors to test")
           ->check(CLI::IsMember(validFunctors, CLI::ignore_case))
           ->delimiter(',');

    app.add_option("-m,--mode", config.targetModes, "Comma-separated list of modes to test")
       ->check(CLI::IsMember({"AoS", "SoASingle", "SoAPair", "SoATriple", "all"}, CLI::ignore_case))
       ->delimiter(',');

    app.add_flag("--n3,!--no-n3", config.newton3, "Enable/Disable Newton3 (enabled by default)");
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
    // Add the version info to the JSON metadata
    benchmark::AddCustomContext("autopas_branch", AUTOPAS_BRANCH);
    benchmark::AddCustomContext("autopas_commit", AUTOPAS_COMMIT);
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

    registerFunctors(config, registry);

    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
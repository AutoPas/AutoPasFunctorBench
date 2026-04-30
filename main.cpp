#include <iostream>

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h>
#include "benchmark/benchmark.h"

#include <autopas/cells/FullParticleCell.h>
#include <random>
#include <CLI/CLI.hpp>

// type aliases for ease of use
using Particle = mdLib::MoleculeLJ;
using Cell = autopas::FullParticleCell<Particle>;

struct BenchmarkConfig {
    int64_t minParticles = 1;
    int64_t maxParticles = 512;
    int64_t cellSize = 3;
    int64_t cutoff = 3;
    std::string targetFunctor = "all";
    std::string targetMode = "all";
    bool newton3 = true;
};

enum FunctorMode
{
    AOS,
    SOASINGLE,
    SOAPAIR,
    SOATRIPLE
};

double distSquared(std::array<double, 3> a, std::array<double, 3> b)
{
    using autopas::utils::ArrayMath::sub;
    using autopas::utils::ArrayMath::dot;
    const auto c = sub(a, b); // 3 FLOPS
    return dot(c, c); // 3+2=5 FLOPs
}

template <typename FunctorType>
void generateParticles(FunctorType& functor, std::vector<Cell>& cells, const size_t numberOfParticles,
                       const double cellSize, const FunctorMode mode)
{
    // generate randomly distributed particles
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, cellSize);

    auto fillCellWithParticles = [&](Cell& cell, double xShift, double yShift, double zShift)
    {
        for (size_t particleId = 0; particleId < numberOfParticles; ++particleId)
        {
            Particle p{
                {dis(gen) + xShift, dis(gen) + yShift, dis(gen)},
                {0., 0., 0.,},
                particleId,
                0
            };
            cell.addParticle(p); // Add particle to the current cell
        }
    };
    switch (mode)
    {
    case AOS:
        fillCellWithParticles(cells[0], 0., 0., 0.);
        break;
    case SOATRIPLE:
        fillCellWithParticles(cells[2], 0., cellSize, 0.);
        functor.SoALoader(cells[2], cells[2]._particleSoABuffer, 0, false);
    case SOAPAIR:
        fillCellWithParticles(cells[1], cellSize, 0., 0.);
        functor.SoALoader(cells[1], cells[1]._particleSoABuffer, 0, false);
    case SOASINGLE:
        fillCellWithParticles(cells[0], 0., 0., 0.);
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


template <typename FunctorType>
static void BM_Functor(benchmark::State& state, FunctorMode functorMode, bool newton3)
{
    const auto numParticles = static_cast<size_t>(state.range(0));
    const auto cellSize = static_cast<double>(state.range(1));
    const auto cutoff = static_cast<double>(state.range(2));
    const double nu{1.0};

    FunctorType functor{cutoff};
    functor.setParticleProperties(nu);

    std::size_t calcsDistTotal = 0;
    std::size_t calcsForceTotal = 0;
    const size_t poolSize = 1000;
    std::vector<std::vector<Cell>> cellPool(poolSize, std::vector<Cell>{3});

    for (auto& cells : cellPool) {
        generateParticles(functor, cells, numParticles, cellSize, functorMode);
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
    // Count interactions for first 5 or less cells
    for (auto i = 0; i < avg; i++)
    {
        const auto [calcsDist, calcsForce] = countInteractions(cellPool[i], cutoff, functorMode);
        calcsDistTotal += calcsDist;
        calcsForceTotal += calcsForce;
    }

    // Per-iteration averages and hit rate as user counters.
    const double avgDist = static_cast<double>(calcsDistTotal) / avg;
    const double avgForce = static_cast<double>(calcsForceTotal) / avg;
    const double hitRate = avgForce / avgDist * 100.0;

    using benchmark::Counter;
    auto roundToPrecision = [](double x, unsigned int precision)
    {
        return std::round(x * std::pow(10, precision)) / std::pow(10, precision);
    };
    state.counters["HitRate [%]"] = roundToPrecision(hitRate, 2);
    state.counters["Time per Triplet"] = Counter(avgDist, Counter::kIsRate | Counter::kInvert, Counter::OneK::kIs1000);
    state.counters["Time per Interaction"] = Counter(avgForce, Counter::kIsRate | Counter::kInvert,
                                                     Counter::OneK::kIs1000);;
    state.counters["# of Triplets"] = avgDist;
}

template <typename FunctorType>
void registerOneBenchmark(const std::string& functorName, const std::string& modeName, FunctorMode functorMode,
                          const BenchmarkConfig& config)
{
    benchmark::RegisterBenchmark(
            "BM_" + functorName + "_" + modeName,
            [=](benchmark::State& state) { BM_Functor<FunctorType>(state, functorMode, config.newton3); })
        ->RangeMultiplier(2)->Ranges({{config.minParticles, config.maxParticles}, {config.cellSize, config.cellSize}, {config.cutoff, config.cutoff}});
}


void RegisterFunctorBenchmarks(const BenchmarkConfig& config)
{

    std::cout << "==========================================" << std::endl;
    std::cout << "AutoPas Functor Benchmark" << std::endl;
    std::cout << "AutoPas Branch: " << AUTOPAS_BRANCH    << std::endl;
    std::cout << "AutoPas Commit: " << AUTOPAS_COMMIT    << std::endl;
    std::cout << "==========================================" << std::endl;

    // some constants that define the benchmark
    constexpr bool mixing{false};
    constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
    constexpr bool globals{false};

    using ATM = mdLib::AxilrodTellerMutoFunctor<Particle, mixing, functorN3Modes, globals>;
    using ATMGlobals = mdLib::AxilrodTellerMutoFunctor<Particle, mixing, functorN3Modes, true>;

    const std::array modes = {
        std::make_pair("AoS", AOS),
        std::make_pair("SoASingle", SOASINGLE),
        std::make_pair("SoAPair", SOAPAIR),
        std::make_pair("SoATriple", SOATRIPLE)
    };

    for (const auto& [modeName, mode] : modes)
    {
        if (config.targetMode != "all" && config.targetMode != modeName)
        {
            continue;
        }
        if (config.targetFunctor == "all" || config.targetFunctor == "ATM")
        {
            registerOneBenchmark<ATM>("ATM", modeName, mode, config);
        }
        if (config.targetFunctor == "all" || config.targetFunctor == "ATMGlobals") {
            registerOneBenchmark<ATMGlobals>("ATMGlobals", modeName, mode, config);
        }
    }
}

void setupCLI(CLI::App& app, BenchmarkConfig& config) {
    app.add_option("--min", config.minParticles, "Minimum number of particles")->default_val(1);
    app.add_option("--max", config.maxParticles, "Maximum number of particles")->default_val(512);
    app.add_option("-c,--cell-size", config.cellSize, "Size of the simulation cell")->default_val(3);
    app.add_option("-r,--cutoff", config.cutoff, "Cutoff radius for interactions")->default_val(3);

    app.add_option("-f,--functor", config.targetFunctor, "Which functor to test")
       ->check(CLI::IsMember({"ATM", "ATMGlobals", "all"}, CLI::ignore_case))
       ->default_val("all");

    app.add_option("-m,--mode", config.targetMode, "Which data layout mode to test")
       ->check(CLI::IsMember({"AoS", "SoASingle", "SoAPair", "SoATriple", "all"}, CLI::ignore_case))
       ->default_val("all");

    app.add_flag("--n3,!--no-n3", config.newton3, "Enable/Disable Newton3 (enabled by default)");
}

bool handleHelpFlag(int argc, char** argv, const CLI::App& app) {
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


    // Read CLI arguments
    CLI::App app{"AutoPas 3-Body Functor Benchmark"};
    BenchmarkConfig config;
    setupCLI(app, config);

    if (handleHelpFlag(argc, argv, app)) {
        benchmark::Initialize(&argc, argv);
        return 0;
    }

    benchmark::Initialize(&argc, argv);
    CLI11_PARSE(app, argc, argv);

    RegisterFunctorBenchmarks(config);

    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
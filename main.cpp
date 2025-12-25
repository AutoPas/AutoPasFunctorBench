#include <iostream>

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h>
#include "benchmark/benchmark.h"

#include <autopas/cells/FullParticleCell.h>
#include <random>

// type aliases for ease of use
using Particle = mdLib::MoleculeLJ;
using Cell = autopas::FullParticleCell<Particle>;

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

    for (auto _ : state)
    {
        state.PauseTiming();
        std::vector<Cell> cells{3};
        generateParticles(functor, cells, numParticles, cellSize, functorMode);
        state.ResumeTiming();

        applyFunctorOnParticles(functor, cells, functorMode, newton3);

        state.PauseTiming();
        const auto [calcsDist, calcsForce] = countInteractions(cells, cutoff, functorMode);
        calcsDistTotal += calcsDist;
        calcsForceTotal += calcsForce;
        state.ResumeTiming();
    }
    // Per-iteration averages and hit rate as user counters.
    const auto iters = static_cast<double>(state.iterations());
    const double avgDist = static_cast<double>(calcsDistTotal) / iters;
    const double avgForce = static_cast<double>(calcsForceTotal) / iters;
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
                          bool newton3)
{
    benchmark::RegisterBenchmark(
            "BM_" + functorName + "_" + modeName,
            [=](benchmark::State& state) { BM_Functor<FunctorType>(state, functorMode, newton3); })
        ->RangeMultiplier(2)->Ranges({{1, 512}, {3, 3}, {3, 3}});
}


void RegisterFunctorBenchmarks()
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
    constexpr bool newton3{true};

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
        registerOneBenchmark<ATM>("ATM", modeName, mode, newton3);
        registerOneBenchmark<ATMGlobals>("ATMGlobals", modeName, mode, newton3);
    }
}

// before BENCHMARK_MAIN()
static bool registerAll = (RegisterFunctorBenchmarks(), true);
int main(int argc, char** argv) {
    // Add the version info to the JSON metadata
    benchmark::AddCustomContext("autopas_branch", AUTOPAS_BRANCH);
    benchmark::AddCustomContext("autopas_commit", AUTOPAS_COMMIT);

    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
#include <iostream>

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h>
#include "benchmark/benchmark.h"

#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/Timer.h>
#include <fstream>
#include <random>

// type aliases for ease of use
using Particle = mdLib::MoleculeLJ;
using Cell = autopas::FullParticleCell<Particle>;

// some constants that define the benchmark
constexpr bool mixing{false};
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool newton3{true};
constexpr bool globals{false};

// using ATM = mdLib::AxilrodTellerFunctor<Particle, false, functorN3Modes, globals>;
using ATM = mdLib::AxilrodTellerMutoFunctor<Particle, false, functorN3Modes, globals>;

enum FunctorMode {
    AOS,
    SOASINGLE,
    SOAPAIR,
    SOATRIPLE
};

double distSquared(std::array<double, 3> a, std::array<double, 3> b) {
    using autopas::utils::ArrayMath::sub;
    using autopas::utils::ArrayMath::dot;
    const auto c = sub(a, b);   // 3 FLOPS
    return dot(c, c);           // 3+2=5 FLOPs
}

void generateParticles(ATM &functor, std::vector<Cell> &cells, const size_t numberOfParticles,
    const double cellSize, const FunctorMode mode){
    // generate randomly distributed particles
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, cellSize);

    auto fillCellWithParticles = [&] (Cell& cell, double xShift, double yShift, double zShift) {
        for (size_t particleId = 0; particleId < numberOfParticles; ++particleId) {
            Particle p{
                { dis(gen) + xShift, dis(gen) + yShift, dis(gen)},
                {0., 0., 0.,},
                particleId,
                0};
            cell.addParticle(p); // Add particle to the current cell
        }
    };

    switch (mode) {
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

void applyAoSFunctor(ATM &functor, Cell &cell) {
    for(std::size_t i = 0; i < cell.size(); ++i) {
        for (std::size_t j = i + 1; j < cell.size(); ++j) {
            for (std::size_t k = j + 1; k < cell.size(); ++k) {
                // timer.at("AoSFunctor on Particles").start();
                functor.AoSFunctor(cell[i], cell[j], cell[k], newton3);
                // timer.at("AoSFunctor on Particles").stop();
            }
        }
    }
}

void applyFunctorOnParticles(ATM &functor, std::vector<Cell> &cells, const FunctorMode mode) {
    switch (mode) {
        case AOS:
            applyAoSFunctor(functor, cells[0]);
            break;
        case SOASINGLE:
            functor.SoAFunctorSingle(cells[0]._particleSoABuffer, newton3);
            break;
        case SOAPAIR:
            functor.SoAFunctorPair(cells[0]._particleSoABuffer, cells[1]._particleSoABuffer, newton3);
            break;
        case SOATRIPLE:
            functor.SoAFunctorTriple(cells[0]._particleSoABuffer, cells[1]._particleSoABuffer, cells[2]._particleSoABuffer, newton3);
            break;
    }
}

// void csvOutput(ATM &functor, std::vector<Cell> &cells) {
//     timer.at("Output").start();
//     std::ofstream csvFile("particles.csv");
//     if (not csvFile.is_open()) {
//         throw std::runtime_error("FILE NOT OPEN!");
//     }
//     for (size_t cellId = 0; cellId < cells.size(); ++cellId) {
//         functor.SoAExtractor(cells[cellId], cells[cellId]._particleSoABuffer, 0);
//         csvFile << "ParticleId,rX,rY,rZ,fX,fY,fZ\n";
//         for (auto p : cells[0]) {
//             using autopas::utils::ArrayUtils::to_string;
//             csvFile << p.getID() << ","
//                     << to_string(p.getR(), ",", {"", ""}) << ","
//                     << to_string(p.getF(), ",", {"", ""})
//                     << "\n";
//         }
//         csvFile.close();
//         timer.at("Output").stop();
//     }
// }

std::tuple<size_t, size_t> countInteractions(std::vector<Cell> &cells, const double cutoff, const FunctorMode mode) {
    size_t calcsDist{0};
    size_t calcsForce{0};
    const auto cutoffSquared{cutoff * cutoff};

    switch (mode) {
        case AOS:
        case SOASINGLE:
            for (size_t i = 0; i < cells[0].size(); i++) {
                for (size_t j = i + 1; j < cells[0].size(); j++) {
                    if (distSquared(cells[0][i].getR(), cells[0][j].getR()) > cutoffSquared) continue;
                    for (size_t k = j + 1; k < cells[0].size(); k++) {
                        if (distSquared(cells[0][i].getR(), cells[0][k].getR()) > cutoffSquared) continue;
                        if (distSquared(cells[0][j].getR(), cells[0][k].getR()) > cutoffSquared) continue;
                        ++calcsForce;
                    }
                }
            }
            calcsDist = cells[0].size() * (cells[0].size() - 1) * (cells[0].size() - 2) / 6;
            break;
        case SOAPAIR:
            for (size_t i = 0; i < cells[0].size(); i++) {
                for (size_t j = i + 1; j < cells[0].size(); j++) {
                    if (distSquared(cells[0][i].getR(), cells[0][j].getR()) > cutoffSquared) continue;
                    for (size_t k = 0; k < cells[1].size(); k++) {
                        if (distSquared(cells[0][i].getR(), cells[1][k].getR()) > cutoffSquared) continue;
                        if (distSquared(cells[0][j].getR(), cells[1][k].getR()) > cutoffSquared) continue;
                        ++calcsForce;
                    }
                }
                for (size_t j = 0; j < cells[1].size(); j++) {
                    if (distSquared(cells[0][i].getR(), cells[1][j].getR()) > cutoffSquared) continue;
                    for (size_t k = j + 1; k < cells[1].size(); k++) {
                        if (distSquared(cells[0][i].getR(), cells[1][k].getR()) > cutoffSquared) continue;
                        if (distSquared(cells[1][j].getR(), cells[1][k].getR()) > cutoffSquared) continue;
                        ++calcsForce;
                    }
                }
            }
            calcsDist = cells[0].size() * cells[1].size() * (cells[0].size() + cells[1].size() - 2) / 2;
            break;
        case SOATRIPLE:
            for (size_t i = 0; i < cells[0].size(); i++) {
                for (size_t j = 0; j < cells[1].size(); j++) {
                    if (distSquared(cells[0][i].getR(), cells[1][j].getR()) > cutoffSquared) continue;
                    for (size_t k = 0; k < cells[2].size(); k++) {
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

class FunctorTestFixture : public benchmark::Fixture
{
public:
    void SetUp(const benchmark::State&) override
    {

    }
    void TearDown(const benchmark::State&) override
    {

    }

    const double cellSize{3.};
    const double cutoff{3.};
    const double nu{1.0};

    const size_t numParticles{100};
    const size_t iterations{1};
    const std::array<FunctorMode, 4> functorsToTest = {
        AOS,
        SOASINGLE,
        SOAPAIR,
        SOATRIPLE
    };

    void runBenchmarkForFunctorMode(benchmark::State& state, FunctorMode functorMode) const
    {
        ATM functor{cutoff};
        functor.setParticleProperties(nu);

        std::size_t calcsDistTotal = 0;
        std::size_t calcsForceTotal = 0;

        for (auto _ : state)
        {
            state.PauseTiming();
            std::vector<Cell> cells{3};
            generateParticles(functor, cells, numParticles, cellSize, functorMode);
            state.ResumeTiming();

            applyFunctorOnParticles(functor, cells, functorMode);

            state.PauseTiming();
            const auto [calcsDist, calcsForce] = countInteractions(cells, cutoff, functorMode);
            calcsDistTotal += calcsDist;
            calcsForceTotal += calcsForce;
            state.ResumeTiming();
        }
        // Per-iteration averages and hit rate as user counters.
        const double iters = static_cast<double>(state.iterations());
        const double avgDist = static_cast<double>(calcsDistTotal) / iters;
        const double avgForce = static_cast<double>(calcsForceTotal) / iters;
        const double hitRate = avgForce / avgDist * 100.0;

        using benchmark::Counter;
        state.counters["HitRate [%]"] = hitRate;
        state.counters["Time per Triplet"] = Counter(avgDist, Counter::kIsRate | Counter::kInvert, Counter::OneK::kIs1000);
        state.counters["Time per Interaction"] = Counter(avgForce, Counter::kIsRate | Counter::kInvert, Counter::OneK::kIs1000);;
    }
};

BENCHMARK_F(FunctorTestFixture, AoS)(benchmark::State& state) {
    runBenchmarkForFunctorMode(state, FunctorMode::AOS);
}

BENCHMARK_F(FunctorTestFixture, SoASingle)(benchmark::State& state) {
    runBenchmarkForFunctorMode(state, FunctorMode::SOASINGLE);
}

BENCHMARK_F(FunctorTestFixture, SoAPair)(benchmark::State& state) {
    runBenchmarkForFunctorMode(state, FunctorMode::SOAPAIR);
}

BENCHMARK_F(FunctorTestFixture, SoATriple)(benchmark::State& state) {
    runBenchmarkForFunctorMode(state, FunctorMode::SOATRIPLE);
}

BENCHMARK_MAIN();

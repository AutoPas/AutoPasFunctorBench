#include <iostream>

#ifdef ENABLE_FAPP
#include "fj_tool/fapp.h"
#endif

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <molecularDynamicsLibrary/ArgonFunctor.h>
#include <molecularDynamicsLibrary/AxilrodTellerFunctor.h>

#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/Timer.h>
#include <fstream>

// type aliases for ease of use
using Particle = mdLib::MoleculeLJ;
// some constants that define the benchmark
constexpr bool shift{false};
constexpr bool mixing{false};
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool newton3{true};
constexpr bool globals{false};

using ArgonFunctor = mdLib::ArgonFunctor<Particle, functorN3Modes, globals>;
using ATMFunctor = mdLib::AxilrodTellerFunctor<Particle, mixing, functorN3Modes, globals>;

double distSquared(std::array<double, 3> a, std::array<double, 3> b) {
    using autopas::utils::ArrayMath::sub;
    using autopas::utils::ArrayMath::dot;
    const auto c = sub(a, b);   // 3 FLOPS
    return dot(c, c);           // 3+2=5 FLOPs
}

std::map<std::string, autopas::utils::Timer> timer{
        {"FunctorArgon",                 autopas::utils::Timer()},
        {"FunctorATM",                 autopas::utils::Timer()}
};

void resetTimer() {
    for (const auto &[name, t]: timer) {
        timer[name].reset();
    }
}

void printTimer() {
    for (const auto &[name, t]: timer) {
        std::cout
                << std::setw(18)
                << std::left
                << name
                << " : "
                << std::setprecision(3)
                << std::setw(8)
                << static_cast<double>(timer[name].getTotalTime()) * 1e-9
                << " [s]\n";
    }
}

std::vector<Particle> generateParticles(const size_t numberOfParticles, double cutoff){
    std::vector<Particle> particles;
    for (size_t particleId = 0; particleId < numberOfParticles; ++particleId) {
        Particle p{
                {
                        rand() / static_cast<double>(RAND_MAX) * cutoff,
                        rand() / static_cast<double>(RAND_MAX) * cutoff,
                        rand() / static_cast<double>(RAND_MAX) * cutoff,
                },
                {0., 0., 0.,},
                // every cell gets its own id space
                particleId,
                0};
        particles.push_back(p);
    }
    return particles;
}

void applyFunctorOnTriplet(ArgonFunctor &functor, Particle &i, Particle &j, Particle &k) {
    timer.at("FunctorArgon").start();
    functor.AoSFunctor(i, j, k, newton3);
    timer.at("FunctorArgon").stop();
}

void applyFunctorOnTriplet(ATMFunctor &functor, Particle &i, Particle &j, Particle &k) {
    timer.at("FunctorATM").start();
    functor.AoSFunctor(i, j, k, newton3);
    timer.at("FunctorATM").stop();
}

template<class FUNCTOR>
void applyFunctorOnParticleSet(FUNCTOR &functor, std::vector<Particle> &particles) {
    for(std::size_t i = 0; i < particles.size(); ++i) {
        for (std::size_t j = i + 1; j < particles.size(); ++j) {
            for (std::size_t k = j + 1; k < particles.size(); ++k) {
                applyFunctorOnTriplet(functor, particles[i], particles[j], particles[k]);
            }
        }
    }
}

std::tuple<size_t, size_t> countInteractions(std::vector<Particle> &particles, double cutoff) {
    size_t calcsDist{0};
    size_t calcsForce{0};
    const auto cutoffSquared{cutoff * cutoff};
    for (const auto &p0: particles) {
        for (const auto &p1: particles) {
            for (const auto &p2: particles) {
                if (p0 != p1 && p0 != p2 && p1 != p2) {
                    ++calcsDist;
                    if (distSquared(p0.getR(), p1.getR()) <= cutoffSquared) {
                        ++calcsForce;
                    }
                }
            }
        }
    }
    return {calcsDist, calcsForce};
}

/**
 * Mini benchmark tool to estimate the inner most kernel performance of AutoPas
 * @return
 */
int main() {

    // Open the file
    std::ofstream file("../Benchmark_Cutoff.csv");

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cout << "Error opening file!" << std::endl;
        return 1;
    }

    // Write headers
    file << "CutoffRadius,Particles,FunctionCalls,Interactions,HitRate,FunctorArgonRuntime,FunctorATMRuntime" << std::endl;

    constexpr double cutoff{3.};


    // define scenario
    constexpr size_t numParticles{100};
    // generate particles at random positions
    std::vector<Particle> particles{generateParticles(numParticles, cutoff)};

    const std::vector<double> radiuses{{1, 1.35, 1.6, 1.82, 2.03, 2.23, 2.43, 2.7, 3., 5.}};
    for (auto i=0 ; i<100; i++) {
        for (auto radius : radiuses) {
            // choose functor based on available architecture
            ArgonFunctor functorArgon{radius};
            ATMFunctor functorATM{radius};

            // actual benchmark
            applyFunctorOnParticleSet(functorArgon, particles);
            applyFunctorOnParticleSet(functorATM, particles);
            // gather data for analysis
            const auto [calcsDistTotal, calcsForceTotal] = countInteractions(particles, radius);

            //using autopas::utils::ArrayUtils::operator<<;

            file
                    << radius << ","
                    << numParticles << ","
                    << calcsDistTotal << ","
                    << calcsForceTotal << ","
                    << (static_cast<double>(calcsForceTotal) / calcsDistTotal) << ","
                    << static_cast<double>(timer["FunctorArgon"].getTotalTime()) * 1e-9 << ","
                    << static_cast<double>(timer["FunctorATM"].getTotalTime()) * 1e-9<< std::endl;

            resetTimer();
        }
    }

    file.close();

}
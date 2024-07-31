#include <iostream>

// Include the functor(s) to test
#include <molecularDynamicsLibrary/KryptonPairFunctor.h>

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/Timer.h>
#include <autopas/utils/ArrayMath.h>
#include <fstream>

// type aliases for ease of use
using Particle = mdLib::MoleculeLJ;
using Cell = autopas::FullParticleCell<mdLib::MoleculeLJ>;
// some constants that define the benchmark
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool newton3{true};
constexpr bool globals{true};

// Chose the Functor
using Functor = mdLib::KryptonPairFunctor<Particle, false, functorN3Modes, globals>;
using FunctorPreCompute = mdLib::KryptonPairFunctor<Particle, true, functorN3Modes, globals>;

/**
 * Mini validation setup that writes functor results to a CSV
 * @return
 */
int main(int argc, char* argv[]) {

    const double cutoff{300.};
    Functor functor{cutoff};
    FunctorPreCompute functorPreCompute{cutoff};

    auto particle1 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
    auto particle2 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);

    auto particle3 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
    auto particle4 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);

    // Open the file
    std::ofstream file("kryptonPair2.csv");

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cout << "Error opening file!" << std::endl;
        return 1;
    }

    // Write headers
    file << "Distance [nm],Force [K/nm],Energy [K]" << std::endl;

    double epot{};
    double virial{};
    size_t resolution = 1000;
    double increment = 7.0 / resolution;
//    double distance = 1.10798137169225 * 0.35;
    double distance = 3.0;
    for (auto i = 0; i < resolution; i++) {
        // move second particle along x-Axis
        particle1.setF({0.,0.,0.});
        particle2.setF({0.,0.,0.});
        particle2.setR({distance, 0., 0.});
        // Compute the interactions
        functor.initTraversal();
        functor.AoSFunctor(particle1, particle2, newton3);
        functor.endTraversal(newton3);

        // Results
        auto force = particle1.getF();
        epot = functor.getPotentialEnergy();
        virial = functor.getVirial();

        // std::cout << "Distance,   " << distance << std::endl;
        // std::cout << "Force,      " << force[0] << std::endl;
        // std::cout << "Potential Energy:" << epot << std::endl;

        // Write to CSV
        file << std::fixed << std::setprecision(15) << distance << "," << force[0] << "," << epot << "," << virial << std::endl;

        // Compute the interactions 2
        particle3.setF({0.,0.,0.});
        particle4.setF({0.,0.,0.});
        particle4.setR({distance, 0., 0.});

        functorPreCompute.initTraversal();
        functorPreCompute.AoSFunctor(particle3, particle4, newton3);
        functorPreCompute.endTraversal(newton3);
        auto forcePreCompute = particle3.getF();

        if (force != forcePreCompute) {
            std::cout << "Different Force for r = " << distance <<
            ", F1 = " << force[0] << ", F2 = " << forcePreCompute[0] <<
            " => DIFF: " << std::abs(force[0] - forcePreCompute[0]) << std::endl;
        }
        distance += increment;
        break;
    }

    file.close();
}

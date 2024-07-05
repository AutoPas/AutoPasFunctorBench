#include <iostream>


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

using Functor = mdLib::KryptonPairFunctor<Particle, false, functorN3Modes, globals>;

void printHelp() {
    std::cout << "Usage: \n"
        << "./AutoPasFunctorBench <functor type> <repetitions> <iterations> <numParticles> <hit rate> <outfile>\n"
        << "Possible values:\n"
        << "functor type: \"aos\", \"single\", \"pair\", \"verlet\"\n"
        << "repetitions, iterations, numParticles: int\n"
        << "hit rate: double\n"
        << "outfile: std::string"
        << std::endl;
}

/**
 * Mini benchmark tool to estimate the inner most kernel performance of AutoPas
 * @return
 */
int main(int argc, char* argv[]) {

    const double cutoff{300.};
    Functor functor{cutoff};
//    std::vector<Particle> {
//        Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0),
//        Particle({0.3, 0., 0.}, {0., 0., 0.}, 1, 0)
//    };
    auto particle1 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
    auto particle2 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);

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
    double increment = 0.7 / resolution;
    double distance = 1.10798137169225 * 0.35;
    while (distance < 15.51) {
        epot = 0.0;
        particle2.setF({0.,0.,0.});
        particle2.setR({distance, 0., 0.});
        functor.initTraversal();
        functor.AoSFunctor(particle1, particle2, newton3);
        functor.endTraversal(newton3);
        auto force = particle1.getF();
        epot = functor.getPotentialEnergy();
        virial = functor.getVirial();
        // std::cout << "Distance,   " << distance << std::endl;
        // std::cout << "Force,      " << force[0] << std::endl;
        // std::cout << "Potential Energy:" << epot << std::endl;

        file << std::fixed << std::setprecision(15) << distance << "," << force[0] << "," << epot << "," << virial << std::endl;
        distance += 0.1;
        break;
    }

    file.close();
}

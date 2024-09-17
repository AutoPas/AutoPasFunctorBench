#include <iostream>
#include <unistd.h> // for getopt
#include <fstream>
#include <variant>

// Include the functor(s) to test
#include <molecularDynamicsLibrary/AxilrodTellerFunctor.h>
#include <molecularDynamicsLibrary/KryptonExtendedATMFunctor.h>
#include <molecularDynamicsLibrary/KryptonPairFunctor.h>
#include <molecularDynamicsLibrary/LJFunctor.h>

#include <molecularDynamicsLibrary/MoleculeLJ.h>
#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/Timer.h>
#include <autopas/utils/ArrayMath.h>

// type aliases for ease of use
using Particle = mdLib::MoleculeLJ;
using Cell = autopas::FullParticleCell<mdLib::MoleculeLJ>;
// some constants that define the experiments
constexpr autopas::FunctorN3Modes functorN3Modes{autopas::FunctorN3Modes::Both};
constexpr bool newton3{true};
constexpr bool globals{true};

// Chose the Functors
using KrEATMFunctor = mdLib::KryptonExtendedATMFunctor<Particle, functorN3Modes, globals>;
using ATMFunctor = mdLib::AxilrodTellerFunctor<Particle, false, functorN3Modes, globals>;
using KrPairFunctor = mdLib::KryptonPairFunctor<Particle, functorN3Modes, globals>;
using LJFunctor = mdLib::LJFunctor<mdLib::MoleculeLJ, true, false, functorN3Modes, true, true>;

enum class FunctorOption {
    LJFunctor,
    KrPairFunctor,
    ATMFunctor,
    KrEATMFunctor
};

enum class GeometryOption {
    linear,
    equilateral,
    isosceles
};

std::string getCSVFileName(FunctorOption functorOption, GeometryOption geometryOption) {
    std::stringstream filenameStream;
    switch (functorOption) {
        case FunctorOption::LJFunctor:
            std::cout << "Using the Lennard-Jones Functor\n";
            filenameStream << "kr-2b-lj";
        break;
        case FunctorOption::KrPairFunctor:
            std::cout << "Using the Ab-initio pairwise Krypton Functor\n";
            filenameStream << "kr-2b-abinitio";
        break;
        case FunctorOption::ATMFunctor:
            std::cout << "Using the Axilrod-Teller-Muto Functor\n";
            filenameStream << "kr-3b-atm";
        break;
        case FunctorOption::KrEATMFunctor:
            std::cout << "Using the extended Axilrod-Teller Functor\n";
            filenameStream << "kr-3b-eatm";
        break;
        default:
            return "";
    }

    switch (geometryOption) {
        case GeometryOption::linear:
            std::cout << "Using a linear geometry\n";
            filenameStream << "-linear";
        break;
        case GeometryOption::equilateral:
            std::cout << "Using a equilateral triangle geometry\n";
            filenameStream << "-equilateral";
        break;
        case GeometryOption::isosceles:
            std::cout << "Using a isosceles triangle geometry\n";
            filenameStream << "-isosceles";
        break;
        default:
            return "";
    }
    filenameStream << ".csv";
    return filenameStream.str();
}

void setEquilateral(Particle &p2, Particle &p3, double distance) {
    p2.setR({0., distance, 0.});
    p3.setR({std::sqrt(3)/2 * distance, 0.5 * distance, 0.});
}

void setLinear(Particle &p1, Particle &p2, Particle &p3, double distance) {
    p1.setR({10.0 + distance, 0., 0.});
    p2.setR({10.0, 0., 0.});
    p3.setR({9.8, 0., 0.});
}

void setIsosceles(Particle& p1, Particle &p2, Particle &p3, double distance) {
    p1.setR({distance, 2.0, 0.0});
    p2.setR({0.0, 0.0, 0.0});
    p3.setR({0.0, 4.0, 0.0});
}

/**
 * Mini validation setup that writes functor results to a CSV
 * @return
 */
int main(int argc, char* argv[]) {

    FunctorOption functorOption{};
    GeometryOption geometryOption{};

    if (argc < 2) {
        std::cerr << "Usage: -f <functor> [-g <geometry>]" << std::endl;
        std::cerr << "functor: 0 : Lennard-Jones (default); 1 : Ab-Initio Krypton Pair; 2 : Axilrod-Teller; 3 : extended Axilrod-Teller\n";
        std::cerr << "Geometry: 0 : linear (default); 1 : equilateral; 2 : isosceles\n";
        return 1;
    }

    int option;
    while ((option = getopt(argc, argv, "f:g:")) != -1) {
        switch (option) {
            case 'f':
                functorOption = static_cast<FunctorOption>(std::stoi(optarg));
            break;
            case 'g':
                geometryOption = static_cast<GeometryOption>(std::stoi(optarg));
            break;
            default:
                std::cerr << "Usage: -f <functor> [-g <geometry>]" << std::endl;
                std::cerr << "functor: 0 : Lennard-Jones (default); 1 : Ab-Initio Krypton Pair; 2 : Axilrod-Teller; 3 : extended Axilrod-Teller\n";
                std::cerr << "Geometry: 0 : linear (default); 1 : equilateral; 2 : isosceles\n";
            return 1;
        }
    }

    constexpr double cutoff{300.};
    constexpr double epsilon = 200.8753;
    constexpr double r_eps = 4.015802;
    const double sigma = r_eps / (std::pow(2.0, 1.0 / 6.0));
    constexpr double nu = 1.61525e6;
    std::cout << "Sigma: " << sigma << std::endl;
    auto particle1 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
    auto particle2 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);
    auto particle3 = Particle({0., 0., 0.}, {0., 0., 0.}, 0, 0);

    std::string filename = getCSVFileName(functorOption, geometryOption);

    // Create the selected functor as a std::variant
    auto createFunctor = [&] () -> std::variant<LJFunctor, KrPairFunctor, ATMFunctor, KrEATMFunctor> {
        switch (functorOption) {
        case FunctorOption::LJFunctor: {
            LJFunctor ljfunctor{cutoff};
            ljfunctor.setParticleProperties(24 * epsilon, sigma * sigma);
            return ljfunctor;
        }
        case FunctorOption::KrPairFunctor:
            return KrPairFunctor{cutoff};
        case FunctorOption::ATMFunctor: {
            ATMFunctor atmfunctor{cutoff};
            atmfunctor.setParticleProperties(nu);
            return atmfunctor;
        }
        case FunctorOption::KrEATMFunctor:
            return KrEATMFunctor{cutoff};
        }
        throw std::invalid_argument("Unknown functor option!");
    };

    auto functor = createFunctor();

    // Open the file
    std::ofstream file(filename);

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cout << "Error opening file!" << std::endl;
        return 1;
    }

    // Write headers
    file << "Distance [A],Energy [K],Force_x [K/A],Force_y [K/A],Force_z [K/A]" << std::endl;

    double distance = 2.0;
    constexpr double maxDistance = 15.0;
    constexpr double increment = 0.02;
    double epot{};

    while (distance <= maxDistance) {
        // Reset forces
        particle1.setF({0.,0.,0.});
        particle2.setF({0.,0.,0.});
        particle3.setF({0.,0.,0.});

        // Move particles
        switch (geometryOption) {
            case GeometryOption::linear: {
                setLinear(particle1, particle2, particle3, distance);
                break;
            }
            case GeometryOption::equilateral: {
                setEquilateral(particle2, particle3, distance);
                break;
            }
            case GeometryOption::isosceles: {
                setIsosceles(particle1, particle2, particle3, distance);
                break;
            }
        }

        // Do one functor interaction (AoS)
        std::visit([&](auto&& f) {
            f.initTraversal();

            using T = std::decay_t<decltype(f)>;
            if constexpr (std::is_same_v<T, LJFunctor> || std::is_same_v<T, KrPairFunctor>) {
                // Call AoSFunctor with two particles for LJFunctor and KrPairFunctor
                f.AoSFunctor(particle1, particle2, newton3);
            } else if constexpr (std::is_same_v<T, ATMFunctor> || std::is_same_v<T, KrEATMFunctor>) {
                // Call AoSFunctor with three particles for ATMFunctor and KrEATMFunctor
                f.AoSFunctor(particle1, particle2, particle3, newton3);
            }

            f.endTraversal(newton3);
            epot = f.getPotentialEnergy();
        }, functor);

        // Results
        auto force = particle1.getF();

        // Write to CSV
        file << distance << std::fixed << std::setprecision(15) << "," << epot << "," << force[0] << "," << force[1] << "," << force[2]  << std::endl;
        distance += increment;
    }

    file.close();
    return 0;
}

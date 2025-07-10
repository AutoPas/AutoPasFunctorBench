#include <iostream>

#include <filesystem>
// Include the functor(s) to test
#include <molecularDynamicsLibrary/AxilrodTellerFunctor.h>
#include <molecularDynamicsLibrary/LJFunctor.h>
//#include <molecularDynamicsLibrary/ArgonFunctor.h>
#include <molecularDynamicsLibrary/KryptonExtendedATMFunctor.h>
#include <molecularDynamicsLibrary/LUT3B.h>

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
constexpr bool flops{true};
constexpr bool useLUTglob{false};
namespace mdFlexibleTypeDefs {
/**
 * If AutoPas is compiled with FLOP logging enabled, use functors with FLOP counting enabled.
 */
    constexpr bool countFLOPs =

            false;


/**
 * If md-flexible is compiled with globals calculations enabled, use functors which calculate globals.
 */
    constexpr bool calcGlobals =

            true;

}  // namespace mdFlexibleTypeDefs
// Chose the Functor
//using Functor = mdLib::AxilrodTellerFunctor<Particle, false, functorN3Modes, globals>;
//using Functor = mdLib::AxilrodTellerFunctor<Particle, false,true, functorN3Modes,globals, flops>


//WORKS FOR feat/3xa/lut
//using Functor = mdLib::AxilrodTellerFunctor<mdLib::MoleculeLJ, false, false, autopas::FunctorN3Modes::Both,mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;


//for auto-personal
//using Functor = mdLib::AxilrodTellerFunctor<mdLib::MoleculeLJ, false, true,false, autopas::FunctorN3Modes::Both,mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
//using Functor = mdLib::LJFunctor<mdLib::MoleculeLJ, false, false,true, useLUTglob, autopas::FunctorN3Modes::Both, mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs, true>;


//for auto-personal argon
//using Functor  =  mdLib::ArgonFunctor<mdLib::MoleculeLJ, autopas::FunctorN3Modes::Both,true, false,mdFlexibleTypeDefs::calcGlobals >;


//for auto-personal krypton
using Functor =  mdLib::KryptonExtendedATMFunctor<mdLib::MoleculeLJ, autopas::FunctorN3Modes::Both, true, false, mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
using FunctorNoLUT =  mdLib::KryptonExtendedATMFunctor<mdLib::MoleculeLJ, autopas::FunctorN3Modes::Both, false, false, mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;

// noble gas krypton
//using Functor =  mdLib::KryptonExtendedATMFunctor<mdLib::MoleculeLJ, autopas::FunctorN3Modes::Both, mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;



//for Force_LUT
//using Functor = mdLib::AxilrodTellerFunctor<Particle, true,false, autopas::FunctorN3Modes::Both, false>
//using Functor = mdLib::AxilrodTellerFunctor<Particle, false, true, autopas::FunctorN3Modes::Both, true>;
void setEquilateral(Particle &p2, Particle &p3, double distance) {
    p2.
    setR({0., distance, 0.});
    p3.setR({std::sqrt(3)/2 * distance, 0.5 * distance, 0.});
}

void setLinear(Particle &p1, Particle &p2, Particle &p3, double distance) {
    p1.setR({3.0 + distance, 0., 0.});
    p2.setR({3.0, 0., 0.});
    p3.setR({0., 0., 0.});
}

void setIsosceles(Particle& p1, Particle &p2, Particle &p3, double distance) {
    p1.setR({distance, 2.0, 0.0});
    p2.setR({0.0, 0.0, 0.0});
    p3.setR({0.0, 4.0, 0.0});}

/**
 * Mini validation setup that writes functor results to a CSV
 * @return
 */
int main(int argc, char* argv[]) {
    const auto cutoff = 5.0;
    auto lut3B = mdLib::LUT3B(0, cutoff * cutoff);

    Functor functor{cutoff, &lut3B};
    FunctorNoLUT functorNoLut{cutoff};

    // Use high resolution to check if LUT converges
    lut3B.setResolution(400);
    lut3B.fill_plain(functor);

    // Particle Setup
    auto particle1 = Particle({0., 0., 0.5}, {0., 0., 0.}, 0, 0);
    auto particle2 = Particle({2.,3., 0.}, {0., 0., 0.}, 0, 0);
    auto particle3 = Particle({2., 0., 0.}, {0., 0., 0.}, 0, 0);

    // Particle Setup for reference
    auto particle1NoLut = Particle(particle1);
    auto particle2NoLut = Particle(particle2);
    auto particle3NoLut = Particle(particle3);

    for (auto i = 0; i < 1; i++) {
        // set Forces to 0
        particle1.setF({0.,0.,0.});
        particle2.setF({0.,0.,0.});
        particle3.setF({0.,0.,0.});

        particle1NoLut.setF({0.,0.,0.});
        particle2NoLut.setF({0.,0.,0.});
        particle3NoLut.setF({0.,0.,0.});

        // Compute the interactions
        functor.initTraversal();
        functorNoLut.initTraversal();

        // Change order as needed
        functor.AoSFunctor(particle2, particle1, particle3, newton3);
        functorNoLut.AoSFunctor(particle2NoLut, particle1NoLut, particle3NoLut, newton3);

        functor.endTraversal(newton3);
        functorNoLut.endTraversal(newton3);

        std::cout << "\nParticle1 Force:      " << particle1.getF()[0] << " " << particle1.getF()[1] << " " << particle1.getF()[2] << std::endl;
        std::cout << "Particle1NoLut Force: " << particle1NoLut.getF()[0] << " " << particle1NoLut.getF()[1] << " " << particle1NoLut.getF()[2] << std::endl;

        std::cout << "\nParticle2 Force:      " << particle2.getF()[0] << " " << particle2.getF()[1] << " " << particle2.getF()[2] << std::endl;
        std::cout << "Particle2NoLut Force: " << particle2NoLut.getF()[0] << " " << particle2NoLut.getF()[1] << " " << particle2NoLut.getF()[2] << std::endl;

        std::cout << "\nParticle3 Force:      " << particle3.getF()[0] << " " << particle3.getF()[1] << " " << particle3.getF()[2] << std::endl;
        std::cout << "Particle3NoLut Force: " << particle3NoLut.getF()[0] << " " << particle3NoLut.getF()[1] << " " << particle3NoLut.getF()[2] << std::endl;
    }

}

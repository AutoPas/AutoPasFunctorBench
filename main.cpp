#include <iostream>

#include <molecularDynamicsLibrary/ArgonInclude/RepulsiveTerm.h>
#include <molecularDynamicsLibrary/ArgonInclude/DispersionTerm.h>

#include <autopas/cells/FullParticleCell.h>
#include <autopas/utils/ArrayMath.h>
#include <fstream>
//Atomic units
/*
static constexpr std::array<double, 23> A{
        {-0.170806873130E+01,  // 000
         -0.316818997395E+02,  // 001
         -0.571545817498E+05,  // 011
         0.848780677578E+02,   // 111
         0.163923794220E+07,   // 002
         0.380809830366E+02,   // 012
         -0.217403993198E+03,  // 112
         0.244997545538E+03,   // 022
         0.128926029735E+03,   // 122
         0.815601247450E+02,   // 222
         0.409987725022E+02,   // 003
         -0.978512983041E+06,  // 013
         0.104383189893E+07,   // 113
         -0.383820796134E+02,  // 023
         0.143934125087E+03,   // 123
         0.102161665959E+04,   // 033
         -0.569593762549E+02,  // 004
         0.178356269769E+04,   // 014
         0.242202158097E+02,   // 114
         -0.279617357863E+01,  // 024
         -0.324585542907E+02,  // 005
         -0.963264559888E-01,  // 015
         -0.898942588279E+05}  // 006
};

static constexpr std::array<double, 23> alpha{
        {0.428132039316E+00,  // 000
         0.503934786518E+00,  // 001
         0.104706730543E+01,  // 011
         0.456769339560E+00,  // 111
         0.131047310452E+01,  // 002
         0.444052360076E+00,  // 012
         0.480469535570E+00,  // 112
         0.737327026170E+00,  // 022
         0.496177745527E+00,  // 122
         0.424365319847E+00,  // 222
         0.428946186456E+00,  // 003
         0.117979281352E+01,  // 013
         0.119534448663E+01,  // 113
         0.416753172892E+00,  // 023
         0.507114743788E+00,  // 123
         0.764351644551E+00,  // 033
         0.422619330972E+00,  // 004
         0.757543022081E+00,  // 014
         0.482734248672E+00,  // 114
         0.419340374650E+00,  // 024
         0.635761316281E+00,  // 005
         0.375600311119E+00,  // 015
         0.130334333132E+01}  // 006
};

static constexpr std::array<double, 5> Z{
        {0.273486414323E+03,   // 111
         -0.213475877256E+05,  // 112
         0.108226781130E+07,   // 122
         -0.213710093072E+07,  // 222
         0.364515182541E+06}   // 113
};

static constexpr std::array<double, 5> beta{
        {0.211602562917E+02,  // 111
         0.149623190559E+01,  // 112
         0.132161541056E+01,  // 122
         0.208199482789E+01,  // 222
         0.179870559008E+01}  // 113
};
*/

//SI units (nm and K)
static constexpr std::array<long double, 23> A{
        {-5.39365445993314E+05,  // 000
         -1.00043526760807E+07,  // 001
         -1.80479894697093E+10,  // 011
         2.68023739515408E+07,   // 111
         5.17630401857978E+11,   // 002
         1.20250233629457E+07,   // 012
         -6.86507513446023E+07,  // 112
         7.73641060191982E+07,   // 022
         4.07116202374599E+07,   // 122
         2.57546504143754E+07,   // 222
         1.29463884038186E+07,   // 003
         -3.08989961490876E+11,  // 013
         3.29616043775900E+11,   // 113
         -1.21201021419532E+07,  // 023
         4.54508019194995E+07,   // 123
         3.22601026022283E+08,   // 033
         -1.79863484497154E+07,  // 004
         5.63204555102674E+08,   // 014
         7.64813924806795E+06,   // 114
         -8.82961781148373E+05,  // 024
         -1.02496007862500E+07,  // 005
         -3.04174890291515E+04,  // 015
         -2.83863618111236E+10}  // 006
};

static constexpr std::array<long double, 23> alpha{
        {8.09052299484753E+00,  // 000
         9.52298731190775E+00,  // 001
         1.97867044131258E+01,  // 011
         8.63168953894591E+00,  // 111
         2.47643526123088E+02,  // 002
         8.39137345537347E+00,  // 012
         9.07955833453440E+00,  // 112
         1.39334614374608E+01,  // 022
         9.37640048180289E+00,  // 122
         8.01934231300047E+00,  // 222
         8.10590814604500E+00,  // 003
         2.22948530135448E+01,  // 013
         2.25887370431209E+01,  // 113
         7.87549358334693E+00,  // 023
         9.58307979519088E+00,  // 123
         1.44441527110870E+01,  // 033
         7.98634790509653E+00,  // 004
         1.43154883935442E+01,  // 014
         9.12235520967071E+00,  // 114
         7.92438461086463E+00,  // 024
         1.20141476840267E+01,  // 005
         7.09781720339141E+00,  // 015
         2.46296194255218E+01}  // 006
};

static constexpr std::array<double, 5> Z{
        {2.81011266190959E-04,   // 111
         -6.14241347348619E-05,  // 112
         8.72021611550457E-06,   // 122
         -4.82191783511564E-08,  // 222
         2.93702828075611E-06}   // 113
};

static constexpr std::array<double, 5> beta{
        {3.99870891182023E+02,  // 111
         2.82746852049202E+01,  // 112
         2.49749116804324E+01,  // 122
         3.93440001759947E+01,  // 222
         3.39906094408459E+01}  // 113
};

// type aliases for ease of use
    using DisplacementHandle = autopas::utils::ArrayMath::Argon::DisplacementHandle;

    namespace autopas::utils::ArrayMath{
        std::array<double, 3> rescale_array(std::array<double, 3> array, double s) {
            return array * s;
        }

        std::array<double, 3> sum_arrays(std::array<double, 3> array1, std::array<double, 3> array2) {
            return array1 + array2;
        }
    }


/**
 * Mini validation setup that writes functor results to a CSV
 * @return*
 **/
int main(int argc, char *argv[]) {

    enum ID {
        I, J, K
    };

    // equilateral triangle
    std::array<double, 3> positionI_unit{{0, 0, 0}};
    std::array<double, 3> positionJ_unit{{0, 1, 0}};
    std::array<double, 3> positionK_unit{{std::sqrt(3) / 2, 0.5, 0}};

    // symmetric linear geometry
    /*std::array<double, 3> positionI_unit{{0, 0, 0}};
    std::array<double, 3> positionJ_unit{{0.5, 0, 0}};
    std::array<double, 3> positionK_unit{{1, 0, 0}};*/

    // Open the file
    std::ofstream file("/Users/irene/TUM/thesis/AutoPas_Thesis/AutoPasFunctorBench/ArgonPotential_EquilateralGeometry.csv");

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cout << "Error opening file!" << std::endl;
        return 1;
    }

    // Write headers
    file << "Distance [nm],Potential [K],Dispersion Potential [K],Repulsive Potential [K]" << std::endl;

    const double R_min = .3;
    const double R_max = 1;
    const double step_size = 0.01;
    for (double r = R_min; r <= R_max; r += step_size) {
        // resize the triangle
        auto positionI = autopas::utils::ArrayMath::rescale_array(positionI_unit, r);
        auto positionJ = autopas::utils::ArrayMath::rescale_array(positionJ_unit, r);
        auto positionK = autopas::utils::ArrayMath::rescale_array(positionK_unit, r);

        auto displacementHandleIJ = DisplacementHandle(positionI, positionJ, I, J);
        auto displacementHandleJK = DisplacementHandle(positionJ, positionK, J, K);
        auto displacementHandleKI = DisplacementHandle(positionK, positionI, K, I);

        const auto dispersionPotential = U_dispersive(Z, beta, displacementHandleIJ, displacementHandleJK,
                                                      displacementHandleKI);
        const auto repulsivePotential = U_repulsive(A, alpha, displacementHandleIJ, displacementHandleJK,
                                                    displacementHandleKI);
        const auto totalPotential = dispersionPotential + repulsivePotential;

        /*const auto dispersiveForceI = autopas::utils::ArrayMath::Argon::F_dispersive<I>(Z, beta, displacementHandleIJ, displacementHandleJK,
                                                      displacementHandleKI);
        const auto repulsiveForceI = autopas::utils::ArrayMath::Argon::F_repulsive<I>(A, alpha, displacementHandleIJ, displacementHandleJK,
                                                    displacementHandleKI);
        const auto totalForceI = autopas::utils::ArrayMath::sum_arrays(dispersiveForceI, repulsiveForceI);
*/
        // std::cout << "Distance,   " << distance << std::endl;
        // std::cout << "Force,      " << force[0] << std::endl;
        // std::cout << "Potential Energy:" << epot << std::endl;

        // Write to CSV
        file << std::fixed << std::setprecision(15) << r << "," << totalPotential << ","
             << dispersionPotential << ',' << repulsivePotential << std::endl;
    }

    file.close();
    return 0;
}


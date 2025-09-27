// Copyright (c) 2025 yoshinarikamakura
#include <forward_list>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include "Definitions.h"
#include "FBMC.h"
using namespace std;

int main() {
    std::mt19937 sharedRng(42);

    const double TEMPERATURE = 300.0;  // (K)
    const double TIMESTEP = 1.0e-15;  // (s)
    const double THERMAL_ENERGY = BOLTZMANN * TEMPERATURE
                                  / ELEMENTARY_CHARGE;  // (eV)
    const double START_MEASUREMENT = 1.0e-12;  // (s)

    Vector3 efield;

    FBMC electron("Si_ELECTRON_PARAM.txt",
                  "Si_ECB.txt",
                  "Si_PH.txt",
                  TEMPERATURE,
                  sharedRng);

    FBMC hole("Si_HOLE_PARAM.txt",
              "Si_EVB.txt",
              "Si_PH.txt",
              TEMPERATURE,
              sharedRng);

    std::vector<double> dos_electron = electron.getDOSTable();
    std::vector<double> dos_hole = hole.getDOSTable();

    efield = Vector3(-1.0e5, 0.0, 0.0);  // (V/m)

    std::cout << "# Electric Field (V/m), "
              << "  Drift Velocity (m/sec), "
              << "  Mean Energy (eV), "
              << "  Alpha (1/m)\n";

    for (int i = 0; i < 10; ++i) {
        double time_integration = 0.0;
        double drift_length = 0.0;
        double energy_integration = 0.0;
        double number_of_impact_ionization = 0.0;
        State state =
        electron.selectStateOnIsoEnergySurface(THERMAL_ENERGY,
                                               Vector3(0.0, 0.0, 0.0),
                                               -1.0);

        for (int step = 0; step < 2000000; ++step) {
            state = electron.flightFree(TIMESTEP, efield, state);

            int scattering_mechanism;
            std::array<double, 2> eh_pair_energies;
            state = electron.scatter(TIMESTEP,
                                     state,
                                     dos_hole,
                                     scattering_mechanism,
                                     eh_pair_energies);

            if (time_integration > START_MEASUREMENT) {
                Vector3 v = electron.getVelocity(state);
                double energy = electron.getEnergy(state);
                drift_length += v.x * TIMESTEP;
                energy_integration += energy * TIMESTEP;
                if (scattering_mechanism == IMPACT_IONIZATION) {
                    number_of_impact_ionization += 1.0;
                }
            }

            time_integration += TIMESTEP;
        }

        std::cout
        << -efield.x << ' '
        << drift_length / (time_integration - START_MEASUREMENT) << ' '
        << energy_integration / (time_integration - START_MEASUREMENT) << ' '
        << number_of_impact_ionization / drift_length
        << std::endl;

        efield.x *= 2.0;
    }

    efield = Vector3(1.0e5, 0.0, 0.0);  // (V/m)
    for (int i = 0; i < 10; ++i) {
        double time_integration = 0.0;
        double drift_length = 0.0;
        double energy_integration = 0.0;
        double number_of_impact_ionization = 0.0;
        State state =
        hole.selectStateOnIsoEnergySurface(THERMAL_ENERGY,
                                           Vector3(0.0, 0.0, 0.0),
                                           1.0);

        for (int step = 0; step < 2000000; ++step) {
            state = hole.flightFree(TIMESTEP, efield, state);

            int scattering_mechanism;
            std::array<double, 2> eh_pair_energies;
            state = hole.scatter(TIMESTEP,
                                 state,
                                 dos_electron,
                                 scattering_mechanism,
                                 eh_pair_energies);

            if (time_integration > START_MEASUREMENT) {
                Vector3 v = hole.getVelocity(state);
                double energy = hole.getEnergy(state);
                drift_length += v.x * TIMESTEP;
                energy_integration += energy * TIMESTEP;
                if (scattering_mechanism == IMPACT_IONIZATION) {
                    number_of_impact_ionization += 1.0;
                }
            }

            time_integration += TIMESTEP;
        }

        std::cout
        << efield.x << ' '
        << drift_length / (time_integration - START_MEASUREMENT) << ' '
        << energy_integration / (time_integration - START_MEASUREMENT) << ' '
        << number_of_impact_ionization / drift_length << ' '
        << number_of_impact_ionization
        << std::endl;

        efield.x *= 2.0;
    }

    return 0;
}

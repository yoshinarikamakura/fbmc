// Copyright (c) 2025 yoshinarikamakura
#include <forward_list>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include "Definitions.h"
#include "FBMC.h"

// Charge of carriers (ELEMENTARY_CHARGE)
const double CHARGE_OF_ELECTRON = -1.0;
const double CHARGE_OF_HOLE = 1.0;

// Lattice temperature (K)
const double TEMPERATURE = 300.0;

// Dielectric constant of Si ()
const double DIELECTRIC_CONSTANT = 11.7;

// Domain of the device with a rectangular shape (m)
const Vector3 RMIN(0.0, 0.0, 0.0);
const Vector3 RMAX(1.0e-3, 1.0e-6, 1.0e-6);

// Doping density in the drift layer (1/m/m/m)
const double ND = 1.0e22;

// Maximum number of carriers
const int NUMBER_OF_MAX_CARRIERS = 128;

// Electric field vector in a one-sided junction (V/m)
inline Vector3 getElectricField(const Vector3 r, const double WX) {
    // Electric field at x = 0 (V/m)
    const double EMAX = -ELEMENTARY_CHARGE * ND * WX
                        / (DIELECTRIC_CONSTANT * EPSILON_0);

    if (r.x > RMIN.x && r.x < RMIN.x + WX) {
        return Vector3((1.0 - (r.x - RMIN.x) / WX) * EMAX, 0.0, 0.0);
    } else {
        return Vector3(0.0, 0.0, 0.0);
    }
}

// Save the carrier positions to files
void takeSnapShot(const std::forward_list <State> carrier_list,
                  const int label) {
    std::string FILENAME_e = "./outputs/electrons"
                             + std::to_string(label) + ".txt";
    std::string FILENAME_h = "./outputs/holes"
                             + std::to_string(label) + ".txt";
    std::ofstream fout_e(FILENAME_e);
    std::ofstream fout_h(FILENAME_h);
    if (fout_e.is_open() && fout_h.is_open()
        && label < 2000) {
        for (auto it = carrier_list.begin(); it != carrier_list.end(); ++it) {
            if ((*it).charge < 0.0) {
                // Electron
                fout_e << (*it).r.x << ' '
                       << (*it).r.y << ' '
                       << (*it).r.z << std::endl;
            } else {
                // Hole
                fout_h << (*it).r.x << ' '
                       << (*it).r.y << ' '
                       << (*it).r.z << std::endl;
            }
        }
    } else {
        std::cerr << "Cannot open file.\n";
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char*argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: ./a.out "
                  << "<depletion layer voltage> "
                  << "<initial energy> "
                  << "<initial position> "
                  << "<initial charge> "
                  << "<seed>\n";
        return 1;
    }

    double argv_VW = atof(argv[1]);
    std::cout << "# Depletion layer voltage (V) = "
              << argv_VW
              << std::endl;

    // Deplection layer width (m)
    const double WX = sqrt(2.0 * DIELECTRIC_CONSTANT * EPSILON_0 * argv_VW
                           / (ELEMENTARY_CHARGE * ND));

    std::cout << "# Depletion layer width (m) = "
              << WX
              << std::endl;

    double argv_INITIAL_CARRIER_ENERGY = atof(argv[2]);
    const double THERMAL_ENERGY = TEMPERATURE * BOLTZMANN / ELEMENTARY_CHARGE;
    if (argv_INITIAL_CARRIER_ENERGY < 1.5 * THERMAL_ENERGY) {
        argv_INITIAL_CARRIER_ENERGY = 1.5 * THERMAL_ENERGY;
    }
    std::cout << "# Initial carrier energy (eV) = "
              << argv_INITIAL_CARRIER_ENERGY
              << std::endl;

    double argv_INITIAL_CARRIER_POSITION = atof(argv[3]);
    std::cout << "# Initial carrier position (%) = "
              << argv_INITIAL_CARRIER_POSITION
              << std::endl;

    double argv_INITIAL_CARRIER_CHARGE = atof(argv[4]);
    std::cout << "# Initial carrier charge (ELEMENTARY_CHARGE) = "
              << argv_INITIAL_CARRIER_CHARGE
              << std::endl;

    int argv_seed = atoi(argv[5]);
    std::cout << "# Random number seed = "
              << argv_seed
              << std::endl;

    std::cout << "# Max. electric field (V/m) = "
              << getElectricField(Vector3(1.0e-10, 0.0, 0.0), WX).x
              << std::endl;

    // Initializing random number generator
    std::mt19937 sharedRng(argv_seed);

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

    const int NUMBER_OF_SAMPLES = 32;

    int n_broken = 0;

    // Loop for samples
    for (int sample = 0; sample < NUMBER_OF_SAMPLES; ++sample) {
        // Initial electron position (m)
        const double initial_x = (WX - RMIN.x) * argv_INITIAL_CARRIER_POSITION
                                 + RMIN.x;
        const double initial_y = (RMIN.y + RMAX.y) * 0.5;
        const double initial_z = (RMIN.z + RMAX.z) * 0.5;

        std::forward_list <State> carrier_list;

        State state;

        if (argv_INITIAL_CARRIER_CHARGE < 0.0) {
            state = electron.selectStateOnIsoEnergySurface(
                        argv_INITIAL_CARRIER_ENERGY,
                        Vector3(initial_x, initial_y, initial_z),
                        argv_INITIAL_CARRIER_CHARGE);
        } else {
            state = hole.selectStateOnIsoEnergySurface(
                        argv_INITIAL_CARRIER_ENERGY,
                        Vector3(initial_x, initial_y, initial_z),
                        argv_INITIAL_CARRIER_CHARGE);
        }

        carrier_list.push_front(state);

        // Time step (s)
        const double TIMESTEP = 1.0e-15;

        const int NUMBER_OF_STEPS = 2000000;

        double time = 0.0;
        int label_for_snapshot = 0;

        // Loop for time
        for (int step = 0; step < NUMBER_OF_STEPS; ++step) {
            if (carrier_list.empty()) {
                std::cout << "# " << sample << ": No carrier.\n";
                break;
            }

            int number_of_electrons = 0;
            int number_of_holes = 0;

            auto prev = carrier_list.before_begin();
            auto curr = carrier_list.begin();

            // Loop for carriers
            while (curr != carrier_list.end()) {
                // Carrier state from the list
                State state = *curr;

                Vector3 efield = getElectricField(state.r, WX);

                // proceed a step
                int scattering_mechanism;
                std::array<double, 2> eh_pair_energies;
                if (state.charge < 0.0) {
                    // Electron
                    ++number_of_electrons;
                    state = electron.flightFree(TIMESTEP, efield, state);
                    state = electron.scatter(TIMESTEP,
                                             state,
                                             scattering_mechanism,
                                             eh_pair_energies);
                } else {
                    // Hole
                    ++number_of_holes;
                    state = hole.flightFree(TIMESTEP, efield, state);
                    state = hole.scatter(TIMESTEP,
                                         state,
                                         scattering_mechanism,
                                         eh_pair_energies);
                }

                // At x-boundaries
                if (state.r.x < RMIN.x) {
                    if (state.charge < 0.0) {
                        state.r.x = 2.0 * RMIN.x - state.r.x;
                        state.k.x = -state.k.x;
                    } else {
                        curr = carrier_list.erase_after(prev);
                        continue;
                    }
                } else if (state.r.x > WX || state.r.x > RMAX.x) {
                    if (state.charge < 0.0) {
                        curr = carrier_list.erase_after(prev);
                        continue;
                    } else {
                        state.r.x = 2.0 * WX - state.r.x;
                        state.k.x = -state.k.x;
                    }
                }

                // Mirror reflextion at y-boundaries
                if (state.r.y < RMIN.y) {
                    state.r.y = 2.0 * RMIN.y - state.r.y;
                    state.k.y = -state.k.y;
                } else if (state.r.y > RMAX.y) {
                    state.r.y = 2.0 * RMAX.y - state.r.y;
                    state.k.y = -state.k.y;
                }

                // Mirror reflextion at z-boundaries
                if (state.r.z < RMIN.z) {
                    state.r.z = 2.0 * RMIN.z - state.r.z;
                    state.k.z = -state.k.z;
                } else if (state.r.z > RMAX.z) {
                    state.r.z = 2.0 * RMAX.z - state.r.z;
                    state.k.z = -state.k.z;
                }

                // Update the state
                *curr = state;

                // Create an e-h pair
                if (scattering_mechanism == IMPACT_IONIZATION) {
                    State new_carrier0, new_carrier1;
                    if (state.charge < 0.0) {
                        // I.I. initiated by electron
                        new_carrier0 = hole.selectStateOnIsoEnergySurface(
                                           eh_pair_energies[0],
                                           state.r,
                                           CHARGE_OF_HOLE);
                        new_carrier1 = electron.selectStateOnIsoEnergySurface(
                                           eh_pair_energies[1],
                                           state.r,
                                           CHARGE_OF_ELECTRON);
                    } else {
                        // I.I. initiated by hole
                        new_carrier0 = electron.selectStateOnIsoEnergySurface(
                                           eh_pair_energies[0],
                                           state.r,
                                           CHARGE_OF_ELECTRON);
                        new_carrier1 = hole.selectStateOnIsoEnergySurface(
                                           eh_pair_energies[1],
                                           state.r,
                                           CHARGE_OF_HOLE);
                    }

                    carrier_list.push_front(new_carrier0);
                    carrier_list.push_front(new_carrier1);
                    ++number_of_electrons;
                    ++number_of_holes;
                }

                prev = curr;
                ++curr;
            }

            if (step % 1000 == 0) {
                std::cout << time << ' '
                          << number_of_electrons << ' '
                          << number_of_holes << ' '
                          << std::endl;
                takeSnapShot(carrier_list, label_for_snapshot++);
            }

            if (number_of_electrons + number_of_holes > NUMBER_OF_MAX_CARRIERS) {
                std::cout << "# " << sample << ": Too many carriers.\n";
                n_broken++;
                break;
            }

            time += TIMESTEP;
        }
    }

    std::cout << static_cast<double>(n_broken) / NUMBER_OF_SAMPLES << std::endl;

    return 0;
}

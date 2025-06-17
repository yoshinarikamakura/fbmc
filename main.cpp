#include <forward_list>
#include <fstream>
#include <iostream>
#include <random>
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
const int NUMBER_OF_MAX_CARRIERS = 1000;

// Electric field vector in a one-sided junction (V/m)
inline Vector3 getElectricField(const Vector3 r, const double WX) {
    // Electric field at x = 0 (V/m)
    const double EMAX = ELEMENTARY_CHARGE * ND * WX / (DIELECTRIC_CONSTANT * EPSILON_0);

    if (r.x > RMIN.x && r.x < RMIN.x + WX) {
        return Vector3((1.0 - (r.x - RMIN.x) / WX) * EMAX, 0.0, 0.0);
    }
    else {
        return 0.0;
    }
}

// Save the carrier positions to files
void takeSnapShot(const std::forward_list <State> carrier_list, const int label) {
    std::string FILENAME_e = "./outputs/electrons" + std::to_string(label) + ".txt";
    std::string FILENAME_h = "./outputs/holes" + std::to_string(label) + ".txt";
    std::ofstream fout_e(FILENAME_e);
    std::ofstream fout_h(FILENAME_h);
    if (fout_e.is_open() && fout_h.is_open()) {
        for (auto it = carrier_list.begin(); it != carrier_list.end(); ++it) {
            if ((*it).charge < 0.0) {
                // Electron
                fout_e << (*it).r.x << ' ' << (*it).r.y << ' ' << (*it).r.z << std::endl;
	    }
            else {
                // Hole
                fout_h << (*it).r.x << ' ' << (*it).r.y << ' ' << (*it).r.z << std::endl;
	    }
        }
    }
    else {
        std::cerr << "Cannot open file.\n";
	exit(EXIT_FAILURE);
    }
}



int main(int argc, char*argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: ./a.out <depletion layer voltage> <initial energy> <seed>\n";
	return 1;
    }

    double argv_VW = atof(argv[1]);
    std::cout << "# Depletion layer voltage (V) = " << argv_VW << std::endl;

    // Deplection layer width (m)
    const double WX = sqrt(2.0 * DIELECTRIC_CONSTANT * EPSILON_0 * argv_VW / (ELEMENTARY_CHARGE * ND));
    std::cout << "# Depletion layer width (m) = " << WX << std::endl;

    double argv_INITIAL_ELECTRON_ENERGY = atof(argv[2]);
    const double THERMAL_ENERGY = TEMPERATURE * BOLTZMANN / ELEMENTARY_CHARGE;
    if (argv_INITIAL_ELECTRON_ENERGY < 1.5 * THERMAL_ENERGY) {
        argv_INITIAL_ELECTRON_ENERGY = 1.5 * THERMAL_ENERGY;
    }
    std::cout << "# Initial electron energy (eV) = " << argv_INITIAL_ELECTRON_ENERGY << std::endl;

    int argv_seed = atoi(argv[3]);
    std::cout << "# Random number seed = " << argv_seed << std::endl;

    // Initializing random number generator
    std::mt19937 sharedRng(argv_seed);

    FBMC electron("Si_ELECTRON_PARAM.txt", "Si_ECB.txt", "Si_PH.txt", TEMPERATURE, sharedRng);
    FBMC hole("Si_HOLE_PARAM.txt", "Si_EVB.txt", "Si_PH.txt", TEMPERATURE, sharedRng);

    // Initial electron position (m)
    const double initial_x = (WX - RMIN.x) * 0.05 + RMIN.x;
    const double initial_y = (RMIN.y + RMAX.y) * 0.5;
    const double initial_z = (RMIN.z + RMAX.z) * 0.5;
    
    State state = electron.selectStateOnIsoEnergySurface(argv_INITIAL_ELECTRON_ENERGY,
                                                         Vector3(initial_x, initial_y, initial_z),
                                                         CHARGE_OF_ELECTRON);
    std::forward_list <State> carrier_list;
    carrier_list.push_front(state);

    // Time step (s)
    const double TIMESTEP = 1.0e-15;

    double time = 0.0;
    int label_for_snapshot = 0;

    // Loop for time
    for (int step = 0; step < 100000; ++step) {
	if (carrier_list.empty()) {
            std::cout << "# No carriers.\n";
            break;
	}

	int number_of_electrons = 0;
	int number_of_holes = 0;
	auto it = carrier_list.before_begin();

        // Loop for carriers
	while (next(it) != carrier_list.end()) {

	    // Carrier state from the list
            State state = *next(it); 

            Vector3 efield = getElectricField(state.r, WX);

	    // proceed a step
            int scattering_mechanism;
            std::array<double, 2> eh_pair_energies;
	    if (state.charge < 0.0) {
                // Electron
                ++number_of_electrons;
                state = electron.flightFree(TIMESTEP, efield, state);
                state = electron.scatter(TIMESTEP, state, scattering_mechanism, eh_pair_energies);
	    }
	    else {
                // Hole
                ++number_of_holes;
                state = hole.flightFree(TIMESTEP, efield, state);
                state = hole.scatter(TIMESTEP, state, scattering_mechanism, eh_pair_energies);
	    }

	    // Annihilation at x-boundaries
	    if (state.r.x < RMIN.x) {
                carrier_list.erase_after(it);
		continue;
	    }
	    else if (state.r.x > WX || state.r.x > RMAX.x) {
                carrier_list.erase_after(it);
		continue;
	    }

	    // Mirror reflextion at y-boundaries
            if (state.r.y < RMIN.y) {
                state.r.y = 2.0 * RMIN.y - state.r.y;
                state.k.y = -state.k.y;
            }
            else if (state.r.y > RMAX.y) {
                state.r.y = 2.0 * RMAX.y - state.r.y;
                state.k.y = -state.k.y;
            }

	    // Mirror reflextion at z-boundaries
            if (state.r.z < RMIN.z) {
                state.r.z = 2.0 * RMIN.z - state.r.z;
                state.k.z = -state.k.z;
            }
            else if (state.r.z > RMAX.z) {
                state.r.z = 2.0 * RMAX.z - state.r.z;
                state.k.z = -state.k.z;
            }

	    // Create an e-h pair
	    if (scattering_mechanism == IMPACT_IONIZATION) {
                State new_carrier0, new_carrier1;
                if (state.charge < 0.0) {
                    // I.I. initiated by electron
                    new_carrier0 = hole.selectStateOnIsoEnergySurface(eh_pair_energies[0],
                                                                      state.r,
                                                                      CHARGE_OF_HOLE);
                    new_carrier1 = electron.selectStateOnIsoEnergySurface(eh_pair_energies[1],
                                                                          state.r,
                                                                          CHARGE_OF_ELECTRON);
		}
		else {
                    // I.I. initiated by hole
                    new_carrier0 = electron.selectStateOnIsoEnergySurface(eh_pair_energies[0],
                                                                          state.r,
                                                                          CHARGE_OF_ELECTRON);
                    new_carrier1 = hole.selectStateOnIsoEnergySurface(eh_pair_energies[1],
                                                                      state.r,
                                                                      CHARGE_OF_HOLE);
		}
                carrier_list.push_front(new_carrier0);
                carrier_list.push_front(new_carrier1);
            }

	    // Update the list
	    *next(it) = state;

	    ++it;
	}

        if (step % 100 == 0) {
            std::cout << time << ' '
                      << number_of_electrons << ' '
		      << number_of_holes << ' '
                      << std::endl;
            takeSnapShot(carrier_list, label_for_snapshot++);
        }

	if (number_of_electrons + number_of_holes > NUMBER_OF_MAX_CARRIERS) {
            std::cout << "# Too many carriers.\n";
            break;
	}

        time += TIMESTEP;
    }

    return 0;
}




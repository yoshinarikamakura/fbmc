#include <iostream>
#include <random>
#include "Definitions.h"
#include "FBMC.h"
using namespace std;

int main(void) {

    mt19937 sharedRng(42);

    const double TEMPERATURE = 300.0; // (K)
    const double TIMESTEP = 1.0e-15; // (s)
    const double THERMAL_ENERGY = BOLTZMANN * TEMPERATURE / ELEMENTARY_CHARGE; // (eV)

    FBMC electron("Si_PARAM.txt", "Si_ECB.txt", "Si_PH.txt", TEMPERATURE, sharedRng);

    /*
    for (int loop = 0; loop < 1000; ++loop) {
        State s = electron.selectStateOnIsoEnergySurface(0.04);
        double energy = electron.getEnergy(s);
        cout << s.k.x << ' ' << s.k.y << ' ' << s.k.z << ' ' << s.n << ' ' << energy << endl;
    }
    exit(-1);
    */
    /*
    for (double energy = 0.001; energy < 10.0; energy *= 1.01) {
        State s = electron.selectStateOnIsoEnergySurface(energy);
        cout << energy << ' ';
        for (int kind = 0; kind < 6; ++kind) {
            double abs = electron.getPhononScatteringRate(kind, true, s);
            double emi = electron.getPhononScatteringRate(kind, false, s);
            cout << abs << ' ' << emi << ' ';
        }
        cout << endl;
    }
    exit(-1);
    */


    Vector3 efield(100.0e5, 0.0, 0.0); // (V/m)
    int scattering_mechanism;
    double time = 0.0;
    State state = electron.selectStateOnIsoEnergySurface(THERMAL_ENERGY);
    for (int step = 0; step < 1000000; ++step) {
        state = electron.flightFree(TIMESTEP, efield, state);
        state = electron.scatter(TIMESTEP, state, scattering_mechanism);
	if (step % 1000 == 0) {
	    cout << time << ' ' << state.r.x << ' ' << state.r.y << ' ' << state.r.z << endl;
	}
        time += TIMESTEP;
    }
    return 0;


    /*
    Vector3 efield(1.0e5, 0.0, 0.0); // (V/m)

    const int START_MEASUREMENT = 1.0e-12; // (s)

    cout << "# Electric Field (V/m), Drift Velocity (m/sec), Mean Energy (eV), Alpha (1/m)\n";
    for (int i = 0; i < 10; ++i) {

        double time_integration = 0.0;
        double drift_length = 0.0;
        double energy_integration = 0.0;
        double number_of_impact_ionization = 0.0;

        State state = electron.selectStateOnIsoEnergySurface(THERMAL_ENERGY);
	int scattering_mechanism = SELF_SCATTERING;

        for (int step = 0; step < 100000; ++step) {

            state = electron.flightFree(TIMESTEP, efield, state);

            state = electron.scatter(TIMESTEP, state, scattering_mechanism);

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

        cout << efield.x << ' '
             << drift_length / (time_integration - START_MEASUREMENT) << ' '
	     << energy_integration / (time_integration - START_MEASUREMENT) << ' ' 
	     << number_of_impact_ionization / drift_length
	     << endl;

	efield.x *= 2.0;
    }	
    */

    return 0;
}

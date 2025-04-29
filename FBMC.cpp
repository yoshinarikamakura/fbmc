#include "FBMC.h"
#include <iostream>
#include "EnergyBand.h"
using namespace std;

FBMC::FBMC(string param_filename, string band_filename, string phonon_filename, double temperature, mt19937& rng)
: mt_(rng),
  carrier_(TABLE_DIRECTORY_NAME_ + band_filename, rng),
  phonon_(TABLE_DIRECTORY_NAME_ + phonon_filename, rng) {

    constructDOSTables(band_filename);

    constructPhononScatteringTables(param_filename, band_filename, temperature);

    cout << "# Look-up tables have been successfully created.\n";
}



State FBMC::flightFree(const double dt, const Vector3 efield, State initial_state) {
    State state = initial_state;
    state.k = state.k + efield * FACTOR_FREE_FLIGHT_ * dt;
    state.k = state.k.reduce_FCC();
    Vector3 v = carrier_.getVelocity(state);
    state.r = state.r + v * dt;
    return state;
}



State FBMC::scatter(const double dt, State initial_state, int& scattering_mechanism) {
    uniform_real_distribution<double> urand(0.0, 1.0);

    const int NETA = phonon_.getNumberOfBands();
    double rgamma = urand(mt_) / dt;
    
    bool isAbsorption;
    double sum = 0.0;
    for (int eta = 0; eta < NETA; ++eta) {

        isAbsorption = false;
        // Phonon Emission Scattering
        sum += getPhononScatteringRate(eta, isAbsorption, initial_state);
        if (sum > rgamma) {
            State final_state = selectStateAfterPhononScattering(eta, isAbsorption, initial_state);
	    scattering_mechanism = PHONON_EMISSION + eta;
            return final_state;
        }

        isAbsorption = true;
        // Phonon Absorption Scattering
        sum += getPhononScatteringRate(eta, isAbsorption, initial_state);
        if (sum > rgamma) {
            State final_state = selectStateAfterPhononScattering(eta, isAbsorption, initial_state);
	    scattering_mechanism = PHONON_ABSORPTION + eta;
            return final_state;
        }
    }

    sum += getImpactIonizationRate(initial_state);
    if (sum > rgamma) {
        State final_state = selectStateAfterImpactIonization(initial_state);
	scattering_mechanism = IMPACT_IONIZATION;
        return final_state;
    }

    // Self-scattering
    scattering_mechanism = SELF_SCATTERING; 
    return initial_state;
}



State FBMC::selectStateOnIsoEnergySurface(const double energy) {

    uniform_real_distribution<double> urand(0.0, 1.0);

    const int ie = static_cast<int>(energy / DELTAE_MINMAX_);
    const int ie_dosmax = static_cast<int>(energy / DELTAE_DOSMAX_);
    const int ie1_min = (ie > NE41_MINMAX_ - 1) ? ie - (NE41_MINMAX_ - 1) : 0;
    const int ie1_max = (ie < NE1_MINMAX_) ? ie : NE1_MINMAX_ - 1;

    vector<pair<int, int>> candidates;
    size_t total = 0;
    for (int ie1 = ie1_min; ie1 <= ie1_max; ++ie1) {
        int ie41_min = ie - ie1;
        for (int ie41 = ie41_min; ie41 < NE41_MINMAX_; ++ie41) {
            auto& elements = minmax_[ie1][ie41];
            total += elements.size();
        }
    }

    candidates.reserve(total);
    for (int ie1 = ie1_min; ie1 <= ie1_max; ++ie1) {
        int ie41_min = ie - ie1;
        for (int ie41 = ie41_min; ie41 < NE41_MINMAX_; ++ie41) {
            for (const auto& p : minmax_[ie1][ie41]) {
                candidates.emplace_back(p);
            }
        }
    }

    const int num_candidates = candidates.size();

    if (num_candidates == 0) {
        cout << "# Error in Table::selectStateOnIsoEnergySurface(),\n";
        cout << "#   ===> No candidates.\n";
        cout << "# Energy = " << energy << endl;
    }

    bool isSuccess = false;
    int it = 0, ib = 0;
    int loop = 0;
    do {
        int ir = static_cast<int>(urand(mt_) * num_candidates);
        it = candidates[ir].first;
        ib = candidates[ir].second;
        if (carrier_.getTetrahedronDOS(it, ib, energy) > dosmax_[ie_dosmax] * urand(mt_)) {
            isSuccess = true;
        }
	else {
            isSuccess = false;
	}
	loop += 1;
    } while (loop < INFINITE_LOOP_THRESHOLD_ && !isSuccess);

    if (!isSuccess) {
        cerr << "# Error in FBMC::selectStateOnIsoEnergySurface(),\n";
        cerr << "#     ===> Infinite loop.\n";
        cerr << "# Energy [eV] = " <<  energy << endl;
	exit(EXIT_FAILURE);
    }

    return carrier_.getStateInTetrahedron(energy, it, ib);
}

#include "FBMC.h"
#include <fstream>
#include <iostream>
#include <map>
#include <set>
using namespace std;

// Table of electron-phonon scattering rate (1/s)
static vector<std::vector<double>> phonon_absorption_scattering_rate_table;
static vector<std::vector<double>> phonon_emission_scattering_rate_table;

// Table of Maximum Matrix Elements of Electron-Phonon Scattering (???)
static vector<std::vector<double>> matrix_element_max;

static double BETA; 
static double FACTOR_PHONON_SCATTERING_RATE; 
static double ENERGY_GEOMETRIC_RATIO;
static vector<double> hwmax;

void FBMC::constructPhononScatteringTables(string param_filename, string band_filename, double temperature) {
    loadMaterialParameter(param_filename);

    cout << "# Material parameters have been loaded.\n";

    const double emax = carrier_.getMaximumEnergy();
    double k_unit = carrier_.getUnitOfWaveVector();
    double factor_dos = carrier_.getFactorForDensityOfStates();
    hwmax = phonon_.getMaximumEnergyPerBand();
    const int NETA = phonon_.getNumberOfBands();

    ENERGY_GEOMETRIC_RATIO = pow(emax / MINIMUM_ENERGY_FOR_SCATTERING_RATE_TABLE_, 1.0 / SIZE_OF_SCATTERING_RATE_TABLE_);
    BETA = ELEMENTARY_CHARGE / (BOLTZMANN * temperature);
    FACTOR_PHONON_SCATTERING_RATE = 0.5 * factor_dos * k_unit * k_unit * M_PI * HBAR / DENSITY_; 
    FACTOR_FREE_FLIGHT_ = ELEMENTARY_CHARGE / (HBAR * k_unit);

    matrix_element_max.resize(NETA);
    for (int eta = 0; eta < NETA; ++eta) {
        matrix_element_max[eta].resize(2);
    }
    matrix_element_max = getMaxToyMatrixElement(emax);

    phonon_absorption_scattering_rate_table.resize(SIZE_OF_SCATTERING_RATE_TABLE_);
    phonon_emission_scattering_rate_table.resize(SIZE_OF_SCATTERING_RATE_TABLE_);
    for (int isc = 0; isc < SIZE_OF_SCATTERING_RATE_TABLE_; ++isc) {
        phonon_absorption_scattering_rate_table[isc].resize(NETA);
        phonon_emission_scattering_rate_table[isc].resize(NETA);
    }
    loadPhononScatteringRate(band_filename);

    cout << "# phonon scattering rate table has been loaded.\n";
} 

void FBMC::loadMaterialParameter(std::string param_filename) {

    ifstream fin_param(TABLE_DIRECTORY_NAME_ + param_filename, ios::in);
    if (!fin_param) {
        cerr << "Cannot open: " << param_filename << endl;
        exit(EXIT_FAILURE);
    }

    map<string, double*> param_map = {
        {"DENSITY", &DENSITY_},
        {"ADAC", &ADAC_},
        {"BDAC", &BDAC_},
        {"ADKOP", &ADKOP_},
        {"BDKOP", &BDKOP_},
        {"ETHII", &ETHII_},
        {"SII", &SII_},
        {"POWII", &POWII_},
        {"AII1", &AII1_},
        {"BII1", &BII1_},
        {"AII2", &AII2_},
        {"BII2", &BII2_}
    };
    string name;
    double value;
    set<string> loaded_params;
    while (fin_param >> name >> value) {
        auto it = param_map.find(name);
        if (it != param_map.end()) {
            *(it->second) = value;
            loaded_params.insert(name);
        } else {
            cerr << "# Invalid material parameter name: " << name << endl;
            exit(EXIT_FAILURE);
        }
    }

    for (const auto& pair : param_map) {
        if (loaded_params.find(pair.first) == loaded_params.end()) {
            cerr << "# Missing material parameter: " << pair.first << endl;
            exit(EXIT_FAILURE);
        }
    }
    fin_param.close();
}

void FBMC::loadPhononScatteringRate(std::string band_filename) {
    string phonon_scattering_rate_table_filename = TABLE_DIRECTORY_NAME_ + "Phonon_scattering_rate_" + band_filename;
    ifstream fin_ephrate(phonon_scattering_rate_table_filename, ios::in);
    const int NETA = phonon_.getNumberOfBands();

    if (fin_ephrate.is_open()) {
        cout << "# File: " + phonon_scattering_rate_table_filename + " is found.\n";

        for (int isc = 0; isc < SIZE_OF_SCATTERING_RATE_TABLE_; ++isc) {
            int isc_scan;
            fin_ephrate >> isc_scan;
	    if (isc_scan == isc) {
                for (int eta = 0; eta < NETA; ++eta) {
                    double r_abs, r_emi;
                    fin_ephrate >> r_abs >> r_emi;
                    phonon_absorption_scattering_rate_table[isc][eta] = r_abs;
                    phonon_emission_scattering_rate_table[isc][eta] = r_emi;
		}
	    }
	    else {
                cerr << "ERROR in making phonon scattering rate table.\n";
                exit(EXIT_FAILURE);
	    }
        }
        fin_ephrate.close();
    }
    else {
        cout << "# Not found phonon scattering rate table.\n";
        cout << "# Now, making phonon scattering rate table. Please wait for a while.\n";

        makePhononScatteringRateTable();
        ofstream fout_ephrate(phonon_scattering_rate_table_filename, ios::out);
        for (int isc = 0; isc < SIZE_OF_SCATTERING_RATE_TABLE_; ++isc) {
            fout_ephrate << isc << endl;
            for (int eta = 0; eta < NETA; ++eta) {
                fout_ephrate << phonon_absorption_scattering_rate_table[isc][eta] << ' '
                             << phonon_emission_scattering_rate_table[isc][eta] << endl;
            }
        }
        fout_ephrate.close();
    }
}

State FBMC::selectStateAfterPhononScattering(const int eta, const bool isAbsorption, const State initial_state) {

    uniform_real_distribution<double> urand(0.0, 1.0);

    const double initial_energy = carrier_.getEnergy(initial_state); // Initial State Energy (eV)

    int ief_min = 0, ief_max = 0;
    if (isAbsorption) { // Phonon Absorption
        ief_min = static_cast<int>(initial_energy / DELTAE_MINMAX_);
        ief_max = static_cast<int>((initial_energy + hwmax[eta]) / DELTAE_MINMAX_);
    }
    else { // Phonon Emission
        ief_min = static_cast<int>((initial_energy - hwmax[eta]) / DELTAE_MINMAX_);
        ief_max = static_cast<int>(initial_energy / DELTAE_MINMAX_);
    }

    int ie1_min = ief_min - (NE41_MINMAX_ - 1);
    if (ie1_min < 0) ie1_min = 0;

    int ie1_max = ief_max;
    if (ie1_max > NE1_MINMAX_ - 1) ie1_max = NE1_MINMAX_ - 1;

    vector<pair<int, int>> candidates;

    size_t total = 0;
    for (int ie1 = ie1_min; ie1 <= ie1_max; ++ie1) {
        int ie41_min = std::max(0, ief_min - ie1);
        for (int ie41 = ie41_min; ie41 < NE41_MINMAX_; ++ie41) {
            auto& elements = minmax_[ie1][ie41];
            total += elements.size();
        }
    }

    candidates.reserve(total);

    for (int ie1 = ie1_min; ie1 <= ie1_max; ++ie1) {
        int ie41_min = max(0, ief_min - ie1);
        for (int ie41 = ie41_min; ie41 < NE41_MINMAX_; ++ie41) {
            for (const auto& p : minmax_[ie1][ie41]) {
                candidates.emplace_back(p);
	    }
        }
    }

    const int num_candidates = candidates.size();

    if (num_candidates == 0) {
        cerr << "# Error in selectStateAfterPhononScattering(),\n";
        cerr << "#   ===> No candidates\n";
        cerr << "# eta = " << eta << endl;
        exit(EXIT_FAILURE);
    }

    int ief_dosmax_min = 0, ief_dosmax_max = 0;
    if (isAbsorption) { // Phonon absorption
        ief_dosmax_min = static_cast<int>(initial_energy / DELTAE_DOSMAX_);
        ief_dosmax_max = static_cast<int>((initial_energy + hwmax[eta]) / DELTAE_DOSMAX_);
    }
    else { // Phonon emission
        ief_dosmax_min = static_cast<int>((initial_energy - hwmax[eta]) / DELTAE_DOSMAX_);
        if (ief_dosmax_min < 0) {
            ief_dosmax_min = 0;
        }
        ief_dosmax_max = static_cast<int>(initial_energy / DELTAE_DOSMAX_);
    }

    double dosmax_max = 0.0;
    for (int ief = ief_dosmax_min; ief <= ief_dosmax_max; ++ief) {
        if (dosmax_[ief] > dosmax_max) {
            dosmax_max = dosmax_[ief];
        }
    }
    double me_max = isAbsorption ? matrix_element_max[eta][1] : matrix_element_max[eta][0];

    double final_energy = 0.0;

    // Loop for final states
    int loop = 0; 
    bool isSuccess = false;
    int it = 0, ib = 0;
    do {
        const int ir = static_cast<int>(urand(mt_) * num_candidates);
        it = candidates[ir].first;
        ib = candidates[ir].second;

	// Change of wave vector before and after e-ph scattering
	Vector3 vq = carrier_.getWaveVectorDifference(initial_state, it);

	// Final state energy
	double sign = isAbsorption ? 1.0 : -1.0;	
        final_energy = initial_energy + sign * phonon_.getEnergy(State{vq, eta});

	// DOS in the final tetrahedron candidate
	double dos = carrier_.getTetrahedronDOS(it, ib, final_energy);

	// Matrix element of e-ph scattering
        double me = getToyMatrixElementForSilicon(State{vq, eta}, isAbsorption, initial_energy);

        isSuccess = (dos * me > dosmax_max * me_max * urand(mt_)) ? true : false;
	loop += 1;
    } while (loop < INFINITE_LOOP_THRESHOLD_ && !isSuccess);

    if (!isSuccess) {
        cerr << "# Warning in selectStateAfterPhononScattering(),\n";
        cerr << "#   ===> Infinite loop\n";
        cerr << "# Initial energy [eV] = " << initial_energy << endl;
        cerr << "# eta = " << eta << endl;
        cerr << "# isAbsorption = " << isAbsorption << endl;
        cerr << "# ki = " << initial_state.k.x << ' ' << initial_state.k.y << ' ' << initial_state.k.z << endl;
//	exit(EXIT_FAILURE);
        return initial_state;
    }

    State final_state = carrier_.getStateInTetrahedron(final_energy, it, ib);
    final_state.r = initial_state.r;

    return final_state;
}

double FBMC::calculatePhononScatteringRate(const int eta, const bool isAbsorption, const State initial_state) {
    const double initial_energy = carrier_.getEnergy(initial_state); // Initial State Energy (eV)
								    //
    int ief_min = 0, ief_max = 0;
    if (isAbsorption == true) { // Phonon absorption
        ief_min = static_cast<int>(initial_energy / DELTAE_MINMAX_);
        ief_max = static_cast<int>((initial_energy + hwmax[eta]) / DELTAE_MINMAX_);
    }
    else { // Phonon emission
        ief_min = static_cast<int>((initial_energy - hwmax[eta]) / DELTAE_MINMAX_);
        ief_max = static_cast<int>(initial_energy / DELTAE_MINMAX_);
    }

    int ie1_min = ief_min - (NE41_MINMAX_ - 1);
    if (ie1_min < 0) ie1_min = 0;

    int ie1_max = ief_max;
    if (ie1_max > NE1_MINMAX_ - 1) ie1_max = NE1_MINMAX_ - 1;

    vector<pair<int, int>> candidates;

    for (int ie1 = ie1_min; ie1 < ie1_max+1; ++ie1) {
        const int ie41_min = (ief_min - ie1 > 0) ? ief_min - ie1 : 0;
        for (int ie41 = ie41_min; ie41 < NE41_MINMAX_; ++ie41) {
            candidates.insert(candidates.end(), minmax_[ie1][ie41].begin(), minmax_[ie1][ie41].end());
        }
    }

    const int num_candidates = candidates.size();

    if (num_candidates == 0) {
        cerr << "# Error in selectStateAfterPhononScattering(),\n";
        cerr << "#   ===> No candidates.\n";
        cerr << "# eta = " << eta << endl;
        exit(EXIT_FAILURE);
    }

    // Loop for final states
    double sum = 0.0;
    for (int i = 0; i < num_candidates; ++i) {
        int it = candidates[i].first;
        int ib = candidates[i].second;

	// Change of wave vector before and after e-ph scattering
        Vector3 vq = carrier_.getWaveVectorDifference(initial_state, it);	

	// Final state energy
	double sign = isAbsorption ? 1.0 : -1.0;	
        double final_energy = initial_energy + sign * phonon_.getEnergy(State{vq, eta});

	// DOS in the final tetrahedron candidate
	double dos = carrier_.getTetrahedronDOS(it, ib, final_energy);

	// Matrix element of e-ph scattering
        double me = getToyMatrixElementForSilicon(State{vq, eta}, isAbsorption, initial_energy);

	sum += me * dos;
    }

    return sum * FACTOR_PHONON_SCATTERING_RATE;
}

void FBMC::makePhononScatteringRateTable(void) {
    const int NETA = phonon_.getNumberOfBands();

    for (int isc = 0; isc < SIZE_OF_SCATTERING_RATE_TABLE_; ++isc) {
        double initial_energy = MINIMUM_ENERGY_FOR_SCATTERING_RATE_TABLE_ * pow(ENERGY_GEOMETRIC_RATIO, isc);

        // Randomly select an initial state
        State initial_state = selectStateOnIsoEnergySurface(initial_energy);
        for (int eta = 0; eta < NETA; ++eta) {
            phonon_absorption_scattering_rate_table[isc][eta] = calculatePhononScatteringRate(eta, true, initial_state);
            phonon_emission_scattering_rate_table[isc][eta] = calculatePhononScatteringRate(eta, false, initial_state);
	}
    }
}

double FBMC::getPhononScatteringRate(const int eta, const bool isAbsorption, const State initial_state) {
    const double initial_energy = carrier_.getEnergy(initial_state);
    const int isc = static_cast<int>(log(initial_energy / MINIMUM_ENERGY_FOR_SCATTERING_RATE_TABLE_) / log(ENERGY_GEOMETRIC_RATIO));

    double rate = 0.0;
    if (isc < 0) {
        if (isAbsorption) {
            rate = phonon_absorption_scattering_rate_table[0][eta];
        }
        else {
            rate = phonon_emission_scattering_rate_table[0][eta];
	}
    }
    else if (isc < SIZE_OF_SCATTERING_RATE_TABLE_ - 1) {
	double e0 = MINIMUM_ENERGY_FOR_SCATTERING_RATE_TABLE_ * pow(ENERGY_GEOMETRIC_RATIO, isc);
	double e1 = e0 * ENERGY_GEOMETRIC_RATIO;
	double y0 = 0.0, y1 = 0.0;
        if (isAbsorption) {
            y0 = phonon_absorption_scattering_rate_table[isc][eta];
            y1 = phonon_absorption_scattering_rate_table[isc + 1][eta];
	}
        else {
            y0 = phonon_emission_scattering_rate_table[isc][eta];
            y1 = phonon_emission_scattering_rate_table[isc + 1][eta];
	}
        rate = y0 + (y1 - y0) * (initial_energy - e0) / (e1 - e0);
    }
    else {
         cerr << "# Error in getPhononScatteringRate(),\n";
         cerr << "#   ===> Too large initial energy.\n";
         cerr << "# Energy (eV) = " << initial_energy;
         exit(EXIT_FAILURE);
    }

    return rate;
}

vector<vector<double>> FBMC::getMaxToyMatrixElement(const double eimax) {
    constexpr int NETA = 6; // for getToyMatrixElementForSilicon()
    vector<vector<double>> memax;

    memax.resize(NETA);
    for (int eta = 0; eta < NETA; ++eta) {
        memax[eta].resize(2);
    }

    Vector3 vqmax = phonon_.getMaximumWaveVector();
    for (int eta = 0; eta < NETA; ++eta) {
        memax[eta][0] = getToyMatrixElementForSilicon(State{vqmax, eta}, false, eimax);
        memax[eta][1] = getToyMatrixElementForSilicon(State{vqmax, eta}, true, eimax);
    }

    return memax;
}

double FBMC::getToyMatrixElementForSilicon(State phonon_state, const bool isAbsorption, const double carrier_energy) {
    const double hw = phonon_.getEnergy(phonon_state);
    double matrix_element = 0.0;

    if (hw > 1.0e-4) {
        double bose = 1.0 / (exp(hw * BETA) - 1.0);
        double q2 = phonon_state.k.squared();
	double k_unit = carrier_.getUnitOfWaveVector();

        double dp2;
        switch (phonon_state.n) {
            case 0: // LO
                dp2 = (ADKOP_ * carrier_energy + BDKOP_) / (k_unit * k_unit);
                break;
            case 1: // TO
                dp2 = (ADKOP_ * carrier_energy + BDKOP_) / (k_unit * k_unit);
                break;
            case 2: // TO
                dp2 = (ADKOP_ * carrier_energy + BDKOP_) / (k_unit * k_unit);
                break;
            case 3: // LA
                dp2 = (ADAC_ * carrier_energy + BDAC_) * q2;
                break;
            case 4: // TA
                dp2 = (ADAC_ * carrier_energy + BDAC_) * q2;
                break;
            case 5: // TA
                dp2 = (ADAC_ * carrier_energy + BDAC_) * q2;
                break;
            default:
                cerr << "# Error in Phonon::getToyMatrixElement()\n";
                exit(EXIT_FAILURE);
                break;
        }

	if (isAbsorption) {
            matrix_element = dp2 * bose / hw;
	}
	else {
            matrix_element = dp2 * (bose + 1.0) / hw;
	}
    }

    return matrix_element;
}

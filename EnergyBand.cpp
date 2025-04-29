#include "EnergyBand.h"
#include <algorithm>
#include <iostream>
#include <fstream>
using namespace std;

EnergyBand::EnergyBand(string band_filename, mt19937& rng) : mt(rng) {

// Load energy band structure from a file
    loadBandFile(band_filename);

    K_UNIT = 2.0 * M_PI / LATTICE_CONSTANT;

// Number of grid points
    NG = NK * NK * NK;

// Number of tetrahedra
    NT = NK * NK * NK * 6;

    FACTOR_DOS = (KMAX.x - KMIN.x) * (KMAX.y - KMIN.y) * (KMAX.z - KMIN.z) * K_UNIT * K_UNIT * K_UNIT
                 / (8.0 * M_PI * M_PI * M_PI * NK * NK * NK);

    FACTOR_GROUP_VELOCITY = Vector3(1.0, 1.0, 1.0) * ELEMENTARY_CHARGE * NK / ((KMAX - KMIN) * K_UNIT * HBAR);

// Grid
    NG = (NK + 1) * (NK + 1) * (NK + 1);
    grid_x.resize(NG);
    grid_y.resize(NG);
    grid_z.resize(NG);
    grid_number.resize(NK+1);
    for (int ix = 0; ix < NK+1; ++ix) {
        grid_number[ix].resize(NK+1);
        for (int iy = 0; iy < NK+1; ++iy) {
            grid_number[ix][iy].resize(NK+1);
        }
    }

    int ig = 0;
    for (int ix = 0; ix < NK+1; ++ix) {
    for (int iy = 0; iy < NK+1; ++iy) {
    for (int iz = 0; iz < NK+1; ++iz) {
        grid_x[ig] = ix;
        grid_y[ig] = iy;
        grid_z[ig] = iz;
        grid_number[ix][iy][iz] = ig;
	ig += 1;
    }
    }
    }

// Tetrahedra
    tetrahedron_vertex1.resize(NT);
    tetrahedron_vertex2.resize(NT);
    tetrahedron_vertex3.resize(NT);
    tetrahedron_vertex4.resize(NT);
    cube_x.resize(NT);
    cube_y.resize(NT);
    cube_z.resize(NT);

    int it = 0;
    for (int ix = 0; ix < NK; ++ix) {
    for (int iy = 0; iy < NK; ++iy) {
    for (int iz = 0; iz < NK; ++iz) {
        const int v1 = grid_number[ix  ][iy  ][iz  ]; 
        const int v2 = grid_number[ix+1][iy  ][iz  ]; 
        const int v3 = grid_number[ix  ][iy+1][iz  ]; 
        const int v4 = grid_number[ix+1][iy+1][iz  ]; 
        const int v5 = grid_number[ix  ][iy  ][iz+1]; 
        const int v6 = grid_number[ix+1][iy  ][iz+1]; 
        const int v7 = grid_number[ix  ][iy+1][iz+1]; 
        const int v8 = grid_number[ix+1][iy+1][iz+1]; 

        tetrahedron_vertex1[it] = v1;
        tetrahedron_vertex2[it] = v2;
        tetrahedron_vertex3[it] = v3;
        tetrahedron_vertex4[it] = v5;
        cube_x[it] = ix;
        cube_y[it] = iy;
        cube_z[it] = iz;
        it += 1;

        tetrahedron_vertex1[it] = v2;
        tetrahedron_vertex2[it] = v5;
        tetrahedron_vertex3[it] = v6;
        tetrahedron_vertex4[it] = v7;
        cube_x[it] = ix;
        cube_y[it] = iy;
        cube_z[it] = iz;
        it += 1;

        tetrahedron_vertex1[it] = v2;
        tetrahedron_vertex2[it] = v3;
        tetrahedron_vertex3[it] = v5;
        tetrahedron_vertex4[it] = v7;
        cube_x[it] = ix;
        cube_y[it] = iy;
        cube_z[it] = iz;
        it += 1;

        tetrahedron_vertex1[it] = v2;
        tetrahedron_vertex2[it] = v4;
        tetrahedron_vertex3[it] = v6;
        tetrahedron_vertex4[it] = v7;
        cube_x[it] = ix;
        cube_y[it] = iy;
        cube_z[it] = iz;
        it += 1;

        tetrahedron_vertex1[it] = v2;
        tetrahedron_vertex2[it] = v3;
        tetrahedron_vertex3[it] = v4;
        tetrahedron_vertex4[it] = v7;
        cube_x[it] = ix;
        cube_y[it] = iy;
        cube_z[it] = iz;
        it += 1;

        tetrahedron_vertex1[it] = v4;
        tetrahedron_vertex2[it] = v6;
        tetrahedron_vertex3[it] = v7;
        tetrahedron_vertex4[it] = v8;
        cube_x[it] = ix;
        cube_y[it] = iy;
        cube_z[it] = iz;
        it += 1;
    }
    }
    }

    cout << "# Look-up tables have been successfully created.\n";
}



double EnergyBand::loadBandFile(std::string band_filename) {

    ifstream fin_band(band_filename, ios::in);
    if (!fin_band) {
        cerr << "# Cannot open: " << band_filename << endl;
        exit(EXIT_FAILURE);
    }

// Read division in k-space
    fin_band >> LATTICE_CONSTANT;
    cout << "# Lattice Constant (m) = " << LATTICE_CONSTANT << endl;

// Read division in k-space
    fin_band >> NK;
    cout << "# Number of k-space divisions = " << NK << endl;

// Read number of bands
    fin_band >> NB;
    cout << "# Number of bands = " << NB << endl;

    fin_band >> KMIN.x >> KMIN.y >> KMIN.z;
    cout << "# Domain of k-space; KMIN = " << KMIN.x << ' ' << KMIN.y << ' ' << KMIN.z << endl;

    fin_band >> KMAX.x >> KMAX.y >> KMAX.z;
    cout << "# Domain of k-space; KMAX = " << KMAX.x << ' ' << KMAX.y << ' ' << KMAX.z << endl;

    grid_energy.resize(NK + 1);
    for (int ix = 0; ix < NK + 1; ++ix) {
        grid_energy[ix].resize(NK + 1);
        for (int iy = 0; iy < NK + 1; ++iy) {
            grid_energy[ix][iy].resize(NK + 1);
            for (int iz = 0; iz < NK + 1; ++iz) {
                grid_energy[ix][iy][iz].resize(NB);
            }
        }
    }

    const int IKMIN_X = KMIN.x / (KMAX.x - KMIN.x) * NK;
    const int IKMIN_Y = KMIN.y / (KMAX.y - KMIN.y) * NK;
    const int IKMIN_Z = KMIN.z / (KMAX.z - KMIN.z) * NK;

    emax = -1.0e10;
    emax_per_band.resize(NB, -1.0e10);
    double emin = 1.0e10;
    double k2max = 0.0;
    for (int ix = 0; ix < NK + 1; ++ix) {
    for (int iy = 0; iy < NK + 1; ++iy) {
    for (int iz = 0; iz < NK + 1; ++iz) {
        int ikx, iky, ikz;
        fin_band >> ikx >> iky >> ikz;

        Vector3 vk(
            static_cast<double>(ix) / NK,
            static_cast<double>(iy) / NK,
            static_cast<double>(iz) / NK
        );
        vk = vk.reduce_FCC();
        double k2 = vk.squared();
        if (k2 > k2max) {
            k2max = k2;
            vkmax = vk;
        }

        for (int n = 0; n < NB; ++n) {
            double e;
            fin_band >> e;
            grid_energy[ikx - IKMIN_X][iky - IKMIN_Y][ikz - IKMIN_Z][n] = e;
	    if (e > emax) emax = e;
            if (e > emax_per_band[n]) emax_per_band[n] = e;
	    if (e < emin) emin = e;
        }
    }
    }
    }
    cout << "# Maximum energy (eV) = " << emax << endl;
    cout << "# Minimum energy (eV) = " << emin << endl;

    fin_band.close();

    return emax;
}



double EnergyBand::getEnergy(const State s) {

    int ix = static_cast<int>((s.k.x - KMIN.x) * NK / (KMAX.x - KMIN.x));
    int iy = static_cast<int>((s.k.y - KMIN.y) * NK / (KMAX.y - KMIN.y));
    int iz = static_cast<int>((s.k.z - KMIN.z) * NK / (KMAX.z - KMIN.z));

    if (ix == NK) ix = NK - 1;
    if (iy == NK) iy = NK - 1;
    if (iz == NK) iz = NK - 1;

    if (ix < 0 || iy < 0 || iz < 0 || ix > NK - 1 || iy > NK - 1 || iz > NK - 1) {
        cerr << "# Error in getEnergy():\n";
        cerr << "# kx, ky, kz = " << s.k.x << ", " << s.k.y << ", " << s.k.z << endl;
        cerr << "# ix, iy, iz = " << ix << ", " << iy << ", " << iz << endl;
        exit(EXIT_FAILURE);
    }

    const double x = (s.k.x - KMIN.x) * NK / (KMAX.x - KMIN.x) - ix;
    const double y = (s.k.y - KMIN.y) * NK / (KMAX.y - KMIN.y) - iy;
    const double z = (s.k.z - KMIN.z) * NK / (KMAX.z - KMIN.z) - iz;

    double a, b, c, d;

    if (x + y < 1.0) {
        if (x + z > 1.0) {
            // type 2
            const double e100 = grid_energy[ix+1][iy  ][iz  ][s.n];
            const double e001 = grid_energy[ix  ][iy  ][iz+1][s.n];
            const double e101 = grid_energy[ix+1][iy  ][iz+1][s.n];
            const double e011 = grid_energy[ix  ][iy+1][iz+1][s.n];
            a = e101 - e001;
            b = e011 - e001;
            c = e101 - e100;
            d = e100 + e001 - e101;
        }
        else {
            if (x + y + z < 1.0) {
                // type1
                const double e000 = grid_energy[ix  ][iy  ][iz  ][s.n];
                const double e100 = grid_energy[ix+1][iy  ][iz  ][s.n];
                const double e010 = grid_energy[ix  ][iy+1][iz  ][s.n];
                const double e001 = grid_energy[ix  ][iy  ][iz+1][s.n];
                a = e100 - e000;
                b = e010 - e000;
                c = e001 - e000;
                d = e000;
            }
            else {
                // type3
                const double e100 = grid_energy[ix+1][iy  ][iz  ][s.n];
                const double e010 = grid_energy[ix  ][iy+1][iz  ][s.n];
                const double e001 = grid_energy[ix  ][iy  ][iz+1][s.n];
                const double e011 = grid_energy[ix  ][iy+1][iz+1][s.n];
                a = e100 + e011 - e010 - e001;
                b = e011 - e001;
                c = e011 - e010;
                d = e010 + e001 - e011;
            }
        }
    }
        else {
        if (x + z < 1.0) {
            // type5
            const double e100 = grid_energy[ix+1][iy  ][iz  ][s.n];
            const double e010 = grid_energy[ix  ][iy+1][iz  ][s.n];
            const double e110 = grid_energy[ix+1][iy+1][iz  ][s.n];
            const double e011 = grid_energy[ix  ][iy+1][iz+1][s.n];
            a = e110 - e010;
            b = e110 - e100;
            c = e011 - e010;
            d = e100 + e010 - e110;
        }
        else {
            if (x + y + z < 2.0) {
                // type4
                const double e011 = grid_energy[ix  ][iy+1][iz+1][s.n];
                const double e100 = grid_energy[ix+1][iy  ][iz  ][s.n];
                const double e101 = grid_energy[ix+1][iy  ][iz+1][s.n];
                const double e110 = grid_energy[ix+1][iy+1][iz  ][s.n];
                a = e110 + e101 - e100 - e011;
                b = e110 - e100;
                c = e101 - e100;
                d = 2.0 * e100 + e011 - e110 - e101;
            }
            else {
                // type6
                const double e110 = grid_energy[ix+1][iy+1][iz  ][s.n];
                const double e101 = grid_energy[ix+1][iy  ][iz+1][s.n];
                const double e011 = grid_energy[ix  ][iy+1][iz+1][s.n];
                const double e111 = grid_energy[ix+1][iy+1][iz+1][s.n];
                a = e111 - e011;
                b = e111 - e101;
                c = e111 - e110;
                d = e110 + e101 + e011 - 2.0 * e111;
            }
        }
    }

    return a * x + b * y + c * z + d;
}



Vector3 EnergyBand::getVelocity(const State s) {

    int ix = static_cast<int>((s.k.x - KMIN.x) * NK / (KMAX.x - KMIN.x));
    int iy = static_cast<int>((s.k.y - KMIN.y) * NK / (KMAX.y - KMIN.y));
    int iz = static_cast<int>((s.k.z - KMIN.z) * NK / (KMAX.z - KMIN.z));

    if (ix == NK) ix = NK - 1;
    if (iy == NK) iy = NK - 1;
    if (iz == NK) iz = NK - 1;

    if (ix < 0 || iy < 0 || iz < 0 || ix > NK - 1 || iy > NK - 1 || iz > NK - 1) {
        cerr << "# Error in getVelocity():\n";
        cerr << "# kx, ky, kz = " << s.k.x << ", " << s.k.y << ", " << s.k.z << endl;
        cerr << "# ix, iy, iz = " << ix << ", " << iy << ", " << iz << endl;
        exit(EXIT_FAILURE);
    }

    const double x = (s.k.x - KMIN.x) * NK / (KMAX.x - KMIN.x) - ix;
    const double y = (s.k.y - KMIN.y) * NK / (KMAX.y - KMIN.y) - iy;
    const double z = (s.k.z - KMIN.z) * NK / (KMAX.z - KMIN.z) - iz;

    Vector3 v;

    if (x + y < 1.0) {
        if (x + z > 1.0) {
            // type 2
            const double e100 = grid_energy[ix+1][iy  ][iz  ][s.n];
            const double e001 = grid_energy[ix  ][iy  ][iz+1][s.n];
            const double e101 = grid_energy[ix+1][iy  ][iz+1][s.n];
            const double e011 = grid_energy[ix  ][iy+1][iz+1][s.n];
            v.x = e101 - e001;
            v.y = e011 - e001;
            v.z = e101 - e100;
        }
        else {
            if (x + y + z < 1.0) {
                // type1
                const double e000 = grid_energy[ix  ][iy  ][iz  ][s.n];
                const double e100 = grid_energy[ix+1][iy  ][iz  ][s.n];
                const double e010 = grid_energy[ix  ][iy+1][iz  ][s.n];
                const double e001 = grid_energy[ix  ][iy  ][iz+1][s.n];
                v.x = e100 - e000;
                v.y = e010 - e000;
                v.z = e001 - e000;
            }
            else {
                // type3
                const double e100 = grid_energy[ix+1][iy  ][iz  ][s.n];
                const double e010 = grid_energy[ix  ][iy+1][iz  ][s.n];
                const double e001 = grid_energy[ix  ][iy  ][iz+1][s.n];
                const double e011 = grid_energy[ix  ][iy+1][iz+1][s.n];
                v.x = e100 + e011 - e010 - e001;
                v.y = e011 - e001;
                v.z = e011 - e010;
            }
        }
    }
    else {
        if (x + z < 1.0) {
            // type5
            const double e100 = grid_energy[ix+1][iy  ][iz  ][s.n];
            const double e010 = grid_energy[ix  ][iy+1][iz  ][s.n];
            const double e110 = grid_energy[ix+1][iy+1][iz  ][s.n];
            const double e011 = grid_energy[ix  ][iy+1][iz+1][s.n];
            v.x = e110 - e010;
            v.y = e110 - e100;
            v.z = e011 - e010;
        }
        else {
            if (x + y + z < 2.0) {
                // type4
                const double e011 = grid_energy[ix  ][iy+1][iz+1][s.n];
                const double e100 = grid_energy[ix+1][iy  ][iz  ][s.n];
                const double e101 = grid_energy[ix+1][iy  ][iz+1][s.n];
                const double e110 = grid_energy[ix+1][iy+1][iz  ][s.n];
                v.x = e110 + e101 - e100 - e011;
                v.y = e110 - e100;
                v.z = e101 - e100;
            }
            else {
                // type6
                const double e110 = grid_energy[ix+1][iy+1][iz  ][s.n];
                const double e101 = grid_energy[ix+1][iy  ][iz+1][s.n];
                const double e011 = grid_energy[ix  ][iy+1][iz+1][s.n];
                const double e111 = grid_energy[ix+1][iy+1][iz+1][s.n];
                v.x = e111 - e011;
                v.y = e111 - e101;
                v.z = e111 - e110;
            }
        }
    }

    return v * FACTOR_GROUP_VELOCITY;
}



double EnergyBand::getDOS(const double e) {

    double sum = 0.0;
    for (int nt = 0; nt < NT; ++nt) {
        for (int n = 0; n < NB; ++n) {
            sum += getTetrahedronDOS(nt, n, e);
	}
    }

    return sum * FACTOR_DOS;
}



double EnergyBand::getTetrahedronDOS(const int nt, const int n, const double e) {

    const int v1 = tetrahedron_vertex1[nt];
    const int v2 = tetrahedron_vertex2[nt];
    const int v3 = tetrahedron_vertex3[nt];
    const int v4 = tetrahedron_vertex4[nt];
    double e1 = grid_energy[grid_x[v1]][grid_y[v1]][grid_z[v1]][n];
    double e2 = grid_energy[grid_x[v2]][grid_y[v2]][grid_z[v2]][n];
    double e3 = grid_energy[grid_x[v3]][grid_y[v3]][grid_z[v3]][n];
    double e4 = grid_energy[grid_x[v4]][grid_y[v4]][grid_z[v4]][n];

    array<double, 4> arr = { e1, e2, e3, e4 };

    // Sort: e1 < e2 < e3 < e4
    sort(arr.begin(), arr.end());

    e1 = arr[0];
    e2 = arr[1];
    e3 = arr[2];
    e4 = arr[3];

    double dos;
    if (e < e1) {
        dos = 0.0;
    }
    else if (e < e2) {
        dos = 0.5 * (e - e1) * (e - e1) / ((e2 - e1) * (e3 - e1) * (e4 - e1));
    }
    else if (e < e3) {
        const double a = e - e1;
        const double a2 = e2 - e1;
        const double a3 = e3 - e1;
        const double a4 = e4 - e1;
        dos = 0.5 * ((a2 - a3 - a4) * a * a + 2.0 * a3 * a4 * a - a2 * a3 * a4)
              / (a3 * a4 * (a3 - a2) * (a4 - a2));
    }
    else if (e < e4) {
        dos = 0.5 * (e4 - e) * (e4 - e) / ((e4 - e1) * (e4 - e2) * (e4 - e3));
    }
    else {
        dos = 0.0;
    }

    return dos;
}



array<double, 4> EnergyBand::getTetrahedronVertexEnergies(int nt, int nb) {

    const int v1 = tetrahedron_vertex1[nt];
    const int v2 = tetrahedron_vertex2[nt];
    const int v3 = tetrahedron_vertex3[nt];
    const int v4 = tetrahedron_vertex4[nt];

    array<double, 4> e;

    e[0] = grid_energy[grid_x[v1]][grid_y[v1]][grid_z[v1]][nb];
    e[1] = grid_energy[grid_x[v2]][grid_y[v2]][grid_z[v2]][nb];
    e[2] = grid_energy[grid_x[v3]][grid_y[v3]][grid_z[v3]][nb];
    e[3] = grid_energy[grid_x[v4]][grid_y[v4]][grid_z[v4]][nb];

    sort(e.begin(), e.end());

    return e;
}



State EnergyBand::getStateInTetrahedron(const double energy, const int nt, const int nb) {

    uniform_real_distribution<double> urand(0.0, 1.0);

    int v1 = tetrahedron_vertex1[nt];
    int v2 = tetrahedron_vertex2[nt];
    int v3 = tetrahedron_vertex3[nt];
    int v4 = tetrahedron_vertex4[nt];

    double e1 = grid_energy[grid_x[v1]][grid_y[v1]][grid_z[v1]][nb];
    double e2 = grid_energy[grid_x[v2]][grid_y[v2]][grid_z[v2]][nb];
    double e3 = grid_energy[grid_x[v3]][grid_y[v3]][grid_z[v3]][nb];
    double e4 = grid_energy[grid_x[v4]][grid_y[v4]][grid_z[v4]][nb];

    array<pair<double, int>, 4> pairs = {
        make_pair(e1, v1),
        make_pair(e2, v2),
        make_pair(e3, v3),
        make_pair(e4, v4)
    };

    // Sort by e
    sort(pairs.begin(), pairs.end());

    e1 = pairs[0].first; v1 = pairs[0].second;
    e2 = pairs[1].first; v2 = pairs[1].second;
    e3 = pairs[2].first; v3 = pairs[2].second;
    e4 = pairs[3].first; v4 = pairs[3].second;

    const Vector3 k21(
        grid_x[v2] - grid_x[v1],
        grid_y[v2] - grid_y[v1],
        grid_z[v2] - grid_z[v1]
    );

    const Vector3 k31(
        grid_x[v3] - grid_x[v1],
        grid_y[v3] - grid_y[v1],
        grid_z[v3] - grid_z[v1]
    );

    const Vector3 k41(
        grid_x[v4] - grid_x[v1],
        grid_y[v4] - grid_y[v1],
        grid_z[v4] - grid_z[v1]
    );

    Vector3 k;

    if (energy < e1) {
        cerr << "# Error in EnergyBand::getStateInTetrahedron()\n";
        cerr << "#   ===> Selected triangle does not contain given energy (energy < e1).\n";
        exit(EXIT_FAILURE);
    }
    else if (energy < e2) {
        const double phi21 = (energy - e1) / (e2 - e1);
        const double phi31 = (energy - e1) / (e3 - e1);
        const double phi41 = (energy - e1) / (e4 - e1);
        const Vector3 a = k21 * phi21;
        const Vector3 b = k31 * phi31;
        const Vector3 c = k41 * phi41;
        double r1 = urand(mt);
        double r2 = urand(mt);
        if (r1 + r2 > 1.0) {
            r1 = 1.0 - r1;
            r2 = 1.0 - r2;
        }
        k = a * (1.0 - r1 - r2) + b * r1 + c * r2;
    }
    else if (energy < e3) {
        const double phi42 = (energy - e2) / (e4 - e2);
        const double phi32 = (energy - e2) / (e3 - e2);
        const double phi31 = (energy - e1) / (e3 - e1);
        const double phi41 = (energy - e1) / (e4 - e1);

        const Vector3 a = k21 * (1.0 - phi42) + k41 * phi42;
        const Vector3 b = k21 * (1.0 - phi32) + k31 * phi32;
        const Vector3 c = k31 * phi31;
        const Vector3 d = k41 * phi41;

        const Vector3 ad = d - a;
        const Vector3 ab = b - a;
        const Vector3 cd = d - c;
        const Vector3 cb = b - c;

        const double abd2 = ad.squared() * ab.squared() - ad.squared() * ad.squared();
        const double bcd2 = cd.squared() * cb.squared() - cd.squared() * cd.squared();
        const double r0 = urand(mt);
        double r1 = urand(mt);
        double r2 = urand(mt);
        if (r1 + r2 > 1.0) {
            r1 = 1.0 - r1;
            r2 = 1.0 - r2;
        }
        const double pr_abd2 = sqrt(abd2) / (sqrt(abd2) + sqrt(bcd2));
        if (r0 < pr_abd2) {
            k = a * (1.0 - r1 - r2) + b * r1 + d * r2;
        }
        else {
            k = b * (1.0 - r1 - r2) + c * r1 + d * r2;
        }
    }
    else if (energy < e4) {
        const double phi42 = (energy - e2) / (e4 - e2);
        const double phi43 = (energy - e3) / (e4 - e3);
        const double phi41 = (energy - e1) / (e4 - e1);
        const Vector3 a = k21 * (1.0 - phi42) + k41 * phi42;
        const Vector3 b = k31 * (1.0 - phi43) + k41 * phi43;
        const Vector3 c = k41 * phi41;
        double r1 = urand(mt);
        double r2 = urand(mt);
        if (r1 + r2 > 1.0) {
            r1 = 1.0 - r1;
            r2 = 1.0 - r2;
        }
        k = a * (1.0 - r1 - r2) + b * r1 + c * r2;
    }
    else {
        cerr << "# Error in EnergyBand::getStateInTetrahedron()\n";
        cerr << "#   ===> Selected triangle does not contain given energy (energy > e3).\n";
        exit(EXIT_FAILURE);
    }

    State s;
    s.k.x = (k.x + grid_x[v1]) * (KMAX.x - KMIN.x) / NK + KMIN.x;
    s.k.y = (k.y + grid_y[v1]) * (KMAX.y - KMIN.y) / NK + KMIN.y;
    s.k.z = (k.z + grid_z[v1]) * (KMAX.z - KMIN.z) / NK + KMIN.z;
    s.n = nb;

    return s;
}



Vector3 EnergyBand::getWaveVectorDifference(const State initial_state, const int it) {

    Vector3 initial_cube = (initial_state.k - KMIN) / (KMAX - KMIN);
    initial_cube = (initial_cube * NK).floor();

    Vector3 final_cube(cube_x[it], cube_y[it], cube_z[it]);
    Vector3 delta_k = (final_cube - initial_cube) / NK;
    delta_k = delta_k.reduce_FCC();

    return delta_k;
}




#pragma once
#include <cmath>

// Physical constants from NIST
// http://physics.nist.gov/cuu/Constants/

// Speed of light in vacuum (m/sec)
constexpr double SPEED_OF_LIGHT = 299792458.0e+00;

// Elementary charge (C)
constexpr double ELEMENTARY_CHARGE = 1.602176462e-19;

// Plank constant over 2 pi (J*sec)
constexpr double HBAR = 1.054571596e-34;

// Electron mass (kg)
constexpr double ELECTRON_MASS = 9.10938188e-31;

// Boltzmann constant (J/K)
constexpr double BOLTZMANN = 1.3806503e-23;

// Electric constant (F/m)
constexpr double EPSILON_0 = 8.854187817e-12;

constexpr int SELF_SCATTERING = -1;
constexpr int PHONON_EMISSION = 0;
constexpr int PHONON_ABSORPTION = 100;
constexpr int IMPACT_IONIZATION = 200;

struct Vector3 {
    double x, y, z;

    constexpr Vector3(double vx = 0.0, double vy = 0.0, double vz = 0.0) : x(vx), y(vy), z(vz) {}

    Vector3 operator+(const Vector3& other) const {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    Vector3 operator-(const Vector3& other) const {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    Vector3 operator*(double scalar) const {
        return Vector3(x * scalar, y * scalar, z * scalar);
    }

    // Hadamard Product
    Vector3 operator*(const Vector3& other) const {
        return Vector3(x * other.x, y * other.y, z * other.z);
    }

    Vector3 operator/(double scalar) const {
        return Vector3(x / scalar, y / scalar, z / scalar);
    }

    // Hadamard Product
    Vector3 operator/(const Vector3& other) const {
        return Vector3(x / other.x, y / other.y, z / other.z);
    }

    double squared() const {
        return x * x + y * y + z * z;
    }

    Vector3 floor() const {
        return Vector3(std::floor(x), std::floor(y), std::floor(z));
    }

    Vector3 reduce() const {
        return Vector3(x - std::floor(x), y - std::floor(y), z - std::floor(z));
    }

    Vector3 reduce_FCC() const {
	double rx = x;
	double ry = y;
	double rz = z;
        double absx = fabs(x);
        double absy = fabs(y);
        double absz = fabs(z);
        if (absx > 1.0 || absy > 1.0 || absz > 1.0 || absx + absy + absz > 1.5) {
            // Out of 1st Brillouin Zone
            double gx = round(x * 0.5) * 2.0;
            double gy = round(y * 0.5) * 2.0;
            double gz = round(z * 0.5) * 2.0;
            if (fabs(x - gx) + fabs(y - gy) + fabs(z - gz) > 1.5) {
                gx = round(x * 0.5 - 0.5) * 2.0 + 1.0;
                gy = round(y * 0.5 - 0.5) * 2.0 + 1.0;
                gz = round(z * 0.5 - 0.5) * 2.0 + 1.0;
            }
	    rx -= gx; ry -= gy; rz -= gz;
        }
        return Vector3(rx, ry, rz);
    }
};

struct State {
    Vector3 k = Vector3(0.0, 0.0, 0.0);
    int n = 0;
    Vector3 r = Vector3(0.0, 0.0, 0.0);
};

#include "Fluid.h"

Fluid::Fluid(size_t nx, size_t ny, size_t nz, double rho, double mul) 
	: rho_(rho), mul_(mul), rhoField(nx, ny, nz, rho){};

const double& Fluid::getDensity(size_t x, size_t y, size_t z) {
	return rhoField(x, y, z);
};

const double& Fluid::getGamma() const {
	return mul_;
};
#include "VelocityField.h"

VelocityField::VelocityField(size_t nx, size_t ny, size_t nz, double value, vectorDouble direction)
	:Field(nx, ny, nz, value), direction_(direction) {};
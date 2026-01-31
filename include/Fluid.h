#pragma once
#include "Matrix.h"
#include "Field.h"
#include "typedef.h"

class Fluid : public Field {
public:
	Fluid() = default;

	Fluid(size_t nx, size_t ny, size_t nz, double rho, double mul);

	//返回特定单元的密度
	const double& getDensity(size_t x, size_t y, size_t z);

	//返回扩散系数Γ
	const double& getGamma() const;

private:
	double rho_;
	double mul_;     //此处mul指计算中的扩散系数gamma(Γ)
	Matrix rhoField;
};
#pragma once
#include "Matrix.h"
#include "typedef.h"

class CoefMatrix {
public:
	CoefMatrix() = default;

	CoefMatrix(size_t dim, size_t nx, size_t ny, size_t nz, double value = 1);

	//返回具体单元具体位置的离散系数
	double& getCoeff(size_t x, size_t y, size_t z, Coeffient_direction dir);

	//设定具体单元具体位置的离散系数
	void setCoeff(size_t x, size_t y, size_t z, Coeffient_direction dir, const double& value);

private:
	size_t dim_;             //维数
	size_t num_coef;         //离散方程的系数
	size_t n_x;              
	size_t n_y;
	size_t n_z;
	vectorMatrix coef_mat;   //系数矩阵

};
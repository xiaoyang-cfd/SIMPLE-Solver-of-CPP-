#include "CoefficientMatrix.h"

CoefMatrix::CoefMatrix(size_t dim, size_t nx, size_t ny, size_t nz, double value) 
	: n_x(nx), n_y(ny), n_z(nz) {
	//设定系数的个数
	if (dim == 2) {
		num_coef = 6;  //针对二维情况，系数分别为aP,aE,aW,aN,aS,Su
	} else if (dim == 3) {
		num_coef = 8;  //针对三维情况，系数分别为aP,aE,aW,aN,aS,Su,aT,aB
	} else {
		throw std::invalid_argument("dim must be two or three");
	}

	//初始化系数矩阵
	Matrix mat(nx, ny, nz, value);
	for (size_t i = 0; i < num_coef; ++i) {
		coef_mat.push_back(mat);
	}
};

//返回具体单元具体位置的离散系数
double& CoefMatrix::getCoeff(size_t x, size_t y, size_t z, Coeffient_direction dir) {
	return coef_mat[dir](x, y, z);
};

void CoefMatrix::setCoeff(size_t x, size_t y, size_t z, Coeffient_direction dir, const double& value) {
	coef_mat[dir](x, y, z) = value;
};
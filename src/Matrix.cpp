#include"Matrix.h"

Matrix::Matrix(size_t num_x, size_t num_y, size_t num_z, double value)
	          : n_x(num_x), n_y(num_y), n_z(num_z), mat(num_x* num_y* num_z, value) { }

double& Matrix::operator()(size_t x, size_t y, size_t z){
	size_t index = (z * n_y * n_x) + (y * n_x) + x;

	if (x >= n_x || y >= n_y || z >= n_z || index >= mat.size()) {
		throw std::out_of_range("The given element index is out of range");
	}
	return mat[index];
}

const double& Matrix::operator()(size_t x, size_t y, size_t z) const {
	size_t index = (z * n_y * n_x) + (y * n_x) + x;

	if (x >= n_x || y >= n_y || z >= n_z || index >= mat.size()) {
		throw std::out_of_range("The given element index is out of range");
	}
	return mat[index];
}

void Matrix::fill(double value) {
	std::fill(mat.begin(), mat.end(), value);
};
#include"Field.h"

Field::Field(size_t nx, size_t ny, size_t nz, double value) 
	: n_x(nx), n_y(ny), n_z(nz), mat_field(nx, ny, nz, value) {
}

std::unique_ptr<Field> Field::createField(size_t nx, size_t ny, size_t nz, double value) {
	return std::make_unique<Field>(nx, ny, nz, value);
};

void Field::printField() {
	std::cout << "[\n";
	for (size_t i = 0; i < n_x; ++i) {
		for (size_t j = 0; j < n_y; ++j) {
			std::cout << mat_field(i, j, 0) << " ";
		}
		std::cout << "\n";
	}
	std::cout << "]" << std::endl;
};

size_t Field::getX() {
	return n_x;
};

size_t Field::getY() {
	return n_y;
};

size_t Field::getZ() {
	return n_z;
};

double& Field::getMat(size_t x, size_t y, size_t z) {
	return mat_field(x, y, z);
};

void Field::copyMat(const Field& mat) {
	this->mat_field = mat.mat_field;
};

void Field::setMat(size_t x, size_t y, size_t z, const double& value) {
	mat_field(x, y, z) = value;
};

void Field::initial(double value) {
	mat_field.fill(value);
};
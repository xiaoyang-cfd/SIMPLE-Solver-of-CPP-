#include "ConstructMesh.h"

//三维形式的构造函数
ConstructMesh::ConstructMesh(vectorDouble range_x, vectorDouble range_y, vectorDouble range_z,
							 size_t nx, size_t ny, size_t nz): 
							n_x(nx), n_y(ny), n_z(nz), nf_x(nx+1), nf_y(ny+1), nf_z(nz+1), 
							cellCenter_x(nx, 0), cellCenter_y(ny, 0), cellCenter_z(nz, 0), 
							cellPoint_x(nx + 1, 0), cellPoint_y(ny + 1, 0), cellPoint_z(nz + 1, 0) {
	if (size(range_x) != 2 || size(range_y) != 2 || size(range_z) != 2) {
		throw std::invalid_argument("range of x/y/z must have specified maximum and minimum values");
	}

	//初始化网格各方向的最大值和最小值
	max_x = std::max(range_x[0], range_x[1]);
	min_x = std::min(range_x[0], range_x[1]);
	max_y = std::max(range_y[0], range_y[1]);
	min_y = std::min(range_y[0], range_y[1]);
	max_z = std::max(range_z[0], range_z[1]);
	min_z = std::min(range_z[0], range_z[1]);

	//计算网格间距
	dx = std::abs(range_x[0] - range_x[1]) / nx;
	dy = std::abs(range_y[0] - range_y[1]) / ny;
	dz = std::abs(range_z[0] - range_z[1]) / nz;

	//初始化网格节点坐标
	for (size_t i = 0; i < nx + 1; ++i) {
		cellPoint_x[i] = min_x + i * dx;
	}
	for (size_t j = 0; j < ny + 1; ++j) {
		cellPoint_y[j] = min_y + j * dy;
	}
	for (size_t k = 0; k < nz + 1; ++k) {
		cellPoint_z[k] = min_z + k * dz;
	}

	//初始化网格中心点坐标
	for (size_t i = 0; i < nx; ++i) {
		cellCenter_x[i] = 0.5 * (cellPoint_x[i] + cellPoint_x[i + 1]);
	}
	for (size_t j = 0; j < ny; ++j) {
		cellCenter_y[j] = 0.5 * (cellPoint_y[j] + cellPoint_y[j + 1]);
	}
	for (size_t k = 0; k < nz; ++k) {
		cellCenter_z[k] = 0.5 * (cellPoint_z[k] + cellPoint_z[k + 1]);
	}

}

//返回x方向的节点坐标
const double& ConstructMesh::getCellPointX(size_t numx) const {
	if (numx >= cellPoint_x.size()) {
		throw std::out_of_range("Index out of range");
	}
	return cellPoint_x[numx];
};

//返回y方向的节点坐标
const double& ConstructMesh::getCellPointY(size_t numy) const {
	if (numy >= cellPoint_y.size()) {
		throw std::out_of_range("Index out of range");
	}
	return cellPoint_y[numy];
};

//返回z方向的节点坐标
const double& ConstructMesh::getCellPointZ(size_t numz) const {
	if (numz >= cellPoint_z.size()) {
		throw std::out_of_range("Index out of range");
	}
	return cellPoint_z[numz];
};

const double& ConstructMesh::getMin(Direction dir) const {
	if (dir == X) {
		return min_x;
	} else if (dir == Y)
	{
		return min_y;
	} else if (dir == Z) {
		return min_z;
	} else {
		throw std::invalid_argument("must input correct direction");
	}
};

const double& ConstructMesh::getMax(Direction dir) const {
	if (dir == X) {
		return max_x;
	}
	else if (dir == Y)
	{
		return max_y;
	}
	else if (dir == Z) {
		return max_z;
	}
	else {
		throw std::invalid_argument("must input correct direction");
	}
};

const size_t& ConstructMesh::getNx() const {
	return n_x;
}
const size_t& ConstructMesh::getNy() const {
	return n_y;
}
const size_t& ConstructMesh::getNz() const {
	return n_z;
}

const double& ConstructMesh::getDx() const {
	return dx;
}
const double& ConstructMesh::getDy() const {
	return dy;
}
const double& ConstructMesh::getDz() const {
	return dz;
}

//打印网格信息
void ConstructMesh::print() {
	std::cout << "max_x:" << max_x << " min_x:" << min_x << std::endl;
	std::cout << "max_y:" << max_y << " min_y:" << min_y << std::endl;
	std::cout << "max_z:" << max_z << " min_z:" << min_z << std::endl;

	std::cout << "n_x:" << n_x << " n_y:" << n_y << " n_z:" << n_z << std::endl;

	std::cout << "nf_x:" << nf_x << " nf_y:" << nf_y << " nf_z:" << nf_z << std::endl;

	std::cout << "dx:" << dx << " dy:" << dy << " dz:" << dz << std::endl;

	std::cout << "coordinate_cell_center_x:";
	for (size_t i = 0; i < cellCenter_x.size(); ++i) {
		std::cout << cellCenter_x[i] << " ";
	}

	std::cout << "\ncoordinate_cell_center_y:";
	for (size_t i = 0; i < cellCenter_y.size(); ++i) {
		std::cout << cellCenter_y[i] << " ";
	}

	std::cout << "\ncoordinate_cell_center_z:";
	for (size_t i = 0; i < cellCenter_z.size(); ++i) {
		std::cout << cellCenter_z[i] << " ";
	}

	std::cout << "\ncoordinate_cell_point_x:";
	for (size_t i = 0; i < cellPoint_x.size(); ++i) {
		std::cout << cellPoint_x[i] << " ";
	}

	std::cout << "\ncoordinate_cell_point_y:";
	for (size_t i = 0; i < cellPoint_y.size(); ++i) {
		std::cout << cellPoint_y[i] << " ";
	}

	std::cout << "\ncoordinate_cell_point_z:";
	for (size_t i = 0; i < cellPoint_z.size(); ++i) {
		std::cout << cellPoint_z[i] << " ";
	}
};
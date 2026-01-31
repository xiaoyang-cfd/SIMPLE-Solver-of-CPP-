#pragma once
#include<iostream>
#include<vector>
#include "typedef.h"
#include "Matrix.h"

class ConstructMesh {
public:
	ConstructMesh() = default;

	ConstructMesh(vectorDouble range_x, vectorDouble range_y, vectorDouble range_z,
		          size_t nx, size_t ny, size_t nz);

	//打印网格信息
	void print();

	//返回节点坐标
	const double& getCellPointX(size_t numx) const; //x方向
	const double& getCellPointY(size_t numy) const; //y方向
	const double& getCellPointZ(size_t numz) const; //z方向

	//返回网格范围
	const double& getMin(Direction dir) const;   //最小值
	const double& getMax(Direction dir) const;   //最大值

	//返回网格单元数
	const size_t& getNx() const; //x方向
	const size_t& getNy() const; //y方向
	const size_t& getNz() const; //z方向

	//返回网格间距
	const double& getDx() const; //x方向
	const double& getDy() const; //y方向
	const double& getDz() const; //z方向

private:
	size_t dim;

	//网格范围
	double max_x;
	double max_y;
	double max_z;
	double min_x;
	double min_y;
	double min_z;

	//网格个数
	size_t n_x;
	size_t n_y;
	size_t n_z;

	//网格面个数
	size_t nf_x;
	size_t nf_y;
	size_t nf_z;

	//网格间距
	double dx;
	double dy;
	double dz;

	//网格中心坐标
	vectorDouble cellCenter_x;
	vectorDouble cellCenter_y;
	vectorDouble cellCenter_z;

	//网格节点坐标
	vectorDouble cellPoint_x;
	vectorDouble cellPoint_y;
	vectorDouble cellPoint_z;

};

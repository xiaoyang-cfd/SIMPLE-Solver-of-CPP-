#pragma once
#include"Matrix.h"

class Field {
public:
	//默认构造函数
	Field() = default;

	//参数构造函数
	Field(size_t nx, size_t ny, size_t nz, double value = 0);

	//虚析构函数
	virtual ~Field() = default;

	//创建该类型对象
	virtual std::unique_ptr<Field> createField(size_t nx, size_t ny, size_t nz, double value = 0);

	//打印场量(只打印xy平面)
	void printField();
	
	//返回nx
	size_t getX();

	//返回ny
	size_t getY();

	//返回nz
	size_t getZ();

	//返回特定位置的矩阵元素
	double& getMat(size_t x, size_t y, size_t z);

	//设定特定位置的矩阵元素
	void setMat(size_t x, size_t y, size_t z, const double& value);

	//拷贝函数，将特定场拷贝给该场
	void copyMat(const Field& mat);

	//重新初始化场，元素设置为value
	void initial(double value);

private:
	size_t n_x;
	size_t n_y;
	size_t n_z;
	Matrix mat_field;

};
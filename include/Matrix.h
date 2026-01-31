#pragma once
#include<vector>
#include<iostream>
#include<memory>
#include"typedef.h"

class Matrix {
public:
	Matrix() = default;

	Matrix(size_t num_x, size_t num_y, size_t num_z, double value = 0);

	double& operator()(size_t x, size_t y, size_t z);

	const double& operator()(size_t x, size_t y, size_t z) const;

	Matrix& operator=(const Matrix& mat) = default;

	//将矩阵mat中元素全部替换为value
	void fill(double value);

private:
	size_t n_x;   
	size_t n_y;
	size_t n_z;
	vectorDouble mat;

};

typedef std::vector<Matrix> vectorMatrix;
#pragma once
#include "Field.h"
#include "typedef.h"

class VelocityField : public Field {
public:
	//默认构造函数
	VelocityField() = default;

	//带参数的构造函数
	VelocityField(size_t nx, size_t ny, size_t nz, double value = 0, vectorDouble direction = {1, 0, 0});


private:
	vectorDouble direction_;


};
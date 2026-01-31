#pragma once
#include "Field.h"

class PressureField :public Field {
public:
	//默认构造函数
	PressureField() = default;

	//带参数的构造函数
	PressureField(size_t nx, size_t ny, size_t nz, double value = 0);

private:


};
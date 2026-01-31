#pragma once
#include <iostream>
#include <fstream>
#include <string>

//写线性迭代求解器的开头标签
void writeHeader_linSolver(std::string file_name);

//用于存储后处理结果的类
class PostProcess {
public:
	//默认初始化将所有成员变量初始化为0
	PostProcess();

	//带参数构造函数
	PostProcess(double norm2_u_, double norm2_v_, double norm2_w_, double norm2_p_,
		double max_norm2_u_, double max_norm2_v_, double max_norm2_w_, double max_norm2_p_);

	//返回最大的二阶范数
	double& get_MaxNorm2_u();
	double& get_MaxNorm2_v();
	double& get_MaxNorm2_w();
	double& get_MaxNorm2_p();
	double& get_MaxNorm2_pp();

	//返回当前二阶范数
	double& get_Norm2_u();
	double& get_Norm2_v();
	double& get_Norm2_w();
	double& get_Norm2_p();
	double& get_Norm2_pp();

	//设置最大的二阶范数
	void set_MaxNorm2_u(const double& value);
	void set_MaxNorm2_v(const double& value);
	void set_MaxNorm2_w(const double& value);
	void set_MaxNorm2_p(const double& value);
	void set_MaxNorm2_pp(const double& value);

	//设置当前二阶范数
	void set_Norm2_u(const double& value);
	void set_Norm2_v(const double& value);
	void set_Norm2_w(const double& value);
	void set_Norm2_p(const double& value);
	void set_Norm2_pp(const double& value);

private:
	double norm2_u;     //当前迭代的速度u二阶范数
	double norm2_v;     //当前迭代的速度v二阶范数
	double norm2_w;     //当前迭代的速度w二阶范数
	double norm2_p;     //当前迭代的压力p二阶范数
	double norm2_pp;    //当前迭代的压力修正值pp二阶范数

	double max_norm2_u; //迄今为止的最大速度u二阶范数
	double max_norm2_v; //迄今为止的最大速度v二阶范数
	double max_norm2_w; //迄今为止的最大速度w二阶范数
	double max_norm2_p; //迄今为止的最大压力p二阶范数
	double max_norm2_pp;//迄今为止的最大压力修正值pp二阶范数
};
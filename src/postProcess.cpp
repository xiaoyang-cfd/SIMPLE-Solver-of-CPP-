#include "postProcess.h"

void writeHeader_linSolver(std::string file_name) {
	//写入文件的标签
	std::ofstream lin_solution(file_name);
	if (!lin_solution.is_open()) {
		std::cerr << "Cannot open file: " << file_name << std::endl;
	}
	else {
		lin_solution << "No.iter_nonlin | No.iter_lin | norm2 | norm_max | res" << std::endl;
	}
}

PostProcess::PostProcess() :
	norm2_u(0.0), norm2_v(0.0), norm2_w(0.0), norm2_p(0.0),
	max_norm2_u(0.0), max_norm2_v(0.0), max_norm2_w(0.0), max_norm2_p(0.0) {};

PostProcess::PostProcess(double norm2_u_, double norm2_v_, double norm2_w_, double norm2_p_,
	double max_norm2_u_, double max_norm2_v_, double max_norm2_w_, double max_norm2_p_) :
	norm2_u(norm2_u_), norm2_v(norm2_v_), norm2_w(norm2_w_), norm2_p(norm2_p_),
	max_norm2_u(max_norm2_u_), max_norm2_v(max_norm2_v_), max_norm2_w(max_norm2_w_), max_norm2_p(max_norm2_p_) {};

//返回最大的二阶范数
double& PostProcess::get_MaxNorm2_u() { return max_norm2_u; }
double& PostProcess::get_MaxNorm2_v() { return max_norm2_v; }
double& PostProcess::get_MaxNorm2_w() { return max_norm2_w; }
double& PostProcess::get_MaxNorm2_p() { return max_norm2_p; }
double& PostProcess::get_MaxNorm2_pp() { return max_norm2_pp; }

//返回当前二阶范数
double& PostProcess::get_Norm2_u() { return norm2_u; }
double& PostProcess::get_Norm2_v() { return norm2_v; }
double& PostProcess::get_Norm2_w() { return norm2_w; }
double& PostProcess::get_Norm2_p() { return norm2_p; }
double& PostProcess::get_Norm2_pp() { return norm2_pp; }

//设置最大的二阶范数
void PostProcess::set_MaxNorm2_u(const double& value) { max_norm2_u = value; }
void PostProcess::set_MaxNorm2_v(const double& value) { max_norm2_v = value; }
void PostProcess::set_MaxNorm2_w(const double& value) { max_norm2_w = value; }
void PostProcess::set_MaxNorm2_p(const double& value) { max_norm2_p = value; }
void PostProcess::set_MaxNorm2_pp(const double& value) { max_norm2_pp = value; }

//设置当前二阶范数
void PostProcess::set_Norm2_u(const double& value) { norm2_u = value; }
void PostProcess::set_Norm2_v(const double& value) { norm2_v = value; }
void PostProcess::set_Norm2_w(const double& value) { norm2_w = value; }
void PostProcess::set_Norm2_p(const double& value) { norm2_p = value; }
void PostProcess::set_Norm2_pp(const double& value) { norm2_pp = value; }
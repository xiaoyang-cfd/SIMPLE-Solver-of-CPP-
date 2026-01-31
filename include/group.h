#pragma once
#include "velocityField.h"
#include "PressureField.h"
#include "CoefficientMatrix.h"

//用于存储所有速度场的结构体
struct VelocityFieldGroup {
	
	VelocityFieldGroup(
		VelocityField& u, VelocityField& v, VelocityField& w,
		VelocityField& uf_, VelocityField& vf_, VelocityField& wf_,
		VelocityField& u0, VelocityField& v0, VelocityField& w0)
		: cu(u),cv(v),cw(w), uf(uf_), vf(vf_), wf(wf_), cu0(u0), cv0(v0), cw0(w0) {}

	//禁用拷贝，防止浅拷贝
	VelocityFieldGroup(const VelocityFieldGroup&) = delete;
	VelocityFieldGroup& operator=(const VelocityFieldGroup&) = delete;

	VelocityField& cu;        //体心速度u
	VelocityField& cv;	      //体心速度v
	VelocityField& cw;	      //体心速度w
	VelocityField& uf;        //面心速度u
	VelocityField& vf;        //面心速度v
	VelocityField& wf;        //面心速度w
	VelocityField& cu0;       //旧的体心速度场u
	VelocityField& cv0;       //旧的体心速度场v
	VelocityField& cw0;       //旧的体心速度场w
};

//用于存储所有压力场的结构体
struct PressureFieldGroup {

	PressureFieldGroup(PressureField& p_, PressureField& p0_, PressureField& pp_, PressureField& pp0_)
		: p(p_), p0(p0_), pp(pp_), pp0(pp0_){}

	//禁用拷贝
	PressureFieldGroup(const PressureFieldGroup&) = delete;
	PressureFieldGroup& operator=(const PressureFieldGroup&) = delete;

	PressureField& p;         //体心压力场
	PressureField& p0;        //旧的体心压力场
	PressureField& pp;        //压力修正场
	PressureField& pp0;	      //旧的压力修正场 
};

//用于存储所有离散方程系数的结构体
struct CoefficientGroup {

	CoefficientGroup(CoefMatrix& coef_u_, CoefMatrix& coef_v_, CoefMatrix& coef_w_, CoefMatrix& coef_p_ )
		: coef_u(coef_u_), coef_v(coef_v_), coef_w(coef_w_), coef_p(coef_p_){}

	//禁用拷贝
	CoefficientGroup(const CoefficientGroup&) = delete;
	CoefficientGroup& operator=(const CoefficientGroup&) = delete;

	CoefMatrix& coef_u;       //速度u的动量方程系数
	CoefMatrix& coef_v;	      //速度v的动量方程系数
	CoefMatrix& coef_w;       //速度w的动量方程系数
	CoefMatrix& coef_p;       //压力修正方程系数
};

//用于存储收敛判定相关参数的结构体
struct ConverSetGroup {

	ConverSetGroup(double lin_tolerance_u_ = 1e-4, double lin_tolerance_p_ = 1e-6, double nonlin_tolerance_ = 1e-6)
		: lin_tolerance_u(lin_tolerance_u_), lin_tolerance_p(lin_tolerance_p_), nonlin_tolerance(nonlin_tolerance_){}

	double lin_tolerance_u;  //速度的线性迭代容差
	double lin_tolerance_p;  //压力的线性迭代容差
	double nonlin_tolerance; //非线性迭代的容差
};

//用于存储求解器相关参数的结构体
struct SolverConfigGroup {

	SolverConfigGroup(bool stop_switch_,
		Conv_scheme conv_scheme_, Iter_method iter_method_,
		size_t nonlin_iter_, size_t numIter_u_, size_t numIter_v_, size_t numIter_w_, size_t numIter_p_,
		double alpha_v_, double alpha_p_) 
		: stop_switch(stop_switch_), conv_scheme(conv_scheme_), iter_method(iter_method_),
		nonlin_iter(nonlin_iter_), numIter_u(numIter_u_), numIter_v(numIter_v_), numIter_w(numIter_w_), numIter_p(numIter_p_),
		alpha_v(alpha_v_), alpha_p(alpha_p_){}

	bool stop_switch;        //跳出循环的开关

	Conv_scheme conv_scheme; //对流项离散格式
	Iter_method iter_method; //线性迭代方法

	size_t nonlin_iter;      //非线性迭代次数
	size_t numIter_u;        //x方向动量方程线性迭代次数
	size_t numIter_v;        //y方向动量方程线性迭代次数
	size_t numIter_w;        //z方向动量方程线性迭代次数
	size_t numIter_p;        //压力修正方程线性迭代次数

	double alpha_v;          //速度松弛因子
	double alpha_p;          //压力松弛因子	
};
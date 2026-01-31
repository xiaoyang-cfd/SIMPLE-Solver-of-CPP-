#pragma once
#include "ConstructMesh.h"
#include "VelocityField.h"
#include "PressureField.h"
#include "CoefficientMatrix.h"
#include "Fluid.h"
#include "BoundaryCondition.h"
#include "group.h"
#include "iterative_solver.h"
#include "postProcess.h"
#include <cmath>
#include <fstream>
#include <string>
#include <functional>
#include <memory>


//定义一个simple求解器类，专门用于定义关于求解器的函数
class SIMPLE_Solver {
public:
	SIMPLE_Solver() = default;

	//构造函数
	SIMPLE_Solver(
		ConstructMesh& mesh,              //网格
		VelocityFieldGroup& velocity,     //速度场
		PressureFieldGroup& pressure,     //压力场
		CoefficientGroup& coefficient,    //离散方程系数
		Fluid& fluid,                     //流体属性
		BcField& bc,                      //边界条件
		ConverSetGroup& converSetting,    //收敛设置
		SolverConfigGroup& solverSetting) //求解器设置
		: m_mesh(mesh), m_velocity(velocity), m_pressure(pressure), m_coefficient(coefficient),
		m_fluid(fluid), m_bc(bc), m_converSetting(converSetting), m_solverSetting(solverSetting),
		m_postProcess(){}

	//存储旧的变量场
	void store_OldField();

	//SIMPLE求解器主函数
	void Solve(size_t num_iter, size_t fre_termiPrint, size_t nstep, std::vector<std::string> lin_fileName);

	//RC插值：体心速度计算面心速度(含迭代轮次判断)
	void InterRhieChow_FaceVel(size_t num_iter);

	//RC插值：体心速度计算面心速度
	void InterRhieChow_FaceVel2();

	//设置边界条件：设置面心速度
	void SetBC_FaceVel();

	//计算动量方程的离散系数(考虑边界条件)――总函数
	void CalcuCoeff_momentEqn(Conv_scheme conv_scheme);

	//计算动量方程的离散系数――子函数1：不考虑边界条件的对流扩散方程
	void CalcuCoeff_momentEqn_ConvAndDiff(Conv_scheme conv_scheme);

	//计算动量方程的离散系数――子函数2：考虑x方向动量方程的边界条件
	void CalcuCoeff_momentEqn_BC_u();

	//计算动量方程的离散系数――子函数3：考虑x方向动量方程的边界条件
	void CalcuCoeff_momentEqn_BC_v();

	//计算动量方程的离散系数――子函数4：考虑x方向动量方程的边界条件
	void CalcuCoeff_momentEqn_BC_w();

	//计算动量方程的离散系数――子函数5：考虑压力梯度项
	void CalcuCoeff_momentEqn_gradP();

	//线性迭代求解器，用于求解离散后的动量方程;参数说明：var―求解的变量类型(U、V、W)
	void IterativeSolver_vel(size_t iter_nonlin, size_t nstep, Variable var, size_t fre_termiPrint, std::string file_name_);

	//线性迭代求解器，用于求解离散后的压力修正方程;
	void IterativeSolver_p(size_t iter_nonlin, size_t nstep, size_t fre_termiPrint, std::string file_name_);

	//计算压力修正方程的离散系数
	void CalcuCoeff_continuityEqn();
	
	//更新压力场
	void update_Pressure();

	//更新体心速度
	void update_CenterVelocity();

	//更新面心速度场(仅更新内部面)
	void update_FaceVelocity();

	//计算残差
	void calculate_residual();

	//判断是否收敛并输出计算过程信息，如残差
	void check_convergence(size_t num_iter);

	//判断是否收敛的开关
	bool check_stop();

	//保存计算结果为VTK文件
	void saveVTK_Vel_Result(std::string filename);      //保存速度场
	void saveVTK_Pressure_Result(std::string filename); //保存压力场

	//计算残差范数
	void calculate_U_norm2();                           //速度场u二范数
	void calculate_V_norm2();                           //速度场v二范数
	void calculate_W_norm2();                           //速度场w二范数
	void calculate_P_norm2();                           //压力场二范数
	void calculate_Pp_norm2();                          //压力修正场二范数

	//===============================================
	//一些特定功能的子函数
	//参数说明：sign_F用于控制F系数的正负
	double a_nb(Conv_scheme conv_scheme, const double& area, const double& rho,
		const double& vel, const double& distance, double sign_F);                    //用于计算对流扩散方程中的单个离散系数

	//===============================================
	//一些别名函数，让主程序更具可阅读性
	double& Uf(size_t x, size_t y, size_t z);   //面心速度u
	double& Vf(size_t x, size_t y, size_t z);   //面心速度u
	double& Wf(size_t x, size_t y, size_t z);   //面心速度u

	double& Uc(size_t x, size_t y, size_t z);   //体心速度u
	double& Vc(size_t x, size_t y, size_t z);   //体心速度u
	double& Wc(size_t x, size_t y, size_t z);   //体心速度u

	double& Uc0(size_t x, size_t y, size_t z);  //旧的体心速度u
	double& Vc0(size_t x, size_t y, size_t z);  //旧的体心速度u
	double& Wc0(size_t x, size_t y, size_t z);  //旧的体心速度u

	double& P(size_t x, size_t y, size_t z);    //体心压力p
	double& P0(size_t x, size_t y, size_t z);   //旧的体心压力p0
	double& Pp(size_t x, size_t y, size_t z);   //体心压力修正值pp
	double& Pp0(size_t x, size_t y, size_t z);  //旧的体心压力修正值pp

	double& Coeff_U(size_t x, size_t y, size_t z, Coeffient_direction coeff_dir);   //U动量方程的离散系数
	double& Coeff_V(size_t x, size_t y, size_t z, Coeffient_direction coeff_dir);   //V动量方程的离散系数
	double& Coeff_W(size_t x, size_t y, size_t z, Coeffient_direction coeff_dir);   //W动量方程的离散系数

	double& L2_U();                             //速度u的二范数
	double& L2_V();                             //速度v的二范数
	double& L2_W();                             //速度w的二范数
	double& L2_P();                             //压力p的二范数
	double& L2_Pp();                            //压力修正值pp的二范数

	double& max_L2_U();                         //速度u的最大二范数
	double& max_L2_V();                         //速度v的最大二范数
	double& max_L2_W();                         //速度w的最大二范数
	double& max_L2_P();                         //压力p的最大二范数
	double& max_L2_Pp();                        //压力修正值pp的最大二范数

private:
	//以下私有成员使用引用，避免重复拷贝

	ConstructMesh& m_mesh;             //网格

	VelocityFieldGroup& m_velocity;    //速度场

	PressureFieldGroup& m_pressure;    //压力场

	CoefficientGroup& m_coefficient;   //离散方程系数

	Fluid& m_fluid;                    //流体属性

	BcField& m_bc;                     //边界条件

	ConverSetGroup& m_converSetting;   //收敛设置

	SolverConfigGroup& m_solverSetting;//求解器设置

	//特有成员
	PostProcess m_postProcess;         //后处理结果存储类

};
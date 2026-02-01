//作者：xiaoyang-cfd
//bilibili名称：奋斗的小阳同学
//编码格式：utf-8

#include<iostream>
#include <chrono>
#include "ConstructMesh.h"
#include "Matrix.h"
#include "Field.h"
#include "VelocityField.h"
#include "PressureField.h"
#include "CoefficientMatrix.h"
#include "Fluid.h"
#include "BoundaryCondition.h"
#include "SIMPLE_Solver.h"
#include "postProcess.h"


int main() {

	//设置基本参数=========================================
	auto start_time = std::chrono::high_resolution_clock::now(); //记录程序开始时间
	
	size_t dim = 3;                           //维度


	//构建网格信息
	vectorDouble range_x{ -0.0625,0.0625 };   //x方向网格范围
	vectorDouble range_y{ -0.0125,0.0125 };   //y方向网格范围
	vectorDouble range_z{ -0.0125,0.0125 };   //z方向网格范围

	//注意网格个数各方向均需大于1，受面插值算法的限制
	size_t nx = 16;                           //x方向网格个数
	size_t ny = 8;                            //y方向网格个数
	size_t nz = 8;                            //z方向网格个数
	

	//实例化网格信息
	ConstructMesh mesh(range_x, range_y, range_z, nx, ny, nz);  
	//mesh.print();   //打印网格信息


	//创建物理场
	//体心速度场，形参说明：x向单元数、y向单元数、z向单元数、元素初始值、速度方向
	VelocityField cu(nx, ny, nz, 0, { 1, 0, 0 });  
	VelocityField cv(nx, ny, nz, 0, { 0, 1, 0 });
	VelocityField cw(nx, ny, nz, 0, { 0, 0, 1 });
	//u.printField();

	//面心速度场
	//注意此处体心与面心的变化
	VelocityField uf(nx+1, ny, nz, 0, { 1, 0, 0 });
	VelocityField vf(nx, ny+1, nz, 0, { 0, 1, 0 });
	VelocityField wf(nx, ny, nz+1, 0, { 0, 0, 1 });
	//vf.printField();

	//旧的体心速度场
	VelocityField cu0(nx, ny, nz, 0, { 1, 0, 0 });
	VelocityField cv0(nx, ny, nz, 0, { 0, 1, 0 });
	VelocityField cw0(nx, ny, nz, 0, { 0, 0, 1 });

	//存储速度场
	VelocityFieldGroup velocityGroup(cu, cv, cw, uf, vf, wf, cu0, cv0, cw0);

	//压力场，形参说明：x向单元数、y向单元数、z向单元数、元素初始值
	PressureField p(nx, ny, nz, 0);
	//p.printField();

	//旧的压力场
	PressureField p0(nx, ny, nz, 0);

	//压力修正场
	PressureField pp(nx, ny, nz, 0);

	//旧的压力修正场
	PressureField pp0(nx, ny, nz, 0);

	//存储压力场
	PressureFieldGroup pressureGroup(p, p0, pp, pp0);


	//创建系数矩阵
	//动量方程系数矩阵，形参说明：维数、x向单元数、y向单元数、z向单元数、元素初始值
	CoefMatrix coef_u(dim, nx, ny, nz, 0);   //x方向
	CoefMatrix coef_v(dim, nx, ny, nz, 0);   //y方向
	CoefMatrix coef_w(dim, nx, ny, nz, 0);   //z方向

	//压力修正方程系数矩阵
	CoefMatrix coef_p(dim, nx, ny, nz, 0);

	//存储离散系数
	CoefficientGroup coeffGroup(coef_u, coef_v, coef_w, coef_p);

	
	//设置求解器参数
	size_t nonlin_iter = 800;                           //非线性迭代次数
	//double dt = 0.01;                                 //时间步长
	bool stop_switch = false;                           //跳出循环的开关
	Conv_scheme conv_scheme = Conv_scheme::UPWIND;      //对流项离散格式
	Iter_method iter_method = Iter_method::GAUSS_SEIDEL;//线性迭代方法
	size_t numIter_u = 5;                               //x方向动量方程线性迭代次数
	size_t numIter_v = 5;                               //y方向动量方程线性迭代次数
	size_t numIter_w = 5;                               //z方向动量方程线性迭代次数
	size_t numIter_p = 200;                             //压力修正方程线性迭代次数
	double alpha_v = 0.85;                              //速度松弛因子
	double alpha_p = 0.15;                              //压力松弛因子
	double lin_tolerance_u = 1e-4;                      //速度的线性迭代容差
	double lin_tolerance_p = 1e-6;                      //压力的线性迭代容差
	double nonlin_tolerance = 1e-6;                     //非线性迭代的容差
	//各个物理量的二阶范数(待补充)

	//存储收敛容差设置
	ConverSetGroup converSetGroup(lin_tolerance_u, lin_tolerance_p, nonlin_tolerance);
	
	//存储求解器设置
	SolverConfigGroup solverConfigGroup(stop_switch, conv_scheme, iter_method,
		nonlin_iter, numIter_u, numIter_v, numIter_w, numIter_p, alpha_v, alpha_p);

	//设置流体属性
	double rho = 1.0;                                   //密度
	double mul = 0.01;                                  //粘度
	Fluid fluid(nx, ny, nz, rho, mul);                  //初始化流体


	//设置几何边界条件
	double value = 0.0;
	double p_outlet = 0.0;                                          //压力的出口定值
	double pp_outlet = 0.0;                                         //压力修正值的出口定值          
	BcField bc(nx, ny, nz, mesh, value, p_outlet, pp_outlet);       //初始化边界，自动设置几何边界条件
	//bc.printType(0, 0, 0, 2);

	//设置物理边界条件，具体而言指将特定几何边界类型与特定物理边界类型关联
	Phy_BC xmin_phy = Phy_BC::INLET;
	Phy_BC xmax_phy = Phy_BC::OUTLET;
	Phy_BC ymin_phy = Phy_BC::WALL;
	Phy_BC ymax_phy = Phy_BC::WALL;
	Phy_BC zmin_phy = Phy_BC::WALL;
	Phy_BC zmax_phy = Phy_BC::WALL;
	//Phy_BC zmin_phy = Phy_BC::NONE;
	//Phy_BC zmax_phy = Phy_BC::NONE;
	bc.setPhyBC(xmin_phy, xmax_phy, ymin_phy, ymax_phy, zmin_phy, zmax_phy);
	
	//设置速度边界条件，具体而言指为所有单元的各个方向面设置速度边界条件
	std::vector<Num_BC> xmax_bc{ Num_BC::NONE, Num_BC::NONE, Num_BC::NONE };             //参数说明：xmax面的u、v、w三个方向的边界类型均为none
	std::vector<Num_BC> xmin_bc{ Num_BC::CONSTANT, Num_BC::CONSTANT, Num_BC::CONSTANT }; //该类型主要用于温度,后续可拓展温度场求解
	std::vector<Num_BC> ymax_bc{ Num_BC::CONSTANT, Num_BC::CONSTANT, Num_BC::CONSTANT };
	std::vector<Num_BC> ymin_bc{ Num_BC::CONSTANT, Num_BC::CONSTANT, Num_BC::CONSTANT };
	std::vector<Num_BC> zmax_bc{ Num_BC::CONSTANT, Num_BC::CONSTANT, Num_BC::CONSTANT };
	std::vector<Num_BC> zmin_bc{ Num_BC::CONSTANT, Num_BC::CONSTANT, Num_BC::CONSTANT };
	std::vector<double> xmax_vel{ 0, 0, 0 };
	std::vector<double> xmin_vel{ 1, 0, 0 };  //参数说明：分别表示u、v、w三个方向的速度，若该方向边界为constant则表示速度值，若该方向边界为flux则表示速度梯度
	std::vector<double> ymax_vel{ 0, 0, 0 };
	std::vector<double> ymin_vel{ 0, 0, 0 };
	std::vector<double> zmax_vel{ 0, 0, 0 };
	std::vector<double> zmin_vel{ 0, 0, 0 };
	bc.setVelTotalBC(xmin_bc, xmin_vel, xmax_bc, xmax_vel, 
					 ymin_bc, ymin_vel, ymax_bc, ymax_vel,
				     zmin_bc, zmin_vel, zmax_bc, zmax_vel);

	//设定输出相关参数
	size_t fre_termiPrint = 30;              //终端中打印迭代信息的频率
	size_t fre_saveVTK = 100;                //保存VTK文件的频率

	//线性迭代求解器结果文件的初始化
	std::string file_name_u = "lin_res_u";   //u动量方程求解的结果文件名
	writeHeader_linSolver(file_name_u);
	std::string file_name_v = "lin_res_v";   //v动量方程求解的结果文件名
	writeHeader_linSolver(file_name_v);
	std::string file_name_w = "lin_res_w";   //w动量方程求解的结果文件名
	writeHeader_linSolver(file_name_w);
	std::string file_name_p = "lin_res_p";   //p动量方程求解的结果文件名
	writeHeader_linSolver(file_name_p);
	std::vector<std::string>file_name_lin{ file_name_u, file_name_v, file_name_w, file_name_p };  //存储线性求解器的文件名


	//开始SIMPLE主循环===================================
	SIMPLE_Solver simple_solver(mesh, velocityGroup, pressureGroup, coeffGroup, fluid, bc, converSetGroup, solverConfigGroup);

	for (size_t iter = 1; iter < nonlin_iter + 1; ++iter) {
		if (iter == nonlin_iter - 1 || iter % fre_termiPrint == 0) {
			std::cout << "================= iteration: " << iter << " =================" << std::endl;
		}
		
		//存储旧的变量场
		simple_solver.store_OldField();

		//求解器部分
		simple_solver.Solve(iter, fre_termiPrint, nonlin_iter, file_name_lin);

		//保存结果文件
		if (iter == nonlin_iter || iter % fre_saveVTK == 0) {
			//保存速度场
			std::string vel_filename = "VelocityField_iter_" + std::to_string(iter) + ".vtk";
			simple_solver.saveVTK_Vel_Result(vel_filename);
			//保存压力场
			std::string pres_filename = "PressureField_iter_" + std::to_string(iter) + ".vtk";
			simple_solver.saveVTK_Pressure_Result(pres_filename);
		}

		//判断是否收敛,若满足精度提前退出循环
		if (simple_solver.check_stop()) {
			break;
		}
	}

	//记录程序结束时间
	auto end_time = std::chrono::high_resolution_clock::now();

	//计算并输出程序运行时间
	auto duration = std::chrono::duration<double>(end_time - start_time).count();
	std::cout << "Total computation time: " << duration << " seconds." << std::endl;

	return 0;
}


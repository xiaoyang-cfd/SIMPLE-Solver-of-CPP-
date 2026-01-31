#include "SIMPLE_Solver.h"

void SIMPLE_Solver::Solve(size_t num_iter, size_t fre_termiPrint, size_t nstep, std::vector<std::string> lin_fileName) {

	/*
	1.设置边界条件
	根据边界条件，调整边界处面心的物理量
	*/
	SetBC_FaceVel();

	/*
	2.计算面心速度
	此处需要进行判断：
	①若为第一轮，则直接采用中心差分，计算面心速度
	②若为后续轮，则采用RC插值，计算面心速度
	*/
	InterRhieChow_FaceVel(num_iter);

	/*
	3.计算动量方程系数
	①先不考虑边界条件，计算连接系数aP,aE,aW,aN,aS(aT,aB),Su
	第一步中，注意对流项的离散格式
	②考虑边界条件，调整计算系数aP,aE,aW,aN,aS(aT,aB),Su
	第二步中，应根据x、y(z)各方向考虑边界条件
	③考虑压力梯度项，调整系数Su
	第三步中，注意压力项离散后对于x、y(z)各方向的源项不同
	*/
	Conv_scheme conv_scheme = m_solverSetting.conv_scheme;      //读取对流离散格式
	CalcuCoeff_momentEqn(conv_scheme);

	/*
	4.求解动量方程，计算体心速度u*、v*(w*)
	属于线性迭代部分，可选择雅各比迭代、高斯赛德尔迭代
	应针对x、y(z)方向分别迭代求解
	*/
	IterativeSolver_vel(num_iter, nstep, Variable::U, fre_termiPrint, lin_fileName[0]);  //求解u方向动量方程
	IterativeSolver_vel(num_iter, nstep, Variable::V, fre_termiPrint, lin_fileName[1]);  //求解v方向动量方程
	IterativeSolver_vel(num_iter, nstep, Variable::W, fre_termiPrint, lin_fileName[2]);  //求解w方向动量方程

	/*
	5.计算面心速度
	使用RC插值方法计算
	*/
	InterRhieChow_FaceVel2();

	/*
	6.更新边界条件
	由于体心速度发生边界，需对面心速度根据边界条件进行修正
	*/
	SetBC_FaceVel();

	/*
	7.计算压力修正方程的系数
	包括apP、apE、apW、apN、apS、apT、apB
	*/
	CalcuCoeff_continuityEqn();

	/*
	8.求解压力修正方程，计算压力修正值pp
	属于线性迭代部分，可选择雅各比迭代、高斯赛德尔迭代
	*/
	IterativeSolver_p(num_iter, nstep, fre_termiPrint, lin_fileName[3]);   //求解压力修正方程

	/*
	9.更新速度场和压力场
	①更新压力场，显式松弛
	②更新体心速度场，显式松弛
	③更新面心速度场，显式松弛
	*/
	update_Pressure();
	update_CenterVelocity();
	update_FaceVelocity();

	/*
	10.判断是否收敛
	①判断误差是否满足要求，若满足停止循环，若不满足返回第一步
	②根据输出频率，终端中输出计算结果的误差
	*/
	calculate_residual();

	check_convergence(num_iter);

};

void SIMPLE_Solver::store_OldField() {
	//存储旧的压力场
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				P0(i, j, k) = P(i, j, k);
				Pp0(i, j, k) = Pp(i, j, k);
			}
		}
	}

	//存储旧的速度场
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				Uc0(i, j, k) = Uc(i, j, k);
				Vc0(i, j, k) = Vc(i, j, k);
				Wc0(i, j, k) = Wc(i, j, k);
			}
		}
	}
};

double& SIMPLE_Solver::Uf(size_t x, size_t y, size_t z) { return m_velocity.uf.getMat(x, y, z); }   //面心速度u
double& SIMPLE_Solver::Vf(size_t x, size_t y, size_t z) { return m_velocity.vf.getMat(x, y, z); }   //面心速度u
double& SIMPLE_Solver::Wf(size_t x, size_t y, size_t z) { return m_velocity.wf.getMat(x, y, z); }   //面心速度u

double& SIMPLE_Solver::Uc(size_t x, size_t y, size_t z) { return m_velocity.cu.getMat(x, y, z); }   //体心速度u
double& SIMPLE_Solver::Vc(size_t x, size_t y, size_t z) { return m_velocity.cv.getMat(x, y, z); }   //体心速度u
double& SIMPLE_Solver::Wc(size_t x, size_t y, size_t z) { return m_velocity.cw.getMat(x, y, z); }   //体心速度u

double& SIMPLE_Solver::Uc0(size_t x, size_t y, size_t z) { return m_velocity.cu0.getMat(x, y, z); } //旧的体心速度u
double& SIMPLE_Solver::Vc0(size_t x, size_t y, size_t z) { return m_velocity.cv0.getMat(x, y, z); } //旧的体心速度u
double& SIMPLE_Solver::Wc0(size_t x, size_t y, size_t z) { return m_velocity.cw0.getMat(x, y, z); } //旧的体心速度u

double& SIMPLE_Solver::P(size_t x, size_t y, size_t z) { return m_pressure.p.getMat(x, y, z); }     //体心压力p
double& SIMPLE_Solver::P0(size_t x, size_t y, size_t z) { return m_pressure.p0.getMat(x, y, z); }   //旧的体心压力p0
double& SIMPLE_Solver::Pp(size_t x, size_t y, size_t z) { return m_pressure.pp.getMat(x, y, z); }   //体心压力修正值pp
double& SIMPLE_Solver::Pp0(size_t x, size_t y, size_t z) { return m_pressure.pp0.getMat(x, y, z); } //体心压力修正值pp0

double& SIMPLE_Solver::Coeff_U(size_t x, size_t y, size_t z, Coeffient_direction coeff_dir) {
	return m_coefficient.coef_u.getCoeff(x, y, z, coeff_dir);
};
double& SIMPLE_Solver::Coeff_V(size_t x, size_t y, size_t z, Coeffient_direction coeff_dir) {
	return m_coefficient.coef_v.getCoeff(x, y, z, coeff_dir);
};
double& SIMPLE_Solver::Coeff_W(size_t x, size_t y, size_t z, Coeffient_direction coeff_dir) {
	return m_coefficient.coef_w.getCoeff(x, y, z, coeff_dir);
};

double& SIMPLE_Solver::L2_U() { return m_postProcess.get_Norm2_u(); }
double& SIMPLE_Solver::L2_V() { return m_postProcess.get_Norm2_v(); }
double& SIMPLE_Solver::L2_W() { return m_postProcess.get_Norm2_w(); }
double& SIMPLE_Solver::L2_P() { return m_postProcess.get_Norm2_p(); }
double& SIMPLE_Solver::L2_Pp() { return m_postProcess.get_Norm2_pp(); }

double& SIMPLE_Solver::max_L2_U() { return m_postProcess.get_MaxNorm2_u(); }
double& SIMPLE_Solver::max_L2_V() { return m_postProcess.get_MaxNorm2_v(); }
double& SIMPLE_Solver::max_L2_W() { return m_postProcess.get_MaxNorm2_w(); }
double& SIMPLE_Solver::max_L2_P() { return m_postProcess.get_MaxNorm2_p(); }
double& SIMPLE_Solver::max_L2_Pp() { return m_postProcess.get_MaxNorm2_pp(); }

void SIMPLE_Solver::InterRhieChow_FaceVel(size_t num_iter) {
	//若为第一次迭代，则使用中心差分法计算面速度
	if (num_iter == 1) {
		//x方向面速度
		for (size_t i = 1; i < m_mesh.getNx(); ++i) {
			for (size_t j = 0; j < m_mesh.getNy(); ++j) {
				for (size_t k = 0; k < m_mesh.getNz(); ++k) {
					Uf(i, j, k) = 0.5 * ( Uc(i, j, k) + Uc(i-1, j, k) );
					//std::cout << "CD vel_u" << "(" << i << ", " << j << ", " << k << ")" << std::endl;
				}
			}
		}

		//y方向面速度
		for (size_t i = 0; i < m_mesh.getNx(); ++i) {
			for (size_t j = 1; j < m_mesh.getNy(); ++j) {
				for (size_t k = 0; k < m_mesh.getNz(); ++k) {
					Vf(i, j, k) = 0.5 * ( Vc(i, j, k) + Vc(i, j - 1, k) );
					//std::cout << "CD vel_v" << "(" << i << ", " << j << ", " << k << ")" << std::endl;
				}
			}
		}
		
		//z方向面速度
		for (size_t i = 0; i < m_mesh.getNx(); ++i) {
			for (size_t j = 0; j < m_mesh.getNy(); ++j) {
				for (size_t k = 1; k < m_mesh.getNz(); ++k) {
					Wf(i, j, k) = 0.5 * (Wc(i, j, k) + Wc(i, j, k - 1));
					//std::cout << "CD vel_w" << "(" << i << ", " << j << ", " << k << ")" << std::endl;
				}
			}
		}
	}

	//非第一次迭代，则使用RC插值计算面心速度
	else {
		//x方向面速度
		/*二维计算公式为uw = 0.5 * (uW + uP) + 0.5 * (pP - pWW) * dy / (2.0 * aW) + 0.5 * (pE - pW) * dy / (2.0 * aP) - 0.5 * (1.0 / aW + 1.0 / aP) * (pP - pW) * dy
		三维计算公式需将dy->dy*dz
		//针对公式为便于编程理解，重新定义一下变量如下：
		变换关系为：uW->uW,uP->uE,pW->pW,pP->pE,pWW->pWW,pE->pEE,uw->uf,aW->aW,aP->aE
		几何关系如图：uWW-uW-uf-uE-uEE
		即该面左边的单元叫W，右边的单元叫E，左左边的单元叫WW，右右边的单元叫EE。
		计算公式变为：uf = 0.5 * (uW + uE) + 0.5 * (pE - pWW) * dy * dz / (2.0 * aW) + 0.5 * (pEE - pW) * dy * dz / (2.0 * aE) - 0.5 * (1.0 / aW + 1.0 / aE) * (pE - pW) * dy * dz
		*/
		double initial_pWW = 0.0;
		double initial_pEE = 0.0;
		const double& dx = m_mesh.getDx();
		const double& dy = m_mesh.getDy();
		const double& dz = m_mesh.getDz();
		for (size_t i = 1; i < m_mesh.getNx(); ++i) {
			for (size_t j = 0; j < m_mesh.getNy(); ++j) {
				for (size_t k = 0; k < m_mesh.getNz(); ++k) {
					double& uW = Uc(i - 1, j, k);
					double& uE = Uc(i, j, k);
					double& pE = P(i, j, k);
					double& pW = P(i - 1, j, k);

					//pWW和pEE有越界的风险，因此需要针对边界特殊处理
					double& pWW = initial_pWW;   //在判定前需提前定义
					if (i == 1) {                //左侧第一个内部面
						//pWW = initial_pWW;     //虽然注释掉，但实际效果在上方实现  
					}
					else {
						pWW = P(i - 2, j, k);
					}
					double& pEE = initial_pEE;
					if (i == m_mesh.getNx() - 1) { //右侧最后一个内部面
						//pEE = initial_pEE;       //虽然注释掉，但实际效果在上方实现  
					}
					else {
						pEE = P(i + 1, j, k);
					}

					//根据边界条件修改pWW和pEE
					const Geo_BC& geo_bc_w = m_bc.getGeoType(i - 1, j, k, fW);//当前面左侧的面的几何边界类型 
					const Geo_BC& geo_bc_e = m_bc.getGeoType(i, j, k, fE);    //当前面右侧的面的几何边界类型
					const Phy_BC& phy_bc_w = m_bc.getPhyType(geo_bc_w);       //当前面左侧的面的物理边界类型
					const Phy_BC& phy_bc_e = m_bc.getPhyType(geo_bc_e);       //当前面右侧的面的物理边界类型
					//修正pWW
					if (phy_bc_w == Phy_BC::INLET || phy_bc_w == Phy_BC::WALL) {
						pWW = P(i - 1, j, k);
					}
					else if (phy_bc_w == Phy_BC::OUTLET) {
						pWW = m_bc.getP_Outlet();
					}
					//修正pEE
					if (phy_bc_e == Phy_BC::INLET || phy_bc_e == Phy_BC::WALL) {
						pEE = P(i, j, k);
					}
					else if (phy_bc_e == Phy_BC::OUTLET) {
						pEE = m_bc.getP_Outlet();
					}

					//动量方程离散的系数
					const double& Coeff_W = m_coefficient.coef_u.getCoeff(i - 1, j, k, aP);  //面左侧单元系数ap
					const double& Coeff_E = m_coefficient.coef_u.getCoeff(i, j, k, aP);      //面右侧单元系数ap

					//RC插值主函数
					Uf(i, j, k) = 0.5 * (uW + uE)
						+ 0.5 * (pE - pWW) * dy * dz / (2.0 * Coeff_W)
						+ 0.5 * (pEE - pW) * dy * dz / (2.0 * Coeff_E)
						- 0.5 * (1.0 / Coeff_W + 1.0 / Coeff_E) * (pE - pW) * dy * dz;

					//std::cout << "RC vel_u" << "(" << i << ", " << j << ", " << k << ")" << std::endl;
				}
			}
		}

		//y方向面速度
		/*
		几何关系如图：uWW-uW-uf-uE-uEE
		即该面下边的单元叫W，上边的单元叫E，下下边的单元叫WW，上上边的单元叫EE。
		计算公式变为：uf = 0.5 * (uW + uE) + 0.5 * (pE - pWW) * dx * dz / (2.0 * aW) + 0.5 * (pEE - pW) * dx * dz / (2.0 * aE) - 0.5 * (1.0 / aW + 1.0 / aE) * (pE - pW) * dx * dz
		*/
		for (size_t i = 0; i < m_mesh.getNx(); ++i) {
			for (size_t j = 1; j < m_mesh.getNy(); ++j) {
				for (size_t k = 0; k < m_mesh.getNz(); ++k) {
					double& uW = Vc(i, j - 1, k);
					double& uE = Vc(i, j, k);
					double& pE = P(i, j, k);
					double& pW = P(i, j - 1, k);

					//pWW和pEE有越界的风险，因此需要针对边界特殊处理
					double& pWW = initial_pWW;   //在判定前需提前定义
					if (j == 1) {                //下侧第一个内部面
						//pWW = initial_pWW;     //虽然注释掉，但实际效果在上方实现  
					} else {
						pWW = P(i, j - 2, k);
					}	
					double& pEE = initial_pEE;
					if (j == m_mesh.getNy() - 1) { //上侧最后一个内部面
						//pEE = initial_pEE;       //虽然注释掉，但实际效果在上方实现  
					}
					else {
						pEE = P(i, j + 1, k);
					}

					//根据边界条件修改pWW和pEE
					const Geo_BC& geo_bc_w = m_bc.getGeoType(i, j - 1, k, fS);//当前面下侧的面的几何边界类型 
					const Geo_BC& geo_bc_e = m_bc.getGeoType(i, j, k, fN);    //当前面上侧的面的几何边界类型
					const Phy_BC& phy_bc_w = m_bc.getPhyType(geo_bc_w);       //当前面下侧的面的物理边界类型
					const Phy_BC& phy_bc_e = m_bc.getPhyType(geo_bc_e);       //当前面上侧的面的物理边界类型
					//修正pWW
					if (phy_bc_w == Phy_BC::INLET || phy_bc_w == Phy_BC::WALL) {
						pWW = P(i, j - 1, k);
					}
					else if (phy_bc_w == Phy_BC::OUTLET) {
						pWW = m_bc.getP_Outlet();
					}
					//修正pEE
					if (phy_bc_e == Phy_BC::INLET || phy_bc_e == Phy_BC::WALL) {
						pEE = P(i, j, k);
					}
					else if (phy_bc_e == Phy_BC::OUTLET) {
						pEE = m_bc.getP_Outlet();
					}

					//动量方程离散的系数
					const double& Coeff_W = m_coefficient.coef_v.getCoeff(i, j - 1, k, aP);  //面下侧单元系数ap
					const double& Coeff_E = m_coefficient.coef_v.getCoeff(i, j, k, aP);      //面上侧单元系数ap

					//RC插值主函数
					Vf(i, j, k) = 0.5 * (uW + uE)
						        + 0.5 * (pE - pWW) * dx * dz / (2.0 * Coeff_W)
						        + 0.5 * (pEE - pW) * dx * dz / (2.0 * Coeff_E)
						        - 0.5 * (1.0 / Coeff_W + 1.0 / Coeff_E) * (pE - pW) * dx * dz;

					//std::cout << "RC vel_v" << "(" << i << ", " << j << ", " << k << ")" << std::endl;
				}
			}
		}

		//z方向面速度
		/*
		几何关系如图：uWW-uW-uf-uE-uEE
		即该面底边的单元叫W，顶边的单元叫E，底底边的单元叫WW，顶顶边的单元叫EE。
		计算公式变为：uf = 0.5 * (uW + uE) + 0.5 * (pE - pWW) * dx * dy / (2.0 * aW) + 0.5 * (pEE - pW) * dx * dy / (2.0 * aE) - 0.5 * (1.0 / aW + 1.0 / aE) * (pE - pW) * dx * dy
		*/
		for (size_t i = 0; i < m_mesh.getNx(); ++i) {
			for (size_t j = 0; j < m_mesh.getNy(); ++j) {
				for (size_t k = 1; k < m_mesh.getNz(); ++k) {
					double& uW = Wc(i, j, k - 1);
					double& uE = Wc(i, j, k);
					double& pE = P(i, j, k);
					double& pW = P(i, j, k - 1);

					//pWW和pEE有越界的风险，因此需要针对边界特殊处理
					double& pWW = initial_pWW;   //在判定前需提前定义
					if (k == 1) {                //下侧第一个内部面
						//pWW = initial_pWW;     //虽然注释掉，但实际效果在上方实现  
					}
					else {
						pWW = P(i, j, k - 2);
					}
					double& pEE = initial_pEE;
					if (k == m_mesh.getNz() - 1) { //上侧最后一个内部面
						//pEE = initial_pEE;       //虽然注释掉，但实际效果在上方实现  
					}
					else {
						pEE = P(i, j, k + 1);
					}

					//根据边界条件修改pWW和pEE
					const Geo_BC& geo_bc_w = m_bc.getGeoType(i, j, k - 1, fB);//当前面底侧的面的几何边界类型 
					const Geo_BC& geo_bc_e = m_bc.getGeoType(i, j, k, fT);    //当前面顶侧的面的几何边界类型
					const Phy_BC& phy_bc_w = m_bc.getPhyType(geo_bc_w);       //当前面底侧的面的物理边界类型
					const Phy_BC& phy_bc_e = m_bc.getPhyType(geo_bc_e);       //当前面顶侧的面的物理边界类型
					//修正pWW
					if (phy_bc_w == Phy_BC::INLET || phy_bc_w == Phy_BC::WALL) {
						pWW = P(i, j, k - 1);
					}
					else if (phy_bc_w == Phy_BC::OUTLET) {
						pWW = m_bc.getP_Outlet();
					}
					//修正pEE
					if (phy_bc_e == Phy_BC::INLET || phy_bc_e == Phy_BC::WALL) {
						pEE = P(i, j, k);
					}
					else if (phy_bc_e == Phy_BC::OUTLET) {
						pEE = m_bc.getP_Outlet();
					}

					//动量方程离散的系数
					const double& Coeff_W = m_coefficient.coef_w.getCoeff(i, j, k - 1, aP);  //面底侧单元系数ap
					const double& Coeff_E = m_coefficient.coef_w.getCoeff(i, j, k, aP);      //面顶侧单元系数ap

					//RC插值主函数
					Wf(i, j, k) = 0.5 * (uW + uE)
						+ 0.5 * (pE - pWW) * dx * dy / (2.0 * Coeff_W)
						+ 0.5 * (pEE - pW) * dx * dy / (2.0 * Coeff_E)
						- 0.5 * (1.0 / Coeff_W + 1.0 / Coeff_E) * (pE - pW) * dx * dy;

					//std::cout << "RC vel_w" << "(" << i << ", " << j << ", " << k << ")" << std::endl;
				}
			}
		}
	}
};

void SIMPLE_Solver::InterRhieChow_FaceVel2() {
	//使用RC插值计算面心速度
 
	//x方向面速度
	/*二维计算公式为uw = 0.5 * (uW + uP) + 0.5 * (pP - pWW) * dy / (2.0 * aW) + 0.5 * (pE - pW) * dy / (2.0 * aP) - 0.5 * (1.0 / aW + 1.0 / aP) * (pP - pW) * dy
	三维计算公式需将dy->dy*dz
	//针对公式为便于编程理解，重新定义一下变量如下：
	变换关系为：uW->uW,uP->uE,pW->pW,pP->pE,pWW->pWW,pE->pEE,uw->uf,aW->aW,aP->aE
	几何关系如图：uWW-uW-uf-uE-uEE
	即该面左边的单元叫W，右边的单元叫E，左左边的单元叫WW，右右边的单元叫EE。
	计算公式变为：uf = 0.5 * (uW + uE) + 0.5 * (pE - pWW) * dy * dz / (2.0 * aW) + 0.5 * (pEE - pW) * dy * dz / (2.0 * aE) - 0.5 * (1.0 / aW + 1.0 / aE) * (pE - pW) * dy * dz
	*/
	double initial_pWW = 0.0;
	double initial_pEE = 0.0;
	const double& dx = m_mesh.getDx();
	const double& dy = m_mesh.getDy();
	const double& dz = m_mesh.getDz();
	for (size_t i = 1; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				double& uW = Uc(i - 1, j, k);
				double& uE = Uc(i, j, k);
				double& pE = P(i, j, k);
				double& pW = P(i - 1, j, k);

				//pWW和pEE有越界的风险，因此需要针对边界特殊处理
				double& pWW = initial_pWW;   //在判定前需提前定义
				if (i == 1) {                //左侧第一个内部面
					//pWW = initial_pWW;     //虽然注释掉，但实际效果在上方实现  
				}
				else {
					pWW = P(i - 2, j, k);
				}
				double& pEE = initial_pEE;
				if (i == m_mesh.getNx() - 1) { //右侧最后一个内部面
					//pEE = initial_pEE;       //虽然注释掉，但实际效果在上方实现  
				}
				else {
					pEE = P(i + 1, j, k);
				}

				//根据边界条件修改pWW和pEE
				const Geo_BC& geo_bc_w = m_bc.getGeoType(i - 1, j, k, fW);//当前面左侧的面的几何边界类型 
				const Geo_BC& geo_bc_e = m_bc.getGeoType(i, j, k, fE);    //当前面右侧的面的几何边界类型
				const Phy_BC& phy_bc_w = m_bc.getPhyType(geo_bc_w);       //当前面左侧的面的物理边界类型
				const Phy_BC& phy_bc_e = m_bc.getPhyType(geo_bc_e);       //当前面右侧的面的物理边界类型
				//修正pWW
				if (phy_bc_w == Phy_BC::INLET || phy_bc_w == Phy_BC::WALL) {
					pWW = P(i - 1, j, k);
				}
				else if (phy_bc_w == Phy_BC::OUTLET) {
					pWW = m_bc.getP_Outlet();
				}
				//修正pEE
				if (phy_bc_e == Phy_BC::INLET || phy_bc_e == Phy_BC::WALL) {
					pEE = P(i, j, k);
				}
				else if (phy_bc_e == Phy_BC::OUTLET) {
					pEE = m_bc.getP_Outlet();
				}

				//动量方程离散的系数
				const double& Coeff_W = m_coefficient.coef_u.getCoeff(i - 1, j, k, aP);  //面左侧单元系数ap
				const double& Coeff_E = m_coefficient.coef_u.getCoeff(i, j, k, aP);      //面右侧单元系数ap

				//RC插值主函数
				Uf(i, j, k) = 0.5 * (uW + uE)
					+ 0.5 * (pE - pWW) * dy * dz / (2.0 * Coeff_W)
					+ 0.5 * (pEE - pW) * dy * dz / (2.0 * Coeff_E)
					- 0.5 * (1.0 / Coeff_W + 1.0 / Coeff_E) * (pE - pW) * dy * dz;

				//std::cout << "RC vel_u" << "(" << i << ", " << j << ", " << k << ")" << std::endl;
			}
		}
	}

	//y方向面速度
	/*
	几何关系如图：uWW-uW-uf-uE-uEE
	即该面下边的单元叫W，上边的单元叫E，下下边的单元叫WW，上上边的单元叫EE。
	计算公式变为：uf = 0.5 * (uW + uE) + 0.5 * (pE - pWW) * dx * dz / (2.0 * aW) + 0.5 * (pEE - pW) * dx * dz / (2.0 * aE) - 0.5 * (1.0 / aW + 1.0 / aE) * (pE - pW) * dx * dz
	*/
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 1; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				double& uW = Vc(i, j - 1, k);
				double& uE = Vc(i, j, k);
				double& pE = P(i, j, k);
				double& pW = P(i, j - 1, k);

				//pWW和pEE有越界的风险，因此需要针对边界特殊处理
				double& pWW = initial_pWW;   //在判定前需提前定义
				if (j == 1) {                //下侧第一个内部面
					//pWW = initial_pWW;     //虽然注释掉，但实际效果在上方实现  
				}
				else {
					pWW = P(i, j - 2, k);
				}
				double& pEE = initial_pEE;
				if (j == m_mesh.getNy() - 1) { //上侧最后一个内部面
					//pEE = initial_pEE;       //虽然注释掉，但实际效果在上方实现  
				}
				else {
					pEE = P(i, j + 1, k);
				}

				//根据边界条件修改pWW和pEE
				const Geo_BC& geo_bc_w = m_bc.getGeoType(i, j - 1, k, fS);//当前面下侧的面的几何边界类型 
				const Geo_BC& geo_bc_e = m_bc.getGeoType(i, j, k, fN);    //当前面上侧的面的几何边界类型
				const Phy_BC& phy_bc_w = m_bc.getPhyType(geo_bc_w);       //当前面下侧的面的物理边界类型
				const Phy_BC& phy_bc_e = m_bc.getPhyType(geo_bc_e);       //当前面上侧的面的物理边界类型
				//修正pWW
				if (phy_bc_w == Phy_BC::INLET || phy_bc_w == Phy_BC::WALL) {
					pWW = P(i, j - 1, k);
				}
				else if (phy_bc_w == Phy_BC::OUTLET) {
					pWW = m_bc.getP_Outlet();
				}
				//修正pEE
				if (phy_bc_e == Phy_BC::INLET || phy_bc_e == Phy_BC::WALL) {
					pEE = P(i, j, k);
				}
				else if (phy_bc_e == Phy_BC::OUTLET) {
					pEE = m_bc.getP_Outlet();
				}

				//动量方程离散的系数
				const double& Coeff_W = m_coefficient.coef_v.getCoeff(i, j - 1, k, aP);  //面下侧单元系数ap
				const double& Coeff_E = m_coefficient.coef_v.getCoeff(i, j, k, aP);      //面上侧单元系数ap

				//RC插值主函数
				Vf(i, j, k) = 0.5 * (uW + uE)
					+ 0.5 * (pE - pWW) * dx * dz / (2.0 * Coeff_W)
					+ 0.5 * (pEE - pW) * dx * dz / (2.0 * Coeff_E)
					- 0.5 * (1.0 / Coeff_W + 1.0 / Coeff_E) * (pE - pW) * dx * dz;

				//std::cout << "RC vel_v" << "(" << i << ", " << j << ", " << k << ")" << std::endl;
			}
		}
	}

	//z方向面速度
	/*
	几何关系如图：uWW-uW-uf-uE-uEE
	即该面底边的单元叫W，顶边的单元叫E，底底边的单元叫WW，顶顶边的单元叫EE。
	计算公式变为：uf = 0.5 * (uW + uE) + 0.5 * (pE - pWW) * dx * dy / (2.0 * aW) + 0.5 * (pEE - pW) * dx * dy / (2.0 * aE) - 0.5 * (1.0 / aW + 1.0 / aE) * (pE - pW) * dx * dy
	*/
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 1; k < m_mesh.getNz(); ++k) {
				double& uW = Wc(i, j, k - 1);
				double& uE = Wc(i, j, k);
				double& pE = P(i, j, k);
				double& pW = P(i, j, k - 1);

				//pWW和pEE有越界的风险，因此需要针对边界特殊处理
				double& pWW = initial_pWW;   //在判定前需提前定义
				if (k == 1) {                //下侧第一个内部面
					//pWW = initial_pWW;     //虽然注释掉，但实际效果在上方实现  
				}
				else {
					pWW = P(i, j, k - 2);
				}
				double& pEE = initial_pEE;
				if (k == m_mesh.getNz() - 1) { //上侧最后一个内部面
					//pEE = initial_pEE;       //虽然注释掉，但实际效果在上方实现  
				}
				else {
					pEE = P(i, j, k + 1);
				}

				//根据边界条件修改pWW和pEE
				const Geo_BC& geo_bc_w = m_bc.getGeoType(i, j, k - 1, fB);//当前面底侧的面的几何边界类型 
				const Geo_BC& geo_bc_e = m_bc.getGeoType(i, j, k, fT);    //当前面顶侧的面的几何边界类型
				const Phy_BC& phy_bc_w = m_bc.getPhyType(geo_bc_w);       //当前面底侧的面的物理边界类型
				const Phy_BC& phy_bc_e = m_bc.getPhyType(geo_bc_e);       //当前面顶侧的面的物理边界类型
				//修正pWW
				if (phy_bc_w == Phy_BC::INLET || phy_bc_w == Phy_BC::WALL) {
					pWW = P(i, j, k - 1);
				}
				else if (phy_bc_w == Phy_BC::OUTLET) {
					pWW = m_bc.getP_Outlet();
				}
				//修正pEE
				if (phy_bc_e == Phy_BC::INLET || phy_bc_e == Phy_BC::WALL) {
					pEE = P(i, j, k);
				}
				else if (phy_bc_e == Phy_BC::OUTLET) {
					pEE = m_bc.getP_Outlet();
				}

				//动量方程离散的系数
				const double& Coeff_W = m_coefficient.coef_w.getCoeff(i, j, k - 1, aP);  //面底侧单元系数ap
				const double& Coeff_E = m_coefficient.coef_w.getCoeff(i, j, k, aP);      //面顶侧单元系数ap

				//RC插值主函数
				Wf(i, j, k) = 0.5 * (uW + uE)
					+ 0.5 * (pE - pWW) * dx * dy / (2.0 * Coeff_W)
					+ 0.5 * (pEE - pW) * dx * dy / (2.0 * Coeff_E)
					- 0.5 * (1.0 / Coeff_W + 1.0 / Coeff_E) * (pE - pW) * dx * dy;

				//std::cout << "RC vel_w" << "(" << i << ", " << j << ", " << k << ")" << std::endl;
			}
		}
	}
};

void SIMPLE_Solver::SetBC_FaceVel() {
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				//左面边界设置
				Face_direction dir = fW;
				const Geo_BC& geo_bc_w = m_bc.getGeoType(i, j, k, dir);//当前面左侧的面的几何边界类型 
				const Phy_BC& phy_bc_w = m_bc.getPhyType(geo_bc_w);    //当前面左侧的面的物理边界类型
				if (phy_bc_w == Phy_BC::INLET) {
					Uf(i, j, k) = m_bc.getBC_U(i, j, k, dir);
				}
				else if (phy_bc_w == Phy_BC::OUTLET) {
					Uf(i, j, k) = Uc(i, j, k);
				}
				else if (phy_bc_w == Phy_BC::WALL) {
					Uf(i, j, k) = m_bc.getBC_U(i, j, k, dir);
				}

				//右面边界设置
				dir = fE;
				const Geo_BC& geo_bc_e = m_bc.getGeoType(i, j, k, dir);//当前面右侧的面的几何边界类型 
				const Phy_BC& phy_bc_e = m_bc.getPhyType(geo_bc_e);    //当前面右侧的面的物理边界类型
				if (phy_bc_e == Phy_BC::INLET) {
					Uf(i + 1, j, k) = m_bc.getBC_U(i, j, k, dir);
				}
				else if (phy_bc_e == Phy_BC::OUTLET) {
					Uf(i + 1, j, k) = Uc(i, j, k);
				}
				else if (phy_bc_e == Phy_BC::WALL) {
					Uf(i + 1, j, k) = m_bc.getBC_U(i, j, k, dir);
				}

				//下面边界设置
				dir = fS;
				const Geo_BC& geo_bc_s = m_bc.getGeoType(i, j, k, dir);//当前面下侧的面的几何边界类型 
				const Phy_BC& phy_bc_s = m_bc.getPhyType(geo_bc_s);    //当前面下侧的面的物理边界类型
				if (phy_bc_s == Phy_BC::INLET) {
					Vf(i, j, k) = m_bc.getBC_V(i, j, k, dir);
				}
				else if (phy_bc_s == Phy_BC::OUTLET) {
					Vf(i, j, k) = Vc(i, j, k);
				}
				else if (phy_bc_s == Phy_BC::WALL) {
					Vf(i, j, k) = m_bc.getBC_V(i, j, k, dir);
				}

				//上面边界设置
				dir = fN;
				const Geo_BC& geo_bc_n = m_bc.getGeoType(i, j, k, dir);//当前面上侧的面的几何边界类型 
				const Phy_BC& phy_bc_n = m_bc.getPhyType(geo_bc_n);    //当前面上侧的面的物理边界类型
				if (phy_bc_n == Phy_BC::INLET) {
					Vf(i, j + 1, k) = m_bc.getBC_V(i, j, k, dir);
				}
				else if (phy_bc_n == Phy_BC::OUTLET) {
					Vf(i, j + 1, k) = Vc(i, j, k);
				}
				else if (phy_bc_n == Phy_BC::WALL) {
					Vf(i, j + 1, k) = m_bc.getBC_V(i, j, k, dir);
				}

				//底面边界设置
				dir = fB;
				const Geo_BC& geo_bc_b = m_bc.getGeoType(i, j, k, dir);//当前面下侧的面的几何边界类型 
				const Phy_BC& phy_bc_b = m_bc.getPhyType(geo_bc_b);    //当前面下侧的面的物理边界类型
				if (phy_bc_b == Phy_BC::INLET) {
					Wf(i, j, k) = m_bc.getBC_W(i, j, k, dir);
				}
				else if (phy_bc_b == Phy_BC::OUTLET) {
					Wf(i, j, k) = Wc(i, j, k);
				}
				else if (phy_bc_b == Phy_BC::WALL) {
					Wf(i, j, k) = m_bc.getBC_W(i, j, k, dir);
				}
				
				//顶面边界设置
				dir = fT;
				const Geo_BC& geo_bc_t = m_bc.getGeoType(i, j, k, dir);//当前面上侧的面的几何边界类型 
				const Phy_BC& phy_bc_t = m_bc.getPhyType(geo_bc_t);    //当前面上侧的面的物理边界类型
				if (phy_bc_t == Phy_BC::INLET) {
					Wf(i, j, k + 1) = m_bc.getBC_W(i, j, k, dir);
				}
				else if (phy_bc_t == Phy_BC::OUTLET) {
					Wf(i, j, k + 1) = Wc(i, j, k);
				}
				else if (phy_bc_t == Phy_BC::WALL) {
					Wf(i, j, k + 1) = m_bc.getBC_W(i, j, k, dir);
				}
			}
		}
	}
};

double SIMPLE_Solver::a_nb(Conv_scheme conv_scheme, const double& area, const double& rho, 
	const double& vel, const double& distance, double sign_F) {
	
	double D = m_fluid.getGamma() * area / distance;
	double F = rho * vel * area;

	//对流项CD格式 
	if (conv_scheme == Conv_scheme::CD) {
		return (D + 0.5 * sign_F * F);
	}
	//对流项upwind格式
	else if (conv_scheme == Conv_scheme::UPWIND) {
		return (D + std::max(0.0, sign_F * F));
	}
	//对流项hybrid格式
	else if (conv_scheme == Conv_scheme::HYBRID) {
		return std::max(sign_F * F, std::max(0.0, D + 0.5 * sign_F * F));
	}
	//对流项powerlaw格式
	else if (conv_scheme == Conv_scheme::POWERLAW) {
		return D * std::max(0.0, pow((1.0 - 0.1 * (F / D)), 5)) + std::max(0.0, sign_F * F);
	}
	else {
		throw std::invalid_argument("you must choose truly method with conv_scheme");
	}
};

void SIMPLE_Solver::CalcuCoeff_momentEqn_ConvAndDiff(Conv_scheme conv_scheme) {
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {

				//====计算aE=====================================

				double area = m_mesh.getDy() * m_mesh.getDz();                                      //面积A
				double rho;                                                                         //密度
				if (i == m_mesh.getNx() - 1) { //最右侧单元需特殊处理
					rho = m_fluid.getDensity(i, j, k);
				}
				else {
					rho = 0.5 * (m_fluid.getDensity(i, j, k) + m_fluid.getDensity(i + 1, j, k));    //密度rho，使用左右两侧单元密度平均处理
				}
				double vel = m_velocity.uf.getMat(i + 1, j, k);                                     //面速度uf
				double distance = m_mesh.getDx();                                                   //x方向网格间距
				double sign_F = -1.0;                                                               //扩散项常数F的正负号,E负、W正、N负、S正、T负、B正 
				double Fe = rho * vel * area;                                                       //E面扩散项常数Fe

				double coeff_aE = a_nb(conv_scheme, area, rho, vel, distance, sign_F);

				//存储aE系数
				Coeffient_direction dir = aE;
				m_coefficient.coef_u.setCoeff(i, j, k, dir, coeff_aE);
				m_coefficient.coef_v.setCoeff(i, j, k, dir, coeff_aE);
				m_coefficient.coef_w.setCoeff(i, j, k, dir, coeff_aE);



				//====计算aW=====================================
				area = m_mesh.getDy() * m_mesh.getDz();
				if (i == 0) {                    //最左侧单元需特殊处理
					rho = m_fluid.getDensity(i, j, k);
				}
				else {
					rho = 0.5 * (m_fluid.getDensity(i, j, k) + m_fluid.getDensity(i - 1, j, k));    //密度rho，使用左右两侧单元密度平均处理
				}
				vel = m_velocity.uf.getMat(i, j, k);
				distance = m_mesh.getDx();                                                          //x方向网格间距
				sign_F = 1.0;                                                                       //扩散项常数F的正负号,E负、W正、N负、S正、T负、B正 
				double Fw = rho * vel * area;                                                       //W面扩散项常数Fw

				double coeff_aW = a_nb(conv_scheme, area, rho, vel, distance, sign_F);

				//存储aW系数
				dir = aW;
				m_coefficient.coef_u.setCoeff(i, j, k, dir, coeff_aW);
				m_coefficient.coef_v.setCoeff(i, j, k, dir, coeff_aW);
				m_coefficient.coef_w.setCoeff(i, j, k, dir, coeff_aW);



				//====计算aN=====================================
				area = m_mesh.getDx() * m_mesh.getDz();
				if (j == m_mesh.getNy() - 1) {   //最上侧单元需特殊处理
					rho = m_fluid.getDensity(i, j, k);
				}
				else {
					rho = 0.5 * (m_fluid.getDensity(i, j, k) + m_fluid.getDensity(i, j + 1, k));    //密度rho，使用上下两侧单元密度平均处理
				}
				vel = m_velocity.vf.getMat(i, j + 1, k);
				distance = m_mesh.getDy();                                                          //y方向网格间距
				sign_F = -1.0;                                                                      //扩散项常数F的正负号,E负、W正、N负、S正、T负、B正 
				double Fn = rho * vel * area;                                                       //N面扩散项常数Fn

				double coeff_aN = a_nb(conv_scheme, area, rho, vel, distance, sign_F);

				//存储aN系数
				dir = aN;
				m_coefficient.coef_u.setCoeff(i, j, k, dir, coeff_aN);
				m_coefficient.coef_v.setCoeff(i, j, k, dir, coeff_aN);
				m_coefficient.coef_w.setCoeff(i, j, k, dir, coeff_aN);



				//====计算aS=====================================
				area = m_mesh.getDx() * m_mesh.getDz();
				if (j == 0) {                    //最下侧单元需特殊处理
					rho = m_fluid.getDensity(i, j, k);
				}
				else {
					rho = 0.5 * (m_fluid.getDensity(i, j, k) + m_fluid.getDensity(i, j - 1, k));    //密度rho，使用上下两侧单元密度平均处理
				}
				vel = m_velocity.vf.getMat(i, j, k);
				distance = m_mesh.getDy();                                                          //y方向网格间距
				sign_F = 1.0;                                                                       //扩散项常数F的正负号,E负、W正、N负、S正、T负、B正 
				double Fs = rho * vel * area;                                                       //S面扩散项常数Fs

				double coeff_aS = a_nb(conv_scheme, area, rho, vel, distance, sign_F);

				//存储aS系数
				dir = aS;
				m_coefficient.coef_u.setCoeff(i, j, k, dir, coeff_aS);
				m_coefficient.coef_v.setCoeff(i, j, k, dir, coeff_aS);
				m_coefficient.coef_w.setCoeff(i, j, k, dir, coeff_aS);



				//====计算aT=====================================
				area = m_mesh.getDx() * m_mesh.getDy();
				if (k == m_mesh.getNz() - 1) {   //最顶侧单元需特殊处理
					rho = m_fluid.getDensity(i, j, k);
				}
				else {
					rho = 0.5 * (m_fluid.getDensity(i, j, k) + m_fluid.getDensity(i, j, k + 1));    //密度rho，使用顶底两侧单元密度平均处理
				}
				vel = m_velocity.wf.getMat(i, j, k + 1);
				distance = m_mesh.getDz();                                                          //z方向网格间距
				sign_F = -1.0;                                                                      //扩散项常数F的正负号,E负、W正、N负、S正、T负、B正 
				double Ft = rho * vel * area;                                                       //T面扩散项常数Ft

				double coeff_aT = a_nb(conv_scheme, area, rho, vel, distance, sign_F);

				//存储aT系数
				dir = aT;
				m_coefficient.coef_u.setCoeff(i, j, k, dir, coeff_aT);
				m_coefficient.coef_v.setCoeff(i, j, k, dir, coeff_aT);
				m_coefficient.coef_w.setCoeff(i, j, k, dir, coeff_aT);



				//====计算aB=====================================
				area = m_mesh.getDx() * m_mesh.getDy();
				if (k == 0) {                   //最底侧单元需特殊处理
					rho = m_fluid.getDensity(i, j, k);
				}
				else {
					rho = 0.5 * (m_fluid.getDensity(i, j, k) + m_fluid.getDensity(i, j, k - 1));    //密度rho，使用顶底两侧单元密度平均处理
				}
				vel = m_velocity.wf.getMat(i, j, k);
				distance = m_mesh.getDz();                                                          //z方向网格间距
				sign_F = 1.0;                                                                       //扩散项常数F的正负号,E负、W正、N负、S正、T负、B正 
				double Fb = rho * vel * area;                                                       //B面扩散项常数Fb

				double coeff_aB = a_nb(conv_scheme, area, rho, vel, distance, sign_F);

				//存储aB系数
				dir = aB;
				m_coefficient.coef_u.setCoeff(i, j, k, dir, coeff_aB);
				m_coefficient.coef_v.setCoeff(i, j, k, dir, coeff_aB);
				m_coefficient.coef_w.setCoeff(i, j, k, dir, coeff_aB);



				//====计算aP=====================================
				double coeff_aP = coeff_aE + coeff_aW + coeff_aN + coeff_aS + coeff_aT + coeff_aB
					+ Fe - Fw + Fn - Fs + Ft - Fb;

				//存储aP系数
				dir = aP;
				m_coefficient.coef_u.setCoeff(i, j, k, dir, coeff_aP);
				m_coefficient.coef_v.setCoeff(i, j, k, dir, coeff_aP);
				m_coefficient.coef_w.setCoeff(i, j, k, dir, coeff_aP);

				//====设置Su=====================================
				dir = Su;
				double coeff_Su = 0.0;
				m_coefficient.coef_u.setCoeff(i, j, k, dir, coeff_Su);
				m_coefficient.coef_v.setCoeff(i, j, k, dir, coeff_Su);
				m_coefficient.coef_w.setCoeff(i, j, k, dir, coeff_Su);
			}
		}
	}
};

void SIMPLE_Solver::CalcuCoeff_momentEqn_BC_u() {
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				double Coeff_a_nb;   //初始化变量，用于处理边界处相邻单元的系数
				double Coeff_a_p;    //初始化变量，用于处理单元形心的系数
				double Coeff_Su;     //初始化变量，用于处理源项系数Su

				//调整E面边界的系数
				Face_direction face_dir = fE;
				Coeffient_direction coeff_dir = aE;
				Geo_BC geo = m_bc.getGeoType(i, j, k, face_dir);
				Phy_BC phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aE
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aE-Sp=aP+aE,因为Sp等于-2aE
					Coeff_Su = m_coefficient.coef_u.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_U(i, j, k, face_dir);               //计算Su，即Su=Su0+2*aE*uE

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aE
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aE
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aE
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aE
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aE
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aE-Sp=aP+aE,因为Sp等于-2aE
					Coeff_Su = m_coefficient.coef_u.getCoeff(i, j, k, Su) +
						2 * Coeff_a_nb * m_bc.getBC_U(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aE*uE

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aE
				}

				//调整W面边界的系数
				face_dir = fW;
				coeff_dir = aW;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aW
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aW-Sp=aP+aW,因为Sp等于-2aW
					Coeff_Su = m_coefficient.coef_u.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_U(i, j, k, face_dir);               //计算Su，即Su=Su0+2*aW*uW

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aW
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aW
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aW
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aW
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aW
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aW-Sp=aP+aW,因为Sp等于-2aW
					Coeff_Su = m_coefficient.coef_u.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_U(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aW*uW

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aW
				}

				//调整N面边界的系数
				face_dir = fN;
				coeff_dir = aN;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aN
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aN-Sp=aP+aN,因为Sp等于-2aN
					Coeff_Su = m_coefficient.coef_u.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_U(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aN*uN

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aN
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aN
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aN
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aN
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aN
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aN-Sp=aP+aN,因为Sp等于-2aN
					Coeff_Su = m_coefficient.coef_u.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_U(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aN*uN

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aN
				}

				//调整S面边界的系数
				face_dir = fS;
				coeff_dir = aS;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aS
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aS-Sp=aP+aS,因为Sp等于-2aS
					Coeff_Su = m_coefficient.coef_u.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_U(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aS*uS

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aS
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aS
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aS
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aS
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aS
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aS-Sp=aP+aS,因为Sp等于-2aS
					Coeff_Su = m_coefficient.coef_u.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_U(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aS*uS

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aS
				}

				//调整T面边界的系数
				face_dir = fT;
				coeff_dir = aT;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aT
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aT-Sp=aP+aT,因为Sp等于-2aT
					Coeff_Su = m_coefficient.coef_u.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_U(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aT*uT

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aT
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aT
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aT
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aT
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aT
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aT-Sp=aP+aT,因为Sp等于-2aT
					Coeff_Su = m_coefficient.coef_u.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_U(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aT*uT

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aT
				}

				//调整B面边界的系数
				face_dir = fB;
				coeff_dir = aB;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aB
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aB-Sp=aP+aB,因为Sp等于-2aB
					Coeff_Su = m_coefficient.coef_u.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_U(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aB*uB

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aB
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aB
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aB
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aB
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_u.getCoeff(i, j, k, coeff_dir);       //读取aB
					Coeff_a_p = m_coefficient.coef_u.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aB-Sp=aP+aB,因为Sp等于-2aB
					Coeff_Su = m_coefficient.coef_u.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_U(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aB*uB

					m_coefficient.coef_u.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_u.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_u.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aB
				}
			}
		}
	}
};

void SIMPLE_Solver::CalcuCoeff_momentEqn_BC_v() {
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				double Coeff_a_nb;   //初始化变量，用于处理边界处相邻单元的系数
				double Coeff_a_p;    //初始化变量，用于处理单元形心的系数
				double Coeff_Su;     //初始化变量，用于处理源项系数Su

				//调整E面边界的系数
				Face_direction face_dir = fE;
				Coeffient_direction coeff_dir = aE;
				Geo_BC geo = m_bc.getGeoType(i, j, k, face_dir);
				Phy_BC phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aE
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aE-Sp=aP+aE,因为Sp等于-2aE
					Coeff_Su = m_coefficient.coef_v.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_V(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aE*uE

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aE
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aE
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aE
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aE
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aE
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aE-Sp=aP+aE,因为Sp等于-2aE
					Coeff_Su = m_coefficient.coef_v.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_V(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aE*uE

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aE
				}

				//调整W面边界的系数
				face_dir = fW;
				coeff_dir = aW;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aW
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aW-Sp=aP+aW,因为Sp等于-2aW
					Coeff_Su = m_coefficient.coef_v.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_V(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aW*uW

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aW
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aW
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aW
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aW
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aW
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aW-Sp=aP+aW,因为Sp等于-2aW
					Coeff_Su = m_coefficient.coef_v.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_V(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aW*uW

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aW
				}

				//调整N面边界的系数
				face_dir = fN;
				coeff_dir = aN;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aN
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aN-Sp=aP+aN,因为Sp等于-2aN
					Coeff_Su = m_coefficient.coef_v.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_V(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aN*uN

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aN
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aN
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aN
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aN
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aN
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aN-Sp=aP+aN,因为Sp等于-2aN
					Coeff_Su = m_coefficient.coef_v.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_V(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aN*uN

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aN
				}

				//调整S面边界的系数
				face_dir = fS;
				coeff_dir = aS;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aS
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aS-Sp=aP+aS,因为Sp等于-2aS
					Coeff_Su = m_coefficient.coef_v.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_V(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aS*uS

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aS
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aS
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aS
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aS
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aS
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aS-Sp=aP+aS,因为Sp等于-2aS
					Coeff_Su = m_coefficient.coef_v.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_V(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aS*uS

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aS
				}

				//调整T面边界的系数
				face_dir = fT;
				coeff_dir = aT;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aT
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aT-Sp=aP+aT,因为Sp等于-2aT
					Coeff_Su = m_coefficient.coef_v.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_V(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aT*uT

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aT
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aT
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aT
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aT
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aT
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aT-Sp=aP+aT,因为Sp等于-2aT
					Coeff_Su = m_coefficient.coef_v.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_V(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aT*uT

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aT
				}

				//调整B面边界的系数
				face_dir = fB;
				coeff_dir = aB;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aB
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aB-Sp=aP+aB,因为Sp等于-2aB
					Coeff_Su = m_coefficient.coef_v.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_V(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aB*uB

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aB
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aB
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aB
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aB
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_v.getCoeff(i, j, k, coeff_dir);       //读取aB
					Coeff_a_p = m_coefficient.coef_v.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aB-Sp=aP+aB,因为Sp等于-2aB
					Coeff_Su = m_coefficient.coef_v.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_V(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aB*uB

					m_coefficient.coef_v.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_v.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_v.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aB
				}
			}
		}
	}
};

void SIMPLE_Solver::CalcuCoeff_momentEqn_BC_w() {
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				double Coeff_a_nb;   //初始化变量，用于处理边界处相邻单元的系数
				double Coeff_a_p;    //初始化变量，用于处理单元形心的系数
				double Coeff_Su;     //初始化变量，用于处理源项系数Su

				//调整E面边界的系数
				Face_direction face_dir = fE;
				Coeffient_direction coeff_dir = aE;
				Geo_BC geo = m_bc.getGeoType(i, j, k, face_dir);
				Phy_BC phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aE
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aE-Sp=aP+aE,因为Sp等于-2aE
					Coeff_Su = m_coefficient.coef_w.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_W(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aE*uE

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aE
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aE
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aE
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aE
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aE
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aE-Sp=aP+aE,因为Sp等于-2aE
					Coeff_Su = m_coefficient.coef_w.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_W(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aE*uE

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aE
				}

				//调整W面边界的系数
				face_dir = fW;
				coeff_dir = aW;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aW
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aW-Sp=aP+aW,因为Sp等于-2aW
					Coeff_Su = m_coefficient.coef_w.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_W(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aW*uW

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aW
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aW
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aW
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aW
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aW
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aW-Sp=aP+aW,因为Sp等于-2aW
					Coeff_Su = m_coefficient.coef_w.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_W(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aW*uW

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aW
				}

				//调整N面边界的系数
				face_dir = fN;
				coeff_dir = aN;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aN
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aN-Sp=aP+aN,因为Sp等于-2aN
					Coeff_Su = m_coefficient.coef_w.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_W(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aN*uN

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aN
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aN
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aN
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aN
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aN
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aN-Sp=aP+aN,因为Sp等于-2aN
					Coeff_Su = m_coefficient.coef_w.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_W(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aN*uN

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aN
				}

				//调整S面边界的系数
				face_dir = fS;
				coeff_dir = aS;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aS
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aS-Sp=aP+aS,因为Sp等于-2aS
					Coeff_Su = m_coefficient.coef_w.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_W(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aS*uS

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aS
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aS
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aS
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aS
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aS
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aS-Sp=aP+aS,因为Sp等于-2aS
					Coeff_Su = m_coefficient.coef_w.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_W(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aS*uS

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aS
				}

				//调整T面边界的系数
				face_dir = fT;
				coeff_dir = aT;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aT
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aT-Sp=aP+aT,因为Sp等于-2aT
					Coeff_Su = m_coefficient.coef_w.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_W(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aT*uT

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aT
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aT
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aT
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aT
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aT
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aT-Sp=aP+aT,因为Sp等于-2aT
					Coeff_Su = m_coefficient.coef_w.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_W(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aT*uT

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aT
				}

				//调整B面边界的系数
				face_dir = fB;
				coeff_dir = aB;
				geo = m_bc.getGeoType(i, j, k, face_dir);
				phyBC = m_bc.getPhyType(geo);
				if (phyBC == Phy_BC::INLET) {  //入口边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aB
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aB-Sp=aP+aB,因为Sp等于-2aB
					Coeff_Su = m_coefficient.coef_w.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_W(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aB*uB

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aB
				}
				else if (phyBC == Phy_BC::OUTLET) {  //出口边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aB
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) - Coeff_a_nb;  //计算aP，即aP-aB
					//Coeff_Su = 0.0;                                                       //计算Su，即Su=0

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					//m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aB
				}
				else if (phyBC == Phy_BC::WALL) {  //壁面边界
					Coeff_a_nb = m_coefficient.coef_w.getCoeff(i, j, k, coeff_dir);       //读取aB
					Coeff_a_p = m_coefficient.coef_w.getCoeff(i, j, k, aP) + Coeff_a_nb;  //计算aP，即aP-aB-Sp=aP+aB,因为Sp等于-2aB
					Coeff_Su = m_coefficient.coef_w.getCoeff(i, j, k, Su) +
						2.0 * Coeff_a_nb * m_bc.getBC_W(i, j, k, face_dir);                 //计算Su，即Su=Su0+2*aB*uB

					m_coefficient.coef_w.setCoeff(i, j, k, aP, Coeff_a_p);                //调整系数aP
					m_coefficient.coef_w.setCoeff(i, j, k, Su, Coeff_Su);                 //调整系数Su
					m_coefficient.coef_w.setCoeff(i, j, k, coeff_dir, 0.0);               //置零系数aB
				}
			}
		}
	}
};

void SIMPLE_Solver::CalcuCoeff_momentEqn_gradP() {
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				//初始化两侧单元的压力变量，用于计算压力梯度离散项
				double p_w;        //u方向指Pw、v方向指Ps、w方向指Pb
				double p_e;        //u方向指Pe、v方向指Pn、w方向指Pt

				//初始化离散方程中的源项
				double Su_0;        //用于存储没考虑压力梯度项前的源项
				double Su_P;         //用于存储考虑压力梯度项后的源项
				
				//u方向的动量方程的压力梯度项====================================
				Face_direction face_dir_w = fW;
				Face_direction face_dir_e = fE;
				Geo_BC geoBC_w = m_bc.getGeoType(i, j, k, face_dir_w);
				Geo_BC geoBC_e = m_bc.getGeoType(i, j, k, face_dir_e);
				Phy_BC phyBC_w = m_bc.getPhyType(geoBC_w);
				Phy_BC phyBC_e = m_bc.getPhyType(geoBC_e);

				if (phyBC_w == Phy_BC::NONE && phyBC_e == Phy_BC::NONE) {
					p_w = m_pressure.p.getMat(i - 1, j, k);  //P_W
					p_e = m_pressure.p.getMat(i + 1, j, k);  //P_E
				}
				else { //考虑边界条件的情况
					//考虑W面
					if (phyBC_w == Phy_BC::INLET || phyBC_w == Phy_BC::WALL) {
						p_w = m_pressure.p.getMat(i, j, k);      //P_W
						p_e = m_pressure.p.getMat(i + 1, j, k);  //P_E
					}
					else if (phyBC_w == Phy_BC::OUTLET) {
						p_w = m_bc.getP_Outlet();                //P_W
						p_e = m_pressure.p.getMat(i + 1, j, k);  //P_E
					}

					//考虑E面
					if (phyBC_e == Phy_BC::INLET || phyBC_e == Phy_BC::WALL) {
						p_w = m_pressure.p.getMat(i - 1, j, k);  //P_W
						p_e = m_pressure.p.getMat(i, j, k);      //P_E
					}
					else if (phyBC_e == Phy_BC::OUTLET) {
						p_w = m_pressure.p.getMat(i - 1, j, k);  //P_W
						p_e = m_bc.getP_Outlet();                //P_E
					}
				}

				//设置压力梯度项的离散系数
				Su_0 = m_coefficient.coef_u.getCoeff(i, j, k, Su);
				Su_P = Su_0 + 0.5 * (p_w - p_e) * m_mesh.getDy() * m_mesh.getDz();  //计算公式为Su=Su0+0.5*(P_W-P_E)*dy*dz
				m_coefficient.coef_u.setCoeff(i, j, k, Su, Su_P);
				

				//v方向的动量方程的压力梯度项====================================
				face_dir_w = fS;
				face_dir_e = fN;
				geoBC_w = m_bc.getGeoType(i, j, k, face_dir_w);
				geoBC_e = m_bc.getGeoType(i, j, k, face_dir_e);
				phyBC_w = m_bc.getPhyType(geoBC_w);
				phyBC_e = m_bc.getPhyType(geoBC_e);

				if (phyBC_w == Phy_BC::NONE && phyBC_e == Phy_BC::NONE) {
					p_w = m_pressure.p.getMat(i, j - 1, k);  //P_S
					p_e = m_pressure.p.getMat(i, j + 1, k);  //P_N
				}
				else { //考虑边界条件的情况
					//考虑S面
					if (phyBC_w == Phy_BC::INLET || phyBC_w == Phy_BC::WALL) {
						p_w = m_pressure.p.getMat(i, j, k);      //P_S
						p_e = m_pressure.p.getMat(i, j + 1, k);  //P_N
					}
					else if (phyBC_w == Phy_BC::OUTLET) {
						p_w = m_bc.getP_Outlet();                //P_S
						p_e = m_pressure.p.getMat(i, j + 1, k);  //P_N
					}

					//考虑N面
					if (phyBC_e == Phy_BC::INLET || phyBC_e == Phy_BC::WALL) {
						p_w = m_pressure.p.getMat(i, j - 1, k);  //P_S
						p_e = m_pressure.p.getMat(i, j, k);      //P_N
					}
					else if (phyBC_e == Phy_BC::OUTLET) {
						p_w = m_pressure.p.getMat(i, j - 1, k);  //P_S
						p_e = m_bc.getP_Outlet();                //P_N
					}
				}

				//设置压力梯度项的离散系数
				Su_0 = m_coefficient.coef_v.getCoeff(i, j, k, Su);
				Su_P = Su_0 + 0.5 * (p_w - p_e) * m_mesh.getDx() * m_mesh.getDz();  //计算公式为Su=Su0+0.5*(P_S-P_N)*dx*dz
				m_coefficient.coef_v.setCoeff(i, j, k, Su, Su_P);


				//w方向的动量方程的压力梯度项====================================
				face_dir_w = fB;
				face_dir_e = fT;
				geoBC_w = m_bc.getGeoType(i, j, k, face_dir_w);
				geoBC_e = m_bc.getGeoType(i, j, k, face_dir_e);
				phyBC_w = m_bc.getPhyType(geoBC_w);
				phyBC_e = m_bc.getPhyType(geoBC_e);

				if (phyBC_w == Phy_BC::NONE && phyBC_e == Phy_BC::NONE) {
					p_w = m_pressure.p.getMat(i, j, k - 1);  //P_B
					p_e = m_pressure.p.getMat(i, j, k + 1);  //P_T
				}
				else { //考虑边界条件的情况
					//考虑B面
					if (phyBC_w == Phy_BC::INLET || phyBC_w == Phy_BC::WALL) {
						p_w = m_pressure.p.getMat(i, j, k);      //P_B
						p_e = m_pressure.p.getMat(i, j, k + 1);  //P_T
					}
					else if (phyBC_w == Phy_BC::OUTLET) {
						p_w = m_bc.getP_Outlet();                //P_B
						p_e = m_pressure.p.getMat(i, j, k + 1);  //P_T
					}

					//考虑T面
					if (phyBC_e == Phy_BC::INLET || phyBC_e == Phy_BC::WALL) {
						p_w = m_pressure.p.getMat(i, j, k - 1);  //P_B
						p_e = m_pressure.p.getMat(i, j, k);      //P_T
					}
					else if (phyBC_e == Phy_BC::OUTLET) {
						p_w = m_pressure.p.getMat(i, j, k - 1);  //P_B
						p_e = m_bc.getP_Outlet();                //P_T
					}
				}

				//设置压力梯度项的离散系数
				Su_0 = m_coefficient.coef_w.getCoeff(i, j, k, Su);
				Su_P = Su_0 + 0.5 * (p_w - p_e) * m_mesh.getDx() * m_mesh.getDy();  //计算公式为Su=Su0+0.5*(P_B-P_T)*dx*dy
				m_coefficient.coef_w.setCoeff(i, j, k, Su, Su_P);
			}
		}
	}
};

void SIMPLE_Solver::IterativeSolver_vel(size_t iter_nonlin, size_t nstep, Variable var, size_t fre_termiPrint, std::string file_name_) {
	Iter_method method = m_solverSetting.iter_method; //读取线性迭代方法，如高斯赛德尔、雅各比

	size_t iter_num;                          //定义变量，用于指定不同变量线性迭代次数
	if (var == Variable::U) {
		iter_num = m_solverSetting.numIter_u; //u动量方程
	}
	else if (var == Variable::V) {
		iter_num = m_solverSetting.numIter_v; //v动量方程
	}
	else if (var == Variable::W) {
		iter_num = m_solverSetting.numIter_w; //w动量方程
	}
	else {
		throw std::invalid_argument("The variables of the linear iterative equation can only be specified as U, V, W.");
	}

	//定义变量场的引用
	//以下为错误使用示例，引用对象初始化后就被绑定，无法更改指向
	//VelocityField& field_value = m_velocity.cu;     
	//VelocityField& field_value0 = m_velocity.cu0;
	//CoefMatrix& field_coeff = m_coefficient.coef_u;
	//double alpha = m_solverSetting.alpha_v; //速度松弛系数
	//if (var == Variable::U) {
	//	field_value = m_velocity.cu;        //u单元中心速度值
	//	field_value0 = m_velocity.cu0;      //旧的u单元中心速度值
	//	field_coeff = m_coefficient.coef_u; //u动量方程系数矩阵
	//}
	//else if (var == Variable::V) {
	//	field_value = m_velocity.cv;        //v单元中心速度值
	//	field_value0 = m_velocity.cv0;      //旧的v单元中心速度值
	//	field_coeff = m_coefficient.coef_v; //v动量方程系数矩阵
	//}
	//else if (var == Variable::W) {
	//	field_value = m_velocity.cw;        //w单元中心速度值
	//	field_value0 = m_velocity.cw0;      //旧的w单元中心速度值
	//	field_coeff = m_coefficient.coef_w; //w动量方程系数矩阵
	//}
	//else {
	//	throw std::invalid_argument("The variables of the field can only be specified as U, V, W.");
	//}
	double alpha = m_solverSetting.alpha_v; //速度松弛系数

	// 定义一个简单的结构体来返回多个引用
	struct FieldRefs {
		VelocityField& value;
		VelocityField& value0;
		CoefMatrix& coeff;
	};

	auto get_field_refs = [&](Variable var) -> FieldRefs {
		if (var == Variable::U) {
			return { m_velocity.cu, m_velocity.cu0, m_coefficient.coef_u };
		}
		else if (var == Variable::V) {
			return { m_velocity.cv, m_velocity.cv0, m_coefficient.coef_v };
		}
		else if (var == Variable::W) {
			return { m_velocity.cw, m_velocity.cw0, m_coefficient.coef_w };
		}
		else {
			throw std::invalid_argument("The variables of the field can only be specified as U, V, W.");
		}
	};
	auto refs = get_field_refs(var);
	VelocityField& field_value = refs.value;        //速度单元中心值
	VelocityField& field_value0 = refs.value0;      //旧的速度单元中心值
	//field_value0.initial(0.0);                      //初始化旧的速度场
	CoefMatrix& field_coeff = refs.coeff;           //动量方程系数矩阵

	//定义用于迭代的变量
	double value_w;    //左侧单元质心的变量值，可指代U\V\W
	double value_e;    //右侧单元质心的变量值
	double value_s;    //下侧单元质心的变量值
	double value_n;    //上侧单元质心的变量值
	double value_b;    //底侧单元质心的变量值
	double value_t;    //顶侧单元质心的变量值
	double value_p;    //中心单元质心的变量值
	double value_p0;   //旧的中心单元质心的变量值

	//定义用于残差计算的变量
	double norm;
	double norm2;
	double norm_max = 0.0;
	double res;


	if (method == Iter_method::JACOBI) { //雅各比迭代
		//总迭代循环

		for (size_t iter = 1; iter < iter_num + 1; ++iter) {

			//将当前物理场拷贝给旧的物理场
			field_value0.copyMat(field_value);

			//初始化范数，用于残差的计算
			norm = 0.0;

			//逐个单元迭代
			for (size_t i = 0; i < m_mesh.getNx(); ++i) {
				for (size_t j = 0; j < m_mesh.getNy(); ++j) {
					for (size_t k = 0; k < m_mesh.getNz(); ++k) {
						//W侧的单元物理量，可指代U\V\W
						if (i == 0) {
							value_w = 0.0;
						}
						else {
							value_w = field_value0.getMat(i - 1, j, k);
						}

						//E侧的单元物理量，可指代U\V\W
						if (i == m_mesh.getNx() - 1) {
							value_e = 0.0;
						}
						else {
							value_e = field_value0.getMat(i + 1, j, k);
						}

						//S侧的单元物理量，可指代U\V\W
						if (j == 0) {
							value_s = 0.0;
						}
						else {
							value_s = field_value0.getMat(i, j - 1, k);
						}

						//N侧的单元物理量，可指代U\V\W
						if (j == m_mesh.getNy() - 1) {
							value_n = 0.0;
						}
						else {
							value_n = field_value0.getMat(i, j + 1, k);
						}

						//B侧的单元物理量，可指代U\V\W
						if (k == 0) {
							value_b = 0.0;
						}
						else {
							value_b = field_value0.getMat(i, j, k - 1);
						}

						//T侧的单元物理量，可指代U\V\W
						if (k == m_mesh.getNz() - 1) {
							value_t = 0.0;
						}
						else {
							value_t = field_value0.getMat(i, j, k + 1);
						}

						//核心部分，计算公式为ΦP=(1/aP)*(aWΦW+aEΦE+aSΦS+aNΦN+aBΦB+aTΦT+Su)
						value_p = (1.0 / field_coeff.getCoeff(i, j, k, aP)) *
							     (field_coeff.getCoeff(i, j, k, aW) * value_w
								+ field_coeff.getCoeff(i, j, k, aE) * value_e
								+ field_coeff.getCoeff(i, j, k, aS) * value_s
								+ field_coeff.getCoeff(i, j, k, aN) * value_n
								+ field_coeff.getCoeff(i, j, k, aB) * value_b
								+ field_coeff.getCoeff(i, j, k, aT) * value_t
								+ field_coeff.getCoeff(i, j, k, Su));

						//ap检查
						double aP_coef = field_coeff.getCoeff(i, j, k, aP);  // 这里 field_coeff = m_coefficient.coef_p
						if (std::fabs(aP_coef) < 1e-12) {
							std::cout << "WARNING: coef_vel aP too small at ("
								<< i << "," << j << "," << k
								<< "), aP=" << aP_coef << std::endl;
							continue;  // 或者把 aP_coef 设成一个小正数再用
						}

						//读取旧的单元中心物理场
						value_p0 = field_value0.getMat(i, j, k);

						//更新单元中心物理场，并使用显式松弛
						value_p = value_p0 + alpha * (value_p - value_p0);
						field_value.setMat(i, j, k, value_p);
						norm += std::pow((alpha * (value_p - value_p0)), 2);
					}
				}
			}

			//计算过程的信息输出=======================================================
			
			//当前迭代步的二范数
			size_t num_cell = m_mesh.getNx() * m_mesh.getNy() * m_mesh.getNz(); //单元数量
			norm2 = std::pow(norm / num_cell, 0.5);                             //计算公式为x2=(sum(norm^2)/ncell)^0.5

			//迭代中最大的二范数
			norm_max = std::max(norm2, norm_max) + 1e-20; //加上极小数防止分母取零

			//相对残差
			res = norm2 / norm_max;

			//判断停止迭代以及输出结果的步骤
			if (res < m_converSetting.lin_tolerance_u || iter == iter_num - 1) {
				//根据非线性迭轮次判断是否输出计算过程残差信息
				if (iter_nonlin % fre_termiPrint == 0 || iter_nonlin == 1 || iter_nonlin == nstep) {
					//终端输出
					std::cout << "No.iter_nonlin:" << iter_nonlin
						<< "|No.iter_lin:" << iter
						<< "|norm2:" << norm2
						<< "|norm_max:" << norm_max
						<< "|res:" << res
						<< std::endl;
					//文件保存
					std::ofstream lin_result(file_name_, std::ios::app);
					if (lin_result.is_open()) {
						lin_result << iter_nonlin << " "
							<< iter << " "
							<< norm2 << " "
							<< norm_max << " "
							<< res << std::endl;
					}
					else {
						std::cerr << "Cannot open file:" << file_name_ << std::endl;
					}
				}
				break;
			}
		}	
	}

	else if (method == Iter_method::GAUSS_SEIDEL) {  //高斯赛德尔迭代
		//总迭代循环
		for (size_t iter = 1; iter < iter_num + 1; ++iter) {

			//初始化范数，用于残差的计算
			norm = 0.0;

			//逐个单元迭代
			for (size_t i = 0; i < m_mesh.getNx(); ++i) {
				for (size_t j = 0; j < m_mesh.getNy(); ++j) {
					for (size_t k = 0; k < m_mesh.getNz(); ++k) {
						//W侧的单元物理量，可指代U\V\W
						if (i == 0) {
							value_w = 0.0;
						}
						else {
							value_w = field_value.getMat(i - 1, j, k);
						}

						//E侧的单元物理量，可指代U\V\W
						if (i == m_mesh.getNx() - 1) {
							value_e = 0.0;
						}
						else {
							value_e = field_value.getMat(i + 1, j, k);
						}

						//S侧的单元物理量，可指代U\V\W
						if (j == 0) {
							value_s = 0.0;
						}
						else {
							value_s = field_value.getMat(i, j - 1, k);
						}

						//N侧的单元物理量，可指代U\V\W
						if (j == m_mesh.getNy() - 1) {
							value_n = 0.0;
						}
						else {
							value_n = field_value.getMat(i, j + 1, k);
						}

						//B侧的单元物理量，可指代U\V\W
						if (k == 0) {
							value_b = 0.0;
						}
						else {
							value_b = field_value.getMat(i, j, k - 1);
						}

						//T侧的单元物理量，可指代U\V\W
						if (k == m_mesh.getNz() - 1) {
							value_t = 0.0;
						}
						else {
							value_t = field_value.getMat(i, j, k + 1);
						}

						//核心部分，计算公式为ΦP=(1/aP)*(aWΦW+aEΦE+aSΦS+aNΦN+aBΦB+aTΦT+Su)
						value_p = (1.0 / field_coeff.getCoeff(i, j, k, aP)) *
								(field_coeff.getCoeff(i, j, k, aW) * value_w
								+ field_coeff.getCoeff(i, j, k, aE) * value_e
								+ field_coeff.getCoeff(i, j, k, aS) * value_s
								+ field_coeff.getCoeff(i, j, k, aN) * value_n
								+ field_coeff.getCoeff(i, j, k, aB) * value_b
								+ field_coeff.getCoeff(i, j, k, aT) * value_t
								+ field_coeff.getCoeff(i, j, k, Su));

						//ap检查
						double aP_coef = field_coeff.getCoeff(i, j, k, aP);  // 这里 field_coeff = m_coefficient.coef_p
						if (std::fabs(aP_coef) < 1e-12) {
							std::cout << "WARNING: coef_vel aP too small at ("
								<< i << "," << j << "," << k
								<< "), aP=" << aP_coef << std::endl;
							continue;  // 或者把 aP_coef 设成一个小正数再用
						}

						//读取旧的单元中心物理场
						value_p0 = field_value.getMat(i, j, k);

						//更新单元中心物理场，并使用显式松弛
						value_p = value_p0 + alpha * (value_p - value_p0);
						field_value.setMat(i, j, k, value_p);
						norm += std::pow((alpha * (value_p - value_p0)), 2);
					}
				}
			}

			//计算过程的信息输出=======================================================

			//当前迭代步的二范数
			size_t num_cell = m_mesh.getNx() * m_mesh.getNy() * m_mesh.getNz(); //单元数量
			norm2 = std::pow(norm / num_cell, 0.5);                             //计算公式为x2=(sum(norm^2)/ncell)^0.5

			//迭代中最大的二范数
			norm_max = std::max(norm2, norm_max) + 1e-20; //加上极小数防止分母取零

			//相对残差
			res = norm2 / norm_max;

			//判断停止迭代以及输出结果的步骤
			if (res < m_converSetting.lin_tolerance_u || iter == iter_num - 1) {
				//根据非线性迭轮次判断是否输出计算过程残差信息
				if (iter_nonlin % fre_termiPrint == 0 || iter_nonlin == 1 || iter_nonlin == nstep) {
					//终端输出
					std::cout << "No.iter_nonlin:" << iter_nonlin
						<< "|No.iter_lin:" << iter
						<< "|norm2:" << norm2
						<< "|norm_max:" << norm_max
						<< "|res:" << res
						<< std::endl;
					//文件保存
					std::ofstream lin_result(file_name_, std::ios::app);
					if (lin_result.is_open()) {
						lin_result << iter_nonlin << " "
							<< iter << " "
							<< norm2 << " "
							<< norm_max << " "
							<< res << std::endl;
					}
					else {
						std::cerr << "Cannot open file:" << file_name_ << std::endl;
					}
				}
				break;
			}
		}
	}

	else {
		throw std::invalid_argument("The linear iterator method can currently only use Jacobi and Gauss-Seidel.");
	}
};

void SIMPLE_Solver::IterativeSolver_p(size_t iter_nonlin, size_t nstep, size_t fre_termiPrint, std::string file_name_) {
	Iter_method method = m_solverSetting.iter_method; //读取线性迭代方法，如高斯赛德尔、雅各比

	size_t iter_num;                                  //定义变量，用于指定不同变量线性迭代次数
	iter_num = m_solverSetting.numIter_p;             //压力修正方程

	//定义变量场的引用
	PressureField& field_value = m_pressure.pp;       //单元中心的压力修正值
	field_value.initial(0.0);                         //初始化压力修正场（注意！！！）
	
	//临时定义一个矩阵，用于存储旧的压力修正值
	Field field_value0{ m_pressure.pp.getX(), m_pressure.pp.getY(), m_pressure.pp.getZ(), 0.0 };
	//PressureField& field_value0 = m_pressure.pp0;   //不直接使用类内成员变量的原因为下一步旧场置零，会导致收敛判断时不便于判定
	field_value0.initial(0.0);                        //旧的初始化压力修正场

	CoefMatrix& field_coeff = m_coefficient.coef_p;
	double alpha = m_solverSetting.alpha_p;           //压力松弛系数

	//定义用于迭代的变量
	double value_w;    //左侧单元质心的变量值，可指代P
	double value_e;    //右侧单元质心的变量值
	double value_s;    //下侧单元质心的变量值
	double value_n;    //上侧单元质心的变量值
	double value_b;    //底侧单元质心的变量值
	double value_t;    //顶侧单元质心的变量值
	double value_p;    //中心单元质心的变量值
	double value_p0;   //旧的中心单元质心的变量值

	//定义用于残差计算的变量
	double norm;
	double norm2;
	double norm_max = 0.0;
	double res;


	if (method == Iter_method::JACOBI) { //雅各比迭代
		//总迭代循环

		for (size_t iter = 1; iter < iter_num + 1; ++iter) {

			//将当前物理场拷贝给旧的物理场
			field_value0.copyMat(field_value);

			//初始化范数，用于残差的计算
			norm = 0.0;

			//逐个单元迭代
			for (size_t i = 0; i < m_mesh.getNx(); ++i) {
				for (size_t j = 0; j < m_mesh.getNy(); ++j) {
					for (size_t k = 0; k < m_mesh.getNz(); ++k) {
						//W侧的单元物理量，可指代P
						if (i == 0) {
							value_w = 0.0;
						}
						else {
							value_w = field_value0.getMat(i - 1, j, k);
						}

						//E侧的单元物理量，可指代P
						if (i == m_mesh.getNx() - 1) {
							value_e = 0.0;
						}
						else {
							value_e = field_value0.getMat(i + 1, j, k);
						}

						//S侧的单元物理量，可指代P
						if (j == 0) {
							value_s = 0.0;
						}
						else {
							value_s = field_value0.getMat(i, j - 1, k);
						}

						//N侧的单元物理量，可指代P
						if (j == m_mesh.getNy() - 1) {
							value_n = 0.0;
						}
						else {
							value_n = field_value0.getMat(i, j + 1, k);
						}

						//B侧的单元物理量，可指代P
						if (k == 0) {
							value_b = 0.0;
						}
						else {
							value_b = field_value0.getMat(i, j, k - 1);
						}

						//T侧的单元物理量，可指代P
						if (k == m_mesh.getNz() - 1) {
							value_t = 0.0;
						}
						else {
							value_t = field_value0.getMat(i, j, k + 1);
						}

						//核心部分，计算公式为ΦP=(1/aP)*(aWΦW+aEΦE+aSΦS+aNΦN+aBΦB+aTΦT+Su)
						value_p = (1.0 / field_coeff.getCoeff(i, j, k, aP)) *
								(field_coeff.getCoeff(i, j, k, aW) * value_w
								+ field_coeff.getCoeff(i, j, k, aE) * value_e
								+ field_coeff.getCoeff(i, j, k, aS) * value_s
								+ field_coeff.getCoeff(i, j, k, aN) * value_n
								+ field_coeff.getCoeff(i, j, k, aB) * value_b
								+ field_coeff.getCoeff(i, j, k, aT) * value_t
								+ field_coeff.getCoeff(i, j, k, Su));

						//ap检查
						//double aP_coef = field_coeff.getCoeff(i, j, k, aP);  // 这里 field_coeff = m_coefficient.coef_p
						//if (std::fabs(aP_coef) < 1e-12) {
						//	std::cout << "WARNING: coef_p aP too small at ("
						//		<< i << "," << j << "," << k
						//		<< "), aP=" << aP_coef << std::endl;
						//	continue;  // 或者把 aP_coef 设成一个小正数再用
						//}

						//读取旧的单元中心物理场
						value_p0 = field_value0.getMat(i, j, k);

						//更新单元中心物理场，并使用显式松弛
						value_p = value_p0 + alpha * (value_p - value_p0);
						field_value.setMat(i, j, k, value_p);
						norm += std::pow((alpha * (value_p - value_p0)), 2);
					}
				}
			}

			//计算过程的信息输出=======================================================

			//当前迭代步的二范数
			size_t num_cell = m_mesh.getNx() * m_mesh.getNy() * m_mesh.getNz(); //单元数量
			norm2 = std::pow(norm / num_cell, 0.5);                             //计算公式为x2=(sum(norm^2)/ncell)^0.5

			//迭代中最大的二范数
			norm_max = std::max(norm2, norm_max) + 1e-20; //加上极小数防止分母取零

			//相对残差
			res = norm2 / norm_max;

			//判断停止迭代以及输出结果的步骤
			if (res < m_converSetting.lin_tolerance_p || iter == iter_num - 1) {
				//根据非线性迭轮次判断是否输出计算过程残差信息
				if (iter_nonlin % fre_termiPrint == 0 || iter_nonlin == 1 || iter_nonlin == nstep) {
					//终端输出
					std::cout << "No.iter_nonlin:" << iter_nonlin
						<< "|No.iter_lin:" << iter
						<< "|norm2:" << norm2
						<< "|norm_max:" << norm_max
						<< "|res:" << res
						<< std::endl;
					//文件保存
					std::ofstream lin_result(file_name_, std::ios::app);
					if (lin_result.is_open()) {
						lin_result << iter_nonlin << " "
							<< iter << " "
							<< norm2 << " "
							<< norm_max << " "
							<< res << std::endl;
					}
					else {
						std::cerr << "Cannot open file:" << file_name_ << std::endl;
					}
				}
				break;
			}
		}
	}

	else if (method == Iter_method::GAUSS_SEIDEL) {  //高斯赛德尔迭代
		//总迭代循环
		for (size_t iter = 1; iter < iter_num + 1; ++iter) {

			//初始化范数，用于残差的计算
			norm = 0.0;

			//逐个单元迭代
			for (size_t i = 0; i < m_mesh.getNx(); ++i) {
				for (size_t j = 0; j < m_mesh.getNy(); ++j) {
					for (size_t k = 0; k < m_mesh.getNz(); ++k) {
						//W侧的单元物理量，可指代P
						if (i == 0) {
							value_w = 0.0;
						}
						else {
							value_w = field_value.getMat(i - 1, j, k);
						}

						//E侧的单元物理量，可指代P
						if (i == m_mesh.getNx() - 1) {
							value_e = 0.0;
						}
						else {
							value_e = field_value.getMat(i + 1, j, k);
						}

						//S侧的单元物理量，可指代P
						if (j == 0) {
							value_s = 0.0;
						}
						else {
							value_s = field_value.getMat(i, j - 1, k);
						}

						//N侧的单元物理量，可指代P
						if (j == m_mesh.getNy() - 1) {
							value_n = 0.0;
						}
						else {
							value_n = field_value.getMat(i, j + 1, k);
						}

						//B侧的单元物理量，可指代P
						if (k == 0) {
							value_b = 0.0;
						}
						else {
							value_b = field_value.getMat(i, j, k - 1);
						}

						//T侧的单元物理量，可指代P
						if (k == m_mesh.getNz() - 1) {
							value_t = 0.0;
						}
						else {
							value_t = field_value.getMat(i, j, k + 1);
						}

						//核心部分，计算公式为ΦP=(1/aP)*(aWΦW+aEΦE+aSΦS+aNΦN+aBΦB+aTΦT+Su)
						value_p = (1.0 / field_coeff.getCoeff(i, j, k, aP)) *
							(field_coeff.getCoeff(i, j, k, aW) * value_w
								+ field_coeff.getCoeff(i, j, k, aE) * value_e
								+ field_coeff.getCoeff(i, j, k, aS) * value_s
								+ field_coeff.getCoeff(i, j, k, aN) * value_n
								+ field_coeff.getCoeff(i, j, k, aB) * value_b
								+ field_coeff.getCoeff(i, j, k, aT) * value_t
								+ field_coeff.getCoeff(i, j, k, Su));

						//ap检查
						//double aP_coef = field_coeff.getCoeff(i, j, k, aP);  // 这里 field_coeff = m_coefficient.coef_p
						//if (std::fabs(aP_coef) < 1e-12) {
						//	std::cout << "WARNING: coef_p aP too small at ("
						//		<< i << "," << j << "," << k
						//		<< "), aP=" << aP_coef << std::endl;
						//	continue;  // 或者把 aP_coef 设成一个小正数再用
						//}

						//读取旧的单元中心物理场
						value_p0 = field_value.getMat(i, j, k);

						//更新单元中心物理场，并使用显式松弛
						value_p = value_p0 + alpha * (value_p - value_p0);
						field_value.setMat(i, j, k, value_p);
						norm += std::pow((alpha * (value_p - value_p0)), 2);
					}
				}
			}

			//计算过程的信息输出=======================================================

			//当前迭代步的二范数
			size_t num_cell = m_mesh.getNx() * m_mesh.getNy() * m_mesh.getNz(); //单元数量
			norm2 = std::pow(norm / num_cell, 0.5);                             //计算公式为x2=(sum(norm^2)/ncell)^0.5

			//迭代中最大的二范数
			norm_max = std::max(norm2, norm_max) + 1e-20; //加上极小数防止分母取零

			//相对残差
			res = norm2 / norm_max;

			//判断停止迭代以及输出结果的步骤
			if (res < m_converSetting.lin_tolerance_p || iter == iter_num - 1) {
				//根据非线性迭轮次判断是否输出计算过程残差信息
				if (iter_nonlin % fre_termiPrint == 0 || iter_nonlin == 1 || iter_nonlin == nstep) {
					//终端输出
					std::cout << "No.iter_nonlin:" << iter_nonlin
						<< "|No.iter_lin:" << iter
						<< "|norm2:" << norm2
						<< "|norm_max:" << norm_max
						<< "|res:" << res
						<< std::endl;
					//文件保存
					std::ofstream lin_result(file_name_, std::ios::app);
					if (lin_result.is_open()) {
						lin_result << iter_nonlin << " "
							<< iter << " "
							<< norm2 << " "
							<< norm_max << " "
							<< res << std::endl;
					}
					else {
						std::cerr << "Cannot open file:" << file_name_ << std::endl;
					}
				}
				break;
			}
		}
	}

	else {
		throw std::invalid_argument("The linear iterator method can currently only use Jacobi and Gauss-Seidel.");
	}
};

void SIMPLE_Solver::CalcuCoeff_momentEqn(Conv_scheme conv_scheme) {
	
	//求解对流扩散方程的离散系数(不考虑边界条件)
	CalcuCoeff_momentEqn_ConvAndDiff(conv_scheme);  

	//计算动量方程的离散系数――子函数2：考虑x方向动量方程的边界条件
	CalcuCoeff_momentEqn_BC_u();

	//计算动量方程的离散系数――子函数3：考虑y方向动量方程的边界条件
	CalcuCoeff_momentEqn_BC_v();

	//计算动量方程的离散系数――子函数4：考虑z方向动量方程的边界条件
	CalcuCoeff_momentEqn_BC_w();

	//计算动量方程的离散系数――子函数5：考虑压力梯度项
	CalcuCoeff_momentEqn_gradP();

};

void SIMPLE_Solver::CalcuCoeff_continuityEqn() {
	//逐单元循环
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				double rho_e;                                         //e界面处的密度
				double rho_w;                                         //w界面处的密度
				double rho_n;                                         //n界面处的密度
				double rho_s;                                         //s界面处的密度
				double rho_t;                                         //t界面处的密度
				double rho_b;                                         //b界面处的密度

				double area_x = m_mesh.getDy() * m_mesh.getDz();      //x方向界面的截面积
				double area_y = m_mesh.getDx() * m_mesh.getDz();      //y方向界面的截面积
				double area_z = m_mesh.getDy() * m_mesh.getDx();      //z方向界面的截面积

				//计算aE
				Face_direction face_dir = fE;                         //E面
				Coeffient_direction coeff_dir = aE;
				Geo_BC geoBC = m_bc.getGeoType(i, j, k, face_dir);    //读取相应几何边界
				double Coeff_aE;                                      //定义计算aE的对象
				Phy_BC phyBC = m_bc.getPhyType(geoBC);                //读取相应物理边界
				if (phyBC != Phy_BC::NONE || i == m_mesh.getNx() - 1) {
					rho_e = m_fluid.getDensity(i, j, k);
					Coeff_aE = 0.0;
					m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aE);
				}
				else {
					rho_e = 0.5 * (m_fluid.getDensity(i, j, k) + m_fluid.getDensity(i + 1, j, k));
					//计算公式为aE=0.5*rho*dx*dx*dy*dy*(1/aP+1/aE)
					Coeff_aE = 0.5 * rho_e * area_x * area_x
						* ( 1.0 / Coeff_U(i, j, k, aP) + 1.0 / Coeff_U(i + 1, j, k, aP));
					m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aE);
				}

				//计算aW
				face_dir = fW;                                        //W面
				coeff_dir = aW;
				geoBC = m_bc.getGeoType(i, j, k, face_dir);           //读取相应几何边界
				double Coeff_aW;                                      //定义计算aW的对象
				phyBC = m_bc.getPhyType(geoBC);                       //读取相应物理边界
				if (phyBC != Phy_BC::NONE || i == 0) {
					rho_w = m_fluid.getDensity(i, j, k);
					Coeff_aW = 0.0;
					m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aW);
				}
				else {
					rho_w = 0.5 * (m_fluid.getDensity(i, j, k) + m_fluid.getDensity(i - 1, j, k));
					Coeff_aW = 0.5 * rho_w * area_x * area_x
						* (1.0 / Coeff_U(i, j, k, aP) + 1.0 / Coeff_U(i - 1, j, k, aP));
					m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aW);
				}

				//计算aN
				face_dir = fN;                                        //N面
				coeff_dir = aN;
				geoBC = m_bc.getGeoType(i, j, k, face_dir);           //读取相应几何边界
				double Coeff_aN;                                      //定义计算aN的对象
				phyBC = m_bc.getPhyType(geoBC);                       //读取相应物理边界
				if (phyBC != Phy_BC::NONE || j == m_mesh.getNy() - 1) {
					rho_n = m_fluid.getDensity(i, j, k);
					Coeff_aN = 0.0;
					m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aN);
				}
				else {
					rho_n = 0.5 * (m_fluid.getDensity(i, j, k) + m_fluid.getDensity(i, j + 1, k));
					Coeff_aN = 0.5 * rho_n * area_y * area_y
						* (1.0 / Coeff_V(i, j, k, aP) + 1.0 / Coeff_V(i, j + 1, k, aP));
					m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aN);
				}

				//计算aS
				face_dir = fS;                                        //S面
				coeff_dir = aS;
				geoBC = m_bc.getGeoType(i, j, k, face_dir);           //读取相应几何边界
				double Coeff_aS;                                      //定义计算aS的对象
				phyBC = m_bc.getPhyType(geoBC);                       //读取相应物理边界
				if (phyBC != Phy_BC::NONE || j == 0) {
					rho_s = m_fluid.getDensity(i, j, k);
					Coeff_aS = 0.0;
					m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aS);
				}
				else {
					rho_s = 0.5 * (m_fluid.getDensity(i, j, k) + m_fluid.getDensity(i, j - 1, k));
					Coeff_aS = 0.5 * rho_s * area_y * area_y
						* (1.0 / Coeff_V(i, j, k, aP) + 1.0 / Coeff_V(i, j - 1, k, aP));
					m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aS);
				}

				//计算aT 
				face_dir = fT;                                        //T面
				coeff_dir = aT;
				geoBC = m_bc.getGeoType(i, j, k, face_dir);           //读取相应几何边界
				double Coeff_aT;                                      //定义计算aT的对象
				phyBC = m_bc.getPhyType(geoBC);                       //读取相应物理边界
				if (phyBC != Phy_BC::NONE || k == m_mesh.getNz() - 1) {
					rho_t = m_fluid.getDensity(i, j, k);
					Coeff_aT = 0.0;
					m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aT);
				}
				else {
					rho_t = 0.5 * (m_fluid.getDensity(i, j, k) + m_fluid.getDensity(i, j, k + 1));
					Coeff_aT = 0.5 * rho_t * area_z * area_z
						* (1.0 / Coeff_W(i, j, k, aP) + 1.0 / Coeff_W(i, j, k + 1, aP));
					m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aT);
				}

				//计算aB
				face_dir = fB;                                        //B面
				coeff_dir = aB;
				geoBC = m_bc.getGeoType(i, j, k, face_dir);           //读取相应几何边界
				double Coeff_aB;                                      //定义计算aB的对象
				phyBC = m_bc.getPhyType(geoBC);                       //读取相应物理边界
				if (phyBC != Phy_BC::NONE || k == 0) {
					rho_b = m_fluid.getDensity(i, j, k);
					Coeff_aB = 0.0;
					m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aB);
				}
				else {
					rho_b = 0.5 * (m_fluid.getDensity(i, j, k) + m_fluid.getDensity(i, j, k - 1));
					Coeff_aB = 0.5 * rho_b * area_z * area_z
						* (1.0 / Coeff_W(i, j, k, aP) + 1.0 / Coeff_W(i, j, k - 1, aP));
					m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aB);
				}

				//计算aP
				coeff_dir = aP;
				double Coeff_aP = (Coeff_aE + Coeff_aW + Coeff_aN + Coeff_aS + Coeff_aT + Coeff_aB);
				m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_aP);

				//计算Su
				coeff_dir = Su;
				double Coeff_Su;                                      //定义计算Su的对象
				Coeff_Su = -(rho_e * Uf(i + 1, j, k) * area_x - rho_w * Uf(i, j, k) * area_x
					+ rho_n * Vf(i, j + 1, k) * area_y - rho_s * Vf(i, j, k) * area_y
					+ rho_t * Wf(i, j, k + 1) * area_z - rho_b * Wf(i, j, k) * area_z);
				m_coefficient.coef_p.setCoeff(i, j, k, coeff_dir, Coeff_Su);
			}
		}
	}
};

void SIMPLE_Solver::update_Pressure() {
	double pressure;                          //定义变量，用于计算单元压力的更新后的值
	double alpha = m_solverSetting.alpha_p;   //松弛系数
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				pressure = P(i, j, k) + alpha * Pp(i, j, k);
				m_pressure.p.setMat(i, j, k, pressure);
			}
		}
	}
};

void SIMPLE_Solver::update_CenterVelocity() {
	//定义用于计算的压力修正变量
	double pp_L;                             //计算u时为P_W,计算v时为P_S,计算w时为P_B
	double pp_R;                             //计算u时为P_E,计算v时为P_N,计算w时为P_T

	//定义变量
	double vel;                              //用于计算单元速度更新后的值
	double alpha = m_solverSetting.alpha_v;  //松弛系数

	double area_x = m_mesh.getDy() * m_mesh.getDz();      //x方向界面的截面积
	double area_y = m_mesh.getDx() * m_mesh.getDz();      //y方向界面的截面积
	double area_z = m_mesh.getDy() * m_mesh.getDx();      //z方向界面的截面积
	
	//修正速度u
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				//根据几何位置设置压力
				if (i == 0) {                   //左侧面
					pp_L = Pp(i, j, k);
				}
				else {
					pp_L = Pp(i - 1, j, k);
				}

				if (i == m_mesh.getNx() - 1) {  //右侧面
					pp_R = Pp(i, j, k);
				}
				else {
					pp_R = Pp(i + 1, j, k);
				}

				//根据边界条件设置压力
				Face_direction face_dir = fW;                         //W面
				Coeffient_direction coeff_dir = aW;
				Geo_BC geoBC = m_bc.getGeoType(i, j, k, face_dir);    //读取相应几何边界
				Phy_BC phyBC = m_bc.getPhyType(geoBC);                //读取相应物理边界
				if (phyBC == Phy_BC::INLET || phyBC == Phy_BC::WALL) {
					pp_L = Pp(i, j, k);
					pp_R = Pp(i + 1, j, k);
				}
				else if (phyBC == Phy_BC::OUTLET) {
					pp_L = m_bc.getPP_Outlet();
					pp_R = Pp(i + 1, j, k);
				}

				face_dir = fE;                                 //E面
				coeff_dir = aE;
				geoBC = m_bc.getGeoType(i, j, k, face_dir);    //读取相应几何边界
				phyBC = m_bc.getPhyType(geoBC);                //读取相应物理边界
				if (phyBC == Phy_BC::INLET || phyBC == Phy_BC::WALL) {
					pp_L = Pp(i - 1, j, k);
					pp_R = Pp(i, j, k);
				}
				else if (phyBC == Phy_BC::OUTLET) {
					pp_L = Pp(i - 1, j, k);
					pp_R = m_bc.getPP_Outlet();
				}

				//计算并更新体心速度u
				//vel = Uc(i, j, k) + alpha * 0.5 * (pp_L - pp_R) * area_x / Coeff_U(i, j, k, aP);
				vel = Uc(i, j, k) + 0.5 * (pp_L - pp_R) * area_x / Coeff_U(i, j, k, aP);            //不使用松弛
				m_velocity.cu.setMat(i, j, k, vel);
			}
		}
	}

	//修正速度v
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				//根据几何位置设置压力
				if (j == 0) {                   //下侧面
					pp_L = Pp(i, j, k);
				}
				else {
					pp_L = Pp(i, j - 1, k);
				}

				if (j == m_mesh.getNy() - 1) {  //上侧面
					pp_R = Pp(i, j, k);
				}
				else {
					pp_R = Pp(i, j + 1, k);
				}

				//根据边界条件设置压力
				Face_direction face_dir = fS;                         //S面
				Coeffient_direction coeff_dir = aS;
				Geo_BC geoBC = m_bc.getGeoType(i, j, k, face_dir);    //读取相应几何边界
				Phy_BC phyBC = m_bc.getPhyType(geoBC);                //读取相应物理边界
				if (phyBC == Phy_BC::INLET || phyBC == Phy_BC::WALL) {
					pp_L = Pp(i, j, k);
					pp_R = Pp(i, j + 1, k);
				}
				else if (phyBC == Phy_BC::OUTLET) {
					pp_L = m_bc.getPP_Outlet();
					pp_R = Pp(i, j + 1, k);
				}

				face_dir = fN;                                 //N面
				coeff_dir = aN;
				geoBC = m_bc.getGeoType(i, j, k, face_dir);    //读取相应几何边界
				phyBC = m_bc.getPhyType(geoBC);                //读取相应物理边界
				if (phyBC == Phy_BC::INLET || phyBC == Phy_BC::WALL) {
					pp_L = Pp(i, j - 1, k);
					pp_R = Pp(i, j, k);
				}
				else if (phyBC == Phy_BC::OUTLET) {
					pp_L = Pp(i, j - 1, k);
					pp_R = m_bc.getPP_Outlet();
				}

				//计算并更新体心速度v
				//vel = Vc(i, j, k) + alpha * 0.5 * (pp_L - pp_R) * area_y / Coeff_V(i, j, k, aP);
				vel = Vc(i, j, k) + 0.5 * (pp_L - pp_R) * area_y / Coeff_V(i, j, k, aP);			//不使用松弛
				m_velocity.cv.setMat(i, j, k, vel);
			}
		}
	}

	//修正速度w
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				//根据几何位置设置压力
				if (k == 0) {                   //底侧面
					pp_L = Pp(i, j, k);
				}
				else {
					pp_L = Pp(i, j, k - 1);
				}

				if (k == m_mesh.getNz() - 1) {  //顶侧面
					pp_R = Pp(i, j, k);
				}
				else {
					pp_R = Pp(i, j, k + 1);
				}

				//根据边界条件设置压力
				Face_direction face_dir = fB;                         //B面
				Coeffient_direction coeff_dir = aB;
				Geo_BC geoBC = m_bc.getGeoType(i, j, k, face_dir);    //读取相应几何边界
				Phy_BC phyBC = m_bc.getPhyType(geoBC);                //读取相应物理边界
				if (phyBC == Phy_BC::INLET || phyBC == Phy_BC::WALL) {
					pp_L = Pp(i, j, k);
					pp_R = Pp(i, j, k + 1);
				}
				else if (phyBC == Phy_BC::OUTLET) {
					pp_L = m_bc.getPP_Outlet();
					pp_R = Pp(i, j, k + 1);
				}

				face_dir = fT;                                 //T面
				coeff_dir = aT;
				geoBC = m_bc.getGeoType(i, j, k, face_dir);    //读取相应几何边界
				phyBC = m_bc.getPhyType(geoBC);                //读取相应物理边界
				if (phyBC == Phy_BC::INLET || phyBC == Phy_BC::WALL) {
					pp_L = Pp(i, j, k - 1);
					pp_R = Pp(i, j, k);
				}
				else if (phyBC == Phy_BC::OUTLET) {
					pp_L = Pp(i, j, k - 1);
					pp_R = m_bc.getPP_Outlet();
				}

				//计算并更新体心速度w
				//vel = Wc(i, j, k) + alpha * 0.5 * (pp_L - pp_R) * area_z / Coeff_W(i, j, k, aP);
				vel = Wc(i, j, k) + 0.5 * (pp_L - pp_R) * area_z / Coeff_W(i, j, k, aP);			//不使用松弛
				m_velocity.cw.setMat(i, j, k, vel);
			}
		}
	}
};

void SIMPLE_Solver::update_FaceVelocity() {
	//定义变量
	double vel;                              //用于计算单元速度更新后的值
	double alpha = m_solverSetting.alpha_v;  //松弛系数

	double area_x = m_mesh.getDy() * m_mesh.getDz();      //x方向界面的截面积
	double area_y = m_mesh.getDx() * m_mesh.getDz();      //y方向界面的截面积
	double area_z = m_mesh.getDy() * m_mesh.getDx();      //z方向界面的截面积

	//更新内部面的uf
	for (size_t i = 1; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				//计算uf
				/*vel = Uf(i, j, k)
					+ alpha * 0.5
					* (1.0 / Coeff_U(i - 1, j, k, aP) + 1.0 / Coeff_U(i, j, k, aP))
					* (Pp(i - 1, j, k) - Pp(i, j, k))
					* area_x;*/
				vel = Uf(i, j, k)      //不使用松弛
					+ 0.5
					* (1.0 / Coeff_U(i - 1, j, k, aP) + 1.0 / Coeff_U(i, j, k, aP))
					* (Pp(i - 1, j, k) - Pp(i, j, k))
					* area_x;
				//更新uf
				m_velocity.uf.setMat(i, j, k, vel);
			}
		}
	}

	//更新内部面的vf
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 1; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				//计算vf
				/*vel = Vf(i, j, k)
					+ alpha * 0.5
					* (1.0 / Coeff_V(i, j - 1, k, aP) + 1.0 / Coeff_V(i, j, k, aP))
					* (Pp(i, j - 1, k) - Pp(i, j, k))
					* area_y;*/
				vel = Vf(i, j, k)       //不使用松弛
					+ 0.5
					* (1.0 / Coeff_V(i, j - 1, k, aP) + 1.0 / Coeff_V(i, j, k, aP))
					* (Pp(i, j - 1, k) - Pp(i, j, k))
					* area_y;
				//更新vf
				m_velocity.vf.setMat(i, j, k, vel);
			}
		}
	}

	//更新内部面的wf
	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 1; k < m_mesh.getNz(); ++k) {
				//计算wf
				/*vel = Wf(i, j, k)    
					+ alpha * 0.5
					* (1.0 / Coeff_W(i, j, k - 1, aP) + 1.0 / Coeff_W(i, j, k, aP))
					* (Pp(i, j, k - 1) - Pp(i, j, k))
					* area_z;*/
				vel = Wf(i, j, k)
					+ 0.5              //不使用松弛
					* (1.0 / Coeff_W(i, j, k - 1, aP) + 1.0 / Coeff_W(i, j, k, aP))
					* (Pp(i, j, k - 1) - Pp(i, j, k))
					* area_z;
				//更新wf
				m_velocity.wf.setMat(i, j, k, vel);
			}
		}
	}
};

void SIMPLE_Solver::saveVTK_Vel_Result(std::string filename) {
	std::ofstream vtk_fid(filename);
	
	size_t nx = m_mesh.getNx() + 1;
	size_t ny = m_mesh.getNy() + 1;
	size_t nz = m_mesh.getNz() + 1;
	size_t ncell = (nx - 1) * (ny - 1) * (nz - 1);

	// 写入标头
	vtk_fid << "# vtk DataFile Version 3.0\n"
		<< "flash 3d grid and solution\n"
		<< "ASCII\n"
		<< "DATASET RECTILINEAR_GRID\n";

	// 写入网格信息
	vtk_fid << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n"
		<< "X_COORDINATES " << nx << " float\n";

	for (size_t i = 0; i < nx; ++i) { vtk_fid << m_mesh.getCellPointX(i) << " "; }
	vtk_fid << "\n";

	vtk_fid << "Y_COORDINATES " << ny << " float\n";
	for (size_t j = 0; j < ny; ++j) { vtk_fid << m_mesh.getCellPointY(j) << " "; }
	vtk_fid << "\n";

	vtk_fid << "Z_COORDINATES " << nz << " float\n";
	for (size_t k = 0; k < nz; ++k) { vtk_fid << m_mesh.getCellPointZ(k) << " "; }
	vtk_fid << "\n";

	// 写入单元格数据
	vtk_fid << "CELL_DATA " << ncell << "\n"
		<< "FIELD FieldData 3\n";

	// 写入速度数据
	// 写入速度u
	vtk_fid << "u" << " 1 " << ncell << " float\n";
	for (size_t k = 0; k < nz - 1; ++k) {
		for (size_t j = 0; j < ny - 1; ++j) {
			for (size_t i = 0; i < nx - 1; ++i) {
				vtk_fid << Uc(i, j, k) << " ";
			}
		}
	}

	// 写入速度v
	vtk_fid << "v" << " 1 " << ncell << " float\n";
	for (size_t k = 0; k < nz - 1; ++k) {
		for (size_t j = 0; j < ny - 1; ++j) {
			for (size_t i = 0; i < nx - 1; ++i) {
				vtk_fid << Vc(i, j, k) << " ";
			}
		}
	}

	// 写入速度w
	vtk_fid << "w" << " 1 " << ncell << " float\n";
	for (size_t k = 0; k < nz - 1; ++k) {
		for (size_t j = 0; j < ny - 1; ++j) {
			for (size_t i = 0; i < nx - 1; ++i) {
				vtk_fid << Wc(i, j, k) << " ";
			}
		}
	}

	vtk_fid.close();
};   

void SIMPLE_Solver::saveVTK_Pressure_Result(std::string filename) {
	std::ofstream vtk_fid(filename);

	size_t nx = m_mesh.getNx() + 1;
	size_t ny = m_mesh.getNy() + 1;
	size_t nz = m_mesh.getNz() + 1;
	size_t ncell = (nx - 1) * (ny - 1) * (nz - 1);

	// 写入标头
	vtk_fid << "# vtk DataFile Version 3.0\n"
		<< "flash 3d grid and solution\n"
		<< "ASCII\n"
		<< "DATASET RECTILINEAR_GRID\n";

	// 写入网格信息
	vtk_fid << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n"
		<< "X_COORDINATES " << nx << " float\n";

	for (size_t i = 0; i < nx; ++i) { vtk_fid << m_mesh.getCellPointX(i) << " "; }
	vtk_fid << "\n";

	vtk_fid << "Y_COORDINATES " << ny << " float\n";
	for (size_t j = 0; j < ny; ++j) { vtk_fid << m_mesh.getCellPointY(j) << " "; }
	vtk_fid << "\n";

	vtk_fid << "Z_COORDINATES " << nz << " float\n";
	for (size_t k = 0; k < nz; ++k) { vtk_fid << m_mesh.getCellPointZ(k) << " "; }
	vtk_fid << "\n";

	// 写入单元格数据
	vtk_fid << "CELL_DATA " << ncell << "\n"
		<< "FIELD FieldData 2\n";

	// 写入压力数据
	vtk_fid << "p" << " 1 " << ncell << " float\n";
	for (size_t k = 0; k < nz - 1; ++k) {
		for (size_t j = 0; j < ny - 1; ++j) {
			for (size_t i = 0; i < nx - 1; ++i) {
				vtk_fid << P(i, j, k) << " ";
			}
		}
	}

	// 写入压力修正数据
	vtk_fid << "pp" << " 1 " << ncell << " float\n";
	for (size_t k = 0; k < nz - 1; ++k) {
		for (size_t j = 0; j < ny - 1; ++j) {
			for (size_t i = 0; i < nx - 1; ++i) {
				vtk_fid << Pp(i, j, k) << " ";
			}
		}
	}

	vtk_fid.close();

};

void SIMPLE_Solver::calculate_U_norm2() {
	//定义变量
	double error{ 0 };       //用于计算各单元更新物理量与旧物理量的差值
	double norm2_U{ 0 };     //用于计算速度u的二范数
	double max_norm2_U{ 0 }; //用于记录目前最大的二范数

	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				error = Uc(i, j, k) - Uc0(i, j, k);
				norm2_U += error * error;
			}
		}
	}

	norm2_U = std::sqrt(norm2_U / (m_mesh.getNx() * m_mesh.getNy() * m_mesh.getNz()));
	
	//目前最大的二范数
	max_norm2_U = std::max(norm2_U, m_postProcess.get_MaxNorm2_u());

	//设定二范数
	m_postProcess.set_Norm2_u(norm2_U);
	m_postProcess.set_MaxNorm2_u(max_norm2_U);
};

void SIMPLE_Solver::calculate_V_norm2() {
	//定义变量
	double error{ 0 };       //用于计算各单元更新物理量与旧物理量的差值
	double norm2_V{ 0 };     //用于计算速度v的二范数
	double max_norm2_V{ 0 }; //用于记录目前最大的二范数

	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				error = Vc(i, j, k) - Vc0(i, j, k);
				norm2_V += error * error;
			}
		}
	}

	norm2_V = std::sqrt(norm2_V / (m_mesh.getNx() * m_mesh.getNy() * m_mesh.getNz()));

	//目前最大的二范数
	max_norm2_V = std::max(norm2_V, m_postProcess.get_MaxNorm2_v());

	//设定二范数
	m_postProcess.set_Norm2_v(norm2_V);
	m_postProcess.set_MaxNorm2_v(max_norm2_V);
};

void SIMPLE_Solver::calculate_W_norm2() {
	//定义变量
	double error{ 0 };       //用于计算各单元更新物理量与旧物理量的差值
	double norm2_W{ 0 };     //用于计算速度w的二范数
	double max_norm2_W{ 0 }; //用于记录目前最大的二范数

	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				error = Wc(i, j, k) - Wc0(i, j, k);
				norm2_W += error * error;
			}
		}
	}

	norm2_W = std::sqrt(norm2_W / (m_mesh.getNx() * m_mesh.getNy() * m_mesh.getNz()));

	//目前最大的二范数
	max_norm2_W = std::max(norm2_W, m_postProcess.get_MaxNorm2_w());
	
	//设定二范数
	m_postProcess.set_Norm2_w(norm2_W);
	m_postProcess.set_MaxNorm2_w(max_norm2_W);
};

void SIMPLE_Solver::calculate_P_norm2() {
	//定义变量
	double error{ 0 };       //用于计算各单元更新物理量与旧物理量的差值
	double norm2_P{ 0 };     //用于计算压力p的二范数

	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				error = P(i, j, k) - P0(i, j, k);
				norm2_P += error * error;
			}
		}
	}

	norm2_P = std::sqrt(norm2_P / (m_mesh.getNx() * m_mesh.getNy() * m_mesh.getNz()));

	//目前最大的二范数
	double max_norm2_P = std::max(norm2_P, m_postProcess.get_MaxNorm2_p());

	//设定二范数
	m_postProcess.set_Norm2_p(norm2_P);
	m_postProcess.set_MaxNorm2_p(max_norm2_P);
};

void SIMPLE_Solver::calculate_Pp_norm2() {
	//定义变量
	double error{ 0 };       //用于计算各单元更新物理量与旧物理量的差值
	double norm2_Pp{ 0 };    //用于计算压力修正值pp的二范数

	for (size_t i = 0; i < m_mesh.getNx(); ++i) {
		for (size_t j = 0; j < m_mesh.getNy(); ++j) {
			for (size_t k = 0; k < m_mesh.getNz(); ++k) {
				error = Pp(i, j, k) - Pp0(i, j, k);
				norm2_Pp += error * error;
			}
		}
	}

	norm2_Pp = std::sqrt(norm2_Pp / (m_mesh.getNx() * m_mesh.getNy() * m_mesh.getNz()));

	//目前最大的二范数
	double max_norm2_Pp = std::max(norm2_Pp, m_postProcess.get_MaxNorm2_pp());

	//设定二范数
	m_postProcess.set_Norm2_pp(norm2_Pp);
	m_postProcess.set_MaxNorm2_pp(max_norm2_Pp);
};

void SIMPLE_Solver::calculate_residual() {
	//计算残差
	calculate_U_norm2();      //计算速度u的二阶范数
	calculate_V_norm2();      //计算速度v的二阶范数
	calculate_W_norm2();      //计算速度w的二阶范数
	calculate_P_norm2();	  //计算压力p的二阶范数
	calculate_Pp_norm2();	  //计算压力pp的二阶范数
};

void SIMPLE_Solver::check_convergence(size_t num_iter) {
	//读取非线性收敛容差
	double nonlin_tol = m_converSetting.nonlin_tolerance;

	if ((L2_U() / max_L2_U() < nonlin_tol &&
		L2_V() / max_L2_V() < nonlin_tol &&
		L2_W() / max_L2_W() < nonlin_tol &&
		L2_P() / max_L2_P() < nonlin_tol) ||
		L2_Pp() / max_L2_Pp() < nonlin_tol) {
		m_solverSetting.stop_switch = true;    //达到收敛条件，停止迭代

		//输出收敛信息
		std::cout << " " << std::endl;
		std::cout << "******************************************" << std::endl;
		std::cout << "Non-linear convergence achieved!" << std::endl;
		std::cout << "******************************************" << std::endl;
		std::cout << " " << std::endl;
		std::cout << "Total iteration number: " << num_iter << std::endl;
		std::cout << "L2_U / max_L2_U: " << L2_U() / max_L2_U() << std::endl;
		std::cout << "L2_V / max_L2_V: " << L2_V() / max_L2_V() << std::endl;
		std::cout << "L2_W / max_L2_W: " << L2_W() / max_L2_W() << std::endl;
		std::cout << "L2_P / max_L2_P: " << L2_P() / max_L2_P() << std::endl;
		std::cout << "L2_Pp / max_L2_Pp: " << L2_Pp() / max_L2_Pp() << std::endl;
		std::cout << " " << std::endl;
	}
};

bool SIMPLE_Solver::check_stop() {
	return m_solverSetting.stop_switch;
};
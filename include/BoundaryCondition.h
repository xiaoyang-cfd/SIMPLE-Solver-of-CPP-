#pragma once
#include "Field.h"
#include "typedef.h"
#include "ConstructMesh.h"
#include <unordered_map>


//存储速度数值边界条件的类
class BC_vel {
public:
	BC_vel();
	
	//设置速度
	void setU(double u1);
	void setV(double v1);
	void setW(double w1);

	//设置速度通量
	void setFluxU(double flux_u1);
	void setFluxV(double flux_v1);
	void setFluxW(double flux_w1);

	//设置边界类型
	void setTypeU(Num_BC type_u1);
	void setTypeV(Num_BC type_v1);
	void setTypeW(Num_BC type_w1);

	//返回速度
	double getU();
	double getV();
	double getW();

	//返回速度通量
	double getFluxU();
	double getFluxV();
	double getFluxW();

	//返回边界类型
	Num_BC getTypeU();
	Num_BC getTypeV();
	Num_BC getTypeW();

private:
	//边界上的速度值
	double u;
	double v;
	double w;

	//边界上的速度通量值
	double flux_u;
	double flux_v;
	double flux_w;

	//边界上的类型,包括内部边none、第一类边界constant、第二类边界flux
	Num_BC type_u;
	Num_BC type_v;
	Num_BC type_w;
};

typedef std::vector<std::vector<BC_vel>> vVectorVelBC;         //速度数值边界类型

class BcField : public Field {
public:
	BcField() = default;

	//构造函数
	BcField(size_t nx, size_t ny, size_t nz, const ConstructMesh& mesh, double value = 0, double p_outlet_ = 0, double pp_outlet_ = 0);
	
	//打印特定边界类型
	//a=1时打印几何边界,a=2时打印物理边界
	void printType(size_t nx, size_t ny, size_t nz, int a = 0);

	//设置特定单元特定面的几何边界类型
	void setGeoType(size_t x, size_t y, size_t z, Face_direction dir, Geo_BC bc);

	//设置几何边界与物理边界的对应关系
	void setPhyBC(Phy_BC xmin_phy, Phy_BC xmax_phy, Phy_BC ymin_phy, Phy_BC ymax_phy, Phy_BC zmin_phy, Phy_BC zmax_phy);

	//设置特定单元特定面的速度边界类型
	void setUBC(size_t x, size_t y, size_t z, Face_direction dir, BC_vel bc);  //设定速度u边界
	void setVBC(size_t x, size_t y, size_t z, Face_direction dir, BC_vel bc);  //设定速度v边界
	void setWBC(size_t x, size_t y, size_t z, Face_direction dir, BC_vel bc);  //设定速度w边界

	//设置所有单元速度数值边界条件
	void setVelTotalBC(std::vector<Num_BC> xmin_bc, std::vector<double> xmin_vel, 
				  std::vector<Num_BC> xmax_bc, std::vector<double> xmax_vel,
		          std::vector<Num_BC> ymin_bc, std::vector<double> ymin_vel,
		          std::vector<Num_BC> ymax_bc, std::vector<double> ymax_vel,
		          std::vector<Num_BC> zmin_bc, std::vector<double> zmin_vel,
		          std::vector<Num_BC> zmax_bc, std::vector<double> zmax_vel );

	//返回特定单元特定面的几何边界类型
	Geo_BC getGeoType(size_t x, size_t y, size_t z, Face_direction dir);

	//返回特定几何边界的物理边界类型
	Phy_BC getPhyType(Geo_BC bc);

	//返回特定单元特定面的速度边界数值
	const double& getBC_U(size_t x, size_t y, size_t z, Face_direction dir);  //速度U
	const double& getBC_V(size_t x, size_t y, size_t z, Face_direction dir);  //速度V
	const double& getBC_W(size_t x, size_t y, size_t z, Face_direction dir);  //速度W

	//返回压力出口的数值
	const double& getP_Outlet();

	//返回压力修正值出口的数值
	const double& getPP_Outlet();

	//返回特定单元特定面的特定边界类型，参数说明：坐标x、y、z，面的方向dir
	//Phy_BC getPhyType(size_t x, size_t y, size_t z, Face_direction dir);

private:
	size_t num_dir;      //每个单元的边的个数
	vVectorGeoBC geo_bc; //几何边界条件，参见Geo_BC类
	//vVectorPhyBC phy_bc; 
	std::unordered_map<Geo_BC, Phy_BC> phy_bc;  //物理边界条件，参见Phy_BC类
	vVectorVelBC vel_bc; //速度数值边界条件，参见BC_cel类
	double p_outlet;     //压力出口边界条件
	double pp_outlet;    //压力修正值出口边界条件
};






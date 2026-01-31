#include"BoundaryCondition.h"

BcField::BcField(size_t nx, size_t ny, size_t nz, const ConstructMesh& mesh, double value, double p_outlet_, double pp_outlet_)
	: Field(nx, ny, nz, value), num_dir(6),
	geo_bc(nx * ny * nz, std::vector<Geo_BC>(num_dir, Geo_BC::NONE)),
	vel_bc(nx * ny * nz, std::vector<BC_vel>(num_dir, BC_vel())),
	p_outlet(p_outlet_), pp_outlet(pp_outlet_){
	phy_bc[Geo_BC::NONE] = Phy_BC::NONE;

	//设定几何的边界条件
	double eps = 1e-12;
	for (size_t i = 0; i < nx; ++i) {
		for (size_t j = 0; j < ny; ++j) {
			for (size_t k = 0; k < nz; ++k) {
				//分别对x、y、z方向判别
				//x方向
				if (std::abs(mesh.getCellPointX(i) - mesh.getMin(X)) < eps) {
					setGeoType(i, j, k, fW, Geo_BC::XMIN);
				}
				if (std::abs(mesh.getCellPointX(i + 1) - mesh.getMax(X)) < eps) {
					setGeoType(i, j, k, fE, Geo_BC::XMAX);
				}

				//y方向
				if (std::abs(mesh.getCellPointY(j) - mesh.getMin(Y)) < eps) {
					setGeoType(i, j, k, fS, Geo_BC::YMIN);
				}
				if (std::abs(mesh.getCellPointY(j + 1) - mesh.getMax(Y)) < eps) {
					setGeoType(i, j, k, fN, Geo_BC::YMAX);
				}

				//z方向
				if (std::abs(mesh.getCellPointZ(k) - mesh.getMin(Z)) < eps) {
					setGeoType(i, j, k, fB, Geo_BC::ZMIN);
				}
				if (std::abs(mesh.getCellPointZ(k + 1) - mesh.getMax(Z)) < eps) {
					setGeoType(i, j, k, fT, Geo_BC::ZMAX);
				}
			}
		}
	}
}

void BcField::setGeoType(size_t x, size_t y, size_t z, Face_direction dir, Geo_BC bc) {
	size_t index = (z * getY() * getX()) + (y * getX()) + x;

	if (x >= getX() || y >= getY() || z >= getZ()) {
		throw std::out_of_range("The given element index is out of range");
	}

	geo_bc[index][dir] = bc;
};

void BcField::setPhyBC(Phy_BC xmin_phy, Phy_BC xmax_phy, Phy_BC ymin_phy, Phy_BC ymax_phy, Phy_BC zmin_phy, Phy_BC zmax_phy) {
	phy_bc[Geo_BC::XMIN] = xmin_phy;
	phy_bc[Geo_BC::XMAX] = xmax_phy;
	phy_bc[Geo_BC::YMIN] = ymin_phy;
	phy_bc[Geo_BC::YMAX] = ymax_phy;
	phy_bc[Geo_BC::ZMIN] = zmin_phy;
	phy_bc[Geo_BC::ZMAX] = zmax_phy;
}

void BcField::setUBC(size_t x, size_t y, size_t z, Face_direction dir, BC_vel bc) {
	size_t index = (z * getY() * getX()) + (y * getX()) + x;

	if (x >= getX() || y >= getY() || z >= getZ()) {
		throw std::out_of_range("The given element index is out of range");
	}

	//将bc的边界赋值给特定单元的特定面
	vel_bc[index][dir].setU(bc.getU());

	vel_bc[index][dir].setFluxU(bc.getFluxU());

	vel_bc[index][dir].setTypeU(bc.getTypeU());
};

void BcField::setVBC(size_t x, size_t y, size_t z, Face_direction dir, BC_vel bc) {
	size_t index = (z * getY() * getX()) + (y * getX()) + x;

	if (x >= getX() || y >= getY() || z >= getZ()) {
		throw std::out_of_range("The given element index is out of range");
	}

	//将bc的边界赋值给特定单元的特定面
	vel_bc[index][dir].setV(bc.getV());

	vel_bc[index][dir].setFluxV(bc.getFluxV());

	vel_bc[index][dir].setTypeV(bc.getTypeV());
};

void BcField::setWBC(size_t x, size_t y, size_t z, Face_direction dir, BC_vel bc) {
	size_t index = (z * getY() * getX()) + (y * getX()) + x;

	if (x >= getX() || y >= getY() || z >= getZ()) {
		throw std::out_of_range("The given element index is out of range");
	}

	//将bc的边界赋值给特定单元的特定面
	vel_bc[index][dir].setW(bc.getW());

	vel_bc[index][dir].setFluxW(bc.getFluxW());

	vel_bc[index][dir].setTypeW(bc.getTypeW());
};

void BcField::setVelTotalBC(std::vector<Num_BC> xmin_bc, std::vector<double> xmin_vel,
							std::vector<Num_BC> xmax_bc, std::vector<double> xmax_vel,
							std::vector<Num_BC> ymin_bc, std::vector<double> ymin_vel,
							std::vector<Num_BC> ymax_bc, std::vector<double> ymax_vel,
							std::vector<Num_BC> zmin_bc, std::vector<double> zmin_vel,
							std::vector<Num_BC> zmax_bc, std::vector<double> zmax_vel) {
	for (size_t i = 0; i < this->getX(); ++i){
		for (size_t j = 0; j < this->getY(); ++j) {
			for (size_t k = 0; k < this->getZ();++k) {
				//循环遍历每个单元
				//归纳总体功能：将类型为constant的边界，该面的速度改为定值，类型改为constant;
				//将类型为flux的边界，该面的速度通量改为定值，类型改为flux;

				//设置xmax面速度边界
				Face_direction dir = fE;
				if (getGeoType(i, j, k, dir) == Geo_BC::XMAX) {
					//设置u速度
					Num_BC type = xmax_bc[0];
					double vel = xmax_vel[0];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeU(type);   
						velBC.setU(vel);
						this->setUBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeU(type);
						velBC.setFluxU(vel);
						this->setUBC(i, j, k, dir, velBC);
					}

					//设置v速度
					type = xmax_bc[1];
					vel = xmax_vel[1];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeV(type);
						velBC.setV(vel);
						this->setVBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeV(type);
						velBC.setFluxV(vel);
						this->setVBC(i, j, k, dir, velBC);
					}

					//设置w速度
					type = xmax_bc[2];
					vel = xmax_vel[2];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeW(type);
						velBC.setW(vel);
						this->setWBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeW(type);
						velBC.setFluxW(vel);
						this->setWBC(i, j, k, dir, velBC);
					}
				}


				//设置xmin面速度边界
				dir = fW;
				if (getGeoType(i, j, k, dir) == Geo_BC::XMIN) {
					//设置u速度
					Num_BC type = xmin_bc[0];
					double vel = xmin_vel[0];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeU(type);
						velBC.setU(vel);
						this->setUBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeU(type);
						velBC.setFluxU(vel);
						this->setUBC(i, j, k, dir, velBC);
					}

					//设置v速度
					type = xmin_bc[1];
					vel = xmin_vel[1];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeV(type);
						velBC.setV(vel);
						this->setVBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeV(type);
						velBC.setFluxV(vel);
						this->setVBC(i, j, k, dir, velBC);
					}

					//设置w速度
					type = xmin_bc[2];
					vel = xmin_vel[2];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeW(type);
						velBC.setW(vel);
						this->setWBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeW(type);
						velBC.setFluxW(vel);
						this->setWBC(i, j, k, dir, velBC);
					}
				}

				//设置ymax面速度边界
				dir = fN;
				if (getGeoType(i, j, k, dir) == Geo_BC::YMAX) {
					//设置u速度
					Num_BC type = ymax_bc[0];
					double vel = ymax_vel[0];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeU(type);
						velBC.setU(vel);
						this->setUBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeU(type);
						velBC.setFluxU(vel);
						this->setUBC(i, j, k, dir, velBC);
					}

					//设置v速度
					type = ymax_bc[1];
					vel = ymax_vel[1];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeV(type);
						velBC.setV(vel);
						this->setVBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeV(type);
						velBC.setFluxV(vel);
						this->setVBC(i, j, k, dir, velBC);
					}

					//设置w速度
					type = ymax_bc[2];
					vel = ymax_vel[2];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeW(type);
						velBC.setW(vel);
						this->setWBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeW(type);
						velBC.setFluxW(vel);
						this->setWBC(i, j, k, dir, velBC);
					}
				}

				//设置ymin面速度边界
				dir = fS;
				if (getGeoType(i, j, k, dir) == Geo_BC::YMIN) {
					//设置u速度
					Num_BC type = ymin_bc[0];
					double vel = ymin_vel[0];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeU(type);
						velBC.setU(vel);
						this->setUBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeU(type);
						velBC.setFluxU(vel);
						this->setUBC(i, j, k, dir, velBC);
					}
					
					//设置v速度
					type = ymin_bc[1];
					vel = ymin_vel[1];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeV(type);
						velBC.setV(vel);
						this->setVBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeV(type);
						velBC.setFluxV(vel);
						this->setVBC(i, j, k, dir, velBC);
					}

					//设置w速度
					type = ymin_bc[2];
					vel = ymin_vel[2];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeW(type);
						velBC.setW(vel);
						this->setWBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeW(type);
						velBC.setFluxW(vel);
						this->setWBC(i, j, k, dir, velBC);
					}
				}

				//设置zmax面速度边界
				dir = fT;
				if (getGeoType(i, j, k, dir) == Geo_BC::ZMAX) {
					//设置u速度
					Num_BC type = zmax_bc[0];
					double vel = zmax_vel[0];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeU(type);
						velBC.setU(vel);
						this->setUBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeU(type);
						velBC.setFluxU(vel);
						this->setUBC(i, j, k, dir, velBC);
					}

					//设置v速度
					type = zmax_bc[1];
					vel = zmax_vel[1];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeV(type);
						velBC.setV(vel);
						this->setVBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeV(type);
						velBC.setFluxV(vel);
						this->setVBC(i, j, k, dir, velBC);
					}

					//设置w速度
					type = zmax_bc[2];
					vel = zmax_vel[2];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeW(type);
						velBC.setW(vel);
						this->setWBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeW(type);
						velBC.setFluxW(vel);
						this->setWBC(i, j, k, dir, velBC);
					}
				}

				//设置zmin面速度边界
				dir = fB;
				if (getGeoType(i, j, k, dir) == Geo_BC::ZMIN) {
					//设置u速度
					Num_BC type = zmin_bc[0];
					double vel = zmin_vel[0];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeU(type);
						velBC.setU(vel);
						this->setUBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeU(type);
						velBC.setFluxU(vel);
						this->setUBC(i, j, k, dir, velBC);
					}

					//设置v速度
					type = zmin_bc[1];
					vel = zmin_vel[1];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeV(type);
						velBC.setV(vel);
						this->setVBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeV(type);
						velBC.setFluxV(vel);
						this->setVBC(i, j, k, dir, velBC);
					}

					//设置w速度
					type = zmin_bc[2];
					vel = zmin_vel[2];
					if (type == Num_BC::CONSTANT) {
						//设置对应面的速度边界
						BC_vel velBC;
						velBC.setTypeW(type);
						velBC.setW(vel);
						this->setWBC(i, j, k, dir, velBC);
					}
					if (type == Num_BC::FLUX) {
						BC_vel velBC;
						velBC.setTypeW(type);
						velBC.setFluxW(vel);
						this->setWBC(i, j, k, dir, velBC);
					}
				}
			}
		}
	}
};

Geo_BC BcField::getGeoType(size_t x, size_t y, size_t z, Face_direction dir) {
	size_t index = (z * getY() * getX()) + (y * getX()) + x;

	if (x >= getX() || y >= getY() || z >= getZ()) {
		throw std::out_of_range("The given element index is out of range");
	}

	return geo_bc[index][dir];
}

Phy_BC BcField::getPhyType(Geo_BC bc) {
	auto it = phy_bc.find(bc);
	if (it != phy_bc.end()) {
		return it->second;
	}
	throw std::invalid_argument("Boundary not found");
};

const double& BcField::getBC_U(size_t x, size_t y, size_t z, Face_direction dir) {
	size_t index = (z * getY() * getX()) + (y * getX()) + x;

	if (x >= getX() || y >= getY() || z >= getZ()) {
		throw std::out_of_range("The given element index is out of range");
	}

	return vel_bc[index][dir].getU();
};

const double& BcField::getBC_V(size_t x, size_t y, size_t z, Face_direction dir) {
	size_t index = (z * getY() * getX()) + (y * getX()) + x;

	if (x >= getX() || y >= getY() || z >= getZ()) {
		throw std::out_of_range("The given element index is out of range");
	}

	return vel_bc[index][dir].getV();
};

const double& BcField::getBC_W(size_t x, size_t y, size_t z, Face_direction dir) {
	size_t index = (z * getY() * getX()) + (y * getX()) + x;

	if (x >= getX() || y >= getY() || z >= getZ()) {
		throw std::out_of_range("The given element index is out of range");
	}

	return vel_bc[index][dir].getW();
};

const double& BcField::getP_Outlet() {
	return p_outlet;
};

const double& BcField::getPP_Outlet() {
	return pp_outlet;
};

//Phy_BC BcField::getPhyType(size_t x, size_t y, size_t z, Face_direction dir) {
//	size_t index = (z * getY() * getX()) + (y * getX()) + x;
//
//	if (x >= getX() || y >= getY() || z >= getZ()) {
//		throw std::out_of_range("The given element index is out of range");
//	}
//
//	return phy_bc[index][dir];
//};

void BcField::printType(size_t nx, size_t ny, size_t nz, int a) {
	if (a == 1) {
		std::cout << "BCid: [ " << nx << "," << ny << "," << nz << " ]";
		std::cout << "\ngeometry_BC:";
		std::cout << "[ E:" << static_cast<int>(getGeoType(nx, ny, nz, fE));
		std::cout << " W:" << static_cast<int>(getGeoType(nx, ny, nz, fW));
		std::cout << " N:" << static_cast<int>(getGeoType(nx, ny, nz, fN));
		std::cout << " S:" << static_cast<int>(getGeoType(nx, ny, nz, fS));
		std::cout << " T:" << static_cast<int>(getGeoType(nx, ny, nz, fT));
		std::cout << " B:" << static_cast<int>(getGeoType(nx, ny, nz, fB));
		std::cout << " ]"  << std::endl;
	}
	/*if (a == 2) {
		std::cout << "BCid: [ " << nx << "," << ny << "," << nz << " ]";
		std::cout << "\nphyics_BC:";
		std::cout << "[ E:" << static_cast<int>(getPhyType(nx, ny, nz, fE));
		std::cout << " W:" << static_cast<int>(getPhyType(nx, ny, nz, fW));
		std::cout << " N:" << static_cast<int>(getPhyType(nx, ny, nz, fN));
		std::cout << " S:" << static_cast<int>(getPhyType(nx, ny, nz, fS));
		std::cout << " T:" << static_cast<int>(getPhyType(nx, ny, nz, fT));
		std::cout << " B:" << static_cast<int>(getPhyType(nx, ny, nz, fB));
		std::cout << " ]" << std::endl;
	}*/
};

BC_vel::BC_vel(): u(0), v(0), w(0), flux_u(0), flux_v(0), flux_w(0),
	type_u(Num_BC::NONE), type_v(Num_BC::NONE), type_w(Num_BC::NONE){}

void BC_vel::setU(double u1) {
	u = u1;
}
void BC_vel::setV(double v1) {
	v = v1;
}
void BC_vel::setW(double w1) {
	w = w1;
}

void BC_vel::setFluxU(double flux_u1) {
	flux_u = flux_u1;
}
void BC_vel::setFluxV(double flux_v1) {
	flux_v = flux_v1;
}
void BC_vel::setFluxW(double flux_w1) {
	flux_w = flux_w1;
}

void BC_vel::setTypeU(Num_BC type_u1) {
	type_u = type_u1;
}
void BC_vel::setTypeV(Num_BC type_v1) {
	type_v = type_v1;
}
void BC_vel::setTypeW(Num_BC type_w1) {
	type_w = type_w1;
}

double BC_vel::getU() {
	return u;
}
double BC_vel::getV() {
	return v;
}
double BC_vel::getW() {
	return w;
}

double BC_vel::getFluxU() {
	return flux_u;
}
double BC_vel::getFluxV() {
	return flux_v;
}
double BC_vel::getFluxW() {
	return flux_w;
}

Num_BC BC_vel::getTypeU() {
	return type_u;
}
Num_BC BC_vel::getTypeV() {
	return type_v;
}
Num_BC BC_vel::getTypeW() {
	return type_w;
}
#pragma once
#include<vector>

//特定类型别名
typedef std::vector<double> vectorDouble;


//记录特定含义枚举类
enum Direction { X = 0, Y = 1, Z = 2 }; //用于区分x、y、z方向

enum Coeffient_direction { aP = 0, aE = 1, aW = 2, aN = 3, aS = 4, Su = 5, aT = 6, aB = 7 }; //存储离散方程系数

enum class Conv_scheme { CD, UPWIND, HYBRID, POWERLAW};        //对流项离散格式

enum class Iter_method { JACOBI, GAUSS_SEIDEL};                //线性迭代方法

enum Face_direction { fE = 0, fW = 1, fN = 2, fS = 3, fT = 4, fB = 5}; //特定单元上各面的方位

enum class Geo_BC { NONE, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX}; //几何边界类型
typedef std::vector<std::vector<Geo_BC>> vVectorGeoBC;         //用于存放各单元各边几何边界类型的类

enum class Phy_BC { NONE, INLET, OUTLET, WALL };               //物理边界类型
typedef std::vector<std::vector<Phy_BC>> vVectorPhyBC;         //用于存放各单元各边物理边界类型的类

enum class Num_BC { NONE, CONSTANT, FLUX };                    //数值边界类型

enum class Variable { U, V, W, P };                            //求解的变量类型


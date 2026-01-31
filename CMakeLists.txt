cmake_minimum_required(VERSION 3.10)                        # 指定构建此项目所需的最低CMake版本
project(SIMPLE_Solver VERSION 1.0.0 LANGUAGES CXX)          # 定义项目的基本信息如下：
                                                            # 项目名称SIMPLE_Solver，后续可通过${PROJECT_NAME}引用
                                                            # VERSION 1.0.0：项目版本号
                                                            # LANGUAGES CXX：项目使用C++语言

# 指定C++语言标准
set(CMAKE_CXX_STANDARD 17)                                  # 要求使用C++17标准
set(CMAKE_CXX_STANDARD_REQUIRED ON)                         # 设置为ON表示必须使用C++17，如果编译器不支持会报错

# 直接设置输出到项目根目录的lib和bin文件夹，其中${CMAKE_SOURCE_DIR}为CMakeLists.txt所在的根目录路径
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin) # 可执行文件（.exe等）的输出目录
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin) # 动态库（.dll/.so）的输出目录
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib) # 静态库（.lib/.a）的输出目录

# 确保目录存在(可选)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

# 包含头文件
include_directories(include)                                # 添加头文件的搜索路径

# 添加src子目录
add_subdirectory(src)                                       # CMake会进入src目录，并执行其中的CMakeLists.txt文件

# 创建可执行文件
add_executable(${PROJECT_NAME} main.cpp)                    # 可执行文件配置，可执行文件的名称会定义为${PROJECT_NAME}.exe

# 链接库 - 使用正确的目标名
target_link_libraries(${PROJECT_NAME} Solver_source)        # 将${PROJECT_NAME}.exe可执行文件与名为Solver_source的库链接

# 设置可执行文件的包含目录
target_include_directories(${PROJECT_NAME} PRIVATE include) # 设置头文件的包含目录，即在include文件夹下搜索头文件
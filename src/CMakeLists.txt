# 自动读取该文件夹下的所有.cpp文件
file(GLOB LIB_SOURCES "*.cpp")

# 添加所有源文件
add_library(Solver_source SHARED ${LIB_SOURCES})                 # 创建动态库，库的目标名称为Solver_source

# 重要：设置Windows下的库属性
set_target_properties(Solver_source PROPERTIES
        # Windows下需要显式设置输出名称
        OUTPUT_NAME "Solver_source"
        # 移除lib前缀（Windows通常不需要）
        PREFIX ""
        # 确保生成导入库（.lib）
        WINDOWS_EXPORT_ALL_SYMBOLS ON                            # 自动导出所有符号，即在生成动态库（.dll）的同时，还额外生成导入库（.lib），供链接使用
)

# 设置包含目录
target_include_directories(Solver_source PUBLIC ../include)      # 为库目标设置头文件的搜索路径，PUBLIC表示这个目录对库本身和使用这个库的目标都可见

# 设置输出目录
set_target_properties(Solver_source PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib         # 设置静态库/导入库（.lib）的输出目录
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin         # 设置动态库（.dll）的输出目录
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin         # 设置可执行文件的输出目录
)
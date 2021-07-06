cmake_minimum_required(VERSION 3.0)

# Build the matrix library
include_directories(matrix)

# Build the dkm library
include_directories(dkm/include)

# Build the split library
include_directories(split)
add_library(split split/split.cpp)

# Build the config library
include_directories(config)
add_library(config config/Config.cpp)

# Build the mathfunc library
include_directories(mathfunc)
add_library(mathfunc mathfunc/MathFunc.cpp)

# Build the pca library
include_directories(pca)
add_library(pca pca/PCA.cpp)

# Build the threadpool library
include_directories(threadpool)
add_library(threadpool threadpool/ThreadPool.cpp)

# Build the mydefine library
include_directories(mydefine)
add_library(mydefine mydefine/MyDefine.cpp)
target_link_libraries(mydefine config mathfunc hetcaller inputparser threadpool)

# Build the inputparser library
include_directories(inputparser)
add_library(inputparser inputparser/InputParser.cpp)
target_link_libraries(inputparser mydefine)

# Build the hetcaller library
include_directories(hetcaller)
add_library(hetcaller hetcaller/HetCaller.cpp)
target_link_libraries(hetcaller mydefine split pca)

cmake_minimum_required(VERSION 3.0)
project(SCClone)

SET(CMAKE_BUILD_TYPE Release)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++11")

add_subdirectory(lib)

include_directories (
	${SCClone_SOURCE_DIR}/lib/config
	${SCClone_SOURCE_DIR}/lib/dkm/include
	${SCClone_SOURCE_DIR}/lib/hetcaller
	${SCClone_SOURCE_DIR}/lib/inputparser
	${SCClone_SOURCE_DIR}/lib/mathfunc
	${SCClone_SOURCE_DIR}/lib/matrix
	${SCClone_SOURCE_DIR}/lib/mydefine
	${SCClone_SOURCE_DIR}/lib/pca
	${SCClone_SOURCE_DIR}/lib/split
	${SCClone_SOURCE_DIR}/lib/threadpool
)

add_executable(scclone src/scclone.cpp)

target_link_libraries(scclone mydefine)

SET_TARGET_PROPERTIES(scclone PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SCClone_BINARY_DIR}/bin)

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
	set_property(TARGET scclone PROPERTY COMPILE_OPTIONS "-pthread")
	set_property(TARGET scclone PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
	target_link_libraries(scclone "${CMAKE_THREAD_LIBS_INIT}")
endif()

INSTALL(TARGETS scclone
	RUNTIME DESTINATION bin
)

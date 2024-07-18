file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/harmonicResistance/*.hpp")

add_library(harmonicResistance DynFullLaplacianInverseSolver.cpp StGreedy.cpp RobustnessGreedy.cpp ${HEADER_LIST})

target_include_directories(
	harmonicResistance
	PUBLIC ../include
)

target_link_libraries(harmonicResistance PUBLIC OpenMP::OpenMP_CXX PRIVATE networkit)

source_group(
	TREE "${PROJECT_SOURCE_DIR}/include"
	PREFIX "Header Files"
	FILES ${HEADER_LIST})

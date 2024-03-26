project(Viewer)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES Main ViewVolume ViewSurface Face FaceContainer Point3D Element 
	Hexahedron Prism Pyramid PyramidBase5 Tetrahedron)

add_executable(viewer ${SOURCE_FILES})

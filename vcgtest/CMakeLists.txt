cmake_minimum_required(VERSION 3.13)
project ( vcgtest )

#set(CMAKE_CXX_FLAGS "-lboost_program_options")
include_directories( ${CMAKE_SOURCE_DIR}/src "${CMAKE_SOURCE_DIR}/src/include" "${CMAKE_SOURCE_DIR}/src/CycleOptimization" "${CMAKE_SOURCE_DIR}/src/CycleOptimization/Annotation" "${CMAKE_SOURCE_DIR}/ann/") 
link_directories( ${CMAKE_SOURCE_DIR}/src/lib) 
# ADD library directory for ANN lib
link_directories( ${CMAKE_SOURCE_DIR}/ANN/lib)

# Boost
find_package( Boost 1.40 COMPONENTS program_options system REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories(${Boost_LIBRARY_DIRS})

# CGAL
find_package( CGAL REQUIRED )
include( ${CGAL_USE_FILE} )

# Usa FetchContent per scaricare la libreria ANN
include(FetchContent)

FetchContent_Declare(
    ann
    GIT_REPOSITORY https://github.com/michaelquigley/ann.git
    GIT_TAG master  # Puoi specificare un tag specifico se necessario
)

# Scarica e costruisci la libreria ANN
FetchContent_MakeAvailable(ann)

# Usa FetchContent per scaricare la libreria VCG
FetchContent_Declare(
    vcglib
    GIT_REPOSITORY https://github.com/cnr-isti-vclab/vcglib.git
    GIT_TAG devel  # Puoi specificare un tag specifico se necessario
)

# Scarica e costruisci la libreria VCG
FetchContent_MakeAvailable(vcglib)

add_executable( vcgtest 
				main.cpp 
				../src/CycleOptimization/Annotation/Annotated_polyhedron_3.h
				../src/CycleOptimization/Annotation/Annotated_polyhedron_items_3.h
				../src/CycleOptimization/Annotation/AnnotationComputation.h
				../src/CycleOptimization/Annotation/AnnotationComputation.cpp
				../src/CycleOptimization/Annotation/HalfedgeDS_annotated_face_base.h
				../src/CycleOptimization/Annotation/HalfedgeDS_annotated_halfedge_base.h
				../src/CycleOptimization/Annotation/HalfedgeDS_annotated_vertex_base.h
				../src/CycleOptimization/Annotation/Polyhedron_annotator_3.h
				../src/CycleOptimization/ShortCycles.cpp
				../src/CycleOptimization/SimpleGraph.h
				../src/CycleOptimization/canonical_loops.h 
				../src/CycleOptimization/canonical_loops.cpp
				../src/CycleOptimization/DijkstraAlgorithm.h 
				../src/CycleOptimization/DijkstraAlgorithm.cpp
				../src/CycleOptimization/edge_annotations_gauss.h 
				../src/CycleOptimization/edge_annotations_gauss.cpp
				../src/CycleOptimization/FibonacciHeap.h 
				../src/CycleOptimization/FibonacciHeap.cpp
				../src/AdjListGraph.cpp
				../src/AdjListGraph.h
				../src/ANNSearch.h
				../src/ANNSearch.cpp
				../src/FilesOutputForOptimalCycles.cpp
				../src/FilesOutputForOptimalCycles.h
				../src/Graph_max_tree.h
				../src/IntersectionOnMesh.h
				../src/IntersectionOnMesh.cpp
				../src/InverseLinkNumberMatrix.h
				../src/InverseLinkNumberMatrix.cpp
				../src/LinkNumberComputing.cpp
				../src/LinkNumberComputing.h
				../src/MapLoopsBackToMeshLevelSetAndArc.cpp
				../src/MapLoopsBackToMeshLevelSetAndArc.h
				../src/NonOverlappedLevelCycleAndArc.cpp
				../src/NonOverlappedLevelCycleAndArc.h
				../src/NonOverlappingCycles.cpp
				../src/NonOverlappingCycles.h
				../src/psbmReebGraph.cpp
				../src/psbmReebGraph.h
				../src/psbmReebGraphElements.h
				../src/ReebGraphPairing.cpp
				../src/ReebGraphPairing.h
				../src/SegmentIntersection.h
				../src/SegmentIntersection.cpp
				../src/SegmentPairIntersection.cpp
				../src/SimpleMesh.cpp
				../src/SimpleMesh.h 
				../src/UniqueDirection.cpp
				../src/UniqueDirection.h 
				../src/RenderVector3.h
				../build/_deps/vcglib-src/wrap/ply/plylib.cpp
				)

target_link_libraries(
	vcgtest
	PUBLIC
		${Boost_LIBRARIES} 
		boost_system 
		boost_program_options 
		ann
		vcglib  
	)
include_directories(../build/_deps/ann-src)
include_directories(../build/_deps/vcglib-src)
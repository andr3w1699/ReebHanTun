#set(CMAKE_CXX_FLAGS "-lboost_program_options")
include_directories( ${CMAKE_SOURCE_DIR}/src "${CMAKE_SOURCE_DIR}/src/include" "${CMAKE_SOURCE_DIR}/src/CycleOptimization" "${CMAKE_SOURCE_DIR}/src/CycleOptimization/Annotation" "${CMAKE_SOURCE_DIR}/ann/") 
# ADD library directory for ANN lib
link_directories(${CMAKE_BINARY_DIR}/_deps/ann-build)

set (CMAKE_CXX_STANDARD 17)


# Boost
find_package( Boost 1.40 COMPONENTS program_options system REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories(${Boost_LIBRARY_DIRS})

# CGAL
find_package( CGAL REQUIRED )
include_directories( ${CGAL_INCLUDE_DIRS} )
link_directories(${CGAL_LIBRARY_DIRS})
# include( ${CGAL_USE_FILE} )
message(STATUS "CGAL_INCLUDE_DIRS: ${CGAL_INCLUDE_DIRS}")

include_directories(${CMAKE_BINARY_DIR}/_deps/ann-src)

add_executable( ReebHanTun
				ComputeReebGraph.cpp
				CycleOptimization/Annotation/Annotated_polyhedron_3.h
				CycleOptimization/Annotation/Annotated_polyhedron_items_3.h
				CycleOptimization/Annotation/AnnotationComputation.h
				CycleOptimization/Annotation/AnnotationComputation.cpp
				CycleOptimization/Annotation/HalfedgeDS_annotated_face_base.h
				CycleOptimization/Annotation/HalfedgeDS_annotated_halfedge_base.h
				CycleOptimization/Annotation/HalfedgeDS_annotated_vertex_base.h
				CycleOptimization/Annotation/Polyhedron_annotator_3.h
				CycleOptimization/ShortCycles.cpp
				CycleOptimization/SimpleGraph.h
				CycleOptimization/canonical_loops.h 
				CycleOptimization/canonical_loops.cpp
				CycleOptimization/DijkstraAlgorithm.h 
				CycleOptimization/DijkstraAlgorithm.cpp
				CycleOptimization/edge_annotations_gauss.h 
				CycleOptimization/edge_annotations_gauss.cpp
				CycleOptimization/FibonacciHeap.h 
				CycleOptimization/FibonacciHeap.cpp
				AdjListGraph.cpp
				AdjListGraph.h
				ANNSearch.h
				ANNSearch.cpp
				FilesOutputForOptimalCycles.cpp
				FilesOutputForOptimalCycles.h
				Graph_max_tree.h
				IntersectionOnMesh.h
				IntersectionOnMesh.cpp
				InverseLinkNumberMatrix.h
				InverseLinkNumberMatrix.cpp
				LinkNumberComputing.cpp
				LinkNumberComputing.h
				MapLoopsBackToMeshLevelSetAndArc.cpp
				MapLoopsBackToMeshLevelSetAndArc.h
				NonOverlappedLevelCycleAndArc.cpp
				NonOverlappedLevelCycleAndArc.h
				NonOverlappingCycles.cpp
				NonOverlappingCycles.h
				psbmReebGraph.cpp
				psbmReebGraph.h
				psbmReebGraphElements.h
				ReebGraphPairing.cpp
				ReebGraphPairing.h
				SegmentIntersection.h
				SegmentIntersection.cpp
				SegmentPairIntersection.cpp
				SimpleMesh.cpp
				SimpleMesh.h 
				UniqueDirection.cpp
				UniqueDirection.h 
				RenderVector3.h
				)
target_link_libraries( ReebHanTun ${Boost_LIBRARIES} boost_system boost_program_options ann CGAL) #ANN.lib libANN.a

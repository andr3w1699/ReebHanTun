/******* HEADER INCLUSION *******/
#include<stdio.h>
#include<vcg/complex/complex.h>
#include <vcg/space/point3.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/bounding.h>
// input output
#include <wrap/io_trimesh/import.h>    // reads mesh files (e.g. .off, .ply formats)
#include<wrap/io_trimesh/export_off.h> // exports a mesh in .off format 
#include<wrap/io_trimesh/export_ply.h> // exports a mesh in .ply format
// mesh ReebHanTun
#include<SimpleMesh.h> //  Defines the _SimpleMesh structure used in ReebHanTun.

// include ComputeReebGraph.cpp
//

// horribile hack to avoid redefinition of main and allowing the use of the functions defined inside
// ComputeReebGraph.cpp

#define main __main
#include "../src/ComputeReebGraph.cpp"
#undef main

#include <time.h>
#include <sstream>
//#include <boost/progress.hpp>
#include <boost/timer/progress_display.hpp>
#include <boost/program_options.hpp>

#include "reebhantun_wrapper.h"
/******* END HEADER INCLUSION *******/

// Using VCG and Standard Namespace
// This avoids needing vcg:: and std:: prefixes for their functions 
using namespace vcg;
using namespace std;

/**** Defining Mesh Types ****/

// Declares three custom mesh components: vertices, edges, and faces.
class MyVertex;
class MyEdge;
class MyFace;

// MyUsedTypes defines the types of elements used in the mesh:
// It specifies that MyVertex, MyEdge, and MyFace will be used.
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>::AsVertexType, vcg::Use<MyEdge>::AsEdgeType,   vcg::Use<MyFace>::AsFaceType>{};

// Defines MyVertex, MyEdge, and MyFace as VCG-compatible mesh elements. Each class inherits from a corresponding
// VCG class (vcg::Vertex, vcg::Edge, vcg::Face) and includes various properties:
// MyVertex: Coord3f --> Stores 3D coordinates, Normal3f --> Stores normal vectors, Color4b --> Stores color (RGBA).
// BitFlags --> Stores extra data (e.g., selection state), VEAdj  --> Stores adjacency information between vertices and edges
// MyEdge: VertexRef: References vertices.
// MyFace: VertexRef: References vertices, Normal3f: Stores normal vectors, Color4b: Stores color (RGBA)
// BitFlags: Stores extra data (e.g., selection state), EFAdj: Stores adjacency information between edges and faces.
class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b, vcg::vertex::BitFlags, vcg::vertex::VEAdj >{};
class MyEdge : public vcg::Edge<MyUsedTypes,vcg::edge::VertexRef> {};
class MyFace  : public vcg::Face < MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::Color4b, vcg::face::BitFlags, vcg::face::EFAdj > {};

// MyMesh class represents the entire mesh. It inherits from vcg::tri::TriMesh and uses std::vector 
// to store collections of MyVertex, MyEdge, and MyFace elements. 
class MyMesh  : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyEdge>, std::vector<MyFace> > {};


/* 
MeshConverter is function to convert a VCG mesh type MyMesh in a ReebHanTun mesh _SimpleMesh. 
It also calculates the bounding box (minBd, maxBd), normals (meshNormal), and triangle orientations (OrientTriangles). 
The fEnlargeFactor is used to scale the vertices.
*/
void MeshConverter (_SimpleMeshVertex &minBd, _SimpleMeshVertex &maxBd, const MyMesh & vcg_mesh,  _SimpleMesh & rht_mesh, std::vector<Vector3> &meshNormal, std::vector<int> &OrientTriangles, const float fEnlargeFactor ) {
// map to store edges as pair of vertex indices
std::map<std::pair<int, int>, int, myPairCompare> edgeMapping;
// vertex conversion 
// reserves space in the OrientTriangles vector to avoid reallocations.
OrientTriangles.reserve(vcg_mesh.face.size());
// resize vecVertex of _SimpleMesh to store as many vertex as Mymesh
rht_mesh.vecVertex.reserve(vcg_mesh.vert.size());
// vertex conversion loop 
// This loop iterates over all vertices in vcg_mesh, 
// scales their coordinates by fEnlargeFactor, and adds them to rht_mesh.vecVertex.
for(size_t i = 0; i < vcg_mesh.vert.size(); ++i) {
	_SimpleMeshVertex tmpVer;
    const MyVertex &v = vcg_mesh.vert[i];
	tmpVer.x = v.P().X() * fEnlargeFactor;
    tmpVer.y = v.P().Y() * fEnlargeFactor;
	tmpVer.z = v.P().Z() * fEnlargeFactor;
	rht_mesh.vecVertex.push_back(tmpVer);
    // Updates Bounding Box based on the current vertex  
    if (rht_mesh.vecVertex.size() == 1) {
       minBd =  tmpVer;
       maxBd = tmpVer;
    } else {
            minBd.x = tmpVer.x < minBd.x ? tmpVer.x : minBd.x;
            minBd.y = tmpVer.y < minBd.y ? tmpVer.y : minBd.y;
            minBd.z = tmpVer.z < minBd.z ? tmpVer.z : minBd.z;
            //
            maxBd.x = tmpVer.x > maxBd.x ? tmpVer.x : maxBd.x;
            maxBd.y = tmpVer.y > maxBd.y ? tmpVer.y : maxBd.y;
            maxBd.z = tmpVer.z > maxBd.z ? tmpVer.z : maxBd.z;
    }

}

// edge conversion 
/*
std::map<std::pair<int, int>, int> edgeMap; 
for(size_t i = 0; i < vcg_mesh.face.size(); ++i) {
    const MyFace &f = vcg_mesh.face[i];
    for(int j = 0; j < 3; ++j) {
        int v0 = vcg::tri::Index(vcg_mesh, f.V(j));
        int v1 = vcg::tri::Index(vcg_mesh, f.V((j+1)%3));
        if (v0 > v1) std::swap(v0, v1);

        std::pair<int, int> edgeKey = std::make_pair(v0, v1);
        if(edgeMap.find(edgeKey)== edgeMap.end()) { 
            edgeMap[edgeKey] = rht_mesh.vecEdge.size();
            rht_mesh.vecEdge.push_back(_SimpleMeshEdge(v0,v1));
        }
    }
}
    */

// triangle conversion 
// reserves space in the vecTriangle vector of rht_mesh to store all triangles from vcg_mesh
rht_mesh.vecTriangle.reserve(vcg_mesh.face.size());
// Looping through faces 
// This loop iterates over all faces in vcg_mesh. For each face, it creates temporary triangle (tmpTri)
// It retrieves the vertex indices (v0, v1, v2) for the current face.
for(size_t i = 0; i < vcg_mesh.face.size(); ++i) {
    _SimpleMeshTriangle tmpTri;
    _SimpleMeshEdge tmpEdge;
    const MyFace &f = vcg_mesh.face[i];
    int v0 = vcg::tri::Index(vcg_mesh, f.V(0));
    int v1 = vcg::tri::Index(vcg_mesh, f.V(1));
    int v2 = vcg::tri::Index(vcg_mesh, f.V(2)); 
    tmpTri.v0 = v0;
	tmpTri.v1 = v1;
	tmpTri.v2 = v2;
	OrientTriangles.push_back(tmpTri.v0);
    OrientTriangles.push_back(tmpTri.v1);
    OrientTriangles.push_back(tmpTri.v2);
	// check the existence of tree edges
    // Calculate Edge Vectors & Normals 
    // Calculates the normal vector for the face using the cross product of two edge vectors and normalizes it.
    Vector3 leftVec, rightVec;
    leftVec[0] = rht_mesh.vecVertex[tmpTri.v2].x - rht_mesh.vecVertex[tmpTri.v1].x;
    leftVec[1] = rht_mesh.vecVertex[tmpTri.v2].y - rht_mesh.vecVertex[tmpTri.v1].y;
    leftVec[2] = rht_mesh.vecVertex[tmpTri.v2].z - rht_mesh.vecVertex[tmpTri.v1].z;

    rightVec[0] = rht_mesh.vecVertex[tmpTri.v0].x - rht_mesh.vecVertex[tmpTri.v1].x;
    rightVec[1] = rht_mesh.vecVertex[tmpTri.v0].y - rht_mesh.vecVertex[tmpTri.v1].y;
    rightVec[2] = rht_mesh.vecVertex[tmpTri.v0].z - rht_mesh.vecVertex[tmpTri.v1].z;
    leftVec = leftVec ^ rightVec;
    unitize(leftVec);
    //leftVec = leftVec / norm(leftVec);
    meshNormal.push_back(leftVec);
    // 
    tmpTri.sortVertices();
	
    // Edge Mapping and Assignement
    // tmpEdgePair: Creates a pair of vertex indices (v0, v1) representing an edge.
	std::pair<int, int> tmpEdgePair(tmpTri.v0, tmpTri.v1);
    // mIter: Declares an iterator for the edgeMapping map.
    std::map<std::pair<int, int>, int, myPairCompare>::iterator mIter;
    // Searches for the edge in the edgeMapping map
    mIter = edgeMapping.find(tmpEdgePair);
    // if the edge is not found in edgeMapping
    if (mIter == edgeMapping.end()) {// new edge
        // sets the vertices of tmpEdge
        tmpEdge.v0 = tmpEdgePair.first;
        tmpEdge.v1 = tmpEdgePair.second;
        // sets the adjacency information 
        // assigns the index of the current triangle to the first position in the AdjTri array of tmpEdge
        tmpEdge.AdjTri[0] = rht_mesh.vecTriangle.size();
        tmpEdge.AdjTriNum = 1;
        // Adds tmpEdge to rht_mesh.vecEdge
        rht_mesh.vecEdge.push_back(tmpEdge);
        // Updates tmpTri.e01 with the index of the new edge
        tmpTri.e01 = rht_mesh.vecEdge.size() - 1;
        // Adds the edge to edgeMapping.
        edgeMapping[tmpEdgePair] = rht_mesh.vecEdge.size() - 1;
    	} else {// existed already
            // Increments the number of triangles adjacent to this edge by 1
            rht_mesh.vecEdge[mIter->second].AdjTriNum++;
            // assigns the index of the current triangle to the second position in the AdjTri array, 
            // indicating that this edge is now adjacent to two triangles.
            rht_mesh.vecEdge[mIter->second].AdjTri[1] = rht_mesh.vecTriangle.size();
            // assigns the index of the existing edge to the triangle's edge index (e01)
            tmpTri.e01 = mIter->second;
            }
            //
    // Repeat for other edges (v1,v2)
    tmpEdgePair.first = tmpTri.v1;
    tmpEdgePair.second = tmpTri.v2;
    mIter = edgeMapping.find(tmpEdgePair);
    if (mIter == edgeMapping.end()) {// new edge
        tmpEdge.v0 = tmpEdgePair.first;
        tmpEdge.v1 = tmpEdgePair.second;
        tmpEdge.AdjTri[0] = rht_mesh.vecTriangle.size();
        tmpEdge.AdjTriNum = 1;
        //
        rht_mesh.vecEdge.push_back(tmpEdge);
        tmpTri.e12 = rht_mesh.vecEdge.size() - 1;
        edgeMapping[tmpEdgePair] = rht_mesh.vecEdge.size() - 1;
    } else {// existed already
        rht_mesh.vecEdge[mIter->second].AdjTriNum++;
        rht_mesh.vecEdge[mIter->second].AdjTri[1] = rht_mesh.vecTriangle.size();
        //
        tmpTri.e12 = mIter->second;
    }
    //
    // Repeat for other edges v0,v2
    tmpEdgePair.first = tmpTri.v0;
    tmpEdgePair.second = tmpTri.v2;
    mIter = edgeMapping.find(tmpEdgePair);
    if (mIter == edgeMapping.end()) {// new edge
        tmpEdge.v0 = tmpEdgePair.first;
        tmpEdge.v1 = tmpEdgePair.second;
        tmpEdge.AdjTri[0] = rht_mesh.vecTriangle.size();
        tmpEdge.AdjTriNum = 1;
        //
        rht_mesh.vecEdge.push_back(tmpEdge);
        tmpTri.e02 = rht_mesh.vecEdge.size() - 1;
        edgeMapping[tmpEdgePair] = rht_mesh.vecEdge.size() - 1;
    } else {// existed already
        rht_mesh.vecEdge[mIter->second].AdjTriNum++;
        rht_mesh.vecEdge[mIter->second].AdjTri[1] = rht_mesh.vecTriangle.size();
        //
        tmpTri.e02 = mIter->second;
    }
        //
    rht_mesh.vecTriangle.push_back(tmpTri);
}                    
	// assign incident edges information to vertex
    // this loop iterates over all edges in rht_mesh.vecEdge
    for (int i = 0; i < int(rht_mesh.vecEdge.size()); i++) {
        // adds the index of the current edge (i) to the adjEdges vector of the vertices.
        rht_mesh.vecVertex[rht_mesh.vecEdge[i].v0].adjEdges.push_back(i);
        rht_mesh.vecVertex[rht_mesh.vecEdge[i].v1].adjEdges.push_back(i);
    }
    std::cout << "Done... " << vcg_mesh.vert.size() << " " << vcg_mesh.face.size() << std::endl;
    std::cout << "ver... " << rht_mesh.vecVertex.size() << " tri " << rht_mesh.vecTriangle.size()
              << " edge" << rht_mesh.vecEdge.size() << std::endl;
    //
    edgeMapping.clear();
    return;
	/*
    // v0 < v1 < v2
    if (v0 > v1) std::swap(v0, v1);
    if (v1 > v2) std::swap(v1, v2);
    if (v0 > v1) std::swap(v0, v1);

    int e01 = edgeMap[std::make_pair(v0, v1)];
    int e12 = edgeMap[std::make_pair(v1, v2)];
    int e02 = edgeMap[std::make_pair(v0, v2)];

    rht_mesh.vecTriangle[i] = _SimpleMeshTriangle(v0, v1, v2, e01, e12, e02);
	*/
}


/* Function to performe the inverse of MeshConverter operation 
* @param : _SimpleMesh &input_mesh reference to a ReebHanThun mesh 
* @param : MyMesh &output_mesh reference to a VCG lib mesh to store the result of the conversion
* @result : output_mesh will be modified to store input_mesh converted in a VCG mesh
            do not modify input_mesh
*/
void ReverseMeshConverter(const _SimpleMesh &input_mesh, MyMesh &output_mesh){

    // Add vertices to the mesh
    // iterates over each vertex in input_mesh.vecVertex
    for(auto& p : input_mesh.vecVertex) {
        MyMesh::VertexType v;
        // Converts the coordinates of the vertex (p.x, p.y, p.z) to a vcg::Point3<float> object
        vcg::Point3<float> point = vcg::Point3<float>(p.x, p.y, p.z);
        v.P() = point;
        // add vertex
        output_mesh.vert.push_back(v);
    }
    // update vertex count 
    output_mesh.vn = output_mesh.vert.size();

    // Add faces to the mesh 
    // loops through triangles of input_mesh.vecTriangle
    for(auto& tri: input_mesh.vecTriangle) {
        // Creates a new face (f) of type MyMesh::FaceType
        MyMesh::FaceType f;
        // Assigns the vertices of the face using the indices from the triangle (tri.v0, tri.v1, tri.v2)
        f.V(0) = &output_mesh.vert[tri.v0];
        f.V(1) = &output_mesh.vert[tri.v1];
        f.V(2) = &output_mesh.vert[tri.v2];
        // add face
        output_mesh.face.push_back(f);
    }
    // update face count in output_mesh
    output_mesh.fn = output_mesh.face.size();

    // Update the normals (for face and for vertex) and bounding box
    tri::UpdateNormal<MyMesh>::PerFaceNormalized(output_mesh);
    tri::UpdateNormal<MyMesh>::PerVertexNormalized(output_mesh);
    tri::UpdateBounding<MyMesh>::Box(output_mesh);

}

void PrintSimpleMesh(const _SimpleMesh &mesh) {
    
    std::cout << "Vertices: " << mesh.vecVertex.size() << std::endl; /*
    for (size_t i = 0; i < mesh.vecVertex.size(); ++i) {
        const _SimpleMeshVertex &v = mesh.vecVertex[i];
        std::cout << "Vertex " << i << ": (" << v.x << ", " << v.y << ", " << v.z << ")" << std::endl;
    }
    */
    std::cout << "Edges: " << mesh.vecEdge.size() << std::endl; /*
    for (size_t i = 0; i < mesh.vecEdge.size(); ++i) {
        const _SimpleMeshEdge &e = mesh.vecEdge[i];
        std::cout << "Edge " << i << ": (" << e.v0 << ", " << e.v1 << ")" << std::endl;
    }
    */
    
    std::cout << "Triangles: " << mesh.vecTriangle.size() << std::endl; /*
    for (size_t i = 0; i < mesh.vecTriangle.size(); ++i) {
        const _SimpleMeshTriangle &t = mesh.vecTriangle[i];
        std::cout << "Triangle " << i << ": (" << t.v0 << ", " << t.v1 << ", " << t.v2 << ") with edges (" << t.e01 << ", " << t.e12 << ", " << t.e02 << ")" << std::endl;
    }
    */
}


MyMesh& ComputeBasis(MyMesh & m_vcg) {
    Vector3 distinctDirection;
    const float fEnlargeFactor = 10000.f;
    // vectors of sets that will store basis loops for handles and tunnels
	std::vector<std::set<int> > v_basis_loops;
	std::vector<std::set<int> > h_basis_loops;
    psbmReebGraph reebGraph;
    std::vector<Vector3>   meshNormal;
    std::vector<int> OrientTriangles;
	std::set<int> extraVertices;
	int nOrgTriangleSize = 0;
	double BoundingBoxRadius;
    int genus = 0;
	_SimpleMeshVertex minBd;
    _SimpleMeshVertex maxBd;

    // declaration of an object of type _SimpleMesh -- ReebHanTun 
    _SimpleMesh m_rht;

    // convert the vcg mesh in reebhantun _SimpleMesh
    MeshConverter(minBd, maxBd, m_vcg, m_rht, meshNormal, OrientTriangles, fEnlargeFactor);

    // The mesh is in the right format 
    // Start ReebHanTun business

    m_rht.SetMeshNormalPtr(&meshNormal);

    // triangles in TRIS are in the same order as triangles in mesh.vecTriangle;
	//
	nOrgTriangleSize = m_rht.vecTriangle.size();
	//
    std::vector<std::vector<std::pair<int, int> > > meshBoundaries;
	if (CheckBoundaries(m_rht, meshBoundaries))
	{// it is a mesh with bondary
		int EulerCharacteristic = m_rht.vecVertex.size() + m_rht.vecTriangle.size() - m_rht.vecEdge.size();
		genus = 1 - (EulerCharacteristic + meshBoundaries.size()) / 2;
		if (genus)
		{
			CloseHoles(m_rht, meshBoundaries, meshNormal, OrientTriangles, extraVertices);
		}
	}
	else
	{// it is a closed mesh
		int EulerCharacteristic = m_rht.vecVertex.size() + m_rht.vecTriangle.size() - m_rht.vecEdge.size();
		genus = 1 - EulerCharacteristic / 2;
	}
	if (!genus)
	{
		std::cout << "NOTHING IS COMPUTED : "  << std::endl;
		std::cout << " ---- MESH HAS GENUS 0!" << std::endl;
		exit(1);
	}
	std::cout << "Mesh has genus : " << genus << std::endl;
	//

    // compute the bounding box
    RandomUniqueDirection(m_rht, distinctDirection);
    reebGraph.ReserveSpaceForEdges(m_rht.vecEdge.size());
    double* scalarField = new double[m_rht.vecVertex.size()];
    int perDir = 0;
    int rayDir = 0;
    for (unsigned int i = 0; i < m_rht.vecVertex.size(); i++)
        {
	        scalarField[i] = m_rht.vecVertex[i].x * distinctDirection[0] +
		    m_rht.vecVertex[i].y * distinctDirection[1] +
		    m_rht.vecVertex[i].z * distinctDirection[2];
        }
    reebGraph.SetHeightDirection(distinctDirection);
    reebGraph.AssignData(&m_rht, scalarField);
    reebGraph.scalarDir = 1;
    std::cout << std::endl;
    {
	    // boost::progress_timer t;
	    //
	    reebGraph.ComputeReebGraph();
	    //
    }
    {
		// boost::progress_timer t;
		////
		std::cout << "Time for mapping and linking :" << std::endl;
		//reebGraph.ComputeCycleAndPairing();
		reebGraph.ComputingCycle_max_tree();
		////
		//std::cout << "mapping" << std::endl;
		//// computing the cycle on surface
		reebGraph.compute_path_on_mesh_for_each_simplified_arc();//pathArcOnMesh, offsetPathArcOnMesh);
		////
		//std::cout << "embed" << std::endl;
		reebGraph.EmbedCycleAsEdgePathOnMesh();
		//std::cout << "linking" << std::endl;
		reebGraph.LinkNumberMatrixComputing();
	}
	//

	//
	if (extraVertices.empty())
		CycleLocalOptimization(m_rht, reebGraph, OrientTriangles, v_basis_loops, h_basis_loops, 1.f / fEnlargeFactor);
	else
		CycleLocalOptimization_bdry(m_rht, reebGraph, OrientTriangles, v_basis_loops, h_basis_loops, extraVertices, 1.f / fEnlargeFactor);
	//
		

	///////////////////////////////////////////////////
    std::cout << "Handle and tunnel loops written in files :  \n";
	FilesOutputForOptimalCycles files_out_op;
    std::string OutputFileName = "provaaaaaa";
	files_out_op.InitMeshPtr(&m_rht);
	files_out_op.WriteCyclesInformation(OutputFileName.c_str(), v_basis_loops, h_basis_loops);
	int orgVertexSize = m_rht.vecVertex.size();
	if (!extraVertices.empty())
		orgVertexSize = *extraVertices.begin();
	files_out_op.WriteGeomviewListFormat(OutputFileName.c_str(), v_basis_loops, h_basis_loops, OrientTriangles, orgVertexSize, nOrgTriangleSize, 1.f / fEnlargeFactor);

    //MyMesh test;
    //ReverseMeshConverter(m_rht, test);
    //vcg::tri::io::ExporterOFF<MyMesh>::Save(test,"test.off");
    //std::cout << "num edges: " << m_rht.vecEdge.size();
	//std::cout << std::endl; 

	for (const auto& s : v_basis_loops) {
        for (const auto& elem : s) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;  // Nuova riga dopo ogni set
    }

	for (const auto& s : h_basis_loops) {
        for (const auto& elem : s) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;  // Nuova riga dopo ogni set
    }
    
    
    
    MyMesh loops;
    for (const auto& s : v_basis_loops) {
        for (const auto& elem : s) {
            if(elem > 0 && elem < m_rht.vecEdge.size())
            {
                int v0 = m_rht.vecEdge[elem].v0;
                int v1 = m_rht.vecEdge[elem].v1;         
                tri::Allocator<MyMesh>::AddEdge(loops, m_vcg.vert[v0].P(),m_vcg.vert[v1].P());
            }
        }
    }
    
    
    for (const auto& s : h_basis_loops) {
        for (const auto& elem : s) {
            if(elem > 0 && elem < m_rht.vecEdge.size())
            {
                int v0 = m_rht.vecEdge[elem].v0;
                int v1 = m_rht.vecEdge[elem].v1;         
                tri::Allocator<MyMesh>::AddEdge(loops, m_vcg.vert[v0].P(),m_vcg.vert[v1].P());
            }
        }
    }
    
    return loops;
}


int main(int argc, char **argv)  {
    // declaration of an object of type MyMesh -- VCG lib
    MyMesh m_vcg;
    
    // load a mesh with VCG lib
    if(vcg::tri::io::ImporterOFF<MyMesh>::Open(m_vcg,argv[1])!=vcg::tri::io::ImporterOFF<MyMesh>::NoError)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }
  
    // Calculate the number of vertices, edges and faces
    int numVertices = m_vcg.VN();
    int numEdges = m_vcg.EN();
    int numFaces = m_vcg.FN();
    

    
    // Print the results
    printf("Number of vertices: %d\n", numVertices);
    printf("Number of edges: %d\n", numEdges);
    printf("Number of faces: %d\n", numFaces);
   
    // Now I've loaded the vcg mesh MyMesh

    ReebHanTunWrapper<MyMesh> wrapper(m_vcg);
    MyMesh loops;
    wrapper.ComputeBasis(loops);
    //wrapper.SetBaseMesh(m_vcg);

    printf("Number of handle loops found: %d\n", wrapper.GetNumHandleLoops());
    printf("Number of Tunnel loops found: %d\n", wrapper.GetNumHandleLoops());

    MyMesh h_loops;
    MyMesh v_loops;
    
    

    wrapper.GetHandleLoops(h_loops, -1);
    wrapper.GetTunnelLoops(v_loops, -1);
    vcg::tri::io::ExporterPLY<MyMesh>::Save(h_loops,"handle_loops.ply", vcg::tri::io::Mask::IOM_EDGEINDEX); 
    vcg::tri::io::ExporterPLY<MyMesh>::Save(v_loops,"tunnel_loops.ply", vcg::tri::io::Mask::IOM_EDGEINDEX);
    vcg::tri::io::ExporterPLY<MyMesh>::Save(loops,"loops.ply", vcg::tri::io::Mask::IOM_EDGEINDEX);
    //vcg::tri::io::ExporterOFF<MyMesh>::Save(loops,"loops.off");
    return 0;
}

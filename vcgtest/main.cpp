#include<stdio.h>
#include<vcg/complex/complex.h>
#include <vcg/space/point3.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/bounding.h>
// input output
#include <wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export_off.h>
#include<wrap/io_trimesh/export_ply.h>
// mesh ReebHanTun
#include<SimpleMesh.h>

// include ComputeReebGraph.cpp
//
#include "psbmReebGraph.h"
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include "SimpleMesh.h"
#include "FilesOutputForOptimalCycles.h"

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
//

using namespace vcg;
using namespace std;

class MyVertex;
class MyEdge;
class MyFace;

struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>::AsVertexType, vcg::Use<MyEdge>::AsEdgeType,   vcg::Use<MyFace>::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b, vcg::vertex::BitFlags, vcg::vertex::VEAdj >{};
class MyEdge : public vcg::Edge<MyUsedTypes> {};
class MyFace  : public vcg::Face < MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::Color4b, vcg::face::BitFlags, vcg::face::EFAdj > {};
class MyMesh  : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyEdge>, std::vector<MyFace> > {};



// function to convert a vcg mesh type _SimpleMesh in a reebhantun mesh _SimpleMesh

void MeshConverter (_SimpleMeshVertex &minBd,
                                      _SimpleMeshVertex &maxBd, const MyMesh & vcg_mesh,  _SimpleMesh & rht_mesh, std::vector<Vector3> &meshNormal, std::vector<int> &OrientTriangles, const float fEnlargeFactor ) {
std::map<std::pair<int, int>, int, myPairCompare> edgeMapping;
// vertex conversion 
OrientTriangles.reserve(vcg_mesh.face.size());
// resize vecVertex of _SimpleMesh to store as many vertex as Mymesh
rht_mesh.vecVertex.reserve(vcg_mesh.vert.size());
for(size_t i = 0; i < vcg_mesh.vert.size(); ++i) {
	_SimpleMeshVertex tmpVer;
    const MyVertex &v = vcg_mesh.vert[i];
	tmpVer.x = v.P().X() * fEnlargeFactor;
    tmpVer.y = v.P().Y() * fEnlargeFactor;
	tmpVer.z = v.P().Z() * fEnlargeFactor;
	rht_mesh.vecVertex.push_back(tmpVer);
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
rht_mesh.vecTriangle.reserve(vcg_mesh.face.size());
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
    //
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
	
	std::pair<int, int> tmpEdgePair(tmpTri.v0, tmpTri.v1);
    std::map<std::pair<int, int>, int, myPairCompare>::iterator mIter;
    mIter = edgeMapping.find(tmpEdgePair);
    if (mIter == edgeMapping.end()) {// new edge
        tmpEdge.v0 = tmpEdgePair.first;
        tmpEdge.v1 = tmpEdgePair.second;
        tmpEdge.AdjTri[0] = rht_mesh.vecTriangle.size();
        tmpEdge.AdjTriNum = 1;
        //
        rht_mesh.vecEdge.push_back(tmpEdge);
        tmpTri.e01 = rht_mesh.vecEdge.size() - 1;
        edgeMapping[tmpEdgePair] = rht_mesh.vecEdge.size() - 1;
    	} else {// existed already
            rht_mesh.vecEdge[mIter->second].AdjTriNum++;
            rht_mesh.vecEdge[mIter->second].AdjTri[1] = rht_mesh.vecTriangle.size();
            //
            tmpTri.e01 = mIter->second;
            }
            //
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
    for (int i = 0; i < int(rht_mesh.vecEdge.size()); i++) {
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
* @param : MyMesh &output_mesh reference to a vcg lib mesh to store the result of the conversion
* @result : output_mesh will be modified to store input_mesh converted in a vcg mesh
            do not modify input_mesh
*/
void ReverseMeshConverter(const _SimpleMesh &input_mesh, MyMesh &output_mesh){

    // Add vertices to the mesh
    for(auto& p : input_mesh.vecVertex) {
        MyMesh::VertexType v;
        vcg::Point3<float> point = vcg::Point3<float>(p.x, p.y, p.z);
        v.P() = point;
        output_mesh.vert.push_back(v);
    }
    output_mesh.vn = output_mesh.vert.size();

    // Add faces to the mesh 
    for(auto& tri: input_mesh.vecTriangle) {
        MyMesh::FaceType f;
        f.V(0) = &output_mesh.vert[tri.v0];
        f.V(1) = &output_mesh.vert[tri.v1];
        f.V(2) = &output_mesh.vert[tri.v2];
        output_mesh.face.push_back(f);
    }
    output_mesh.fn = output_mesh.face.size();

    // Update the normals and bounding box
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


int main(int argc, char **argv)  {
    Vector3 distinctDirection;
    const float fEnlargeFactor = 10000.f;
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

    printf("prova\n");
    // definition of type MyMesh -- VCG lib
    MyMesh m_vcg;
    // definition of type _SimpleMesh -- ReebHanTun 
    _SimpleMesh m_rht;
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
  // convert the vcg mesh in reebhantun _SimpleMesh

   MeshConverter(minBd, maxBd, m_vcg, m_rht, meshNormal, OrientTriangles, fEnlargeFactor);
   

    
   // The mesh is in the right format 
   // Start ReebHanTun business

   // function load data!!

  // m_rht.LoadMeshInOFFformat(minBd, maxBd, meshNormal, argv[1], OrientTriangles, fEnlargeFactor);//"E:\\RProgramming\\Models\\torus.off");// "D:\\MeshModels\\OFF-models\\HighGenusCubeHC1.off");//

   m_rht.SetMeshNormalPtr(&meshNormal);
   //PrintSimpleMesh(m_rht);
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

// end function load data
// begin main ComputeReebGraph.cpp
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
	
	

    MyMesh test;
    ReverseMeshConverter(m_rht, test);

  //vcg::tri::io::ExporterOFF<MyMesh>::Save(test,"test.off");

    std::cout << "num edges: " << m_rht.vecEdge.size();
	std::cout << std::endl; 

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
			if(elem > 0 && elem < m_rht.vecEdge.size()) {
			_SimpleMeshEdge edge = m_rht.vecEdge[elem];
			int v0 = edge.v0;
			int v1 = edge.v1;
			loops.vert.push_back(m_vcg.vert[v0]);
			loops.vert.push_back(m_vcg.vert[v1]);
			}                  
    }
  }

   for (const auto& s : h_basis_loops) {
        for (const auto& elem : s) {
			if(elem > 0 && elem < m_rht.vecEdge.size()) {
			_SimpleMeshEdge edge = m_rht.vecEdge[elem];
			int v0 = edge.v0;
			int v1 = edge.v1;
			loops.vert.push_back(m_vcg.vert[v0]);
			loops.vert.push_back(m_vcg.vert[v1]);
			}                  
    }
  }
  loops.vn = loops.vert.size();
 // vcg::tri::io::ExporterPLY<MyMesh>::Save(loops,"wallMesh.ply", vcg::tri::io::Mask::IOM_VERTCOLOR); 
 vcg::tri::io::ExporterOFF<MyMesh>::Save(loops,"loops.off");
  return 0;
}

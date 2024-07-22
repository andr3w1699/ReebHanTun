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



bool CheckBoundaries(_SimpleMesh &inMesh, std::vector<std::vector<std::pair<int,int> > > &boundries)
{
	bool ret = false;
	std::vector<bool> edge_flag(inMesh.vecEdge.size(), false);
	for (unsigned int i = 0; i < inMesh.vecEdge.size(); i++)
	{
		if (!edge_flag[i] && inMesh.vecEdge[i].AdjTriNum == 1)
		{// this is a boundary edge
			ret = true;
			std::vector<std::pair<int, int> > bdrEdges;
			//
			int terminate_vertex = inMesh.vecEdge[i].v0;
			int current_vertex = inMesh.vecEdge[i].v1;
			int current_edge = i;
			int loop_counter = 0;
			while (current_vertex != terminate_vertex)
			{
				bdrEdges.push_back(std::pair<int, int>(current_edge, inMesh.vecEdge[current_edge].AdjTri[0]));
				edge_flag[current_edge] = true;
				//
				unsigned int evid = 0;
				for (;evid < inMesh.vecVertex[current_vertex].adjEdges.size(); evid++)
				{
					int rot_edge_id = inMesh.vecVertex[current_vertex].adjEdges[evid];
					if (rot_edge_id != current_edge && inMesh.vecEdge[rot_edge_id].AdjTriNum == 1)
						break;
				}
				if (evid < inMesh.vecVertex[current_vertex].adjEdges.size())
				{
					current_edge = inMesh.vecVertex[current_vertex].adjEdges[evid];
					current_vertex = inMesh.vecEdge[current_edge].v0 + inMesh.vecEdge[current_edge].v1 - current_vertex;
				}
				else
				{
					std::cout << "READING MESH MODEL ERROR : " << std::endl;
					std::cout << "--- CAN NOT FIND ADJACENT BOUDNARY EDGE !" << std::endl;
					std::cout << "--- PLEASE CHECK THE MESH MODEL" << std::endl;
					exit(0);
				}
				loop_counter++;
				if (loop_counter > inMesh.vecEdge.size())
				{
					std::cout << "READING MESH MODEL ERROR : " << std::endl;
					std::cout << "--- PLEASE CHECK THE MESH MODEL" << std::endl;
					exit(0);
				}
			}
			edge_flag[current_edge] = true;
			bdrEdges.push_back(std::pair<int, int>(current_edge, inMesh.vecEdge[current_edge].AdjTri[0]));
			//
			boundries.push_back(bdrEdges);
		}
		edge_flag[i] = true;
	}
	return ret;
}
void FindEdgeOrientation(_SimpleMesh &inMesh, const int fid, std::vector<int> &OrientTriangles, std::pair<int, int> &EdgeOrientation)
{
	int start_index = 3 * fid;
	std::vector<int> TriOrientation(OrientTriangles.begin() + start_index, OrientTriangles.begin() + start_index + 3);
	TriOrientation.push_back(TriOrientation.front());
	//
	for (unsigned int i = 0; i < 3; i++)
	{
		if ((EdgeOrientation.first == TriOrientation[i] && EdgeOrientation.second == TriOrientation[i + 1]) ||
			(EdgeOrientation.second == TriOrientation[i] && EdgeOrientation.first == TriOrientation[i + 1]) )
		{
			EdgeOrientation.first = TriOrientation[i + 1];
			EdgeOrientation.second = TriOrientation[i];
		}
	}
	return;
}
void CloseHoles(_SimpleMesh &inMesh, std::vector<std::vector<std::pair<int, int> > > &meshBoundaries,
	std::vector<Vector3> &inMeshNormal,
	std::vector<int> &OrientTriangles, std::set<int> &extraVertices)
{
	for (unsigned int b = 0; b < meshBoundaries.size(); b++)
	{
		Vector3 centroid(0.0, 0.0, 0.0);
		int curEdge = 0;
		int nextEdge = 0;
		int curVertex = 0;
		for (unsigned int i = 0; i < meshBoundaries[b].size(); i++)
		{// the edges on the boundary is ordered
			curEdge = meshBoundaries[b][i].first;
			nextEdge = meshBoundaries[b][0].first;
			if (i < meshBoundaries[b].size() - 1)
				nextEdge = meshBoundaries[b][i + 1].first;
			// shared vertex
			curVertex = inMesh.vecEdge[curEdge].v0;
			if (curVertex != inMesh.vecEdge[nextEdge].v0 &&
				curVertex != inMesh.vecEdge[nextEdge].v1)
				curVertex = inMesh.vecEdge[curEdge].v1;
			//
			centroid = centroid + Vector3(	inMesh.vecVertex[curVertex].x,
											inMesh.vecVertex[curVertex].y,
											inMesh.vecVertex[curVertex].z);
			//
		}
		centroid = centroid / meshBoundaries[b].size();
		//
		// add vertex
		int centroid_index = inMesh.vecVertex.size();
		extraVertices.insert(centroid_index);
		_SimpleMeshVertex tmpVer;
		tmpVer.x = centroid[0] ;
		tmpVer.y = centroid[1] ;
		tmpVer.z = centroid[2] ;
		//
		inMesh.vecVertex.push_back(tmpVer);
		//
		// add edges and triangles
		std::map<int, int> edge_mapping;
		for (unsigned int i = 0; i < meshBoundaries[b].size(); i++)
		{
			// Find orientation of this edge
			curEdge = meshBoundaries[b][i].first;
			std::pair<int, int> EdgeOrientation(inMesh.vecEdge[curEdge].v0, inMesh.vecEdge[curEdge].v1);
			FindEdgeOrientation(inMesh, meshBoundaries[b][i].second, OrientTriangles, EdgeOrientation);
			//
			OrientTriangles.push_back(EdgeOrientation.first);
			OrientTriangles.push_back(EdgeOrientation.second);
			OrientTriangles.push_back(centroid_index);
			//
			_SimpleMeshTriangle tmpTri;
			_SimpleMeshEdge tmpEdge;
			//
			tmpTri.v0 = EdgeOrientation.first;
			tmpTri.v1 = EdgeOrientation.second;
			tmpTri.v2 = centroid_index;
			//
			//
			Vector3 leftVec, rightVec;
			leftVec[0] = inMesh.vecVertex[tmpTri.v2].x - inMesh.vecVertex[tmpTri.v1].x;
			leftVec[1] = inMesh.vecVertex[tmpTri.v2].y - inMesh.vecVertex[tmpTri.v1].y;
			leftVec[2] = inMesh.vecVertex[tmpTri.v2].z - inMesh.vecVertex[tmpTri.v1].z;

			rightVec[0] = inMesh.vecVertex[tmpTri.v0].x - inMesh.vecVertex[tmpTri.v1].x;
			rightVec[1] = inMesh.vecVertex[tmpTri.v0].y - inMesh.vecVertex[tmpTri.v1].y;
			rightVec[2] = inMesh.vecVertex[tmpTri.v0].z - inMesh.vecVertex[tmpTri.v1].z;
			leftVec = leftVec^rightVec;
			unitize(leftVec);
			//leftVec = leftVec / norm(leftVec);
			inMeshNormal.push_back(leftVec);
			//
			tmpTri.sortVertices();
			// [v0 v1] existing edge
			inMesh.vecEdge[curEdge].AdjTriNum++;
			inMesh.vecEdge[curEdge].AdjTri[1]=inMesh.vecTriangle.size();
			//
			tmpTri.e01 = curEdge;
			//
			std::map<int, int>::iterator mIter;
			mIter = edge_mapping.find(tmpTri.v1);
			if (mIter == edge_mapping.end())
			{// new edge
				tmpEdge.v0 = tmpTri.v1;
				tmpEdge.v1 = tmpTri.v2;
				tmpEdge.AdjTri[0] = inMesh.vecTriangle.size();
				tmpEdge.AdjTriNum = 1;
				//
				inMesh.vecEdge.push_back(tmpEdge);
				tmpTri.e12 = inMesh.vecEdge.size() - 1;
				edge_mapping[tmpTri.v1] = tmpTri.e12;
				//
				inMesh.vecVertex[tmpEdge.v0].adjEdges.push_back(tmpTri.e12);
				inMesh.vecVertex[tmpEdge.v1].adjEdges.push_back(tmpTri.e12);
			}
			else
			{// existed edge
				inMesh.vecEdge[mIter->second].AdjTriNum++;
				inMesh.vecEdge[mIter->second].AdjTri[1] = inMesh.vecTriangle.size();
				//
				tmpTri.e12 = mIter->second;
			}
			mIter = edge_mapping.find(tmpTri.v0);
			if (mIter == edge_mapping.end())
			{
				// new edge
				tmpEdge.v0 = tmpTri.v0;
				tmpEdge.v1 = tmpTri.v2;
				tmpEdge.AdjTri[0] = inMesh.vecTriangle.size();
				tmpEdge.AdjTriNum = 1;
				//
				inMesh.vecEdge.push_back(tmpEdge);
				tmpTri.e02 = inMesh.vecEdge.size() - 1;
				edge_mapping[tmpTri.v0] = tmpTri.e02;
				//
				inMesh.vecVertex[tmpEdge.v0].adjEdges.push_back(tmpTri.e02);
				inMesh.vecVertex[tmpEdge.v1].adjEdges.push_back(tmpTri.e02);
			}
			else
			{// existed already
				inMesh.vecEdge[mIter->second].AdjTriNum++;
				inMesh.vecEdge[mIter->second].AdjTri[1] = inMesh.vecTriangle.size();
				//
				tmpTri.e02 = mIter->second;
			}
			//
			inMesh.vecTriangle.push_back(tmpTri);
		}
	}
	return;
}

void RandomUniqueDirection(_SimpleMesh &mesh, Vector3 &uniDirection)
{
	//
	std::vector<Vector3> vecEdgeDirections (2 * mesh.vecEdge.size());
	// all vectors are pointed to positive Z direction
	for (unsigned int i = 0; i < mesh.vecEdge.size(); i++)
	{
		vecEdgeDirections[i][0] = mesh.vecVertex[mesh.vecEdge[i].v0].x - mesh.vecVertex[mesh.vecEdge[i].v1].x;
		vecEdgeDirections[i][1] = mesh.vecVertex[mesh.vecEdge[i].v0].y - mesh.vecVertex[mesh.vecEdge[i].v1].y;
		vecEdgeDirections[i][2] = mesh.vecVertex[mesh.vecEdge[i].v0].z - mesh.vecVertex[mesh.vecEdge[i].v1].z;
		//
		if (vecEdgeDirections[i][2] < 0)
		{
			for (int j = 0; j < 3; j++)
				vecEdgeDirections[i][j] = - vecEdgeDirections[i][j];
		}
		//
		unitize(vecEdgeDirections[i]);
	}
	int mesh_edge_size = mesh.vecEdge.size();
	int opp_v_0;
	int opp_v_1;
	for (unsigned int i = 0; i < mesh.vecEdge.size(); i++)
	{
		if (mesh.vecEdge[i].AdjTriNum == 2)
		{
			opp_v_0 = mesh.vecEdge[i].AdjTri[0];
			if (mesh.vecTriangle[opp_v_0].v0 != mesh.vecEdge[i].v0 &&
				mesh.vecTriangle[opp_v_0].v0 != mesh.vecEdge[i].v1)
			{
				opp_v_0 = mesh.vecTriangle[opp_v_0].v0;
			}
			else
			{
				if (mesh.vecTriangle[opp_v_0].v1 != mesh.vecEdge[i].v0 &&
					mesh.vecTriangle[opp_v_0].v1 != mesh.vecEdge[i].v1)
				{
					opp_v_0 = mesh.vecTriangle[opp_v_0].v1;
				}
				else
				{
					opp_v_0 = mesh.vecTriangle[opp_v_0].v2;
				}
			}
			//
			opp_v_1 = mesh.vecEdge[i].AdjTri[1];
			if (mesh.vecTriangle[opp_v_1].v0 != mesh.vecEdge[i].v0 &&
				mesh.vecTriangle[opp_v_1].v0 != mesh.vecEdge[i].v1)
			{
				opp_v_1 = mesh.vecTriangle[opp_v_1].v0;
			}
			else
			{
				if (mesh.vecTriangle[opp_v_1].v1 != mesh.vecEdge[i].v0 &&
					mesh.vecTriangle[opp_v_1].v1 != mesh.vecEdge[i].v1)
				{
					opp_v_1 = mesh.vecTriangle[opp_v_1].v1;
				}
				else
				{
					opp_v_1 = mesh.vecTriangle[opp_v_1].v2;
				}
			}
			//
			vecEdgeDirections[i + mesh_edge_size][0] = mesh.vecVertex[opp_v_0].x - mesh.vecVertex[opp_v_1].x;
			vecEdgeDirections[i + mesh_edge_size][1] = mesh.vecVertex[opp_v_0].y - mesh.vecVertex[opp_v_1].y;
			vecEdgeDirections[i + mesh_edge_size][2] = mesh.vecVertex[opp_v_0].z - mesh.vecVertex[opp_v_1].z;
			//
			if (vecEdgeDirections[i + mesh_edge_size][2] < 0)
			{
				for (int j = 0; j < 3; j++)
					vecEdgeDirections[i + mesh_edge_size][j] = - vecEdgeDirections[i + mesh_edge_size][j];
			}
			//
			unitize(vecEdgeDirections[i + mesh_edge_size]);
		}
	}
	//
	srand(time(NULL));
	//
	int runTimesCounter = 0;
	int halfRandMax = RAND_MAX >> 2;
	bool bFindDirection = false;
	double minAangleValue = 1.0;
	double error_precision = 1e-7;
	double max_running_counter = 2000;
	while (!bFindDirection)
	{
		runTimesCounter++;
		for (int i = 0 ; i < 3; i++)
		{
			uniDirection[i] = rand() - halfRandMax;
		}
		if (uniDirection[2] < 0)
		{
			for (int i = 0; i < 3; i++)
				uniDirection[i] = - uniDirection[i];
		}
		//
		unitize(uniDirection);
		//
		bFindDirection = true;
		minAangleValue = 1.0;
		for (unsigned int i = 0; i < vecEdgeDirections.size(); i++)
		{
			double angleValue = uniDirection * vecEdgeDirections[i];
			if (angleValue < error_precision && angleValue > -error_precision)
			{
				bFindDirection = false;
				break;
			}
			if (std::abs(minAangleValue) > std::abs(angleValue))
				minAangleValue = angleValue;
		}
		if (runTimesCounter > max_running_counter)
			break;
	}
	//
	//std::cout << "run times : " << runTimesCounter << std::endl;
	//std::cout << "min angle : " << minAangleValue << std::endl;
	//
	return;
}
//
void CycleLocalOptimization(_SimpleMesh &locMesh, psbmReebGraph &reebgraph, std::vector<int> &OrientTriangles, std::vector<std::set<int> >& v_basis,
	std::vector<std::set<int> >& h_basis,  const float fEnlargeFactor);
void CycleLocalOptimization_bdry(_SimpleMesh &locMesh, psbmReebGraph &reebgraph, std::vector<int> &OrientTriangles,
							std::vector<std::set<int> >& out_v_basis_loops,
							std::vector<std::set<int> >& out_h_basis_loops,
							std::set<int> extraVertices, const float fEnlargeFactor);
////

// function to convert a vcg mesh type _SimpleMesh in a reebhantun mesh _SimpleMesh

void MeshConverter (const MyMesh & vcg_mesh,  _SimpleMesh & rht_mesh) {

// vertex conversion 

// resize vecVertex of _SimpleMesh to store as many vertex as Mymesh
rht_mesh.vecVertex.resize(vcg_mesh.vert.size());
for(size_t i = 0; i < vcg_mesh.vert.size(); ++i) {
    const MyVertex &v = vcg_mesh.vert[i];
    rht_mesh.vecVertex[i] = _SimpleMeshVertex(v.P().X(), v.P().Y(), v.P().Z());
}

// edge conversion 

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

// triangle conversion 
rht_mesh.vecTriangle.resize(vcg_mesh.face.size());
for(size_t i = 0; i < vcg_mesh.face.size(); ++i) {
    const MyFace &f = vcg_mesh.face[i];
    int v0 = vcg::tri::Index(vcg_mesh, f.V(0));
    int v1 = vcg::tri::Index(vcg_mesh, f.V(1));
    int v2 = vcg::tri::Index(vcg_mesh, f.V(2)); 

    // v0 < v1 < v2
    if (v0 > v1) std::swap(v0, v1);
    if (v1 > v2) std::swap(v1, v2);
    if (v0 > v1) std::swap(v0, v1);

    int e01 = edgeMap[std::make_pair(v0, v1)];
    int e12 = edgeMap[std::make_pair(v1, v2)];
    int e02 = edgeMap[std::make_pair(v0, v2)];

    rht_mesh.vecTriangle[i] = _SimpleMeshTriangle(v0, v1, v2, e01, e12, e02);
}
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

  // MeshConverter(m_vcg, m_rht);
   

    
   // The mesh is in the right format 
   // Start ReebHanTun business

   // function load data!!
   _SimpleMeshVertex minBd;
	_SimpleMeshVertex maxBd;
   m_rht.LoadMeshInOFFformat(minBd, maxBd, meshNormal, argv[1], OrientTriangles, fEnlargeFactor);//"E:\\RProgramming\\Models\\torus.off");// "D:\\MeshModels\\OFF-models\\HighGenusCubeHC1.off");//

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
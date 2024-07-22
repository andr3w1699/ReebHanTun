#include<stdio.h>
#include<vcg/complex/complex.h>
#include <vcg/space/point3.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/bounding.h>
// input output
#include <wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export_off.h>
// mesh ReebHanTun
#include<SimpleMesh.h>

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
    printf("prova\n");
    // definition of type MyMesh -- VCG lib
    MyMesh m_vcg;
    // definition of type _SimpleMesh -- ReebHanTun 
    _SimpleMesh m_rht;
    std::vector<Vector3>   meshNormal;
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

   MeshConverter(m_vcg, m_rht);
   
   // The mesh is in the right format 
   // Start ReebHanTun business
   m_rht.SetMeshNormalPtr(&meshNormal);
   //PrintSimpleMesh(m_rht);


   MyMesh test;
   ReverseMeshConverter(m_rht, test);
 
  vcg::tri::io::ExporterOFF<MyMesh>::Save(test,"test.off");
  return 0;
}
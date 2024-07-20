#include<stdio.h>
#include<vcg/complex/complex.h>
// input output
#include <wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export_off.h>
// mesh ReebHanTun
#include<SimpleMesh.h>

class MyFace;
class MyVertex;

struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>::AsVertexType,    vcg::Use<MyFace>::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face < MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::Color4b, vcg::face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};



_SimpleMesh MeshConverter (MyMesh vcg_mesh) {
_SimpleMesh rht_mesh;

return rht_mesh;
}


int main(int argc, char **argv)  {
    printf("prova");
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
  
    // Calculate the number of vertices and faces
    int numVertices = m_vcg.VN();
    int numFaces = m_vcg.FN();
    
    // Print the results
    printf("Number of vertices: %d\n", numVertices);
    printf("Number of faces: %d\n", numFaces);
   
  // Now I've loaded the vcg mesh MyMesh
  // convert the vcg mesh in reebhantun _SimpleMesh




  vcg::tri::io::ExporterOFF<MyMesh>::Save(m_vcg,"test.off");
  return 0;
}
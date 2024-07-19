#include<stdio.h>
#include<vcg/complex/complex.h>
// input output
#include <wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export_off.h>


class MyFace;
class MyVertex;

struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>::AsVertexType,    vcg::Use<MyFace>::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face < MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::Color4b, vcg::face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};
int main(int argc, char **argv)  {
    printf("prova");
    // definition of type MyMesh
    MyMesh m;
    if(vcg::tri::io::ImporterOFF<MyMesh>::Open(m,argv[1])!=vcg::tri::io::ImporterOFF<MyMesh>::NoError)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }
  
    // Calculate the number of vertices, edges, and faces
    int numVertices = m.VN();
    int numFaces = m.FN();
    

    // Print the results
    printf("Number of vertices: %d\n", numVertices);
    printf("Number of faces: %d\n", numFaces);
   

  vcg::tri::io::ExporterOFF<MyMesh>::Save(m,"test.off");
  return 0;
}
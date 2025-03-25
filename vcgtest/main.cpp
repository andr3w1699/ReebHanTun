#include<vcg/complex/complex.h>
#include <wrap/io_trimesh/import.h>    // reads mesh files (e.g. .off, .ply formats)
#include<wrap/io_trimesh/export_off.h> // exports a mesh in .off format 
#include<wrap/io_trimesh/export_ply.h> // exports a mesh in .ply format


#include "reebhantun_wrapper.h"

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

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b, vcg::vertex::Qualityf, vcg::vertex::BitFlags, vcg::vertex::VEAdj >{};
class MyEdge : public vcg::Edge<MyUsedTypes,vcg::edge::VertexRef> {};
class MyFace  : public vcg::Face < MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::Color4b, vcg::face::BitFlags, vcg::face::EFAdj > {};

// MyMesh class represents the entire mesh. It inherits from vcg::tri::TriMesh and uses std::vector 
// to store collections of MyVertex, MyEdge, and MyFace elements. 
class MyMesh  : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyEdge>, std::vector<MyFace> > {};




int main(int argc, char **argv)  {
    // declaration of an object of type MyMesh -- VCG lib
    MyMesh m_vcg;

    // load a mesh with VCG lib
    if(vcg::tri::io::ImporterOFF<MyMesh>::Open(m_vcg,argv[1])!=vcg::tri::io::ImporterOFF<MyMesh>::NoError)
      {
        printf("Error reading file  %s\n",argv[1]);
        exit(0);
      }
  
    // initialize ReebHanTun 
    // instatiate an object ReebHanTunWrapper by passing 
    // the mesh you want to compute the basis of loops 
    vcg::tri::ReebHanTunWrapper<MyMesh> wrapper;
    // declare a mesh to store the computed basis 
    MyMesh loops;
    // compute the basis and store the results in the mesh loops (all the basis tunnel+handle is computed)
    wrapper.ComputeBasis(m_vcg, loops);

    // use the provided getter functions to print the number of computed loops 
    printf("Number of handle loops found: %d\n", wrapper.GetNumHandleLoops());
    printf("Number of Tunnel loops found: %d\n", wrapper.GetNumTunnelLoops());


    vcg::Point3<RealTypeForVector3> reebGraphDirection = wrapper.GetStoredDirection();
    std::cout << "Direction used to compute ReebGraph: (" << reebGraphDirection.X() << ", " << reebGraphDirection.Y() << ", " << reebGraphDirection.Z() << ")" << std::endl;

    // prepare two mesh to store the handle and tunnel loops 
    MyMesh h_loops;
    MyMesh v_loops;
    
    // get the last handle loop in the basis 
    // or you can specify an index of the loop you want
    // provided it belongs to the permitted range   
    wrapper.GetHandleLoops(m_vcg, h_loops, wrapper.GetNumHandleLoops()-1);
    // get all the basis of tunnel loops by specifying -1 as index 
    wrapper.GetTunnelLoops(m_vcg, v_loops, -1);

    MyMesh rht_mesh;
    wrapper.GetColoredMeshByScalarField(rht_mesh);

    vcg::tri::io::ExporterPLY<MyMesh>::Save(h_loops,"handle_loops.ply", vcg::tri::io::Mask::IOM_EDGEINDEX); 
    vcg::tri::io::ExporterPLY<MyMesh>::Save(v_loops,"tunnel_loops.ply", vcg::tri::io::Mask::IOM_EDGEINDEX);
    vcg::tri::io::ExporterPLY<MyMesh>::Save(loops,"loops.ply", vcg::tri::io::Mask::IOM_EDGEINDEX);
    vcg::tri::io::ExporterPLY<MyMesh>::Save(rht_mesh, "ReebHanTun_Mesh_with_color.ply", vcg::tri::io::Mask::IOM_VERTCOLOR);
    return 0;
}

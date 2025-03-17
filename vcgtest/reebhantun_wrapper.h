#ifndef REEBHANTUN_WRAPPER_H
#define REEBHANTUN_WRAPPER_H

#include <psbmReebGraph.h>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include "SimpleMesh.h"
#include <FilesOutputForOptimalCycles.h>


///
/** \addtogroup trimesh */
/*@{*/
/// Wrapper to use ReebHanTun code inside the vcg library. 
/// Non-static class parametric w.r.t. to MeshType that acts 
/// as a wrapper to run ReebHanTun code inside the vcg library.
/// 1. Invoke the constructor by passing a reference to a mesh MeshType
///    for which we want to compute the basis
/// 2. Invoke the method ComputeBasis by passing a reference to a mesh MeshType 
///    where you want the loop to be saved
/// 
template <class MeshType>
class ReebHanTunWrapper
{
    

private:
    // typedef of parametric type VertexType & FaceType
    typedef typename MeshType::VertexType           VertexType;  
    typedef typename MeshType::FaceType             FaceType;

    // mesh for which we want to compute the basis -- VCG format 
    // reference to an object of parametric type MeshType
    MeshType& vcg_mesh;

    // mesh for which we want to compute the basis in ReebHanTun format, object of type _SimpleMesh 
    _SimpleMesh m_rht;

    // private variables where the tunnel & handle loops are saved 
    // after being computed bu the function ComputeBasis
    std::vector<std::set<int> > v_basis_loops;
    std::vector<std::set<int> > h_basis_loops;

    // flag set to 1 when the loops has been computed for the current mesh 
    // initially false
    bool is_basis_computed = false;  


    // define the typedef for struct Params 
    typedef struct Params
    {
        // enlarge/scaling factor used internally by ReebHanTun
        const float enlarge_factor;

        // enlarge_factor is set to 10000
        Params() : enlarge_factor(10000.0f) {}
    } Params;

    // declare the private variable of type Params
    Params _params;

    // define the typedef for struct Stats
    typedef struct Stats
    {
        // num of handle_loops
        int handle_num;
        // num of tunnel_loops
        int tunnel_num;

        // constructor method to initialize the fields to -1
        Stats() : handle_num(-1), tunnel_num(-1) {}
    } StatsType;

    // declare the private variable of type Stats
    StatsType _stats;

/* 
MeshConverter is function to convert a VCG mesh type MyMesh in a ReebHanTun mesh _SimpleMesh. 
It also calculates the bounding box (minBd, maxBd), normals (meshNormal), and triangle orientations (OrientTriangles). 
The fEnlargeFactor is used to scale the vertices.
*/
void MeshConverter (_SimpleMeshVertex &minBd, _SimpleMeshVertex &maxBd, const MeshType & vcg_mesh,  _SimpleMesh & rht_mesh, std::vector<Vector3> &meshNormal, std::vector<int> &OrientTriangles, const float fEnlargeFactor ) {
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
        const VertexType &v = vcg_mesh.vert[i];
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
        const FaceType &f = vcg_mesh.face[i];
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
    
    public:

    // Constructor that accepts a reference to a MeshType object for which the loops are computed
    ReebHanTunWrapper(MeshType &meshRef) : vcg_mesh(meshRef) {}
        
    /**
     * @brief SetBaseMesh
     * @param mesh
     * define the base mesh for which the loops are computed 
     * 
     */
    void SetBaseMesh(MeshType &mesh) {
        // reinitialize internal structures 
        // preparing to compute a basis for a new mesh 
        this->~ReebHanTunWrapper(); // // Manually call destructor
        new (this) ReebHanTunWrapper(mesh);  // Placement new: reconstruct the object
    }
    


    /**
     * @brief ComputeBasis
     * @param MeshType &loops : reference to a mesh where the computed basis is going to be saved 
     * Function to invoke to compute a basis of handle and tunnel loops on the current mesh 
     */
   void ComputeBasis(MeshType &loops) {
        Vector3 distinctDirection;
        // vectors of sets that will store basis loops for handles and tunnels
        psbmReebGraph reebGraph;
        std::vector<Vector3>   meshNormal;
        std::vector<int> OrientTriangles;
        std::set<int> extraVertices;
        int nOrgTriangleSize = 0;
        double BoundingBoxRadius;
        int genus = 0;
        _SimpleMeshVertex minBd;
        _SimpleMeshVertex maxBd;
    
        // convert the vcg mesh in reebhantun _SimpleMesh
        MeshConverter(minBd, maxBd, vcg_mesh, m_rht, meshNormal, OrientTriangles, _params.enlarge_factor);
    
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
            CycleLocalOptimization(m_rht, reebGraph, OrientTriangles, v_basis_loops, h_basis_loops, 1.f / _params.enlarge_factor);
        else
            CycleLocalOptimization_bdry(m_rht, reebGraph, OrientTriangles, v_basis_loops, h_basis_loops, extraVertices, 1.f / _params.enlarge_factor);
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
        files_out_op.WriteGeomviewListFormat(OutputFileName.c_str(), v_basis_loops, h_basis_loops, OrientTriangles, orgVertexSize, nOrgTriangleSize, 1.f / _params.enlarge_factor);
    
        //MyMesh test;
        //ReverseMeshConverter(m_rht, test);
        //vcg::tri::io::ExporterOFF<MyMesh>::Save(test,"test.off");
        //std::cout << "num edges: " << m_rht.vecEdge.size();
        //std::cout << std::endl; 
    
        for (const auto& s : v_basis_loops) {
            for (const auto& elem : s) {
                std::cout << elem << " ";
            }
            std::cout << std::endl;  // new line after each set
        }
    
        for (const auto& s : h_basis_loops) {
            for (const auto& elem : s) {
                std::cout << elem << " ";
            }
            std::cout << std::endl;  // new line after each set 
        }
        
        // now the basis has been computed for the current mesh
        is_basis_computed = true;
        // Update the Stats struct with the number of loops
        _stats.handle_num = h_basis_loops.size();
        _stats.tunnel_num = v_basis_loops.size();
        
        for (const auto& s : v_basis_loops) {
            for (const auto& elem : s) {
                if(elem > 0 && elem < m_rht.vecEdge.size())
                {
                    int v0 = m_rht.vecEdge[elem].v0;
                    int v1 = m_rht.vecEdge[elem].v1;         
                    vcg::tri::Allocator<MeshType>::AddEdge(loops, vcg_mesh.vert[v0].P(),vcg_mesh.vert[v1].P());
                }
            }
        }
        
        
        for (const auto& s : h_basis_loops) {
            for (const auto& elem : s) {
                if(elem > 0 && elem < m_rht.vecEdge.size())
                {
                    int v0 = m_rht.vecEdge[elem].v0;
                    int v1 = m_rht.vecEdge[elem].v1;         
                    vcg::tri::Allocator<MeshType>::AddEdge(loops, vcg_mesh.vert[v0].P(),vcg_mesh.vert[v1].P());
                }
            }
        }
        
        
    }

    
    /**
     * @brief GetNumHandleLoops
     * getter method to retrieve the number of tunnel loops computed for the current mesh 
     * return -1 if the basis hasn't been computed yet 
     */
    int GetNumHandleLoops() {
        return _stats.handle_num;
    }

    /**
     * @brief GetNumTunnelLoops
     * getter method to retrieve the number of tunnel loops computed for the current mesh
     * return -1 if the basis hasn't been computed yet 
     */
    int GetNumTunnelLoops() {
        return _stats.tunnel_num;
    }


    /**
     * @brief GetHandleLoops
     * @param loop the mesh to contain the loop(s)
     * @param index the index of the required loop (default is -1 for all loops in a single mesh)
     * 
     * Get the computed loop
     */
    void GetHandleLoops(MeshType &loop, int index=-1)
    {
        if(is_basis_computed==false)
            std::cout << "BASIS NOT COMPUTED FOR THE CURRENT MESH: INVOKE COMPUTEBASIS"  << std::endl;
        else if(index==-1) {
        for (const auto& s : h_basis_loops) {
            for (const auto& elem : s) {
                if(elem > 0 && elem < m_rht.vecEdge.size())
                {
                    int v0 = m_rht.vecEdge[elem].v0;
                    int v1 = m_rht.vecEdge[elem].v1;         
                    vcg::tri::Allocator<MeshType>::AddEdge(loop, vcg_mesh.vert[v0].P(),vcg_mesh.vert[v1].P());
                }
            }
        }
      }
      else if(index>= 0 && index<h_basis_loops.size()) {
        std::set<int> s = h_basis_loops[index];
        for (const auto& elem : s) {
            if(elem > 0 && elem < m_rht.vecEdge.size())
            {
                int v0 = m_rht.vecEdge[elem].v0;
                int v1 = m_rht.vecEdge[elem].v1;         
                vcg::tri::Allocator<MeshType>::AddEdge(loop, vcg_mesh.vert[v0].P(),vcg_mesh.vert[v1].P());
            }
        }
    }
    else std::cout << "LOOP INDEX OUT OF BOUND"  << std::endl;
    }
    
    void GetTunnelLoops(MeshType &loop, int index=-1)
    {   
        if(is_basis_computed==false)
        std::cout << "BASIS NOT COMPUTED FOR THE CURRENT MESH: INVOKE COMPUTEBASIS"  << std::endl;
        else if(index==-1) {
            for (const auto& s : v_basis_loops) {
                for (const auto& elem : s) {
                    if(elem > 0 && elem < m_rht.vecEdge.size())
                    {
                        int v0 = m_rht.vecEdge[elem].v0;
                        int v1 = m_rht.vecEdge[elem].v1;         
                        vcg::tri::Allocator<MeshType>::AddEdge(loop, vcg_mesh.vert[v0].P(),vcg_mesh.vert[v1].P());
                    }
                }
            }         
          }
        else if(index>= 0 && index<v_basis_loops.size()) {
            std::set<int> s = v_basis_loops[index];
            for (const auto& elem : s) {
                if(elem > 0 && elem < m_rht.vecEdge.size())
                {
                    int v0 = m_rht.vecEdge[elem].v0;
                    int v1 = m_rht.vecEdge[elem].v1;         
                    vcg::tri::Allocator<MeshType>::AddEdge(loop, vcg_mesh.vert[v0].P(),vcg_mesh.vert[v1].P());
                }
            }
        }
        else std::cout << "LOOP INDEX OUT OF BOUND"  << std::endl;
    }
};



#endif // REEBHANTUN_WRAPPER_H

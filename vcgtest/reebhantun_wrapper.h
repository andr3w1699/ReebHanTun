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
/// 
template <class MeshType>
class ReebHanTunWrapper
{
    

private:
    _SimpleMesh _mesh;
    struct Params
    {
        float enlarge_factor=10000.0f;
    };

    struct Stats
    {
        int handle_num;
        int tunnel_num;        
    };
    
    public:
    typedef typename MeshType::VertexType           VertexType;
        
        /**
     * @brief SetBaseMesh
     * @param mesh
     * define the base mesh for which the loops are computed 
     * 
     */
    void SetBaseMesh(MeshType &mesh) {
        _mesh = &mesh;
    }
    
    /**
     * @brief ComputeBasis
     * @param mesh
     * define the base mesh for which the loops are computed
     */
    void ComputeBasis(Stats &s)
    {
        
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
        
    }
    
    void GetTunnelLoops()
    {
        
    }
};


#endif // REEBHANTUN_WRAPPER_H

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
from numba import jit, njit



def sphere_delaunay(r_sphere):
    """
        Cretates a triangulation of the spherical surface 
        by using the Convex Hull.
        
        Input:  
               r_sphere  = np.array[n_vertices, 3]. Coordinates of the points
               
        Output:  
               tri  = np.array[n_tris, 3]. Triangulation array
               
    """
    hull = ConvexHull(r_sphere)
    indices = hull.simplices

    ntri = indices.size//3
    tri = np.zeros((ntri, 3), dtype=np.int32)

    for i in range(0,ntri):
        tri[i,:] = indices[i] 

    return tri



def triplot(r0, tri):
    """
        Data reppresentation for easy visualisation of a triangulated surface.
        
        Input:  
               r0  = np.array[n_vertices, 3]. Coordinates of the points
               tri  = np.array[n_tris, 3]. Triangulation array
        
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_trisurf(r0[:,0], r0[:,1], r0[:,2], 
                    triangles=tri, 
                    color=(0.5,0.5,0.5,1.0), 
                    edgecolor=(0.0,0.0,0.0, 1.0),
                    linewidth=0.1,
                    antialiased=True,
                    shade=False)
    plt.show()



@njit
def topo_defects(tri, Np):
    """
        Calculates coordianation number and from it
        the number of five-fold and seven fold defects
        as well as the total toplogical charge
        
        Input:  
               tri  = np.array[n_tris, 3]. Triangulation array
               Np   = int. Total number of particles defining the surface 
    """
    N_tri = len(tri)
    cn = np.zeros(Np)
    for t in range(0, N_tri):
        v1 = tri[t,0]
        v2 = tri[t,1]
        v3 = tri[t,2]
        cn[v1] += 1
        cn[v2] += 1 
        cn[v3] += 1

    # Five-fold and seven fold defects
    n5=0
    n7=0
    q=0
    for i in range(0,Np):
        if cn[i] == 5:
            n5 += 1
        if cn[i] == 7:
            n7 += 1
        q += 6-cn[i]
    q=int(q)
    topology = [n5, n7, q]
    
    return topology



@jit
def area_variability(r0, tri):
    """
        Calculate the area of each a triangle and the resulting  
        dispersion in area per triangle
        
        Input:  
               tri  = np.array[n_tris, 3]. Triangulation array
               Np   = int. Total number of particles defining the surface 
               
        Output:
                std(Area_triangles)/mean(Area_triangles)
    """
    
    N_tri = len(tri)
    areas = np.zeros(N_tri)
    for t in range(0, N_tri):
        v1 = tri[t,0]
        v2 = tri[t,1]
        v3 = tri[t,2]
        r21 = r0[v2,:] - r0[v1,:]
        r31 = r0[v3,:] - r0[v1,:]
        areas[t] = np.linalg.norm(np.cross(r21, r31))/2
    return np.std(areas)/np.mean(areas)



@jit
def check_orient(r0, tri):
    """
        For a convex polyhedron, we can change the orientation
        such that we have homogeneous orientation everywhere on
        the surface.
        
        Input:  
               r0  = np.array[n_vertices, 3]. Coordinates of the points
               tri  = np.array[n_tris, 3]. Triangulation array
               
        Ouput:
            tri  = np.array[n_tris, 3]. Triangulation array with coherent orientation
            orientation = int. Indiation of the correct orientation in the original mesh
        
    """
    N_tri = len(tri)
    Np = len(r0)

    # Calculate normal vectors to the surface
    oriented = 1
    DELTA = 0.3	
    nv = np.zeros((N_tri,3))
    for t in range(0, N_tri):
        v1 = tri[t,0]
        v2 = tri[t,1]
        v3 = tri[t,2]
        r21 = r0[v2,:] - r0[v1,:]
        r31 = r0[v3,:] - r0[v1,:]
        tnv = np.cross(r21, r31)
        nv[t,:] = tnv
        nv[t,:] /= np.linalg.norm(nv[t,:])
        
        # Move along the normal direction and check change in radius
        mp = (r0[v1,:] + r0[v2,:] + r0[v3,:])/3.0
        tmp = mp + DELTA*nv[t,:]        
        tr = np.linalg.norm(tmp)

        if tr >= 1.0:
            oriented = 0
            tri[t,1] = tri[t,1] + tri[t,2]
            tri[t,2] = tri[t,1] - tri[t,2]
            tri[t,1] = tri[t,1] - tri[t,2]
        
    return tri, oriented


def export_obj(r, tri):
    """
        Export the coordinates and the mesh as obj fileobj
        
         Input:  
            r0  = np.array[n_vertices, 3]. Coordinates of the points
            tri  = np.array[n_tris, 3]. Triangulation array
    """
    fileobj = open("./outfiles/init_coords.obj","w")
    for i in range(0, len(r)):
        fileobj.writelines("v\t" + str(r[i,0]) + "\t" + str(r[i,1]) + "\t" + str(r[i,2]) + "\n" ) 
    for t in range(0, len(tri)):
        fileobj.writelines("f\t" + str(tri[t,0]+1) + "\t" + str(tri[t,1]+1) + "\t" + str(tri[t,2]+1) + "\n" ) 
        
    fileobj.close() 
    
    
if __name__ == "__main__":
    
    SPHERE_COORD_FNAME = "./outfiles/init_coords.dat"
    MESH_FNAME = "./outfiles/mesh.dat"
    
    r0 = np.loadtxt(SPHERE_COORD_FNAME, delimiter="\t") # Coordinates minim
    tri = sphere_delaunay(r0)
    
    topology = topo_defects(tri, len(r0))
    relserr = area_variability(r0, tri)
    tri_oriented, oriented = check_orient(r0, tri)
    np.savetxt(MESH_FNAME, tri_oriented.astype(int), fmt='%i', delimiter="\t")
    export_obj(r0, tri_oriented)
    
    print(" - Five-fold defects, n5 = " + str(topology[0]))
    print(" - Seven-fold defects, n7 = " + str(topology[1]))
    print(" - Total topological charge, q = " + str(topology[2]) )
    print(" - Triangle area variability = " + str(relserr*100) + "%")
    
    if oriented == 0:
        print(" - Surface initially not oriented. Corrected.")
        
    if( topology[2] != 12 or relserr > 0.05):
        print(" - Recommended increase in the number of Monte Carlo steps.")
    
    triplot(r0,tri)


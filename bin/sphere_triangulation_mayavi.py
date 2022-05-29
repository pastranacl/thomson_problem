import numpy as np
import mayavi.mlab as mlab

SPHERE_COORD_FNAME = "rsphere_eq.dat"
MESH_FNAME = "mesh.dat"


# Reppresentation of the data
def sphere_plot(r0, tri):
        
    x=r0[:,0]
    y=r0[:,1]
    z=r0[:,2]
    s=r0[:,1]
    Np = len(r0)
    s= np.random.rand(Np)
    #mlab.triangular_mesh(x, y, z, triangles, scalars=s)
    
    mlab.triangular_mesh(x, y, z, tri, scalars=s, transparent=True)  
    mlab.triangular_mesh(x, y, z, tri, color=(0,0,0), representation='fancymesh', transparent=True)  
    mlab.orientation_axes()
    mlab.show() 




if __name__ == "__main__":
    r0 = np.loadtxt(SPHERE_COORD_FNAME, delimiter="\t") # Coordinates minim
    tri = np.loadtxt(MESH_FNAME, delimiter="\t") # Triangles
    sphere_plot(r0,tri)

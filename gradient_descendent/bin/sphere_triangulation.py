# Cesar L. Pastrana, 2021
# Triangulation of the sphere using the Convex hull and reppresentation

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull

SPHERE_COORD_FNAME = "minim_coords.dat"
MESH_FNAME = "mesh.dat"


# Triangulation and data saving
def sphere_delaunay():
    r_sphere = np.loadtxt(SPHERE_COORD_FNAME, delimiter="\t") # Coordinates minim
    hull = ConvexHull(r_sphere)
    indices = hull.simplices

    ntri = indices.size/3
    tri = np.zeros((ntri, 3), dtype = int)

    for i in range(0,ntri):
        tri[i,:] = indices[i] 

    np.savetxt(MESH_FNAME, tri.astype(int), fmt='%i', delimiter="\t")


# Reppresentation of the data
def sphere_plot():

    tri = np.loadtxt(MESH_FNAME, delimiter="\t") # Triangles
    r0 = np.loadtxt(SPHERE_COORD_FNAME, delimiter="\t") # Coordinates minim
    x0=r0[:,0]
    y0=r0[:,1]
    z0=r0[:,2]
    fig = plt.figure()
    ax = fig.gca(projection='3d')


    ax.plot_trisurf(x0, y0, z0, triangles=tri, 
                    color=(0.5,0.5,0.5,1.0), 
                    edgecolor=(0.0,0.0,0.0, 1.0),
                    linewidth=0.1,
                    antialiased=True,
                    shade=False)
    plt.show()

# Calculates the coordination number
def topo_defects():
    tri = np.loadtxt(MESH_FNAME, delimiter="\t") 
    r0 = np.loadtxt(SPHERE_COORD_FNAME, delimiter="\t") 
    N_tri = len(tri)
    Np = len(r0)
    cn = np.zeros(Np)
    for t in range(0, N_tri):
        v1 = int(tri[t,0])
        v2 = int(tri[t,1])
        v3 = int(tri[t,2])
        cn[v1] = cn[v1]+1
        cn[v2] = cn[v2]+1
        cn[v3] = cn[v3]+1

    # Five-fold and seven fold defects
    n5=0
    n7=0
    for i in range(0,Np):
        if cn[i] == 5:
            n5=n5+1
        if cn[i] == 7:
            n7=n7+1

    topology = [n5, n7, n5 - n7]
    return topology	

    print "Five-fold defects: " + str(n5)
    print "Seven-fold defects: " + str(n7)

# Calculate the area of each a triangle and the resulting dispersion in area
def area_variability():
    tri = np.loadtxt(MESH_FNAME, delimiter="\t") 
    r0 = np.loadtxt(SPHERE_COORD_FNAME, delimiter="\t") 
    N_tri = len(tri)
    areas = np.zeros(N_tri)
    for t in range(0, N_tri):
        v1 = int(tri[t,0])
        v2 = int(tri[t,1])
        v3 = int(tri[t,2])
        r21 = r0[v2,:] - r0[v1,:]
        r31 = r0[v3,:] - r0[v1,:]
        areas[t] = np.linalg.norm(np.cross(r21, r31))/2
        
    return np.std(areas)/np.mean(areas)

if __name__ == "__main__":
    sphere_delaunay()
    topolgy = topo_defects();
    relserr = area_variability()

    print "Five-fold defects, n5 = " + str(topolgy[0])
    print "Seven-fold defects, n7 = " + str(topolgy[1])
    print "Total topological charge, q = " + str(topolgy[2])
    print "Triangle area variability = " + str(relserr*100) + "%"

    if( topolgy[2] != 12 or relserr > 0.05):
        print "Recommended increase in the number of Monte Carlo steps "
    sphere_plot()



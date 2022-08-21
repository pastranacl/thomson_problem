# Cesar L. Pastrana, 2021
# Triangulation of the sphere using the Convex hull and reppresentation

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull



# Triangulation and data saving
def sphere_delaunay(r_sphere):

    hull = ConvexHull(r_sphere)
    indices = hull.simplices

    ntri = indices.size//3
    tri = np.zeros((ntri, 3))

    for i in range(0,ntri):
        tri[i,:] = indices[i] 

    return tri


# Reppresentation of hte data
def sphere_plot(r0, tri):
        
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
def topo_defects(r0, tri):

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
            n5 += 1
        if cn[i] == 7:
            n7 += 1

    topology = [n5, n7, n5 - n7]
    return topology	

    print("Five-fold defects: " + str(n5))
    print("Seven-fold defects: " + str(n7))



# Calculate the area of each a triangle and the resulting dispersion in area
def area_variability(r0, tri):
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



# Calculate orientation of the surface and exchange the order to have homogeneous inward orientation
def check_orient(r0, tri):
    
    N_tri = len(tri)
    Np = len(r0)

    # Scale to unity
    radius = 0
    for i in range(0,Np):
        radius += np.linalg.norm(r0[i,:])
    radius /= Np
    r0 /= radius
    
    # Calculate normal vectors to the surface
    oriented = 1
    DELTA = 0.3	
    nv = np.zeros((N_tri,3))
    for t in range(0, N_tri):
        v1 = int(tri[t,0])
        v2 = int(tri[t,1])
        v3 = int(tri[t,2])
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
    
    if oriented == 0:
        np.savetxt(MESH_FNAME, tri.astype(int), fmt='%i', delimiter="\t")
        
    return oriented


if __name__ == "__main__":
    
    SPHERE_COORD_FNAME = "rsphere_eq.dat"
    MESH_FNAME = "mesh.dat"
    
    r0 = np.loadtxt(SPHERE_COORD_FNAME, delimiter="\t") # Coordinates minim
    tri = sphere_delaunay(r0)
    np.savetxt(MESH_FNAME, tri.astype(int), fmt='%i', delimiter="\t")

    topolgy = topo_defects(r0, tri)
    relserr = area_variability(r0, tri)
    oriented = check_orient(r0, tri)

    print("Five-fold defects, n5 = " + str(topolgy[0]))
    print("Seven-fold defects, n7 = " + str(topolgy[1]))
    print("Total topological charge, q = " + str(topolgy[2]) )
    print("Triangle area variability = " + str(relserr*100) + "%")
    
    if oriented == 0:
        print("Surface initially not oriented. Corrected.")
    else:
        print("Surface is oriented")
        
    if( topolgy[2] != 12 or relserr > 0.05):
        print("Recommended new minimisation")
    #sphere_plot(r0,tri)
    r0 = np.loadtxt(SPHERE_COORD_FNAME, delimiter="\t") # Coordinates minim
    tri = np.loadtxt(MESH_FNAME, delimiter="\t") # Triangles
    nv = np.zeros((len(tri), 3))
    cm_tri = np.zeros((len(tri), 3))
    for t in range(0, len(tri)):
        v1= int(tri[t,0])
        v2= int(tri[t,1])
        v3= int(tri[t,2])
        r21 = r0[v2,:]-r0[v3,:]
        r31 = r0[v3,:]-r0[v1,:]
        nv[t,:] = np.cross(r21,r31);
        nv[t,:] /= np.linalg.norm(nv[t,:])
        cm_tri[t,:] = (r0[v1,:] + r0[v2,:]+ r0[v3,:])/3
        
    # Plot
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_trisurf(r0[:,0],
                    r0[:,1], 
                    r0[:,2],
                    triangles=tri, 
                    color=(0.5,0.5,0.5,0.8), 
                    edgecolor=(0.0,0.0,0.0, 1.0),
                    linewidth=0.2,
                    antialiased=True,
                    shade=False)

    plt.axis('off')
    plt.savefig("n" + str(len(r0)) + ".pdf", dpi=100, bbox_inches='tight')
    plt.show()

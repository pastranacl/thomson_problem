import numpy as np
import jax
import jax.numpy as jnp
import tqdm

@jax.jit
def energy(r, i_idx, j_idx):
    """
        NxN energy of the particles to be executed on the GPU
    """
    locE = jnp.where(i_idx==j_idx, 0,  1./jnp.linalg.norm(r[:,None] - r, axis=2) )
    return jnp.sum(locE)/2



def mcsphere(npart, drad, mc_steps, init_temp, anneal_time):
    """
        Generates a collection of particles homogeneously 
        distributed on the surface of a spherical shell.
        
        Input:
            npart       = int. Number of particles
            drad        = float. Maximum angular displacement of the Monte Carlo moves.
            mc_steps    = int. Number of Monte carlo steps to executed
            init_temp   = float. Initial temperature
            annel_time  = int. Annealing time (steps between temperature changes)
        
        Output:
            r = np.array [npart,3]. Cartesian coordinates of the    
                vertices defining the particles on the sphercal shell

    """
    kB = 1.38064852e-23
    T = init_temp;
    dT = (T-1.0)/(mc_steps/anneal_time)
    BETA = 1.0/(kB*T)

    r = np.zeros((npart, 3), dtype=np.float32)
    tr = np.zeros(3, dtype=np.float32)
    azi = np.zeros(npart, dtype=np.float32)
    polar = np.zeros(npart, dtype=np.float32)

    # We need these for the calculation of the energy (pair-pair)
    i_idx = jnp.arange(npart)[:,None]
    j_idx = jnp.arange(npart)[:]
    jax.vmap(energy)
             
    # Initial random positioning of the particles in the surface of the sphere 
    for i in range(0, npart):
        azi[i] = 2*np.pi*np.random.rand()
        polar[i] = np.pi*np.random.rand()
        r[i,0] = np.sin(polar[i])*np.cos(azi[i])
        r[i,1] = np.sin(polar[i])*np.sin(azi[i])
        r[i,2] = jnp.cos(polar[i])
        
    # Initial configuration energy
    E_curr = energy(r, i_idx, j_idx)

    # MAIN LOOP
    for n in range(1, mc_steps):

        if n == anneal_time:
            T -= dT
            BETA = 1.0/(kB*T);

        # Select a random particle and atempt a random move 
        tp = np.random.randint(0, npart)
        
        # Attempted displacement
        tpolar = polar[tp]
        tazi   = azi[tp]
        tr = np.copy(r[tp,:])

        polar[tp] += drad*(np.random.rand()-0.5);
        azi[tp] += drad*(np.random.rand()-0.5);
        r[tp,0] = np.sin(polar[tp])*np.cos(azi[tp])
        r[tp,1] = np.sin(polar[tp])*np.sin(azi[tp])
        r[tp,2] = jnp.cos(polar[tp])

        # Calculates the energy of the new configuration 
        # applied Metropolis criterion 
        tE = energy(r, i_idx, j_idx)
        DE = tE - E_curr;
        if np.random.rand() <= np.exp(-DE*BETA):
            E_curr = tE;
        else:
            r[tp,:] = np.copy(tr)
            azi[tp] = tazi
            polar[tp] = tpolar
        
        # Adjust the steps
        if mc_steps//2:
            drad /= 2
        if 3*mc_steps//4:
            drad /= 2
        
    return r




if __name__ == "__main__":
    NP = 8_000                      # Number of particles
    DRAD = 2.0*np.pi/180.0          # Maximum anguar displacement
    MC_STEPS = 6_000_000            # Number of Monte-Carlo steps
    T0 = 300                        # Initial temperature
    ANNEAL_TIME = 50                # Annealing time (steps between temperature changes)
    
    r = mcsphere(NP, DRAD, MC_STEPS, T0, ANNEAL_TIME) 
    np.savetxt("init_coords.dat", r, delimiter="\t")
   
   
   
   

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import matplotlib.animation as animation
import matplotlib as mpl
import math

# Hydro Modules
import SPHAdvance as hadv
import SPHSetup as hset


# Advance the simulation
def advance(n, r, v, gtime, tstep, t, lim, counter, press, dens, temps, allxs, allys, allvmags, alldens, alltemps, framenum):
    press, dens, temps, r, v, counter = hadv.hydrosim(r, v, temps, gtime, tstep, counter)
    t = t + gtime

    # Storing the data needed for generating an animation
    for i in range(n):
        allxs[framenum, i] = r[0,i]
        allys[framenum, i] = r[1,i]
        alldens[framenum, i] = dens[i]
        alltemps[framenum, i] = temps[i]
        allvmags[framenum, i] = math.sqrt((v[0,i]**2) + (v[1,i]**2) + (v[2,i]**2))

    return r, v, t, counter, press, dens, temps, allxs, allys, allvmags, alldens, alltemps
    
    

# Setting up the simulation
print("Enter the number of particles to be simulated:")
n = int(input())

t = 0
counter = 1
dframes = []
tframes = []
vframes = []

print("Enter the total amount of time that should be simulated (in seconds):")
lim = float(input())
print("Enter the amount of simulated time between successive graphs (in seconds):")
gtime = float(input())
print("Enter the timestep length (in seconds):")
tstep = float(input())
print("Enter an initial temperature (in Kelvin):")
temps = float(input())
temps = np.resize(temps, n)

print("""Current Initial Configurations:
1 - Sphere (zero initial velocity for all particles)
2 - Sphere (with turbulent velocity field)
3 - Cuboid (zero initial velocity for all particles)
""")
print("Enter the ID of an initial configuration from the list above:")
boxform = int(input())

if boxform==1 or boxform==2:
    print("Enter the radius of the sphere (in metres):")
    rad = float(input())
    print("Generating initial conditions...")
    if boxform==1:
        r, v = hset.sprsetup(n, 3, 0, rad)
    else:
        r, v = hset.sprsetup(n, 3, 1, rad)

elif boxform==3: # Using elif to allow more boxforms to be added easily in the future
    print("Enter the size of the cuboid along the x-axis (in metres):")
    xlim = float(input())
    print("Enter the size of the cuboid along the y-axis (in metres):")
    ylim = float(input())
    print("Enter the size of the cuboid along the z-axis (in metres):")
    zlim = float(input())
    print("Generating initial conditions...")
    r = hset.cubsetup(n, lim=[xlim, ylim, zlim], ndims=3)
    v = np.zeros(shape=np.shape(r))

print("Initial conditions generated.")
nframes = int(lim/gtime) # Assuming there is not a frame for the initial state at t=0
                         # This assumption is subject to change in the future
                         # For now there is no frame at t=0 because the densities
                         # are not calculated during the setup

allxs = np.zeros(shape=(nframes, n))
allys = np.zeros(shape=(nframes, n))
allvmags = np.zeros(shape=(nframes, n))
alldens = np.zeros(shape=(nframes, n))
alltemps = np.zeros(shape=(nframes, n))
press = np.zeros(shape=n)
dens = np.zeros(shape=n)

print("Beginning Simulation")
for framenum in range(nframes):
    r, v, t, counter, press, dens, temps, allxs, allys, allvmags, alldens, alltemps = advance(n, r, v, gtime, tstep, t, lim, counter, press, dens, temps,
                                                                                              allxs, allys, allvmags, alldens, alltemps, framenum)
    print("Simulation continuing.", ((t/lim)*100), "% complete")

print("Generating Animated Figures")

# Animation
dfig, dax = plt.subplots()
tfig, tax = plt.subplots()
vfig, vax = plt.subplots()
dmin = np.min(alldens)
dmax = np.max(alldens)
tmin = np.min(alltemps)
tmax = np.max(alltemps)
vmin = np.min(allvmags)
vmax = np.max(allvmags)
dnorm = cols.Normalize(dmin, dmax)
tnorm = cols.Normalize(tmin, tmax)
vnorm = cols.Normalize(vmin, vmax)
cmap = mpl.colormaps['viridis']

for i in range(nframes):
    dframe = dax.scatter(allxs[i,:], allys[i,:], c=alldens[i,:], norm=dnorm, cmap=cmap.resampled(256))
    tframe = tax.scatter(allxs[i,:], allys[i,:], c=alltemps[i,:], norm=tnorm, cmap=cmap.resampled(256))
    vframe = vax.scatter(allxs[i,:], allys[i,:], c=allvmags[i,:], norm=vnorm, cmap=cmap.resampled(256))
    dframes.append([dframe])
    tframes.append([tframe])
    vframes.append([vframe])
    
dfig.colorbar(dframe, ax=dax, label=r'Density [kg $m^{-3}$]')
tfig.colorbar(tframe, ax=tax, label='Temperature [K]')
vfig.colorbar(vframe, ax=vax, label=r'Magnitude of Velocity [m $s^{-2}$]')
dax.set_xlabel('Position along x-axis')
dax.set_ylabel('Position along y-axis')
tax.set_xlabel('Position along x-axis')
tax.set_ylabel('Position along y-axis')
vax.set_xlabel('Position along x-axis')
vax.set_ylabel('Position along y-axis')
dfig.suptitle('Density Map')
tfig.suptitle('Temperature Map')
vfig.suptitle('Velocity Map')
anidens = animation.ArtistAnimation(fig=dfig, artists=dframes, repeat=False)
anitemp = animation.ArtistAnimation(fig=tfig, artists=tframes, repeat=False)
anivels = animation.ArtistAnimation(fig=vfig, artists=vframes, repeat=False)
anidens.save(filename="DensityMap.gif", writer="pillow")
anitemp.save(filename="TemperatureMap.gif", writer="pillow")
anivels.save(filename="VelocityMap.gif", writer="pillow")
plt.show()
    


    



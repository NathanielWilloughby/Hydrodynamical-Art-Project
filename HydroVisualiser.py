import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as cols
import matplotlib.animation as animation
import math

import SPHAdvance as hadv
import SPHSetup as hset

def advance(n, r, v, gtime, tstep, t, lim, counter, press, dens, temps,
            xs, ys, vxs, vys, dfig, tfig, dax, tax, dframes, tframes):
    # Note that the only reason xs, ys, vxs, vys are used as parameters
    # is so that they don't need to be fully recreated on each execution of this
    # subroutine. Their actual contents upon being passed into the subroutine
    # are not important, as long as they are of the correct shape
    #press, dens, temps, r, v, counter = hadv.hydrosim(r, v, temps, gtime, tstep, counter, n, 3)
    press, dens, temps, r, v, counter = hadv.hydrosim(r, v, temps, gtime, tstep, counter)

    t = t + gtime

    # Placing r, v components into individual arrays
    for i in range(n):
        xs[i] = r[0,i]
        ys[i] = r[1,i]
        vxs[i] = v[0,i]
        vys[i] = v[1,i]

    # Finding minimum and maximum values
    # Note that only dmin, dmax, tmin and tmax are strictly required
    # xmin, xmax, ymin and ymax are included for the sake of completeness
    # There are no variables zmin or zmax as the z-data is not plotted
    xmin = np.min(xs)
    xmax = np.max(xs)
    ymin = np.min(ys)
    ymax = np.max(ys)
    dmin = np.min(dens)
    dmax = np.max(dens)
    tmin = np.min(temps)
    tmax = np.max(temps)

    # Generating animation frames
    dnorm = cols.Normalize(dmin, dmax)
    tnorm = cols.Normalize(tmin, tmax)
    cmap = cm.get_cmap()
    dframe = dax.scatter(xs, ys, c=dens, norm=dnorm, cmap=cmap.resampled(256))
    tframe = tax.scatter(xs, ys, c=temps, norm=tnorm, cmap=cmap.resampled(256))
    dframes.append([dframe])
    tframes.append([tframe])
    return r, v, t, counter, press, dens, temps, dframes, tframes
        

# This goes before arrays are initialised because array sizes depend on n
print("Enter the number of particles to be simulated:")
n = int(input())

# Initialising arrays, variables and figures
t = 0
counter = 1
dframes = []
tframes = []
press = np.zeros(shape=n)
dens = np.zeros(shape=n)
xs = np.zeros(shape=n)
ys = np.zeros(shape=n)
vxs = np.zeros(shape=n)
vys = np.zeros(shape=n)
vinit = np.zeros(shape=n)

print("Enter the total amount of time that should be simulated (in seconds):")
lim = float(input())
print("Enter the amount of simulated time between successive graphs (in seconds):")
gtime = float(input())

ngraphs = int(lim/gtime)
print(ngraphs)
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
    if boxform==1:
        r, v = hset.sprsetup(n, 3, 0, rad)
    else:
        r, v = hset.sprsetup(n, 3, 1, rad)

elif boxform==3:
    print("Enter the size of the cuboid along the x-axis (in metres):")
    xlim = float(input())
    print("Enter the size of the cuboid along the y-axis (in metres):")
    ylim = float(input())
    print("Enter the size of the cuboid along the z-axis (in metres):")
    zlim = float(input())
    r = hset.cubsetup(n, lim=[xlim, ylim, zlim], ndims=3)
    v = np.zeros(shape=np.shape(r))

print("Initial Conditions Generated.")

dfig, dax = plt.subplots()
tfig, tax = plt.subplots()

# Creating figures for the initial state - these will be the first frames of the animations
for i in range(n):
    xs[i] = r[0,i]
    ys[i] = r[1,i]
    vxs[i] = v[0,i]
    vys[i] = v[1,i]
    vinit[i] = math.sqrt((v[0,i]**2) + (v[1,i]**2) + (v[2,i]**2))


xmin = np.min(xs)
xmax = np.max(xs)
ymin = np.min(ys)
ymax = np.max(ys)
dmin = np.min(dens)
dmax = np.max(dens)
tmin = np.min(temps)
tmax = np.max(temps)
vinitmin = np.min(vinit)
vinitmax = np.max(vinit)

cmap = cm.get_cmap()
vnorm = cols.Normalize(vinitmin, vinitmax)
plt.scatter(xs, ys, c=vinit, norm=vnorm, cmap=cmap.resampled(256))
plt.colorbar()
plt.show()

dnorm = cols.Normalize(dmin, dmax)
tnorm = cols.Normalize(tmin, tmax)
dframe = dax.scatter(xs, ys, c=dens, norm=dnorm, cmap=cmap.resampled(256))
tframe = tax.scatter(xs, ys, c=temps, norm=tnorm, cmap=cmap.resampled(256))
dframes.append([dframe])
tframes.append([tframe])
# End of Initialisation
print("Beginning Simulation")
for i in range(ngraphs):
    r, v, t, counter, press, dens, temps, dframes, tframes = advance(n, r, v, gtime, tstep, t, lim, counter, press, dens, temps,
                                                                     xs, ys, vxs, vys, dfig, tfig, dax, tax, dframes, tframes)
    print("Simulation Continuing.", ((t/lim)*100), "% complete")

print("Simulation Completed")
anidens = animation.ArtistAnimation(fig=dfig, artists=dframes, interval=gtime*1000, repeat=False)
anitemp = animation.ArtistAnimation(fig=tfig, artists=tframes, interval=gtime*1000, repeat=False)
anidens.save(filename="DensityMap.gif", writer="pillow")
anitemp.save(filename="TemperatureMap.gif", writer="pillow")
plt.show()



        


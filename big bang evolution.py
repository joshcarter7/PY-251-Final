The text below is python code from my final project for an undergraduate class. The assignment I chose was to model the evolution of the early universe by starting with particles in 
random positions as a rough representation of the distribution of matter soon after the Big Bang. The particles were then allowed to evolve due to the force of gravity exerted by the 
other particles with varying time steps and total amounts of time elapsed. The end result showed “filament” structures that resemble what we expect structure of the early universe to 
look like.



###################################################################################################



# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 19:01:55 2019

"""

import numpy as np
import pylab as py
import time
import matplotlib
import matplotlib.animation as animation
import matplotlib.image as mgimg
matplotlib.use('Agg')

py.close('all')

t = time.time()

G = 1.0
dt = .01
n = 80
timeSteps = 200
def Force(particles):
    numParticles = np.size(particles,0) # Finds the number of rows
    F3D = np.zeros([numParticles, 2, numParticles])
    for p1 in range(numParticles - 1): # Don't need to calculate anything for last particle
        for p2 in range(p1+1,numParticles):

            # Particle 1
            m1 = particles[p1,0]
            x1 = particles[p1,1]
            y1 = particles[p1,2]
            # Particle 2
            m2 = particles[p2,0]
            x2 = particles[p2,1]
            y2 = particles[p2,2]
            
            Fx, Fy = Fg(m1, m2, x1, y1, x2, y2)
            
            F3D[p2,0,p1] = Fx
            F3D[p2,1,p1] = Fy
            # Newtons Third law
            F3D[p1,0,p2] = -Fx
            F3D[p1,1,p2] = -Fy
    # Sum the forces
    F = np.zeros([numParticles,2])
    for i in range(numParticles):
        F[i,0] = np.sum(F3D[:,0,i])
        F[i,1] = np.sum(F3D[:,1,i])
    return F
        
def Fg(m1, m2, x1, y1, x2, y2):
    # Calculates the Force of gravity vector betweeen 2 objects
    k = G*m1*m2/((x2 - x1)**2 + (y2 - y1)**2)**1.5
    Fx = k*(x2-x1)
    Fy = k*(y2-y1)
    return Fx, Fy 

#Updates particles
def update_particles(particles, forces):
    for i in range(0, (n**2)):
        m = particles[i,0]
        particles[i,1] += particles[i,3]*dt + forces[i,0]*.5*dt**2/m
        particles[i,2] += particles[i,4]*dt + forces[i,1]*.5*dt**2/m
        x = particles[i,1]
        y = particles[i,2]
        particles[i,3] += forces[i,0]*dt/m
        particles[i,4] += forces[i,1]*dt/m
        # Periodic boundary conditions
        
        if x > n:
            x = x%n
        if x < 0:
            x = n - np.abs(x)%n
        if y > n:
            y = y%n
        if y < 0:
            y = n - np.abs(y)%n

        particles[i,1] = x
        particles[i,2] = y
    return particles

# Update mesh
def update_mesh(particles):
    n = np.sqrt(len(particles))
    n = n.astype(int)
    grid = np.zeros([n,n])
    for p in range(len(particles)):
        x = np.floor(particles[p,1])
        y = np.floor(particles[p,2])
        x = x.astype(int)
        y = y.astype(int)
        grid[y,x] += particles[p,0]
    return grid



#sets up grid with random data

x = np.arange(1,n+1)
y = np.arange(1,n+1)
X,Y = np.meshgrid(x,y)

grid = np.random.rand(n,n)




# Intitialize our particles
particles = np.zeros([n**2,5])
ip = 0
# determine the intitial conditions 
for i in range(n):
    for j in range(n):
        # mass
        particles[ip,0] = grid[j,i]
        # x
        particles[ip,1] = i + 0.5
        # y
        particles[ip,2] = j + 0.5
        # update index
        ip += 1




maxD = 3
fig = py.pcolor(X,Y,grid,cmap = 'inferno')

py.colorbar(fig)
py.clim(0,maxD)
filename = "frame0.png"
# save the image
py.savefig(filename)
py.close('all')
myimages = []
for trial in range(timeSteps):
    forces = Force(particles)
    particles = update_particles(particles, forces)
    filename = "frame" + str(trial+1) + ".png" 
    #grid = update_meshrid(particles)
    grid = update_mesh(particles)
    py.figure(trial + 2)
    fig = py.pcolor(X,Y,grid,cmap = 'inferno')
    py.colorbar(fig)
    py.clim(0,maxD)
    # save the image
    py.savefig(filename)
    py.close('all')

    # save the image name for the animation
    img = mgimg.imread(filename)
    imgplot = py.imshow(img)
    # append AxesImage object to the list
    myimages.append([imgplot])
    
fig1 = py.figure()
## create an instance of animation
my_anim = animation.ArtistAnimation(fig1, myimages, interval=1000, blit=True, repeat_delay=1000)
#my_anim.save("animation.mp4")



elapsed = time.time() - t
print("It took "+str(elapsed) +" seconds to run")
        
        
        






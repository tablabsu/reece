#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 14:54:44 2021

@author: reecekeller
"""

# -*- coding: utf-8 -*-

#save_results_to = '/Users/reece/.spyder-py3/rotationPlots/'
save_results_to = '/Users/reecekeller/Documents/tabLab/rotationPlots/'

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches as mpatches
import random
import time

start_time = time.time()
jet= plt.get_cmap('jet')
colors = iter(jet(np.linspace(0,1,10)))

def arrow(self, x, y, dx, dy, **kwargs):
    kwargs.setdefault('arrowstyle', 'simple, head_width=10, head_length=10')
    kwargs.setdefault('fc', 'black')
    x = self.convert_xunits(x)
    y = self.convert_yunits(y)
    dx = self.convert_xunits(dx)
    dy = self.convert_yunits(dy)
    posA = x, y
    posB = x+dx, y+dy
    a = mpatches.FancyArrowPatch(posA=posA, posB=posB, **kwargs)
    self.add_artist(a)
    return a

def brownianRK4(x, y, vx, vy, theta, gamma, m, Fact, L):     
    Fx = np.random.normal(0, np.sqrt(2*gamma*T))
    Fy = np.random.normal(0, np.sqrt(2*gamma*T))
    
    k1x=vx*dt
    k2x=(vx+k1x/2)*dt
    k3x=(vx+k2x/2)*dt
    k4x=(vx+k3x)*dt
    
    k1vx=-(gamma/m)*vx*dt
    k2vx=-(gamma/m)*(vx+k1vx/2)*dt
    k3vx=-(gamma/m)*(vx+k2vx/2)*dt
    k4vx=-(gamma/m)*(vx+k3vx)*dt
    
    k1y=vy*dt
    k2y=(vy+k1y/2)*dt
    k3y=(vy+k2y/2)*dt
    k4y=(vy+k3y)*dt
    
    k1vy=-(gamma/m)*vy*dt
    k2vy=-(gamma/m)*(vy+k1vy/2)*dt
    k3vy=-(gamma/m)*(vy+k2vy/2)*dt
    k4vy=-(gamma/m)*(vy+k3vy)*dt
    
    x = (x+(1/6)*(k1x+2*k2x+2*k3x+k4x))%L
    y = (y+(1/6)*(k1y+2*k2y+2*k3y+k4y))%L
    
    vx = vx+(1/6)*(k1vx+2*k2vx+2*k3vx+k4vx)+(dt/m)*Fx+Fact*np.cos(theta)
    vy = vy+(1/6)*(k1vy+2*k2vy+2*k3vy+k4vy)+(dt/m)*Fy+Fact*np.sin(theta)
    return x, y, vx, vy

def sector(theta, psi, phi):
    #phi half cone width
    start = psi-phi
    end = psi+phi
    if (theta >= start and theta <= end):
        #return start, end, theta;
        return True
    else:
        #return 0, 0, -1
        return False

def sector2(theta, psi, phi):
    check1 = np.abs(theta-psi) <= phi;
    check2 = np.abs(theta+2*np.pi-psi) <= phi;
    check3 = np.abs(theta-2*np.pi-psi) <= phi;
    if (check1 or check2 or check3):
        return True
    else:
        return False
    
def angle(x, y):
    
    #90, 180, and 270 degree axis conditions
    if (y > 0 and (x==0 or np.abs(x) < 10**-2)):
        return np.pi/2
    elif (y < 0 and (x==0 or np.abs(x) < 10**-2)):
        return 3*np.pi/2
    elif (y == 0 and x < 0):
        return np.pi
    
    #2nd and 3rd quadrant
    elif ((x < 0 and y > 0) or (x < 0 and y < 0)):
        return np.arctan(y/x) + np.pi
    
    #4th quadrant
    elif (x > 0 and y < 0):
        return np.arctan(y/x) % (2*np.pi)
    else:
        return np.arctan(y/x)
    
def plotCones(x, y, theta, phi, r):
    for k in range(n-1):
        fig, ax = plt.subplots(figsize=(13, 10))
        #ax = plt.gca()
        for p in range(nParticles):
            z = np.e**(1j*theta[k, p])
            zl = np.e**(1j*phi)*z
            zr = np.e**(-1j*phi)*z
            
            zy=np.imag(z)
            zx=np.real(z)
            
            zly = np.imag(zl)
            zlx = np.real(zl)
            
            zry = np.imag(zr)
            zrx = np.real(zr)
            #print(zlx**2+zly**2)
            
            if p==0:
                cl = 'r'
            elif p==1:
                cl = 'b'
            elif p==2:
                cl = 'g'
            elif p==3:
                cl = 'm'
            else:
                cl = 'c'
            plt.quiver(x[k, p], y[k, p], zx, zy, color='k', angles = 'xy')
            # plt.quiver(x[k, p], y[k, p], zlx, zly, color=cl, angles = 'xy', )
            # plt.quiver(x[k, p], y[k, p], zrx, zry, color=cl, angles = 'xy')
            #arrow(ax, x[k, p], y[k, p], r*zrx, r*zry, color = cl)
            #arrow(ax, x[k, p], y[k, p], r*zlx, r*zly, color = cl)
            plt.plot([x[k, p], x[k, p]+r*zlx], [y[k, p], y[k, p]+r*zly], color = cl, label = 'p%d = %d'%(p, hits[k, p]))
            plt.plot([x[k, p], x[k, p]+r*zrx], [y[k, p], y[k, p]+r*zry], color = cl)
            t = np.linspace(theta[k, p] - phi , theta[k, p]+phi , 100) # angular range of sector
            xr = r* np.cos(t) + x[k, p] # sector x coords
            yr = r* np.sin(t) + y[k, p] # sector y coords
            plt.plot(xr,yr, color = cl, ls = '-', lw = 2) # plot the sector
        plt.ion()
        # if k<25:
        #     for p in range(nParticles):
        #         plt.plot(x[:k+1, p], y[:k+1, p], 'k.', markersize=2)
        # else:
        #     for p in range(nParticles):
        #         plt.plot(x[k-24:k+1, p], y[k-24:k+1, p], 'k.', markersize=2)
        #plt.title('D_r=%1.3f, F_Active=%1.3f, t=%d'%(D_r, Fact, k))
        #plt.title("theta = %1.2f, r = %d, N_p = %d" %(2*phi*180/np.pi, r, nParticles))
        plt.title('Collision Counting with 5 Particles')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', fontsize = 15)
        plt.grid(alpha=.5)
        #plt.style.use('dark_background') #background   
        plt.axis("square")
        plt.xlim([0, L])
        plt.ylim([0, L])
        plt.savefig(save_results_to + 'image %d'%k)
        plt.close(fig)
        plt.pause(0.01)
        
def plotTrajectory(x, y, theta):
    start = time.time()
    for k in range(n-1):
        fig, ax = plt.subplots(figsize=(10, 8))
        #ax = plt.gca()
        for p in range(nParticles):
            z = np.e**(1j*theta[k, p])
            zy=np.imag(z)
            zx=np.real(z)
            if (p==0):
                cl = 'r'
                zl = np.e**(1j*phi)*z
                zr = np.e**(-1j*phi)*z
                
                zy=np.imag(z)
                zx=np.real(z)
                
                zly = np.imag(zl)
                zlx = np.real(zl)
                
                zry = np.imag(zr)
                zrx = np.real(zr)
                plt.plot([x[k, p], x[k, p]+r*zlx], [y[k, p], y[k, p]+r*zly], color = cl, label = 'p0=%d'%hits[k, 0])
                plt.plot([x[k, p], x[k, p]+r*zrx], [y[k, p], y[k, p]+r*zry], color = cl)
                t = np.linspace(theta[k, p] - phi , theta[k, p]+phi , 100) # angular range of sector
                xr = r* np.cos(t) + x[k, p] # sector x coords
                yr = r* np.sin(t) + y[k, p] # sector y coords
                plt.plot(xr,yr, color = cl, ls = '-', lw = 2) # plot the sector
                plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', fontsize = 15)
            else:
                cl = 'k'
            plt.quiver(x[k, p], y[k, p], zx, zy, color=cl, angles = 'xy')   
        plt.ion()
        # if k<25:
        #     for p in range(nParticles):
        #         plt.plot(x[:k+1, p], y[:k+1, p], 'k.', markersize=2)
        # else:
        #     for p in range(nParticles):
        #         plt.plot(x[k-24:k+1, p], y[k-24:k+1, p], 'k.', markersize=2)
        #plt.title('D_r=%1.3f, F_Active=%1.3f, t=%d'%(D_r, Fact, k))
        plt.title("theta = %1.2f, r = %d, N_p = %d" %(2*phi*180/np.pi, r, nParticles))
        plt.grid(alpha=.5)
        #plt.style.use('dark_background') #background   
        plt.axis("square")
        plt.xlim([0, L])
        plt.ylim([0, L])
        plt.savefig(save_results_to + 'image %d'%k)
        plt.close(fig)
        plt.pause(0.01)
    end = time.time()
    t = end - start
    return t


def intCircle(n, dt, gamma, m, T, D_r, Fact, nParticles, L, r):
    x=np.zeros([n, nParticles]); y=np.zeros([n, nParticles]); 
    vx=np.zeros([n, nParticles]); vy=np.zeros([n, nParticles]); 
    theta=np.zeros([n, nParticles]); 
    startX=np.random.randint(L/2, size=(1, nParticles))
    startY=np.random.randint(L/2, size=(1, nParticles))
    x[0, :]=startX
    y[0, :]=startY
    
    startTheta=random.sample(range(0, 360), nParticles)
    theta[0, :] = np.asarray(startTheta)*np.pi/180
        
    hits = np.zeros([n, nParticles])
    
    v_a = np.zeros([n-1, 1])
    sigma = np.sqrt(2*D_r*dt)
    for i in range(n-1):
        for p in range(nParticles):
            dTheta=np.random.normal(0, sigma)
            neighbors=[]
            count  = 0
            for k in range(nParticles):
                inCircle = (x[i, k]-x[i, p])**2+(y[i, k]-y[i, p])**2 <= r**2
                if (inCircle):
                    neighbors.append(k)      
                    count +=1
            if (len(neighbors) == 0):
                thetaNeighbors=theta[i, p]
            else:
                thetaNeighbors = [theta[i, :][j] for j in neighbors]
            thetaMean_r = np.mean(thetaNeighbors)
            theta[i+1, p]=thetaMean_r + dTheta
            hits[i, p] = count
            x[i+1, p], y[i+1, p], vx[i+1, p], vy[i+1, p] = brownianRK4(x[i, p], y[i, p], vx[i, p], vy[i, p], theta[i, p], gamma, m, Fact, L)
            
            vxSum = sum(vx[i+1, :])
            vySum = sum(vy[i+1, :])
            vAvg = sum(np.sqrt(vx[i+1, :]**2 + vy[i+1, :]**2))
            v_a[i] = np.sqrt(vxSum**2+vySum**2)/vAvg
            
    return x, y, theta, v_a, hits

def intCone(n, dt, gamma, m, T, D_r, Fact, nParticles, L, phi, r):
    start = time.time();
    x=np.zeros([n, nParticles]); y=np.zeros([n, nParticles]); 
    vx=np.zeros([n, nParticles]); vy=np.zeros([n, nParticles]); 
    theta=np.zeros([n, nParticles]); 
    startX=np.random.randint(L, size=(1, nParticles))
    startY=np.random.randint(L, size=(1, nParticles))
    x[0, :]=startX
    y[0, :]=startY
    
    startTheta=random.sample(range(0, 360), nParticles)
    theta[0, :] = np.asarray(startTheta)*np.pi/180
    #theta[0, :] = [np.pi/2 for i in range(nParticles)]
    
    hits = np.zeros([n, nParticles])
    
    v_a = np.zeros([n-1, 1])
    sigma = np.sqrt(2*D_r*dt)
    for i in range(n-1):
        for p in range(nParticles):
            dTheta=np.random.normal(0, sigma)
            neighbors=[]
            psi = theta[i, p]
            counter = 0
            for k in range(nParticles):
                # p'th particle is the center 
                #checking all k particles
                inCircle = (x[i, k]-x[i, p])**2+(y[i, k]-y[i, p])**2 < r**2
                if (inCircle):
                    if (x[i, k]!=x[i, p] or y[i, k]!=y[i, p]):
                        vec = [x[i, k]-x[i, p], y[i, k]-y[i, p]]
                        delta = angle(vec[0], vec[1])
                        if (sector2(delta, psi, phi)):
                            neighbors.append(k)
                            counter+=1
                    # elif (x[i, k] == x[i, p] and y[i, k]==y[i, p]):
                    #     neighbors.append(k)
                    #     counter+=1
            if (len(neighbors) == 0):
                thetaNeighbors = theta[i, p]
            else:
                thetaNeighbors = [theta[i, :][j] for j in neighbors]
            thetaMean_r = np.mean(thetaNeighbors)
            theta[i+1, p]= (thetaMean_r + dTheta) % (2*np.pi)
            hits[i, p] = counter                     
 
            x[i+1, p], y[i+1, p], vx[i+1, p], vy[i+1, p] = brownianRK4(x[i, p], y[i, p], vx[i, p], vy[i, p], theta[i, p], gamma, m, Fact, L)
            
            vxSum = sum(vx[i+1, :])
            vySum = sum(vy[i+1, :])
            vAvg = sum(np.sqrt(vx[i+1, :]**2 + vy[i+1, :]**2))
            v_a[i] = np.sqrt(vxSum**2+vySum**2)/vAvg
            end = time.time()
            t = end - start
    return x, y, theta, v_a, hits, t


def error(k, steadyState, nParticles, L):
    data = []
    for i in range(k):
        _, _, _, v_a = intCircle(n, dt, gamma, m, T, D_r, Fact, nParticles, L)
        data.append(v_a[steadyState])
    sx = np.std(data)
    return sx

def avg(hits, nParticles):
    sums=[]
    for p in range(nParticles):
        sums.append(sum(hits[:, p]))
    return np.mean(sums)

def singular(hits, phi):
    plt.plot(hits, label = 'theta = %1.2f' % (2*phi*180/np.pi))
    plt.title('Collisions for p0')
    plt.xlabel('t')
    plt.ylabel('collisions')
    plt.legend()
    plt.grid()
    
def normalized(hits, phi):
    plt.plot(hits/(2*phi), label = 'theta = %1.2f' % (2*phi*180/np.pi))
    plt.title('Normalized Plot')
    plt.xlabel('t')
    plt.ylabel('collisions per degree')
    plt.legend()
    plt.grid()
    
def v2arr(v_a):
    new = []
    for i in range(np.shape(v_a)[0]):
        new.append(v_a[i][0])
    return new
dt=0.1; D_r=0.3; Fact=1
n = 501; gamma=1; m=1; T=1
nParticles = 40; L=20
phi = np.pi/3
r=1

c=(L/20)*dt*np.sqrt(n)*Fact

#phis = np.asarray([np.pi/12, np.pi/4, np.pi/3, np.pi/2]);
phis = np.linspace(0, np.pi, 20)
Hits = []
nP = [20, 40, 80]
L = [2.2*np.sqrt(4), 3.1*np.sqrt(4), 4.47*np.sqrt(4)]
master=[]
for i in range(3):
    order=[]
    for phi in phis:
        x, y, theta, v_a, hits, tRun = intCone(n, dt, gamma, m, T, D_r, Fact, nP[i], L[i], phi, r)
        #v = v2arr(v_a)
        v = np.mean(v_a[60:])
        order.append(v)
        #Hits.append(hits[:, 0])
    master.append(order)
#Hits = np.asarray(Hits)

    
plt.plot(phis, master[0], '-or', phis, master[1], '-ob', phis, master[2], '-og')
plt.xlabel('View Angle')
plt.ylabel('Order Param')


#tPlot = plotTrajectory(x, y, theta)

# angles = np.asarray([15, 30, 60, 90, 150, 230, 360])*np.pi/180

# for phi in angles: 
#     _, _, _, v_a, _ = intCone(n, dt, gamma, m, T, D_r, Fact, nParticles, L, phi, r)
#     plt.plot(v_a, label = 'phi = %d'%(phi*180/np.pi))
#     plt.legend()


# D_r = np.linspace(0, 60, 60)

# nP = [20, 40, 80, 150]
# angles = np.asarray([15, 30, 60, 90, 150, 230, 360])*np.pi/180
# # L = [2.2, 3.1, 4.47, 5.4];
# ls2=[]
# barsTot = []
# va=[]
# for i in range(len(nP)):
#     ls=[]
#     bars = []
#     v=[]
#     for j in range(len(D_r)):
#         _, _, _, v_a, _ = intCone(n, dt, gamma, m, T, D_r[j], Fact, nP[i], L, phi, r)
#         data = v_a[60:]
#         ls.append(np.mean(data))
#         bars.append(np.std(data))
#         v.append(v_a)
#     barsTot.append(bars)
#     ls2.append(ls)
#     va.append(v)

# plt.fill_between(D_r, np.asarray(ls2[0])-np.asarray(barsTot[0]), np.asarray(ls2[0])+np.asarray(barsTot[0]), color = 'b', alpha = 0.3)
# plt.fill_between(D_r, np.asarray(ls2[1])-np.asarray(barsTot[1]), np.asarray(ls2[1])+np.asarray(barsTot[1]), color = 'g', alpha = 0.3)
# plt.fill_between(D_r, np.asarray(ls2[2])-np.asarray(barsTot[2]), np.asarray(ls2[2])+np.asarray(barsTot[2]), color = 'm', alpha = 0.3)
# plt.fill_between(D_r, np.asarray(ls2[3])-np.asarray(barsTot[3]), np.asarray(ls2[3])+np.asarray(barsTot[3]), color = 'r', alpha = 0.3)

# plt.title('Vision Cone Analysis')
# plt.ylabel('v_a')
# plt.xlabel('D_r')

# plt.plot(D_r, ls2[0], 'b-x', label = 'N=20')    
# plt.plot(D_r, ls2[1], 'g-x', label = 'N=40')
# plt.plot(D_r, ls2[2], 'm-x', label = 'N=80')
# plt.plot(D_r, ls2[3], 'r-x', label = 'N=150')
# plt.legend()


# # density = [i/20 for i in nP]
# # plt.plot(density, ls2[0], 'm-x')
# # plt.fill_between(density, np.asarray(ls2[0])-np.asarray(barsTot[0]), np.asarray(ls2[0])+np.asarray(barsTot[0]), color = 'm', alpha = 0.3)
# # plt.title('L = 20, D_r = 2')
# # plt.xlabel('Density')
# # plt.ylabel('v_a')
# print("--- %s seconds ---" % (time.time() - start_time))

###########################

# steadyState = 60
# bars1=[]; bars2=[]; bars3=[]; v_ss1 = []; v_ss2=[]; v_aa=[]; v_ss3=[];
# xx=[]; yy=[]
# fig = plt.figure(figsize=(10, 10))
# ax1 = fig.add_subplot(2, 1, 1)
# k=20;
# for p in nP:
#     sx1 = error(k, steadyState, p, L)
#     sx2 = error(k, steadyState*2, p, L)
#     sx3 = error(k, steadyState*3, p, L)
#     bars1.append(sx1); bars2.append(sx2); bars3.append(sx3)
#     x, y, theta, v_a = accum(n, dt, gamma, m, T, sigma, Fact, p, L)
#     v_ss1.append(v_a[steadyState])
#     v_ss2.append(v_a[2*steadyState])
#     v_ss3.append(v_a[3*steadyState])
#     rho = p/boxSize
#     print('\n Particles=%d\n'%p, 'Start theta', theta[0, :])
#     print('\n', 'Start x', x[0, :])
#     print('\n', 'Start y', y[0, :])
#     plt.plot(v_a, label='N=%d'%p)
# plt.legend()
# plt.grid()
# plt.xlabel('t'); plt.ylabel('v_a')
# plt.title('# of avgs: %d, t_s = %d' % (k, steadyState))
# ax2 = fig.add_subplot(2, 1, 2)
# plt.plot(nP, v_ss1, 'ro', nP, v_ss2, 'ko', nP, v_ss3, 'go')
# plt.errorbar(nP, v_ss1, yerr=bars1, ls='none', ecolor='r')
# plt.errorbar(nP, v_ss2, yerr=bars2, ls='none', ecolor='k')
# plt.errorbar(nP, v_ss3, yerr=bars3, ls='none', ecolor='g')
# plt.xlabel('N'); plt.ylabel('v_a')
# plt.legend()
# plt.grid()

 ################################

# for k
#  v_mid = [np.cos(theta[i, p]), np.sin(theta[i, p])]
#                 vL = [np.cos(phi)*v_mid[0]-np.sin(phi)*v_mid[1], np.sin(phi)*v_mid[0]+np.cos(phi)*v_mid[1]]
#                 vR = [np.cos(2*np.pi-phi)*v_mid[0]-np.sin(2*np.pi-phi)*v_mid[1], np.sin(2*np.pi-phi)*v_mid[0]+np.cos(2*np.pi-phi)*v_mid[1]]
#                 vL = r*np.asarray(vL)/np.linalg.norm(vL)
#                 vR = r*np.asarray(vR)/np.linalg.norm(vR)
#                 n_vR = np.asarray([-vR[1], vR[0]])
#                 n_vL = np.asarray([-vL[1], vL[0]])
#                 for k in range(nParticles):
#                     if ((x[i, k]-x[i, p])**2+(y[i, k]-y[i, p])**2<=r**2):
#                         #v = [x[i, k], y[i, k]]
#                         v = [np.cos(theta[i, k]), np.sin(theta[i, k])]
#                         v = r*np.asarray(v)/np.linalg.norm(v)
#                         left = np.dot(v, n_vL)
#                         right = np.dot(v, n_vR)
#                         if (left <= 0 and right >= 0):
#                             neighbors.append(k)       


    
    
    
    
    
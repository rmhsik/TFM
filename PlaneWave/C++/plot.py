import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import copy

dataDir = 'Data/PhiG/'
imageDir = 'Images/PhiG/'
mcmap = copy.copy(matplotlib.cm.get_cmap('gray'))
mcmap.set_bad((0,0,0))

for filename in os.listdir(dataDir):
    data = np.genfromtxt(dataDir+'/'+filename,dtype=np.double,delimiter=' ')
    matplotlib.use('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #plot = ax.imshow(data,norm=matplotlib.colors.LogNorm(vmin = 1E-6, vmax = 1.0),origin='lower',cmap=mcmap,aspect=0.5)
    plot = ax.imshow(data,vmin = 1E-6, vmax = 1.0,origin='lower',cmap=mcmap,aspect=0.5)
    fig.colorbar(plot)
    plt.savefig(imageDir+'/'+filename[:-4]+'.png')
    plt.close()


dataDir = 'Data/PhiE/'
imageDir = 'Images/PhiE/'

for filename in os.listdir(dataDir):
    data = np.genfromtxt(dataDir+'/'+filename,dtype=np.double,delimiter=' ')
    matplotlib.use('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #plot = ax.imshow(data,norm=matplotlib.colors.LogNorm(vmin = 1E-6, vmax = 1.0),origin='lower',cmap=mcmap,aspect=0.5)
    plot = ax.imshow(data,vmin = 1E-6, vmax = 1.0,origin='lower',cmap=mcmap,aspect=0.5)
    fig.colorbar(plot)
    plt.savefig(imageDir+'/'+filename[:-4]+'.png')
    plt.close()

dataDir = 'Data/PhiT/'
imageDir = 'Images/PhiT/'

for filename in os.listdir(dataDir):
    data = np.genfromtxt(dataDir+'/'+filename,dtype=np.double,delimiter=' ')
    matplotlib.use('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #plot = ax.imshow(data,norm=matplotlib.colors.LogNorm(vmin = 1E-6, vmax = 1.0),origin='lower',cmap=mcmap,aspect=0.5)
    plot = ax.imshow(data,vmin = 1E-6, vmax = 1.0,origin='lower',cmap=mcmap,aspect=0.5)
    fig.colorbar(plot)
    plt.savefig(imageDir+'/'+filename[:-4]+'.png')
    plt.close()

dataDir = 'Data/Potential.dat'
imageDir = 'Images/'

data = np.genfromtxt(dataDir,dtype=np.double,delimiter=' ')
matplotlib.use('Agg')
fig = plt.figure()
ax = fig.add_subplot(111)
plot = ax.imshow(data,origin='lower',cmap='RdGy',aspect=0.5)
fig.colorbar(plot)
plt.savefig(imageDir+'/'+'Potential.png')
plt.close()
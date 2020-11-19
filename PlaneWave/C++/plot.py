import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os

dataDir = 'Data/PhiG/'
imageDir = 'Images/PhiG/'

for filename in os.listdir(dataDir):
    data = np.genfromtxt(dataDir+'/'+filename,dtype=np.double,delimiter=' ')
    matplotlib.use('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot = ax.imshow(data,origin='lower',cmap='RdGy',aspect=2.0,vmin=0.1, vmax=1.0)
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
    plot = ax.imshow(data,origin='lower',cmap='RdGy',aspect=2.0,vmin=0.1, vmax=1.0)
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
    plot = ax.imshow(data,origin='lower',cmap='RdGy',aspect=2.0,vmin=0.1, vmax=1.0)
    fig.colorbar(plot)
    plt.savefig(imageDir+'/'+filename[:-4]+'.png')
    plt.close()
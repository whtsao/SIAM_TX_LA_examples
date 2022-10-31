import numpy as np
import csv
import os
import matplotlib.pyplot as plt

file_vof = 'column_gauges.csv'

def readProbeFile(filename):
    with open (filename, 'r') as csvfile:
        #read data
        data=np.loadtxt(csvfile, delimiter=",", skiprows=1)
        time=data[:,0]
        data = data[:,1:]
        #parse header
        csvfile.seek(0)
        header = csvfile.readline()
        header = header.replace("time","")
        header = header.replace("[","")
        header = header.replace("]","")
        header = header.replace(","," ")
        header = header.split()
        probeType = []
        probex = []
        probey = []
        probez = []        
        for ii in range(0,len(header),4):
            probeType.append(header[ii])
            probex.append(float(header[ii+1]))
            probey.append(float(header[ii+2]))
            probez.append(float(header[ii+3]))
        probeCoord = np.column_stack([np.array(probex),np.array(probey),np.array(probez)])
        datalist = [probeType,probeCoord,time,data]
        return datalist

data_vof = readProbeFile(file_vof)

#####################################################################################

# Extracting probes
time = data_vof[2]
vof = data_vof[3]

ETA = []
tank_dim = [37.7, 0.75]
waterLevel = 0.4

# gauge_x = [6.0,
#            11.0,
#            12.0,
#            13.0,
#            14.0,
#            15.0,
#            16.0,
#            17.0]

gauge_x = [6.0,
           10.8,
           12.8,
           13.8,
           14.8,
           16.0,
           17.6]

for j in range(data_vof[1].shape[0]//2):
    #eta = []
    #for i in range(vof.shape[0]-1):
    #    eta.append(tank_dim[1]-vof[1:,j][i]-waterLevel)
    ETA.append(tank_dim[1]-vof[1:,j]-waterLevel)
#ETA = np.array(ETA)
ETA = np.vstack(ETA)

#####################################################################################

# Plotting the probes
fig = plt.figure(figsize=(25,15))
ax = ['' for x in range(data_vof[1].shape[0]//2)]
for i in range(data_vof[1].shape[0]//2):
    ax[i] = fig.add_subplot(4,3,i+1)
    ax[i].plot(time[1:], ETA[i], 'r')
    ax[i].set_ylim([-0.02,0.02])
    ax[i].tick_params(labelsize=10)
    ax[i].set_title('Eta [m] against time [sec] at x='+str(gauge_x[i]), color='b', fontsize=12)
    ax[i].grid()
plt.tight_layout() 
plt.savefig('eta.png')   
#plt.show()

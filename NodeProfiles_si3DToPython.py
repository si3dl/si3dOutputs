# -------------------------------------------------------------------------
# This script uses the text files generated from the si3D simulations on
# single nodes and saves the results in a python format as a dictionary. The code is based of previous versions created by Alicia Cortes and others authors related to the use of this numerical model.

# Version 1:
# - This version transform the text file result from SI3D into a python dictionary for the file name specified.
# - This script also creates a datetime variable for plotting purposes and time of simulation.
# - This script also considers the condition when tracers are used in the simulations. The code will work by utlizing the ntracer variable.

# USER GUIDE:
# Output = NodeProfles_si3DToPython(file,TimeEndSim,ipt,ntracer,startDate):
# VariableName = Format. Description of variable.
# file = String. Name of node profile from si3D. The file name follows 'tfim_jm.txt'
# TimeEndSim = Float. Seconds of end of simulation. It is obtained from si3d input file 'si3d_inp.txt' as tl
# ipt = Float. Number of time steps between consecutive outputs of resulting vertical node profiles from si3D simulations. It is obtained from si3d input file 'si3d_inp.txt' as ipt
# ntracer = Float. Number of tracers used in the numerical simulations. It is obtained from si3d input file 'si3d_inp.txt' as ntr
# StartDate = Datetime. Date of start of numerical simulations.

# NOTES:
# 1. The resulting structure will have on the columns the time steps and on the rows the depths.

# Copy Right Sergio A. Valbuena 2021
# Author: Sergio Valbuena
# Date: 05-12-2020

# Library Import
import datetime as Dt
import numpy as np
import os

# Definition of variable
def is_eof(f):
    cur = f.tell()    # save current position
    f.seek(0, os.SEEK_END)
    end = f.tell()    # find the size of file
    f.seek(cur, os.SEEK_SET)
    return cur == end

def NodeProfles_si3DToPython(file,TimeEndSim,dt,ipt,ntracer,startDate):
    TimeEndSimhrs = TimeEndSim/60/60
    fid = open(file,'r+')
    headerline = fid.readline()
    runNumber = fid.readline()
    line0 = fid.readline()
    im = int(float(line0[4:7]))
    jm = int(float(line0[13:16]))
    km = int(float(line0[23:26]))
    # dt = int(float(line0[33:38]))
    dtout = dt*ipt/3600
    # hs = float(line0[49:55])
    # ddz = float(line0[64:69])
    line1 = fid.readline()
    line2 = fid.readline()
    line3 = fid.readline()
    line4 = fid.readline()
    ntime = 0

    del line1, line2, line3, line4, headerline, runNumber, im, jm

    timesim = np.empty(int(TimeEndSim/dt/ipt)+1)
    h = np.empty(int(TimeEndSim/dt/ipt)+1)
    z = np.empty((km,int(TimeEndSim/dt/ipt)+1))
    u = np.empty((km,int(TimeEndSim/dt/ipt)+1))
    v = np.empty((km,int(TimeEndSim/dt/ipt)+1))
    w = np.empty((km,int(TimeEndSim/dt/ipt)+1))
    Av = np.empty((km,int(TimeEndSim/dt/ipt)+1))
    Dv = np.empty((km,int(TimeEndSim/dt/ipt)+1))
    s = np.empty((km,int(TimeEndSim/dt/ipt)+1))

    if ntracer == 0:
        # Read first line for surface node at horizontal location (i,j)
        data1 = np.fromfile(fid, count = 10, sep=' ', dtype=np.float32)
        timesim[ntime] = data1[0]
        h[ntime] = data1[2]
        z[0,ntime] = data1[3]
        u[0,ntime] = data1[4]
        v[0,ntime] = data1[5]
        w[0,ntime] = data1[6]
        Av[0,ntime] = data1[7]
        Dv[0,ntime] = data1[8]
        s[0,ntime] = data1[9]
        # Read reast of rows for initial time t =0
        data2 = np.fromfile(fid,count = 7*(km-1), sep = ' ', dtype=np.float32)
        z[1:km,ntime] = data2[0:7*(km-1):7]
        u[1:km,ntime] = data2[1:7*(km-1):7]
        v[1:km,ntime] = data2[2:7*(km-1):7]
        w[1:km,ntime] = data2[3:7*(km-1):7]
        Av[1:km,ntime] = data2[4:7*(km-1):7]
        Dv[1:km,ntime] = data2[5:7*(km-1):7]
        s[1:km,ntime] = data2[6:7*(km-1):7]

        del data1, data2
        while (timesim[ntime] +dtout) <= TimeEndSimhrs and is_eof(fid) == False:
            ntime += 1
            data1 = np.fromfile(fid, count = 10, sep=' ', dtype=np.float32)
            timesim[ntime] = data1[0]
            h[ntime] = data1[2]
            z[0,ntime] = data1[3]
            u[0,ntime] = data1[4]
            v[0,ntime] = data1[5]
            w[0,ntime] = data1[6]
            Av[0,ntime] = data1[7]
            Dv[0,ntime] = data1[8]
            s[0,ntime] = data1[9]

            data2 = np.fromfile(fid,count = 7*(km-1), sep = ' ', dtype=np.float32)
            z[1:km,ntime] = data2[0:7*(km-1):7]
            u[1:km,ntime] = data2[1:7*(km-1):7]
            v[1:km,ntime] = data2[2:7*(km-1):7]
            w[1:km,ntime] = data2[3:7*(km-1):7]
            Av[1:km,ntime] = data2[4:7*(km-1):7]
            Dv[1:km,ntime] = data2[5:7*(km-1):7]
            s[1:km,ntime] = data2[6:7*(km-1):7]

            del data1, data2
    elif ntracer != 0:
        print('Needs to be completed')

    fid.close()

    # Replacing -99 to NaN values
    z[z == -99] = np.nan
    u[u == -99] = np.nan
    v[v == -99] = np.nan
    w[w == -99] = np.nan
    Av[Av == -99] = np.nan
    Dv[Dv == -99] = np.nan
    s[s == -99] = np.nan
    h[h == -99] = np.nan

    # Creating the dictionary for the node.
    dummy = {}
    dummy['TimeSimHrs'] = timesim
    Deltasec = TimeEndSim
    EndSim = int(Deltasec + dt*ipt)
    step = int(dt*ipt)
    dummy['TimeDateLocal'] = np.array([startDate + Dt.timedelta(0,t) for t in range(0,EndSim,step)])
    dummy['z'] = z
    dummy['u'] = u
    dummy['v'] = v
    dummy['w'] = w
    dummy['Av'] = Av
    dummy['Dv'] = Dv
    dummy['T'] = s
    dummy['h'] = h
    dummy['comments'] = [['Sim_Time_hrs','Time in hours of the simulation run'],['TimeDateLocal','Date time variable of simulation'],['u','[cm/s] U component of velocity'],['v','[cm/s] V component of velocity'],['w','[cm/s] W component of velocity'],['Av','[cm2/s}'],['Dv','[Cm2/s]'], ['T','[C] Temperature'],['h',' [cm] Surface Level']]

    return dummy

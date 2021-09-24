# -------------------------------------------------------------------------
# This script uses the binary files generated from the SI3D simulations for horizontal planes (H-Planes) and saves the data as python dictionaries. The code is based of previous versions created by Alicia Cortes and others authors related to the use of this numerical model.

# Version 1:
# - This version transform the binary file result from SI3D into a python dictionary for a plane considered in the simulation.
# - This script also creates DateTime variable for plotting purposes.

# USER GUIDE:
# output = HPlane_Si3dToPython(Plane,dx,SimPath):
# variable = Format. Description.
# Plane = Integer. Number of the plane/s obtained from the Si3D numerical simulations. It is obtained from si3d input file 'si3d_inp.txt' as plane x
# dx = Float. Delta x value from numerical simulation. Grid resolution 'm'. 
# NOTES: 1. The resulting structure will have 3D matrices where the third dimensions is for the time steps. Columns are x and rows are y.

# Author: Sergio Valbuena
# Date: 05-12-2020

# Library Import
import os
import numpy as np
import pandas as pd

# Definition of variable
def is_eof(f):
    cur = f.tell()    # save current position
    f.seek(0, os.SEEK_END)
    end = f.tell()    # find the size of file
    f.seek(cur, os.SEEK_SET)
    return cur == end

def HPlane_Si3dToPython(Plane,dx):
    Hplane = 'plane_'+str(Plane)
    fid = open(Hplane,"r+")
    dum1 = np.fromfile(fid, count = 1, dtype = np.int32)
    n_frames = np.fromfile(fid, count = 1, dtype = np.int32)
    n_frames = n_frames[0]
    dum2 = np.fromfile(fid, count = 1, dtype = np.int32)
    dum3 = np.fromfile(fid, count = 1, dtype = np.int32)
    ipoints = np.fromfile(fid, count = 1, dtype = np.int32)
    ipoints = ipoints[0]
    dum4 = np.fromfile(fid, count = 1, dtype = np.int32)

    year = np.empty(n_frames+1)
    month = np.empty(n_frames+1)
    day = np.empty(n_frames+1)
    hour = np.empty(n_frames+1)
    u = np.empty((ipoints,n_frames+1))
    v = np.empty((ipoints,n_frames+1))
    w = np.empty((ipoints,n_frames+1))
    T = np.empty((ipoints,n_frames+1))
    Az = np.empty((ipoints,n_frames+1))
    l = np.empty((ipoints,n_frames+1))

    for count in range(0,n_frames+1):
        dum5 = np.fromfile(fid, count = 1, dtype = np.int32)
        st = is_eof(fid)
        if st == False:
            istep = np.fromfile(fid, count = 1, dtype = np.int32)
            year[count] = np.fromfile(fid, count = 1, dtype = np.int32)
            month[count] = np.fromfile(fid, count = 1, dtype = np.int32)
            day[count] = np.fromfile(fid, count = 1, dtype = np.int32)
            hour[count] = np.fromfile(fid, count = 1, dtype = np.float32)
            # ... Read all data for present time slice
            if count == 0:
                # ... GEOMETRY only in first step
                out_array = np.fromfile(fid,count = 8*ipoints, dtype = np.float32)
                x = out_array[0:len(out_array)-7:8]
                y = out_array[1:len(out_array)-6:8]
                u[:,count] = out_array[2:len(out_array)-5:8]
                v[:,count] = out_array[3:len(out_array)-4:8]
                w[:,count] = out_array[4:len(out_array)-3:8]
                T[:,count] = out_array[5:len(out_array)-2:8]
                Az[:,count] = out_array[6:len(out_array)-1:8]
                l[:,count] = out_array[7:len(out_array):8]
            else:
                # ... GEOMETRY only in first step
                out_array = np.fromfile(fid,count = 6*ipoints, dtype = np.float32)
                u[:,count] = out_array[0:len(out_array)-5:6]
                v[:,count] = out_array[1:len(out_array)-4:6]
                w[:,count] = out_array[2:len(out_array)-3:6]
                T[:,count] = out_array[3:len(out_array)-2:6]
                Az[:,count] = out_array[4:len(out_array)-1:6]
                l[:,count] = out_array[5:len(out_array):6]
            dum6 = np.fromfile(fid,count = 1, dtype = np.int32)
            del out_array
        else:
            n_frames = count - 1
            year = year[0:n_frames+1]
            month = month[0:n_frames+1]
            day = day[0:n_frames+1]
            hour = hour[0:n_frames+1]
            u = u[:,0:n_frames+1]
            v = v[:,0:n_frames+1]
            w = w[:,0:n_frames+1]
            T = T[:,0:n_frames+1]
            Az = Az[:,0:n_frames+1]
            l = l[:,0:n_frames+1]
            break
    # print('No of Frames = ',n_frames)
    years = year.astype(int)
    months = month.astype(int)
    days = day.astype(int)
    hour /= 100
    minsec = hour%1
    hours = hour.astype(int)
    minsec = hours%1
    minut = minsec.astype(int)
    sec = minsec*60
    df = pd.DataFrame({'year': years, 'month': months, 'day': days, 'hour': hours, 'minute': minut, 'second': sec})
    Time = pd.to_datetime(df[["year","month","day","hour","minute","second"]])

    x = x.astype(int)
    y = y.astype(int)

    xv = np.arange(1,max(x),1)
    yv = np.arange(1,max(y),1)
    [xg,yg] = np.meshgrid(xv,yv)
    xg = xg*dx
    yg = yg*dx
    [m1,m2] = np.shape(xg)

    ug = np.empty((m1,m2,n_frames+1))
    ug[:,:] = np.nan
    vg = np.empty((m1,m2,n_frames+1))
    vg[:,:] = np.nan
    wg = np.empty((m1,m2,n_frames+1))
    wg[:,:] = np.nan
    Tg = np.empty((m1,m2,n_frames+1))
    Tg[:,:] = np.nan
    Azg = np.empty((m1,m2,n_frames+1))
    Azg[:,:] = np.nan
    lg = np.empty((m1,m2,n_frames+1))
    lg[:,:] = np.nan

    for frame in range(0,n_frames+1):
        # ------- ISOLATE RECORDS FOR a GIVEN TIME STEP
        # ........ EW_velocity
        u_vect = u[:,frame]             # [m/s]
        # ......... NS_velocity
        v_vect = v[:,frame]             # [m/s]
        # # ......... Vertical_velocity
        w_vect = w[:,frame]             # [m/s]
        # # ......... Scalar field
        T_vect = T[:,frame]             # [C]
        # # ......... TKE field
        Az_vect = Az[:,frame]
        # # ......... TKE field
        l_vect = l[:,frame]             # [m]
        # ------- ARRANGE RECORDS IN SPATIALLY MEANINGFUL DOMAIN
        for j in range(int(np.min(y))-2,int(np.max(y)-1)):
            dum = np.where(y == j+2)
            # ncols = len(dum)
            ug[j,x[dum]-2,frame] = u_vect[dum]
            vg[j,x[dum]-2,frame] = v_vect[dum]
            wg[j,x[dum]-2,frame] = w_vect[dum]
            Tg[j,x[dum]-2,frame] = T_vect[dum]
            Azg[j,x[dum]-2,frame] = Az_vect[dum]
            lg[j,x[dum]-2,frame] = l_vect[dum]
            del dum

    del dum1, dum2, dum3, dum4, dum5, dum6

    output = {}
    output['x'] = xg
    output['y'] = yg
    output['u'] = ug
    output['v'] = vg
    output['w'] = wg
    output['T'] = Tg
    output['Az'] = Azg
    output['l'] = lg
    output['Time'] = Time
    output['n_frames'] = n_frames

    fid.close()

    return output

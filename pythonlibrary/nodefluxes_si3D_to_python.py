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
import pandas as pd
import os

# Definition of variable


def is_eof(f):
    cur = f.tell()    # save current position
    f.seek(0, os.SEEK_END)
    end = f.tell()    # find the size of file
    f.seek(cur, os.SEEK_SET)
    return cur == end


def nodefluxes_si3D_to_python(file, TimeEndSim, dt, ipt, startDate):
    TimeEndSimhrs = TimeEndSim / 60 / 60
    fid = open(file, 'r+')
    headerline = fid.readline()
    runNumber = fid.readline()
    line0 = fid.readline()
    im = int(float(line0[4:7]))
    jm = int(float(line0[13:16]))
    km = int(float(line0[23:26]))
    km += 1
    dtout = dt * ipt / 3600
    line1 = fid.readline()
    line2 = fid.readline()
    keys = fid.readline()
    keys = keys.strip().split()
    keys1 = keys[2:]
    del keys
    keys = keys1
    units = fid.readline()
    units = units.strip().split()
    units = units[3:]
    ntime = 0

    del line1, line2, headerline, runNumber, im, jm

    output = {}
    output['time'] = np.empty(int(TimeEndSim / dt / ipt) + 1) * np.nan
    output['step'] = np.empty(int(TimeEndSim / dt / ipt) + 1) * np.nan
    for k in keys:
        output[k] = np.empty((km, int(TimeEndSim / dt / ipt) + 1)) * np.nan

    data1 = np.fromfile(fid, count=37, sep=' ', dtype=np.float32)
    data2 = np.fromfile(fid, count=35 * (km - 1), sep=' ', dtype=np.float32)
    output['time'][ntime] = data1[0]
    output['step'][ntime] = data1[1]

    idx = 0
    for k in keys:
        if k == 'time' or k == 'step':
            continue
        output[k][0, ntime] = data1[idx + 2]
        output[k][1:km, ntime] = data2[idx:35 * (km - 1):35]
        idx += 1
    del data1, data2
    while (output['time'][ntime] + dtout) <= TimeEndSimhrs and is_eof(fid) is False:
        ntime += 1
        data1 = np.fromfile(fid, count=37, sep=' ', dtype=np.float32)
        data2 = np.fromfile(fid, count=35 * (km - 1), sep=' ', dtype=np.float32)
        output['time'][ntime] = data1[0]
        output['step'][ntime] = data1[1]
        idx = 0
        for k in keys:
            if k == 'time' or k == 'step':
                continue
            output[k][0, ntime] = data1[idx + 2]
            output[k][1:km, ntime] = data2[idx:35 * (km - 1):35]
            idx += 1

        del data1, data2
    fid.close()

    Deltasec = TimeEndSim
    EndSim = int(Deltasec + dt * ipt)
    step = int(dt * ipt)
    output['dateLocal'] = [startDate + Dt.timedelta(0, t) for t in range(0, EndSim, step)]
    output['dateLocal'] = pd.to_datetime(output['dateLocal'])
    output['units'] = units
    output['z'] = output['depth']
    keys = keys[1:]
    output['keys'] = keys
    del output['depth']

    return output

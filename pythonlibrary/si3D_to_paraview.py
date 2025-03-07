# -------------------------------------------------------------------------
# This script uses the binary files generated from the SI3D simulations for
# 3D outputs and saves the data as vts files to visualize the results in
# ParaView. It was used previous codes created by Alicia Cortes
# and Francisco Rueda and the script uses gridToVTK.

# NOTES:

# 1. The paraview files will be named with the format 'si3D_xxx.vts', where
# the xxx will be numerical values representing the hours after the
# simulation started.

# 2. It is imperative that Si3D was run using the same values for iht and
# ipxml. It is also imperative that the run was done using the same itspfh
# and itspf parameters from the input file.

# 3. The function exports:
#   a. Individual files for each time step saved in the 3D binary file
#   b. File to load into paraview with the time information, in seconds, of
#       each time step. The file name is si3d.pvd
#   c. File as reference for dates and time in seconds. The file name is ParaviewRef.txt

# Author: Sergio Valbuena
# Date: 03-10-2024

import numpy as np
import pandas as pd
import os
import time
from datetime import datetime, timedelta
import vtk
from vtk.util.numpy_support import numpy_to_vtk

def is_eof(f):
    cur = f.tell()    # save current position
    f.seek(0, os.SEEK_END)
    end = f.tell()    # find the size of file
    f.seek(cur, os.SEEK_SET)
    return cur == end


def si3D_to_paraview(pathfile, pathsave, startdate, deltaZ, dx, dz, dt, iTurb, itspf, nTracer, concTr):
    # Function to convert SI3D binary files to vtk for ParaView visualization

    # Constants
    FileName3D = 'si3d_3D'
    PlaneName = 'plane_2'
    outputFile = 'si3d'
    FileNameZ = 'si3d_layer.txt'
    fileTracer = 'tracer_'

    # Move to working directory
    os.chdir(pathfile)

    beg1 = time.time()

    # Read Depth information
    if deltaZ:
        M = np.loadtxt(FileNameZ, skiprows=5)
        layer = M[:, 0]
        depth = M[:, 1]
        depth[depth == -100] = 0  # Adjust as needed
        ddz = dz
    else:
        ddz = dz
    del dz

    # Creation of the .pvd file to add time to the series of paraview files
    os.chdir(pathsave)
    # Create and open the .pvd file
    fidPV = open(outputFile + '.pvd', 'wt+')
    # Write the initial lines to the file
    fidPV.write('%s\n' % '<?xml version="1.0"?>')
    fidPV.write('%s\n' % '<VTKFile type="Collection" version="1.0" byte_order="LittleEndian">')
    fidPV.write('%s\n' % '\t<Collection>')

    # Reading binary files
    os.chdir(pathfile)
    fid3D = open(FileName3D, 'rb')
    fidPL = open(PlaneName, 'rb')

    _ = np.fromfile(fid3D, count=1, dtype='int32')
    _ = np.fromfile(fidPL, count=1, dtype='int32')

    # Read the number of frames
    n_frames3D = np.fromfile(fid3D, count=1, dtype='int32')[0]
    n_framesPL = np.fromfile(fidPL, count=1, dtype='int32')[0]

    _ = np.fromfile(fidPL, count=1, dtype='int32')
    _ = np.fromfile(fidPL, count=1, dtype='int32')
    ipointsPL = np.fromfile(fidPL, count=1, dtype='int32')[0]
    _ = np.fromfile(fidPL, count=1, dtype='int32')

    _ = np.fromfile(fid3D, count=1, dtype='int32')
    _ = np.fromfile(fid3D, count=1, dtype='int32')
    ipoints3D = np.fromfile(fid3D, count=1, dtype='int32')[0]
    _ = np.fromfile(fid3D, count=1, dtype='int32')

    if nTracer > 0:
        fidTr = []
        n_framesTr = np.full(nTracer, np.nan)
        ipointsTr = np.full(nTracer, np.nan)
        for tr in range(0, nTracer):
            filenameTr = fileTracer + str(tr + 1)
            fidTr.append(open(filenameTr, 'rb'))
            _ = np.fromfile(fidTr[tr], count=1, dtype='int32')
            n_framesTr[tr] = np.fromfile(fidTr[tr], count=1, dtype='int32')[0]
            if (n_frames3D != n_framesPL) or (n_frames3D != n_framesTr[tr]):
                raise ValueError('The number of frames between the surface plane and 3D file are not the same')
            else:
                n_frames = n_frames3D
                _ = np.fromfile(fidTr[tr], count=1, dtype='int32')
                _ = np.fromfile(fidTr[tr], count=1, dtype='int32')
                ipointsTr[tr] = int(np.fromfile(fidTr[tr], count=1, dtype='int32')[0])
                _ = np.fromfile(fidTr[tr], count=1, dtype='int32')
    else:
        if (n_frames3D != n_framesPL):
            raise ValueError('The number of frames between the surface plane and 3D file are not the same')
        else:
            n_frames = n_frames3D

    if itspf == 0:
        n_frames += 1
    elif itspf != 0:
        n_frames += 2

    istep = np.full(n_frames, 0)
    year1 = np.full(n_frames, np.nan)
    month1 = np.full(n_frames, np.nan)
    day1 = np.full(n_frames, np.nan)
    hour1 = np.full(n_frames, np.nan)

    for n in range(0, n_frames):
        beg2 = time.time()
        _ = np.fromfile(fid3D, count=1, dtype='int32')
        _ = np.fromfile(fidPL, count=1, dtype='int32')

        st3d = is_eof(fid3D)
        stpl = is_eof(fidPL)
        sttr = np.full(nTracer, np.nan)

        if nTracer > 0:
            for tr in range(0, nTracer):
                _ = np.fromfile(fidTr[tr], count=1, dtype='int32')
                sttr[tr] = is_eof(fidTr[tr])
        if (st3d == 0) or (stpl == 0) or (np.sum(sttr) == 0):
            istep[n] = np.fromfile(fid3D, count=1, dtype='int32')
            year1[n] = np.fromfile(fid3D, count=1, dtype='int32')
            month1[n] = np.fromfile(fid3D, count=1, dtype='int32')
            day1[n] = np.fromfile(fid3D, count=1, dtype='int32')
            hour1[n] = np.fromfile(fid3D, count=1, dtype='int32')

            _ = np.fromfile(fidPL, count=1, dtype='int32')
            _ = np.fromfile(fidPL, count=1, dtype='int32')
            _ = np.fromfile(fidPL, count=1, dtype='int32')
            _ = np.fromfile(fidPL, count=1, dtype='int32')
            _ = np.fromfile(fidPL, count=1, dtype='int32')

            if nTracer > 0:
                for tr in range(0, nTracer):
                    _ = np.fromfile(fidTr[tr], count=1, dtype='int32')
                    _ = np.fromfile(fidTr[tr], count=1, dtype='int32')
                    _ = np.fromfile(fidTr[tr], count=1, dtype='int32')
                    _ = np.fromfile(fidTr[tr], count=1, dtype='int32')
                    _ = np.fromfile(fidTr[tr], count=1, dtype='int32')
            if n == 0:
                if iTurb == 1:
                    outarr3D = np.fromfile(fid3D, count=13 * ipoints3D, dtype='float32')
                    x = outarr3D[0:len(outarr3D) - 12:13].astype(int)
                    y = outarr3D[1:len(outarr3D) - 11:13].astype(int)
                    z = outarr3D[2:len(outarr3D) - 10:13].astype(int)
                    # h = outarr3D[3:len(outarr3D) - 9:13]
                    u = outarr3D[4:len(outarr3D) - 8:13]
                    v = outarr3D[5:len(outarr3D) - 7:13]
                    w = outarr3D[6:len(outarr3D) - 6:13]
                    T = outarr3D[7:len(outarr3D) - 5:13]
                    Dv = outarr3D[8:len(outarr3D) - 4:13]
                    q2 = outarr3D[9:len(outarr3D) - 3:13]
                    q2l = outarr3D[10:len(outarr3D) - 2:13]
                    kh = outarr3D[11:len(outarr3D) - 1:13]
                    Av = outarr3D[12:len(outarr3D):13]
                else:
                    outarr3D = np.fromfile(fid3D, count=8 * ipoints3D, dtype='float32')
                    x = outarr3D[0:len(outarr3D) - 7:8].astype(int)
                    y = outarr3D[1:len(outarr3D) - 6:8].astype(int)
                    z = outarr3D[2:len(outarr3D) - 5:8].astype(int)
                    # h = outarr3D[3:len(outarr3D) - 4:8]
                    u = outarr3D[4:len(outarr3D) - 3:8]
                    v = outarr3D[5:len(outarr3D) - 2:8]
                    w = outarr3D[6:len(outarr3D) - 1:8]
                    T = outarr3D[7:len(outarr3D):8]

                outarrPL = np.fromfile(fidPL, count=8 * ipointsPL, dtype='float32')
                dumPLx = outarrPL[0:len(outarrPL) - 7:8].astype(int)
                dumPLy = outarrPL[1:len(outarrPL) - 6:8].astype(int)
                s = outarrPL[7:len(outarrPL):8]

                if nTracer > 0:
                    tracerc = np.full((int(ipointsTr[0]), nTracer), np.nan)
                    for tr in range(0, nTracer):
                        outarrTr = np.fromfile(fidTr[tr], count=4 * int(ipointsTr[tr]), dtype='float32')
                        # dumTrx = outarrTr[0:len(outarrTr) - 3:4].astype(int)
                        # dumTry = outarrTr[1:len(outarrTr) - 2:4].astype(int)
                        # dumTrz = outarrTr[2:len(outarrTr) - 1:4].astype(int)
                        tracerc[:, tr] = outarrTr[3:len(outarrTr):4]
            else:
                if iTurb == 1:
                    outarr3D = np.fromfile(fid3D, count=10 * ipoints3D, dtype='float32')
                    # h = outarr3D[0:len(outarr3D) - 9:10]
                    u = outarr3D[1:len(outarr3D) - 8:10]
                    v = outarr3D[2:len(outarr3D) - 7:10]
                    w = outarr3D[3:len(outarr3D) - 6:10]
                    T = outarr3D[4:len(outarr3D) - 5:10]
                    Dv = outarr3D[5:len(outarr3D) - 4:10]
                    q2 = outarr3D[6:len(outarr3D) - 3:10]
                    q2l = outarr3D[7:len(outarr3D) - 2:10]
                    kh = outarr3D[8:len(outarr3D) - 1:10]
                    Av = outarr3D[9:len(outarr3D):10]
                else:
                    outarr3D = np.fromfile(fid3D, count=5 * ipoints3D, dtype='float32')
                    # h = outarr3D[0:len(outarr3D) - 4:5]
                    u = outarr3D[1:len(outarr3D) - 3:5]
                    v = outarr3D[2:len(outarr3D) - 2:5]
                    w = outarr3D[3:len(outarr3D) - 1:5]
                    T = outarr3D[4:len(outarr3D):5]

                outarrPL = np.fromfile(fidPL, count=6 * ipointsPL, dtype='float32')
                s = outarrPL[5:len(outarrPL):6]

                if nTracer > 0:
                    tracerc = np.full((int(ipointsTr[0]), nTracer), np.nan)
                    for tr in range(0, nTracer):
                        outarrTr = np.fromfile(fidTr[tr], count=1 * int(ipointsTr[tr]), dtype='float32')
                        tracerc[:, tr] = outarrTr[0:len(outarrTr):1]

            # Code Generator
            xmin = np.min(x)
            xmax = np.max(x)
            ymin = np.min(y)
            ymax = np.max(y)
            zmin = np.min(z)
            zp = np.unique(z)

            # Adding results of Grad(u,v) = 0
            iy = y == ymin
            filtx = x[iy]
            filtz = z[iy]
            idatasum = np.sum(iy)
            ynew = np.concatenate((y, (ymin - 1) * np.ones(idatasum)))
            znew = np.concatenate((z, filtz))
            xnew = np.concatenate((x, filtx))
            filtu = u[iy]
            filtv = v[iy]
            filtw = w[iy]
            filtT = T[iy]
            # filth = h[iy]
            # h = np.concatenate((h, filth))
            u = np.concatenate((u, filtu))
            v = np.concatenate((v, filtv))
            w = np.concatenate((w, filtw))
            T = np.concatenate((T, filtT))

            if nTracer > 0:
                tracerc1 = np.full((len(xnew), nTracer), np.nan)
                for tr in range(0, nTracer):
                    filtTr = tracerc[iy, tr]
                    tracerc1[:, tr] = np.concatenate((tracerc[:, tr], filtTr))
                tracerc = tracerc1
                del tracerc1

            if iTurb == 1:
                filtDV = Dv[iy]
                filtq2 = q2[iy]
                filtq2l = q2l[iy]
                filtkh = kh[iy]
                filtAv = Av[iy]
                Dv = np.concatenate((Dv, filtDV))
                q2 = np.concatenate((q2, filtq2))
                q2l = np.concatenate((q2l, filtq2l))
                kh = np.concatenate((kh, filtkh))
                Av = np.concatenate((Av, filtAv))

            ix = xnew == xmin
            filty = ynew[ix]
            filtz = znew[ix]
            idatasum = np.sum(ix)
            xnew = np.concatenate((xnew, (xmin - 1) * np.ones(idatasum)))
            znew = np.concatenate((znew, filtz))
            ynew = np.concatenate((ynew, filty))
            filtu = u[ix]
            filtv = v[ix]
            filtw = w[ix]
            filtT = T[ix]
            # filth = h[ix]
            # h = np.concatenate((h, filth))
            u = np.concatenate((u, filtu))
            v = np.concatenate((v, filtv))
            w = np.concatenate((w, filtw))
            T = np.concatenate((T, filtT))

            if nTracer > 0:
                tracerc1 = np.full((len(xnew), nTracer), np.nan)
                for tr in range(0, nTracer):
                    filtTr = tracerc[ix, tr]
                    tracerc1[:, tr] = np.concatenate((tracerc[:, tr], filtTr))
                tracerc = tracerc1
                del tracerc1
            if iTurb == 1:
                filtDV = Dv[ix]
                filtq2 = q2[ix]
                filtq2l = q2l[ix]
                filtkh = kh[ix]
                filtAv = Av[ix]
                Dv = np.concatenate((Dv, filtDV))
                q2 = np.concatenate((q2, filtq2))
                q2l = np.concatenate((q2l, filtq2l))
                kh = np.concatenate((kh, filtkh))
                Av = np.concatenate((Av, filtAv))

            iz = znew == zmin
            filty = ynew[iz]
            filtx = xnew[iz]
            idatasum = np.sum(iz)
            xnew = np.concatenate((xnew, filtx))
            ynew = np.concatenate((ynew, filty))
            znew = np.concatenate((znew, (zmin - 1) * np.ones(idatasum)))
            filtu = u[iz]
            filtv = v[iz]
            filtw = w[iz]
            filtT = T[iz]
            # filth = h[iz]
            # h = np.concatenate((h, filth))
            u = np.concatenate((u, filtu))
            v = np.concatenate((v, filtv))
            w = np.concatenate((w, filtw))
            T = np.concatenate((T, filtT))

            if nTracer > 0:
                tracerc1 = np.full((len(xnew), nTracer), np.nan)
                for tr in range(0, nTracer):
                    filtTr = tracerc[iz, tr]
                    tracerc1[:, tr] = np.concatenate((tracerc[:, tr], filtTr))
                tracerc = tracerc1
                del tracerc1
            if iTurb == 1:
                filtDV = Dv[iz]
                filtq2 = q2[iz]
                filtq2l = q2l[iz]
                filtkh = kh[iz]
                filtAv = Av[iz]
                Dv = np.concatenate((Dv, filtDV))
                q2 = np.concatenate((q2, filtq2))
                q2l = np.concatenate((q2l, filtq2l))
                kh = np.concatenate((kh, filtkh))
                Av = np.concatenate((Av, filtAv))

            if n == 0:
                xnew = xnew.astype(int)
                ynew = ynew.astype(int)
                znew = znew.astype(int)
                xmin = int(np.min(xnew))
                xmax = int(np.max(xnew))
                ymin = int(np.min(ynew))
                ymax = int(np.max(ymax))
                zp = np.unique(znew).astype(int)

                # To create 3D domain for which there is a numerical solution
                xgrid = np.arange(xmin, xmax + 1)
                ygrid = np.arange(ymin, ymax + 1)
                xg, yg, zg = np.meshgrid(xgrid, ygrid, zp)
                xv = xg.flatten()
                yv = yg.flatten()
                zv = zg.flatten()

                # To define the actual dimension of the matrix that describes the domain
                # of the numerical solution including and extra level for the bottom.
                if deltaZ:
                    layer = np.concatenate(([1], layer[:]))
                    depth = np.concatenate((depth[:], [depth[-1] + ddz]))
                    zp = -depth
                else:
                    zp = -(zp - 1) * ddz
                xp = (xg[0, :, 0] - 1) * dx
                yp = (yg[:, 0, 0] - 1) * dx

                # To create the 3D domain with the real dimension of the structured grid
                xgf, ygf, zgf = np.meshgrid(xp, yp, zp)
                nx, ny, nz = np.shape(xgf)
                xgf = xgf.ravel(order='C')
                ygf = ygf.ravel(order='C')
                zgf = zgf.ravel(order='C')

                del xp, yp, zp, xg, yg, zg

                # To create code for the surface plane
                coord_simPL = np.column_stack((dumPLx, dumPLy))
                codePL = np.array(['A{}A{}A2'.format(row[0], row[1]) for row in coord_simPL])

                coord_sim = np.column_stack((xnew, ynew, znew))
                codex = np.array(['A{}A{}A{}'.format(row[0], row[1], row[2]) for row in coord_sim])
                coord_dom = np.column_stack((xv, yv, zv))
                codev = np.array(['A{}A{}A{}'.format(row[0], row[1], row[2]) for row in coord_dom])

                _, icodev, icodex = np.intersect1d(codev, codex, return_indices=True)
                _, icodevPL, icodexPL = np.intersect1d(codev, codePL, return_indices=True)

                if (len(icodev) != len(xnew)) or (len(icodex) != len(xnew)):
                    raise ValueError('The code generator is not working properly. The length of the sum of idata must be the same length as the vectors with the solution of si3D')

            if iTurb == 1:
                TKE = q2 / 2
                ml = q2l / q2

            _ = np.fromfile(fid3D, count=1, dtype='int32')
            _ = np.fromfile(fidPL, count=1, dtype='int32')
            if nTracer > 0:
                for tr in range(0, nTracer):
                    _ = np.fromfile(fidTr[tr], count=1, dtype='int32')

            del outarrPL, outarr3D
            if nTracer > 0:
                del outarrTr

            # Paraview file create .vtk
            lv = np.full(len(xgf.ravel()), np.nan)
            sv = np.full(len(xgf.ravel()), np.nan)
            uv = np.full(len(xgf.ravel()), np.nan)
            vv = np.full(len(xgf.ravel()), np.nan)
            wv = np.full(len(xgf.ravel()), np.nan)
            Tv = np.full(len(xgf.ravel()), np.nan)

            lv[icodev] = sv[icodex]
            lv[icodevPL] = s[icodexPL]
            uv[icodev] = u[icodex]
            vv[icodev] = v[icodex]
            wv[icodev] = w[icodex]
            Tv[icodev] = T[icodex]

            if iTurb == 1:
                Dvv = np.full(len(xgf.ravel()), np.nan)
                TKEv = np.full(len(xgf.ravel()), np.nan)
                mlv = np.full(len(xgf.ravel()), np.nan)
                khv = np.full(len(xgf.ravel()), np.nan)
                Avv = np.full(len(xgf.ravel()), np.nan)

                Dvv[icodev] = Dv[icodex] / (100**2)
                TKEv[icodev] = TKE[icodex]
                mlv[icodev] = ml[icodex]
                khv[icodev] = kh[icodex]
                Avv[icodev] = Av[icodex]

            # To convert into float32
            uv = uv.astype('float32')
            vv = vv.astype('float32')
            wv = wv.astype('float32')
            Tv = Tv.astype('float32')
            xgf = xgf.astype('float32')
            ygf = ygf.astype('float32')
            zgf = zgf.astype('float32')

            data = {'xgf': xgf,
                    'ygf': ygf,
                    'zgf': zgf}
            df = pd.DataFrame(data)
            df['Tv'] = Tv
            df['uv'] = uv
            df['vv'] = vv
            df['wv'] = wv

            if iTurb == 1:
                Dvv = Dvv.astype('float32')
                TKEv = TKEv.astype('float32')
                mlv = mlv.astype('float32')
                khv = khv.astype('float32')
                Avv = Avv.astype('float32')
                df['Dvv'] = Dvv
                df['TKEv'] = TKEv
                df['mlv'] = mlv
                df['khv'] = khv
                df['Avv'] = Avv
            if nTracer > 0:
                for tr in range(0, nTracer):
                    title = concTr[tr]
                    conc = np.full(len(xgf.ravel()), np.nan)
                    conc[icodev] = tracerc[icodex, tr]
                    df[title] = conc
            
            df = df.sort_values(by=['zgf','xgf','ygf'])
            df = df.reset_index(drop=True)

            xg = df['xgf'].values
            yg = df['ygf'].values
            zg = df['zgf'].values
            Tv = df['Tv'].values
            uv = df['uv'].values
            vv = df['vv'].values
            wv = df['wv'].values
            # To save data into .vts file for visualization in paraview
            os.chdir(pathsave)
            outputname = outputFile + '_' + str(round(istep[n] * dt / 3600, 2))
            
            fidPV.write('%s' % '\t\t<DataSet timestep="' + str(istep[n] * dt) + '" file="' + outputname + '.vts"/>\n')

            dims = (nx, ny, nz)
            points = np.c_[xg, yg, zg]
            vtk_points = vtk.vtkPoints()
            vtk_points.SetData(numpy_to_vtk(points, deep=True))
            grid = vtk.vtkStructuredGrid()
            grid.SetDimensions(nx, ny, nz)
            grid.SetPoints(vtk_points)

            T_vtk = numpy_to_vtk(Tv, deep=True)
            T_vtk.SetName("T(C)")
            grid.GetPointData().AddArray(T_vtk)

            V = np.column_stack((uv, vv, wv))
            V = V.flatten()
            vel = numpy_to_vtk(V, deep=True, array_type=vtk.VTK_FLOAT)
            V_vtk = vtk.vtkFloatArray()
            V_vtk.SetNumberOfComponents(3)
            V_vtk.SetName('u(m/s)')
            V_vtk.SetArray(vel, len(V), vtk.VTK_FLOAT)
            grid.GetPointData().AddArray(V_vtk)

            if iTurb == 1:
                Dv_vtk = numpy_to_vtk(df['Dvv'].values, deep=True)
                Dv_vtk.SetName('Dv(m2/s)')
                grid.GetPointData().AddArray(Dv_vtk)
                TKE_vtk = numpy_to_vtk(df['TKEv'].values, deep=True)
                TKE_vtk.SetName('TKE(m2/s2)')
                grid.GetPointData().AddArray(TKE_vtk)
                ml_vtk = numpy_to_vtk(df['mlv'].values, deep=True)
                ml_vtk.SetName('ml(m)')
                grid.GetPointData().AddArray(ml_vtk)
                kh_vtk = numpy_to_vtk(df['khv'].values, deep=True)
                kh_vtk.SetName('kh(m2/s)')
                grid.GetPointData().AddArray(kh_vtk)
                Av_vtk = numpy_to_vtk(df['Avv'].values, deep=True)
                Av_vtk.SetName('Av(m2/s)')
                grid.GetPointData().AddArray(Av_vtk)
            if nTracer > 0:
                for tr in range(0, nTracer):
                    title = concTr[tr]
                    conc = df[title].values
                    tr_vtk = numpy_to_vtk(conc, deep=True)
                    tr_vtk.SetName(title)
                    grid.GetPointData().AddArray(tr_vtk)

            outputname += '.vts'
            writer = vtk.vtkXMLStructuredGridWriter()
            writer.SetFileName(outputname)
            writer.SetInputData(grid)
            writer.Write()

            del df, data

            # To create files using a different library

            # vectors_pt = {'u(m/s)': (uv, vv, wv)}
            # vectors_cell = {'u(m/s)': (uc, vc, wc)}
            # scalars_pt = {'T(C)': Tv, 'l(m)': lv}
            # scalars_cell = {'T(C)': Tc, 'l(m)': lc}

            # if nTracer > 0:
            #     for tr in range(0, nTracer):
            #         conc = np.full(len(xgf.ravel()), np.nan)
            #         title = concTr[tr]
            #         conc[icodev] = tracerc[icodex, tr]
            #         tr_C = conc.reshape(np.shape(xgf))
            #         conc_c = np.full((nx - 1, ny - 1, nz - 1), np.nan)
            #         conc_c[:, :, :] = tr_C[1:, 1:, 1:]
            #         scalars_pt.update({title: tr_C})
            #         scalars_cell.update({title: conc_c})
            # if iTurb == 1:
            #     scalars_pt.update({'Dv(m2/s)': Dvv, 'TKE(m2/s)': TKEv, 'ml(m)': mlv, 'kh(m2/s)': khv, 'Av(m2/s)': Avv})
            #     scalars_cell.update({'Dv(m2/s)': Dvc, 'TKE(m2/s)': TKEc, 'ml(m)': mlc, 'kh(m2/s)': khc, 'Av(m2/s)': Avc})

            # pointinfo = {}
            # cellinfo = {}
            # pointinfo.update(vectors_pt)
            # pointinfo.update(scalars_pt)
            # cellinfo.update(vectors_cell)
            # cellinfo.update(scalars_cell)
            # gridToVTK(outputname, xgf, ygf, zgf, pointData=pointinfo, cellData=cellinfo)
            # gridToVTK(outputname, xgf, ygf, zgf, cellData=cellinfo)
            # gridToVTK(outputname, xgf, ygf, zgf, pointData=pointinfo)

        else:
            n_frames = n - 1

        end2 = time.time()
        print('Time to create vtk file for time frame ' + str(n) + ' is ' + str(round(end2 - beg2, 3)) + ' seconds')

    fid3D.close()
    fidPL.close()
    if nTracer > 0:
        for tr in range(0, nTracer):
            fidTr[tr].close()

    fidPV.write('%s\n' % '\t</Collection>')
    fidPV.write('%s' % '</VTKFile>')

    # Creation of paraview reference file
    os.chdir(pathsave)
    startdate = datetime.strptime(startdate, '%Y-%m-%d %H:%M:%S')
    Date = [startdate + timedelta(seconds=i * dt) for i in istep]
    datestr = [date.strftime('%Y-%m-%d %H:%M:%S') for date in Date]
    Id = list(range(1, len(Date) + 1))

    ParaviewRef = pd.DataFrame({
        'Id': Id,
        'seconds': [i * dt for i in istep],
        'Date': Date,
        'DateStr': datestr
    })

    ParaviewRef.to_csv('ParaviewRef.txt', sep='\t', index=False)

    end1 = time.time()
    print('Time needed to run all the code is ' + str(round((end1 - beg1) / 60, 2)) + ' minutes')

    return n_frames

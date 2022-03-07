% NEED TO PUT TEXT TO EXPLAIN THE PROCESSOR.

clc;
close all;
clear variables;

% root = 'S:\si3D\';
% FileSim = 'si3dTest';
root = 'G:\My Drive\Lake_Tahoe\Projects\TahoeWarningWebsite\';
FileSims = {'Turb_0_itspf_0'};

%% Code section
for i = 1:length(FileSims)
    FileSim = char(FileSims(i));
    PathFile = [root,FileSim];
    PathSave = [root,FileSim,'\Paraview'];

    % To obtain simulation information like dt,dz,dx and date
    cd(PathFile)
    inputFile = readlines('si3d_inp.txt');
    year = inputFile(6);
    month = inputFile(7);
    day = inputFile(8);
    hour = inputFile(9);
    year = char(regexp(year,'\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+','match'));
    month = char(regexp(month,'\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+','match'));
    day = char(regexp(day,'\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+','match'));
    hour = char(regexp(hour,'\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+','match'));
    minute = hour(3:4);
    hour = hour(1:2);
    StartDate = [year,'-',month,'-',day,' ',hour,':',minute,':00'];
    
    % Obtain dx from si3d_inp.txt
    dx = inputFile(17);
    dx = char(dx);
    dx = str2double(dx(15:34));
    % Obtain dz from si3d_inp.txt
    dz = inputFile(19);
    dz = char(dz);
    dz = str2double(dz(15:34));
    % Obtain dt from si3d_inp.txt
    dt = inputFile(23);
    dt = char(dt);
    dt = str2double(dt(15:34));
    
    % Obtain idz from si3d_inp.txt and know if si3d_layer.txt is necessary
    idz = inputFile(24);
    idz = char(idz);
    idz = idz(15:34);
    idz = str2double(idz);
    if idz == 0
        DeltaZ = 'constant';
    elseif idz == -1
        DeltaZ = 'variable';
    end
    % Obtain iTurb from si3d_inp.txt. This indicates whether ptrack file
    % has turbulent parameters
    iTurb = inputFile(76);
    iTurb = char(iTurb);
    iTurb = str2double(iTurb(15:34));
    % Obtain itspf from si3d_inp.txt. This indicates whether time steps
    % saved started from the start of the simulation or after a certain
    % time step
    itspf = inputFile(75);
    itspf = char(itspf);
    itspf = str2double(itspf(15:34));
    
    %% To generate Paraview files
    if isfile('ptrack_hydro.bnr')
        % File exists.
        n_frames = Si3DtoParaview(PathFile,PathSave,StartDate,DeltaZ,dx,dz,dt,iTurb,itspf);
%         cd(PathFile)
%         delete 'ptrack_hydro.bnr'
    else
        disp('Error 3D file was not found')
    end
end
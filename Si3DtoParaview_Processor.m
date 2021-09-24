lines
% To process all the simulations ran but not previously processed.

clc;
close all;
clear variables;

% root = 'S:\si3D\';
% FileSim = 'si3dTest';
root = 'S:\si3D\si3D_EcoModeling\';
FileSims = {'EcoModelSpring_R14'};
iTurb = 0;

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

    dx = inputFile(17);
    dx = char(dx);
    dx = str2double(dx(15:34));
    dz = inputFile(19);
    dz = char(dz);
    dz = str2double(dz(15:34));
    dt = inputFile(23);
    dt = char(dt);
    dt = str2double(dt(15:34));

    idz = inputFile(24);
    idz = char(idz);
    idz = idz(15:34);
    idz = str2double(idz);
    if idz == 0
        DeltaZ = 'constant';
    elseif idz == -1
        DeltaZ = 'variable';
    end
    %% To generate Paraview files
    if isfile('ptrack_hydro.bnr')
        % File exists.
        n_frames = Si3DtoParaview(PathFile,PathSave,StartDate,DeltaZ,dx,dz,dt,iTurb);
        cd(PathFile)
        delete 'ptrack_hydro.bnr'
    else
        disp('Error 3D file was not found')
    end
end

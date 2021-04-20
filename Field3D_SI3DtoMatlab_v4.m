% -------------------------------------------------------------------------
% This script uses the binary files generated from the SI3D simulations for
% 3D outputs and saves the data as matlab structures. The script also has
% the capability of creating .vtk files to visualize the results in
% ParaView, and it finally includes a flag in case it is desired to have 3D
% plots in Matlab as well. It was used previous codes created by Alicia Cortes
% and Francisco Rueda and the script uses a vtkwrite function created following
% submission by  William Thielicke and David Gingras.

% Version 1: - This version transform the binary file result from SI3D into
% a matlab variables that can be saved as matlab data for the 3D domain
% considered in the simulation.
% - This script also creates a DateTime
% variable for plotting purposes.
% Version 2: - This version includes the creation of time step files for
% paraview visualization of the 3D numerical simulations.
% Version 3: Structure modified for RAM memory limitations and to make the
% code less computationally expensive. It was also added the creation of a
% file to load in paraview that will have the time of the simulation in
% seconds. This was done to be able to use the Langrarian Particle Tracking
% filter in Paraview. The Paraview ref file was also modified to save the new time.
% Version 4: The code was modified for the newest version of si3D where the
% 3D output has the turbulent variables. This version not longer considers
% the plot of matlab outputs and focuses only in created files.

% NOTES:

% 1. vtkwrite function is needed IF Paraview output is desired. This script
% will run the section where .vtk files are created

% 2. The resulting matlab data will be squeezed vectors that are the
% elements of a 3D domain of size m1,m2,m3.

% 3. The paraview files will be named with the format 'si3D_xxx.vtk', where
% the xxx will be numerical values representing the hours after the
% simulation started.

% 4. Please consider that the more positive flags the user activate, the
% code will require more RAM and it is possible that it crashes due to
% memory limitations. If RAM memory is limited, it is recommended to run
% one flag at a time.

% 5. It is recommended to run the code from a text editor for less memory
% consumption.

% Author: Sergio Valbuena
% Date: 05-12-2020

clear variables
close all
clc
ticT = tic;

% ----------------------- USER SECTION START-------------------------------
%% Load data and define variables
root = 'S:\';
SimName = 'L50B50_H100_W0.5_N0.006694_f40_0.25Ti';
Pathfilelayer = [root,'si3D\',SimName];
Pathfile = [root,'si3D\',SimName];
PathSave = [root,'\si3D\',SimName,'\Paraview\'];
Path_Bathy = 'G:\My Drive\Lake_Tahoe\Projects\Upwelling_3DModel\Bathymetry';
FileName3D = 'ptrack_hydro.bnr';
PlaneName = 'plane_2';
StartDate = '2018-01-01 00:00:00';
outputFile = 'si3D';

% Please define dx for the numerical simulation
dx = 200; % [m]
% Please define depth of Lake H
H = 100.5;
% Please define if dz is constant or variable in the simulation
DeltaZ = 'constant';
dz = 1;                         % [m] Only needed if DeltaZ is constant
FileNameZ = 'si3d_layer.txt';   % [m] Only needed if DeltaZ is variable

% Please define the dt that the simulation used
% 3D files with lags from the starting of the simulation is not developed yet.
dt = 10;        % [s] Time step in seconds

% ---------------------- USER SECTION END --------------------------------
% --------------------  CODE SECTION START -------------------------------
%% Definition of depths and horizontal dimensions
switch DeltaZ
    case 'variable'
        cd(Pathfilelayer)
        M = dlmread(FileNameZ,'',5,0);
        Layer = M(:,1);
        Depth = M(:,2);
        Depth(Depth==-100) = 0;
        clear M
    case 'constant'
        ddz = dz;
end
clear dz

%% Creation of .pvd file to add time to the series of paraview files
cd(PathSave)
fidPV = fopen([outputFile,'.vtk.series'],'wt+');
fprintf(fidPV,'%s\n','{');
fprintf(fidPV,'%s\n','	"file-series-version" : "1.0",');
fprintf(fidPV,'%s\n',' 	"files" : [');

%% Read binary file from SI3D output
cd(Pathfile)
fid3D = fopen(FileName3D);
fidPL = fopen(PlaneName);
dum3D1 = fread(fid3D,1,'int32');
dumPL1 = fread(fidPL,1,'int32');
n_frames3D = fread(fid3D,1,'int32');
n_framesPL = fread(fidPL,1,'int32');
if n_frames3D ~= n_framesPL
    error('The number of frames between the surface plane and 3D file are not the same')
elseif n_frames3D == n_framesPL
    n_frames = n_frames3D;
end
dumPL2 = fread(fidPL,1,'int32');
dumPL3 = fread(fidPL,1,'int32');
ipointsPL = fread(fidPL,1,'int32');
dumPL4 = fread(fidPL,1,'int32');

dum3D2 = fread(fid3D,1,'int32');
dum3D3 = fread(fid3D,1,'int32');
ipoints3D = fread(fid3D,1,'int32');
dum3D4 = fread(fid3D,1,'int32');

istep = zeros(n_frames,1);
year1 = zeros(n_frames,1);
month1 = zeros(n_frames,1);
day1 = zeros(n_frames,1);
hour1 = zeros(n_frames,1);

h = zeros(ipoints3D,1);
u = zeros(ipoints3D,1);
v = zeros(ipoints3D,1);
w = zeros(ipoints3D,1);
Dv = zeros(ipoints3D,1);
T = zeros(ipoints3D,1);
q2 = zeros(ipoints3D,1);
kh = zeros(ipoints3D,1);
q2l = zeros(ipoints3D,1);
s = zeros(ipointsPL,1);


tic;
for count = 1:n_frames+1
    dum3D5 = fread(fid3D,1,'int32');
    dumPL5 = fread(fidPL,1,'int32');
    st = feof(fid3D);
    if (st == 0)
        istep(count)=fread(fid3D,1,'int32');
        year1(count)=fread(fid3D,1,'int32');
        month1(count)=fread(fid3D,1,'int32');
        day1(count)=fread(fid3D,1,'int32');
        hour1(count)=fread(fid3D,1,'float32');
        dumPL6=fread(fidPL,1,'int32');
        dumPL7=fread(fidPL,1,'int32');
        dumPL8=fread(fidPL,1,'int32');
        dumPL9=fread(fidPL,1,'int32');
        dumPL10=fread(fidPL,1,'float32');
        % ... Read all data for present time slice
        if (count == 1)
            out_array3D=fread(fid3D,13*ipoints3D,'float32');
            x = out_array3D(1:13:length(out_array3D)-12);
            y = out_array3D(2:13:length(out_array3D)-11);
            z = out_array3D(3:13:length(out_array3D)-10);
            h = out_array3D(4:13:length(out_array3D)-9);
            u = out_array3D(5:13:length(out_array3D)-8);
            v = out_array3D(6:13:length(out_array3D)-7);
            w = out_array3D(7:13:length(out_array3D)-6);
            Dv = out_array3D(8:13:length(out_array3D)-5);
            T = out_array3D(9:13:length(out_array3D)-4);
            q2 = out_array3D(10:13:length(out_array3D)-3);
            q2l = out_array3D(11:13:length(out_array3D)-2);
            kh = out_array3D(12:13:length(out_array3D)-1);
            Av = out_array3D(13:13:length(out_array3D));
            
            out_arrayPL=fread(fidPL,8*ipointsPL,'float32');
            dumPLx = out_arrayPL(1:8:length(out_arrayPL)-7);
            dumPLy = out_arrayPL(2:8:length(out_arrayPL)-6);
            dumPLu = out_arrayPL(3:8:length(out_arrayPL)-5);
            dumPLv = out_arrayPL(4:8:length(out_arrayPL)-4);
            dumPLw = out_arrayPL(5:8:length(out_arrayPL)-3);
            dumPLT = out_arrayPL(6:8:length(out_arrayPL)-2);
            dumPLAz = out_arrayPL(7:8:length(out_arrayPL)-1);
            s = out_arrayPL(8:8:length(out_arrayPL));
            
            %% Code generator
            % To select the unique values of the array results that compound
            % the matrix of the numerical solution
            xmin = min(x);
            xmax = max(x);
            ymin = min(y);
            ymax = max(y);
            zp = unique(z);
            
            % To create the 3D domain for which there is a numerical solution
            [xg,yg,zg] = meshgrid(xmin:xmax,ymin:ymax,zp);
            xv = xg(:);
            yv = yg(:);
            zv = zg(:);
            % To define the actual dimensions of the matrix that describes the
            % domain of the numerical solution including an extra level for bottom.
            switch DeltaZ
                case 'variable'
                    if length(Layer)~= length(zp)
                        idata = ~ismember(Layer,zp);
                        Layer(idata) = [];
                        Depth(idata) = [];
                    end
                    idata = zp == Layer;
                    zp = -Depth(idata);
                case 'constant'
                    zp = -(zp-2)*ddz;
            end
            xp = xg(1,:,1)'*dx;
            yp = yg(:,1,1)*dx;
            % To create the 3D domain with the real dimensions of the
            % structured grid
            [xgf,ygf,zgf] = meshgrid(xp,yp,zp);
            [m1,m2,m3] = size(xgf);
            clear xp yp zp xg yg zg
            
            
            % To create code for the surface plane
            [xp,yp] = meshgrid(min(dumPLx):max(dumPLx),min(dumPLy):max(dumPLy));
            [m1p,m2p] = size(xp);
            coord_simPL = [dumPLx,dumPLy];
            codePL = num2str(coord_simPL);
            codecellPL = cellstr(codePL);
            clearvars codePL
            codePL = regexprep(codecellPL,'\s+','A');
            codePL = strcat(codePL,'A2');
            
            tic;
            % To create a code that idenfies the indexes of the numerical
            % simulations within the whole domain.
            % The numerical values are converted into strings for code
            % identification.
            coord_sim = [x,y,z];
            coord_dom = [xv,yv,zv];
            codex = num2str(coord_sim);
            codev = num2str(coord_dom);
            
            % clearvars xv yv zv
            codecellx = cellstr(codex);
            codecellv = cellstr(codev);
            clearvars codex codev
            codex = regexprep(codecellx,'\s+','A');
            codev = regexprep(codecellv,'\s+','A');
            
            [a,icodev,icodex] = intersect(codev,codex,'sorted');
            [b,icodevPL,icodexPL] = intersect(codev,codePL,'sorted');
            sum1 = length(icodev);
            sum2 = length(icodex);
            clearvars codecellx codecellv codev codex
            disp(['Time needed to generate the code for indexing data is ',num2str(toc),' seconds'])
            
            if sum1 ~= length(x) || sum2 ~= length(x)
                error('The code generator is not working properly. The lenght of the sum of idata must be the same length as the vectors with the solution of SI3D')
            end
        else
            out_array3D=fread(fid3D,10*ipoints3D,'float32');
            h = out_array3D(1:10:length(out_array3D)-9);
            u = out_array3D(2:10:length(out_array3D)-8);
            v = out_array3D(3:10:length(out_array3D)-7);
            w = out_array3D(4:10:length(out_array3D)-6);
            Dv = out_array3D(5:10:length(out_array3D)-5);
            T = out_array3D(6:10:length(out_array3D)-4);
            q2 = out_array3D(7:10:length(out_array3D)-3);
            q2l = out_array3D(8:10:length(out_array3D)-2);
            kh = out_array3D(9:10:length(out_array3D)-1);
            Av = out_array3D(10:10:length(out_array3D));
            
            out_arrayPL=fread(fidPL,6*ipointsPL,'float32');
            dumPLu = out_arrayPL(1:6:length(out_arrayPL)-5);
            dumPLv = out_arrayPL(2:6:length(out_arrayPL)-4);
            dumPLw = out_arrayPL(3:6:length(out_arrayPL)-3);
            dumPLT = out_arrayPL(4:6:length(out_arrayPL)-2);
            dumPLAz = out_arrayPL(5:6:length(out_arrayPL)-1);
            s = out_arrayPL(6:6:length(out_arrayPL));
        end
        TKE = q2./2; % Turbulent Kinectic Energy
        ml = q2l./q2; % Mixing macroscale
        
        dum3D6 = fread(fid3D,1,'int32');
        dumPL11 = fread(fidPL,1,'int32');
        clear out_arrayPL out_array3D;
        
        %% Paraview file creation .vtk
        lv = NaN(length(xgf(:)),1);
        sv = zeros(length(xgf(:)),1);
        uv = NaN(length(xgf(:)),1);
        vv = NaN(length(xgf(:)),1);
        wv = NaN(length(xgf(:)),1);
        Tv = NaN(length(xgf(:)),1);
        Dvv = NaN(length(xgf(:)),1);
        TKEv = NaN(length(xgf(:)),1);
        mlv = NaN(length(xgf(:)),1);
        khv = NaN(length(xgf(:)),1);
        Avv = NaN(length(xgf(:)),1);
        
        lv(icodev) = sv(icodex);
        lv(icodevPL) = s(icodexPL);
        uv(icodev) = u(icodex);
        vv(icodev) = v(icodex);
        wv(icodev) = w(icodex);
        Tv(icodev) = T(icodex);
        Dvv(icodev) = Dv(icodex)./(100^2);
        TKEv(icodev) = TKE(icodex);
        mlv(icodev) = ml(icodex);
        khv(icodev) = kh(icodex);
        Avv(icodev) = Av(icodex);
        
        tic1 = tic;
        cd(PathSave)
        vtkwriteSV_v1([outputFile,'_',num2str((istep(count))*dt/3600),'.vtk'],...
            'structured_grid',xgf-dx,ygf-dx,zgf,'vectors','U(m/s)',uv,vv,wv,...
            'scalars','T(C)',Tv,'scalars','l(m)',lv,'scalars','Dv(m2/s)',Dvv,...
            'scalars','TKE(m2/s2)',TKEv,'scalars','ml(m)',mlv,...
            'scalars','kh(m2/s)',khv,'scalars','Av(m2/s)',Avv,'precision','6','binary')
        disp(['Time frame ',num2str(count),' is ',num2str(toc(tic1)),' seconds'])      
        if count == n_frames + 1
            fprintf(fidPV,'%s\n',['		{ "name" : "',outputFile,'_',num2str((istep(count))*dt/3600),'.vtk",',' "time" : ',num2str((istep(count))*dt),'}']);
        else
            fprintf(fidPV,'%s\n',['		{ "name" : "',outputFile,'_',num2str((istep(count))*dt/3600),'.vtk",',' "time" : ',num2str((istep(count))*dt),'},']);
        end
    else
        n_frames = count-1;
        break
    end
    
end
fclose(fid3D);
fclose(fidPL);

%% Creation of paraview reference file

cd(PathSave)
StartDate = datetime(StartDate,'Format','yyy-MM-dd HH:mm:ss');
Date = StartDate + istep.*dt/3600/24;
DateStr = datestr(Date);
Id = transpose(1:length(Date));
Id = Id - 1;
ParaviewRef = table(Id,istep.*dt,Date,DateStr);
ParaviewRef.Properties.VariableNames = {'Id','Seconds','Date','DateStr'};
writetable(ParaviewRef,'ParaviewRef.txt','Delimiter','tab')

fprintf(fidPV,'%s\n','	]');
fprintf(fidPV,'%s\n','}');
fclose(fidPV);


disp(['Time needed to read the binary file from SI3D and create paraview files is ',num2str(toc),' seconds'])

disp(['Total time needed to create all the files requested is ',num2str(toc(ticT)/60),' minutes'])

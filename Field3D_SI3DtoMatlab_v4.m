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
% SimName = 'psi3D_SVFR_R6';
SimName = 'LakeTahoe_R1';
root = 'S:\';
Path_file_layer = [root,'si3D\',SimName];
Path_file = [root,'si3D\',SimName];
Path_save_struct = [root,'si3D\',SimName,'\MatlabOutputs\'];
Path_save_Paraview = [root,'\si3D\',SimName,'\Paraview\'];
Path_Bathy = 'G:\My Drive\Lake_Tahoe\Projects\Upwelling_3DModel\Bathymetry';
FileName3D = 'ptrack_hydro.bnr';
StartDate = '2018-05-26 00:00:00';
outputFile = 'si3D';

% Please define if Matlab structure is going to be saved
saveMatlabFlag = 'NO';
% Please define if paraview function is going to be used
saveParaviewFlag = 'YES';
% Please define if vorticity will be calculate
vorticityFlag = 'NO';
% Please define type of bathymetry 'Lake' 'Rectangular'
BathyFlag = 'Lake';

% Please define dx for the numerical simulation
dx = 100; % [m]
% Please define depth of Lake H
H = 501;
% Please define if dz is constant or variable in the simulation
DeltaZ = 'variable';
dz = 0.25;                         % [m] Only needed if DeltaZ is constant
FileNameZ = 'si3d_layer.txt';   % [m] Only needed if DeltaZ is variable

% Please define the dt that the simulation used
% 3D files with lags from the starting of the simulation is not developed yet.
dt = 10;        % [s] Time step in seconds

% ---------------------- USER SECTION END --------------------------------
% --------------------  CODE SECTION START -------------------------------
%% Definition of depths and horizontal dimensions
switch DeltaZ
    case 'variable'
        cd(Path_file_layer)
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
switch saveParaviewFlag
    case 'YES'
        cd(Path_save_Paraview)
        fid_PV = fopen([outputFile,'.vtk.series'],'wt+');
        fprintf(fid_PV,'%s\n','{');
        fprintf(fid_PV,'%s\n','	"file-series-version" : "1.0",');
        fprintf(fid_PV,'%s\n',' 	"files" : [');
    case 'NO'
end

%% Reading of wind field and Bathymetry
cd(Path_file)
surfbc = readtable('surfbc.txt');
uwind = surfbc.Var9;
vwind = surfbc.Var10;
wwind = zeros(length(uwind),1);
twind = surfbc.Var1/24;
StartDate = datetime(StartDate,'Format','yyyy-MM-dd HH:mm:ss');
Timew = StartDate:600/86400:StartDate + twind(end);
wind_Dir = atan2d(uwind,vwind);
wind_Dir = wind_Dir + 180;

Timewind = StartDate:3600/86400:StartDate + twind(end);
idata = ismember(Timew,Timewind);
uwindf = uwind(idata);
vwindf = vwind(idata);
wwindf = wwind(idata);
Dirwindf = wind_Dir(idata);
clearvars surfbc wind_Dir uwind vwind wwind twind Timew idata

switch BathyFlag
    case 'Lake'
        cd(Path_Bathy)
        load(['LakeTahoeGrid',num2str(dx),'m'])
        eval(['File = LakeTahoeGrid',num2str(dx),'m;'])
        Bathys = File.depth(2:end,2:end);
    case 'Rectangular'
end

clearvars File LakeTahoeGrid*
%% Read binary file from SI3D output
cd(Path_file)
fid = fopen(FileName3D);
dum1 = fread(fid,1,'int32');
n_frames = fread(fid,1,'int32');
dum2 = fread(fid,1,'int32');
dum3 = fread(fid,1,'int32');
ipoints = fread(fid,1,'int32');
dum4 = fread(fid,1,'int32');
istep = zeros(n_frames,1);
year1 = zeros(n_frames,1);
month1 = zeros(n_frames,1);
day1 = zeros(n_frames,1);
hour1 = zeros(n_frames,1);

h = zeros(ipoints,1);
u = zeros(ipoints,1);
v = zeros(ipoints,1);
w = zeros(ipoints,1);
Dv = zeros(ipoints,1);
T = zeros(ipoints,1);
q2 = zeros(ipoints,1);
kh = zeros(ipoints,1);
q2l = zeros(ipoints,1);

tic;
for count = 1:n_frames+1
    dum5 = fread(fid,1,'int32');
    st = feof(fid);
    if (st == 0)
        istep(count)=fread(fid,1,'int32');
        year1(count)=fread(fid,1,'int32');
        month1(count)=fread(fid,1,'int32');
        day1(count)=fread(fid,1,'int32');
        hour1(count)=fread(fid,1,'float32');
        % ... Read all data for present time slice
        if (count == 1)
            out_array=fread(fid,13*ipoints,'float32');
            x = out_array(1:13:length(out_array)-12);
            y = out_array(2:13:length(out_array)-11);
            z = out_array(3:13:length(out_array)-10);
            h = out_array(4:13:length(out_array)-9);
            u = out_array(5:13:length(out_array)-8);
            v = out_array(6:13:length(out_array)-7);
            w = out_array(7:13:length(out_array)-6);
            Dv = out_array(8:13:length(out_array)-5);
            T = out_array(9:13:length(out_array)-4);
            q2 = out_array(10:13:length(out_array)-3);
            q2l = out_array(11:13:length(out_array)-2);
            kh = out_array(12:13:length(out_array)-1);
            Av = out_array(13:13:length(out_array));

            %% Code generator
            % To select the unique values of the array results that compound
            % the matrix of the numerical solution
            xmin = min(x);
            xmax = max(x);
            ymin = min(y);
            ymax = max(y);
            zp = unique(z);
%             zp(end+1) = zp(end) + 1;
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
                    zp = -zp*ddz;
            end
            xp = xg(1,:,1)'*dx;
            yp = yg(:,1,1)*dx;
            % To create the 3D domain with the real dimensions of the
            % structured grid
            [xgf,ygf,zgf] = meshgrid(xp,yp,zp);
            [m1,m2,m3] = size(xgf);
            clear xp yp zp xg yg zg

            % Creation of bathymetry
            switch BathyFlag
                case 'Lake'
                    Bathy = NaN(m1,m2,m3);
                    Bathy(:,:,1) = Bathys;
                    Bathy = Bathy(:);
                case 'Rectangular'
            end

            tic

            switch saveMatlabFlag
                case 'YES'
                    data.x = xgf;
                    data.y = ygf;
                    data.z = zgf;
                case 'NO'
            end

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
            sum1 = length(icodev);
            sum2 = length(icodex);
            clearvars codecellx codecellv codev codex
            disp(['Time needed to generate the code for indexing data is ',num2str(toc),' seconds'])

            if sum1 ~= length(x) || sum2 ~= length(x)
                error('The code generator is not working properly. The lenght of the sum of idata must be the same length as the vectors with the solution of SI3D')
            end
        else
            out_array=fread(fid,10*ipoints,'float32');
            h = out_array(1:10:length(out_array)-9);
            u = out_array(2:10:length(out_array)-8);
            v = out_array(3:10:length(out_array)-7);
            w = out_array(4:10:length(out_array)-6);
            Dv = out_array(5:10:length(out_array)-5);
            T = out_array(6:10:length(out_array)-4);
            q2 = out_array(7:10:length(out_array)-3);
            q2l = out_array(8:10:length(out_array)-2);
            kh = out_array(9:10:length(out_array)-1);
            Av = out_array(10:10:length(out_array));
        end
        TKE = q2./2; % Turbulent Kinectic Energy
        ml = q2l./q2; % Mixing macroscale

        dum = fread(fid,1,'int32');
        clear out_array;

        %% Paraview file creation .vtk
        hv = NaN(length(xgf(:)),1);
        uv = NaN(length(xgf(:)),1);
        vv = NaN(length(xgf(:)),1);
        wv = NaN(length(xgf(:)),1);
        Tv = NaN(length(xgf(:)),1);
        Dvv = NaN(length(xgf(:)),1);
        TKEv = NaN(length(xgf(:)),1);
        mlv = NaN(length(xgf(:)),1);
        khv = NaN(length(xgf(:)),1);
        Avv = NaN(length(xgf(:)),1);

        hv(icodev) = h(icodex);
        uv(icodev) = u(icodex);
        vv(icodev) = v(icodex);
        wv(icodev) = w(icodex);
        Tv(icodev) = T(icodex);
        Dvv(icodev) = Dv(icodex)./(100^2);
        TKEv(icodev) = TKE(icodex);
        mlv(icodev) = ml(icodex);
        khv(icodev) = kh(icodex);
        Avv(icodev) = kh(icodex);

        switch saveParaviewFlag
            case 'YES'
                % Start date in date time
                tic1 = tic;
                cd(Path_save_Paraview)
                switch vorticityFlag
                    case 'YES'
                        u3D = reshape(uv,m1,m2,m3);
                        v3D = reshape(vv,m1,m2,m3);
                        w3D = reshape(wv,m1,m2,m3);
                        [wu,wv,ww,cav] = curl(xgf,ygf,zgf,u3D,v3D,w3D);
%                         vtkwriteSV_v1([outputFile,'_',num2str((istep(count))*dt/3600),'.vtk'],...
%                             'structured_grid',xgf,ygf,zgf,'vectors','U(m/s)',uv,vv,wv,'vectors','Vorticity',wu,wv,ww,...
%                             'scalars','T(C)',Tv,'scalars','Dv(m2/s)',Dvv,'scalars','cav',cav,'scalars','TKE(m2/s2)',TKEv,...
%                             'scalars','ml(m)',mlv,'scalars','kh(m2/s2)',khv,'scalars','Av(m2/s)',Avv,...
%                             'precision','6','Binary')
                        vtkwriteSV_v1([outputFile,'_',num2str((istep(count))*dt/3600),'.vtk'],...
                            'structured_grid',xgf,ygf,zgf,'vectors','U(m/s)',uv,vv,wv,'vectors','Vorticity',wu,wv,ww,...
                            'scalars','T(C)',Tv,...
                            'precision','6','Binary')
                        disp(['Time frame ',num2str(count),' is ',num2str(toc(tic1)),' seconds'])
                    case 'NO'
                        switch BathyFlag
                            case 'Lake'
%                                 vtkwriteSV_v1([outputFile,'_',num2str((istep(count))*dt/3600),'.vtk']...
%                                     ,'structured_grid',xgf,ygf,zgf,'vectors','U(m/s)',uv,vv,wv,'scalars','T(C)',Tv,...
%                                     'scalars','Dv(m2/s)',Dvv,'scalars','TKE(m2/s2)',TKEv,'scalars','ml(m)',mlv,...
%                                     'scalars','kh(m2/s)',khv,'scalars','Av(m2/s)',Avv,'scalars','Depth(m)',Bathy,...
%                                     'precision','6','Binary')
                                vtkwriteSV_v1([outputFile,'_',num2str((istep(count))*dt/3600),'.vtk'],...
                                    'structured_grid',xgf,ygf,zgf,'vectors','U(m/s)',uv,vv,wv,...
                                    'scalars','T(C)',Tv,...
                                    'precision','6','binary')
                                disp(['Time frame ',num2str(count),' is ',num2str(toc(tic1)),' seconds'])
                            case 'Rectangular'
%                                 vtkwriteSV_v1([outputFile,'_',num2str((istep(count))*dt/3600),'.vtk']...
%                                     ,'structured_grid',xgf,ygf,zgf,'vectors','U(m/s)',uv,vv,wv,'scalars','T(C)',Tv,...
%                                     'scalars','Dv(m2/s)',Dvv,'scalars','TKE(m2/s2)',TKEv,'scalars','ml(m)',mlv,...
%                                     'scalars','kh(m2/s)',khv,'scalars','Av(m2/s)',Avv,...
%                                     'precision','6','Binary')
                                vtkwriteSV_v1([outputFile,'_',num2str((istep(count))*dt/3600),'.vtk'],...
                                    'structured_grid',xgf,ygf,zgf,'vectors','U(m/s)',uv,vv,wv,...
                                    'scalars','T(C)',Tv,...
                                    'precision','6','Binary')
                                disp(['Time frame ',num2str(count),' is ',num2str(toc(tic1)),' seconds'])
                        end
                end

                if count == n_frames + 1
                    fprintf(fid_PV,'%s\n',['		{ "name" : "',outputFile,'_',num2str((istep(count))*dt/3600),'.vtk",',' "time" : ',num2str((istep(count))*dt),'}']);
                else
                    fprintf(fid_PV,'%s\n',['		{ "name" : "',outputFile,'_',num2str((istep(count))*dt/3600),'.vtk",',' "time" : ',num2str((istep(count))*dt),'},']);
                end

            case 'NO'
                disp('Paraview file was not created')
            otherwise
                error('Did not specify the save paraview flag')
        end
        switch saveMatlabFlag
            case 'YES'
                tic
                data.h(:,count) = hv;
                data.u(:,count) = uv;
                data.v(:,count) = vv;
                data.w(:,count) = wv;
                data.T(:,count) = Tv;
                data.Dv(:,count) = Dvv;
                data.TKE(:,count) = TKE;
                data.ml(:,count) = ml;
                data.kh(:,count) = kh;
                data.Av(:,count) = Avv;
                disp(['Time frame ',num2str(count),' is ',num2str(toc),' seconds'])
            case 'NO'
        end
    else
        n_frames = count-1;
        break
    end

end
fclose(fid);

%% Creation of paraview reference file
switch saveParaviewFlag
    case 'YES'
        cd(Path_save_Paraview)
        StartDate = datetime(StartDate,'Format','yyy-MM-dd HH:mm:ss');
        Date = StartDate + istep.*dt/3600/24;
        DateStr = datestr(Date);
        Id = transpose(1:length(Date));
        Id = Id - 1;
        ParaviewRef = table(Id,istep.*dt,Date,DateStr);
        ParaviewRef.Properties.VariableNames = {'Id','Seconds','Date','DateStr'};
        writetable(ParaviewRef,'ParaviewRef.txt','Delimiter','tab')

        fprintf(fid_PV,'%s\n','	]');
        fprintf(fid_PV,'%s\n','}');
        fclose(fid_PV);
    case 'NO'
end

disp(['Time needed to read the binary file from SI3D and create paraview files is ',num2str(toc),' seconds'])

%% Saving Matlab File data
switch saveMatlabFlag
    case 'YES'
        tic
        eval([outputFile,'=data;'])

        cd(Path_save_struct)
        save(outputFile,outputFile,'-v7.3')

        disp(['Time needed to save the data is ',num2str(toc),' seconds'])
    case 'NO'
        disp('Matlab Structure not saved')
    otherwise
        error('Did not specify the save_flag')
end

disp(['Total time needed to create all the files requested is ',num2str(toc(ticT)/60),' minutes'])

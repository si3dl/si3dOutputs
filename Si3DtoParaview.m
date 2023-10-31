function n_frames = Si3DtoParaview(Pathfile,PathSave,StartDate,DeltaZ,dx,dz,dt,iTurb,itspf,nTracer,concTr)
% -------------------------------------------------------------------------
% This script uses the binary files generated from the SI3D simulations for
% 3D outputs and saves the data as vtk files to visualize the results in
% ParaView. It was used previous codes created by Alicia Cortes
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
% Version 5: The code was modified to include the cases from SI3D when
% turbulence parameters are either saved or not. It also includes the new
% formulation to save 3D files and surface planes from a time step
% different from 0.

% NOTES:

% 1. vtkwrite function is needed IF Paraview output is desired. This script
% will run the section where .vtk files are created

% 2. The paraview files will be named with the format 'si3D_xxx.vtk', where
% the xxx will be numerical values representing the hours after the
% simulation started.

% 3. It is imperative that Si3D was run using the same values for iht and
% ipxml. It is also imperative that the run was done using the same itspfh
% and itspf parameters from the input file.

% 5. It is recommended to run the code from a text editor for less memory
% consumption.

% Author: Sergio Valbuena
% Date: 02-12-2022

% ----------------------- USER SECTION START-------------------------------
FileName3D = 'ptrack_hydro.bnr';
PlaneName = 'plane_2';
outputFile = 'si3D';
FileNameZ = 'si3d_layer.txt';   % [m] Only needed if DeltaZ is variable
fileTracer = 'tracer_';
% ---------------------- USER SECTION END --------------------------------
% --------------------  CODE SECTION START -------------------------------
%% Definition of depths and horizontal dimensions
switch DeltaZ
    case 'variable'
        cd(Pathfile)
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

fread(fid3D,1,'int32');
fread(fidPL,1,'int32');
n_frames3D = fread(fid3D,1,'int32');
n_framesPL = fread(fidPL,1,'int32');

fread(fidPL,1,'int32');
fread(fidPL,1,'int32');
ipointsPL = fread(fidPL,1,'int32');
fread(fidPL,1,'int32');

fread(fid3D,1,'int32');
fread(fid3D,1,'int32');
ipoints3D = fread(fid3D,1,'int32');
fread(fid3D,1,'int32');

if nTracer > 0
    fidTr = NaN(nTracer,1);
    n_framesTr = NaN(nTracer,1);
    ipointsTr = NaN(nTracer,1);
    for tr = 1:nTracer
        fileNameTr = [fileTracer,num2str(tr)];
        fidTr(tr) = fopen(fileNameTr);
        fread(fidTr(tr),1,'int32');
        n_framesTr(tr) = fread(fidTr(tr),1,'int32');
        if n_frames3D ~= n_framesPL || n_frames3D ~= n_framesTr(tr)
            error('The number of frames between the surface plane and 3D file are not the same')
        elseif n_frames3D == n_framesPL && n_frames3D == n_framesTr(tr)
            n_frames = n_frames3D;
        end
        fread(fidTr(tr),1,'int32');
        fread(fidTr(tr),1,'int32');
        ipointsTr(tr) = fread(fidTr(tr),1,'int32');
        fread(fidTr(tr),1,'int32');
    end
else
    if n_frames3D ~= n_framesPL
        error('The number of frames between the surface plane and 3D file are not the same')
    elseif n_frames3D == n_framesPL
        n_frames = n_frames3D;
    end
end

if itspf == 0
    n_frames = n_frames + 1;
elseif itspf ~= 0
    n_frames = n_frames + 2;
end

istep = zeros(n_frames,1);
year1 = zeros(n_frames,1);
month1 = zeros(n_frames,1);
day1 = zeros(n_frames,1);
hour1 = zeros(n_frames,1);

ticT = tic;
for count = 1:n_frames
    fread(fid3D,1,'int32');
    fread(fidPL,1,'int32');
    st3d = feof(fid3D);
    stpl = feof(fidPL);
    sttr = 0;
    if nTracer > 0
        for tr = 1:nTracer
            fread(fidTr(tr),1,'int32');
            sttr = feof(fidTr(tr));
        end
    end
    if (st3d == 0) || (stpl == 0) || (sttr == 0) 
        istep(count)=fread(fid3D,1,'int32');
        year1(count)=fread(fid3D,1,'int32');
        month1(count)=fread(fid3D,1,'int32');
        day1(count)=fread(fid3D,1,'int32');
        hour1(count)=fread(fid3D,1,'float32');
        fread(fidPL,1,'int32');
        fread(fidPL,1,'int32');
        fread(fidPL,1,'int32');
        fread(fidPL,1,'int32');
        fread(fidPL,1,'float32');
        if nTracer > 0
            for tr = 1:nTracer
                fread(fidTr(tr),1,'int32');
                fread(fidTr(tr),1,'int32');
                fread(fidTr(tr),1,'int32');
                fread(fidTr(tr),1,'int32');
                fread(fidTr(tr),1,'int32');
            end
        end
        % ... Read all data for present time slice
        if (count == 1)
            if (iTurb == 1)
                out_array3D=fread(fid3D,13*ipoints3D,'float32');
                x = out_array3D(1:13:length(out_array3D)-12);
                y = out_array3D(2:13:length(out_array3D)-11);
                z = out_array3D(3:13:length(out_array3D)-10);
                %h = out_array3D(4:13:length(out_array3D)-9);
                u = out_array3D(5:13:length(out_array3D)-8);
                v = out_array3D(6:13:length(out_array3D)-7);
                w = out_array3D(7:13:length(out_array3D)-6);
                Dv = out_array3D(8:13:length(out_array3D)-5);
                T = out_array3D(9:13:length(out_array3D)-4);
                q2 = out_array3D(10:13:length(out_array3D)-3);
                q2l = out_array3D(11:13:length(out_array3D)-2);
                kh = out_array3D(12:13:length(out_array3D)-1);
                Av = out_array3D(13:13:length(out_array3D));

                out_arrayPL = fread(fidPL,8*ipointsPL,'float32');
                dumPLx = out_arrayPL(1:8:length(out_arrayPL)-7);
                dumPLy = out_arrayPL(2:8:length(out_arrayPL)-6);
                s = out_arrayPL(8:8:length(out_arrayPL));

                if nTracer > 0
                    tracerc = NaN(ipointsTr,nTracer);
                    for i = tr:nTracer
                        out_arrayTr=fread(fidTr(tr),5*ipoints,'float32');
                        dumTrx = out_arrayTr(1:5:length(out_arrayTr)-4);
                        dumTry = out_arrayTr(2:5:length(out_arrayTr)-3);
                        dumTrz = out_arrayTr(3:5:length(out_arrayTr)-2);
                        %tracerm(:,count)=out_array(4:5:length(out_array)-1);
                        tracerc(:,tr)=out_arrayTr(5:5:length(out_arrayTr));
                    end
                end

            elseif (iTurb == 0)
                out_array3D=fread(fid3D,8*ipoints3D,'float32');
                x = out_array3D(1:8:length(out_array3D)-7);
                y = out_array3D(2:8:length(out_array3D)-6);
                z = out_array3D(3:8:length(out_array3D)-5);
                %h = out_array3D(4:9:length(out_array3D)-4);
                u = out_array3D(5:8:length(out_array3D)-3);
                v = out_array3D(6:8:length(out_array3D)-2);
                w = out_array3D(7:8:length(out_array3D)-1);
                T = out_array3D(8:8:length(out_array3D));

                out_arrayPL = fread(fidPL,8*ipointsPL,'float32');
                dumPLx = out_arrayPL(1:8:length(out_arrayPL)-7);
                dumPLy = out_arrayPL(2:8:length(out_arrayPL)-6);
                s = out_arrayPL(8:8:length(out_arrayPL));

                if nTracer > 0
                    tracerc = NaN(ipoints3D,nTracer);
                    for tr = 1:nTracer
                        out_arrayTr=fread(fidTr(tr),5*ipointsTr(tr),'float32');
                        dumTrx = out_arrayTr(1:5:length(out_arrayTr)-4);
                        dumTry = out_arrayTr(2:5:length(out_arrayTr)-3);
                        dumTrz = out_arrayTr(3:5:length(out_arrayTr)-2);
                        %tracerm(:,count)=out_array(4:5:length(out_array)-1);
                        tracerc(:,tr)=out_arrayTr(5:5:length(out_arrayTr));
                    end
                end
            end

           %% Code generator
            % To select the unique values of the array results that compound
            % the matrix of the numerical solution
            xmin = min(x);
            xmax = max(x);
            ymin = min(y);
            ymax = max(y);
            zp = unique(z);

            %% Adding results of Grad(u,v) = 0
            idata = y == ymin;
            filtx = x(idata);
            filtz = z(idata);
            idatasum = sum(idata);
            ynew = [y;(ymin-1)*ones(idatasum,1)];
            znew = [z;filtz];
            xnew = [x;filtx];
            filtu = u(idata);
            filtv = v(idata);
            filtw = w(idata);
            filtT = T(idata);
            u = [u;filtu];
            v = [v;filtv];
            w = [w;filtw];
            T = [T;filtT];

            if nTracer > 0
                for tr = 1:nTracer
                    filtTr = tracerc(idata,tr);
                    tracerc1(:,tr) = [tracerc(:,tr);filtTr];
                end
                tracerc = tracerc1;
                clearvars tracerc1
            end

            if iTurb
                filtDv = Dv(idata);
                filtq2 = q2(idata);
                filtq2l = q2l(idata);
                filtkh = kh(idata);
                filtAv = Av(idata);
                Dv = [Dv;filtDv];
                Av = [Av;filtAv];
                q2 = [q2;filtq2];
                q2l = [q2l;filtq2l];
                kh = [kh;filtkh];
            end

            idata = xnew == xmin;
            filty = ynew(idata);
            filtz = znew(idata);
            idatasum = sum(idata);
            xnew = [xnew;(xmin-1)*ones(idatasum,1)];
            znew = [znew;filtz];
            ynew = [ynew;filty];
            filtu = u(idata);
            filtv = v(idata);
            filtw = w(idata);
            filtT = T(idata);
            u = [u;filtu];
            v = [v;filtv];
            w = [w;filtw];
            T = [T;filtT];

            if nTracer > 0
                for tr = 1:nTracer
                    filtTr = tracerc(idata,tr);
                    tracerc1(:,tr) = [tracerc(:,tr);filtTr];
                end
                tracerc = tracerc1;
                clearvars tracerc1
            end

            if iTurb
                filtDv = Dv(idata);
                filtq2 = q2(idata);
                filtq2l = q2l(idata);
                filtkh = kh(idata);
                filtAv = Av(idata);
                Dv = [Dv;filtDv];
                Av = [Av;filtAv];
                q2 = [q2;filtq2];
                q2l = [q2l;filtq2l];
                kh = [kh;filtkh];
            end
            iz = znew == 2;
            filtx = xnew(iz);
            filty = ynew(iz);
            filtz = znew(iz);
            filtu = u(iz);
            filtv = v(iz);
            filtw = w(iz);
            filtT = T(iz);
            xnew = [xnew;filtx];
            ynew = [ynew;filty];
            znew = [znew;ones(length(filtz),1)];
            u = [u;filtu];
            v = [v;filtv];
            w = [w;filtw];
            T = [T;filtT];

            if nTracer > 0
                for tr = 1:nTracer
                    filtTr = tracerc(iz,tr);
                    tracerc1(:,tr) = [tracerc(:,tr);filtTr];
                end
                tracerc = tracerc1;
                clearvars tracerc1
            end

            %%
            xmin = min(xnew);
            xmax = max(xnew);
            ymin = min(ynew);
            ymax = max(ynew);
            zp = unique(znew);   
            
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
                        Layer = Layer - 1;
                        idata = ~ismember(Layer,zp);
                        Layer(idata) = [];
                        Depth(idata) = [];
                    end
                    zp = -Depth;
                case 'constant'
                    zp = -(zp-1)*ddz;
%                     zp = [zp;zp(end)-ddz/2];
%                     zp = sort(zp,'descend');
            end
            xp = (xg(1,:,1)-1)'*dx;
            yp = (yg(:,1,1)-1)*dx;

            % To create the 3D domain with the real dimensions of the
            % structured grid
            [xgf,ygf,zgf] = meshgrid(xp,yp,zp);

            clear xp yp zp xg yg zg
            
            % To create code for the surface plane
            coord_simPL = [dumPLx,dumPLy];
            codePL = num2str(coord_simPL);
            codecellPL = cellstr(codePL);
            clearvars codePL
            codePL = regexprep(codecellPL,'\s+','A');
            codePL = strcat(codePL,'A1');
            
            tic;
            % To create a code that idenfies the indexes of the numerical
            % simulations within the whole domain.
            % The numerical values are converted into strings for code
            % identification.
            coord_sim = [xnew,ynew,znew];
            coord_dom = [xv,yv,zv];
            codex = num2str(coord_sim);
            codev = num2str(coord_dom);
            
            % clearvars xv yv zv
            codecellx = cellstr(codex);
            codecellv = cellstr(codev);
            clearvars codex codev
            codex = regexprep(codecellx,'\s+','A');
            codev = regexprep(codecellv,'\s+','A');
            
            [~,icodev,icodex] = intersect(codev,codex,'sorted');
            [~,icodevPL,icodexPL] = intersect(codev,codePL,'sorted');
            sum1 = length(icodev);
            sum2 = length(icodex);
            clearvars codecellx codecellv codev codex
            disp(['Time needed to generate the code for indexing data is ',num2str(toc),' seconds'])
            
            if sum1 ~= length(xnew) || sum2 ~= length(xnew)
                error('The code generator is not working properly. The lenght of the sum of idata must be the same length as the vectors with the solution of SI3D')
            end
        else
            if (iTurb == 1)
                out_array3D=fread(fid3D,10*ipoints3D,'float32');
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
                s = out_arrayPL(6:6:length(out_arrayPL));
                if nTracer > 0
                    for tr = 1:nTracer
                    out_arrayTr=fread(fidTr(tr),2*ipoints3D,'float32');
                    %tracerm(:,count)=out_array(4:5:length(out_array)-1);
                    tracerc(:,tr)=out_arrayTr(2:2:length(out_arrayTr));
                    end
                end

            elseif (iTurb == 0)
                out_array3D=fread(fid3D,5*ipoints3D,'float32');
                u = out_array3D(2:5:length(out_array3D)-3);
                v = out_array3D(3:5:length(out_array3D)-2);
                w = out_array3D(4:5:length(out_array3D)-1);
                T = out_array3D(5:5:length(out_array3D));                
                out_arrayPL=fread(fidPL,6*ipointsPL,'float32');
                s = out_arrayPL(6:6:length(out_arrayPL));
                if nTracer > 0
                    for tr = 1:nTracer
                    out_arrayTr=fread(fidTr(tr),2*ipointsTr(tr),'float32');
                    %tracerm(:,count)=out_array(1:2:length(out_array)-1);
                    tracerc(:,tr)=out_arrayTr(2:2:length(out_arrayTr));
                    end
                end
            end

            %% Adding results of Grad(u,v) = 0
            idata = y == min(y);
            filtx = x(idata);
            filtz = z(idata);
            idatasum = sum(idata);
            ynew = [y;(min(y)-1)*ones(idatasum,1)];
            znew = [z;filtz];
            xnew = [x;filtx];
            filtu = u(idata);
            filtv = v(idata);
            filtw = w(idata);
            filtT = T(idata);
            u = [u;filtu];
            v = [v;filtv];
            w = [w;filtw];
            T = [T;filtT];

            if nTracer > 0
                for tr = 1:nTracer
                    filtTr = tracerc(idata,tr);
                    tracerc1(:,tr) = [tracerc(:,tr);filtTr];
                end
                tracerc = tracerc1;
                clearvars tracerc1
            end

            if iTurb
                filtDv = Dv(idata);
                filtq2 = q2(idata);
                filtq2l = q2l(idata);
                filtkh = kh(idata);
                filtAv = Av(idata);
                Dv = [Dv;filtDv];
                Av = [Av;filtAv];
                q2 = [q2;filtq2];
                q2l = [q2l;filtq2l];
                kh = [kh;filtkh];
            end

            idata = xnew == min(xnew);
            filty = ynew(idata);
            filtz = znew(idata);
            idatasum = sum(idata);
            xnew = [xnew;(min(xnew)-1)*ones(idatasum,1)];
            znew = [znew;filtz];
            ynew = [ynew;filty];
            filtu = u(idata);
            filtv = v(idata);
            filtw = w(idata);
            filtT = T(idata);
            u = [u;filtu];
            v = [v;filtv];
            w = [w;filtw];
            T = [T;filtT];

            if nTracer > 0
                for tr = 1:nTracer
                    filtTr = tracerc(idata,tr);
                    tracerc1(:,tr) = [tracerc(:,tr);filtTr];
                end
                tracerc = tracerc1;
                clearvars tracerc1
            end

            if iTurb
                filtDv = Dv(idata);
                filtq2 = q2(idata);
                filtq2l = q2l(idata);
                filtkh = kh(idata);
                filtAv = Av(idata);
                Dv = [Dv;filtDv];
                Av = [Av;filtAv];
                q2 = [q2;filtq2];
                q2l = [q2l;filtq2l];
                kh = [kh;filtkh];
            end

            iz = znew == 2;
            filtx = xnew(iz);
            filty = ynew(iz);
            filtz = znew(iz);
            filtu = u(iz);
            filtv = v(iz);
            filtw = w(iz);
            filtT = T(iz);
            xnew = [xnew;filtx];
            ynew = [ynew;filty];
            znew = [znew;1*ones(length(filtz),1)];
            u = [u;filtu];
            v = [v;filtv];
            w = [w;filtw];
            T = [T;filtT];
            if nTracer > 0
                for tr = 1:nTracer
                    filtTr = tracerc(iz,tr);
                    tracerc1(:,tr) = [tracerc(:,tr);filtTr];
                end
                tracerc = tracerc1;
                clearvars tracerc1;
            end

        end
        
        if (iTurb == 1)
            TKE = q2./2; % Turbulent Kinectic Energy
            ml = q2l./q2; % Mixing macroscale
        end

        fread(fid3D,1,'int32');
        fread(fidPL,1,'int32');
        if nTracer > 0
            for tr = 1:nTracer
                fread(fidTr(tr),1,'int32');
            end
        end
        clear out_arrayPL out_array3D out_arrayTr
        
        %% Paraview file creation .vtk
        lv = NaN(length(xgf(:)),1);
        sv = zeros(length(xgf(:)),1);
        uv = NaN(length(xgf(:)),1);
        vv = NaN(length(xgf(:)),1);
        wv = NaN(length(xgf(:)),1);
        Tv = NaN(length(xgf(:)),1);
           
        lv(icodev) = sv(icodex);
        lv(icodevPL) = s(icodexPL);
        uv(icodev) = u(icodex);
        vv(icodev) = v(icodex);
        wv(icodev) = w(icodex);
        Tv(icodev) = T(icodex);

        if nTracer > 0 
            conc = NaN(length(xgf(:)),nTracer);
            for tr = 1:nTracer
                conc(icodev,tr) = tracerc(icodex,tr);
            end
        end
    
        if iTurb == 1
            Dvv = NaN(length(xgf(:)),1);
            TKEv = NaN(length(xgf(:)),1);
            mlv = NaN(length(xgf(:)),1);
            khv = NaN(length(xgf(:)),1);
            Avv = NaN(length(xgf(:)),1);
        
            Dvv(icodev) = Dv(icodex)./(100^2);
            TKEv(icodev) = TKE(icodex);
            mlv(icodev) = ml(icodex);
            khv(icodev) = kh(icodex);
            Avv(icodev) = Av(icodex);
        end
        
        tic1 = tic;
        cd(PathSave)

        if iTurb == 1
            scalarType = {   'T(C)',   'l(m)','Dv(m2/s)','TKE(m2/s2)',  'ml(m)','kh(m2/s)','Av(m2/s)'};
            varInp =     {       Tv,       lv,       Dvv,        TKEv,      mlv,       khv,       Avv};
            if nTracer > 0
                for tr = 1:nTracer
                    scalarType(end+1) = concTr(tr);
                    varInp(end+1) = conc(:,tr);
                end
            end
            vtkwriteSV_v1([outputFile,'_',num2str((istep(count))*dt/3600),'.vtk'],...
                'structured_grid',xgf,ygf,zgf,'vectors','U(m/s)',uv,vv,wv,...
                'scalars',scalarType,varInp,'precision','6','binary')
            disp(['Time frame ',num2str(count),' is ',num2str(toc(tic1)),' seconds'])
        elseif iTurb == 0
            scalarType = {   'T(C)',   'l(m)'};
            varInp =     {       Tv,       lv};
            if nTracer > 0
                for tr = 1:nTracer
                    
                    scalarType{end+1} = concTr(tr);
                    varInp{end+1} = conc(:,tr);
                end
            end
            vtkwriteSV_v1([outputFile,'_',num2str((istep(count))*dt/3600),'.vtk'],...
                'structured_grid',xgf,ygf,zgf,'vectors','U(m/s)',uv,vv,wv,...
                'scalars',scalarType,varInp,'precision','6','binary')
            disp(['Time frame ',num2str(count),' is ',num2str(toc(tic1)),' seconds'])
        end
        
        if count == n_frames
            fprintf(fidPV,'%s\n',['		{ "name" : "',outputFile,'_',num2str((istep(count))*dt/3600),'.vtk",',' "time" : ',num2str((istep(count))*dt),'}']);
        else
            fprintf(fidPV,'%s\n',['		{ "name" : "',outputFile,'_',num2str((istep(count))*dt/3600),'.vtk",',' "time" : ',num2str((istep(count))*dt),'},']);
        end
        clearvars tracerc conc
    else
        n_frames = count-1;
        break
    end
    
end
fclose(fid3D);
fclose(fidPL);
if nTracer > 0
    for tr = 1:nTracer
        fclose(fidTr(tr));
    end
end

%%
fprintf(fidPV,'%s\n','	]');
fprintf(fidPV,'%s\n','}');
fclose(fidPV);
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

disp(['Time needed to read the binary file from SI3D and create paraview files is ',num2str(toc),' seconds'])

disp(['Total time needed to create all the files requested is ',num2str(toc(ticT)/60),' minutes'])
end
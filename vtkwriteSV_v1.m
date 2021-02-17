function vtkwriteSV_v1( filename,dataType,varargin)
% vtkwriteSV_v1 Writes 3D Matlab array into VTK file format for paraview
% visualizations

% Way to use vtkwriteSV_v1:
% vtkwrite(filename,'structured_grid',x,y,z,'vectors',title,u,v,w,'scalars',
% title,T). This usage of the function will write a 3D structured grid that
% contains vectors and scalars values. x,y,z,u,v,w,r have to be the same
% size and contain the corresponding positon.

% 'structured grid' is an option to select this type of grid. 'Unstructure
% Grid' is considered to be included in a further version of this function.
% As of today, vtkwriteSV_v1 only works for structured grids.

% 'x','y','z' must be a 3x3 matrix with the correct sizes. This is to make
% sure that the domain is a complete cube and contains all points
% acordingly referenced.

% Each variable to write must be specified as scalars or vectors followed
% by the name of the variable and the components. 3 in the case of a vector
% and 1 in the case of a scalar. 

% The file can be saved using ASCII or Binary formats. Both will be
% recognized by Paraview. Having Binary files will decrease the size of
% these.

% NOTE: vtk files are for a single time step/frame. Therefore, care must be
% used when using this function for a unsteady transient simulation. For
% loop to create all time steps can be used. 

% Version 1.0
% Copyright, Sergio Valbuena, 2020
% Codes are modified from William Thielicke and David Gingras's submission.

fid = fopen(filename, 'w');
% VTK files contain five major parts
% 1. VTK DataFile Version
fprintf(fid, '# vtk DataFile Version 2.0\n');
% 2. Title
fprintf(fid, 'VTK from Matlab\n');

binaryflag = any(strcmpi(varargin, 'BINARY'));
if any(strcmpi(varargin, 'PRECISION'))
    precision = num2str(varargin{find(strcmpi(varargin, 'PRECISION'))+1});
else
    precision = '4';
end

switch upper(dataType)
    case {'STRUCTURED_GRID','UNSTRUCTURED_GRID'}
        % 3. The format of the data. As said before it can be ASCII or Binary
        % fprintf is used for ASCII and fwrite for binary.
        if numel(varargin)<6, error('Not enough input arguments'); end
        setdataformat(fid, binaryflag);
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        if sum(size(x)==size(y) & size(y)==size(z))~=length(size(x))
            error('Input dimesions do not match')
        end
        n_elements = numel(x);
        
        % 4. Type of Dataset Developed to consider unstructed and
        % structured grids. V1 only works for structured grids.
        if strcmpi(dataType,'STRUCTURED_GRID')
            fprintf(fid, 'DATASET STRUCTURED_GRID\n');
            fprintf(fid, 'DIMENSIONS %d %d %d\n', size(x,1), size(x,2), size(x,3));
        else
            fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
        end
        fprintf(fid, ['POINTS ' num2str(n_elements) ' float\n']);
        output = [x(:)'; y(:)'; z(:)'];
        
        if ~binaryflag
            spec = ['%0.', precision, 'f '];
            fprintf(fid, spec, output);
        else
            fwrite(fid, output, 'float', 'b');
        end
        % 5.This final part describe the dataset attributes and begins with the
        % keywords 'POINT_DATA' or 'CELL_DATA', followed by an integer number
        % specifying the number of points or cells. Other keyword/data combination
        % then define the actual dataset attribute values.
        fprintf(fid, ['\nPOINT_DATA ' num2str(n_elements)]);
        % Parse remaining argument.
        vidx = find(strcmpi(varargin,'VECTORS'));
        sidx = find(strcmpi(varargin,'SCALARS'));
        if vidx~=0
            for ii = 1:length(vidx)
                title = varargin{vidx(ii)+1};
                % Data enteries begin with a keyword specifying data type
                % and numeric format.
                fprintf(fid, ['\nVECTORS ', title,' float\n']);
                output = [varargin{ vidx(ii) + 2 }(:)';...
                    varargin{ vidx(ii) + 3 }(:)';...
                    varargin{ vidx(ii) + 4 }(:)'];
                if ~binaryflag
                    spec = ['%0.', precision, 'f '];
                    fprintf(fid, spec, output);
                else
                    fwrite(fid, output, 'float', 'b');
                end
            end
        end
        if sidx~=0
            for ii = 1:length(sidx)
                title = varargin{sidx(ii)+1};
                fprintf(fid, ['\nSCALARS ', title,' float\n']);
                fprintf(fid, 'LOOKUP_TABLE default\n');
                if ~binaryflag
                    spec = ['%0.', precision, 'f '];
                    fprintf(fid, spec, varargin{ sidx(ii) + 2});
                else
                    output = reshape(varargin{sidx(ii)+2},1,n_elements);
                    fwrite(fid, output,'float','b');
                end
            end
        end
end
fclose(fid);
end

% Function of Binary format
function setdataformat(fid, binaryflag)
if ~binaryflag
    fprintf(fid, 'ASCII\n');
else
    a = fprintf(fid, 'BINARY\n');
end
end

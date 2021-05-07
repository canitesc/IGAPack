% MSH_TO_VTU: Export an unstructured mesh to VTK format for plotting.
%
%  msh_to_vtu (pts, sigma, disp, npts, filename, fieldname)
%
% INPUT:
%
%     pts:       points at which the field was computed
%     sigma:     stresses of the field at the selected points
%     disp:      displacements of the field at the selected points
%     npts:      number of points on each element, and for each parametric direction
%     filename:  name of the output file
%     fieldname: how to name the saved variable in the vtk file
%
% OUTPUT:
%
%    a vtk unstructured mesh file named <filename> is produced
%
% Copyright (C) 2009, 2010 Carlo de Falco, Rafael Vazquez
% Copyright (C) 2011, 2012 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Modified by VP Nguyen, Cardiff University, UK.

function msh_to_vtu (pts, sigma, disp, errorU,errorStress,npts, filename)

rdim = size (pts, 1);
ndim = numel (npts);

str1 = cat (2,'<?xml version="1.0"?> \n', ...
    '<VTKFile type="UnstructuredGrid" version="0.1"> \n', ...
    '<UnstructuredGrid> \n', ...
    '<Piece NumberOfPoints="%d" NumberOfCells="%d"> \n', ...
    '<PointData %s="%s">\n', ...
    '<DataArray type="Float32" Name="%s" format="ascii" NumberOfComponents="%d"> \n');

str2 = cat (2,...
    '</PointData> \n', ...
    '<Points> \n', ...
    '<DataArray type="Float32" NumberOfComponents="3"> \n');

str3 = cat (2, '\n', ...
    '</DataArray> \n', ...
    '</Points> \n', ...
    '<Cells> \n', ...
    '<DataArray type="Int32" Name="connectivity">\n');

str4 = cat (2, '\n </DataArray> \n', ...
    '<DataArray type="Int32" Name="offsets"> \n');

str5 = cat (2, '\n', ...
    '</DataArray> \n', ...
    '<DataArray type="Int32" Name="types"> \n');

str6 = cat (2, '\n', ...
    '</DataArray>\n', ...
    '</Cells> \n', ...
    '</Piece> \n', ...
    '</UnstructuredGrid> \n', ...
    '</VTKFile> \n');

% Check whether we are saving scalars or vectors.
% Even for 2D data, everything is saved in 3D
if (numel (size (sigma)) == 2)
    fieldclass = 'Scalars';
    ncomp = 1;
else
    fieldclass = 'Vectors';
    ncomp = 3;
end

if (rdim == 2)
    pts(3,:,:) = 0;
    %     if (ncomp == 3)
    %       values(3,:,:) = 0;
    %     end
end

% Number of points and cells
nel     = size (pts, 3);
npoints = prod (npts) * nel;
ncells  = prod (npts - 1) * nel;

if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.vtu'))
    filename = cat (2, filename, '.vtu');
end

fid = fopen (filename, 'w');
if (fid < 0)
    error ('msh_to_vtu: could not open file %s', filename);
end

% Connectivities for the local cells in one Bezier element
if (ndim == 2)
    xel = repmat (1:npts(1)-1, 1, npts(2)-1)';
    yel = repmat ((1:npts(2)-1), npts(1)-1, 1);
    yel = yel(:);
    local_conn = [sub2ind(npts, xel, yel), ...
        sub2ind(npts, xel+1, yel), ...
        sub2ind(npts, xel+1, yel+1), ...
        sub2ind(npts, xel, yel+1)]';
else
    error ('The unstructured plot is not ready for 3D geometries')
end

fprintf (fid, str1, ...
    npoints, ncells, fieldclass, 'fields', 'stress', ncomp);

fprintf (fid, '%g ', sigma(:));

fprintf (fid, '</DataArray> \n');

fprintf (fid, '<DataArray  type="Float32"  Name="ErrorStress" NumberOfComponents="3" format="ascii"> \n');
fprintf (fid, '%g ', errorStress(:));

fprintf (fid, '</DataArray> \n');


fprintf (fid, '<DataArray  type="Float64"  Name="u" NumberOfComponents="3" format="ascii"> \n');
fprintf (fid, '%g  ', disp(:));
fprintf (fid, '</DataArray> \n');

fprintf (fid, '<DataArray  type="Float64"  Name="ErrorU" NumberOfComponents="2" format="ascii"> \n');
fprintf (fid, '%g  ', errorU(:));
fprintf (fid, '</DataArray> \n');

fprintf (fid, str2);
fprintf (fid, '%g ', pts(:));

fprintf (fid, str3);
for iel = 1:nel
    fprintf (fid, '%d ', local_conn + (iel-1)*prod(npts) - 1);
end

fprintf (fid, str4);
fprintf (fid, '%d ', 4:4:ncells*4);
fprintf (fid, str5);
if (ndim == 2)
    fprintf (fid, '%d ', 9 * ones (1,ncells));
else
    fprintf (fid, '%d ', 12 * ones (1, ncells));
end
fprintf (fid, str6);

fclose (fid);

end

function nrbexport (varargin)

%
% NRBEXPORT: export NURBS geometries to a format compatible with the one used in GeoPDEs (version 0.6).
% 
% Calling Sequence:
% 
%   nrbexport (nurbs, filename);
%   nrbexport (nurbs, interfaces, boundaries, filename);
% 
% INPUT:
% 
%   nurbs    :  NURBS curve, surface or volume, see nrbmak.
%   interfaces: interface information for GeoPDEs (see nrbmultipatch)
%   boundaries: boundary information for GeoPDEs (see nrbmultipatch)
%   filename :  name of the output file.
% 
% 
% Description:
% 
%   The data of the nurbs structure is written in the file, in a format 
%     that can be read by GeoPDEs.
%
%    Copyright (C) 2011, 2014 Rafael Vazquez
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

if (nargin == 2)
  nurbs = varargin{1};
  filename = varargin{2};
  if (numel (nurbs) > 1)
    warning ('Automatically creating the interface information with nrbmultipatch')
    [interfaces, boundaries] = nrbmultipatch (nurbs);
  else
    interfaces = []; boundaries = [];
  end
elseif (nargin == 4)
  nurbs = varargin{1};
  interfaces = varargin{2};
  boundaries = varargin{3};
  filename = varargin{4};
else
  error ('nrbexport: wrong number of input arguments') 
end

fid = fopen (filename, 'w');
if (fid < 0)
  error ('nrbexport: cannot open file %s', filename);
end

ndim = numel (nurbs(1).order);
npatch = numel (nurbs);
fprintf (fid, '%s\n', '# nurbs mesh v.0.7');
fprintf (fid, '%s\n', '#');
fprintf (fid, '%s\n', ['# ' date]);
fprintf (fid, '%s\n', '#');

fprintf (fid, '%4i', ndim, npatch, numel(interfaces), 1);
fprintf (fid, '\n');
for iptc = 1:npatch
  fprintf (fid, '%s %i \n', 'PATCH', iptc);
  fprintf (fid, '%4i', nurbs(iptc).order-1);
  fprintf (fid, '\n');
  fprintf (fid, '%4i', nurbs(iptc).number);
  fprintf (fid, '\n');
  for ii = 1:ndim
    fprintf (fid, '%1.7f   ', nurbs(iptc).knots{ii});
    fprintf (fid, '\n');
  end

  if (ndim == 2)
    for ii = 1:ndim
      fprintf (fid, '%1.15f   ', nurbs(iptc).coefs(ii,:,:));
      fprintf (fid, '\n');
    end
    fprintf (fid, '%1.15f   ', nurbs(iptc).coefs(4,:,:));
    fprintf (fid, '\n');
  elseif (ndim == 3)
    for ii = 1:ndim
      fprintf (fid, '%1.15f   ', nurbs(iptc).coefs(ii,:,:,:));
      fprintf (fid, '\n');
    end
    fprintf (fid, '%1.15f   ', nurbs(iptc).coefs(4,:,:,:));
    fprintf (fid, '\n');
  end
end

for intrfc = 1:numel(interfaces)
  if (isfield (interfaces, 'ref'))
    fprintf (fid, '%s \n', interfaces(intrfc).ref);
  else
    fprintf (fid, '%s %i \n', 'INTERFACE', intrfc);
  end
  fprintf (fid, '%i %i \n', interfaces(intrfc).patch1, interfaces(intrfc).side1);
  fprintf (fid, '%i %i \n', interfaces(intrfc).patch2, interfaces(intrfc).side2);
  if (ndim == 2)
    fprintf (fid, '%i \n', interfaces(intrfc).ornt);
  elseif (ndim == 3)
    fprintf (fid, '%i %i %i \n', interfaces(intrfc).flag, interfaces(intrfc).ornt1, interfaces(intrfc).ornt2);
  end
end

% The subdomain part should be fixed
fprintf (fid, '%s \n', 'SUBDOMAIN 1');
fprintf (fid, '%i ', 1:npatch);
fprintf (fid, '\n');
    

for ibnd = 1:numel (boundaries)
  if (isfield (boundaries, 'name'))
    fprintf (fid, '%s \n', boundaries(ibnd).name);
  else
    fprintf (fid, '%s %i \n', 'BOUNDARY', ibnd);
  end
  fprintf (fid, '%i \n', boundaries(ibnd).nsides);
  for ii = 1:boundaries(ibnd).nsides
    fprintf (fid, '%i %i \n', boundaries(ibnd).patches(ii), boundaries(ibnd).faces(ii));
  end
end

fclose (fid);
end

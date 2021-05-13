/* Copyright (C) 2009, 2020 Carlo de Falco, Rafael Vazquez

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <octave/oct.h>
#include "low_level_functions.h"

DEFUN_DLD(basisfunder, args, nargout,"\n\
 BASISFUNDER:  B-Spline Basis function derivatives\n\
\n\
 Calling Sequence:\n\
\n\
   ders = basisfunder (ii, pl, uu, k, nd)\n\
\n\
    INPUT:\n\
\n\
      ii  - knot span\n\
      pl  - degree of curve\n\
      uu  - parametric points\n\
      k   - knot vector\n\
      nd  - number of derivatives to compute\n\
\n\
    OUTPUT:\n\
\n\
      ders - ders(n, i, :) (i-1)-th derivative at n-th point\n\
\n\
   Adapted from Algorithm A2.3 from 'The NURBS BOOK' pg72. \n\
\n\
")
{
  octave_value_list retval;

  if (nargout != 1 || args.length () != 5)
    print_usage ();

  const NDArray i = args(0).array_value ();
  int pl = args(1).int_value ();
  const NDArray u = args(2).array_value ();
  const RowVector U = args(3).row_vector_value ();
  int nd = args(4).int_value ();

  if (i.numel () != u.numel ())
    print_usage ();
 
  NDArray dersv (dim_vector (i.numel (), nd+1, pl+1), 0.0);
  NDArray ders(dim_vector(nd+1, pl+1), 0.0);
  for ( octave_idx_type jj(0); jj < i.numel (); jj++)
    {
      basisfunder (int (i(jj)), pl, u(jj), U, nd, ders);

      for (octave_idx_type kk(0); kk < nd+1; kk++)
        for (octave_idx_type ll(0); ll < pl+1; ll++)
          {
            dersv(jj, kk, ll) = ders(kk, ll);
          }
    }
  retval(0) = dersv;

  return retval;
}

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

DEFUN_DLD(bspderiv, args, nargout,"\n\
 BSPDERIV:  B-Spline derivative\n\
\n\
\n\
 Calling Sequence:\n\
\n\
          [dc,dk] = bspderiv(d,c,k)\n\
\n\
  INPUT:\n\
 \n\
    d - degree of the B-Spline\n\
    c - control points   double matrix(mc,nc)\n\
    k - knot sequence    double vector(nk)\n\
 \n\
  OUTPUT:\n\
 \n\
    dc - control points of the derivative     double  matrix(mc,nc)\n\
    dk - knot sequence of the derivative      double  vector(nk)\n\
 \n\
  Modified version of Algorithm A3.3 from 'The NURBS BOOK' pg98.\n\
")
{
  //if (bspderiv_bad_arguments(args, nargout)) 
  //  return octave_value_list(); 
  
  if ((nargout != 1 && nargout != 2) || args.length () != 3)
    print_usage ();

  int       d = args(0).int_value();
  const Matrix    c = args(1).matrix_value();
  const RowVector k = args(2).row_vector_value();
  octave_value_list retval;
  octave_idx_type mc = c.rows(), nc = c.cols(), nk = k.numel();
  Matrix dc (mc, nc-1, 0.0);
  RowVector dk(nk-2, 0.0);

  double tmp;
      
  for (octave_idx_type i(0); i<=nc-2; i++)
    {
      tmp = (double)d / (k(i+d+1) - k(i+1));
      for ( octave_idx_type j(0); j<=mc-1; j++)
        dc(j,i) = tmp*(c(j,i+1) - c(j,i));        
    }
      
  for ( octave_idx_type i(1); i <= nk-2; i++)
    dk(i-1) = k(i);
      
  if (nargout>1)
    retval(1) = octave_value(dk);
    retval(0) = octave_value(dc);

  return(retval);
}

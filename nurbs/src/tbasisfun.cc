/* Copyright (C) 2009 Carlo de Falco
   Copyright (C) 2012 Rafael Vazquez

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
#include <iostream>

void onebasisfun__ (double u, octave_idx_type p, RowVector U, double *N)
{
  *N = 0.0;
  if ((u <= U.min ()) || ( u > U.max ()))
    return;
  else if (p == 0)
    {
      *N = 1.0;
      return;
    }
  else if (p == 1) 
    {
      if (u < U(1)) 
        {
          *N = (u - U(0)) / (U(1) - U(0));
          return;
        }
      else 
        {
          *N = (U(2) - u) / (U(2) - U(1));
          return;
        }
    }
  else if (p == 2) 
    {
      double ln = u - U(0);
      double dn = U(3) - u;
      double ld = U(2) - U(0); 
      double dd = U(3) - U(1);
      if (u < U(1)) 
        {
          *N = ln*ln / (ld * (U(1) - U(0)));
          return;
        }
      else if (u > U(2)) 
        {
          *N = dn*dn / (dd * (U(3) - U(2)));
          return;
        }
      else 
        {
          if (ld != 0)            
            *N += ln * (U(2) - u) / ((U(2) - U(1)) * ld);
          
          if (dd != 0)
            *N += dn * (u - U(1)) / ((U(2) - U(1)) * dd);
          return;
        }
    }
 
  double ln = u - U(0);
  double ld = U(U.length () - 2) - U(0);
  if (ld != 0)
    {
      double tmp;
      onebasisfun__ (u, p-1, U.extract (0, U.length () - 2), &tmp);
      *N += ln * tmp / ld; 
    }
  double dn = U(U.length () - 1) - u;
  double dd = U(U.length () - 1) - U(1);
  if (dd != 0)
    {
      double tmp;
      onebasisfun__ (u, p-1, U.extract (1, U.length () - 1), &tmp);
      *N += dn * tmp / dd;
    }
  return;
}

void onebasisfun__ (double u, double p, RowVector U, double *N)
{ onebasisfun__ (u, static_cast<octave_idx_type> (p), U, N); }

void onebasisfunder__ (double u, octave_idx_type p, RowVector U, double *N, double *Nder)
{
  double aux;
  *N = 0.0; *Nder = 0.0;
  if ((u <= U.min ()) || ( u > U.max ()))       
    return;
  else if (p == 0)
    {
      *N = 1.0;
      *Nder = 0.0;
      return;
    }  
  else {    
    double ln = u - U(0);
    double ld = U(U.length () - 2) - U(0);

    if (ld != 0)
      {
        onebasisfun__ (u, p-1, U.extract (0, U.length () - 2), &aux);
        aux = aux / ld;
        *N += ln * aux;
        *Nder += p * aux;
      }
    
    double dn = U(U.length () - 1) - u;
    double dd = U(U.length () - 1) - U(1);

    if (dd != 0)
      { 
        onebasisfun__ (u, p-1, U.extract (1, U.length () - 1), &aux);
        aux = aux / dd;
        *N    += dn *aux;
        *Nder -= p * aux;
      } 
  }
}
   
DEFUN_DLD(tbasisfun, args, nargout,"\
TBASISFUN: Compute a B- or T-Spline basis function, and its derivatives, from its local knot vector.\n\
\n\
 usage:\n\
\n\
 [N, Nder] = tbasisfun (u, p, U)\n\
 [N, Nder] = tbasisfun ([u; v], [p q], {U, V})\n\
 [N, Nder] = tbasisfun ([u; v; w], [p q r], {U, V, W})\n\
 \n\
 INPUT:\n\
  u or [u; v] : points in parameter space where the basis function is to be\n\
  evaluated \n\
  \n\
  U or {U, V} : local knot vector\n\
\n\
  p or [p q] : polynomial order of the basis function\n\
\n\
 OUTPUT:\n\
  N : basis function evaluated at the given parametric points\n\
  Nder : gradient of the basis function evaluated at the given points\n")

{
  
  octave_value_list retval;
  Matrix u = args(0).matrix_value ();

  RowVector N(u.cols ());
  double *Nptr = N.fortran_vec ();

  if (! args(2).is_cell ())
    {

      double p = args(1).idx_type_value ();
      RowVector U = args(2).row_vector_value (true, true);
      assert (U.numel () == p+2);
      
      if (nargout == 1)
        for (octave_idx_type ii = 0; ii < u.numel (); ii++)
          onebasisfun__ (u(ii), p, U, &(Nptr[ii]));

      if (nargout == 2) 
        {
          RowVector Nder(u.cols ());
          double *Nderptr = Nder.fortran_vec ();

          for (octave_idx_type ii=0; ii<u.numel (); ii++)
            onebasisfunder__ (u(ii), p, U, &(Nptr[ii]), &(Nderptr[ii]));

          retval(1) = Nder;
        }      
    } 
  else 
    {
      RowVector p = args(1).row_vector_value ();
      if (p.length() == 2) 
        {
          Cell C = args(2).cell_value ();
          RowVector U = C(0).row_vector_value (true, true);
          RowVector V = C(1).row_vector_value (true, true);
          
          if (nargout == 1) 
            {
              for (octave_idx_type ii=0; ii<u.cols (); ii++)
                {
                  double Nu, Nv;
                  onebasisfun__ (u(0, ii), octave_idx_type(p(0)), U, &Nu);
                  onebasisfun__ (u(1, ii), octave_idx_type(p(1)), V, &Nv);
                  Nptr[ii] = Nu * Nv;
                }
            }
          else if (nargout == 2) 
            {
              double Nu, Nv, Ndu, Ndv;
              Matrix Nder (2, u.cols());
              double *Nderptr = Nder.fortran_vec ();
              for (octave_idx_type ii = 0; ii < u.cols (); ii++)
                {
                  onebasisfunder__ (u(0, ii), octave_idx_type(p(0)), U, &Nu, &Ndu);
                  onebasisfunder__ (u(1, ii), octave_idx_type(p(1)), V, &Nv, &Ndv);
                  Nptr[ii] = Nu * Nv;
                
                  Nderptr[0 + (2 * ii)] = Ndu * Nv;
                  Nderptr[1 + (2 * ii)] = Ndv * Nu;
                }
              retval(1) = Nder;
            }
        
        } 
      else if (p.length() == 3) 
        {
          Cell C = args(2).cell_value ();
          RowVector U = C(0).row_vector_value (true, true);
          RowVector V = C(1).row_vector_value (true, true);
          RowVector W = C(2).row_vector_value (true, true);
        
          if (nargout == 1)
            {
              for (octave_idx_type ii = 0; ii < u.cols (); ii++)
                {
                  double Nu, Nv, Nw;
                  onebasisfun__ (u(0, ii), octave_idx_type(p(0)), U, &Nu);
                  onebasisfun__ (u(1, ii), octave_idx_type(p(1)), V, &Nv);
                  onebasisfun__ (u(2, ii), octave_idx_type(p(2)), W, &Nw);
                  Nptr[ii] = Nu * Nv * Nw;
                }
            }
          else if (nargout == 2) 
            {
              double Nu, Nv, Nw, Ndu, Ndv, Ndw;
              Matrix Nder (3, u.cols());
              double *Nderptr = Nder.fortran_vec ();
              for (octave_idx_type ii=0; ii<u.cols (); ii++)
              {
                onebasisfunder__ (u(0, ii), octave_idx_type(p(0)), U, &Nu, &Ndu);
                onebasisfunder__ (u(1, ii), octave_idx_type(p(1)), V, &Nv, &Ndv);
                onebasisfunder__ (u(2, ii), octave_idx_type(p(2)), W, &Nw, &Ndw);
                Nptr[ii] = Nu * Nv * Nw;
                Nderptr[0 + (3 * ii)] = Ndu * Nv * Nw;
                Nderptr[1 + (3 * ii)] = Ndv * Nu * Nw;
                Nderptr[2 + (3 * ii)] = Ndw * Nu * Nv;
              }
            retval(1) = Nder;
          }
        
      }
    }
  retval(0) = octave_value (N);
  return retval;
}

/*

%!demo
%! U = {[0 0 1/2 1 1], [0 0 0 1 1]};
%! p = [3, 3];
%! [X, Y] = meshgrid (linspace(0, 1, 30));
%! u = [X(:), Y(:)]';
%! N = tbasisfun (u, p, U);
%! surf (X, Y, reshape (N, size(X)))
%! title('Basis function associated to a local knot vector')
%! hold off

%!test
%! U = [0 1/2 1];
%! p = 1;
%! u = [0.3 0.4 0.6 0.7];
%! [N, Nder] = tbasisfun (u, p, U);
%! assert (N, [0.6 0.8 0.8 0.6], 1e-12);
%! assert (Nder, [2 2 -2 -2], 1e-12);

%!test
%! U = {[0 1/2 1] [0 1/2 1]};
%! p = [1 1];
%! u = [0.3 0.4 0.6 0.7; 0.3 0.4 0.6 0.7];
%! [N, Nder] = tbasisfun (u, p, U);
%! assert (N, [0.36 0.64 0.64 0.36], 1e-12);
%! assert (Nder, [1.2 1.6 -1.6 -1.2; 1.2 1.6 -1.6 -1.2], 1e-12);

%!test
%! U = {[0 1/2 1] [0 1/2 1] [0 1/2 1]};
%! p = [1 1 1];
%! u = [0.4 0.4 0.6 0.6; 0.4 0.4 0.6 0.6; 0.4 0.6 0.4 0.6];
%! [N, Nder] = tbasisfun (u, p, U);
%! assert (N, [0.512 0.512 0.512 0.512], 1e-12);
%! assert (Nder, [1.28 1.28 -1.28 -1.28; 1.28 1.28 -1.28 -1.28; 1.28 -1.28 1.28 -1.28], 1e-12);

*/

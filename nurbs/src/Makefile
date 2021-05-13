## Copyright (C) 2009-2017 Carlo de Falco
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


OCTFILES=basisfun.oct \
basisfunder.oct \
bspderiv.oct \
bspeval.oct \
curvederivcpts.oct \
nrb_srf_basisfun__.oct \
nrb_srf_basisfun_der__.oct \
nrbsurfderiveval.oct \
surfderivcpts.oct \
surfderiveval.oct \
tbasisfun.oct 

MKOCTFILE ?= mkoctfile

all: $(OCTFILES)

low_level_functions.o: low_level_functions.cc
	$(MKOCTFILE) -c $<

%.oct:  %.cc low_level_functions.o
	$(MKOCTFILE) $< low_level_functions.o

clean:
	-rm -f *.o core octave-core *.oct *~


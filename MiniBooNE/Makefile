# GLoBES -- General LOng Baseline Experiment Simulator
# (C) 2002 - 2007,  The GLoBES Team
#
# GLoBES is mainly intended for academic purposes. Proper
# credit must be given if you use GLoBES or parts of it. Please
# read the section 'Credit' in the README file.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# The configure script should have replaced anything as required to
# obtain a working Makefile, which may serve as a template for own
# applications.
#
# This simple Makefile is for the GLoBES examples
#
# Compile example files with ``make example1'' etc.
#
# This Makefile assumes dynamic libraries, installed to either
# the default prefix /usr/local/ or to a user-defined directory 
# called ${prefix}.
#
# For linking against a specific version of GLoBES, libglobes.so can be 
# replaced by the respective library, such as libglobes.so.0.0.1



prefix = /usr/local
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

local_CFLAGS = -g -O4

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)



BIN = comb-dec comb-osc nue-osc anue-osc nue-dec anue-dec numu-osc anue-dec-only
OBJ = comb-dec.o comb-osc.o nue-osc.o anue-osc.o nue-dec.o anue-dec.o numu-osc.o anue-dec-only.o myio.o

all: $(BIN)
comb-dec: comb-dec.o myio.o
	gcc comb-dec.o myio.o -o comb-dec $(LDFLAGS) $(local_LDFLAGS)

comb-osc: comb-osc.o myio.o
	gcc comb-osc.o myio.o -o comb-osc $(LDFLAGS) $(local_LDFLAGS)

anue-dec-only: anue-dec-only.o myio.o
	gcc anue-dec-only.o myio.o -o anue-dec-only $(LDFLAGS) $(local_LDFLAGS)

numu-osc: numu-osc.o myio.o
	gcc numu-osc.o myio.o -o numu-osc $(LDFLAGS) $(local_LDFLAGS)

nue-osc: nue-osc.o myio.o
	gcc nue-osc.o myio.o -o nue-osc $(LDFLAGS) $(local_LDFLAGS)

anue-osc: anue-osc.o myio.o
	gcc anue-osc.o myio.o -o anue-osc $(LDFLAGS) $(local_LDFLAGS)

nue-dec: nue-dec.o myio.o
	gcc nue-dec.o myio.o -o nue-dec $(LDFLAGS) $(local_LDFLAGS)

anue-dec: anue-dec.o myio.o
	gcc anue-dec.o myio.o -o anue-dec $(LDFLAGS) $(local_LDFLAGS)

%.o : %.c
	gcc $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)
.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ)

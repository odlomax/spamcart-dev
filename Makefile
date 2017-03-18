# The MIT License (MIT)
# 
# Copyright (c) 2016 odlomax
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

SHELL=/bin/sh

# LIB_DIR=/software/libraries/fftw/3.3.4/gnu-4.4.7/lib
# INC_DIR=/software/libraries/fftw/3.3.4/gnu-4.4.7/include

LIB_DIR=/usr/local/lib/
INC_DIR=/usr/local/include/

FFTW_LIB=fftw3

COMPILER=gfortran

EXECUTABLE=spamcart

OPTIONS=-Wall -std=f2008 -pedantic

OPTIMISE_OPT=-O3 -flto -march=native
DEBUG_OPT=-g3 -Og -fcheck=all -fbacktrace -ffpe-trap=invalid,overflow -Wextra
PARALLEL_OPT=-fopenmp
PROFILE_OPT=-pg

# OPTIONS+= $(PROFILE_OPT)
OPTIONS+= $(OPTIMISE_OPT)
OPTIONS+= $(PARALLEL_OPT)
# OPTIONS+= $(DEBUG_OPT)


SOURCE_DIR=./src/

OBJECTS=parameters.o maths.o io.o string.o options.o atomic_update.o dust.o dust_d03.o triangular_array.o particle.o kernel.o kernel_m4.o line.o ms_star_sed.o source.o source_point_bb.o source_point_custom.o source_point_ms_star.o source_external_bb.o source_external_ps05.o binary_tree_node.o binary_tree.o ray.o image_tree_node.o datacube.o simulation.o main.o

COMPILE_LINE=$(COMPILER) $(OPTIONS) -I$(INC_DIR) -c $(SOURCE_DIR)

spamcart: $(OBJECTS) $(MODULES)
	$(COMPILER) $(OPTIONS) -o $(EXECUTABLE) $(OBJECTS) -L$(LIB_DIR) -l$(FFTW_LIB)

parameters.o: $(SOURCE_DIR)parameters.f90
	$(COMPILE_LINE)parameters.f90

io.o: $(SOURCE_DIR)io.f90
	$(COMPILE_LINE)io.f90

atomic_update.o: $(SOURCE_DIR)atomic_update.f90
	$(COMPILE_LINE)atomic_update.f90

maths.o: $(SOURCE_DIR)maths.f90
	$(COMPILE_LINE)maths.f90

string.o: $(SOURCE_DIR)string.f90
	$(COMPILE_LINE)string.f90
	
options.o: $(SOURCE_DIR)options.f90
	$(COMPILE_LINE)options.f90
	
dust.o: $(SOURCE_DIR)dust.f90
	$(COMPILE_LINE)dust.f90

dust_d03.o: $(SOURCE_DIR)dust_d03.f90
	$(COMPILE_LINE)dust_d03.f90

triangular_array.o: $(SOURCE_DIR)triangular_array.f90
	$(COMPILE_LINE)triangular_array.f90

particle.o: $(SOURCE_DIR)particle.f90
	$(COMPILE_LINE)particle.f90

kernel.o: $(SOURCE_DIR)kernel.f90
	$(COMPILE_LINE)kernel.f90

kernel_m4.o: $(SOURCE_DIR)kernel_m4.f90
	$(COMPILE_LINE)kernel_m4.f90

line.o: $(SOURCE_DIR)line.f90
	$(COMPILE_LINE)line.f90

binary_tree_node.o: $(SOURCE_DIR)binary_tree_node.f90
	$(COMPILE_LINE)binary_tree_node.f90

binary_tree.o: $(SOURCE_DIR)binary_tree.f90
	$(COMPILE_LINE)binary_tree.f90

ray.o: $(SOURCE_DIR)ray.f90
	$(COMPILE_LINE)ray.f90

ms_star_sed.o: $(SOURCE_DIR)ms_star_sed.f90
	$(COMPILE_LINE)ms_star_sed.f90

source.o: $(SOURCE_DIR)source.f90
	$(COMPILE_LINE)source.f90 -Wno-unused-dummy-argument

source_point_bb.o: $(SOURCE_DIR)source_point_bb.f90
	$(COMPILE_LINE)source_point_bb.f90 -Wno-unused-dummy-argument
	
source_point_custom.o: $(SOURCE_DIR)source_point_custom.f90
	$(COMPILE_LINE)source_point_custom.f90 -Wno-unused-dummy-argument

source_point_ms_star.o: $(SOURCE_DIR)source_point_ms_star.f90
	$(COMPILE_LINE)source_point_ms_star.f90 -Wno-unused-dummy-argument

source_external_bb.o: $(SOURCE_DIR)source_external_bb.f90
	$(COMPILE_LINE)source_external_bb.f90 -Wno-unused-dummy-argument
	
source_external_ps05.o: $(SOURCE_DIR)source_external_ps05.f90
	$(COMPILE_LINE)source_external_ps05.f90 -Wno-unused-dummy-argument

image_tree_node.o: $(SOURCE_DIR)image_tree_node.f90
	$(COMPILE_LINE)image_tree_node.f90

datacube.o: $(SOURCE_DIR)datacube.f90
	$(COMPILE_LINE)datacube.f90

simulation.o: $(SOURCE_DIR)simulation.f90
	$(COMPILE_LINE)simulation.f90

main.o: $(SOURCE_DIR)main.f90
	$(COMPILE_LINE)main.f90


clean:
	rm *.o *.mod
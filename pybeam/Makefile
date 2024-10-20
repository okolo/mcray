# Compilers
CC=gcc
CXX=c++

# Python paths
PYINCLUDE=$(shell python3 -m pybind11 --includes)

# CRbeam paths
# -I</path/to/CRbeam/headers/> -L</path/to/compiled/CRbeam/lib -lMCPropag -fopenmp -lgsl
CRBEAM_FLAGS=-I../src/lib/ -I../src/external/ -L../bin/mcray/

all: pybeam test

test: pybeam
	python3 test.py

pybind11:
	pip3 install pybind11

pybeam:
	$(CXX) -O3 -Wall -shared -std=c++11 -fPIC $(PYINCLUDE) $(CRBEAM_FLAGS) main.cpp -o pybeam.so -ldl -lmcray -fopenmp -lgsl -lgfortran
#	python3 setup.py install

clean:
	rm pybeam.so

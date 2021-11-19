#!/bin/bash

#Fortran compilation
#ifort -c mathtools.f90

#Python compilation of fortran code
python3 -m numpy.f2py -c mathtools.f90 -m mathtools --fcompiler=intelem --quiet

#Include in python:
#import mathtools
#mathtools.padeseries() (das nach dem Punkt ist der Name der subroutine, aufpassen mit gross- und kleinschreibung)

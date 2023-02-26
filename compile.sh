#! /bin/bash

execfile=tunelj-8feb2021.x

rm ${execfile}

gfortran -std=f2008 -c global/myglobal.F08
gfortran -std=f2008 -c sysutils/sysutils.F08
gfortran -std=f2008 -c parser/myparser.F08
gfortran -std=f2008 -c grotopo/grotopo.F08
gfortran -std=f2008 myglobal.o myparser.o grotopo.o sysutils.o main.F08 -o ${execfile}

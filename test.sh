#! /bin/bash

clear
rm g_nitrile_hbond examples/*{xvg,log} *xvg *log
make -j 8
if [ ! -f g_nitrile_hbond ] ; then
    exit
    fi
if [ ! -z $1 ] ; then
    ./g_nitrile_hbond -h
else
time ./g_nitrile_hbond -s examples/tpr.tpr -f examples/xtc.xtc -a1 351 -a2 352 -select "resname SOL and same residue as within 0.5 of resname CNC and name NE" -oa -op
    fi


#!/bin/sh
swig -python -Wall designspacetoolbox_interface.i
cp designspacetoolbox_interface_wrap.c dspace/SWIG/designspacetoolbox_wrap.c
cp dspace_interface.py dspace/SWIG/

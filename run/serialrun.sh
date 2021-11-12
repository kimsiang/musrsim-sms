#!/bin/bash
#   i=730044
#   ../bin/Linux-g++/musrSim $i.mac
#   i=730050
#   ../bin/Linux-g++/musrSim $i.mac
  
   i=720004
   while [ "$i" -lt 720016 ]; do
    ../bin/Linux-g++/musrSim $i.mac
   i=$(expr $i + 1)
   done

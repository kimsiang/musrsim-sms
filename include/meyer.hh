/***************************************************************************
 *  musrSim - the program for the simulation of (mainly) muSR instruments. *
 *          More info on http://lmu.web.psi.ch/simulation/index.html .     *
 *          musrSim is based od Geant4 (http://geant4.web.cern.ch/geant4/) *
 *                                                                         *
 *  Copyright (C) 2009 by Paul Scherrer Institut, 5232 Villigen PSI,       *
 *                                                       Switzerland       *
 *                                                                         *
 *  This program is free software; you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation; either version 2 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program; if not, write to the Free Software            *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              *
 ***************************************************************************/

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$//
//
//  &           &&&&&&&&&&                            &&&&&&&         &&&&&&&&
//  &           &                                   &&      &&       &       &&
//  &           &                   &      &       &                &       &&
//  &           &&&&&&&            &      &         &&&&&&         &&&&&&&&  
//  &           &                 &      &&                &      &      &&
//  &           &                &&     &  &     &&      &&      &        &
//  &&&&&&&&&&  &&&&&&&&&&      &  &&&&&    &&    &&&&&&&       &        &&   
//                             &
//                            &
//                           &
//                         &  
//                             MEYER
/*
  fIRST IMPLEMENTATION BY ANLSEM,H. IN FORTRAN
  C++ CONVERSION T.K.PARAISO 04-2005

  !!! IMPORTANT  !!!
 
  Notice: 
  Tables definition changes between FORTRAN and C++:
  1/ Fortran indices start at 1 and C++ indices start at 0
  2/ Tables are defined as table[column][row] in Fortran
  table[row][column] in c++

  usefull reference 
  http://gershwin.ens.fr/vdaniel/Doc-Locale/Langages-Program-Scientific/Fortran/Tutorial/arrays.htm

*/
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$// 


#ifndef meyer_h
#define meyer_h 1

#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <ios>
#include "globals.hh"

class meyer
{
 public:
  meyer();
 ~meyer();

 void GFunctions(double*, double*, double);
 void Get_F_Function_Meyer(double tau, double Ekin, double Z1, double Z2, double m1, double m2);
 void F_Functions_Meyer( double tau,double thetaSchlange,double *f1,double *f2);
 void Get_F_Function(double tau,double theta, double Ekin, double Z1, double Z2, double m1, double m2, double* F);

};

#endif

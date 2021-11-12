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

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$//
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


#include "meyer.hh"
#include <iomanip>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <ios>
// James Lord:
// already defines its own constant Pi so no need to look for M_PI anywhere
// #define _USE_MATH_DEFINES 1
// #include <math.h>
using namespace std;

meyer::meyer()
{;}

meyer::~meyer()
{;}

void meyer::GFunctions(double* g1,double* g2, double tau)
{

  //Diese Routine gibt in Abhaengigkeit von der reduzierten Dicke 'tau'
  //Funktionswerte fuer g1 und g2 zurueck. g1 und g2 sind dabei die von
  //Meyer angegebenen tabellierten Funktionen fuer die Berechnung von Halbwerts-
  //breiten von Streuwinkelverteilungen. (L.Meyer, phys.stat.sol. (b) 44, 253
  //(1971))


  double help;

  int i;


  double tau_[] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
		   2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0,
		   10.0, 12.0, 14.0, 16.0, 18.0, 20.0 };

  double g1_[]  = {0.050,0.115,0.183,0.245,0.305,0.363,0.419,0.473,0.525,0.575,
		   0.689,0.799,0.905,1.010,1.100,1.190,1.370,1.540,1.700,1.850,
		   1.990,2.270,2.540,2.800,3.050,3.290 };

  double g2_[]  = {0.00,1.25,0.91,0.79,0.73,0.69,0.65,0.63,0.61,0.59,
		   0.56,0.53,0.50,0.47,0.45,0.43,0.40,0.37,0.34,0.32,
		   0.30,0.26,0.22,0.18,0.15,0.13 };


  if (tau<tau_[0])// tau_[0] is the lowest in c++; in fortran it is tau_[1]!!! TAO!
    {
      std::cout<<"SUBROUTINE G_Functions:"<<std::endl;
      std::cout<<"  Fehler bei Berechnung der g-Funktionen fuer Winkelaufstreuung:"<<std::endl;
      std::cout<<"  aktuelles tau ist kleiner als kleinster Tabellenwert:"<<std::endl;
      std::cout<<"  tau     = "<< tau<<std::endl;
      std::cout<<"  tau_[0] = "<< tau_[0]<<std::endl;
      return;
    }
  i = 1;

  do
    {
      i = i + 1;
      if (i==26)
	{
	  std::cout<<"SUBROUTINE G_Functions:"<<std::endl;
	  std::cout<<"  Fehler bei Berechnung der g-Funktionen fuer Winkelaufstreuung:"<<std::endl;
	  std::cout<<"  aktuelles tau ist groesser als groesster Tabellenwert:"<<std::endl;
	  std::cout<<"  tau      = "<< tau <<std::endl;
	  std::cout<<"  tau_[26] = "<< tau_[26] <<std::endl;
	  break;
	}
    }while(tau>tau_[i]);


  //lineare Interpolation zwischen Tabellenwerten:

  help = (tau-tau_[i-1])/(tau_[i]-tau_[i-1]);

  *g1 = g1_[i-1] + help*(g1_[i]-g1_[i-1]);
  *g2 = g2_[i-1] + help*(g2_[i]-g2_[i-1]);



}



//==========================================================================================



void meyer::Get_F_Function_Meyer(double tau, double Ekin, double Z1, double Z2, double m1, double m2)
{

  double thetaSchlange,thetaSchlangeMax;
  double theta,thetaMax,thetaStep;
  double f1,f2,F;
  
  
  //---------------------------------
  //- Parameters:
  
  //  double Z1, Z2;  ! die atomaren Nummern von Projektil und Target

  double a0;                  //      ! Bohrscher Radius in cm
  double screeningPar;        //      ! Screeningparameter "a" in cm fuer Teilchen der
                              //      ! Kernladungszahl Z1=1 in Kohlenstoff (Z2 = 6)
                              //      ! bei Streichung von Z1 (vgl. Referenz, S. 268)
  double r0Meyer;             //      ! r0(C) berechnet aus dem screeningParameter "a"
                              //      ! und dem ebenfalls bei Meyer angegebenem
                              //      ! Verhaeltnis a/r0=0.26 (vgl. Referenz, S. 263 oben)
  double eSquare;             //      ! elektrische Ladung zum Quadrat in keV*cm
  double Pi	;             //      ! die Kreiszahl
  
  ///	
  a0 = 5.29E-9;//unit == centimeter
  
  //the screening parameter
  double D= exp(2/3*log(Z1))+exp(2/3*log(Z2));
  double a=0.885*a0/sqrt(D);
  screeningPar=a;      //  screeningPar = 2.5764E-9;
  r0Meyer = 9.909E-9;  
  eSquare = 1.44E-10;
  Pi = 3.141592654;
  
  
  double Meyer_faktor3;
  double Meyer_faktor4;
  double zzz;//	    ! "Hilfsparameter"
  double Meyer_faktor5;
  
  Meyer_faktor3 = (screeningPar/r0Meyer) * (screeningPar/r0Meyer);
  
  Meyer_faktor4 = (m1+m2)/m2/2.;
    //((1./9.+12.)/12.)/2. ;// TAO m1+m2/m2/2.
  // in meyer article, we then have b= mf4/Ekine1= (m1+m2)/(m2*2*m1*v1²)
  zzz = screeningPar / (2.*Z1*Z2*eSquare);
  Meyer_faktor5 = zzz*zzz / (8*Pi*Pi);
  


  //---------------------------------
  //---------------------------------
  //---------------------------------

  
  
  int nBin,nBinMax;
  nBinMax=201;
  double value[201]; //    /0.,nBinMax*0./
  double area[201]  ;  //    /   nBinMax*0./
  double integ[201]; //      /0.,nBinMax*0./
  
  
  //    common /MeyerTable/ value,area,integ,thetaStep,nBin
  
  int i;
  double rHelp;
  
  int   HB_memsize;
  HB_memsize=500000;
  //  double memory[500000];
  //        COMMON /PAWC/ memory


//nur noch fuer Testzwecke:

  double fValues[203];
  double fValuesFolded[203];
  
  int idh;
  idh = 50;

  //INCLUDE "mutrack$sourcedirectory:COM_DIRS.INC"
  //	character filename*20           ! Name der Ausgabe-Dateien
  //	COMMON /filename/ filename
  
//----------------------------------------------------------------------------
  
//Festlegen des maximalen Theta-Wertes sowie der Schrittweite:

  if (tau<0.2) 
    {
      std::cout<< "Subroutine ''Get_F_Function_Meyer'':"<<std::endl;
      std::cout<< "Effektive Dicke ist kleiner als 0.2 => kann ich nicht ... => STOP"<<std::endl;
      return;
    }
  else if (tau<=2.)
    {
      //	    ! => Tabelle A
      thetaSchlangeMax = 4.0;
    }
  else if (tau<=8.)
    {
      //! => Tabelle B
      thetaSchlangeMax = 7.0;
    }
  else if (tau<=20.)
    {
      //! => Tabelle C
      thetaSchlangeMax = 20.0;
    }
  else
    {
      std::cout<< "Subroutine ''Get_F_Function_Meyer'':"<<std::endl;
      std::cout<< "Effektive Dicke ist groesser als 20 => kann ich nicht ... => STOP"<<std::endl;
      return;
    }
  //  std::cout<< "M4:   "<<Meyer_faktor4<<std::endl;
  //  std::cout<< "Ekin: "<<Ekin <<std::endl;
  thetaMax = thetaSchlangeMax / Meyer_faktor4 / Ekin/Pi*180;
  if (thetaMax>50.)
    {
      thetaStep = .5;
    }

  else if (thetaMax>25)
    {
      thetaStep = .25;
    }
  else if (thetaMax>12.5)
    {
      thetaStep = .125;
    }
  else
    {
      thetaStep = .0625;
    }


  //Tabelle der F-Werte erstellen:

  nBin = 0;
  std::cout<<"thetamax  = "<<thetaMax << std::endl;


  theta=thetaStep;
  // begining of do loop
  for( theta = thetaStep; theta<=thetaMax; theta+=thetaStep)
    {
      //      std::cout<<"theta"<<theta << std::endl;
      //	    ! Berechne aus theta das 'reduzierte' thetaSchlange (dabei gleich
      //	    ! noch von degree bei theta in Radiant bei thetaSchlange umrechnen):
      //
      thetaSchlange = Meyer_faktor4 * Ekin * theta  *Pi/180;

      //  ! Auslesen der Tabellenwerte fuer die f-Funktionen:

      F_Functions_Meyer(tau,thetaSchlange,&f1,&f2);
      if (thetaSchlange==-1)
	{
	  //! wir sind jenseits von thetaSchlangeMax
	  goto bigtheta;
	  //    endif
	}

      //	    ! Berechnen der Streuintensitaet:
      F = Meyer_faktor4*Meyer_faktor4 * Ekin*Ekin /2 /Pi * (f1 - Meyer_faktor3*f2);// TAO, Anselm was: Meyer_faktor5 * Ekin*Ekin * (f1 - Meyer_faktor3*f2);

      nBin = nBin + 1;
      if (nBin>nBinMax)
	{
	  std::cout<< "nBin > nBinMax  =>  EXIT";
	  break;
	}

      value[nBin] = sin(theta)*F;

      fValues[nBin+1] = F;                  //    ! fuer Testzwecke
      fValuesFolded[nBin+1] = sin(theta/180*Pi)*F;//    ! fuer Testzwecke


    }// end of do loop


  //Berechnen der Flaecheninhalte der einzelnen Kanaele sowie der Integrale:

 bigtheta:for( i = 1;i<= nBin; i++)
   {
     area[i]  = (value[i]+value[i-1])/2.* thetaStep;
     integ[i] = integ[i-1] + area[i];
   }


  //Normiere totale Flaeche auf 1:

  rHelp = integ[nBin];
  for( i = 1; i<=nBin; i++)
    {
      value[i] = value[i] / rHelp;
      area[i]  = area[i]  / rHelp;
      integ[i] = integ[i] / rHelp;
    }


  //vorerst noch: gib Tabelle in Datei und Histogrammfile aus:

  //! Berechne die Werte fuer theta=0:

  F_Functions_Meyer(tau,0.,&f1,&f2);
  F = Meyer_faktor4*Meyer_faktor4 * Ekin*Ekin /2 /Pi * (f1 - Meyer_faktor3*f2);// TAO, Anselm was: Meyer_faktor5 * Ekin*Ekin *          (f1 - Meyer_faktor3*f2);
  fValues[1]       = F;
  fValuesFolded[1] = 0.;

  //! Gib die Werte in das Tabellenfile aus:

  /*    ofstream Mprint("testmeyer.out");
 theta = thetaStep;
  if (!Mprint.is_open()) exit(8);
   for( i = 1; i<=nBin+1;i++)
     {
       Mprint << theta<< " "<< fValues[i]/fValues[1]<<" " << fValuesFolded[i]<<std::endl;
       theta = theta + thetaStep;
     }

   Mprint.close();
  */
}

//=============================================================================================================================================================//
void meyer:: Get_F_Function(double tau,double theta, double Ekin, double Z1, double Z2, double m1, double m2, double* F)
{

  double thetaSchlange; //thetaSchlangeMax;
  double thetaStep; //thetaMax
  double f1,f2;


  //---------------------------------
  //- Parameters:

  //  double Z1, Z2;  ! die atomaren Nummern von Projektil und Target

  double a0;                  //      ! Bohrscher Radius in cm
  double screeningPar;        //      ! Screeningparameter "a" in cm fuer Teilchen der
                              //      ! Kernladungszahl Z1=1 in Kohlenstoff (Z2 = 6)
                              //      ! bei Streichung von Z1 (vgl. Referenz, S. 268)
  double r0Meyer;             //      ! r0(C) berechnet aus dem screeningParameter "a"
                              //      ! und dem ebenfalls bei Meyer angegebenem
                              //      ! Verhaeltnis a/r0=0.26 (vgl. Referenz, S. 263 oben)
  double eSquare;             //      ! elektrische Ladung zum Quadrat in keV*cm
  double Pi;                  //      ! die Kreiszahl

  ///
  a0 = 5.29E-9;//unit == centimeter
  
  //the screening parameter
  double D= exp(2/3*log(Z1))+exp(2/3*log(Z2));
  double a=0.885*a0/sqrt(D);
  screeningPar=a;      //  screeningPar = 2.5764E-9;
  r0Meyer = 9.909E-9;  
  eSquare = 1.44E-10;
  Pi = 3.141592654;
  
  
  double Meyer_faktor3;
  double Meyer_faktor4;
  double zzz;//	    ! "Hilfsparameter"
  double Meyer_faktor5;
  
  Meyer_faktor3 = (screeningPar/r0Meyer) * (screeningPar/r0Meyer);
  
  Meyer_faktor4 = (m1+m2)/m2/2.;
    //((1./9.+12.)/12.)/2. ;// TAO m1+m2/m2/2.
  // in meyer article, we then have b= mf4/Ekine1= (m1+m2)/(m2*2*m1*v1²)
  zzz = screeningPar / (2.*Z1*Z2*eSquare);
  Meyer_faktor5 = zzz*zzz / (8*Pi*Pi);
  


  //---------------------------------
  //---------------------------------
  //---------------------------------

  
  
  int nBin,nBinMax;
  nBinMax=201;
  double value[201]; //    /0.,nBinMax*0./
  double area[201]  ;  //    /   nBinMax*0./
  double integ[201]; //      /0.,nBinMax*0./
  
  
  //    common /MeyerTable/ value,area,integ,thetaStep,nBin
  
  int i;
  //double rHelp;
  
  int   HB_memsize;
  HB_memsize=500000;
  //  double memory[500000];
  //        COMMON /PAWC/ memory


//nur noch fuer Testzwecke:

  //double fValues[203];
  //double fValuesFolded[203];
  
  int idh;
  idh = 50;

  //INCLUDE "mutrack$sourcedirectory:COM_DIRS.INC"
  //	character filename*20           ! Name der Ausgabe-Dateien
  //	COMMON /filename/ filename
  
//----------------------------------------------------------------------------
  
//Festlegen des maximalen Theta-Wertes sowie der Schrittweite:

  
  //Tabelle der F-Werte erstellen:
  
  nBin = 0;
  
      thetaSchlange = Meyer_faktor4 * Ekin * theta  *Pi/180;
      
      
      F_Functions_Meyer(tau,thetaSchlange,&f1,&f2);
      if (thetaSchlange==-1)
	{
	  goto bigtheta;
	}
      
      *F = Meyer_faktor4*Meyer_faktor4 * Ekin*Ekin /2 /Pi * (f1 - Meyer_faktor3*f2);
  
 bigtheta:for( i = 1;i<= nBin; i++)
   {
     area[i]  = (value[i]+value[i-1])/2.* thetaStep;
     integ[i] = integ[i-1] + area[i];
   }
  

  //Normiere totale Flaeche auf 1:

  //! Berechne die Werte fuer theta=0:
      double norm;
  F_Functions_Meyer(tau,0.,&f1,&f2);
  norm = Meyer_faktor4*Meyer_faktor4 * Ekin*Ekin /2 /Pi * (f1 - Meyer_faktor3*f2);
  *F=*F/norm;
 
}

//===============================================================================================

void meyer:: F_Functions_Meyer( double tau,double thetaSchlange,double *f1,double *f2)
{

  //Diese Routine gibt in Abhaengigkeit von 'thetaSchlange' und 'tau'
  //Funktionswerte fuer f1 und f2 zurueck. f1 und f2 entsprechen dabei den
  //bei Meyer angegebenen Funktion gleichen Namens. Die in dieser Routine
  //verwendeten Tabellen sind eben dieser Referenz entnommen:
  //L.Meyer, phys.stat.sol. (b) 44, 253 (1971)

  double  f1_[2], f2_[2];
  
  int column_,column,row;
  
  int iColumn;
  double weightCol, weightRow;
  
  //----------------------------------------------------------------------------
  
  //die Tabellendaten der Referenz (Tabellen 2 und 3):
  
  int nColumn;
  nColumn=24;
  
  double tau_[25]=    {
    0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0,
    3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10., 12., 14., 16., 18., 20. 
  };
  
  
  int nRowA=24;
  double thetaSchlangeA[25]=
    {
      .00, .05, .10, .15, .20, .25, .30, .35, .40, .45, .50, .60,
      .70, .80, .90, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0 
    };
  
  int nRowB=23;
  double thetaSchlangeB[24]=
    {
      0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.8,
      2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0	  
    };
  
  int nRowC=23;
  double thetaSchlangeC[24]=
    {
      0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0,
      7.0, 8.0, 9.0, 10., 11., 12., 13., 14., 15., 16., 18., 20.
    };
  
  
  double f1_A[25][9]=
    { 
      {1.96E+2,4.55E+1,2.11E+1,1.25E+1,8.48E+0,6.21E+0,4.80E+0,3.86E+0,3.20E+0},
      {9.82E+1,3.72E+1,1.97E+1,1.20E+1,8.27E+0,6.11E+0,4.74E+0,3.83E+0,3.17E+0},
      {3.96E+1,2.58E+1,1.65E+1,1.09E+1,7.73E+0,5.82E+0,4.58E+0,3.72E+0,3.10E+0},
      {1.76E+1,1.58E+1,1.27E+1,9.26E+0,6.93E+0,5.38E+0,4.31E+0,3.55E+0,2.99E+0},
      {8.62E+0,1.01E+1,9.45E+0,7.58E+0,6.02E+0,4.85E+0,3.98E+0,3.33E+0,2.84E+0},
      {4.65E+0,6.55E+0,6.91E+0,6.06E+0,5.11E+0,4.28E+0,3.62E+0,3.08E+0,2.66E+0},
      {2.74E+0,4.45E+0,5.03E+0,4.78E+0,4.27E+0,3.72E+0,3.23E+0,2.82E+0,2.47E+0},
      {1.77E+0,3.02E+0,3.71E+0,3.76E+0,3.53E+0,3.20E+0,2.86E+0,2.55E+0,2.27E+0},
      {1.22E+0,2.19E+0,2.78E+0,2.96E+0,2.91E+0,2.73E+0,2.51E+0,2.28E+0,2.07E+0},
      {8.82E-1,1.59E+0,2.12E+0,2.35E+0,2.39E+0,2.32E+0,2.19E+0,2.03E+0,1.87E+0},
      {6.55E-1,1.20E+0,1.64E+0,1.88E+0,1.97E+0,1.96E+0,1.90E+0,1.79E+0,1.68E+0},
      {3.80E-1,7.15E-1,1.01E+0,1.22E+0,1.35E+0,1.40E+0,1.41E+0,1.39E+0,1.34E+0},
      {2.26E-1,4.45E-1,6.44E-1,8.08E-1,9.28E-1,1.01E+0,1.05E+0,1.06E+0,1.05E+0},
      {1.39E-1,2.80E-1,4.21E-1,5.45E-1,6.46E-1,7.22E-1,7.75E-1,8.07E-1,8.21E-1},
      {8.22E-2,1.76E-1,2.78E-1,3.71E-1,4.53E-1,5.21E-1,5.74E-1,6.12E-1,6.37E-1},
      {5.04E-2,1.11E-1,1.86E-1,2.57E-1,3.22E-1,3.79E-1,4.27E-1,4.65E-1,4.94E-1},
      {2.51E-2,5.60E-2,9.24E-2,1.31E-1,1.69E-1,2.02E-1,2.40E-1,2.71E-1,2.97E-1},
      {1.52E-2,3.20E-2,5.08E-2,7.23E-2,9.51E-2,1.18E-1,1.41E-1,1.63E-1,1.83E-1},
      {1.03E-2,2.05E-2,3.22E-2,4.55E-2,6.01E-2,7.53E-2,9.02E-2,1.05E-1,1.19E-1},
      {8.80E-3,1.48E-2,2.25E-2,3.13E-2,4.01E-2,5.03E-2,6.01E-2,7.01E-2,8.01E-2},
      {6.10E-3,1.15E-2,1.71E-2,2.28E-2,2.89E-2,3.52E-2,4.18E-2,4.86E-2,5.55E-2},
      {0.00   ,0.00   ,0.00   ,0.00   ,0.00   ,1.71E-2,1.98E-2,2.28E-2,2.58E-2},
      {0.00   ,0.00   ,0.00   ,0.00   ,0.00   ,8.90E-3,1.02E-2,1.16E-2,1.31E-2},
      {0.00   ,0.00   ,0.00   ,0.00   ,0.00   ,4.90E-3,5.70E-3,6.40E-3,7.20E-3},
      {0.00   ,0.00   ,0.00   ,0.00   ,0.00   ,2.90E-3,3.40E-3,3.90E-3,4.30E-3}
    };

  double f1_B[24][9]=
    {
      {2.71E+0,1.92E+0,1.46E+0,1.16E+0,9.52E-1,8.03E-1,6.90E-1,5.32E-1,4.28E-1},
      {2.45E+0,1.79E+0,1.39E+0,1.12E+0,9.23E-1,7.82E-1,6.75E-1,5.23E-1,4.23E-1},
      {1.87E+0,1.48E+0,1.20E+0,9.96E-1,8.42E-1,7.24E-1,6.32E-1,4.98E-1,4.07E-1},
      {1.56E+0,1.30E+0,1.09E+0,9.19E-1,7.89E-1,6.86E-1,6.03E-1,4.80E-1,3.95E-1},
      {1.28E+0,1.11E+0,9.62E-1,8.33E-1,7.27E-1,6.40E-1,5.69E-1,4.59E-1,3.81E-1},
      {8.23E-1,7.90E-1,7.29E-1,6.64E-1,6.01E-1,5.44E-1,4.94E-1,4.12E-1,3.49E-1},
      {5.14E-1,5.36E-1,5.29E-1,5.07E-1,4.78E-1,4.47E-1,4.16E-1,3.60E-1,3.13E-1},
      {3.19E-1,3.58E-1,3.76E-1,3.78E-1,3.70E-1,3.57E-1,3.45E-1,3.08E-1,2.76E-1},
      {2.02E-1,2.40E-1,2.64E-1,2.77E-1,2.82E-1,2.80E-1,2.65E-1,2.59E-1,2.39E-1},
      {1.67E-1,1.96E-1,2.20E-1,2.36E-1,2.44E-1,2.47E-1,2.45E-1,2.35E-1,2.21E-1},
      {1.33E-1,1.61E-1,1.85E-1,2.02E-1,2.12E-1,2.18E-1,2.18E-1,2.14E-1,2.03E-1},
      {8.99E-2,1.12E-1,1.32E-1,1.48E-1,1.59E-1,1.67E-1,1.68E-1,1.75E-1,1.72E-1},
      {6.24E-2,7.94E-2,9.50E-2,1.09E-1,1.20E-1,1.29E-1,1.35E-1,1.42E-1,1.43E-1},
      {4.55E-2,5.74E-2,6.98E-2,8.11E-2,9.09E-2,9.92E-2,1.06E-1,1.15E-1,1.19E-1},
      {3.35E-2,4.22E-2,5.19E-2,6.11E-2,6.95E-2,7.69E-2,8.33E-2,9.28E-2,9.85E-2},
      {2.50E-2,3.16E-2,3.92E-2,4.66E-2,5.35E-2,6.00E-2,6.57E-2,7.49E-2,8.13E-2},
      {1.90E-2,2.40E-2,2.99E-2,3.58E-2,4.16E-2,4.70E-2,5.20E-2,6.05E-2,6.70E-2},
      {1.47E-2,1.86E-2,2.32E-2,2.79E-2,3.25E-2,3.70E-2,4.12E-2,4.89E-2,5.51E-2},
      {8.10E-3,1.04E-2,1.30E-2,1.57E-2,1.84E-2,2.12E-2,2.40E-2,2.93E-2,3.42E-2},
      {4.80E-3,6.20E-3,7.70E-3,9.30E-3,1.09E-2,1.26E-2,1.44E-2,1.79E-2,2.14E-2},
      {2.80E-3,3.80E-3,4.70E-3,5.70E-3,6.70E-3,7.50E-3,8.90E-3,1.13E-2,1.36E-2},
      {1.70E-3,2.30E-3,2.90E-3,3.60E-3,4.20E-3,4.90E-3,5.60E-3,7.20E-3,8.80E-3},
      {0.00   ,0.00   ,0.00   ,0.00   ,0.00   ,0.00   ,2.00E-3,2.80E-3,3.50E-3},
      {0.00   ,0.00   ,0.00   ,0.00   ,0.00   ,0.00   ,8.80E-4,1.20E-3,1.60E-3}
    };
  
  double f1_C[24][7]=
    {
      { 3.65E-1,2.62E-1,2.05E-1,1.67E-1,1.41E-1,1.21E-1,1.05E-1},
      { 3.33E-1,2.50E-1,1.95E-1,1.61E-1,1.36E-1,1.18E-1,1.03E-1},
      { 2.75E-1,2.18E-1,1.76E-1,1.48E-1,1.27E-1,1.11E-1,9.80E-2},
      { 2.04E-1,1.75E-1,1.50E-1,1.29E-1,1.13E-1,1.01E-1,9.00E-2},
      { 1.41E-1,1.31E-1,1.19E-1,1.08E-1,9.71E-2,8.88E-2,8.01E-2},
      { 9.32E-2,9.42E-2,9.10E-2,8.75E-2,8.00E-2,7.44E-2,6.91E-2},
      { 5.98E-2,6.52E-2,6.72E-2,6.62E-2,6.40E-2,6.12E-2,5.82E-2},
      { 3.83E-2,4.45E-2,4.80E-2,4.96E-2,4.98E-2,4.90E-2,4.77E-2},
      { 2.46E-2,3.01E-2,3.40E-2,3.65E-2,3.79E-2,3.84E-2,3.83E-2},
      { 1.59E-2,2.03E-2,2.39E-2,2.66E-2,2.85E-2,2.97E-2,3.04E-2},
      { 1.04E-2,1.37E-2,1.66E-2,1.92E-2,2.12E-2,2.27E-2,2.37E-2},
      { 4.39E-3,6.26E-3,8.26E-3,9.96E-3,1.15E-2,1.29E-2,1.41E-2},
      { 2.06E-3,3.02E-3,4.24E-3,5.28E-3,6.32E-3,7.32E-3,8.26E-3},
      { 1.21E-3,1.69E-3,2.24E-3,2.85E-3,3.50E-3,4.16E-3,4.82E-3},
      { 8.50E-4,1.10E-3,1.38E-3,1.65E-3,2.03E-3,2.45E-3,2.88E-3},
      { 5.90E-4,7.40E-4,8.50E-4,9.90E-4,1.23E-3,1.49E-3,1.71E-3},
      { 3.90E-4,4.60E-4,5.20E-4,6.30E-4,7.65E-4,9.65E-4,1.12E-3},
      { 2.40E-4,2.70E-4,3.10E-4,3.98E-4,4.97E-4,6.03E-4,7.18E-4},
      { 1.50E-4,1.70E-4,2.15E-4,2.70E-4,3.35E-4,4.35E-4,5.00E-4},
      { 1.00E-4,1.20E-4,1.46E-4,1.90E-4,2.40E-4,2.88E-4,3.43E-4},
      { 0.00   ,0.00   ,1.04E-4,1.41E-4,1.80E-4,2.10E-4,2.50E-4},
      { 0.00   ,0.00   ,8.20E-5,1.06E-4,1.38E-4,1.58E-4,1.85E-4},
      { 0.00   ,0.00   ,5.40E-5,7.00E-5,8.60E-5,1.03E-4,1.20E-4},
      { 0.00   ,0.00   ,4.20E-5,5.40E-5,6.50E-5,7.70E-5,8.80E-5}
    };
  
  double f2_A[25][9]=
    {
      {3.52E+3, 3.27E+2, 9.08E+1, 3.85E+1, 2.00E+1, 1.18E+1, 7.55E+0, 5.16E+0, 3.71E+0},
      { 2.58E+2, 1.63E+2, 7.30E+1, 3.42E+1, 1.85E+1, 1.11E+1, 7.18E+0, 4.96E+0, 3.59E+0},
      { -1.12E+2, 4.84E+0, 3.56E+1, 2.34E+1, 1.45E+1, 9.33E+0, 6.37E+0, 4.51E+0, 3.32E+0},
      { -5.60E+1,-1.12E+1, 9.87E+0, 1.24E+1, 9.59E+0, 7.01E+0, 5.16E+0, 3.83E+0, 2.91E+0},
      { -2.13E+1,-1.22E+1,-2.23E+0, 3.88E+0, 5.15E+0, 4.65E+0, 3.87E+0, 3.12E+0, 2.45E+0},
      { -8.25E+0,-9.58E+0,-5.59E+0,-1.40E+0, 1.76E+0, 2.71E+0, 2.71E+0, 2.35E+0, 1.95E+0},
      { -3.22E+0,-6.12E+0,-5.28E+0,-2.87E+0,-1.92E-1, 1.32E+0, 1.69E+0, 1.74E+0, 1.48E+0},
      { -1.11E+0,-3.40E+0,-4.12E+0,-3.08E+0,-6.30E-1, 3.60E-1, 9.20E-1, 1.03E+0, 1.04E+0},
      { -2.27E-1,-2.00E+0,-2.93E+0,-2.69E+0,-1.48E+0,-3.14E-1, 2.69E-1, 5.28E-1, 6.09E-1},
      { 1.54E-1,-1.09E+0,-2.10E+0,-2.15E+0,-1.47E+0,-6.77E-1,-1.80E-1, 1.08E-1, 2.70E-1},
      { 3.28E-1,-6.30E-1,-1.50E+0,-1.68E+0,-1.34E+0,-8.43E-1,-4.60E-1,-1.85E-1,-4.67E-3},
      { 3.32E-1,-2.06E-1,-7.32E-1,-9.90E-1,-9.42E-1,-8.20E-1,-6.06E-1,-4.51E-1,-3.01E-1},
      { 2.72E-1,-3.34E-2,-3.49E-1,-5.65E-1,-6.03E-1,-5.79E-1,-5.05E-1,-4.31E-1,-3.45E-1},
      { 2.02E-1, 2.80E-2,-1.54E-1,-3.00E-1,-3.59E-1,-3.76E-1,-4.60E-1,-3.40E-1,-3.08E-1},
      { 1.38E-1, 4.84E-2,-5.56E-2,-1.44E-1,-2.04E-1,-2.39E-1,-2.54E-1,-2.49E-1,-2.48E-1},
      { 9.47E-2, 4.86E-2,-1.08E-2,-6.44E-2,-1.02E-1,-1.34E-1,-1.62E-1,-1.79E-1,-1.87E-1},
      { 5.33E-2, 3.71E-2, 1.85E-2, 1.63E-3,-1.69E-2,-3.69E-2,-5.66E-2,-7.78E-2,-9.33E-2},
      { 3.38E-2, 2.40E-2, 1.62E-2, 9.90E-3, 3.76E-3,-4.93E-3,-1.66E-2,-3.05E-2,-4.22E-2},
      { 2.12E-2, 1.56E-2, 1.05E-2, 7.80E-3, 7.92E-3, 6.30E-3, 3.20E-4,-8.50E-3,-1.66E-2},
      { 1.40E-2, 9.20E-3, 5.30E-3, 4.70E-3, 6.31E-3, 8.40E-3, 5.30E-3, 8.80E-4,-3.30E-3},
      { 9.20E-3, 4.70E-3, 1.70E-3, 2.60E-3, 4.49E-3, 6.60E-3, 6.00E-3, 4.70E-3, 2.80E-3},
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   },
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   },
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   },
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00 }
    };

  double f2_B[24][9]=
    {
      {2.75E+0, 1.94E+0, 9.13E-1, 6.06E-1, 4.26E-1, 3.14E-1, 2.40E-1, 1.51E-1, 1.03E-1},
      { 1.94E+0, 1.16E+0, 7.56E-1, 5.26E-1, 3.81E-1, 2.87E-1, 2.23E-1, 1.43E-1, 9.78E-2},
      { 5.85E-1, 5.04E-1, 4.10E-1, 3.30E-1, 2.69E-1, 2.17E-1, 1.78E-1, 1.22E-1, 8.71E-2},
      { 7.83E-2, 2.00E-1, 2.35E-1, 2.19E-1, 1.97E-1, 1.73E-1, 1.48E-1, 1.08E-1, 7.93E-2},
      { -1.82E-1, 1.56E-2, 1.04E-1, 1.36E-1, 1.38E-1, 1.31E-1, 1.19E-1, 9.46E-2, 7.19E-2},
      { -2.71E-1,-1.66E-1,-7.29E-2,-4.74E-3, 3.60E-2, 5.50E-2, 6.28E-2, 5.98E-2, 5.09E-2},
      { -1.87E-1,-1.58E-1,-1.09E-1,-5.80E-2,-2.03E-2, 2.48E-3, 1.99E-2, 3.36E-2, 3.27E-2},
      { -1.01E-1,-1.05E-1,-8.95E-2,-6.63E-2,-3.93E-2,-2.38E-2,-9.22E-3, 8.47E-3, 1.52E-2},
      { -5.19E-2,-6.47E-2,-6.51E-2,-5.62E-2,-4.51E-2,-3.49E-2,-2.45E-2,-8.19E-3, 2.05E-3},
      { -3.68E-2,-4.89E-2,-5.36E-2,-5.06E-2,-4.27E-2,-3.65E-2,-2.80E-2,-1.33E-2,-3.47E-3},
      { -2.33E-2,-3.69E-2,-4.41E-2,-4.38E-2,-3.97E-2,-3.50E-2,-2.88E-2,-1.60E-2,-6.68E-3},
      { -8.76E-3,-2.07E-2,-2.90E-2,-3.17E-2,-3.09E-2,-2.92E-2,-2.63E-2,-1.79E-2,-1.03E-2},
      { -1.20E-3,-1.11E-2,-1.90E-2,-2.20E-2,-2.32E-2,-2.24E-2,-2.10E-2,-1.66E-2,-1.11E-2},
      { 1.72E-3,-4.82E-3,-1.02E-2,-1.42E-2,-1.65E-2,-1.66E-2,-1.60E-2,-1.39E-2,-1.09E-2},
      { 2.68E-3,-1.18E-3,-5.19E-3,-8.30E-5,-1.01E-2,-1.14E-2,-1.16E-2,-1.16E-2,-9.99E-3},
      { 2.81E-3, 8.21E-4,-1.96E-3,-3.99E-3,-5.89E-3,-7.13E-3,-8.15E-3,-9.05E-3,-8.60E-3},
      { 2.61E-3, 1.35E-3,-2.99E-4,-1.79E-3,-3.12E-3,-4.44E-3,-5.61E-3,-7.01E-3,-7.27E-3},
      { 2.06E-3, 1.45E-3, 4.64E-4,-5.97E-4,-1.71E-3,-2.79E-3,-3.84E-3,-5.29E-3,-5.90E-3},
      { 1.07E-3, 9.39E-4, 8.22E-4, 3.58E-4,-1.15E-4,-6.60E-4,-1.18E-3,-2.15E-3,-2.88E-3},
      { 4.97E-4, 5.46E-4, 6.15E-4, 5.56E-4, 3.14E-4, 9.80E-5,-1.30E-4,-5.98E-4,-1.07E-4},
      { 1.85E-4, 3.11E-4, 4.25E-4, 4.08E-4, 3.63E-4, 3.04E-4, 2.24E-4, 2.80E-5,-2.10E-4},
      { 4.80E-5, 1.48E-4, 2.44E-4, 2.80E-4, 3.01E-4, 3.11E-4, 3.13E-4, 2.40E-4, 1.10E-4},
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 1.39E-4, 1.80E-4, 1.80E-4},
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 4.38E-5, 7.30E-5, 8.40E-5}
    };
  
  double f2_C[24][7]=
    {
      {7.36E-2, 4.21E-2, 2.69E-2, 1.83E-2, 1.34E-2, 1.01E-2, 7.88E-3},
      { 5.79E-2, 3.61E-2, 2.34E-2, 1.64E-2, 1.21E-2, 9.26E-3, 7.28E-3},
      { 2.94E-2, 2.17E-2, 1.60E-2, 1.23E-2, 9.49E-3, 7.45E-3, 5.95E-3},
      { 2.30E-3, 7.07E-3, 7.76E-3, 7.02E-3, 6.13E-3, 5.17E-3, 4.34E-3},
      { -7.50E-3,-2.00E-3, 9.93E-4, 2.36E-3, 2.82E-3, 2.86E-3, 2.72E-3},
      { 8.27E-3,-5.37E-3,-2.58E-3,-7.96E-4, 3.75E-4, 9.71E-4, 1.28E-3},
      { -5.79E-3,-5.12E-3,-3.86E-3,-2.46E-3,-1.20E-3,-3.74E-4, 1.74E-4},
      { -3.26E-3,-3.43E-3,-3.26E-3,-2.68E-3,-1.84E-3,-1.12E-3,-4.54E-4},
      { -1.46E-3,-1.49E-3,-2.20E-3,-2.18E-3,-1.85E-3,-1.40E-3,-8.15E-4},
      { -4.29E-4,-9.44E-4,-1.29E-3,-1.50E-3,-1.51E-3,-1.36E-3,-9.57E-4},
      { -3.30E-5,-3.66E-4,-6.78E-4,-9.38E-4,-1.09E-3,-1.09E-3,-9.56E-4},
      { 1.50E-4, 3.10E-5,-1.38E-4,-3.06E-4,-4.67E-4,-5.48E-4,-6.08E-4},
      { 1.00E-4, 8.50E-5, 2.30E-5,-6.60E-5,-1.58E-4,-2.40E-4,-3.05E-4},
      { 5.40E-5, 6.50E-5, 4.90E-5, 1.20E-5,-3.60E-5,-8.90E-5,-1.31E-4},
      { 2.90E-5, 4.30E-5, 4.40E-5, 2.90E-5, 5.10E-6,-2.20E-5,-4.80E-5},
      { 1.40E-5, 2.40E-5, 2.80E-5, 2.60E-5, 1.90E-5, 7.50E-6,-1.10E-5},
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   },
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   },
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   },
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   },
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   },
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   },
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   },
      { 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00   , 0.00 }
    };
  



  //=============================================================================

  //Bestimme, welche Reihen der Tabellen fuer Interpolation benoetigt werden:
  
  if (tau<tau_[0])
    {
      std::cout<<"tau is less than the lowest tabulated value:"<<std::endl;
      std::cout<<"tau = "<<tau<<std::endl;
      std::cout<<"minimum = "<<tau_[0]<<std::endl;
      return;
    }
  else if (tau>tau_[nColumn])
    {
      std::cout<<"tau is greater than the highest tabulated value:"<<std::endl;
      std::cout<<"tau = "<<tau<<std::endl;
      std::cout<<"maximum = "<<tau_[nColumn]<<std::endl;
      return;
    }
  
  
  column_ = 0;
  do
    {
      if(tau>tau_[column_])// TAO IF LOOP NOT GET THE CORRECT COLUNM INTERPOLATION
	{
	  column_ = column_ + 1;
	}
    }while (tau>tau_[column_]);

#ifdef DEBUGMEYER 
  std::cout<<"column=  " << column_  <<std::endl;
  std::cout<<"tau c  " << tau_[column_]  <<std::endl;
  std::cout<<"tau c-1  " << tau_[column_-1]  <<std::endl;
#endif


  //  ! Das Gewicht der Reihe zu groesserem Tau:
  if(column_==0)
    {
      weightCol=1; 
    }
  else
    {  
      weightCol = (tau-tau_[column_-1]) / (tau_[column_]-tau_[column_-1]);
    }
  
  
  //Besorge fuer gegebenes 'thetaSchlange' die interpolierten f1- und f2 -Werte
  //der beiden relevanten Reihen:
  //iColumn = 1 => Reihe mit hoeherem Index
  //iColumn = 2 => Reihe mit kleinerem Index
  
  
  iColumn = 1;
  
  //  5 continue;  
  do{
    
    if (column_<=8)
      {
	//! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//! Werte aus 1. Tabelle: 0.2 <= tau <= 1.8
	
	column = column_;
	//	    std::cout<<"thetaSchlange = "<<thetaSchlange<<std::endl;
	
	if (thetaSchlange<thetaSchlangeA[0])
	  {
	    std::cout<<"thetaSchlange is less than the lowest tabulated value in table 1:"<<std::endl;
	    std::cout<<"thetaSchlange = "<<thetaSchlange<<std::endl;
	    std::cout<<"minimum       = "<<thetaSchlangeA[0]<<std::endl;
	    return;
	  }
	else if (thetaSchlange>thetaSchlangeA[nRowA])
	  {
	    std::cout<<"thetaSchlange is greater than the highest tabulated value in table 1:"<<std::endl;
	    std::cout<<"thetaSchlange = "<<thetaSchlange<<std::endl;
	    std::cout<<"maximum       = "<<thetaSchlangeA[nRowA]<<std::endl;
	    thetaSchlange = -1.;
	    return;
	  }
	
	row = 0;
	do 
	  {
	    if (thetaSchlange>thetaSchlangeA[row])
	      {
		row = row + 1;
	      }
	  }while (thetaSchlange>thetaSchlangeA[row]);
#ifdef DEBUGMEYER 
	std::cout<<"row=  " << row  <<std::endl;
#endif
	
	
	//! Gewicht des Tabellenwertes zu groesseren ThetaSchlange:
	
	if(row==0)
	  {
	    weightRow=1;
	    f1_[iColumn] = weightRow  * f1_A[row][column];
	    f2_[iColumn] = weightRow  * f2_A[row][column];
	  }
	else
	  {
	    weightRow = (thetaSchlange-thetaSchlangeA[row-1]) / 
	      (thetaSchlangeA[row]-thetaSchlangeA[row-1]);
	    
	    f1_[iColumn] = (1.-weightRow) * f1_A[row-1][column] +
	      weightRow  * f1_A[row][column];
	    f2_[iColumn] = (1.-weightRow) * f2_A[row-1][column] +
	      weightRow  * f2_A[row][column];
	  }
	    

      }    
    
    else if (column_>8&&column_<=17)
      {
	//! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//! Werte aus 2. Tabelle: 2.0 <= tau <= 7.0
	
	column = column_ - 9;
	
	if (thetaSchlange<thetaSchlangeB[0]) 
	  {
	    std::cout<< "thetaSchlange is less than the lowest tabulated value in table 2:"<<std::endl;
	    std::cout<< "thetaSchlange = "<<thetaSchlange<<std::endl;
	    std::cout<< "minimum       = "<<thetaSchlangeB[0]<<std::endl;
	    return;
	  }
	else if (thetaSchlange>thetaSchlangeB[nRowB]) 
	  {
	    std::cout<< "thetaSchlange is greater than the highest tabulated value in table 2:";
	    std::cout<< "thetaSchlange = "<<thetaSchlange;
	    std::cout<< "maximum = "<<thetaSchlangeB[nRowB];
	    //	call exit
	    thetaSchlange = -1.;
	    return;
	  }
	
	row = 0;
	do
	  {
	    if(thetaSchlange>thetaSchlangeB[row])
	      {
		row = row + 1;
	      }
	  } while (thetaSchlange>thetaSchlangeB[row]);
#ifdef DEBUGMEYER 
	std::cout<<"row=  " << row  <<std::endl;
#endif
	
	
	//	    ! Gewicht des Tabellenwertes zu groesseren ThetaSchlange:
	if(row==0)
	  {
	    weightRow=1;
	    f1_[iColumn] = weightRow  * f1_B[row][column];
	    f2_[iColumn] = weightRow  * f2_B[row][column];
	  }
	else
	  {
	    weightRow = (thetaSchlange-thetaSchlangeB[row-1]) / 
	      (thetaSchlangeB[row]-thetaSchlangeB[row-1]);
	    
	    f1_[iColumn] = (1.-weightRow) * f1_B[row-1][column] +
	      weightRow  * f1_B[row][column];
	    f2_[iColumn] = (1.-weightRow) * f2_B[row-1][column] +
	      weightRow  * f2_B[row][column];
	  }
      }
 


    else
      {
	//! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//				     ! Werte aus 3. Tabelle: 8.0 <= tau <= 20.
	
	column = column_ - 18;
	
	if (thetaSchlange<thetaSchlangeC[0]) 
	  {
	    std::cout<< "thetaSchlange is less than the lowest tabulated value in table 3:"<<std::endl;
	    std::cout<< "thetaSchlange = "<<thetaSchlange<<std::endl;
	    std::cout<< "minimum       = "<<thetaSchlangeC[0]<<std::endl;
	    return;
	  }
	else if (thetaSchlange>thetaSchlangeC[nRowC]) 
	  {
	    std::cout<< "thetaSchlange is greater than the highest tabulated value in table 3:";
	    std::cout<< "\n thetaSchlange = "<<thetaSchlange<<std::endl;
	    std::cout<< "\n maximum = "<<thetaSchlangeC[nRowC]<<std::endl;
	    thetaSchlange = -1.;
	    return;
	  }
	
	row = 0;
	do 
	  {	
	    if(thetaSchlange>thetaSchlangeC[row])
	      {
		row = row + 1;
	      }
	  }while (thetaSchlange>thetaSchlangeC[row]);
#ifdef DEBUGMEYER 
	std::cout<<"row=  " << row  <<std::endl;
#endif
	
	
	//	    ! Gewicht des Tabellenwertes zu groesseren ThetaSchlange:
	
	if(row==0)
	  {
	    weightRow=1;
	    f1_[iColumn] = weightRow  * f1_C[row][column];
	    f2_[iColumn] = weightRow  * f2_C[row][column];
	  }
	else
	  {
	    weightRow = (thetaSchlange-thetaSchlangeC[row-1]) / 
	      (thetaSchlangeC[row]-thetaSchlangeC[row-1]);
	    
	    f1_[iColumn] = (1.-weightRow) * f1_C[row-1][column] +
	      weightRow  * f1_C[row][column];
	    f2_[iColumn] = (1.-weightRow) * f2_C[row-1][column] +
	      weightRow  * f2_C[row][column];
	  }
	
      }
    
    
    
#ifdef DEBUGMEYER  
    std::cout<<"f1_[iColumn]=  " << f1_[iColumn]  <<std::endl;
    std::cout<<"f2_[iColumn]=  " << f2_[iColumn]  <<std::endl;
    
    std::cout<<"wc: "<<weightCol<<std::endl;
    std::cout<<"wr: "<<weightRow<<std::endl;
    std::cout<<"icol: "<<iColumn<<std::endl;
#endif
    
    iColumn++ ;
    
  }while(iColumn<=2);	      
  
  
#ifdef DEBUGMEYER  
  std::cout<<"f1: "<<*f1<<std::endl;
  std::cout<<"f2: "<<*f2<<std::endl;
#endif
  
    
  *f1 = weightCol*f1_[1] + (1.-weightCol)*f1_[2];
  *f2 = weightCol*f2_[1] + (1.-weightCol)*f2_[2];
  
}
 

//========================================================================================

  



/*-
	options /extend_source

	subroutine throwMeyerAngle (theta)
c	==================================

	implicit none

	real lowerbound,y1,y2,f,root,radiant,fraction
	integer bin,nBin
	integer nBinMax
	parameter (nBinMax=201)

	real theta,thetaStep
	real value(0:nBinMax)  /0.,nBinMax*0./
	real area(nBinMax)     /   nBinMax*0./
	real integ(0:nBinMax)  /0.,nBinMax*0./
	common /MeyerTable/ value,area,integ,thetaStep,nBin

	real rHelp

	real random
	integer seed
	common /seed/ seed


c bin: Nummer des Bins, innerhalb dessen das Integral den Wert von 
c random erreicht oder ueberschreitet:

	random = ran(seed)

	bin = 1
	do while (random.GT.integ(bin))
	    bin = bin + 1
	    if (bin.GT.nBin) then
		write(*,*) 'error 1'
		call exit
	    endif
	enddo

	fraction = (random-integ(bin-1)) / (integ(bin)-integ(bin-1))
	y1 = value(bin-1)
	y2 = value(bin)
	f = thetaStep / (y2-y1)
	rHelp = y1*f

	radiant = rHelp*rHelp + fraction*thetaStep*(y1+y2)*f
	root = SQRT(radiant)
	lowerBound = real(bin-1)*thetaStep
	if (f.GT.0) then
	    theta = lowerBound - rHelp + root
	else
	    theta = lowerBound - rHelp - root
	endif


	END


c===============================================================================

	options /extend_source

*/

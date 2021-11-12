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

#include "musrParameters.hh"
// #include "musrErrorMessage.hh"  - musrErrorMessage class can not be used inside "musrParameters" constructor, because
//                                   musrErrorMessage is crated later!

#include "CLHEP/Units/SystemOfUnits.h"

musrParameters::musrParameters(G4String steeringFileName) 
{
  pointerToParameters = this;
  boolG4RegionRequested = false;
  mySteeringFileName = steeringFileName;
  myStopFileName = mySteeringFileName;
  myStopFileName.replace(myStopFileName.length()-3,myStopFileName.length()-1,"stop");
  myRandomNumberFileName = mySteeringFileName;
  myRandomNumberFileName.replace(myRandomNumberFileName.length()-3,myRandomNumberFileName.length()-1,"rndm");

  // Read in the parameters, which have to be known before the detector construction is run
  // (and therefore the parameters can not be read in in the musrDetectorConstruction.cc class).

  FILE *fSteeringFile=fopen(steeringFileName.c_str(),"r");
  if (fSteeringFile==NULL) {
    G4cout<<"musrParameters::musrParameters:  steeringFileName=\""<<steeringFileName
	  <<"\" not opened for some reason."<<G4endl;
    G4cout << "S T O P    F O R C E D" << G4endl;
    exit(1);
  }
  G4cout<<"musrParameters::musrParameters:  steeringFileName=\""<<steeringFileName<<"\" opened."<<G4endl;
  char  line[501];

  while (!feof(fSteeringFile)) {
    fgets(line,500,fSteeringFile);
    if ((line[0]!='#')&&(line[0]!='\n')&&(line[0]!='\r')) {
      char tmpString0[100]="Unset";
      sscanf(&line[0],"%s",tmpString0);
      // First find out how many events will be generated (needs to be known at an early stage, if the
      // field is to be set in steps):
      if (strcmp(tmpString0,"/run/beamOn")==0) {
	int nev;
	sscanf(&line[0],"%*s %d", &nev);
	musrParameters::nrOfEventsToBeGenerated = nev;
      }
      if (strncmp(tmpString0,"/gps/",5)==0) {
	musrParameters::boolG4GeneralParticleSource = true;
	G4cout<<"\n========================================================================"<<G4endl;
	G4cout<<"musrParameters.cc:  GPS (General Particle Source) requested in the macro."<<G4endl;
	G4cout<<"                    GPS will be used instead of the primary generator action."<<G4endl;
      }

      // Now find some private parameters that need to be initialised at early stage
      if ( (strcmp(tmpString0,"/musr/ignore")!=0)&&(strcmp(tmpString0,"/musr/command")!=0) ) continue;

      char tmpString1[100]="Unset", tmpString2[100]="Unset", tmpString3[100]="Unset";
      sscanf(&line[0],"%*s %s %s %s",tmpString1,tmpString2,tmpString3);
      //      if (strcmp(tmpString1,"G4GeneralParticleSource")==0){
      //      	if (strcmp(tmpString2,"true")==0){ musrParameters::boolG4GeneralParticleSource = true; }
      //      }
      if (strcmp(tmpString1,"G4OpticalPhotons")==0){
      	if (strcmp(tmpString2,"true")==0){ musrParameters::boolG4OpticalPhotons = true; }
      }
      if (strcmp(tmpString1,"G4OpticalPhotonsUnprocess")==0){
      	if (strcmp(tmpString2,"true")==0){ musrParameters::boolG4OpticalPhotonsUnprocess = true; }
      }
      if (strcmp(tmpString1,"region")==0) {
	boolG4RegionRequested = true;
      } 
      if ( (strcmp(tmpString1,"construct")==0) && boolG4RegionRequested) {
	G4cout<<"musrParameters.cc: User requests to construct a detector volume "<<tmpString3<<G4endl;
	G4cout<<"                   after a previous G4Region definition."<<G4endl;
	G4cout<<"                   Perhaps not a crutial problem, but the execution will be stopped"<<G4endl;
	G4cout<<"                   anyway just to be on the safe side.  Please correct "<<steeringFileName<<" file"<<G4endl;
	G4cout<<"                   such that G4Region is called only after all detector volumes have been created."<<G4endl;
	G4cout<<"                    S T O P     F O R C E D!"<<G4endl;
	exit(1);
      }
    }
  }
  fclose(fSteeringFile);
}

musrParameters::~musrParameters() {}

musrParameters* musrParameters::pointerToParameters=NULL;
musrParameters* musrParameters::GetInstance() {
  return pointerToParameters;
}

G4String musrParameters::mySteeringFileName="Unset";
G4String musrParameters::myStopFileName="Unsetblablabla034tdk40928jfmfnakfh921djf02UNSET";
G4String musrParameters::myRandomNumberFileName="Unsetblablabla034tdk40928jfmfnakfh921djf02UNSED";
G4bool   musrParameters::storeOnlyEventsWithHits=true;
G4int    musrParameters::storeOnlyEventsWithHitInDetID=0;
G4double musrParameters::signalSeparationTime=100*CLHEP::nanosecond;
G4bool   musrParameters::storeOnlyTheFirstTimeHit=false;
G4int    musrParameters::storeOnlyEventsWithMuonsDecayedInDetID=0;
G4bool   musrParameters::field_DecayWithSpin=false;
G4bool   musrParameters::killAllElectrons=false;
G4bool   musrParameters::killAllPositrons=false;
G4bool   musrParameters::killAllGammas=false;
G4bool   musrParameters::killAllNeutrinos=true;
G4bool   musrParameters::boolG4GeneralParticleSource=false;
G4bool   musrParameters::boolG4OpticalPhotons=false;
G4bool   musrParameters::boolG4OpticalPhotonsUnprocess=false;
//cks G4bool   musrParameters::includeMuoniumProcesses =true; // TS
//G4bool   musrParameters::boolG4GeneralParticleSource=true;
G4int    musrParameters::nrOfEventsToBeGenerated=0;
G4bool   musrParameters::finiteRiseTimeInScintillator=false;
G4int    musrParameters::maximumTimePerEvent=60;
G4int    musrParameters::maximumNrOfStepsPerTrack=100000;

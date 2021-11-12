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

#include "musrScintHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include <fstream>
#include <iostream>
#include <iomanip>

G4Allocator<musrScintHit> musrScintHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrScintHit::musrScintHit() {
}

G4int musrScintHit::ScintMultihit=0;
G4int musrScintHit::runIDoldScint=-1;
G4int musrScintHit::eventIDoldScint=-1;
G4int musrScintHit::NIS=0;
G4int musrScintHit::ScintChamberNbold=-1;
G4int musrScintHit::verboseLevel=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrScintHit::~musrScintHit() {
  //save the Tree header. The file will be automatically closed
  //when going out of the function scope
  //  rootTree->Write();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrScintHit::musrScintHit(const musrScintHit& right)
  : G4VHit()
{
  particleName = right.particleName;
  trackID   = right.trackID;
  edep      = right.edep;
  pre_pos       = right.pre_pos;
  pol       = right.pol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const musrScintHit& musrScintHit::operator=(const musrScintHit& right)
{
  particleName = right.particleName;
  trackID   = right.trackID;
  edep      = right.edep;
  pre_pos       = right.pre_pos;
  pol       = right.pol;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int musrScintHit::operator==(const musrScintHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrScintHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pre_pos);
    circle.SetScreenSize(0.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrScintHit::Print()
{
  if (verboseLevel>2) G4cout<<"VERBOSE 3: Kamil:  musrScintHit::Print()"<<G4endl;
  
  G4RunManager* fRunManager = G4RunManager::GetRunManager();
  eventID = fRunManager->GetCurrentEvent()->GetEventID();
  runID   = fRunManager->GetCurrentRun()->GetRunID();
  G4int ScintMultihitSwitch=0;
  
  if (runID != runIDoldScint)  {
      NIS=0;
      ScintMultihit = 0;
      eventIDoldScint = -1;
      ScintChamberNbold = -1;
    }
  
  //cks  if (particleName== "e+" and (eventIDoldScint != eventID or ScintChamberNbold != IBchamberNb)) {
  if (particleName== "e+") {
    if (eventIDoldScint == eventID) {
      ScintMultihit++;
      ScintMultihitSwitch=1;
    }
    NIS++;

    G4FieldManager *fMgr=G4TransportationManager::GetTransportationManager()->GetFieldManager();
    point[0]=0.;
    point[1]=0.;
    point[2]=0.;
    B[2]=0.0;
    if(!fMgr->DoesFieldChangeEnergy())  {          //then we have a magnetic field
      mfield = fMgr->GetDetectorField();
      mfield->GetFieldValue(point,B);
      B[0]=B[0]/CLHEP::tesla;
      B[1]=B[1]/CLHEP::tesla;
      B[2]=B[2]/CLHEP::tesla;
    }
    //      G4cout << "  Segment: " << IBchamberNb << G4endl;
    //      G4cout <<"Position " << pos.x()/cm<<" "<<pos.y()/cm <<" "<< pos.z()/cm <<G4endl;
    //	 G4cout << "Field is "<< B[2]<<G4endl;
    
    std::ofstream posfile1;
    posfile1.open ("scint.dat", std::ios::out | std::ios::app);
    posfile1 << runID << " " << eventID 
	   << " " << logicalVolume
	   << " " << ScintMultihitSwitch 
	   <<" " << edep
	   << " " << fabs(B[2]) 
           <<" "<< pre_pos.x()/CLHEP::cm<<" "<<pre_pos.y()/CLHEP::cm <<" "<< pre_pos.z()/CLHEP::cm  
	   << " " << globalTime/CLHEP::s 
        // << " " << IBchamberNb 
      //	     << " first=" << firstStepInVolume << " last=" << lastStepInVolume
	   << G4endl;
    posfile1.close();
  }
  eventIDoldScint=eventID;
  runIDoldScint = runID;
  //  ScintChamberNbold = IBchamberNb; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


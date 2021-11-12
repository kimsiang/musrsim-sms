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

#include "musrSteppingVerbose.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrSteppingVerbose::musrSteppingVerbose()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrSteppingVerbose::~musrSteppingVerbose()
{}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrSteppingVerbose::StepInfo()
{
  CopyState();
  
  G4int prec = G4cout.precision(3);

  if( verboseLevel >= 1 ){
    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){
      G4cout << G4endl;    
      G4cout << std::setw( 5) << "#Step#"     << " "
	     << std::setw( 6) << "X"          << "    "
	     << std::setw( 6) << "Y"          << "    "  
	     << std::setw( 6) << "Z"          << "    "
	     << std::setw( 9) << "KineE"      << " "
	     << std::setw( 9) << "dEStep"     << " "  
	     << std::setw(10) << "StepLeng"     
	     << std::setw(10) << "TrakLeng" 
	     << std::setw(10) << "Volume"    << "  "
	     << std::setw(10) << "Process"   << G4endl;	          
    }

    G4cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
	<< std::setw(6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
      //	<< std::setw(6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
	<< std::setw(6) << G4BestUnit(fStep->GetDeltaEnergy(),"Energy")
	<< std::setw(6) << G4BestUnit(fStep->GetStepLength(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetTrackLength(),"Length")
	<< "  ";

    // if( fStepStatus != fWorldBoundary){ 
    if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << std::setw(10) << fTrack->GetVolume()->GetName();
    } else {
      G4cout << std::setw(10) << "OutOfWorld";
    }

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << "  " 
        << std::setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	                                ->GetProcessName();
    } else {
      G4cout << "   UserLimit";
    }

    G4cout << G4endl;

    if( verboseLevel == 2 ){
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
	                    fN2ndariesAlongStepDoIt +
	                    fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0){
	G4cout << "    :----- List of 2ndaries - "
	       << "#SpawnInStep=" << std::setw(3) << tN2ndariesTot 
	       << "(Rest="  << std::setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << std::setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << std::setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << std::setw(3) << (*fSecondary).size()
	       << " ---------------"
	       << G4endl;

	for(size_t lp1=(*fSecondary).size()-tN2ndariesTot; 
                        lp1<(*fSecondary).size(); lp1++){
	  G4cout << "    : "
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
		 << std::setw(10)
		 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
	  G4cout << G4endl;
	}
              
	G4cout << "    :-----------------------------"
	       << "----------------------------------"
	       << "-- EndOf2ndaries Info ---------------"
	       << G4endl;
      }
    }
    
  }
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrSteppingVerbose::TrackingStarted()
{

  CopyState();
G4int prec = G4cout.precision(3);
  if( verboseLevel > 0 ){

    G4cout << std::setw( 5) << "Step#"      << " "
           << std::setw( 6) << "X"          << "    "
	   << std::setw( 6) << "Y"          << "    "  
	   << std::setw( 6) << "Z"          << "    "
	   << std::setw( 9) << "KineE"      << " "
	   << std::setw( 9) << "dEStep"     << " "  
	   << std::setw(10) << "StepLeng"  
	   << std::setw(10) << "TrakLeng"
	   << std::setw(10) << "Volume"     << "  "
	   << std::setw(10) << "Process"    << G4endl;	     

    G4cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
	<< std::setw(6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
	<< std::setw(6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
	<< std::setw(6) << G4BestUnit(fStep->GetStepLength(),"Length")
	<< std::setw(6) << G4BestUnit(fTrack->GetTrackLength(),"Length")
	<< "  ";

    if(fTrack->GetNextVolume()){
      G4cout << std::setw(10) << fTrack->GetVolume()->GetName();
    } else {
      G4cout << std::setw(10) << "OutOfWorld";
    }
    G4cout  << "    initStep" << G4endl;
  }
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

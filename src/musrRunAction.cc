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

// Make G4Timer appear first!
#include "G4Timer.hh"
#include "musrRunAction.hh"
#include "musrEventAction.hh"
#include "musrSteppingAction.hh"
#include "musrScintSD.hh"
#include "G4Run.hh"
#include "musrErrorMessage.hh"
#include "F04GlobalField.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrRunAction::musrRunAction() {
  timer = new G4Timer;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrRunAction::~musrRunAction() {
  delete timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrRunAction::BeginOfRunAction(const G4Run* aRun)  {
  timer->Start();
  G4int run_id= aRun->GetRunID();
  //  if (run_id%100 == 0) {
  G4cout << "### Run " << run_id << G4endl;
  musrRootOutput::GetRootInstance()->StoreGeantParameter(6,run_id);
  musrRootOutput::GetRootInstance()->BeginOfRunAction();

  // Initiate global electromagnetic field (if it exists):
  if (F04GlobalField::Exists()) {
    FieldList* fields = F04GlobalField::getObject()->getFields();

    if (fields) {
      if (fields->size()>0) {
	G4int jjj=0;
	G4cout<<"\n------------ The following fields were defined:  ----------------"<<G4endl;
	FieldList::iterator i;
	for (i=fields->begin(); i!=fields->end(); ++i) {
	  (*i)->construct();
	  //	  G4cout<<"\t==> "<<(*i)->GetElementFieldName()<<G4endl;
	  // Get the nominal field value for the given field and store it in the Root output 
	  G4double nomFieldValue = (*i)->GetNominalFieldValue();
	  musrRootOutput::GetRootInstance()->SetFieldNomVal(jjj,nomFieldValue);
	  jjj++;
	}
	G4cout<<"-----------------------------------------------------------------"<<G4endl;
      }
    }

    // Print out the field values at the points user requested to be printed out:
    F04GlobalField::getObject()->PrintFieldAtRequestedPoints();
  }
  
  // initialise the "boolIsVvvInfoRequested" variable in the musrSteppingAction.cc and musrScintSD.cc
  G4bool boolIsVvvInfoRequested = false;
  if ((musrRootOutput::store_det_VvvKine)||(musrRootOutput::store_det_VvvX)||
      (musrRootOutput::store_det_VvvY)||(musrRootOutput::store_det_VvvZ)||
      (musrRootOutput::store_det_VvvVolID)||(musrRootOutput::store_det_VvvProcID)||
      (musrRootOutput::store_det_VvvTrackID)||(musrRootOutput::store_det_VvvParticleID)) {
    boolIsVvvInfoRequested = true;
  }
  musrSteppingAction::GetInstance()->SetVvvInfoRequested(boolIsVvvInfoRequested);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrRunAction::EndOfRunAction(const G4Run* aRun)  { 
  musrScintSD::GetInstance()->EndOfRun();

  musrRootOutput::GetRootInstance()->StoreGeantParameter(5,aRun->GetNumberOfEvent());
  musrRootOutput::GetRootInstance()->EndOfRunAction();

  if (F04GlobalField::Exists()) {
    F04GlobalField* myGlobalField = F04GlobalField::getObject();   //Was causing seg. fault when called here
    if (myGlobalField!=NULL)  {delete myGlobalField;}
  }

  musrErrorMessage::GetInstance()->PrintErrorSummary();

  timer->Stop();
  G4cout << "musrRunAction::EndOfRunAction:"<<G4endl;
  G4cout << "    Number of events    = " << aRun->GetNumberOfEvent()<<G4endl;
  //         << " " << *timer << G4endl;
  G4cout << "    User elapsed time   = "<<timer->GetUserElapsed()/3600<<"h   = "
	 <<timer->GetUserElapsed()/60<<"min   = "<<timer->GetUserElapsed()<<"s."<<G4endl;
  G4cout << "    Real elapsed time   = "<<timer->GetRealElapsed()/3600<<"h   = "
	 <<timer->GetRealElapsed()/60<<"min   = "<<timer->GetRealElapsed()<<"s."<<G4endl;
  G4cout << "    System elapsed time = "<<timer->GetSystemElapsed()/3600<<"h   = "
	 <<timer->GetSystemElapsed()/60<<"min   = "<<timer->GetSystemElapsed()<<"s."<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




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

#include "musrEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "musrParameters.hh"
#include "musrRootOutput.hh"
#include "musrErrorMessage.hh"
#include "musrSteppingAction.hh"
#include "F04GlobalField.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include <sys/stat.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int  musrEventAction::nHowOftenToPrintEvent=10000;
G4double  musrEventAction::maximumRunTimeAllowed=85000;

musrEventAction::musrEventAction() {
  pointer=this;
  timeOfRunStart = -1000;
}
musrEventAction* musrEventAction::pointer=0;
musrEventAction* musrEventAction::GetInstance() {
  return pointer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
musrEventAction::~musrEventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrEventAction::BeginOfEventAction(const G4Event* evt) {
  //  G4cout<<"musrEventAction::BeginOfEventAction:"<<G4endl;
  long thisEventNr = (long) (evt->GetEventID());

  //  if (thisEventNr == 44654) {trackingManager->SetVerboseLevel(2);}
  musrSteppingAction::GetInstance()->DoAtTheBeginningOfEvent();

  if (F04GlobalField::Exists()) {
    F04GlobalField::getObject() -> CheckWhetherAnyNominalFieldValueNeedsToBeChanged(thisEventNr);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void musrEventAction::EndOfEventAction(const G4Event* evt)  {
  //  cout << ":." << flush;
  long thisEventNr = (long) evt->GetEventID();

  //  musrSteppingAction::GetInstance()->DoAtTheEndOfEvent();
  
  // write out the root tree:
  musrRootOutput* myRootOutput = musrRootOutput::GetRootInstance();
  G4RunManager* fRunManager = G4RunManager::GetRunManager();
  myRootOutput->SetRunID(fRunManager->GetCurrentRun()->GetRunID());
  myRootOutput->SetEventID(fRunManager->GetCurrentEvent()->GetEventID());
  myRootOutput->FillEvent();

  // get number of stored trajectories
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  //  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  
  // periodic printing
  //
  //  if (thisEventNr != 0 and thisEventNr%10000 == 0) {
  if (timeOfRunStart == -1000) timeOfRunStart=time(0);   // initiate the time when execution starts once at the beginning of run
  if ((time(0)-timeOfRunStart)>maximumRunTimeAllowed) {
    // Stop the execution of the run - the run took already too long time
    char eMessage[200];
    sprintf(eMessage,"musrEventAction::EndOfEventAction(): Run execution exceeded the allowed maximum time (maximum = %f sec) ==> RUN STOPPED",maximumRunTimeAllowed);
    musrErrorMessage::GetInstance()->musrError(WARNING,eMessage,false);
    G4RunManager* fRunManager = G4RunManager::GetRunManager();
    fRunManager->AbortRun(true);
  }

  // Stop the run, if the user wished that (by creating a file "RUNNUMBER.stop"
  struct stat stFileInfo; 
  if (stat((musrParameters::myStopFileName).c_str(),&stFileInfo)==0) {
    // Stop the execution of the run - the file "RUNNUMBER.stop" was found
    remove((musrParameters::myStopFileName).c_str());
    char eMessage[200];
    sprintf(eMessage,"musrEventAction::EndOfEventAction(): Request to stop the run was found (%s) ==> RUN STOPPED",musrParameters::myStopFileName.c_str());
    musrErrorMessage::GetInstance()->musrError(WARNING,eMessage,false);
    G4RunManager* fRunManager = G4RunManager::GetRunManager();
    fRunManager->AbortRun(true);
  }

  if ((thisEventNr != 0) && (thisEventNr%nHowOftenToPrintEvent == 0)) {
    time_t curr=time(0);
    //char * ctime(const time_t * tp);
    G4cout << ">>> Event " << evt->GetEventID() <<". Running already for "<<curr-timeOfRunStart<<" seconds.  Present time: "<< ctime(&curr);
    G4cout.flush();
    //    G4cout << "                   seed set to "<< CLHEP::HepRandom::getTheSeed();//<< G4endl; 
  }
    
  // extract the trajectories and draw them
  if (G4VVisManager::GetConcreteInstance())  {
    for (G4int i=0; i<n_trajectories; i++) 
      { G4Trajectory* trj = (G4Trajectory*)
	  ((*(evt->GetTrajectoryContainer()))[i]);
	trj->DrawTrajectory(); //removed argument 1000 (JSL)
      }
  }
}

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

#include "musrDetectorMessenger.hh"
#include "musrDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

#include "G4RunManager.hh"   //cks  included in order to be able to change run ID
#include "Randomize.hh"      //cks  included in order to initialise the random nr. generator by time
#include <time.h>            //cks   -----------------------------||-------------------------------
#include "musrEventAction.hh" // cks needed for setting the variable "nHowOftenToPrintEvent"
#include "musrPrimaryGeneratorAction.hh" // cks needed for the initialisation of the random nr. generator by event nr.
#include "musrParameters.hh"
#include <vector>
#include "globals.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


musrDetectorMessenger::musrDetectorMessenger(musrDetectorConstruction* myDet)
  :myDetector(myDet)
{ 
  musrDir = new G4UIdirectory("/musr/");
  musrDir->SetGuidance("UI commands specific to this example.");

  Ignore1Cmd = new G4UIcmdWithAString("/musr/ignore",this);
  Ignore1Cmd->SetGuidance("This command is ignored by the messenger, but used for the detector construction.");
  Ignore1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  Ignore2Cmd = new G4UIcmdWithAString("/musr/command",this);
  Ignore2Cmd->SetGuidance("This command is ignored by the messenger, but used for the detector construction.");
  Ignore2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  runDir = new G4UIdirectory("/musr/run/");
  runDir->SetGuidance("musr run control");

  RunIDSetCmd = new G4UIcmdWithAnInteger("/musr/run/runID",this);
  RunIDSetCmd->SetGuidance("Set the run number");
  RunIDSetCmd->SetParameterName("something",false);
  RunIDSetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RandomOptionCmd = new G4UIcmdWithAnInteger("/musr/run/randomOption",this);
  RandomOptionCmd->SetGuidance("Specify the random number generator initialisation");
  RandomOptionCmd->SetGuidance("   0 ... no initialisation (default)");
  RandomOptionCmd->SetGuidance("   1 ... use actual computer time to initialise now");
  RandomOptionCmd->SetGuidance("   2 ... use event number to initialise at the beginning of each event");
  RandomOptionCmd->SetGuidance("   3 ... read in the random no. initial values for each event from a file");
  RandomOptionCmd->SetGuidance("   4 ... read in the random no. initial values for each event from the file kamil.rndm");
  RandomOptionCmd->SetParameterName("randomOpt",false);
  RandomOptionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  HowOftenToPrintEventCmd = new G4UIcmdWithAnInteger("/musr/run/howOftenToPrintEvent",this);
  HowOftenToPrintEventCmd->SetGuidance("Each n-th event will be notified.  Set _n_ by this command.");
  HowOftenToPrintEventCmd->SetParameterName("howOftenToPrintEv",false);
  HowOftenToPrintEventCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RndmEventToSaveSeedsCmd = new G4UIcmdWithAnInteger("/musr/run/rndmEventToSaveSeeds",this);
  RndmEventToSaveSeedsCmd -> SetGuidance("Save seeds of the random number generators of the given event number.");
  RndmEventToSaveSeedsCmd -> SetParameterName("rndmEventToSaveSe",false);
  RndmEventToSaveSeedsCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);

  detDir = new G4UIdirectory("/musr/det/");
  detDir->SetGuidance("detector control.");

  UpdateCmd = new G4UIcmdWithoutParameter("/musr/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrDetectorMessenger::~musrDetectorMessenger()
{
  delete UpdateCmd;
  delete detDir;
  delete musrDir;
  delete Ignore1Cmd;
  delete Ignore2Cmd;
  delete RunIDSetCmd;
  delete RandomOptionCmd;
  delete HowOftenToPrintEventCmd;
  delete RndmEventToSaveSeedsCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)  { 

  if( command == UpdateCmd )
    { myDetector->UpdateGeometry(); }

  if( command == RunIDSetCmd )
    { (G4RunManager::GetRunManager())->SetRunIDCounter(RunIDSetCmd->GetNewIntValue(newValue));}
  
  if( command == RandomOptionCmd )
    { 
      G4int RandomOption=RandomOptionCmd->GetNewIntValue(newValue);
      if (RandomOption == 1) {
	//	G4long seed=time(0); //returns time in seconds as an integer
	//	HepRandom::setTheSeed(seed);//changes the seed of the random engine
	G4cout << "******************************************" << G4endl;
	G4cout << "*** Random Seed set by the system time ***" << G4endl;
	G4cout << "******************************************" << G4endl;
	long seeds[2];
	time_t systime = time(NULL);
	seeds[0] = (long) systime;
	seeds[1] = (long) (systime*G4UniformRand());
	G4cout << "seed1: " << seeds[0] << "; seed2: " << seeds[1] << G4endl;
	CLHEP::HepRandom::setTheSeeds(seeds);
	CLHEP::HepRandom::showEngineStatus();
      }
      else if (RandomOption == 2) {
	G4cout << "*******************************************" << G4endl;
	G4cout << "*** Random Seed set by the event number ***" << G4endl;
	G4cout << "*******************************************" << G4endl;
	//	musrEventAction::setRandomNrSeedAccordingEventNr=1;
	musrPrimaryGeneratorAction::setRandomNrSeedAccordingEventNr=1;
	//	musrEventAction::setMyEventNr(70);
      }
      else if (RandomOption == 3) {
	G4cout << "*******************************************" << G4endl;
	G4cout << "*** Random Seed set from external file  ***" << G4endl;
	G4cout << "*******************************************" << G4endl;
	//	musrEventAction::setRandomNrSeedFromFile=1;
	musrPrimaryGeneratorAction::setRandomNrSeedFromFile=1;
	std::ifstream indata;
	int num;

	//	indata.open("randomNum.dat"); // opens the file
	indata.open(musrParameters::myRandomNumberFileName);
	if(!indata) { // file couldn't be opened
	  G4cout << "Error: file could not be opened" << G4endl;
	  exit(1);
	}
	std::vector<int> * seedVector = musrPrimaryGeneratorAction::GetPointerToSeedVector();
	std::vector<int> tmpSeedVector;
	indata >> num;
	while ( !indata.eof() ) { // keep reading until end-of-file
	  //	  G4cout << "The next number is " << num << G4endl;
	  tmpSeedVector.push_back(num);
	  //	  musrPrimaryGeneratorAction::lastEventID_in_pointerToSeedVector = num;
	  indata >> num; // sets EOF flag if no value found
	}
	indata.close();
	// save last eventID for further use in musrPrimaryGeneratorAction.cc
	musrPrimaryGeneratorAction::lastEventID_in_pointerToSeedVector = tmpSeedVector.back();
	// save the eventIDs in the reversed order so that they can be popped-back (using pop_back)
	std::vector<int>::reverse_iterator rit;
	for ( rit=tmpSeedVector.rbegin() ; rit < tmpSeedVector.rend(); ++rit ) {
	  seedVector->push_back(*rit);
	}
	tmpSeedVector.clear();
      }
      else if (RandomOption == 4) {
	G4cout << "*********************************************" << G4endl;
	G4cout << "*** Random Seed set from kamil.rndm file  ***" << G4endl;
	G4cout << "*********************************************" << G4endl;
	musrPrimaryGeneratorAction::setRandomNrSeedFromFile_RNDM=1;
      }
    }
  if ( command == HowOftenToPrintEventCmd ) 
    {
      G4int n = HowOftenToPrintEventCmd->GetNewIntValue(newValue);
      musrEventAction::nHowOftenToPrintEvent=n;
    }
  if ( command == RndmEventToSaveSeedsCmd ) 
    {
      G4int n = RndmEventToSaveSeedsCmd->GetNewIntValue(newValue);
      musrPrimaryGeneratorAction::nRndmEventToSaveSeeds=n;
    }
  
}

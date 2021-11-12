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

#ifndef musrSteppingAction_h
#define musrSteppingAction_h 1
#include "G4UserSteppingAction.hh"
#include "G4ProcessManager.hh"
#include "globals.hh"
#include "musrRootOutput.hh"
#include <fstream>

class G4Timer;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class musrSteppingAction : public G4UserSteppingAction
{
  public:

    static musrSteppingAction* GetInstance();
    musrSteppingAction();
    ~musrSteppingAction();

    void UserSteppingAction(const G4Step *theStep);
    void DoAtTheBeginningOfEvent();
    void SetLogicalVolumeAsSpecialSaveVolume(G4String logicName, G4int volumeID);
    void SetVolumeForMuonEventReweighting(G4String logicName, G4int weight);
    G4bool  GetInfoAboutOldTrack(G4int trackID, G4int& parentTrackID, G4int& particleID, G4double& vertexKine,
				 G4ThreeVector& vertexPosition, G4String& vertexLogVol, G4String& vertexProcess);
    G4bool  AreTracksCommingFromSameParent(G4int trackID1, G4int trackID2, G4String volumeName);
    static const G4int maxNumberOfOldTracks=200;
    G4bool IsVvvInfoRequested() {return boolIsVvvInfoRequested;}
    void   SetVvvInfoRequested(G4bool boolvar) {boolIsVvvInfoRequested = boolvar;}
    void   SetCalculationOfFieldIntegralRequested(G4bool decision) {boolCalculateFieldIntegral = decision;}

  private:
    G4Timer* timer;
    time_t realTimeWhenThisEventStarted;
    static musrSteppingAction* pointer;
    musrRootOutput* myRootOutput;

    G4bool muAlreadyWasInTargetInThisEvent;
    G4bool muAlreadyWasInM0InThisEvent;
    G4bool muAlreadyWasInM1InThisEvent;
    G4bool muAlreadyWasInM2InThisEvent;
    G4bool radioactiveElectronAlreadySavedInThisEvent;
    G4bool boolIsAnySpecialSaveVolumeDefined;
    G4bool boolIsVvvInfoRequested;
    std::map<G4String,G4int>  saveVolumeMapping;
    G4String lastActualVolume;
    G4bool boolMuonEventReweighting;
    G4bool boolCalculateFieldIntegral;
    std::map<G4String,G4int>  volumeMuonWeightMapping;

    G4int         indexOfOldTrack;
    std::map<G4int,G4int> myOldTracksMap;
    G4int         particleID_oldTrack[maxNumberOfOldTracks];
    G4int         parentTrackID_oldTrack[maxNumberOfOldTracks];
    G4double      vertexKine_oldTrack[maxNumberOfOldTracks];
    G4ThreeVector vertexPosition_oldTrack[maxNumberOfOldTracks];
    G4String      vertexLogVol_oldTrack[maxNumberOfOldTracks];
    G4String      vertexProcess_oldTrack[maxNumberOfOldTracks];

    // for field integral along the muon path (if its calculation is requested by the user)
    G4double BxIntegral, ByIntegral, BzIntegral;
    G4double BzIntegral1, BzIntegral2, BzIntegral3;
    G4double CoordinateForFieldIntegral[4];
    G4double FieldForFieldIntegral[6];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

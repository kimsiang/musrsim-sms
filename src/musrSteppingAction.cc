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

#include "musrSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"  // needed for the event nr. comparison
#include "G4Run.hh"         // ---------------||------------------
#include "G4MagneticField.hh"          // needed for storing the magnetic field to the Root class
#include "G4FieldManager.hh"           // ---------------||------------------
#include "G4TransportationManager.hh"  // ---------------||------------------
//#include "G4OpBoundaryProcess.hh"      // Optical photon process
#include "musrErrorMessage.hh"
#include "musrParameters.hh"
#include "F04GlobalField.hh"
//#include "TCanvas.h"
//#include "TMath.h"
//#include "TF1.h"
//#include "G4GeometryTolerance.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Double_t poissonf(Double_t* x, Double_t* par) {
//  return par[0]*TMath::Poisson(x[0],par[1]);
//}



musrSteppingAction::musrSteppingAction()  {
  pointer=this;

  boolIsAnySpecialSaveVolumeDefined = false;
  boolIsVvvInfoRequested = false;
  boolMuonEventReweighting = false;
  boolCalculateFieldIntegral = false;
  myRootOutput = musrRootOutput::GetRootInstance();
  if (myRootOutput == NULL) {
    musrErrorMessage::GetInstance()->musrError(FATAL,
					       "musrSteppingAction::musrSteppingAction():  pointer to the musrRootOutput class not found! ==> EXECUTION STOPPED",true);
  }
  lastActualVolume="Unset";
}

musrSteppingAction::~musrSteppingAction()  {
}


 musrSteppingAction* musrSteppingAction::pointer=0;
 musrSteppingAction* musrSteppingAction::GetInstance()
{
  return pointer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void musrSteppingAction::DoAtTheBeginningOfEvent() {
  //  G4cout<<"musrSteppingAction::DoAtTheBeginningOfEvent:  at the beginning"<<G4endl;
  radioactiveElectronAlreadySavedInThisEvent=false;
  muAlreadyWasInTargetInThisEvent=false;
  muAlreadyWasInM0InThisEvent=false;
  muAlreadyWasInM1InThisEvent=false;
  muAlreadyWasInM2InThisEvent=false;
  myOldTracksMap.clear();
  indexOfOldTrack = -1;
  realTimeWhenThisEventStarted=time(0);
  BxIntegral=0;  ByIntegral=0;  BzIntegral=0; BzIntegral1=0;  BzIntegral2=0;  BzIntegral3=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrSteppingAction::UserSteppingAction(const G4Step* aStep)  {

  G4Track* aTrack = aStep->GetTrack();

  // kill the track, if required by user:
  if (aTrack->GetDefinition()) {
    G4String p_name = aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName();
    if      ((musrParameters::killAllPositrons)&&(p_name == "e+")) {aTrack->SetTrackStatus(fStopAndKill); return;}   // suspend the track
    else if ((musrParameters::killAllElectrons)&&(p_name == "e-")) {aTrack->SetTrackStatus(fStopAndKill); return;}
    else if ((musrParameters::killAllGammas)&&(p_name == "gamma")) {aTrack->SetTrackStatus(fStopAndKill); return;}
    else if ((musrParameters::killAllNeutrinos)&&((p_name == "nu_mu")||(p_name == "anti_nu_mu")
		  ||(p_name == "nu_e")||(p_name == "anti_nu_e")))  {aTrack->SetTrackStatus(fStopAndKill); return;}
  }
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  G4ThreeVector preStepPosition = preStepPoint->GetPosition();
  G4ThreeVector postStepPosition = postStepPoint->GetPosition();

  //  suspend the track if too many steps has already happened (relevant at high field)
  if (aTrack->GetCurrentStepNumber()>musrParameters::maximumNrOfStepsPerTrack) {
    char eMessage[200];
    sprintf(eMessage,"musrSteppingAction: Current number of steps for the track > %d ==> TRACK KILLED",
	    musrParameters::maximumNrOfStepsPerTrack);
    musrErrorMessage::GetInstance()->musrError(WARNING,eMessage,true);
    G4double x = postStepPosition.x()/CLHEP::mm;
    G4double y = postStepPosition.y()/CLHEP::mm;
    G4double z = postStepPosition.z()/CLHEP::mm;
    G4double E=aTrack->GetVertexKineticEnergy()/CLHEP::MeV;
    myRootOutput->htest1->Fill(x,y);
    myRootOutput->htest2->Fill(sqrt(x*x+y*y),z);
    myRootOutput->htest3->Fill(E);
    aTrack->SetTrackStatus(fStopAndKill);
  }

  //  abort the event if it takes too long to finish (e.g. more than 60 seconds)
  if ((time(0) - realTimeWhenThisEventStarted)>musrParameters::maximumTimePerEvent) {
    G4RunManager* fRunManager = G4RunManager::GetRunManager();
    G4cout<<"musrSteppingAction: event "<<fRunManager->GetCurrentEvent()->GetEventID()
  	  <<" aborted because calculation took already "<<musrParameters::maximumTimePerEvent<<" seconds."<<G4endl;
    char eMessage[200];
    sprintf(eMessage,"musrSteppingAction:  event aborted because its calculation takes more than %d seconds.",
	    musrParameters::maximumTimePerEvent);
    musrErrorMessage::GetInstance()->musrError(WARNING,eMessage,true);
    myRootOutput->SetEventWeight(0);
    fRunManager->AbortEvent();
  }

  // Temporary fix to avoid crashes caused by particles with unphysically high energies
  // (probably corrupted event?)
  //  if ((aStep->GetPreStepPoint()->GetKineticEnergy()) > (1*GeV)) {
  if ((preStepPoint->GetKineticEnergy()) > (100000*CLHEP::GeV)) {
    musrErrorMessage::GetInstance()->musrError(SERIOUS,
         "musrSteppingAction: kinetic energy of a particle larger than 1GeV!  STRANGE FOR muSR!",false);
    G4RunManager* fRunManager = G4RunManager::GetRunManager();
    G4cout<<"   Event nr.:"<<fRunManager->GetCurrentEvent()->GetEventID()
	  <<", the particle \""<<  aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName()
	  <<"\" has energy of "<<(preStepPoint->GetKineticEnergy())/CLHEP::GeV<<" GeV."<<G4endl;
    G4cout<<"             Deleting the event!"<<G4endl;
    G4cout.flush();
    myRootOutput->SetEventWeight(0);
    fRunManager->AbortEvent();
  }

  if (aTrack->GetDefinition()) {
    G4ParticleDefinition* p_definition = aTrack->GetDynamicParticle()->GetDefinition();
    G4String              p_name       = p_definition->GetParticleName();
    //    G4ProcessManager*     p_manager    = p_definition->GetProcessManager();
    G4LogicalVolume* actualLogicalVolume = aTrack->GetVolume()->GetLogicalVolume();
    G4String actualVolume = actualLogicalVolume->GetName();

    // Delete track if the particle is in the "kill" volume.
    // There is an example how to delete the track in example/novice/N04.
    // It is done in a different way here, because the example/novice/N04 was not doing
    // exactly what I wanted.
    if((actualVolume(0,8)=="log_kill")||(actualVolume(0,8)=="log_Kill")) {
      aTrack->SetTrackStatus(fStopAndKill);   // suspend the track
    }
    if ((p_name=="nu_mu")||(p_name=="anti_nu_mu")||(p_name=="nu_e")||(p_name=="anti_nu_e")) {
      //aTrack->SetTrackStatus(fStopAndKill);   // suspend the tracks of neutrinos
    }


    // Save info about the old tracks, if the user wishes to have Vvv info in the output Root Tree.
    if (boolIsVvvInfoRequested) {
      G4VPhysicalVolume* nextVolume = aTrack->GetNextVolume();
      if (nextVolume!=NULL) {
       	if ((nextVolume->GetLogicalVolume()->GetSensitiveDetector()!=NULL)||(actualLogicalVolume->GetSensitiveDetector()!=NULL)) {
	  G4int trackID = aTrack->GetTrackID();
	  std::map<G4int,G4int>::iterator itr;
	  itr = myOldTracksMap.find(trackID);
	  if (itr==myOldTracksMap.end()) {
	    // track with this trackID has not been found in the map (has not been saved yet) ==> save it
	    indexOfOldTrack++;
	    myOldTracksMap.insert( std::pair<G4int,G4int>(trackID,indexOfOldTrack) );
	    if (indexOfOldTrack<maxNumberOfOldTracks) {
	      particleID_oldTrack[indexOfOldTrack]     = aTrack->GetDefinition()->GetPDGEncoding();
	      parentTrackID_oldTrack[indexOfOldTrack]  = aTrack->GetParentID();
	      vertexKine_oldTrack[indexOfOldTrack]     = aTrack->GetVertexKineticEnergy();
	      vertexPosition_oldTrack[indexOfOldTrack] = aTrack->GetVertexPosition();
	      vertexLogVol_oldTrack[indexOfOldTrack]   = aTrack->GetLogicalVolumeAtVertex()->GetName();
	      if ((aTrack->GetCreatorProcess())!=NULL) { vertexProcess_oldTrack[indexOfOldTrack] = aTrack->GetCreatorProcess()->GetProcessName();}
	      else { vertexProcess_oldTrack[indexOfOldTrack] = "initialParticle";}
	    }
	    else {
	      musrErrorMessage::GetInstance()->musrError(WARNING,
		     "musrSteppingAction:  Maximum number of oldTracks reached ==> det_VvvXXX variables might be affected.",true);
	    }
	  }
	}
      }
    }


    //  This are the data just for the radioactive decay (when using the radioactive source):
    if ((musrParameters::boolG4GeneralParticleSource))  {
      //      &&(!radioactiveElectronAlreadySavedInThisEvent)) {
      if (aTrack->GetTrackID() != 1 ){
	if (aTrack->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay") {
          if (aTrack->GetDefinition()->GetParticleName()=="e-") {
	    if (aTrack->GetCurrentStepNumber()==1) {
	      G4double electron_kinetic_energy=preStepPoint->GetKineticEnergy();
	      myRootOutput->htest4->Fill(electron_kinetic_energy);
	    }
	  }
        }
      }
    }


    // Check if particle comes to the special volume
    if (boolIsAnySpecialSaveVolumeDefined) {
      //    G4bool isFirstStepInVolume=aStep->IsFirstStepInVolume();
      //          This does not work!!! (aStep->IsFirstStepInVolume() is always zero.)  I do not understand why!
      G4bool isFirstStepInVolume=false;
      if (actualVolume!=lastActualVolume) {
	lastActualVolume=actualVolume;
	isFirstStepInVolume=true;
      }

      if (isFirstStepInVolume) {
	G4int tmpVolumeID=saveVolumeMapping[actualVolume];
	if (tmpVolumeID!=0) {
	  G4int particle_id_save=p_definition->GetPDGEncoding();
	  G4double ke_save=preStepPoint->GetKineticEnergy();
	  G4double x_save=preStepPosition.x();
	  G4double y_save=preStepPosition.y();
	  G4double z_save=preStepPosition.z();
	  G4double time_save=preStepPoint->GetGlobalTime();
	  G4double px_save=preStepPoint->GetMomentum().x();
	  G4double py_save=preStepPoint->GetMomentum().y();
	  G4double pz_save=preStepPoint->GetMomentum().z();
	  G4double polx_save=preStepPoint->GetPolarization().x();
	  G4double poly_save=preStepPoint->GetPolarization().y();
	  G4double polz_save=preStepPoint->GetPolarization().z();
	  myRootOutput->SetSaveDetectorInfo(tmpVolumeID,particle_id_save,ke_save,x_save,y_save,z_save,time_save,px_save,py_save,pz_save,polx_save,poly_save,polz_save);
	  //
          //-----------------------------------------------------------------------------------------
          //  Uncoment for iterative musrSim runs (e.g. when searching for a quadrupole triplet focus using a python script)
	  // //  myRootOutput->htest7->Fill(sqrt(x_save*x_save+y_save*y_save));
	  // //  myRootOutput->htest7->Fill(x_save*x_save+y_save*y_save);
	  // myRootOutput->htest7->Fill(x_save);
	  // myRootOutput->htest8->Fill(y_save);
	  //-----------------------------------------------------------------------------------------
	}
      }
    }


    if ((p_name == "mu+")||(p_name == "mu-")||(p_name == "Mu"))  {
      // Store the information about the muon when it enters the target, M0, M1 or M2 for the fist time
      // in a given event (i.e. the code has to be called just once during the event).
      if ((actualVolume=="log_target")||(actualVolume=="log_Target")) {
	if (!muAlreadyWasInTargetInThisEvent) {
	  muAlreadyWasInTargetInThisEvent=true;
	  myRootOutput->SetPolInTarget(aTrack->GetPolarization());
	  myRootOutput->SetTimeInTarget(aTrack->GetGlobalTime());
	  myRootOutput->SetMomentumInTarget(preStepPoint->GetMomentum());
	}
      }
      else if ((actualVolume=="log_M0")||(actualVolume=="log_m0")) {
	if (!muAlreadyWasInM0InThisEvent) {
	  muAlreadyWasInM0InThisEvent=true;
	  myRootOutput->SetPolInM0(aTrack->GetPolarization());
	  myRootOutput->SetTimeInM0(aTrack->GetGlobalTime());
	}
      }
      else if ((actualVolume=="log_M1")||(actualVolume=="log_m1")) {
	if (!muAlreadyWasInM1InThisEvent) {
	  muAlreadyWasInM1InThisEvent=true;
	  myRootOutput->SetPolInM1(aTrack->GetPolarization());
	  myRootOutput->SetTimeInM1(aTrack->GetGlobalTime());
	}
      }
      else if ((actualVolume=="log_M2")||(actualVolume=="log_m2")) {
	if (!muAlreadyWasInM2InThisEvent) {
	  muAlreadyWasInM2InThisEvent=true;
	  myRootOutput->SetPolInM2(aTrack->GetPolarization());
	  myRootOutput->SetTimeInM2(aTrack->GetGlobalTime());
	}
      }


      // Calculate the field integral along the muon path, if requested by the user.
      // (2008.10.20. - idea of Robert is to calculate the distribution of the field integral for different muon paths
      // to see what  (delta_I/I)  is still acceptable for hith field muSR.  To calculate it properly, the user must
      // force Geant to use very small step size.
      if (boolCalculateFieldIntegral) {
	if (F04GlobalField::Exists()) {
	  CoordinateForFieldIntegral[0] = postStepPosition.x();
	  CoordinateForFieldIntegral[1] = postStepPosition.y();
	  CoordinateForFieldIntegral[2] = postStepPosition.z();
	  CoordinateForFieldIntegral[3] = aTrack->GetGlobalTime();
	  F04GlobalField::getObject() -> GetFieldValue(CoordinateForFieldIntegral,FieldForFieldIntegral);
	  //	  G4cout<<"B=("<<FieldForFieldIntegral[0]/tesla<<","<<FieldForFieldIntegral[1]/tesla<<","<<FieldForFieldIntegral[2]/tesla<<")"<<G4endl;
	  G4double stepLength = aStep->GetStepLength();
	  //	  BxIntegral += stepLength;
	  //	  ByIntegral += stepLength;
	  //	  BzIntegral += stepLength;
	  BxIntegral += stepLength * FieldForFieldIntegral[0];
	  ByIntegral += stepLength * FieldForFieldIntegral[1];
	  BzIntegral += stepLength * FieldForFieldIntegral[2];
	  BzIntegral1 += ( postStepPosition.z() - preStepPoint->GetPosition().z() ) * FieldForFieldIntegral[2];
	  BzIntegral2 += stepLength;
	  BzIntegral3 += ( postStepPosition.z() - preStepPoint->GetPosition().z() );
	  //	  G4cout<<"BzIntegral="<<BzIntegral<<"  stepLength="<<stepLength<<"FieldForFieldIntegral[2]="<<FieldForFieldIntegral[2]<<G4endl;
	}
      }

      //     Pick up process "DecayWithSpin":
      const G4VProcess* process = postStepPoint->GetProcessDefinedStep();
      if (process!=NULL) {
	G4String processName = process->GetProcessName();
	if (processName=="DecayWithSpin") {
	  //	  std::cout<<"musrSteppingAction: DecayWithSpin"<<std::endl;
	  musrParameters::field_DecayWithSpin=true;

	  // Test whether the event reweighting is requeseted for this volume
	  // (i.e. user may request reweighting of events depending on the volume,
          // in which the muon stops and decays).
	  if (boolMuonEventReweighting) {
	    G4int weight = volumeMuonWeightMapping[actualVolume];
	    if (weight!=0) {
	      G4double randomNumber = weight * G4UniformRand();
	      if (randomNumber < (weight-1.)) {
		//	      G4cout<<"Event will be aborted"<<G4endl;
		musrErrorMessage::GetInstance()->musrError(INFO,
			        "musrSteppingAction:  event deleted because of the reweighting.",true);
		G4RunManager* fRunManager = G4RunManager::GetRunManager();
		weight=0;
		fRunManager->AbortEvent();
	      }
	      //	    else {
	      //	      G4cout<<"Event will be reweighted"<<G4endl;
	      //	    }
	      myRootOutput->SetEventWeight(weight);
	    }
	  }

	  // Store the information about the decaying muon and store it in the Root tree
	  G4double timeOfDecay_tmp=aTrack->GetGlobalTime();
	  G4double BFieldAtOrigin[6] = {0.,0.,0.,0.,0.,0.};
	  G4double PointOfDecay[4] ={postStepPosition.x(),postStepPosition.y(),postStepPosition.z(),timeOfDecay_tmp};
	  if (F04GlobalField::Exists()) {
	    F04GlobalField* myGlobalField = F04GlobalField::getObject();
	    myGlobalField->GetFieldValue(PointOfDecay,BFieldAtOrigin);
	  }
	  else {
	    G4FieldManager *fMgr=G4TransportationManager::GetTransportationManager()->GetFieldManager();
	    //	    G4cout<<"Debug 1"<<G4endl;
	    if (fMgr!=NULL) {
	      //	      G4cout<<"Debug 2"<<G4endl;
	      if(!fMgr->DoesFieldChangeEnergy())  {            //then we have a magnetic field
		//		G4cout<<"Debug 3"<<G4endl;
		fMgr->GetDetectorField()->GetFieldValue(PointOfDecay,BFieldAtOrigin);
	      }
	    }
	  }
	  myRootOutput->SetDecayTime(timeOfDecay_tmp);
	  myRootOutput->SetDecayPolarisation(aTrack->GetPolarization());
	  myRootOutput->SetDecayPosition(postStepPosition);
	  myRootOutput->SetDecayDetectorID(actualVolume);
	  myRootOutput->SetBField(BFieldAtOrigin);
	  if (boolCalculateFieldIntegral)  {
	    myRootOutput->SetBFieldIntegral(BxIntegral,ByIntegral,BzIntegral,BzIntegral1,BzIntegral2,BzIntegral3);
	  }
	  //	  G4cout<<"==================================================  BzIntegral = "<<BzIntegral<<G4endl;

	  // store the information about the emerging positron
	  const G4TrackVector* secondary = fpSteppingManager->GetSecondary();
	  G4int n_secondaries= (*secondary).size();
	  for (G4int i=0; i<n_secondaries; i++) {
	    if ( ((*secondary)[i]->GetDefinition()->GetParticleName()) == "e+" ) {
	      myRootOutput->SetInitialPositronMomentum((*secondary)[i]->GetMomentum());
	    }
	  }
	}
      }
    }

    else {   // particle is not muon
      // Delete track if the particle is far away from the detector (i.e. in the "shield" volume).
      // There is an example how to delete the track in example/novice/N04.
      // It is done in a different way here, because the example/novice/N04 was not doing
      // exactly what I wanted.
      if((actualVolume(0,10)=="log_shield")||(actualVolume(0,10)=="log_Shield")) {
	aTrack->SetTrackStatus(fStopAndKill);   // suspend the track
      }
    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void musrSteppingAction::SetLogicalVolumeAsSpecialSaveVolume(G4String logicName, G4int volumeID)  {
  boolIsAnySpecialSaveVolumeDefined = true;
  saveVolumeMapping[logicName]=volumeID;
}


void musrSteppingAction::SetVolumeForMuonEventReweighting(G4String logicName, G4int weight) {
  boolMuonEventReweighting = true;
  volumeMuonWeightMapping[logicName]=weight;
}

G4bool  musrSteppingAction::GetInfoAboutOldTrack(G4int trackID, G4int& parentTrackID, G4int& particleID, G4double& vertexKine,
						 G4ThreeVector& vertexPosition, G4String& vertexLogVol, G4String& vertexProcess) {
  //  G4int ind = myOldTracksMap[trackID];
  //  G4cout<<"musrSteppingAction::GetInfoAboutOldTrack: trackID="<<trackID<<"\t myOldTracksMap[trackID]="<<myOldTracksMap[trackID]<<G4endl;
  std::map<G4int,G4int>::iterator itr;
  itr = myOldTracksMap.find(trackID);
  if ( itr==myOldTracksMap.end() ) {
    //  if ((ind==0)||(ind>=maxNumberOfOldTracks)) {
    char eMessage[200];
    sprintf(eMessage,"musrSteppingAction::GetInfoAboutOldTrack: trackID not found in myOldTracksMap, det_VvvXXX variables might be affected");
    musrErrorMessage::GetInstance()->musrError(WARNING,eMessage,false);
    G4cout<<"                                  Requested trackID="<<trackID<<G4endl;
    //    G4cout<<"Saved tracks:"<<G4endl;
    //    for (itr=myOldTracksMap.begin(); itr!=myOldTracksMap.end(); itr++) {
    //      G4cout<<"first="<<itr->first<<"\tsecond="<<itr->second<<G4endl;
    //    }
    return false;
  }
  else {
    G4int ind = itr->second;
    if (ind>=maxNumberOfOldTracks) {
      G4cout<<"musrSteppingAction::GetInfoAboutOldTrack:  Problem!  ind>maxNumberOfOldTracks! ("<<ind<<">"<<maxNumberOfOldTracks<<")"<<G4endl;
      G4cout<<"                                                     itr->first = "<<itr->first<<",  trackID = "<<trackID<<G4endl;
      return false;
    }
    parentTrackID=parentTrackID_oldTrack[ind];
    if (trackID==parentTrackID) {
      G4cout<<"musrSteppingAction::GetInfoAboutOldTrack:  Problem!  trackID==parentTrackID! ("<<trackID<<"=="<<parentTrackID<<")"<<G4endl;
      return false;
    }
    particleID=particleID_oldTrack[ind];
    vertexKine=vertexKine_oldTrack[ind];
    vertexPosition= vertexPosition_oldTrack[ind];
    vertexLogVol=vertexLogVol_oldTrack[ind];
    vertexProcess=vertexProcess_oldTrack[ind];
  }
  //  G4cout<<"GetInfoAboutOldTrack: trackID="<<trackID<<"\t parentTrackID="<<parentTrackID<<"\t particleID="<<particleID;
  //  G4cout<<"\t vertexKine="<<vertexKine<<"\t vertexLogVol="<<vertexLogVol<<"\t vertexProcess="<<vertexProcess<<G4endl;
  return true;
}




G4bool  musrSteppingAction::AreTracksCommingFromSameParent(G4int trackID1, G4int trackID2, G4String volumeName){
  // There are two tracks with different track IDs.  This routine finds the parents of both of them,
  // which were created outside logical volume "volumeID".  If both tracks have the same parent, the
  // functions returns "true".
  std::map<G4int,G4int>::iterator itr;
  G4int ind;

  G4int track1;
  G4int trID = trackID1;
  do {
    track1=trID;
    itr = myOldTracksMap.find(trID);
    if ( itr==myOldTracksMap.end() ) {
      G4cout<<"musrSteppingAction::AreTracksCommingFromSameParent()  Strange, trackID1 ="<<trackID1<<" not found"<<G4endl;
      return false;
    }
    ind   = itr->second;
    trID = parentTrackID_oldTrack[ind];
  } while (vertexLogVol_oldTrack[ind]==volumeName);

  G4int track2;
  trID = trackID2;
  do {
    track2=trID;
    itr = myOldTracksMap.find(trID);
    if ( itr==myOldTracksMap.end() ) {
      G4cout<<"musrSteppingAction::AreTracksCommingFromSameParent()  Strange, trackID2 ="<<trackID2<<" not found"<<G4endl;
      return false;
    }
    ind   = itr->second;
    trID = parentTrackID_oldTrack[ind];
  } while (vertexLogVol_oldTrack[ind]==volumeName);


  if (track1==track2) {return true;}
    //G4cout<<"track1="<<track1<<"\ttrack2="<<track2<<G4endl; return true;}
  //  G4cout<<"\t\t\t\ttrack1="<<track1<<"\ttrack2="<<track2<<G4endl;
  return false;
}

//  //Double_t musrSteppingAction::poissonf(Double_t* x, Double_t* par)                                         
//  //{ 
//Double_t musrSteppingAction::poissonf(Double_t* x, Double_t* par) {
//  //  if (x<0)
//  //    return 0;
//  //  else if (x == 0.0)
//  //    return 1./Math::Exp(par);
//  //  else {
//  //    Double_t lnpoisson = x*log(par)-par-LnGamma(x+1.);
//  //    return Exp(lnpoisson);
//  //  }
//  return par[0]*TMath::Poisson(x[0],par[1]);
//}

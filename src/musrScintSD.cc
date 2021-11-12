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

#include "musrScintSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include <algorithm>   // needed for the sort() function
#include "G4VProcess.hh"  // needed for the degugging message of the process name
#include "G4OpBoundaryProcess.hh"      // Optical photon process
#include "G4RunManager.hh"
#include "musrParameters.hh"
#include "musrErrorMessage.hh"
#include "musrSteppingAction.hh"
//#include "TCanvas.h"
#include "TH1D.h"
#include "TMath.h"
#include "TF1.h"
#include <vector>

//bool myREMOVEfunction (int i,int j) { return (i<j); }
//bool timeOrdering (musrScintHit hit1, musrScintHit hit2) { 
//  return (hit1.GetGlobalTime()<hit2.GetGlobalTime());
//}
//
//bool timeOrdering2 (std::map<int,double>::iterator i1, std::map<int,double>::iterator m2) {
//  return ( (*i1).first()<(*i2).second() );
//}
//
//bool timeOrdering2 (std::pair<int,double> p1, std::pair<int,double> p2) {
//  return ( p1.first()<p2.second() );
//}
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Double_t poissonf(Double_t* x, Double_t* par) {
//  return par[0]*TMath::Poisson(x[0],par[1]);
//}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


const G4double musrScintSD::OPSA_C1_threshold (0.5);
const G4double musrScintSD::OPSA_C2_threshold (1.);
const G4double musrScintSD::OPSA_C3_threshold (2.);
const G4double musrScintSD::OPSA_C4_threshold (3.);
const G4double musrScintSD::OPSA_C5_threshold (4.);
const G4double musrScintSD::OPSA_C6_threshold (5.);
const G4double musrScintSD::OPSA_C7_threshold (6.);
const G4double musrScintSD::OPSA_C8_threshold (7.);
const G4double musrScintSD::OPSA_C9_threshold (8.);
const G4double musrScintSD::OPSA_C10_threshold (9.);
const G4double musrScintSD::OPSA_C11_threshold (10.);
const G4double musrScintSD::OPSA_C12_threshold (11.);
const G4double musrScintSD::OPSA_C13_threshold (12.);
const G4double musrScintSD::OPSA_C14_threshold (13.);
const G4double musrScintSD::OPSA_C15_threshold (14.);
const G4double musrScintSD::OPSA_C16_threshold (15.);
const G4double musrScintSD::OPSA_C17_threshold (16.);
const G4double musrScintSD::OPSA_C18_threshold (17.);
const G4double musrScintSD::OPSA_C19_threshold (18.);
const G4double musrScintSD::OPSA_C20_threshold (19.);
const G4double musrScintSD::OPSA_C21_threshold (20.);
const G4double musrScintSD::OPSA_C22_threshold (22.);
const G4double musrScintSD::OPSA_C23_threshold (24.);
const G4double musrScintSD::OPSA_C24_threshold (26.);
const G4double musrScintSD::OPSA_C25_threshold (28.);
const G4double musrScintSD::OPSA_C26_threshold (30.);
const G4double musrScintSD::OPSA_C27_threshold (35.);
const G4double musrScintSD::OPSA_C28_threshold (40.);
const G4double musrScintSD::OPSA_C29_threshold (50.);
const G4double musrScintSD::OPSA_C30_threshold (60.);
const G4double musrScintSD::OPSA_C31_threshold (70.);
const G4double musrScintSD::OPSA_C32_threshold (80.);
const G4double musrScintSD::OPSA_C33_threshold (90.);
const G4double musrScintSD::OPSA_C34_threshold (100.);
const G4double musrScintSD::OPSA_C35_threshold (125.);
const G4double musrScintSD::OPSA_C36_threshold (150.);
const G4double musrScintSD::OPSA_C37_threshold (170.);
const G4double musrScintSD::OPSA_C38_threshold (200.);
const G4double musrScintSD::OPSA_C39_threshold (250.);
const G4double musrScintSD::OPSA_C40_threshold (300.);

musrScintSD::musrScintSD(G4String name)
:G4VSensitiveDetector(name)
{
  pointer=this;
  G4String HCname;
  collectionName.insert(HCname="scintCollection");
  OPSA_minNrOfDetectedPhotons = 1;
  OPSA_signalSeparationTime   = 10.;
  OPSA_fracA = 0.01;
  OPSA_fracB = 0.05;
  OPSA_C_threshold = 0.5;
  OPSA_D_threshold = 0.5;
  bool_multimapOfEventIDsForOPSAhistosEXISTS = false;
  bool_StoreThisOPSAhistSUMMED = false;
  bool_StoreThisOPSAhistALL = false;
  OPSAhistoNbin = 100;
  OPSAhistoMin =0;
  OPSAhistoMax = 10.;
  OPSAhistoBinWidth=0.1;
  OPSAhistoBinWidth1000=100.;
  bool_pulseShapeExists=false;
  OPSA_CFD_a1 = -0.2;
  OPSA_CFD_delay = 2.;
  OPSA_CFD_timeShiftOffset = 0.;
  APDcell_nx =1; APDcell_ny=10; APDcell_nz=10; 
  APDcell_ax =0.3; APDcell_ay=0.3; APDcell_az=0.3;
  APDcellsEffectRequested=false;
  APDcellsTimeVariationRequested=false;
  APDcrossTalkRequested=false;
  APDcrossTalkProb=0.;
  // verboseLevel = 9; // JSL
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrScintSD::~musrScintSD(){ }

musrScintSD* musrScintSD::pointer=NULL;
musrScintSD* musrScintSD::GetInstance() {return pointer;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrScintSD::Initialize(G4HCofThisEvent* HCE) {
  // This "Initialize" method is called at the beginning of each event
  //   --> perhaps some this code could be executed once per run?
  if (verboseLevel>1) G4cout<<"VERBOSE 2:  musrScintSD::Initialize\n";
  scintCollection = new musrScintHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID<0) { 
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    if (verboseLevel>1) G4cout<<"VERBOSE 2:  musrScintSD::HCID was <0\n, now HCID="<<HCID<<"\n";
  }
  HCE->AddHitsCollection( HCID, scintCollection ); 
  myStoreOnlyEventsWithHits = musrParameters::storeOnlyEventsWithHits;
  mySignalSeparationTime    = musrParameters::signalSeparationTime;
  myStoreOnlyTheFirstTimeHit= musrParameters::storeOnlyTheFirstTimeHit;
  myStoreOnlyEventsWithHitInDetID = musrParameters::storeOnlyEventsWithHitInDetID;
  myStoreOnlyEventsWithMuonsDecayedInDetID = musrParameters::storeOnlyEventsWithMuonsDecayedInDetID;
  musrSteppingAction* myMusrSteppingAction = musrSteppingAction::GetInstance();
  boolIsVvvInfoRequested = myMusrSteppingAction->IsVvvInfoRequested();
  myRootOutput = musrRootOutput::GetRootInstance();

 // In case of optical photons, delete all optHitDetectorMap* from the previous event (if they exist).
  if (musrParameters::boolG4OpticalPhotons) {
    if (!musrParameters::boolG4OpticalPhotonsUnprocess) {
      for (optHitMapType::const_iterator it=optHitMap.begin() ; it != optHitMap.end(); it++ ) {
		  if (verboseLevel>1) G4cout<<"VERBOSE 2 : deleting optHitDetectorMap" <<"\n";
	delete (it->second);
      }
		  if (verboseLevel>1) G4cout<<"VERBOSE 2 : clearing optHitMap" <<"\n";
      optHitMap.clear();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool musrScintSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  if (verboseLevel>1) G4cout<<"VERBOSE 2:  musrScintSD::ProcessHits\n";
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) {
    return false;
  }

  G4Track* aTrack = aStep->GetTrack();
  //  G4LogicalVolume* hitLogicalVolume = aTrack->GetVolume()->GetLogicalVolume();
  //  G4String actualVolume=aTrack->GetVolume()->GetLogicalVolume()->GetName();
  G4String hitLogicalVolumeName = aTrack->GetVolume()->GetLogicalVolume()->GetName();
  G4String particleName=aTrack->GetDefinition()->GetParticleName();
  //  if (particleName=="opticalphoton") {G4cout<<"UFON JE TU:  edep="<<edep<<G4endl; return false;}
  if (particleName=="opticalphoton") {
    //    G4double APDcellsTimeVariation = G4RandGauss::shoot(0,APDcellsTimeVariationSigma);
    //    if (!musrParameters::boolG4OpticalPhotonsUnprocess) ProcessOpticalPhoton(aStep,APDcellsTimeVariation);
    if (!musrParameters::boolG4OpticalPhotonsUnprocess) ProcessOpticalPhoton(aStep);
    return false;
  }

  // If requested, store only the hit that happened first (usefull for some special studies, not for a serious simulation)
  if (myStoreOnlyTheFirstTimeHit) {
    G4int NbHits = scintCollection->entries(); 
    if (NbHits>0) {
      aTrack->SetTrackStatus(fStopAndKill);
      return false;
    }
  }
  
  

  if (verboseLevel>1) G4cout<<"VERBOSE 2 : creating a new musrScintHit" <<"\n";
  musrScintHit* newHit = new musrScintHit();
  //  newHit->SetParticleName (aTrack->GetDefinition()->GetParticleName());
  newHit->SetParticleName(particleName);
  G4int particleID = aTrack->GetDefinition()->GetPDGEncoding();
  newHit->SetParticleID (particleID);
  newHit->SetEdep     (edep);
  newHit->SetPrePos   (aStep->GetPreStepPoint()->GetPosition());
  newHit->SetPostPos  (aStep->GetPostStepPoint()->GetPosition());
  newHit->SetPol      (aTrack->GetPolarization());
  //  G4LogicalVolume* hitLogicalVolume = aTrack->GetVolume()->GetLogicalVolume();
  newHit->SetLogVolName(hitLogicalVolumeName);
  //  newHit->SetLogVolName(hitLogicalVolume->GetName());
  newHit->SetGlobTime(aTrack->GetGlobalTime());
  //  Warning - aStep->IsFirstStepInVolume() only available in Geant version >= 4.8.2 !
  //  newHit->SetFirstStepInVolumeFlag(aStep->IsFirstStepInVolume());
  //  newHit->SetLastStepInVolumeFlag(aStep->IsLastStepInVolume());
  newHit->SetKineticEnergy(aTrack->GetKineticEnergy());
  // newHit->SetKineticEnergy(aTrack->GetKineticEnergy()+edep);
  G4double vertexKineticEnergy = aTrack->GetVertexKineticEnergy();
  newHit->SetVertexKineticEnergy(vertexKineticEnergy);
  G4ThreeVector vertexPosition = aTrack->GetVertexPosition();
  newHit->SetVertexPosition(vertexPosition);
  const G4LogicalVolume* vertexLogicalVolume = aTrack->GetLogicalVolumeAtVertex();
  G4String vertexLogicalVolumeName = vertexLogicalVolume->GetName();
  newHit->SetLogicalVolumeAtVertex(vertexLogicalVolumeName);
  G4String processName;
  if ((aTrack->GetCreatorProcess())!=0) { processName=aTrack->GetCreatorProcess()->GetProcessName(); }
  else {processName="initialParticle";}   //if no process found, the track comes from the generated particle
  newHit->SetCreatorProcessName(processName);
  G4int trackID = aTrack->GetTrackID();
  newHit->SetTrackID  (trackID);
  newHit->SetStepLength   (aStep->GetStepLength());

  scintCollection->insert( newHit );
  //  newHit->Print();
  newHit->Draw();
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void musrScintSD::ProcessOpticalPhoton(G4Step* aStep) {
//void musrScintSD::ProcessOpticalPhoton(G4Step* aStep, G4double APDcellsTimeVariation) {
  //  //Was the photon absorbed by the absorption process ?
  //  const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
  //  G4String processName =  (process) ? process->GetProcessName() : "Unknown";
  //  G4cout<<"processName="<<processName<<G4endl;
  //  if (processName=="OpAbsorption") {
  //    G4cout<<"\n ABSORPTION\n";
  //  }

  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  static G4OpBoundaryProcess* boundaryProc=NULL;
  //find the boundary process only once
  if(!boundaryProc){
    G4ProcessManager* pm = aStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    for(G4int i=0; i<nprocesses; i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
	boundaryProc = (G4OpBoundaryProcess*)(*pv)[i];
	break;
      }
    }
  }
  
  boundaryStatus=boundaryProc->GetStatus();
  
  //Check to see if the partcile was actually at a boundary
  //Otherwise the boundary status may not be valid
  //Prior to Geant4.6.0-p1 this would not have been enough to check
  if(aStep->GetPostStepPoint()->GetStepStatus()==fGeomBoundary){
    //    G4cout<<"      boundaryStatus="<<boundaryStatus<<"   ";
    //    G4String actualVolume = aStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetName();
    G4String actualVolume = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName();
    Int_t detID =  myRootOutput->ConvertVolumeToID(actualVolume);
    //    G4cout<<" detID ="<<detID<<"  actualVolume="<<actualVolume<<G4endl;
    optHitMapType::iterator iter = optHitMap.find(detID);
    if (iter==optHitMap.end()) {  // optHitDetectorMapType does not exist ==> create it
      if (verboseLevel>1) G4cout<<"VERBOSE 2 : creating a new OptHitDetectorMap" <<"\n";
      optHitDetectorMapType* optHitDetectorMapTmp = new optHitDetectorMapType;
      optHitMap.insert(std::pair<Int_t,optHitDetectorMapType*>(detID,optHitDetectorMapTmp));
      iter = optHitMap.find(detID);
    }
    optHitDetectorMapType* optHitDetectorMap = (*iter).second;
    //	  optHitDetectorMapType optHitDetectorMap = optHitMap[detID];
    G4double tmpTime = aStep->GetPostStepPoint()->GetGlobalTime();
    //    G4cout<<"   tmpTime="<<tmpTime;
    //	  optHitDetectorMap->insert(std::pair<G4double,G4int>(tmpTime,boundaryStatus));
    
    if (boundaryStatus!=Detection) {
      char message[200];
      sprintf(message,"musrScintSD.cc::ProcessOpticalPhoton(): Optical photon boundary status is not Detection but %i",boundaryStatus);
      //      musrErrorMessage::GetInstance()->musrError(FATAL,message,false);
      musrErrorMessage::GetInstance()->musrError(WARNING,message,true);
    }

    // Do APD cell variation time if requested.  Even if not requested, generate the random numbers in order to
    // have reproducible simulation.
    G4double APDcellsTimeVariation = G4RandGauss::shoot(0,APDcellsTimeVariationSigma);
    if (APDcellsTimeVariationRequested) tmpTime += APDcellsTimeVariation;
    
    // If the simulation of the effect of the finite number of cells in APD is requested, find out
    // whether a given cell already fired - if so, do not store the photon.
    G4int APDcellID = (APDcellsEffectRequested) ? FindAPDcellID(aStep) : 0;
    FireAPDcell(optHitDetectorMap,APDcellID,tmpTime,1);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrScintSD::EndOfEvent(G4HCofThisEvent*) {
  if (verboseLevel>1) { 
    G4cout<<"VERBOSE 2:  musrScintSD::EndOfEvent"<<G4endl;
    G4int NbHits = scintCollection->entries();
    G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
	   << " hits in the scint chambers: " << G4endl;
  } 
  

  G4int NbHits = scintCollection->entries(); 
  
  if (myStoreOnlyEventsWithHits) {
    if (NbHits<=0) {
      return;
    }
    else if (myStoreOnlyEventsWithHitInDetID!=0) {
      for (G4int i=0; i<NbHits; i++) {
	musrScintHit* aHit = (*scintCollection)[i];
	G4String aHitVolumeName = aHit->GetLogVolName();
	G4int    aHitVolumeID   = myRootOutput->ConvertVolumeToID(aHitVolumeName);
	if (aHitVolumeID==myStoreOnlyEventsWithHitInDetID) break;  // hit in the requested detector was identified
	if (i==(NbHits-1)) return;                              // no hit identified in the requested detector
      }
    }
  }

  if (myStoreOnlyEventsWithMuonsDecayedInDetID) {
    if ((myRootOutput->GetDecayDetectorID())!=myStoreOnlyEventsWithMuonsDecayedInDetID) return;
  }

  //  Sort out hits and fill them into root
  if (NbHits>0) {
    const G4int det_IDmax = musrRootOutput::det_nMax;
    G4double det_edep[det_IDmax];
    G4int    det_nsteps[det_IDmax];
    G4double det_length[det_IDmax];
    G4int    det_ID[det_IDmax]; 
    G4double det_edep_el[det_IDmax];
    G4double det_edep_pos[det_IDmax];
    G4double det_edep_gam[det_IDmax];
    G4double det_edep_mup[det_IDmax];
    G4double det_time_start[det_IDmax];
    G4double det_time_end[det_IDmax];
    G4double det_x[det_IDmax];
    G4double det_y[det_IDmax];
    G4double det_z[det_IDmax];
    G4double det_kine[det_IDmax];
    G4double det_VrtxKine[det_IDmax];
    G4double det_VrtxX[det_IDmax];
    G4double det_VrtxY[det_IDmax];
    G4double det_VrtxZ[det_IDmax];
    G4int    det_VrtxVolID[det_IDmax];
    G4int    det_VrtxProcID[det_IDmax];
    G4int    det_VrtxTrackID[det_IDmax];
    G4int    det_VrtxParticleID[det_IDmax];
    G4int det_VvvTrackSign[det_IDmax];

    //  Sort hits according to the time.  Using std::map is convenient, because it sorts
    //  its entries according to the key (the first variable of the pair).
    std::multimap<G4double,G4int> myHitTimeMapping;    // "map" replaced by "multimap" (needed for radioactive decay,
                                                       // in which case times are huge and due to the limited rounding 
                                                       // precision become equal  --> map ignores the same "keys",
                                                       // multimap does not.
    std::multimap<G4double,G4int>::iterator it;
    for (G4int i=0; i<NbHits; i++) {
      musrScintHit* aHit = (*scintCollection)[i];
      G4double tmptime=aHit->GetGlobalTime();
      //      G4cout<<"Hit nr "<<i<<"  at time="<<tmptime<<"  with edep="<<aHit->GetEdep()/MeV
      //	    <<"   detID="<<myRootOutput->ConvertVolumeToID(aHit->GetLogVolName())<< G4endl;
      myHitTimeMapping.insert ( std::pair<G4double,G4int>(tmptime,i) );
    }

    //  Loop over all hits (which are sorted according to their time):
    G4int nSignals=0;
    for (it=myHitTimeMapping.begin(); it!=myHitTimeMapping.end(); it++) {
      //      G4cout << "Key:" << it->first;
      //      G4cout << "  Value:" << it->second << "\n";
      G4int ii = it->second;      // ii  is the index of the hits, which is sorted according to time
      musrScintHit* aHit = (*scintCollection)[ii];
      G4String aHitVolumeName = aHit->GetLogVolName();
      G4int    aHitVolumeID   = myRootOutput->ConvertVolumeToID(aHitVolumeName);
      G4double aHitTime       = aHit->GetGlobalTime();
      G4int    aHitTrackID    = aHit->GetTrackID();

      // Loop over all already defined signals and check whether the hit falls into any of them
      G4bool signalAssigned=false;
      for (G4int j=0; j<nSignals; j++) {
	if ( (aHitVolumeID==det_ID[j]) && ((aHitTime-det_time_end[j])<mySignalSeparationTime) ) {
	  signalAssigned=true;
	  det_edep[j]                 += aHit->GetEdep();
	  det_nsteps[j]++;
	  det_length[j]               += aHit->GetStepLength();
	  det_time_end[j]              = aHitTime;
	  G4String aParticleName       = aHit->GetParticleName();
	  if (aParticleName=="e-") {
	    det_edep_el[j]            += aHit->GetEdep();
	  }  else if (aParticleName=="e+") {
	    det_edep_pos[j]           += aHit->GetEdep();
	  }  else if (aParticleName=="gamma") {
	    det_edep_gam[j]           += aHit->GetEdep();
	  }  else if ((aParticleName=="mu+")||(aParticleName=="mu-")) {
	    det_edep_mup[j]           += aHit->GetEdep();
	  }  else {
	    char message[200];
	    sprintf(message,"musrScintSD.cc::EndOfEvent(): untreated particle \"%s\" deposited energy.",aParticleName.c_str());
	    musrErrorMessage::GetInstance()->musrError(WARNING,message,true);
	  }
	  // Check whether the signals consits of more then just one hit,  in which case make the track ID negative:
	  if (abs(det_VrtxTrackID[j])!=aHitTrackID) {
	    det_VrtxTrackID[j]=-1*abs(det_VrtxTrackID[j]);
	    if (boolIsVvvInfoRequested) {
	      if (det_VvvTrackSign[j]==1) {
		musrSteppingAction* myMusrSteppingAction = musrSteppingAction::GetInstance();
		if (!(myMusrSteppingAction->AreTracksCommingFromSameParent(aHitTrackID,abs(det_VrtxTrackID[j]),aHitVolumeName))) {det_VvvTrackSign[j]=-1;}
	      }
	    }
	  }
	  break;
	}
      }
      if (!signalAssigned) {    // The hit does not belong to any existing signal --> create a new signal.
	// Check, whether the maximum number of signals was not exceeded:
	if ( nSignals >= (det_IDmax-1) ) {
	  char message[200];
	  sprintf(message,"musrScintSD.cc::EndOfEvent(): number of signals exceeds maximal allowed value.");
	  musrErrorMessage::GetInstance()->musrError(WARNING,message,true);
	}
	else {
	  det_edep[nSignals]                  = aHit->GetEdep();
	  det_nsteps[nSignals]                = 1;
	  det_length[nSignals]                = aHit->GetStepLength();
	  det_ID[nSignals]                    = aHitVolumeID;
	  det_time_start[nSignals]            = aHitTime;
	  det_time_end[nSignals]              = aHitTime;
	  det_edep_el[nSignals]               = 0;
	  det_edep_pos[nSignals]              = 0;
	  det_edep_gam[nSignals]              = 0;
	  det_edep_mup[nSignals]              = 0;
	  G4String aParticleName              = aHit->GetParticleName();
	  if (aParticleName=="e-") {
	    det_edep_el[nSignals]            += aHit->GetEdep();
	  }  else if (aParticleName=="e+") {
	    det_edep_pos[nSignals]           += aHit->GetEdep();
	  }  else if (aParticleName=="gamma") {
	    det_edep_gam[nSignals]           += aHit->GetEdep();
	  }  else if ((aParticleName=="mu+")||(aParticleName=="mu-")) {
	    det_edep_mup[nSignals]           += aHit->GetEdep();
	  }  else {
	    char message[200];
	    sprintf(message,"musrScintSD.cc::EndOfEvent(): UNTREATED PARTICLE \"%s\" deposited energy.",aParticleName.c_str());
	    musrErrorMessage::GetInstance()->musrError(WARNING,message,true);
	  }
	  G4ThreeVector prePos = aHit->GetPrePos();
	  det_x[nSignals]=prePos.x();
	  det_y[nSignals]=prePos.y();
	  det_z[nSignals]=prePos.z();
	  det_kine[nSignals]             = aHit->GetKineticEnergy();
	  det_VrtxKine[nSignals]         = aHit->GetVertexKineticEnergy();
	  G4ThreeVector VrtxPos          = aHit->GetVertexPosition();
	  det_VrtxX[nSignals]            = VrtxPos.x();
	  det_VrtxY[nSignals]            = VrtxPos.y();
	  det_VrtxZ[nSignals]            = VrtxPos.z();
	  G4String logicalVolumeAtVertex = aHit->GetLogicalVolumeAtVertex();
	  det_VrtxVolID[nSignals]        = myRootOutput->ConvertVolumeToID(logicalVolumeAtVertex);
	  G4String creatorProcessName    = aHit->GetCreatorProcessName();
	  det_VrtxProcID[nSignals]       = myRootOutput->ConvertProcessToID(creatorProcessName);
	  det_VrtxTrackID[nSignals]      = aHit->GetTrackID();
	  det_VrtxParticleID[nSignals]   = aHit->GetParticleID();
	  det_VvvTrackSign[nSignals]     = 1;
	  nSignals++;
	}
      }   // end of "if (!signalAssigned)"
    }  // end of the for loop over the hits

    // Sort the signals according to the energy (in decreasing order)
    std::map<G4double,G4int> mySignalMapping;
    std::map<G4double,G4int>::iterator itt;
    for (G4int i=0; i<nSignals; i++) {
      mySignalMapping.insert ( std::pair<G4double,G4int>(-det_edep[i],i) );
    }

    // Write out the signals (sorted according to energy) to the musrRootOutput class:
    G4int jj=-1;
    for (itt=mySignalMapping.begin(); itt!=mySignalMapping.end(); itt++) {
      jj++;
      G4int ii = itt->second;
      myRootOutput->SetDetectorInfo(jj,det_ID[ii],det_VrtxParticleID[ii],det_edep[ii],
                                    det_edep_el[ii],det_edep_pos[ii],
				    det_edep_gam[ii],det_edep_mup[ii],det_nsteps[ii],det_length[ii],
				    det_time_start[ii],det_time_end[ii],det_x[ii],det_y[ii],det_z[ii],
				    det_kine[ii],det_VrtxKine[ii],det_VrtxX[ii],det_VrtxY[ii],det_VrtxZ[ii],
				    det_VrtxVolID[ii],det_VrtxProcID[ii],det_VrtxTrackID[ii] );

      if (boolIsVvvInfoRequested) {
	G4int oldTrackID = abs(det_VrtxTrackID[ii]);
	musrSteppingAction* myMusrSteppingAction = musrSteppingAction::GetInstance();
	G4int vvvParentTrackID= -1; G4int vvvPparticleID; G4double vvvKine; G4ThreeVector vvvPosition; G4String vvvLogVol; G4String vvvProcess;
	G4int vvvLogVolID=-999;
	G4bool oldTrackRetrievedOK;
	do {
	  if (vvvParentTrackID>-1) {oldTrackID=vvvParentTrackID;}
	  oldTrackRetrievedOK = myMusrSteppingAction->GetInfoAboutOldTrack(oldTrackID, 
				     vvvParentTrackID, vvvPparticleID, vvvKine, vvvPosition, vvvLogVol, vvvProcess);
	  if (oldTrackRetrievedOK) {
	    //	  G4cout<<"musrScintSD: Old Track seems to be achieved fine -> Lets save it to Root tree"<<G4endl;
	    vvvLogVolID=myRootOutput->ConvertVolumeToID(vvvLogVol);
	  }
	  else { 
	    G4cout<<" oldTrackRetrievedOK is false"<<G4endl; 
	    oldTrackID       = -999; 
	    vvvParentTrackID = -999;
	    vvvPparticleID   = -999;
	    vvvKine          = -999;
	    vvvPosition      = G4ThreeVector(-999,-999,-999);
	    vvvLogVol        = "undefined";
	    vvvProcess       = "undefined";
	  }  
	} while (oldTrackRetrievedOK && (vvvLogVolID==det_ID[ii]));

	if (oldTrackRetrievedOK) {
	  G4int vvvProcessID=myRootOutput->ConvertProcessToID(vvvProcess);
	  myRootOutput->SetDetectorInfoVvv(jj,vvvKine,vvvPosition.x(),vvvPosition.y(),vvvPosition.z(),
					 vvvLogVolID,vvvProcessID,oldTrackID*det_VvvTrackSign[ii],vvvPparticleID);
	}
      }   //end "boolIsVvvInfoRequested"

    }
  }   //end "if (NbHits>0)"

  // Analyse optical photons if they were produced
  if (musrParameters::boolG4OpticalPhotons) {
    if (!musrParameters::boolG4OpticalPhotonsUnprocess)  EndOfEvent_OptiacalPhotons();
  }
}  



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void musrScintSD::EndOfEvent_OptiacalPhotons() {
  //  G4double  kCarTolerance = G4GeometryTolerance::GetInstance()  ->GetSurfaceTolerance();
  //  G4cout<<"  DEBUG 10:  kCarTolerance="<<kCarTolerance<<G4endl;

  if (optHitMap.empty()) return;
  
  G4RunManager* fRunManager = G4RunManager::GetRunManager();
  G4int eeeventID = fRunManager->GetCurrentEvent()->GetEventID();

  for (optHitMapType::const_iterator it=optHitMap.begin() ; it != optHitMap.end(); it++ ) {
    G4bool boolStoreThisOPSAhist = false;
    //    G4bool boolStoreThisOPSAhistSUMMED = false;
    G4int  OPSA_detID= it->first;
    optHitDetectorMapType* optHitDetectorMap = it->second;

    // Check whether OPSA histograming of times of optical photon detection is required for this eventID.
    if (bool_multimapOfEventIDsForOPSAhistosEXISTS) {
      if (bool_StoreThisOPSAhistALL) boolStoreThisOPSAhist = true;
      else if (multimapOfEventIDsForOPSAhistos.find(eeeventID)!=multimapOfEventIDsForOPSAhistos.end()) {
	//  Now check whether the histogramming is required for the currently analysed detector.
	std::pair<multimapOfEventIDsForOPSAhistos_Type::iterator,multimapOfEventIDsForOPSAhistos_Type::iterator> retOPSAhist;
	multimapOfEventIDsForOPSAhistos_Type::iterator itOPSAhist;
	retOPSAhist = multimapOfEventIDsForOPSAhistos.equal_range(eeeventID);
	for (itOPSAhist = retOPSAhist.first; itOPSAhist!=retOPSAhist.second; itOPSAhist++) {
	  //  Store histograms if the second parameters of eventsForOPSAhistos (i.e. detector ID) 
          //      is set to 0 or corresponds to the currently analysed sensitive detector.
	  if ( (itOPSAhist->second == 0) || (itOPSAhist->second == OPSA_detID) ) boolStoreThisOPSAhist = true;
	}
      }
      //      if (multimapOfEventIDsForOPSAhistos.find(-1)!=multimapOfEventIDsForOPSAhistos.end()) {
      //	// The user requires to store OPSA timing histograms summed up for all events together
      //	boolStoreThisOPSAhistSUMMED = true;
      //      }
    }

    if (optHitDetectorMap->empty()) continue;

    // Do APD cell variation time if requested.  Even if not requested, generate the random numbers in order to
    // have reproducible simulation.
    //    for (optHitDetectorMapType::iterator it2 = optHitDetectorMap->begin(); it2 != optHitDetectorMap->end(); it2++ ) {
    //      G4double APDcellsTimeVariation = G4RandGauss::shoot(0,APDcellsTimeVariationSigma);
    //      if (APDcellsTimeVariationRequested) {
    //	G4int  APDcellID = it2->second;
    //	G4double tmpTime = it2->first + APDcellsTimeVariation;
    //	optHitDetectorMap->erase(it2);
    //	optHitDetectorMap->insert(std::pair<G4double,G4int>(tmpTime,APDcellID));
    //      }
    //    }

    // Simulate cross talk if requested.  Even if not requested, generate the random numbers in order to
    // have reproducible simulation.
    for (optHitDetectorMapType::const_iterator it2 = optHitDetectorMap->begin(); it2 != optHitDetectorMap->end(); it2++ ) {
      G4double rand = G4UniformRand();
      if (APDcrossTalkRequested) {
	if (rand<APDcrossTalkProb) {
	  G4int  APDcellID = (it2->second).GetAPDcellID();
	  G4double tmpTime = it2->first;
	  if (!APDcellsEffectRequested) FireAPDcell(optHitDetectorMap,APDcellID,tmpTime,0);
	  else {
	    if (rand<(0.5*APDcrossTalkProb)) APDcellID--;
	    else                             APDcellID++;
	    if ((APDcellID<0)||(APDcellID>=APDcell_nx*APDcell_ny*APDcell_nz)) continue;  // the cross-talk cell is outside range	  
	    FireAPDcell(optHitDetectorMap,APDcellID,tmpTime,0);
	  }
	}
      }
    }

    optHitDetectorMapType::const_iterator it2_START = optHitDetectorMap->begin();
    optHitDetectorMapType::const_iterator it2_STOP  = it2_START;
    optHitDetectorMapType::const_iterator it2_LAST = optHitDetectorMap->end();  it2_LAST--;
    Double_t time = -1000, lastTime = -1000;
    G4int    OPSA_nPhot       = 0;
    G4int    OPSA_nPhotPrim   = 0;
    G4double OPSA_timeFirst   = -1000000;
    G4double OPSA_timeSecond  = -1000000;
    G4double OPSA_timeThird   = -1000000;
    G4double OPSA_timeA       = -1000000;
    G4double OPSA_timeB       = -1000000;
    G4double OPSA_timeC       = -1000000;
    G4double OPSA_timeD       = -1000000;
    G4double OPSA_timeMean    = -1000000;
    G4double OPSA_timeLast    = -1000000;
    G4double OPSA_CFD_time    = -1000000;
    G4double OPSA_CFD_ampl    = -1000;
    G4int    iHistNr          = -1;
    G4int    iHistNrSUM       = -1;
    //    TF1* poiss = NULL;
    for (optHitDetectorMapType::const_iterator it2 = optHitDetectorMap->begin(); it2 != optHitDetectorMap->end(); it2++ ) {
      //delete      G4cout<<"KAMIL: Nr of optHitDetectorMap->size()="<<optHitDetectorMap->size()<<G4endl;
      OPSA_nPhot++;
      lastTime = time;
      time = it2->first;
      if (OPSA_nPhot==1) lastTime=time;  // First photon does not have a proper "lastTime" defined
      if ( (it2==it2_LAST) || ((time-lastTime) > OPSA_signalSeparationTime)) {
	// The iterator it2 reached last element of the map optHitDetectorMap       or
	// the time difference between two subsequently detected photons is too big ==> divide the signal into more signals
	it2_STOP = it2;
	if (it2==it2_LAST) it2_STOP++;          // if we are at the end of optHitDetectorMap, 
                                                //    the it2_STOP should point just behind the last map element
	else OPSA_nPhot--;                      // if we split the optHitDetectorMap, however, the latest optical photon
	                                        //    should have not be counted, because it belong to the next signal.
	optHitDetectorMapType optHitDetectorSUBmap(it2_START, it2_STOP);
	it2_START = it2_STOP;
	if (OPSA_nPhot >= OPSA_minNrOfDetectedPhotons) {   // ignore hits with too low number of detected photons
	  G4double OPSA_f_nPhot = OPSA_nPhot;
	  G4int NA = int (OPSA_fracA * OPSA_f_nPhot + 0.5);  if (NA<=0) NA=1;
	  G4int NB = int (OPSA_fracB * OPSA_f_nPhot + 0.5);  if (NB<=0) NB=1;

	  Int_t nP=0;
	  // Define OPSA histograms if required for this event
	  if ((boolStoreThisOPSAhist)||(bool_pulseShapeExists)) {
	    iHistNr++;
	    char nameHist[200];  sprintf(nameHist,"OPSAhist_%d_%d_%d",eeeventID,OPSA_detID,iHistNr);
	    char nameHistTitle[200]; sprintf(nameHistTitle,"OPSAhist_%d_%d_%d;time (ns);Nr of photons",eeeventID,OPSA_detID,iHistNr);
            if (verboseLevel>1) G4cout<<"VERBOSE 2 : creating a new TH1D for OPSAhisto" <<"\n";
	    OPSAhisto = new TH1D(nameHist, nameHistTitle, OPSAhistoNbin, OPSAhistoMin, OPSAhistoMax);
	    //	    poiss = new TF1("poiss",poissonf,0.,.5,2); // x in [0;300], 2
	    //	    poiss->SetParameter(0,1);
	    //	    poiss->SetParameter(1,1);
	    if (bool_pulseShapeExists) {
	      sprintf(nameHist,"OPSAshape_%d_%d_%d",eeeventID,OPSA_detID,iHistNr);
	      sprintf(nameHistTitle,"OPSAshape_%d_%d_%d;time (ns);Pulse signal",eeeventID,OPSA_detID,iHistNr);
              if (verboseLevel>1) G4cout<<"VERBOSE 2 : creating a new TH1D for OPSAshape" <<"\n";
	      OPSAshape = new TH1D(nameHist, nameHistTitle, OPSAhistoNbin, OPSAhistoMin, OPSAhistoMax);
	      sprintf(nameHist,"OPSA_CFD_%d_%d_%d",eeeventID,OPSA_detID,iHistNr);
	      sprintf(nameHistTitle,"OPSA_CFD_%d_%d_%d;time (ns);CFD signal",eeeventID,OPSA_detID,iHistNr);
              if (verboseLevel>1) G4cout<<"VERBOSE 2 : creating a new TH1D for OPSA_CFD" <<"\n";
	      OPSA_CFD = new TH1D(nameHist, nameHistTitle, OPSAhistoNbin, OPSAhistoMin, OPSAhistoMax);
	    }
	  }
	  if (bool_StoreThisOPSAhistSUMMED) {
	    iHistNrSUM++;
	    char nameHist[200];  sprintf(nameHist,"OPSAhistSUM_%d_%d",OPSA_detID,iHistNrSUM);
	    char nameHist0[200];  sprintf(nameHist0,"OPSAhistSUM0_%d_%d",OPSA_detID,iHistNrSUM);
	    if (mapOfOPSAsumHistograms.find(nameHist) != mapOfOPSAsumHistograms.end()) {
	      OPSAhistoSUM  = mapOfOPSAsumHistograms[nameHist];
	      OPSAhistoSUM0 = mapOfOPSAsum0Histograms[nameHist0];
	      //      G4cout<<" OPSAhistoSUM histogram found:"<<OPSAhistoSUM<<G4endl;
	    }
	    else {  // create histogram because it does not exist
	      char nameHistTitle[200]; sprintf(nameHistTitle,"OPSAhistSUM_%d_%d;time (ns);Nr of photons",OPSA_detID,iHistNrSUM);
	      char nameHistTitle0[200]; sprintf(nameHistTitle0,"OPSAhistSUM0_%d_%d;time (ns);Nr of photons",OPSA_detID,iHistNrSUM);
              if (verboseLevel>1) G4cout<<"VERBOSE 2 : creating new TH1Ds for OPSAhistoSUM,0" <<"\n";
	      OPSAhistoSUM = new TH1D(nameHist, nameHistTitle, OPSAhistoNbin, OPSAhistoMin, OPSAhistoMax);
	      OPSAhistoSUM0= new TH1D(nameHist0,nameHistTitle0,OPSAhistoNbin, OPSAhistoMin, OPSAhistoMax);
	      mapOfOPSAsumHistograms.insert(std::pair<std::string,TH1D*>(nameHist,OPSAhistoSUM));
	      mapOfOPSAsum0Histograms.insert(std::pair<std::string,TH1D*>(nameHist0,OPSAhistoSUM0));
	      //      G4cout<<" New OPSAhistoSUM histogram created:"<<OPSAhistoSUM<<G4endl;
	    }
	  }
	  
	  for (optHitDetectorMapType::const_iterator it3 = optHitDetectorSUBmap.begin(); it3 != optHitDetectorSUBmap.end(); it3++) {
	    nP++;
	    G4double timePhot = it3->first;
            OPSA_nPhotPrim   += (it3->second).GetAPDcellNphot();
	    //	    if (APDcellsTimeVariationRequested)  timePhot += G4RandGauss::shoot(0,APDcellsTimeVariationSigma); // Shifted above
	    if (nP==1)          OPSA_timeFirst = timePhot;
	    if (nP==2)          OPSA_timeSecond= timePhot;
	    if (nP==3)          OPSA_timeThird = timePhot;
	    if (nP==NA)         OPSA_timeA     = timePhot;
	    if (nP==NB)         OPSA_timeB     = timePhot;
	    if (nP==OPSA_nPhot) OPSA_timeLast  = timePhot;
	    if ((boolStoreThisOPSAhist)||(bool_pulseShapeExists)) {
	      G4double tt = timePhot-OPSA_timeFirst+0.00000000001;
	      OPSAhisto->Fill(tt);
	      if (bool_pulseShapeExists) {  // fill contribution of this detected photon to the overall signal shape
		Int_t iBin = int(tt/OPSAhistoBinWidth);
		Double_t rest = tt - iBin*OPSAhistoBinWidth;
		for (Int_t iBinNew=0; iBinNew<OPSAhistoNbin-iBin; iBinNew++) { // loop over all bins of the shape histogram
		  Int_t iPS = int( (double(iBinNew))*OPSAhistoBinWidth1000+rest*1000 );  // iPS is time in picoseconds 
		                                                                // and also the index of the shape array
		  if (iPS>=iPSmax-1) break;  // out of the pulse shape info ==> leave this loop
		  OPSAshape->Fill((iBin+iBinNew+0.5)*OPSAhistoBinWidth,pulseShapeArray[iPS]);
		}
	      }
	    }
	    if (bool_StoreThisOPSAhistSUMMED) {
	      OPSAhistoSUM ->Fill(timePhot-OPSA_timeFirst+0.00000000001);
	      OPSAhistoSUM0->Fill(timePhot+0.00000000001);
	    }
	  }
	  
	  //	    OPSAhisto->Fit(poiss,"Q");
	  if (boolStoreThisOPSAhist) OPSAhisto->Write();
	  OPSA_timeMean = OPSAhisto->GetMean() + OPSA_timeFirst;

	  // if required, convert the histogram with photon times to histogram of electronic pulse shape
	  G4double timeCFDarray[1000];
	  G4double OPSA_timeC1[41];
	  if (bool_pulseShapeExists) {
	    //	    for (Int_t iBin=1; iBin<=OPSAhistoNbin; iBin++) { // loop over all bins of photon histogram
	    //	      Double_t photonCount = OPSAhisto ->GetBinContent(iBin);
	    //	    //	      if (photonCount==0.) break;  //nothing to be done for bins without any photons
	    // here was error - instead of continue I used break !!!  
	    //        if (photonCount==0.) continue;  //nothing to be done for bins without any photons
	    //	      for (Int_t iBinNew=0; iBinNew<(OPSAhistoNbin-iBin); iBinNew++) { // loop over all bins from the actual one to the end
	    //		Int_t iPS = int( (double(iBinNew)+0.5)*OPSAhistoBinWidth1000 );
	    //		if (iPS>=iPSmax-1) break;  // out of the pulse shape info ==> leave this loop
	    //		//		  OPSAshape->AddBinContent(iBin+iBinNew,photonCount*pulseShapeArray[iPS]*OPSAhistoBinWidth);  // iPS corresponds to time in picoseconds
	    //		OPSAshape->Fill((iBin+iBinNew-0.5)*OPSAhistoBinWidth,photonCount*pulseShapeArray[iPS]*OPSAhistoBinWidth);
	    //	      }
	    //	    }
	    if (boolStoreThisOPSAhist) OPSAshape -> Write();
	    

	    //  If requested, store CFD times for various CFD delays and amplitudes of the reverted signal
	    if (musrRootOutput::store_odet_timeCFDarray) {
	      G4double OPSA_CFD_a1_saveThis    = OPSA_CFD_a1;
	      G4double OPSA_CFD_delay_saveThis = OPSA_CFD_delay;
	      const int nDelay = 13;
	      G4double delay[nDelay]={0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0};
	      const int nAmpl1 = 5;
	      G4double ampl1[nAmpl1] = {-0.05, -0.1, -0.2, -0.3, -0.4};
	      for (Int_t kk=0; kk<nDelay; kk++) {
		for (Int_t ll=0; ll<nAmpl1; ll++) {
		  OPSA_CFD_a1    = ampl1[ll];
		  OPSA_CFD_delay = delay[kk];
		  FindCFDtime(OPSA_CFD_time, OPSA_CFD_ampl, OPSA_timeFirst);
		  int index = (ll+1)*100+kk;
		  if (index<1000) timeCFDarray[(ll+1)*100+kk] = OPSA_CFD_time;
		  else {G4cout<<"FATALNI ERROR in musrScintSD by calculating (ll+1)*100+kk"<<G4endl;}
		  //		G4cout<<"   delay= "<<OPSA_CFD_delay<<"  ampl1= "<<OPSA_CFD_a1
		  //		      <<"  OPSA_CFD_time= "<<OPSA_CFD_time<<"   OPSA_CFD_ampl= "<<OPSA_CFD_ampl<<G4endl;
		}
	      }
	      OPSA_CFD_a1    = OPSA_CFD_a1_saveThis;
	      OPSA_CFD_delay = OPSA_CFD_delay_saveThis;
	    }

	    // Now fill the histogram with the CFD signal
	    FindCFDtime(OPSA_CFD_time, OPSA_CFD_ampl, OPSA_timeFirst);
	    //	    G4cout<<"OPSA_CFD_time = "<<OPSA_CFD_time<<"   OPSA_CFD_ampl = "<<OPSA_CFD_ampl<<G4endl;
	    if (boolStoreThisOPSAhist) OPSA_CFD -> Write();

	    // Find the timeC and timeD from the shape histogram
	    Double_t yValue=-1000;
	    G4bool   OPSA_C_time_found = false;
	    G4bool   OPSA_D_time_found = false;
	    Double_t D_threshold = OPSA_D_threshold * (OPSAshape->GetMaximum());  // covert relative "OPSA_D_threshold" into 
	                                                                          // the absolute threshold using signal amplitude
	    for (Int_t iBin=1; iBin<=OPSAhistoNbin; iBin++) { // loop over all bins of CFD histogram
	      Double_t xValue = (iBin-0.5)*OPSAhistoBinWidth;
	      Double_t oldYvalue = yValue;
	      yValue          = OPSAshape->GetBinContent(iBin);
	      if ( (yValue>=OPSA_C_threshold) && (!OPSA_C_time_found) ) {  // signal just moved above the threshold
		OPSA_C_time_found = true;
		OPSA_timeC = xValue - (yValue-OPSA_C_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		OPSA_timeC += OPSA_timeFirst;
	      }
	      
	      if ( (yValue>=D_threshold) && (!OPSA_D_time_found) ) {  // signal just moved above the threshold
		OPSA_D_time_found = true;
		OPSA_timeD = xValue - (yValue-D_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		OPSA_timeD += OPSA_timeFirst;
	      }
	      
	      if (OPSA_C_time_found && OPSA_D_time_found) break;  // times C and D found, no need to continue looping.
	    }


	    if (musrRootOutput::store_odet_timeC1) {
	      Double_t yValue=-1000;
	      for (Int_t i=1; i<=40; i++) {OPSA_timeC1[i] = -1000000.;}	      
	      G4bool  OPSA_C1_time_found = false;   G4bool  OPSA_C2_time_found = false;
	      G4bool  OPSA_C3_time_found = false;   G4bool  OPSA_C4_time_found = false;
	      G4bool  OPSA_C5_time_found = false;   G4bool  OPSA_C6_time_found = false;
	      G4bool  OPSA_C7_time_found = false;   G4bool  OPSA_C8_time_found = false;
	      G4bool  OPSA_C9_time_found = false;   G4bool  OPSA_C10_time_found = false;
	      G4bool  OPSA_C11_time_found = false;  G4bool  OPSA_C12_time_found = false;
	      G4bool  OPSA_C13_time_found = false;  G4bool  OPSA_C14_time_found = false;
	      G4bool  OPSA_C15_time_found = false;  G4bool  OPSA_C16_time_found = false;
	      G4bool  OPSA_C17_time_found = false;  G4bool  OPSA_C18_time_found = false;
	      G4bool  OPSA_C19_time_found = false;  G4bool  OPSA_C20_time_found = false;
	      G4bool  OPSA_C21_time_found = false;  G4bool  OPSA_C22_time_found = false;
	      G4bool  OPSA_C23_time_found = false;  G4bool  OPSA_C24_time_found = false;
	      G4bool  OPSA_C25_time_found = false;  G4bool  OPSA_C26_time_found = false;
	      G4bool  OPSA_C27_time_found = false;  G4bool  OPSA_C28_time_found = false;
	      G4bool  OPSA_C29_time_found = false;  G4bool  OPSA_C30_time_found = false;
	      G4bool  OPSA_C31_time_found = false;  G4bool  OPSA_C32_time_found = false;
	      G4bool  OPSA_C33_time_found = false;  G4bool  OPSA_C34_time_found = false;
	      G4bool  OPSA_C35_time_found = false;  G4bool  OPSA_C36_time_found = false;
	      G4bool  OPSA_C37_time_found = false;  G4bool  OPSA_C38_time_found = false;
	      G4bool  OPSA_C39_time_found = false;  G4bool  OPSA_C40_time_found = false;
	      Double_t OPSAshape_histoMax = OPSAshape->GetMaximum();
	      for (Int_t iBin=1; iBin<=OPSAhistoNbin; iBin++) { // loop over all bins of CFD histogram
		Double_t xValue = (iBin-0.5)*OPSAhistoBinWidth;
		Double_t oldYvalue = yValue;
		yValue          = OPSAshape->GetBinContent(iBin);
		if ( (yValue>=OPSA_C1_threshold) && (!OPSA_C1_time_found) ) {  // signal just moved above the threshold
		  OPSA_C1_time_found = true;
		  OPSA_timeC1[1] = OPSA_timeFirst + xValue - (yValue-OPSA_C1_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C2_threshold) && (!OPSA_C2_time_found) ) {  // signal just moved above the threshold
		  OPSA_C2_time_found = true;
		  OPSA_timeC1[2] = OPSA_timeFirst + xValue - (yValue-OPSA_C2_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C3_threshold) && (!OPSA_C3_time_found) ) {  // signal just moved above the threshold
		  OPSA_C3_time_found = true;
		  OPSA_timeC1[3] = OPSA_timeFirst + xValue - (yValue-OPSA_C3_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C4_threshold) && (!OPSA_C4_time_found) ) {  // signal just moved above the threshold
		  OPSA_C4_time_found = true;
		  OPSA_timeC1[4] = OPSA_timeFirst + xValue - (yValue-OPSA_C4_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C5_threshold) && (!OPSA_C5_time_found) ) {  // signal just moved above the threshold
		  OPSA_C5_time_found = true;
		  OPSA_timeC1[5] = OPSA_timeFirst + xValue - (yValue-OPSA_C5_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C6_threshold) && (!OPSA_C6_time_found) ) {  // signal just moved above the threshold
		  OPSA_C6_time_found = true;
		  OPSA_timeC1[6] = OPSA_timeFirst + xValue - (yValue-OPSA_C6_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C7_threshold) && (!OPSA_C7_time_found) ) {  // signal just moved above the threshold
		  OPSA_C7_time_found = true;
		  OPSA_timeC1[7] = OPSA_timeFirst + xValue - (yValue-OPSA_C7_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C8_threshold) && (!OPSA_C8_time_found) ) {  // signal just moved above the threshold
		  OPSA_C8_time_found = true;
		  OPSA_timeC1[8] = OPSA_timeFirst + xValue - (yValue-OPSA_C8_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C9_threshold) && (!OPSA_C9_time_found) ) {  // signal just moved above the threshold
		  OPSA_C9_time_found = true;
		  OPSA_timeC1[9] = OPSA_timeFirst + xValue - (yValue-OPSA_C9_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C10_threshold) && (!OPSA_C10_time_found) ) {  // signal just moved above the threshold
		  OPSA_C10_time_found = true;
		  OPSA_timeC1[10] = OPSA_timeFirst + xValue - (yValue-OPSA_C10_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C11_threshold) && (!OPSA_C11_time_found) ) {  // signal just moved above the threshold
		  OPSA_C11_time_found = true;
		  OPSA_timeC1[11] = OPSA_timeFirst + xValue - (yValue-OPSA_C11_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C12_threshold) && (!OPSA_C12_time_found) ) {  // signal just moved above the threshold
		  OPSA_C12_time_found = true;
		  OPSA_timeC1[12] = OPSA_timeFirst + xValue - (yValue-OPSA_C12_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C13_threshold) && (!OPSA_C13_time_found) ) {  // signal just moved above the threshold
		  OPSA_C13_time_found = true;
		  OPSA_timeC1[13] = OPSA_timeFirst + xValue - (yValue-OPSA_C13_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C14_threshold) && (!OPSA_C14_time_found) ) {  // signal just moved above the threshold
		  OPSA_C14_time_found = true;
		  OPSA_timeC1[14] = OPSA_timeFirst + xValue - (yValue-OPSA_C14_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C15_threshold) && (!OPSA_C15_time_found) ) {  // signal just moved above the threshold
		  OPSA_C15_time_found = true;
		  OPSA_timeC1[15] = OPSA_timeFirst + xValue - (yValue-OPSA_C15_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C16_threshold) && (!OPSA_C16_time_found) ) {  // signal just moved above the threshold
		  OPSA_C16_time_found = true;
		  OPSA_timeC1[16] = OPSA_timeFirst + xValue - (yValue-OPSA_C16_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C17_threshold) && (!OPSA_C17_time_found) ) {  // signal just moved above the threshold
		  OPSA_C17_time_found = true;
		  OPSA_timeC1[17] = OPSA_timeFirst + xValue - (yValue-OPSA_C17_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C18_threshold) && (!OPSA_C18_time_found) ) {  // signal just moved above the threshold
		  OPSA_C18_time_found = true;
		  OPSA_timeC1[18] = OPSA_timeFirst + xValue - (yValue-OPSA_C18_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C19_threshold) && (!OPSA_C19_time_found) ) {  // signal just moved above the threshold
		  OPSA_C19_time_found = true;
		  OPSA_timeC1[19] = OPSA_timeFirst + xValue - (yValue-OPSA_C19_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C20_threshold) && (!OPSA_C20_time_found) ) {  // signal just moved above the threshold
		  OPSA_C20_time_found = true;
		  OPSA_timeC1[20] = OPSA_timeFirst + xValue - (yValue-OPSA_C20_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C21_threshold) && (!OPSA_C21_time_found) ) {  // signal just moved above the threshold
		  OPSA_C21_time_found = true;
		  OPSA_timeC1[21] = OPSA_timeFirst + xValue - (yValue-OPSA_C21_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C22_threshold) && (!OPSA_C22_time_found) ) {  // signal just moved above the threshold
		  OPSA_C22_time_found = true;
		  OPSA_timeC1[22] = OPSA_timeFirst + xValue - (yValue-OPSA_C22_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C23_threshold) && (!OPSA_C23_time_found) ) {  // signal just moved above the threshold
		  OPSA_C23_time_found = true;
		  OPSA_timeC1[23] = OPSA_timeFirst + xValue - (yValue-OPSA_C23_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C24_threshold) && (!OPSA_C24_time_found) ) {  // signal just moved above the threshold
		  OPSA_C24_time_found = true;
		  OPSA_timeC1[24] = OPSA_timeFirst + xValue - (yValue-OPSA_C24_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C25_threshold) && (!OPSA_C25_time_found) ) {  // signal just moved above the threshold
		  OPSA_C25_time_found = true;
		  OPSA_timeC1[25] = OPSA_timeFirst + xValue - (yValue-OPSA_C25_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C26_threshold) && (!OPSA_C26_time_found) ) {  // signal just moved above the threshold
		  OPSA_C26_time_found = true;
		  OPSA_timeC1[26] = OPSA_timeFirst + xValue - (yValue-OPSA_C26_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C27_threshold) && (!OPSA_C27_time_found) ) {  // signal just moved above the threshold
		  OPSA_C27_time_found = true;
		  OPSA_timeC1[27] = OPSA_timeFirst + xValue - (yValue-OPSA_C27_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C28_threshold) && (!OPSA_C28_time_found) ) {  // signal just moved above the threshold
		  OPSA_C28_time_found = true;
		  OPSA_timeC1[28] = OPSA_timeFirst + xValue - (yValue-OPSA_C28_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C29_threshold) && (!OPSA_C29_time_found) ) {  // signal just moved above the threshold
		  OPSA_C29_time_found = true;
		  OPSA_timeC1[29] = OPSA_timeFirst + xValue - (yValue-OPSA_C29_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C30_threshold) && (!OPSA_C30_time_found) ) {  // signal just moved above the threshold
		  OPSA_C30_time_found = true;
		  OPSA_timeC1[30] = OPSA_timeFirst + xValue - (yValue-OPSA_C30_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C31_threshold) && (!OPSA_C31_time_found) ) {  // signal just moved above the threshold
		  OPSA_C31_time_found = true;
		  OPSA_timeC1[31] = OPSA_timeFirst + xValue - (yValue-OPSA_C31_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C32_threshold) && (!OPSA_C32_time_found) ) {  // signal just moved above the threshold
		  OPSA_C32_time_found = true;
		  OPSA_timeC1[32] = OPSA_timeFirst + xValue - (yValue-OPSA_C32_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C33_threshold) && (!OPSA_C33_time_found) ) {  // signal just moved above the threshold
		  OPSA_C33_time_found = true;
		  OPSA_timeC1[33] = OPSA_timeFirst + xValue - (yValue-OPSA_C33_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C34_threshold) && (!OPSA_C34_time_found) ) {  // signal just moved above the threshold
		  OPSA_C34_time_found = true;
		  OPSA_timeC1[34] = OPSA_timeFirst + xValue - (yValue-OPSA_C34_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C35_threshold) && (!OPSA_C35_time_found) ) {  // signal just moved above the threshold
		  OPSA_C35_time_found = true;
		  OPSA_timeC1[35] = OPSA_timeFirst + xValue - (yValue-OPSA_C35_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C36_threshold) && (!OPSA_C36_time_found) ) {  // signal just moved above the threshold
		  OPSA_C36_time_found = true;
		  OPSA_timeC1[36] = OPSA_timeFirst + xValue - (yValue-OPSA_C36_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C37_threshold) && (!OPSA_C37_time_found) ) {  // signal just moved above the threshold
		  OPSA_C37_time_found = true;
		  OPSA_timeC1[37] = OPSA_timeFirst + xValue - (yValue-OPSA_C37_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C38_threshold) && (!OPSA_C38_time_found) ) {  // signal just moved above the threshold
		  OPSA_C38_time_found = true;
		  OPSA_timeC1[38] = OPSA_timeFirst + xValue - (yValue-OPSA_C38_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C39_threshold) && (!OPSA_C39_time_found) ) {  // signal just moved above the threshold
		  OPSA_C39_time_found = true;
		  OPSA_timeC1[39] = OPSA_timeFirst + xValue - (yValue-OPSA_C39_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}
		if ( (yValue>=OPSA_C40_threshold) && (!OPSA_C40_time_found) ) {  // signal just moved above the threshold
		  OPSA_C40_time_found = true;
		  OPSA_timeC1[40] = OPSA_timeFirst + xValue - (yValue-OPSA_C40_threshold)/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
		}

		if (OPSA_C40_time_found) {
		  if (OPSA_C1_time_found  && OPSA_C2_time_found  && OPSA_C3_time_found  && OPSA_C4_time_found && 
		      OPSA_C5_time_found  && OPSA_C6_time_found  && OPSA_C7_time_found  && OPSA_C8_time_found && 
		      OPSA_C9_time_found  && OPSA_C10_time_found && OPSA_C11_time_found && OPSA_C12_time_found && 
		      OPSA_C13_time_found && OPSA_C14_time_found && OPSA_C15_time_found && OPSA_C16_time_found && 
		      OPSA_C17_time_found && OPSA_C18_time_found && OPSA_C19_time_found && OPSA_C20_time_found &&
		      OPSA_C21_time_found && OPSA_C22_time_found && OPSA_C23_time_found && OPSA_C24_time_found && 
		      OPSA_C25_time_found && OPSA_C26_time_found && OPSA_C27_time_found && OPSA_C28_time_found &&
		      OPSA_C29_time_found && OPSA_C30_time_found && OPSA_C31_time_found && OPSA_C32_time_found && 
		      OPSA_C33_time_found && OPSA_C34_time_found && OPSA_C35_time_found && OPSA_C36_time_found && 
		      OPSA_C37_time_found && OPSA_C38_time_found && OPSA_C39_time_found) break; // all times found
		}
		
		if (yValue==OPSAshape_histoMax) break;  // no need for further searching, because maximum was achieved.
	      }

	      //	      myRootOutput -> SetTimeC1SpecialInfo(OPSA_timeC1);
	    }

	    //	    G4cout<<" OPSA_CFD_time = "<< OPSA_CFD_time<<",   OPSA_CFD_ampl="<<OPSA_CFD_ampl
	    //		  <<",  OPSA_timeFirst="<<OPSA_timeFirst<<",  OPSA_timeC="<<OPSA_timeC <<",  OPSA_timeD="<<OPSA_timeD<<G4endl;
	  } // end if (pulseShapeExists)

	  // Delete the histograms from memory if they were created
	  if ((boolStoreThisOPSAhist)||(bool_pulseShapeExists)) {
            if (verboseLevel>1) G4cout<<"VERBOSE 2 : deleting OPSAhisto" <<"\n";
	    delete OPSAhisto;
	    if (bool_pulseShapeExists) {
              if (verboseLevel>1) G4cout<<"VERBOSE 2 : deleting OPSAshape and CFD" <<"\n";
	      delete OPSAshape;
	      delete OPSA_CFD;
	    }
	  }

          if (verboseLevel>1) G4cout<<"VERBOSE 2 : creating new signalInfo" <<"\n";
	  signalInfo* mySignalInfo = new signalInfo(OPSA_detID,OPSA_nPhot,OPSA_nPhotPrim,OPSA_timeFirst,OPSA_timeSecond,OPSA_timeThird,
						    OPSA_timeA,OPSA_timeB,OPSA_timeC,OPSA_timeD,OPSA_timeMean,OPSA_timeLast,
						    OPSA_CFD_time,OPSA_CFD_ampl,timeCFDarray,OPSA_timeC1);
	  OPSA_signal_Map.insert(std::pair<G4int,signalInfo*>(OPSA_nPhot,mySignalInfo) );
	}
	OPSA_nPhot       = 0;
	OPSA_nPhotPrim   = 0;
	OPSA_timeFirst   = -1000000;
	OPSA_timeSecond  = -1000000;
	OPSA_timeThird   = -1000000;
	OPSA_timeA       = -1000000;
	OPSA_timeB       = -1000000;
	OPSA_timeC       = -1000000;
	OPSA_timeD       = -1000000;
	OPSA_timeMean    = -1000000;
	OPSA_timeLast    = -1000000;
	OPSA_CFD_time    = -1000000;
	OPSA_CFD_ampl    = -1000;
      }
    }
  }

  //  // Now delete all optHitDetectorMap*
  //  for (optHitMapType::const_iterator it=optHitMap.begin() ; it != optHitMap.end(); it++ ) {
  //    delete (it->second);
  //  }
  //  optHitMap.clear();


  //  Now store the information of all mySignalInfo*  to rootOutputFile
  G4int nn=0;
  for (OPSA_signal_MapType::reverse_iterator ritOPSA = OPSA_signal_Map.rbegin(); ritOPSA != OPSA_signal_Map.rend(); ritOPSA++) {
    (ritOPSA->second) -> transferDataToRoot(myRootOutput,nn);
    nn++;
  }
  //  Now delete all mySignalInfo*  from  OPSA_signal_Map
  for (OPSA_signal_MapType::const_iterator itOPSA = OPSA_signal_Map.begin(); itOPSA != OPSA_signal_Map.end(); itOPSA++) {
    if (verboseLevel>1) G4cout<<"VERBOSE 2 : deleting itOPSA->second" <<"\n";
    delete (itOPSA->second);
  }
  OPSA_signal_Map.clear();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void musrScintSD::EndOfRun() {
  if (musrParameters::boolG4OpticalPhotons) {
    if (!mapOfOPSAsumHistograms.empty()) {
      for (mapOfOPSAsumHistograms_Type::const_iterator it = mapOfOPSAsumHistograms.begin(); it!=mapOfOPSAsumHistograms.end(); it++) {
	(it->second) -> Write();
      }
    }
    if (!mapOfOPSAsum0Histograms.empty()) {
      for (mapOfOPSAsumHistograms_Type::const_iterator it = mapOfOPSAsum0Histograms.begin(); it!=mapOfOPSAsum0Histograms.end(); it++) {
	(it->second) -> Write();
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void musrScintSD::ReadInPulseShapeArray(const char* filename) {
  bool_pulseShapeExists = true;
  std::ifstream file( filename );
  //  FILE* file = fopen(filename,"r");
  if (!(file.is_open())) {
  //  if (file==NULL) {
    G4cout << "musrScintSD::ReadInPulseShapeArray:  requested pulse shape data file \""<< filename << "\" not opened (found) !!!"<<G4endl;
    musrErrorMessage::GetInstance()->musrError(FATAL,"musrScintSD: pulse shape data file not found !",false);
  }
  do {
    file.ignore(256, '\n');
  } while (file.peek() == '%');

  G4int i = -1;
  G4int idummy;
  G4double ampl;
  while (!file.eof()) {
  //  while (!feof(fSteeringFile)) {
    //    if (file.eof()) break;
    file >> idummy >> ampl;
    if (file.eof()) break;
    i++;

    if (i!=idummy) { 
      G4cout<<"musrScintSD::ReadInPulseShapeArray: pulse shape data file corrupted! i<>idummy  (i="<<i<<", idummy="<<idummy<<")"<<G4endl;
      musrErrorMessage::GetInstance()->musrError(FATAL,"musrScintSD: pulse shape data file corrupted !",false);
    }
    if (i>10000) {
      G4cout<<"musrScintSD::ReadInPulseShapeArray: pulse shape data file can not handle such big array ! i>10000  (i="<<i<<")"<<G4endl;
      musrErrorMessage::GetInstance()->musrError(FATAL,"musrScintSD: pulse shape data file problem !",false);
    }
    pulseShapeArray[i] = ampl;
  }

  iPSmax = i+1;
  file.close( );

  //  for (int j=0; j<iPSmax; j++) {G4cout<<j<<"  "<< pulseShapeArray[j]<<G4endl;}
  //  G4cout<<G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int musrScintSD::FindAPDcellID(G4Step* aStep) {
    // TEST:
  G4StepPoint* postPoint = aStep->GetPostStepPoint();
  G4ThreeVector postStepPoint = postPoint->GetPosition();
  const G4AffineTransform* trans = postPoint->GetTouchable()->GetHistory()->GetPtrTopTransform();
  G4ThreeVector postStepPointLocal = trans->TransformPoint(postStepPoint);
  //  G4cout<<"musrScintSD::FindAPDcellID:  postStepPointLocal="<<postStepPointLocal<<G4endl;
  G4int ix = (APDcell_nx>1) ? int ((postStepPointLocal.x() / APDcell_ax + APDcell_nx/2.)) : 0;
  G4int iy = (APDcell_ny>1) ? int ((postStepPointLocal.y() / APDcell_ay + APDcell_ny/2.)) : 0;
  G4int iz = (APDcell_nz>1) ? int ((postStepPointLocal.z() / APDcell_az + APDcell_nz/2.)) : 0;
  G4int in = ix *  APDcell_ny * APDcell_nz  +  iy * APDcell_nz + iz;
  //  G4cout<<"musrScintSD::FindAPDcellID:  in = "<<in<<",   ix = "<<ix<<",  iy = "<<iy<<",  iz = "<<iz<<G4endl;
  if ((ix<0)||(ix>APDcell_nx)) musrErrorMessage::GetInstance()->musrError(FATAL,
		"musrScintSD::FindAPDcellID: APD cell out of range (ix). Wrong dimensions of APD in \"/musr/command OPSA APDcells\" command?",false);
  if ((iy<0)||(iy>APDcell_ny)) musrErrorMessage::GetInstance()->musrError(FATAL,
		"musrScintSD::FindAPDcellID: APD cell out of range (iy). Wrong dimensions of APD in \"/musr/command OPSA APDcells\" command?",false);
  if ((iz<0)||(iz>APDcell_nz)) musrErrorMessage::GetInstance()->musrError(FATAL,
		"musrScintSD::FindAPDcellID: APD cell out of range (iz). Wrong dimensions of APD in \"/musr/command OPSA APDcells\" command?",false);
  return in;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrScintSD::FireAPDcell(optHitDetectorMapType* optHitDetectorMap, G4int APDcellID, G4double time, G4int nTruePhe) {
  if ( (musrRootOutput::store_phot_time) && (nTruePhe>0) ) myRootOutput->SetPhotDetTime(time);
  if (!APDcellsEffectRequested) {
    if (verboseLevel>1) G4cout<<"VERBOSE 2 : storing photon hit (non APD)" <<"\n";
    APDidClass tmpAPDid(0,nTruePhe);
    //    optHitDetectorMap->insert(std::pair<G4double,G4int>(time,0));
    optHitDetectorMap->insert(std::pair<G4double,APDidClass>(time,tmpAPDid));
  }
  else {
    G4bool APDcellAlreadyFired = false;
    if (verboseLevel>1) G4cout<<"VERBOSE 2 : searching for the APD element" <<"\n";
    for (optHitDetectorMapType::iterator it2 = optHitDetectorMap->begin(); it2 != optHitDetectorMap->end(); it2++ ) {
      if ((it2->second).GetAPDcellID()==APDcellID) {  // this cell already fired before ==> check times
        if (verboseLevel>1) G4cout<<"VERBOSE 2 : this cell already fired before" <<"\n";
	APDcellAlreadyFired = true;
	if (time<(it2->first)) {       // the new cell fired before the already saved cell ==> replace it
	  if (verboseLevel>1) G4cout<<"VERBOSE 2 : and the new photon is earlier, replacing" <<"\n";
	  G4int tmpAPDcellNphot = (it2->second).GetAPDcellNphot() + nTruePhe;
	  APDidClass tmpAPDid(APDcellID,tmpAPDcellNphot);
	  optHitDetectorMap->erase(it2);
	  //	  optHitDetectorMap->insert(std::pair<G4double,G4int>(time,APDcellID));
	  optHitDetectorMap->insert(std::pair<G4double,APDidClass>(time,tmpAPDid));
	}
	break;
      }
    }
    if (!APDcellAlreadyFired) {
      //      optHitDetectorMap->insert(std::pair<G4double,G4int>(time,APDcellID));
      if (verboseLevel>1) G4cout<<"VERBOSE 2 : this cell hasn't fired, storing" <<"\n";
      APDidClass tmpAPDid(APDcellID,nTruePhe);
      optHitDetectorMap->insert(std::pair<G4double,APDidClass>(time,tmpAPDid));
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrScintSD::FindCFDtime(G4double& OPSA_CFD_time, G4double& OPSA_CFD_ampl, G4double timeOfFirstPhoton) {
  OPSA_CFD->Reset("M");
  for (Int_t iBin=1; iBin<=OPSAhistoNbin; iBin++) { // loop over all bins of electronic signal histogram
    Double_t signalValue = OPSAshape->GetBinContent(iBin);
    Double_t xValue      = (iBin-0.5)*OPSAhistoBinWidth;
    OPSA_CFD->Fill(xValue,signalValue*OPSA_CFD_a1);
    Double_t xValueShifted = xValue + OPSA_CFD_delay;
    if (xValueShifted>OPSAhistoMax) continue;   //return if out of range;
    OPSA_CFD->Fill(xValueShifted,signalValue);
  } 
	    
  // Find the timeCFD from CFD signal
  Double_t oldYvalue;
  Double_t yValue=-1000;
  OPSA_CFD_ampl = OPSA_CFD->GetMaximum();
  int binmax = OPSA_CFD->GetMaximumBin();
  int binmim = OPSA_CFD->GetMinimumBin();
  for (Int_t iBin=binmim; iBin<=binmax; iBin++) { // loop over bins between min and max of CFD histogram
    Double_t xValue = OPSA_CFD->GetXaxis()->GetBinCenter(iBin);
    oldYvalue       = yValue;
    yValue          = OPSA_CFD->GetBinContent(iBin);
    if (yValue >= 0) {       // signal just crossed the y=0 axis
      OPSA_CFD_time = xValue - yValue/(yValue-oldYvalue)*OPSAhistoBinWidth;   //linear interpolation
      OPSA_CFD_time += timeOfFirstPhoton - OPSA_CFD_delay + OPSA_CFD_timeShiftOffset;       // correct for 
      break;
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

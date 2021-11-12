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

#include "globals.hh"
#include "G4ios.hh"
#include "musrPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"

#include "G4MuonDecayChannel.hh"
#include "G4DecayTable.hh"
//cks Added to have Geant default muon decay with spin
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"
#include "G4RadioactiveDecay.hh"
#include "G4IonConstructor.hh"
//TS Classes which account for Muonium as "particle" and its spin
#include "musrMuonium.hh"
#include "MuDecayChannel.hh"
#include "MuDecayChannelWithSpin.hh"
//
#include "musrParameters.hh"
#include "musrErrorMessage.hh"
// cks 2009-06-08   G4StepLimiter and/or G4UserSpecialCuts are needed to activate the "G4UserLimits"
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"
#include "G4OpticalPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrPhysicsList::musrPhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 0.1*CLHEP::mm;
  cutForGamma = 0.1*CLHEP::mm;
  cutForElectron = 0.1*CLHEP::mm;
  cutForMuon = 0.01*CLHEP::mm;
  SetVerboseLevel(0);

  if (musrParameters::boolG4OpticalPhotons) {
    

    //    // * Optical Physics
    //    G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
    //    RegisterPhysics( opticalPhysics );
    //    
    //    // adjust some parameters for the optical physics
    //    opticalPhysics->SetWLSTimeProfile("delta");
    //    
    //    opticalPhysics->SetScintillationYieldFactor(1.0);
    //    opticalPhysics->SetScintillationExcitationRatio(0.0);
    //    
    //    opticalPhysics->SetMaxNumPhotonsPerStep(100);
    //    opticalPhysics->SetMaxBetaChangePerStep(10.0);
    //    
    //    opticalPhysics->SetTrackSecondariesFirst(true);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrPhysicsList::~musrPhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrPhysicsList::ConstructParticle()
{

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();

  if (musrParameters::boolG4OpticalPhotons) {
    // optical photon
    G4OpticalPhoton::OpticalPhotonDefinition();
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrPhysicsList::ConstructLeptons()
{
  // leptons
  //  e+/-
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  // mu+/-
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  //cks
  //  G4DecayTable* MuonPlusDecayTable = new G4DecayTable();
  //  MuonPlusDecayTable -> Insert(new musrMuonDecayChannel("mu+",1.00));
  //  G4MuonPlus::MuonPlusDefinition() -> SetDecayTable(MuonPlusDecayTable);
  //csk
  //
  // Muonium - TS
  musrMuonium::MuoniumDefinition();
  //
  // nu_e
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  // nu_mu
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  //cks:  Trial to use Geant4 muon decay with spin
  G4DecayTable* MuonPlusDecayTable = new G4DecayTable();
  //  MuonPlusDecayTable -> Insert(new G4MuonDecayChannelWithSpin("mu+",1.00));
  MuonPlusDecayTable -> Insert(new G4MuonDecayChannelWithSpin("mu+",0.986));
  MuonPlusDecayTable -> Insert(new G4MuonRadiativeDecayChannelWithSpin("mu+",0.014));
  G4MuonPlus::MuonPlusDefinition() -> SetDecayTable(MuonPlusDecayTable);
  //
  G4DecayTable* MuonMinusDecayTable = new G4DecayTable();
  //  MuonMinusDecayTable -> Insert(new G4MuonDecayChannelWithSpin("mu-",1.00));
  MuonMinusDecayTable -> Insert(new G4MuonDecayChannelWithSpin("mu-",0.986));
  MuonMinusDecayTable -> Insert(new G4MuonRadiativeDecayChannelWithSpin("mu-",0.014));
  G4MuonMinus::MuonMinusDefinition() -> SetDecayTable(MuonMinusDecayTable);
  //csk
  //
  //TS: Using the muonium decay with and without spin
  G4DecayTable* MuoniumDecayTable = new G4DecayTable();
  MuoniumDecayTable -> Insert(new MuDecayChannel("Mu",0.50));
  MuoniumDecayTable -> Insert(new MuDecayChannelWithSpin("Mu",0.5));
  musrMuonium::MuoniumDefinition() -> SetDecayTable(MuoniumDecayTable);
  //MuoniumDecayTable ->DumpInfo(); // Info on muonium decay channels
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrPhysicsList::ConstructMesons()
{
  //  mesons
  //    light mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrPhysicsList::ConstructBaryons()
{
  //  baryons
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();

  // ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"

// #include "G4MultipleScattering.hh"   // obsolete in Geant4.9.4

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"

#include "G4UserSpecialCuts.hh"

//#include "musrAtRestSpinRotation.hh"

//cks 14.12.2011  // For low energy physics processes:
//cks 14.12.2011  #include "G4LowEnergyCompton.hh"
//cks 14.12.2011  //#include "G4LowEnergyPolarizedCompton.hh"
//cks 14.12.2011  #include "G4LowEnergyGammaConversion.hh"
//cks 14.12.2011  #include "G4LowEnergyPhotoElectric.hh"
//cks 14.12.2011  #include "G4LowEnergyRayleigh.hh"
//cks 14.12.2011  #include "G4LowEnergyBremsstrahlung.hh"
//cks 14.12.2011  #include "G4LowEnergyIonisation.hh"
//cks 14.12.2011  #include "G4hLowEnergyIonisation.hh"

//cks 14.12.2011  // For Penelope processes:
//cks 14.12.2011  #include "G4PenelopeCompton.hh"
//cks 14.12.2011  #include "G4PenelopeGammaConversion.hh"
//cks 14.12.2011  #include "G4PenelopePhotoElectric.hh"
//cks 14.12.2011  #include "G4PenelopeRayleigh.hh"
//cks 14.12.2011  #include "G4PenelopeIonisation.hh"
//cks 14.12.2011  #include "G4PenelopeBremsstrahlung.hh"
//cks 14.12.2011  #include "G4PenelopeAnnihilation.hh"

// For Coulomb scattering instead of multiple scattering
#include "G4CoulombScattering.hh"
//cks 14.12.2011  #include "G4CoulombScatteringModel.hh"

// For Muonium formation in the Carbon foil
#include "musrMuFormation.hh" // includes the yield function Y = Y(E).

// For muon energy loss in thin Carbon foil
#include "musrMuEnergyLossLandau.hh"

// For a simple Muonium "scattering" when Mu hits solid materials
#include "musrMuScatter.hh"

// For optical photons
#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"

// For models
#include "G4ProcessTable.hh"
#include "G4WentzelVIModel.hh"
//#include "G4UrbanMscModel90.hh"
//cks 14.12.2011  #include "G4UrbanMscModel92.hh"
//#include "G4UrbanMscModel93.hh"
#include "G4UrbanMscModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrPhysicsList::ConstructEM()
{
  // cks  2008.08.22. - Adding the possibility to define the processes from the steering file:
  char charSteeringFileName[1000]; strcpy(charSteeringFileName,(musrParameters::mySteeringFileName).c_str());
  FILE *fSteeringFile=fopen(charSteeringFileName,"r");
  if (fSteeringFile==NULL) {
    sprintf(eMessage,"musrPhysicsList::ConstructEM():  Failed to open macro file \"%s\" .",charSteeringFileName);
                     musrErrorMessage::GetInstance()->musrError(FATAL,eMessage,false);
  }
  else {G4cout<<"musrPhysicsList: The Physics processes are being defined:"<<G4endl;}

  char  line[501];
  while (!feof(fSteeringFile)) {
    fgets(line,500,fSteeringFile);
    if ((line[0]!='#')&&(line[0]!='\n')&&(line[0]!='\r')) {
      char tmpString0[100]="Unset", tmpString1[100]="Unset",tmpString2[100]="Unset";
      sscanf(&line[0],"%s %s %s",tmpString0,tmpString1,tmpString2);
      if ( (strcmp(tmpString0,"/musr/ignore")!=0)&&(strcmp(tmpString0,"/musr/command")!=0) ) continue;
      if (strcmp(tmpString1,"process")!=0)  continue; 
      
      float tmpCutValue;
      if ((strcmp(tmpString2,"addProcess")==0)||(strcmp(tmpString2,"addDiscreteProcess")==0)||
	  (strcmp(tmpString2,"addModel")==0)||(strcmp(tmpString2,"SetLossFluctuations_OFF")==0)) {
	char charParticleName[100], charProcessName[100];
	sscanf(&line[0],"%*s %*s %s %s %s",tmpString2,charParticleName,charProcessName);
	G4cout<<"musrPhysicsList: Defining process "<<charProcessName<<" for "<<charParticleName<<G4endl;
	G4String stringProcessName  = charProcessName;
	G4String stringParticleName = charParticleName;
	G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable() -> FindParticle(stringParticleName);
	//      G4cout<<"particleDefinition of "<<stringParticleName<<" = "<<particleDefinition<<G4endl;
	if (particleDefinition==NULL) {
	  if ((stringParticleName=="opticalphoton")&& (!(musrParameters::boolG4OpticalPhotons))) {
	    musrErrorMessage::GetInstance()->musrError(INFO,
						       "musrPhysicsList: Ignoring optical photon process definition because G4OpticalPhotons in *.mac file is set to false",false);
	    continue;
	  }
	  sprintf(eMessage,"musrPhysicsList:  Partile \"%s\" not found in G4ParticleTable when trying to find or assign process  \"%s\".",
		  charParticleName,charProcessName);
	  musrErrorMessage::GetInstance()->musrError(FATAL,eMessage,false);
	}
	G4ProcessManager* pManager = particleDefinition->GetProcessManager();

	if (strcmp(tmpString2,"addDiscreteProcess")==0) {
	  if      (stringProcessName=="G4PhotoElectricEffect")      pManager->AddDiscreteProcess(new G4PhotoElectricEffect);
	  else if (stringProcessName=="G4ComptonScattering")        pManager->AddDiscreteProcess(new G4ComptonScattering);
	  else if (stringProcessName=="G4GammaConversion")          pManager->AddDiscreteProcess(new G4GammaConversion);
	  else if (stringProcessName=="G4RayleighScattering")       pManager->AddDiscreteProcess(new G4RayleighScattering);
//cks 14.12.2011  	  else if (stringProcessName=="G4PenelopePhotoElectric")    pManager->AddDiscreteProcess(new G4PenelopePhotoElectric);
//cks 14.12.2011  	  else if (stringProcessName=="G4PenelopeCompton")          pManager->AddDiscreteProcess(new G4PenelopeCompton);
//cks 14.12.2011  	  else if (stringProcessName=="G4PenelopeGammaConversion")  pManager->AddDiscreteProcess(new G4PenelopeGammaConversion);
//cks 14.12.2011  	  else if (stringProcessName=="G4PenelopeRayleigh")         pManager->AddDiscreteProcess(new G4PenelopeRayleigh);
//cks 14.12.2011  	  else if (stringProcessName=="G4LowEnergyPhotoElectric")   pManager->AddDiscreteProcess(new G4LowEnergyPhotoElectric);
//cks 14.12.2011  	  else if (stringProcessName=="G4LowEnergyCompton")         pManager->AddDiscreteProcess(new G4LowEnergyCompton);
//cks 14.12.2011  	  else if (stringProcessName=="G4LowEnergyGammaConversion") pManager->AddDiscreteProcess(new G4LowEnergyGammaConversion);
//cks 14.12.2011  	  else if (stringProcessName=="G4LowEnergyRayleigh")        pManager->AddDiscreteProcess(new G4LowEnergyRayleigh);

	  else if (stringProcessName=="G4CoulombScattering")        pManager->AddDiscreteProcess(new G4CoulombScattering);
	  
	  else if (stringProcessName=="G4OpAbsorption")             pManager->AddDiscreteProcess(new G4OpAbsorption);
	  else if (stringProcessName=="G4OpRayleigh")               pManager->AddDiscreteProcess(new G4OpRayleigh);
	  else if (stringProcessName=="G4OpBoundaryProcess")        pManager->AddDiscreteProcess(new G4OpBoundaryProcess);
	  else if (stringProcessName=="G4OpWLS")                    pManager->AddDiscreteProcess(new G4OpWLS);
	  else if (stringProcessName=="G4Cerenkov")              {
					G4Cerenkov* myCerenkov = new G4Cerenkov();
					  myCerenkov->SetMaxNumPhotonsPerStep(20); // settings from GEANT example N06
					  myCerenkov->SetMaxBetaChangePerStep(10.0);
					  myCerenkov->SetTrackSecondariesFirst(false); // as for Scintillation, do main high energy particles first
					pManager->AddDiscreteProcess(myCerenkov); // added JSL
	  }
	  else {
	    sprintf(eMessage,"musrPhysicsList:  Process \"%s\" is not implemented in musrPhysicsList.cc for addDiscreteProcess.  It can be easily added.",
		    charProcessName);
	    musrErrorMessage::GetInstance()->musrError(FATAL,eMessage,false);
	  }
	}
	else if (strcmp(tmpString2,"addProcess")==0) {
	  G4int nr1, nr2, nr3;
	  char charRegion1[100]="", charRegion2[100]="", charRegion3[100]="", charControlString[10]="";
	  sscanf(&line[0],"%*s %*s %*s %*s %*s %d %d %d %s %s %s %s",&nr1,&nr2,&nr3,charRegion1,charRegion2,charRegion3,charControlString);
	  //	  if      (stringProcessName=="G4MultipleScattering")         pManager->AddProcess(new G4MultipleScattering,nr1,nr2,nr3); 
	                                      // obsolete in Geant4.9.4
	  if (stringProcessName=="G4eMultipleScattering")        pManager->AddProcess(new G4eMultipleScattering,nr1,nr2,nr3);
	  else if (stringProcessName=="G4eIonisation")                pManager->AddProcess(new G4eIonisation,nr1,nr2,nr3);
	  else if (stringProcessName=="G4eBremsstrahlung")            pManager->AddProcess(new G4eBremsstrahlung,nr1,nr2,nr3);
	  else if (stringProcessName=="G4eplusAnnihilation")          pManager->AddProcess(new G4eplusAnnihilation,nr1,nr2,nr3);
//cks 14.12.2011  	  else if (stringProcessName=="G4PenelopeIonisation")         pManager->AddProcess(new G4PenelopeIonisation,nr1,nr2,nr3);
//cks 14.12.2011  	  else if (stringProcessName=="G4PenelopeBremsstrahlung")     pManager->AddProcess(new G4PenelopeBremsstrahlung,nr1,nr2,nr3);
//cks 14.12.2011  	  else if (stringProcessName=="G4PenelopeAnnihilation")       pManager->AddProcess(new G4PenelopeAnnihilation,nr1,nr2,nr3);
  	//  else if (stringProcessName=="G4LowEnergyIonisation")        pManager->AddProcess(new G4LowEnergyIonisation,nr1,nr2,nr3);
  	//  else if (stringProcessName=="G4LowEnergyBremsstrahlung")    pManager->AddProcess(new G4LowEnergyBremsstrahlung,nr1,nr2,nr3);
          else if (stringProcessName=="G4MuMultipleScattering")       pManager->AddProcess(new G4MuMultipleScattering,nr1,nr2,nr3);
	  else if (stringProcessName=="G4MuIonisation")               pManager->AddProcess(new G4MuIonisation,nr1,nr2,nr3);
	  else if (stringProcessName=="G4MuBremsstrahlung")           pManager->AddProcess(new G4MuBremsstrahlung,nr1,nr2,nr3);
	  else if (stringProcessName=="G4MuPairProduction")           pManager->AddProcess(new G4MuPairProduction,nr1,nr2,nr3);
          // cks 2009-06-08   G4StepLimiter and/or G4UserSpecialCuts are needed to activate the "G4UserLimits"
	  else if (stringProcessName=="G4StepLimiter")                pManager->AddProcess(new G4StepLimiter,nr1,nr2,nr3);
	  else if (stringProcessName=="G4UserSpecialCuts")            pManager->AddProcess(new G4UserSpecialCuts,nr1,nr2,nr3);
//	  else if (stringProcessName=="G4DecayWithSpin")      pManager->AddProcess(new G4DecayWithSpin,nr1,nr2,nr3);
          //  else if (stringProcessName=="G4hIonisation")        pManager->AddProcess(new G4hIonisation,nr1,nr2,nr3);
	//  else if (stringProcessName=="G4hLowEnergyIonisation")      pManager->AddProcess(new G4hLowEnergyIonisation,nr1,nr2,nr3);
	  else if (stringProcessName=="musrMuFormation")              pManager->AddProcess(new musrMuFormation,nr1,nr2,nr3);
	  else if (stringProcessName=="musrMuEnergyLossLandau")       pManager->AddProcess(new musrMuEnergyLossLandau,nr1,nr2,nr3);
	  else if (stringProcessName=="G4Scintillation")              {
        	      G4Scintillation* myScintillation = new G4Scintillation();
		      //		      if (musrParameters::finiteRiseTimeInScintillator) myScintillation->SetFiniteRiseTime(true);
	              pManager->AddProcess(myScintillation,nr1,nr2,nr3);
	  }
	  else if (stringProcessName=="G4Cerenkov")              {
					G4Cerenkov* myCerenkov = new G4Cerenkov();
					  myCerenkov->SetMaxNumPhotonsPerStep(20); // settings from GEANT example N06
					  myCerenkov->SetMaxBetaChangePerStep(10.0);
					  myCerenkov->SetTrackSecondariesFirst(false); // as for Scintillation, do main high energy particles first
					pManager->AddProcess(myCerenkov,nr1,nr2,nr3); // added JSL
					pManager->SetProcessOrdering(myCerenkov,idxPostStep); // from Example N06
	  }
	  // cks:  musrMuScatter could be uncommented here, but testing is needed, because Toni has some strange comments
	  //       in his original "musrPhysicsList.cc about implementing musrMuScatter.
	  //	  else if (stringProcessName=="musrMuScatter")       pManager->AddProcess(new musrMuScatter,nr1,nr2,nr3);
	  //
	  // The G4MultipleScattering is obsolete in Geant4.9.4
	  //	  else if (stringProcessName=="MultipleAndCoulombScattering") {
	  //	    G4MultipleScattering* multScat = new G4MultipleScattering();
	  //	    //	    G4CoulombScattering*  coulScat = new G4CoulombScattering();
	  //	    G4CoulombScatteringModel* coulScatModel = new G4CoulombScatteringModel();
	  //	    if (strcmp(charRegion1,"")!=0) {
	  //	      G4Region* regionForCoulomb = FindG4Region(charRegion1,line);
	  //	      G4cout<<"  Adding Coulomb scattering model to multiple scattering model for region "<<charRegion1<<G4endl;
	  //	      multScat->AddEmModel(0,coulScatModel,regionForCoulomb);
	  //	      //	      multScat->AddEmModel(0,multScat,regionForCoulomb);
	  //	    }
	  //	    if (strcmp(charRegion2,"")!=0) {
	  //	      G4Region* regionForCoulomb = FindG4Region(charRegion2,line);
	  //	      G4cout<<"  Adding Coulomb scattering model to multiple scattering model for region "<<charRegion2<<G4endl;
	  //	      multScat->AddEmModel(0,coulScatModel,regionForCoulomb);
	  //	    }
	  //	    if (strcmp(charRegion3,"")!=0) {
	  //	      G4Region* regionForCoulomb = FindG4Region(charRegion3,line);
	  //	      G4cout<<"  Adding Coulomb scattering model to multiple scattering model for region "<<charRegion3<<G4endl;
	  //	      multScat->AddEmModel(0,coulScatModel,regionForCoulomb);
	  //	    }
	  //	    if (strcmp(charControlString,"")!=0) {
	  //	      G4cout<<"More than 3 regions requested for Coulomb Scattering, but presently only up to 3 such regions are supported."<<G4endl;
	  //	      G4cout<<"Please extend the number of supported regions in musrPhysicsList.cc to higher number."<<G4endl;
	  //	      G4cout<<"The extention of the code to larger number of regions is not very difficult."<<G4endl;
	  //	      G4cout<<"   S T O P      F O R C E D"<<G4endl;
	  //	      exit(1);
	  //	    }
	  //	    pManager->AddProcess(multScat,nr1,nr2,nr3);
	  //	  }
	  else {
	    sprintf(eMessage,"musrPhysicsList:  Process \"%s\" is not implemented in musrPhysicsList.cc for addProcess.  It can be easily added.",
		    charProcessName);
	    musrErrorMessage::GetInstance()->musrError(FATAL,eMessage,false);
	  }
	}
	else if (strcmp(tmpString2,"addModel")==0) {
	  G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();
	  //	  G4VProcess* myProc =  processTable->FindProcess(charProcessName,particleDefinition);
	  // pManager
	  
	  G4int modelPriority;
	  char charModelName[100]="";
	  sscanf(&line[0],"%*s %*s %*s %*s %*s %s %d",charModelName,&modelPriority);
	  G4String stringModelName = charModelName;
	  G4cout<<" Adding model "<<charModelName<<" with priority "<<modelPriority<<" to process "<<charProcessName<<"."<<G4endl;
	  //	  if (myProc==NULL) {
	  //	    sprintf(eMessage,"musrPhysicsList:  Trying to add model \"%s\" , but the process \"%s\" not found in G4ProcessTable.",
	  //		    charModelName,charProcessName);
	  //	    musrErrorMessage::GetInstance()->musrError(FATAL,eMessage,false);
	  //	  }
	  if      ((stringModelName=="G4WentzelVIModel")&&(stringProcessName=="G4MuMultipleScattering")) {
	    G4MuMultipleScattering* mmm = (G4MuMultipleScattering*) processTable->FindProcess("muMsc",particleDefinition);
	    mmm->AddEmModel(modelPriority, new G4WentzelVIModel());
	  }
	  else if ((stringModelName=="G4UrbanMscModel")&&(stringProcessName=="G4MuMultipleScattering")) {
	    G4MuMultipleScattering* mmm = (G4MuMultipleScattering*) processTable->FindProcess("muMsc",particleDefinition);
	    mmm->AddEmModel(modelPriority, new G4UrbanMscModel());
	  }
//cks 14.12.2011  	  else if ((stringModelName=="G4UrbanMscModel92")&&(stringProcessName=="G4MuMultipleScattering")) {
//cks 14.12.2011  	    G4MuMultipleScattering* mmm = (G4MuMultipleScattering*) processTable->FindProcess("muMsc",particleDefinition);
//cks 14.12.2011  	    mmm->AddEmModel(modelPriority, new G4UrbanMscModel92());
//cks 14.12.2011  	  }
//	  else if ((stringModelName=="G4UrbanMscModel93")&&(stringProcessName=="G4MuMultipleScattering")) {
//	    G4MuMultipleScattering* mmm = (G4MuMultipleScattering*) processTable->FindProcess("muMsc",particleDefinition);
//	    mmm->AddEmModel(modelPriority, new G4UrbanMscModel93());
//	  }
	  else {
	    sprintf(eMessage,"musrPhysicsList:  Model \"%s\" is not implemented for \"%s\" in musrPhysicsList.cc for addModel.  It can be easily added.",
		    charModelName,charProcessName);
	    musrErrorMessage::GetInstance()->musrError(FATAL,eMessage,false);
	  }

	}

	else if (strcmp(tmpString2,"SetLossFluctuations_OFF")==0) {
	  G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();
	  if (strcmp(charProcessName,"G4eIonisation")==0) {
	    G4cout<<"musrPhysicsList.cc:  Switching off loss fluctuations for "<<charParticleName<<" in "<<charProcessName<<G4endl;
	    G4eIonisation* mmm = (G4eIonisation*) processTable->FindProcess("eIoni",particleDefinition);
	    mmm->SetLossFluctuations(false);
	  }
	  else if (strcmp(charProcessName,"G4MuIonisation")==0) {
	    G4cout<<"musrPhysicsList.cc:  Switching off loss fluctuations for "<<charParticleName<<" in "<<charProcessName<<G4endl;
	    G4MuIonisation* mmm = (G4MuIonisation*) processTable->FindProcess("muIoni",particleDefinition);
	    mmm->SetLossFluctuations(false);
	  }
	  else {
	    G4cout<<"musrPhysicsList.cc:  Switching off loss fluctuations for "
		  <<charParticleName<<" in "<<charProcessName<<" NOT SUPPORTED YET"<<G4endl;
	  }
	}


      }


      else if (strcmp(tmpString2,"cutForGamma")==0) {
	sscanf(&line[0],"%*s %*s %*s %g",&tmpCutValue);
	cutForGamma = tmpCutValue;
      }

      else if (strcmp(tmpString2,"cutForElectron")==0) {
	sscanf(&line[0],"%*s %*s %*s %g",&tmpCutValue);
	cutForElectron = tmpCutValue;
      }

      //cks UNFORTUNATELY cutForMuon does not work!!!
      //      else if (strcmp(tmpString2,"cutForMuon")==0) {
      //	sscanf(&line[0],"%*s %*s %*s %g",&tmpCutValue);
      //	cutForMuon = tmpCutValue;
      //      }


      else ReportProblemWithProcessDefinition(line);
    }
  }
  fclose(fSteeringFile);

  G4cout<<"\n\n\n\n"<<G4endl;
  // csk  2008.08.22.

  //del G4String myTypeOfProcesses  = musrParameters::GetInstance()->GetMyTypeOfProcesses();
  //del  G4cout<<"musrPhysicsList::ConstructEM():  myTypeOfProcesses="<<myTypeOfProcesses<<G4endl;

  auto theParticleIterator=GetParticleIterator(); //Geant4.10.3
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if ((particleName == "gamma")||(particleName == "e-")||(particleName == "e+")) {
      // do nothing
    }
      
    else if ((particleName=="mu+")||(particleName=="mu-")) {          //muon  
      G4DecayWithSpin* theDecayProcess = new G4DecayWithSpin();
      //      theDecayProcess->SetVerboseLevel(2);
      pmanager->AddProcess(theDecayProcess);
      pmanager ->SetProcessOrderingToLast(theDecayProcess, idxAtRest);
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
    }

    else if (particleName=="Mu") {
      //  TS:
      // Muonium "scattering"    Kamil:  the following 3 lines could be replaced by reading the musrMuScatter 
      //                                 process through the steering file
      G4VProcess* aMuScatt = new musrMuScatter();
      pmanager->AddProcess(aMuScatt);
      pmanager->SetProcessOrdering(aMuScatt, idxPostStep, 1);
      //
      G4Decay* theDecayProcess = new G4Decay();
      //musrDecayWithSpin* theDecayProcess = new musrDecayWithSpin();
      pmanager->AddProcess(theDecayProcess);
      pmanager->SetProcessOrderingToLast(theDecayProcess, idxAtRest);
      pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
    }
      
    
    else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4hMultipleScattering,-1, 1,1);
      //      if (myTypeOfProcesses=="highenergy") {
      //	pmanager->AddProcess(new G4hIonisation,       -1, 2,2);
      //    }
      //    else {
      pmanager->AddProcess(new G4hIonisation,       -1, 2,2);
      //    }
      ///pmanager->AddProcess(new G4UserSpecialCuts,  -1,-1,3);      
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



#include "G4Decay.hh"
void musrPhysicsList::ConstructGeneral()  {
  if (musrParameters::boolG4GeneralParticleSource) {
    G4RadioactiveDecay*  theRadioactiveDecay = new G4RadioactiveDecay();
    G4GenericIon* ion = G4GenericIon::GenericIon();
    
    auto theParticleIterator=GetParticleIterator(); //Geant4.10.3
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      
      if (particle == ion) {
	pmanager->AddProcess(theRadioactiveDecay, 0, -1, 3);
      }
    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

void musrPhysicsList::SetCuts()
{
  //G4VUserPhysicsList::SetCutsWithDefault method sets 
  //the default cut value for all particle types 
  //
  //cks 6.10.2009  SetCutsWithDefault();

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForElectron, "e+");
  //cks - This command does not work for muons:  SetCutValue(cutForMuon, "mu-");
  //cks - This command does not work for muons:  SetCutValue(cutForMuon, "mu+");

  //  G4cout<<"Kamil:   cutForGamma = "<<cutForGamma/mm<<" mm"<<G4endl;
  //  G4cout<<"Kamil:   cutForElectron = "<<cutForElectron/mm<<" mm"<<G4endl;
  //  G4cout<<"Kamil:   cutForMuons = "<<cutForMuon/mm<<" mm"<<G4endl;

    
  if (verboseLevel>0) DumpCutValuesTable();
  DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrPhysicsList::ReportProblemWithProcessDefinition(char myString[501]) {
  G4cout<<"\nE R R O R  in musrPhysicsList.cc:    "
	<<"Unknown keyword requested in the steering (*.mac) file :"<<G4endl;
  G4cout<<"  "<<myString<<G4endl;
  G4cout<<"S T O P     F O R C E D!"<<G4endl;
  exit(1);
}

G4Region* musrPhysicsList::FindG4Region(G4String regionName, char* lineOfSteeringFile) {
  G4Region* myRegion = G4RegionStore::GetInstance()->GetRegion(regionName,false);
  if( myRegion != NULL )  { // G4Region found
    return myRegion;
  }
  else {  // G4Region not found
    G4cout<<"musrPhysicsList:  G4Region "<<regionName<<" not found."<<G4endl;
    G4cout<<"                  The critical command line of the steering file is:"<<G4endl;
    G4cout<<"  "<<lineOfSteeringFile<<G4endl;
    G4cout<<"      S T O P     F O R C E D"<<G4endl;
    exit(1);
  }
}

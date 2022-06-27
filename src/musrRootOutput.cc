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

#include "musrRootOutput.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "musrErrorMessage.hh"
#include "musrParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrRootOutput::musrRootOutput(std::string name) {
    TTree::SetMaxTreeSize(100000000000LL);      // Set maximum size of the tree file
    // to 100 GB (instead of 1.9 GB).
    run_name = name;
    pointerToRoot=this;
    boolIsAnySpecialSaveVolumeDefined=false;
    nFieldNomVal=0;
    strcpy(rootOutputDirectoryName,"data");

    ProcessIDMapping["DecayWithSpin"]=1;
    ProcessIDMapping["eIoni"]=2;
    ProcessIDMapping["eBrem"]=3;
    ProcessIDMapping["annihil"]=4;
    ProcessIDMapping["LowEnCompton"]=5;
    ProcessIDMapping["LowEnConversion"]=6;
    ProcessIDMapping["LowEnBrem"]=7;
    ProcessIDMapping["LowEnergyIoni"]=8;
    ProcessIDMapping["LowEnPhotoElec"]=9;
    ProcessIDMapping["RadioactiveDecay"]=10;
    ProcessIDMapping["muIoni"]=11;
    ProcessIDMapping["MuFormation"]=12;
    ProcessIDMapping["Decay"]=13;
    ProcessIDMapping["conv"]=14;
    ProcessIDMapping["compt"]=15;
    ProcessIDMapping["phot"]=16;
    ProcessIDMapping["dInelastic"] = 17;
    ProcessIDMapping["electronNuclear"] = 18;
    ProcessIDMapping["GammaToMuPair"] = 19;
    ProcessIDMapping["hBertiniCaptureAtRest"] = 20;
    ProcessIDMapping["hIoni"] = 21;
    ProcessIDMapping["lambdaInelastic"] = 22;
    ProcessIDMapping["muMinusCaptureAtRest"] = 23;
    ProcessIDMapping["muonNuclear"] = 24;
    ProcessIDMapping["nCapture"] = 25;
    ProcessIDMapping["neutronInelastic"] = 26;
    ProcessIDMapping["photonNuclear"] = 27;
    ProcessIDMapping["pi-Inelastic"] = 28;
    ProcessIDMapping["pi+Inelastic"] = 29;
    ProcessIDMapping["protonInelastic"] = 30;
    ProcessIDMapping["tInelastic"] = 31;
    ProcessIDMapping["hadElastic"] = 32;
    ProcessIDMapping["initialParticle"]=100;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrRootOutput::~musrRootOutput() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrRootOutput* musrRootOutput::pointerToRoot=0;
musrRootOutput* musrRootOutput::GetRootInstance() {
    return pointerToRoot;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool musrRootOutput::store_runID = true;
G4bool musrRootOutput::store_eventID = true;
G4bool musrRootOutput::store_weight = true;
G4bool musrRootOutput::store_timeToNextEvent = true;
G4bool musrRootOutput::store_BFieldAtDecay = true;
G4bool musrRootOutput::store_muIniTime = true;
G4bool musrRootOutput::store_muIniPosX = true;
G4bool musrRootOutput::store_muIniPosY = true;
G4bool musrRootOutput::store_muIniPosZ = true;
G4bool musrRootOutput::store_muIniMomX = true;
G4bool musrRootOutput::store_muIniMomY = true;
G4bool musrRootOutput::store_muIniMomZ = true;
G4bool musrRootOutput::store_muIniPolX = true;
G4bool musrRootOutput::store_muIniPolY = true;
G4bool musrRootOutput::store_muIniPolZ = true;
G4bool musrRootOutput::store_muDecayDetID= true;
G4bool musrRootOutput::store_muDecayPosX = true;
G4bool musrRootOutput::store_muDecayPosY = true;
G4bool musrRootOutput::store_muDecayPosZ = true;
G4bool musrRootOutput::store_muDecayTime = true;
G4bool musrRootOutput::store_muDecayPolX = true;
G4bool musrRootOutput::store_muDecayPolY = true;
G4bool musrRootOutput::store_muDecayPolZ = true;
G4bool musrRootOutput::store_muTargetTime = false;
G4bool musrRootOutput::store_muTargetPolX = false;
G4bool musrRootOutput::store_muTargetPolY = false;
G4bool musrRootOutput::store_muTargetPolZ = false;
G4bool musrRootOutput::store_muTargetMomX = false;
G4bool musrRootOutput::store_muTargetMomY = false;
G4bool musrRootOutput::store_muTargetMomZ = false;
G4bool musrRootOutput::store_muM0Time = false;
G4bool musrRootOutput::store_muM0PolX = false;
G4bool musrRootOutput::store_muM0PolY = false;
G4bool musrRootOutput::store_muM0PolZ = false;
G4bool musrRootOutput::store_muM1Time = false;
G4bool musrRootOutput::store_muM1PolX = false;
G4bool musrRootOutput::store_muM1PolY = false;
G4bool musrRootOutput::store_muM1PolZ = false;
G4bool musrRootOutput::store_muM2Time = false;
G4bool musrRootOutput::store_muM2PolX = false;
G4bool musrRootOutput::store_muM2PolY = false;
G4bool musrRootOutput::store_muM2PolZ = false;
G4bool musrRootOutput::store_posIniMomX = true;
G4bool musrRootOutput::store_posIniMomY = true;
G4bool musrRootOutput::store_posIniMomZ = true;
G4bool musrRootOutput::store_nOptPhot   = true;
G4bool musrRootOutput::store_nOptPhotDet= true;
G4bool musrRootOutput::store_phot_time  = false;
G4bool musrRootOutput::store_det_ID = true;
G4bool musrRootOutput::store_det_edep = true;
G4bool musrRootOutput::store_det_edep_el = true;
G4bool musrRootOutput::store_det_edep_pos = true;
G4bool musrRootOutput::store_det_edep_gam = true;
G4bool musrRootOutput::store_det_edep_mup = true;
G4bool musrRootOutput::store_det_kine_mup = true;
G4bool musrRootOutput::store_det_edep_mun = true;
G4bool musrRootOutput::store_det_kine_mun = true;
G4bool musrRootOutput::store_det_nsteps = true;
G4bool musrRootOutput::store_det_length = true;
G4bool musrRootOutput::store_det_start = true;
G4bool musrRootOutput::store_det_end = true;
G4bool musrRootOutput::store_det_x = true;
G4bool musrRootOutput::store_det_y = true;
G4bool musrRootOutput::store_det_z = true;
G4bool musrRootOutput::store_det_x_mup = true;
G4bool musrRootOutput::store_det_y_mup = true;
G4bool musrRootOutput::store_det_z_mup = true;
G4bool musrRootOutput::store_det_x_mun = true;
G4bool musrRootOutput::store_det_y_mun = true;
G4bool musrRootOutput::store_det_z_mun = true;
G4bool musrRootOutput::store_det_kine = true;
G4bool musrRootOutput::store_det_VrtxKine = true;
G4bool musrRootOutput::store_det_VrtxX = true;
G4bool musrRootOutput::store_det_VrtxY = true;
G4bool musrRootOutput::store_det_VrtxZ = true;
G4bool musrRootOutput::store_det_VrtxVolID = true;
G4bool musrRootOutput::store_det_VrtxProcID = true;
G4bool musrRootOutput::store_det_VrtxTrackID = true;
G4bool musrRootOutput::store_det_VrtxParentTrackID = true;
G4bool musrRootOutput::store_det_VrtxParticleID = true;
G4bool musrRootOutput::store_det_VvvKine = true;
G4bool musrRootOutput::store_det_VvvX = true;
G4bool musrRootOutput::store_det_VvvY = true;
G4bool musrRootOutput::store_det_VvvZ = true;
G4bool musrRootOutput::store_det_VvvVolID = true;
G4bool musrRootOutput::store_det_VvvProcID = true;
G4bool musrRootOutput::store_det_VvvTrackID = true;
G4bool musrRootOutput::store_det_VvvParticleID = true;
G4bool musrRootOutput::store_fieldNomVal = true;
G4bool musrRootOutput::store_fieldIntegralBx = false;
G4bool musrRootOutput::store_fieldIntegralBy = false;
G4bool musrRootOutput::store_fieldIntegralBz = false;
G4bool musrRootOutput::store_fieldIntegralBz1 = false;
G4bool musrRootOutput::store_fieldIntegralBz2 = false;
G4bool musrRootOutput::store_fieldIntegralBz3 = false;
G4bool musrRootOutput::store_odet_ID = true;
G4bool musrRootOutput::store_odet_nPhot = true;
G4bool musrRootOutput::store_odet_nPhotPrim = true;
G4bool musrRootOutput::store_odet_timeFirst = true;
G4bool musrRootOutput::store_odet_timeSecond = true;
G4bool musrRootOutput::store_odet_timeThird = true;
G4bool musrRootOutput::store_odet_timeA = true;
G4bool musrRootOutput::store_odet_timeB = true;
G4bool musrRootOutput::store_odet_timeC = true;
G4bool musrRootOutput::store_odet_timeD = true;
G4bool musrRootOutput::store_odet_timeMean = true;
G4bool musrRootOutput::store_odet_timeLast = true;
G4bool musrRootOutput::store_odet_timeCFD = true;
G4bool musrRootOutput::store_odet_amplCFD = true;
G4bool musrRootOutput::store_odet_timeCFDarray = false;
G4bool musrRootOutput::store_odet_timeC1 = false;

G4int musrRootOutput::oldEventNumberInG4EqEMFieldWithSpinFunction=-1;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrRootOutput::BeginOfRunAction() {
    G4cout << "musrRootOutput::BeginOfRunAction()  Defining the Root tree and branches:"<<G4endl;
    G4int tmpRunNr=(G4RunManager::GetRunManager())->GetCurrentRun()->GetRunID();
    char RootOutputFileName[200];
    //  sprintf(RootOutputFileName, "data/musr_%i.root", tmpRunNr);
    if(run_name!="") {
        sprintf(RootOutputFileName, "%s/musrSim_%i_%s.root",rootOutputDirectoryName,tmpRunNr,run_name.c_str());
    }
    else {
        sprintf(RootOutputFileName, "%s/musrSim_%i.root",rootOutputDirectoryName,tmpRunNr);
    }

    rootFile=new TFile(RootOutputFileName,"recreate");
    if (rootFile->IsZombie()) {
        char message[200];
        sprintf(message,"musrRootOutput::BeginOfRunAction() Root output file %s can not be created",RootOutputFileName);
        musrErrorMessage::GetInstance()->musrError(FATAL,message,false);
    }
    rootTree=new TTree("t1","a simple Tree with simple variables");
    if (store_runID)        {rootTree->Branch("runID",&runID,"runID/I");}
    if (store_eventID)      {rootTree->Branch("eventID",&eventID,"eventID/I");}
    if (store_weight)       {rootTree->Branch("weight",&weight,"weight/D");}
    if (store_timeToNextEvent){rootTree->Branch("timeToNextEvent",&timeToNextEvent,"timeToNextEvent/D");}
    if (store_BFieldAtDecay) {rootTree->Branch("BFieldAtDecay",&B_t,"Bx/D:By:Bz:B3:B4:B5");}
    if (store_muIniTime)    {rootTree->Branch("muIniTime",&muIniTime,"muIniTime/D");}
    if (store_muIniPosX)    {rootTree->Branch("muIniPosX",&muIniPosX,"muIniPosX/D");}
    if (store_muIniPosY)    {rootTree->Branch("muIniPosY",&muIniPosY,"muIniPosY/D");}
    if (store_muIniPosZ)    {rootTree->Branch("muIniPosZ",&muIniPosZ,"muIniPosZ/D");}
    if (store_muIniMomX)    {rootTree->Branch("muIniMomX",&muIniMomX,"muIniMomX/D");}
    if (store_muIniMomY)    {rootTree->Branch("muIniMomY",&muIniMomY,"muIniMomY/D");}
    if (store_muIniMomZ)    {rootTree->Branch("muIniMomZ",&muIniMomZ,"muIniMomZ/D");}
    if (store_muIniPolX)    {rootTree->Branch("muIniPolX",&muIniPolX,"muIniPolX/D");}
    if (store_muIniPolY)    {rootTree->Branch("muIniPolY",&muIniPolY,"muIniPolY/D");}
    if (store_muIniPolZ)    {rootTree->Branch("muIniPolZ",&muIniPolZ,"muIniPolZ/D");}
    if (store_muDecayDetID) {rootTree->Branch("muDecayDetID",&muDecayDetID,"muDecayDetID/I");}
    if (store_muDecayPosX)  {rootTree->Branch("muDecayPosX",&muDecayPosX,"muDecayPosX/D");}
    if (store_muDecayPosY)  {rootTree->Branch("muDecayPosY",&muDecayPosY,"muDecayPosY/D");}
    if (store_muDecayPosZ)  {rootTree->Branch("muDecayPosZ",&muDecayPosZ,"muDecayPosZ/D");}
    if (store_muDecayTime)  {rootTree->Branch("muDecayTime",&muDecayTime,"muDecayTime/D");}
    if (store_muDecayPolX)  {rootTree->Branch("muDecayPolX",&muDecayPolX,"muDecayPolX/D");}
    if (store_muDecayPolY)  {rootTree->Branch("muDecayPolY",&muDecayPolY,"muDecayPolY/D");}
    if (store_muDecayPolZ)  {rootTree->Branch("muDecayPolZ",&muDecayPolZ,"muDecayPolZ/D");}
    if (store_muTargetTime) {rootTree->Branch("muTargetTime",&muTargetTime,"muTargetTime/D");}
    if (store_muTargetPolX) {rootTree->Branch("muTargetPolX",&muTargetPolX,"muTargetPolX/D");}
    if (store_muTargetPolY) {rootTree->Branch("muTargetPolY",&muTargetPolY,"muTargetPolY/D");}
    if (store_muTargetPolZ) {rootTree->Branch("muTargetPolZ",&muTargetPolZ,"muTargetPolZ/D");}
    if (store_muTargetMomX) {rootTree->Branch("muTargetMomX",&muTargetMomX,"muTargetMomX/D");}
    if (store_muTargetMomY) {rootTree->Branch("muTargetMomY",&muTargetMomY,"muTargetMomY/D");}
    if (store_muTargetMomZ) {rootTree->Branch("muTargetMomZ",&muTargetMomZ,"muTargetMomZ/D");}
    if (store_muM0Time)     {rootTree->Branch("muM0Time",&muM0Time,"muM0Time/D");}
    if (store_muM0PolX)     {rootTree->Branch("muM0PolX",&muM0PolX,"muM0PolX/D");}
    if (store_muM0PolY)     {rootTree->Branch("muM0PolY",&muM0PolY,"muM0PolY/D");}
    if (store_muM0PolZ)     {rootTree->Branch("muM0PolZ",&muM0PolZ,"muM0PolZ/D");}
    if (store_muM1Time)     {rootTree->Branch("muM1Time",&muM1Time,"muM1Time/D");}
    if (store_muM1PolX)     {rootTree->Branch("muM1PolX",&muM1PolX,"muM1PolX/D");}
    if (store_muM1PolY)     {rootTree->Branch("muM1PolY",&muM1PolY,"muM1PolY/D");}
    if (store_muM1PolZ)     {rootTree->Branch("muM1PolZ",&muM1PolZ,"muM1PolZ/D");}
    if (store_muM2Time)     {rootTree->Branch("muM2Time",&muM2Time,"muM2Time/D");}
    if (store_muM2PolX)     {rootTree->Branch("muM2PolX",&muM2PolX,"muM2PolX/D");}
    if (store_muM2PolY)     {rootTree->Branch("muM2PolY",&muM2PolY,"muM2PolY/D");}
    if (store_muM2PolZ)     {rootTree->Branch("muM2PolZ",&muM2PolZ,"muM2PolZ/D");}
    if (store_posIniMomX)   {rootTree->Branch("posIniMomX",&posIniMomx,"posIniMomX/D");}
    if (store_posIniMomY)   {rootTree->Branch("posIniMomY",&posIniMomy,"posIniMomY/D");}
    if (store_posIniMomZ)   {rootTree->Branch("posIniMomZ",&posIniMomz,"posIniMomZ/D");}
    if (musrParameters::boolG4OpticalPhotons) {
        if (store_nOptPhot)     {rootTree->Branch("nOptPhot",&nOptPhot,"nOptPhot/I");}
        if (store_nOptPhotDet)  {rootTree->Branch("nOptPhotDet",&nOptPhotDet,"nOptPhotDet/I");}
        if (store_phot_time)    {rootTree->Branch("phot_time",&phot_time,"phot_time[nOptPhotDet]/D");}
    }
    //  if (store_globalTime)   {rootTree->Branch("globalTime",&globalTime,"globalTime/D");}
    //  if (store_fieldValue)   {rootTree->Branch("fieldValue",&fieldValue,"fieldValue/D");}
    if (store_fieldNomVal)  {
        rootTree->Branch("nFieldNomVal",&nFieldNomVal,"nFieldNomVal/I");
        rootTree->Branch("fieldNomVal",&fieldNomVal,"fieldNomVal[nFieldNomVal]/D");
    }
    if (store_fieldIntegralBx) {rootTree->Branch("BxIntegral",&BxIntegral,"BxIntegral/D");}
    if (store_fieldIntegralBy) {rootTree->Branch("ByIntegral",&ByIntegral,"ByIntegral/D");}
    if (store_fieldIntegralBz) {rootTree->Branch("BzIntegral",&BzIntegral,"BzIntegral/D");}
    if (store_fieldIntegralBz1) {rootTree->Branch("BzIntegral1",&BzIntegral1,"BzIntegral1/D");}
    if (store_fieldIntegralBz2) {rootTree->Branch("BzIntegral2",&BzIntegral2,"BzIntegral2/D");}
    if (store_fieldIntegralBz3) {rootTree->Branch("BzIntegral3",&BzIntegral3,"BzIntegral3/D");}

    rootTree->Branch("det_n",&det_n,"det_n/I");
    if (store_det_ID)       {rootTree->Branch("det_ID",&det_ID,"det_ID[det_n]/I");}
    if (store_det_edep)     {rootTree->Branch("det_edep",&det_edep,"det_edep[det_n]/D");}
    if (store_det_edep_el)  {rootTree->Branch("det_edep_el",&det_edep_el,"det_edep_el[det_n]/D");}
    if (store_det_edep_pos) {rootTree->Branch("det_edep_pos",&det_edep_pos,"det_edep_pos[det_n]/D");}
    if (store_det_edep_gam) {rootTree->Branch("det_edep_gam",&det_edep_gam,"det_edep_gam[det_n]/D");}
    if (store_det_edep_mup) {rootTree->Branch("det_edep_mup",&det_edep_mup,"det_edep_mup[det_n]/D");}
    if (store_det_kine_mup) {rootTree->Branch("det_kine_mup",&det_kine_mup,"det_kine_mup[det_n]/D");}
    if (store_det_edep_mun) {rootTree->Branch("det_edep_mun",&det_edep_mun,"det_edep_mun[det_n]/D");}
    if (store_det_kine_mun) {rootTree->Branch("det_kine_mun",&det_kine_mun,"det_kine_mun[det_n]/D");}
    if (store_det_nsteps)   {rootTree->Branch("det_nsteps",&det_nsteps,"det_nsteps[det_n]/I");}
    if (store_det_length)   {rootTree->Branch("det_length",&det_length,"det_length[det_n]/D");}
    if (store_det_start)    {rootTree->Branch("det_time_start",&det_time_start,"det_time_start[det_n]/D");}
    if (store_det_end)      {rootTree->Branch("det_time_end",&det_time_end,"det_time_end[det_n]/D");}
    if (store_det_x)        {rootTree->Branch("det_x",&det_x,"det_x[det_n]/D");}
    if (store_det_y)        {rootTree->Branch("det_y",&det_y,"det_y[det_n]/D");}
    if (store_det_z)        {rootTree->Branch("det_z",&det_z,"det_z[det_n]/D");}
    if (store_det_x_mup)        {rootTree->Branch("det_x_mup",&det_x_mup,"det_x_mup[det_n]/D");}
    if (store_det_y_mup)        {rootTree->Branch("det_y_mup",&det_y_mup,"det_y_mup[det_n]/D");}
    if (store_det_z_mup)        {rootTree->Branch("det_z_mup",&det_z_mup,"det_z_mup[det_n]/D");}
    if (store_det_x_mun)        {rootTree->Branch("det_x_mun",&det_x_mun,"det_x_mun[det_n]/D");}
    if (store_det_y_mun)        {rootTree->Branch("det_y_mun",&det_y_mun,"det_y_mun[det_n]/D");}
    if (store_det_z_mun)        {rootTree->Branch("det_z_mun",&det_z_mun,"det_z_mun[det_n]/D");}
    if (store_det_kine)     {rootTree->Branch("det_kine",&det_kine,"det_kine[det_n]/D");}
    if (store_det_VrtxKine) {rootTree->Branch("det_VrtxKine",&det_VrtxKine,"det_VrtxKine[det_n]/D");}
    if (store_det_VrtxX)    {rootTree->Branch("det_VrtxX",&det_VrtxX,"det_VrtxX[det_n]/D");}
    if (store_det_VrtxY)    {rootTree->Branch("det_VrtxY",&det_VrtxY,"det_VrtxY[det_n]/D");}
    if (store_det_VrtxZ)    {rootTree->Branch("det_VrtxZ",&det_VrtxZ,"det_VrtxZ[det_n]/D");}
    if (store_det_VrtxVolID){rootTree->Branch("det_VrtxVolID",&det_VrtxVolID,"det_VrtxVolID[det_n]/I");}
    if (store_det_VrtxProcID){rootTree->Branch("det_VrtxProcID",&det_VrtxProcID,"det_VrtxProcID[det_n]/I");}
    if (store_det_VrtxTrackID){rootTree->Branch("det_VrtxTrackID",&det_VrtxTrackID,"det_VrtxTrackID[det_n]/I");}
    if (store_det_VrtxParentTrackID){rootTree->Branch("det_VrtxPrtTrackID",&det_VrtxPrtTrackID,"det_VrtxPrtTrackID[det_n]/I");}
    if (store_det_VrtxParticleID){rootTree->Branch("det_VrtxParticleID",&det_VrtxParticleID,"det_VrtxParticleID[det_n]/I");}
    if (store_det_VvvKine) {rootTree->Branch("det_VvvKine",&det_VvvKine,"det_VvvKine[det_n]/D");}
    if (store_det_VvvX)    {rootTree->Branch("det_VvvX",&det_VvvX,"det_VvvX[det_n]/D");}
    if (store_det_VvvY)    {rootTree->Branch("det_VvvY",&det_VvvY,"det_VvvY[det_n]/D");}
    if (store_det_VvvZ)    {rootTree->Branch("det_VvvZ",&det_VvvZ,"det_VvvZ[det_n]/D");}
    if (store_det_VvvVolID){rootTree->Branch("det_VvvVolID",&det_VvvVolID,"det_VvvVolID[det_n]/I");}
    if (store_det_VvvProcID){rootTree->Branch("det_VvvProcID",&det_VvvProcID,"det_VvvProcID[det_n]/I");}
    if (store_det_VvvTrackID){rootTree->Branch("det_VvvTrackID",&det_VvvTrackID,"det_VvvTrackID[det_n]/I");}
    if (store_det_VvvParticleID){rootTree->Branch("det_VvvParticleID",&det_VvvParticleID,"det_VvvParticleID[det_n]/I");}

    if (boolIsAnySpecialSaveVolumeDefined) {
        rootTree->Branch("save_n",&save_n,"save_n/I");
        rootTree->Branch("save_detID",&save_detID,"save_detID[save_n]/I");
        rootTree->Branch("save_particleID",&save_particleID,"save_particleID[save_n]/I");
        rootTree->Branch("save_ke",&save_ke,"save_ke[save_n]/D");
        rootTree->Branch("save_time",&save_time,"save_time[save_n]/D");
        rootTree->Branch("save_x",&save_x,"save_x[save_n]/D");
        rootTree->Branch("save_y",&save_y,"save_y[save_n]/D");
        rootTree->Branch("save_z",&save_z,"save_z[save_n]/D");
        rootTree->Branch("save_px",&save_px,"save_px[save_n]/D");
        rootTree->Branch("save_py",&save_py,"save_py[save_n]/D");
        rootTree->Branch("save_pz",&save_pz,"save_pz[save_n]/D");
        rootTree->Branch("save_polx",&save_polx,"save_polx[save_n]/D");
        rootTree->Branch("save_poly",&save_poly,"save_poly[save_n]/D");
        rootTree->Branch("save_polz",&save_polz,"save_polz[save_n]/D");
    }

    if (musrParameters::boolG4OpticalPhotons) {
        if (store_odet_ID || store_odet_nPhot || store_odet_nPhotPrim || store_odet_timeFirst || store_odet_timeSecond ||
            store_odet_timeThird || store_odet_timeA || store_odet_timeB || store_odet_timeC || store_odet_timeD ||
            store_odet_timeMean || store_odet_timeLast || store_odet_timeCFD || store_odet_amplCFD || store_odet_timeC1)
        {rootTree->Branch("odet_n",&odet_n,"odet_n/I");}
        if (store_odet_ID)          {rootTree->Branch("odet_ID",&odet_ID,"odet_ID[odet_n]/I");}
        if (store_odet_nPhot)       {rootTree->Branch("odet_nPhot",&odet_nPhot,"odet_nPhot[odet_n]/I");}
        if (store_odet_nPhotPrim)   {rootTree->Branch("odet_nPhotPrim",&odet_nPhotPrim,"odet_nPhotPrim[odet_n]/I");}
        if (store_odet_timeFirst)   {rootTree->Branch("odet_timeFirst",&odet_timeFirst,"odet_timeFirst[odet_n]/D");}
        if (store_odet_timeSecond)  {rootTree->Branch("odet_timeSecond",&odet_timeSecond,"odet_timeSecond[odet_n]/D");}
        if (store_odet_timeThird)   {rootTree->Branch("odet_timeThird",&odet_timeThird,"odet_timeThird[odet_n]/D");}
        if (store_odet_timeA)       {rootTree->Branch("odet_timeA",&odet_timeA,"odet_timeA[odet_n]/D");}
        if (store_odet_timeB)       {rootTree->Branch("odet_timeB",&odet_timeB,"odet_timeB[odet_n]/D");}
        if (store_odet_timeC)       {rootTree->Branch("odet_timeC",&odet_timeC,"odet_timeC[odet_n]/D");}
        if (store_odet_timeD)       {rootTree->Branch("odet_timeD",&odet_timeD,"odet_timeD[odet_n]/D");}
        if (store_odet_timeMean)       {rootTree->Branch("odet_timeMean",&odet_timeMean,"odet_timeMean[odet_n]/D");}
        if (store_odet_timeLast)    {rootTree->Branch("odet_timeLast",&odet_timeLast,"odet_timeLast[odet_n]/D");}
        if (store_odet_timeCFD)     {rootTree->Branch("odet_timeCFD",&odet_timeCFD,"odet_timeCFD[odet_n]/D");}
        if (store_odet_amplCFD)     {rootTree->Branch("odet_amplCFD",&odet_amplCFD,"odet_amplCFD[odet_n]/D");}
        if (store_odet_timeCFDarray)     {
            rootTree->Branch("odet_timeCFD100",&odet_timeCFD100,"odet_timeCFD100[odet_n]/D");
            rootTree->Branch("odet_timeCFD101",&odet_timeCFD101,"odet_timeCFD101[odet_n]/D");
            rootTree->Branch("odet_timeCFD102",&odet_timeCFD102,"odet_timeCFD102[odet_n]/D");
            rootTree->Branch("odet_timeCFD103",&odet_timeCFD103,"odet_timeCFD103[odet_n]/D");
            rootTree->Branch("odet_timeCFD104",&odet_timeCFD104,"odet_timeCFD104[odet_n]/D");
            rootTree->Branch("odet_timeCFD105",&odet_timeCFD105,"odet_timeCFD105[odet_n]/D");
            rootTree->Branch("odet_timeCFD106",&odet_timeCFD106,"odet_timeCFD106[odet_n]/D");
            rootTree->Branch("odet_timeCFD107",&odet_timeCFD107,"odet_timeCFD107[odet_n]/D");
            rootTree->Branch("odet_timeCFD108",&odet_timeCFD108,"odet_timeCFD108[odet_n]/D");
            rootTree->Branch("odet_timeCFD109",&odet_timeCFD109,"odet_timeCFD109[odet_n]/D");
            rootTree->Branch("odet_timeCFD110",&odet_timeCFD110,"odet_timeCFD110[odet_n]/D");
            rootTree->Branch("odet_timeCFD111",&odet_timeCFD111,"odet_timeCFD111[odet_n]/D");
            rootTree->Branch("odet_timeCFD112",&odet_timeCFD112,"odet_timeCFD112[odet_n]/D");
            rootTree->Branch("odet_timeCFD200",&odet_timeCFD200,"odet_timeCFD200[odet_n]/D");
            rootTree->Branch("odet_timeCFD201",&odet_timeCFD201,"odet_timeCFD201[odet_n]/D");
            rootTree->Branch("odet_timeCFD202",&odet_timeCFD202,"odet_timeCFD202[odet_n]/D");
            rootTree->Branch("odet_timeCFD203",&odet_timeCFD203,"odet_timeCFD203[odet_n]/D");
            rootTree->Branch("odet_timeCFD204",&odet_timeCFD204,"odet_timeCFD204[odet_n]/D");
            rootTree->Branch("odet_timeCFD205",&odet_timeCFD205,"odet_timeCFD205[odet_n]/D");
            rootTree->Branch("odet_timeCFD206",&odet_timeCFD206,"odet_timeCFD206[odet_n]/D");
            rootTree->Branch("odet_timeCFD207",&odet_timeCFD207,"odet_timeCFD207[odet_n]/D");
            rootTree->Branch("odet_timeCFD208",&odet_timeCFD208,"odet_timeCFD208[odet_n]/D");
            rootTree->Branch("odet_timeCFD209",&odet_timeCFD209,"odet_timeCFD209[odet_n]/D");
            rootTree->Branch("odet_timeCFD210",&odet_timeCFD210,"odet_timeCFD210[odet_n]/D");
            rootTree->Branch("odet_timeCFD211",&odet_timeCFD211,"odet_timeCFD211[odet_n]/D");
            rootTree->Branch("odet_timeCFD212",&odet_timeCFD212,"odet_timeCFD212[odet_n]/D");
            rootTree->Branch("odet_timeCFD300",&odet_timeCFD300,"odet_timeCFD300[odet_n]/D");
            rootTree->Branch("odet_timeCFD301",&odet_timeCFD301,"odet_timeCFD301[odet_n]/D");
            rootTree->Branch("odet_timeCFD302",&odet_timeCFD302,"odet_timeCFD302[odet_n]/D");
            rootTree->Branch("odet_timeCFD303",&odet_timeCFD303,"odet_timeCFD303[odet_n]/D");
            rootTree->Branch("odet_timeCFD304",&odet_timeCFD304,"odet_timeCFD304[odet_n]/D");
            rootTree->Branch("odet_timeCFD305",&odet_timeCFD305,"odet_timeCFD305[odet_n]/D");
            rootTree->Branch("odet_timeCFD306",&odet_timeCFD306,"odet_timeCFD306[odet_n]/D");
            rootTree->Branch("odet_timeCFD307",&odet_timeCFD307,"odet_timeCFD307[odet_n]/D");
            rootTree->Branch("odet_timeCFD308",&odet_timeCFD308,"odet_timeCFD308[odet_n]/D");
            rootTree->Branch("odet_timeCFD309",&odet_timeCFD309,"odet_timeCFD309[odet_n]/D");
            rootTree->Branch("odet_timeCFD310",&odet_timeCFD310,"odet_timeCFD310[odet_n]/D");
            rootTree->Branch("odet_timeCFD311",&odet_timeCFD311,"odet_timeCFD311[odet_n]/D");
            rootTree->Branch("odet_timeCFD312",&odet_timeCFD312,"odet_timeCFD312[odet_n]/D");
            rootTree->Branch("odet_timeCFD400",&odet_timeCFD400,"odet_timeCFD400[odet_n]/D");
            rootTree->Branch("odet_timeCFD401",&odet_timeCFD401,"odet_timeCFD401[odet_n]/D");
            rootTree->Branch("odet_timeCFD402",&odet_timeCFD402,"odet_timeCFD402[odet_n]/D");
            rootTree->Branch("odet_timeCFD403",&odet_timeCFD403,"odet_timeCFD403[odet_n]/D");
            rootTree->Branch("odet_timeCFD404",&odet_timeCFD404,"odet_timeCFD404[odet_n]/D");
            rootTree->Branch("odet_timeCFD405",&odet_timeCFD405,"odet_timeCFD405[odet_n]/D");
            rootTree->Branch("odet_timeCFD406",&odet_timeCFD406,"odet_timeCFD406[odet_n]/D");
            rootTree->Branch("odet_timeCFD407",&odet_timeCFD407,"odet_timeCFD407[odet_n]/D");
            rootTree->Branch("odet_timeCFD408",&odet_timeCFD408,"odet_timeCFD408[odet_n]/D");
            rootTree->Branch("odet_timeCFD409",&odet_timeCFD409,"odet_timeCFD409[odet_n]/D");
            rootTree->Branch("odet_timeCFD410",&odet_timeCFD410,"odet_timeCFD410[odet_n]/D");
            rootTree->Branch("odet_timeCFD411",&odet_timeCFD411,"odet_timeCFD411[odet_n]/D");
            rootTree->Branch("odet_timeCFD412",&odet_timeCFD412,"odet_timeCFD412[odet_n]/D");
            rootTree->Branch("odet_timeCFD500",&odet_timeCFD500,"odet_timeCFD500[odet_n]/D");
            rootTree->Branch("odet_timeCFD501",&odet_timeCFD501,"odet_timeCFD501[odet_n]/D");
            rootTree->Branch("odet_timeCFD502",&odet_timeCFD502,"odet_timeCFD502[odet_n]/D");
            rootTree->Branch("odet_timeCFD503",&odet_timeCFD503,"odet_timeCFD503[odet_n]/D");
            rootTree->Branch("odet_timeCFD504",&odet_timeCFD504,"odet_timeCFD504[odet_n]/D");
            rootTree->Branch("odet_timeCFD505",&odet_timeCFD505,"odet_timeCFD505[odet_n]/D");
            rootTree->Branch("odet_timeCFD506",&odet_timeCFD506,"odet_timeCFD506[odet_n]/D");
            rootTree->Branch("odet_timeCFD507",&odet_timeCFD507,"odet_timeCFD507[odet_n]/D");
            rootTree->Branch("odet_timeCFD508",&odet_timeCFD508,"odet_timeCFD508[odet_n]/D");
            rootTree->Branch("odet_timeCFD509",&odet_timeCFD509,"odet_timeCFD509[odet_n]/D");
            rootTree->Branch("odet_timeCFD510",&odet_timeCFD510,"odet_timeCFD510[odet_n]/D");
            rootTree->Branch("odet_timeCFD511",&odet_timeCFD511,"odet_timeCFD511[odet_n]/D");
            rootTree->Branch("odet_timeCFD512",&odet_timeCFD512,"odet_timeCFD512[odet_n]/D");
        }
        if (store_odet_timeC1)       {
            rootTree->Branch("odet_timeC1",&odet_timeC1,"odet_timeC1[odet_n]/D");
            rootTree->Branch("odet_timeC2",&odet_timeC2,"odet_timeC2[odet_n]/D");
            rootTree->Branch("odet_timeC3",&odet_timeC3,"odet_timeC3[odet_n]/D");
            rootTree->Branch("odet_timeC4",&odet_timeC4,"odet_timeC4[odet_n]/D");
            rootTree->Branch("odet_timeC5",&odet_timeC5,"odet_timeC5[odet_n]/D");
            rootTree->Branch("odet_timeC6",&odet_timeC6,"odet_timeC6[odet_n]/D");
            rootTree->Branch("odet_timeC7",&odet_timeC7,"odet_timeC7[odet_n]/D");
            rootTree->Branch("odet_timeC8",&odet_timeC8,"odet_timeC8[odet_n]/D");
            rootTree->Branch("odet_timeC9",&odet_timeC9,"odet_timeC9[odet_n]/D");
            rootTree->Branch("odet_timeC10",&odet_timeC10,"odet_timeC10[odet_n]/D");
            rootTree->Branch("odet_timeC11",&odet_timeC11,"odet_timeC11[odet_n]/D");
            rootTree->Branch("odet_timeC12",&odet_timeC12,"odet_timeC12[odet_n]/D");
            rootTree->Branch("odet_timeC13",&odet_timeC13,"odet_timeC13[odet_n]/D");
            rootTree->Branch("odet_timeC14",&odet_timeC14,"odet_timeC14[odet_n]/D");
            rootTree->Branch("odet_timeC15",&odet_timeC15,"odet_timeC15[odet_n]/D");
            rootTree->Branch("odet_timeC16",&odet_timeC16,"odet_timeC16[odet_n]/D");
            rootTree->Branch("odet_timeC17",&odet_timeC17,"odet_timeC17[odet_n]/D");
            rootTree->Branch("odet_timeC18",&odet_timeC18,"odet_timeC18[odet_n]/D");
            rootTree->Branch("odet_timeC19",&odet_timeC19,"odet_timeC19[odet_n]/D");
            rootTree->Branch("odet_timeC20",&odet_timeC20,"odet_timeC20[odet_n]/D");
            rootTree->Branch("odet_timeC21",&odet_timeC21,"odet_timeC21[odet_n]/D");
            rootTree->Branch("odet_timeC22",&odet_timeC22,"odet_timeC22[odet_n]/D");
            rootTree->Branch("odet_timeC23",&odet_timeC23,"odet_timeC23[odet_n]/D");
            rootTree->Branch("odet_timeC24",&odet_timeC24,"odet_timeC24[odet_n]/D");
            rootTree->Branch("odet_timeC25",&odet_timeC25,"odet_timeC25[odet_n]/D");
            rootTree->Branch("odet_timeC26",&odet_timeC26,"odet_timeC26[odet_n]/D");
            rootTree->Branch("odet_timeC27",&odet_timeC27,"odet_timeC27[odet_n]/D");
            rootTree->Branch("odet_timeC28",&odet_timeC28,"odet_timeC28[odet_n]/D");
            rootTree->Branch("odet_timeC29",&odet_timeC29,"odet_timeC29[odet_n]/D");
            rootTree->Branch("odet_timeC30",&odet_timeC30,"odet_timeC30[odet_n]/D");
            rootTree->Branch("odet_timeC31",&odet_timeC31,"odet_timeC31[odet_n]/D");
            rootTree->Branch("odet_timeC32",&odet_timeC32,"odet_timeC32[odet_n]/D");
            rootTree->Branch("odet_timeC33",&odet_timeC33,"odet_timeC33[odet_n]/D");
            rootTree->Branch("odet_timeC34",&odet_timeC34,"odet_timeC34[odet_n]/D");
            rootTree->Branch("odet_timeC35",&odet_timeC35,"odet_timeC35[odet_n]/D");
            rootTree->Branch("odet_timeC36",&odet_timeC36,"odet_timeC36[odet_n]/D");
            rootTree->Branch("odet_timeC37",&odet_timeC37,"odet_timeC37[odet_n]/D");
            rootTree->Branch("odet_timeC38",&odet_timeC38,"odet_timeC38[odet_n]/D");
            rootTree->Branch("odet_timeC39",&odet_timeC39,"odet_timeC39[odet_n]/D");
            rootTree->Branch("odet_timeC40",&odet_timeC40,"odet_timeC40[odet_n]/D");
        }
    }
    //  htest1 = new TH1F("htest1","The debugging histogram 1",50,-4.,4.);
    //  htest2 = new TH1F("htest2","The debugging histogram 2",50,0.,3.142);
    htest1 = new TH2F("htest1","x, y",50,-200.,200.,50,-200.,200.);
    htest2 = new TH2F("htest2","R, z",50,0.,250.,50,-150.,150.);
    htest3 = new TH1F("htest3","Energy in MeV",55,0.,55.);
    htest4 = new TH1F("htest4","Radioactive electron kinetic energy",250,0.,2.5);
    htest5 = new TH1F("htest5","The debugging histogram 5",50,-4.,4.);
    htest6 = new TH1F("htest6","The debugging histogram 6",50,0.,3.142);
    htest7 = new TH1F("htest7","The debugging histogram 7",100000,0.,100.);
    htest8 = new TH1F("htest8","The debugging histogram 8",100000,0.,100.);

    G4cout << "musrRootOutput::BeginOfRunAction()  The Root tree and branches were defined."<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void musrRootOutput::EndOfRunAction() {
    G4cout << "musrRootOutput::EndOfRunAction() - Writing out the Root tree:"<<G4endl;
    rootTree->Write();
    htest1->Write();
    htest2->Write();
    htest3->Write();
    htest4->Write();
    htest5->Write();
    htest6->Write();
    htest7->Write();
    //-----------------------------------------------------------------------------------------
    //  Uncoment for iterative musrSim runs (e.g. when searching for a quadrupole triplet focus using a python script)
    //  //  std::cout<<"DEBUG: FOCUS TEST: sigma="<<(htest7->GetMean())/(htest7->GetEntries())<<std::endl;
    //  //  std::cout<<"DEBUG: FOCUS TEST: sigma="<<(htest7->GetMean())<<std::endl;
    //  Double_t xRMS_tmp = htest7->GetRMS();
    //  Double_t yRMS_tmp = htest8->GetRMS();
    //  std::cout<<"DEBUG: FOCUS TEST: sigma="<< sqrt(xRMS_tmp*xRMS_tmp + yRMS_tmp*yRMS_tmp) <<std::endl;
    //  htest8->Write();
    //-----------------------------------------------------------------------------------------
    //
    // Variables exported from Geant simulation to the Root output
    //  static const Int_t nGeantParamD=10;
    TVectorD TVector_GeantParametersD(maxNGeantParameters);
    for (Int_t i=0; i<maxNGeantParameters; i++) {
        TVector_GeantParametersD[i]=GeantParametersD[i];
    }
    //
    TH1D* hGeantParameters = new TH1D("hGeantParameters","hGeantParameters",maxNGeantParameters,0.,float(maxNGeantParameters));
    TVector_GeantParametersD.Write("geantParametersD");
    for (Int_t i=0; i<maxNGeantParameters; i++) {
        hGeantParameters->SetBinContent(i,GeantParametersD[i]);
    }
    hGeantParameters->Write();
    rootFile->Close();
    G4cout<<"musrRootOutput::EndOfRunAction() - Root tree written out."<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrRootOutput::FillEvent() {
    htest5->Fill(atan2(posIniMomy,posIniMomx));
    htest6->Fill(atan2(sqrt(posIniMomx*posIniMomx+posIniMomy*posIniMomy),posIniMomz));
    if (weight>0.) {
        if ( !((musrParameters::storeOnlyEventsWithHits)&&(det_n<=0)&&(odet_n<=0)) ) {
            rootTree->Fill();
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrRootOutput::ClearAllRootVariables() {
    runID=-1000;
    eventID=-1000;
    weight=1.;
    timeToNextEvent = -1000;
    B_t[0]=-1000.;B_t[1]=-1000.;B_t[2]=-1000.;B_t[3]=-1000.;B_t[4]=-1000.;B_t[5]=-1000.;
    muIniTime=-1000;
    muIniPosX=-1000; muIniPosY=-1000; muIniPosZ=-1000;
    muIniMomX=-1000; muIniMomY=-1000; muIniMomZ=-1000;
    muIniPolX=-1000; muIniPolY=-1000; muIniPolZ=-1000;
    muDecayDetID=-1000;
    muDecayPolX=-1000; muDecayPolY=-1000; muDecayPolZ=-1000;
    muTargetTime=-1000; muTargetPolX=-1000; muTargetPolY=-1000; muTargetPolZ=-1000;
    muTargetMomX=-1000; muTargetMomY=-1000; muTargetMomZ=-1000;
    muM0Time=-1000; muM0PolX=-1000; muM0PolY=-1000; muM0PolZ=-1000;
    muM1Time=-1000; muM1PolX=-1000; muM1PolY=-1000; muM1PolZ=-1000;
    muM2Time=-1000; muM2PolX=-1000; muM2PolY=-1000; muM2PolZ=-1000;
    muDecayPosX=-1000;muDecayPosY=-1000;muDecayPosZ=-1000;
    muDecayTime=-1000;
    posIniMomx=-1000;posIniMomy=-1000;posIniMomz=-1000;
    nOptPhot=0;
    nOptPhotDet=0;
    BxIntegral = -1000; ByIntegral = -1000; BzIntegral = -1000;
    BzIntegral1 = -1000; BzIntegral2 = -1000; BzIntegral3 = -1000;
    det_n=0;
    save_n=0;
    odet_n=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void musrRootOutput::SetVolumeIDMapping(std::string logivol, int volumeID) {
    //  This function assigns a unique number to each sensitive detector name.
    //  The numbers are used in the root tree, as it is easier to work with numbers
    //  rather than with strings.
    if (SensDetectorMapping[logivol]) {
        char message[200];
        sprintf(message,"musrRootOutput::SetVolumeIDMapping: Sensitive volume %s already assigned",logivol.c_str());
        musrErrorMessage::GetInstance()->musrError(FATAL,message,false);
    }
    else{
        SensDetectorMapping[logivol]=volumeID;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int musrRootOutput::ConvertVolumeToID(std::string logivol) {
    G4int volumeID = SensDetectorMapping[logivol];
    if (volumeID==0) {
        char message[200];
        sprintf(message,"musrRootOutput::ConvertVolumeToID: No ID number assigned to sensitive volume %s .",logivol.c_str());
        musrErrorMessage::GetInstance()->musrError(SERIOUS,message,true);
    }
    return volumeID;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int musrRootOutput::ConvertProcessToID(std::string processName) {
    G4int processID = ProcessIDMapping[processName];
    if (processID==0) {
        char message[200];
        sprintf(message,"musrRootOutput::ConvertProcessToID: No ID number assigned to the process \"%s\" .",processName.c_str());
        musrErrorMessage::GetInstance()->musrError(WARNING,message,true);
    }
    return processID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrRootOutput::SetSaveDetectorInfo (G4int ID, G4int particleID, G4double ke,
                                          G4double x, G4double y, G4double z, G4double time,
                                          G4double px, G4double py, G4double pz, G4double polx, G4double poly, G4double polz) {
    if (save_n>=save_nMax) {
        char message[200];
        sprintf(message,"musrRootOutput.cc::SetSaveDetectorInfo(): more \"save\" hits then allowed: save_nMax=%i",save_nMax);
        musrErrorMessage::GetInstance()->musrError(SERIOUS,message,true);
    }
    else {
        save_detID[save_n] = ID;
        save_particleID[save_n] = particleID;
        save_ke[save_n] = ke / CLHEP::MeV;
        save_x[save_n] = x / CLHEP::mm;
        save_y[save_n] = y / CLHEP::mm;
        save_z[save_n] = z / CLHEP::mm;
        save_time[save_n] = time / CLHEP::microsecond;
        save_px[save_n] = px / CLHEP::MeV;
        save_py[save_n] = py / CLHEP::MeV;
        save_pz[save_n] = pz / CLHEP::MeV;
        save_polx[save_n] = polx;
        save_poly[save_n] = poly;
        save_polz[save_n] = polz;
        save_n++;
    }
}

void musrRootOutput::SetFieldNomVal(G4int i, G4double value) {
    if (i<maxNFieldnNominalValues) {
        // cks the following will probably not be correct for electric field,
        //     because the units are tesla.  Should be modified.
        fieldNomVal[i]=value/CLHEP::tesla;
    }
    else {
        char message[200];
        sprintf(message,
                "musrRootOutput.cc::SetFieldNomVal(): more electromagnetic fields then allowed: maxNFieldnNominalValues=%i",
                maxNFieldnNominalValues);
        musrErrorMessage::GetInstance()->musrError(SERIOUS,message,true);
    }
}


void musrRootOutput::SetDetectorInfo (G4int nDetectors, G4int ID, G4int particleID, G4double edep,
                                      G4double edep_el, G4double edep_pos,
                                      G4double edep_gam, G4double edep_mup, G4double kine_mup, G4double edep_mun, G4double kine_mun, G4int nsteps, G4double length, G4double t1,
                                      G4double t2, G4double x, G4double y, G4double z,
                                      G4double x_mup, G4double y_mup, G4double z_mup, G4double x_mun, G4double y_mun, G4double z_mun,
                                      G4double ek, G4double ekVertex, G4double xVertex, G4double yVertex, G4double zVertex,
                                      G4int idVolVertex, G4int idProcVertex, G4int idTrackVertex, G4int idPrtTrackVertex)
{
    if ((nDetectors<0)||(nDetectors>=(det_nMax-1))) {
        char message[200];
        sprintf(message,"musrRootOutput.cc::SetDetectorInfo: nDetectors %i is larger than det_nMax = %i",nDetectors,det_nMax);
        musrErrorMessage::GetInstance()->musrError(SERIOUS,message,false);
        return;
    }
    else {
        det_n=nDetectors+1;
        det_ID[nDetectors]=ID;
        det_edep[nDetectors]=edep/CLHEP::MeV;
        det_edep_el[nDetectors]=edep_el/CLHEP::MeV;
        det_edep_pos[nDetectors]=edep_pos/CLHEP::MeV;
        det_edep_gam[nDetectors]=edep_gam/CLHEP::MeV;
        det_edep_mup[nDetectors]=edep_mup/CLHEP::MeV;
        det_kine_mup[nDetectors]=kine_mup/CLHEP::MeV;
        det_edep_mun[nDetectors]=edep_mun/CLHEP::MeV;
        det_kine_mun[nDetectors]=kine_mun/CLHEP::MeV;
        det_nsteps[nDetectors]=nsteps;
        det_length[nDetectors]=length/CLHEP::mm;
        det_time_start[nDetectors]=t1/CLHEP::microsecond;
        det_time_end[nDetectors]=t2/CLHEP::microsecond;
        det_x[nDetectors]=x/CLHEP::mm;
        det_y[nDetectors]=y/CLHEP::mm;
        det_z[nDetectors]=z/CLHEP::mm;
        det_x_mup[nDetectors]=x_mup/CLHEP::mm;
        det_y_mup[nDetectors]=y_mup/CLHEP::mm;
        det_z_mup[nDetectors]=z_mup/CLHEP::mm;
        det_x_mun[nDetectors]=x_mun/CLHEP::mm;
        det_y_mun[nDetectors]=y_mun/CLHEP::mm;
        det_z_mun[nDetectors]=z_mun/CLHEP::mm;
        det_kine[nDetectors]=ek/CLHEP::MeV;
        det_VrtxKine[nDetectors]=ekVertex/CLHEP::MeV;
        det_VrtxX[nDetectors]=xVertex/CLHEP::mm;
        det_VrtxY[nDetectors]=yVertex/CLHEP::mm;
        det_VrtxZ[nDetectors]=zVertex/CLHEP::mm;
        det_VrtxVolID[nDetectors]=idVolVertex;
        det_VrtxProcID[nDetectors]=idProcVertex;
        det_VrtxPrtTrackID[nDetectors]=idPrtTrackVertex;
        det_VrtxTrackID[nDetectors]=idTrackVertex;
        det_VrtxParticleID[nDetectors]=particleID;
    }
}

void musrRootOutput::SetDetectorInfoVvv (G4int nDetectors,
                                         G4double ekVertex, G4double xVertex, G4double yVertex, G4double zVertex,
                                         G4int idVolVertex, G4int idProcVertex, G4int idTrackVertex, G4int particleID)   {
    if ((nDetectors<0)||(nDetectors>=(det_nMax-1))) {
        char message[200];
        sprintf(message,"musrRootOutput.cc::SetDetectorInfoVvv: nDetectors %i is larger than det_nMax = %i",nDetectors,det_nMax);
        musrErrorMessage::GetInstance()->musrError(SERIOUS,message,false);
        return;
    }
    else {
        det_VvvKine[nDetectors]=ekVertex/CLHEP::MeV;
        det_VvvX[nDetectors]=xVertex/CLHEP::mm;
        det_VvvY[nDetectors]=yVertex/CLHEP::mm;
        det_VvvZ[nDetectors]=zVertex/CLHEP::mm;
        det_VvvVolID[nDetectors]=idVolVertex;
        det_VvvProcID[nDetectors]=idProcVertex;
        det_VvvTrackID[nDetectors]=idTrackVertex;
        det_VvvParticleID[nDetectors]=particleID;
    }
}


void musrRootOutput::SetOPSAinfo    (G4int nDetectors, G4int ID, G4int nPhot, G4int nPhotPrim, G4double timeFirst,
                                     G4double timeSecond, G4double timeThird, G4double timeA, G4double timeB, G4double timeC,
                                     G4double timeD, G4double timeMean, G4double timeLast, G4double timeCFD, G4double amplCFD)
{
    if ((nDetectors<0)||(nDetectors>=(odet_nMax-1))) {
        char message[200];
        sprintf(message,"musrRootOutput.cc::SetOPSAInfo: nDetectors %i is larger than det_nMax = %i",nDetectors,det_nMax);
        musrErrorMessage::GetInstance()->musrError(SERIOUS,message,false);
        return;
    }
    else {
        odet_n=nDetectors+1;
        odet_ID[nDetectors]=ID;
        odet_nPhot[nDetectors]=nPhot;
        odet_nPhotPrim[nDetectors]=nPhotPrim;
        odet_timeFirst[nDetectors]=timeFirst/CLHEP::microsecond;
        odet_timeSecond[nDetectors]=timeSecond/CLHEP::microsecond;
        odet_timeThird[nDetectors]=timeThird/CLHEP::microsecond;
        odet_timeA[nDetectors]=timeA/CLHEP::microsecond;
        odet_timeB[nDetectors]=timeB/CLHEP::microsecond;
        odet_timeC[nDetectors]=timeC/CLHEP::microsecond;
        odet_timeD[nDetectors]=timeD/CLHEP::microsecond;
        odet_timeMean[nDetectors]=timeMean/CLHEP::microsecond;
        odet_timeLast[nDetectors]=timeLast/CLHEP::microsecond;
        odet_timeCFD[nDetectors]=timeCFD/CLHEP::microsecond;
        odet_amplCFD[nDetectors]=amplCFD;
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void musrRootOutput::SetCFDSpecialInfo (G4int n, G4double time) {
    if      (n==100) {odet_timeCFD100[odet_n-1] = time/CLHEP::microsecond;} // G4cout<<"OKKKKK odet_n-1="<<odet_n-1<<G4endl;}
    else if (n==101) {odet_timeCFD101[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==102) {odet_timeCFD102[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==103) {odet_timeCFD103[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==104) {odet_timeCFD104[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==105) {odet_timeCFD105[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==106) {odet_timeCFD106[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==107) {odet_timeCFD107[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==108) {odet_timeCFD108[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==109) {odet_timeCFD109[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==110) {odet_timeCFD110[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==111) {odet_timeCFD111[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==112) {odet_timeCFD112[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==200) {odet_timeCFD200[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==201) {odet_timeCFD201[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==202) {odet_timeCFD202[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==203) {odet_timeCFD203[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==204) {odet_timeCFD204[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==205) {odet_timeCFD205[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==206) {odet_timeCFD206[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==207) {odet_timeCFD207[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==208) {odet_timeCFD208[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==209) {odet_timeCFD209[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==210) {odet_timeCFD210[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==211) {odet_timeCFD211[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==212) {odet_timeCFD212[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==300) {odet_timeCFD300[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==301) {odet_timeCFD301[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==302) {odet_timeCFD302[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==303) {odet_timeCFD303[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==304) {odet_timeCFD304[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==305) {odet_timeCFD305[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==306) {odet_timeCFD306[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==307) {odet_timeCFD307[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==308) {odet_timeCFD308[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==309) {odet_timeCFD309[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==310) {odet_timeCFD310[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==311) {odet_timeCFD311[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==312) {odet_timeCFD312[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==400) {odet_timeCFD400[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==401) {odet_timeCFD401[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==402) {odet_timeCFD402[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==403) {odet_timeCFD403[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==404) {odet_timeCFD404[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==405) {odet_timeCFD405[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==406) {odet_timeCFD406[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==407) {odet_timeCFD407[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==408) {odet_timeCFD408[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==409) {odet_timeCFD409[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==410) {odet_timeCFD410[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==411) {odet_timeCFD411[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==412) {odet_timeCFD412[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==500) {odet_timeCFD500[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==501) {odet_timeCFD501[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==502) {odet_timeCFD502[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==503) {odet_timeCFD503[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==504) {odet_timeCFD504[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==505) {odet_timeCFD505[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==506) {odet_timeCFD506[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==507) {odet_timeCFD507[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==508) {odet_timeCFD508[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==509) {odet_timeCFD509[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==510) {odet_timeCFD510[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==511) {odet_timeCFD511[odet_n-1] = time/CLHEP::microsecond;}
    else if (n==512) {odet_timeCFD512[odet_n-1] = time/CLHEP::microsecond;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void musrRootOutput::SetTimeC1SpecialInfo (G4double* time) {
    odet_timeC1[odet_n-1] = time[1]/CLHEP::microsecond;
    odet_timeC2[odet_n-1] = time[2]/CLHEP::microsecond;
    odet_timeC3[odet_n-1] = time[3]/CLHEP::microsecond;
    odet_timeC4[odet_n-1] = time[4]/CLHEP::microsecond;
    odet_timeC5[odet_n-1] = time[5]/CLHEP::microsecond;
    odet_timeC6[odet_n-1] = time[6]/CLHEP::microsecond;
    odet_timeC7[odet_n-1] = time[7]/CLHEP::microsecond;
    odet_timeC8[odet_n-1] = time[8]/CLHEP::microsecond;
    odet_timeC9[odet_n-1] = time[9]/CLHEP::microsecond;
    odet_timeC10[odet_n-1] = time[10]/CLHEP::microsecond;
    odet_timeC11[odet_n-1] = time[11]/CLHEP::microsecond;
    odet_timeC12[odet_n-1] = time[12]/CLHEP::microsecond;
    odet_timeC13[odet_n-1] = time[13]/CLHEP::microsecond;
    odet_timeC14[odet_n-1] = time[14]/CLHEP::microsecond;
    odet_timeC15[odet_n-1] = time[15]/CLHEP::microsecond;
    odet_timeC16[odet_n-1] = time[16]/CLHEP::microsecond;
    odet_timeC17[odet_n-1] = time[17]/CLHEP::microsecond;
    odet_timeC18[odet_n-1] = time[18]/CLHEP::microsecond;
    odet_timeC19[odet_n-1] = time[19]/CLHEP::microsecond;
    odet_timeC20[odet_n-1] = time[20]/CLHEP::microsecond;
    odet_timeC21[odet_n-1] = time[21]/CLHEP::microsecond;
    odet_timeC22[odet_n-1] = time[22]/CLHEP::microsecond;
    odet_timeC23[odet_n-1] = time[23]/CLHEP::microsecond;
    odet_timeC24[odet_n-1] = time[24]/CLHEP::microsecond;
    odet_timeC25[odet_n-1] = time[25]/CLHEP::microsecond;
    odet_timeC26[odet_n-1] = time[26]/CLHEP::microsecond;
    odet_timeC27[odet_n-1] = time[27]/CLHEP::microsecond;
    odet_timeC28[odet_n-1] = time[28]/CLHEP::microsecond;
    odet_timeC29[odet_n-1] = time[29]/CLHEP::microsecond;
    odet_timeC30[odet_n-1] = time[30]/CLHEP::microsecond;
    odet_timeC31[odet_n-1] = time[31]/CLHEP::microsecond;
    odet_timeC32[odet_n-1] = time[32]/CLHEP::microsecond;
    odet_timeC33[odet_n-1] = time[33]/CLHEP::microsecond;
    odet_timeC34[odet_n-1] = time[34]/CLHEP::microsecond;
    odet_timeC35[odet_n-1] = time[35]/CLHEP::microsecond;
    odet_timeC36[odet_n-1] = time[36]/CLHEP::microsecond;
    odet_timeC37[odet_n-1] = time[37]/CLHEP::microsecond;
    odet_timeC38[odet_n-1] = time[38]/CLHEP::microsecond;
    odet_timeC39[odet_n-1] = time[39]/CLHEP::microsecond;
    odet_timeC40[odet_n-1] = time[40]/CLHEP::microsecond;

    //  std::cout<<"odet_timeC1["<<odet_n-1<<"]="<<odet_timeC1[odet_n-1]<<std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void musrRootOutput::setRootOutputDirectoryName(char dirName[1000]) {
    strcpy(rootOutputDirectoryName,dirName);
    char message[200];
    sprintf(message,"musrRootOutput.cc::setRootOutputDirectoryName: Root output file will be stored in directory %s",dirName);
    musrErrorMessage::GetInstance()->musrError(INFO,message,false);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void musrRootOutput::SetPhotDetTime(G4double time) {
    if (nOptPhotDet<maxNOptPhotDet) {
        phot_time[nOptPhotDet] = time;
        nOptPhotDet++;
    }
    else {
        nOptPhotDet++; // still want to know how many in total (JSL)
        char message[200];
        sprintf(message,"musrRootOutput.cc::SetPhotDetTime: number of individual photons larger than maxNOptPhotDet (=%d)",maxNOptPhotDet);
        musrErrorMessage::GetInstance()->musrError(WARNING,message,true); // had silent=false and printed all messages (JSL)
    }
}

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

#ifndef musrRootOutput_h
#define musrRootOutput_h 1
//#include "G4UserRunAction.hh"
#include <CLHEP/Units/PhysicalConstants.h>
#include "globals.hh"
#include "G4ThreeVector.hh"
//  ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
//
#include <map>

#include "G4ios.hh"
#include <fstream>
#include <vector>
#include <cmath>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class musrRootOutput  {
  public:
    musrRootOutput();
   ~musrRootOutput();
    static musrRootOutput* GetRootInstance();

  public:
    void BeginOfRunAction();
    void EndOfRunAction();
    void FillEvent();
    void ClearAllRootVariables();
    void SetVolumeIDMapping(std::string logivol, int volumeID);
    G4int ConvertVolumeToID(std::string logivol);
    G4int ConvertProcessToID(std::string processName);
    void SetSpecialSaveVolumeDefined() {boolIsAnySpecialSaveVolumeDefined=true;};

  // Getting variables (just for debugging)
    G4double GetDecayPositionZ() {return muDecayPosZ;};
    G4double GetDecayTime()      {return muDecayTime*CLHEP::microsecond;};
    G4double GetTimeInTarget()   {return muTargetTime*CLHEP::microsecond;};

  // Setting variables common to the whole event:
    void SetRunID          (G4int id) {runID = id;};
    void SetEventID        (G4int id) {eventID = id;};
    void SetTimeToNextEvent(G4double deltaT) {timeToNextEvent = deltaT/CLHEP::microsecond;}
    void SetDecayDetectorID (std::string detectorName) {muDecayDetID = SensDetectorMapping[detectorName];};
    Int_t GetDecayDetectorID() {return  muDecayDetID;};
    void SetBField       (G4double F[6]) {B_t[0]=F[0]/CLHEP::tesla; B_t[1]=F[1]/CLHEP::tesla; B_t[2]=F[2]/CLHEP::tesla; 
                                          B_t[3]=F[3]/CLHEP::tesla; B_t[4]=F[4]/CLHEP::tesla; B_t[5]=F[5]/CLHEP::tesla;};
    void SetDecayPolarisation (G4ThreeVector pol) {muDecayPolX=pol.x(); muDecayPolY=pol.y(); muDecayPolZ=pol.z();};
    void SetDecayPosition (G4ThreeVector pos) {muDecayPosX=pos.x()/CLHEP::mm; muDecayPosY=pos.y()/CLHEP::mm; 
                                               muDecayPosZ=pos.z()/CLHEP::mm;};
    void SetEventWeight  (G4double w) {weight *= w;}
    void SetDetectorInfo (G4int nDetectors, G4int ID, G4int particleID, G4double edep, 
                          G4double edep_el, G4double edep_pos, 
			  G4double edep_gam, G4double edep_mup,G4int nsteps, G4double length, G4double t1, 
			  G4double t2, G4double x, G4double y, G4double z,
			  G4double ek, G4double ekVertex, G4double xVertex, G4double yVertex, G4double zVertex, 
			  G4int idVolVertex, G4int idProcVertex, G4int idTrackVertex) ;

    void SetDetectorInfoVvv (G4int nDetectors,
			     G4double ekVertex, G4double xVertex, G4double yVertex, G4double zVertex, 
			     G4int idVolVertex, G4int idProcVertex, G4int idTrackVertex, G4int particleID) ;

    void SetOPSAinfo    (G4int nDetectors, G4int ID, G4int nPhot, G4int nPhotPrim, G4double timeFirst, G4double timeSecond, 
			 G4double timeThird, G4double timeA, G4double timeB, G4double timeC, G4double timeD,
			 G4double timeMean, G4double timeLast, G4double timeCFD, G4double amplCFD);

    void SetCFDSpecialInfo (G4int n, G4double time);

    void SetTimeC1SpecialInfo (G4double* time);

    void SetSaveDetectorInfo (G4int ID, G4int particleID, G4double ke, G4double x, G4double y, G4double z, G4double time, 
			      G4double px, G4double py, G4double pz, G4double polx, G4double poly, G4double polz) ;

    void SetInitialMuonParameters(G4double x, G4double y, G4double z, G4double px, G4double py, G4double pz, 
				  G4double xpolaris, G4double ypolaris, G4double zpolaris, G4double particleTime) {
      muIniTime=particleTime/CLHEP::microsecond;
      muIniPosX=x;  muIniPosY=y;  muIniPosZ=z;
      muIniMomX=px; muIniMomY=py; muIniMomZ=pz;
      muIniPolX=xpolaris; muIniPolY=ypolaris; muIniPolZ=zpolaris; 
    }
    void PrintInitialMuonParameters() {
      G4cout<<"musrRootOutput.hh: Initial muon parameters: x="<<muIniPosX<<", y="<<muIniPosY<<", z="<<muIniPosZ
	    <<", px="<<muIniMomX << ", py="<<muIniMomY<<", pz="<<muIniMomZ<<G4endl;
      G4cout<<"            polx="<<muIniPolX<<", poly="<<muIniPolY<<", polz="<<muIniPolZ<<G4endl;
      G4cout<<"            time at which muon was generated = "<<muIniTime<<G4endl;
      G4cout<<"            numberOfGeneratedEvents = "<<GeantParametersD[7]<<G4endl;
    }

    void SetPolInTarget(G4ThreeVector pol) {muTargetPolX=pol.x(); muTargetPolY=pol.y(); muTargetPolZ=pol.z();}
    void SetTimeInTarget(G4double time) {muTargetTime = time/CLHEP::microsecond;}
    void SetMomentumInTarget(G4ThreeVector mom) {muTargetMomX=(mom.x())/CLHEP::MeV; muTargetMomY=(mom.y())/CLHEP::MeV; muTargetMomZ=(mom.z())/CLHEP::MeV;}
    void SetPolInM0(G4ThreeVector pol) {muM0PolX=pol.x(); muM0PolY=pol.y(); muM0PolZ=pol.z();}
    void SetTimeInM0(G4double time) {muM0Time = time/CLHEP::microsecond;}
    void SetPolInM1(G4ThreeVector pol) {muM1PolX=pol.x(); muM1PolY=pol.y(); muM1PolZ=pol.z();}
    void SetTimeInM1(G4double time) {muM1Time = time/CLHEP::microsecond;}
    void SetPolInM2(G4ThreeVector pol) {muM2PolX=pol.x(); muM2PolY=pol.y(); muM2PolZ=pol.z();}
    void SetTimeInM2(G4double time) {muM2Time = time/CLHEP::microsecond;}
    void SetInitialPositronMomentum(G4ThreeVector mom) {posIniMomx=mom.x();  posIniMomy=mom.y(); posIniMomz=mom.z();}
    void SetNOptPhot(G4int value) {nOptPhot=value;}
    void SetPhotDetTime(G4double time);
    void SetDecayTime(G4double time) {muDecayTime=time/CLHEP::microsecond;}
    void SetNrFieldNomVal(G4int n) {nFieldNomVal = n;}
    void SetFieldNomVal(G4int i, G4double value);
    G4int GetNrOfVolumes() {return det_nMax;}
    void SetBFieldIntegral(G4double BxInt,G4double ByInt,G4double BzInt,G4double BzInt1,G4double BzInt2,G4double BzInt3) {
      BxIntegral=BxInt/CLHEP::m/CLHEP::tesla; ByIntegral=ByInt/CLHEP::m/CLHEP::tesla; BzIntegral=BzInt/CLHEP::m/CLHEP::tesla;
      BzIntegral1=BzInt1/CLHEP::m/CLHEP::tesla;BzIntegral2=BzInt2/CLHEP::mm;BzIntegral3=BzInt3/CLHEP::mm;
    }
    void StoreGeantParameter(Int_t i, Double_t value) {
      if (i<maxNGeantParameters) { GeantParametersD[i]=value; }
      else {G4cout<<"musrRootOutput.hh::StoreGeantParameter:  index="<<i<<" out of range"
		  <<" (maxNGeantParameters=" <<maxNGeantParameters<<")"<<G4endl;}
    };
    void setRootOutputDirectoryName(char dirName[1000]);

    TH2F *htest1, *htest2;
    TH1F *htest3, *htest4, *htest5, *htest6, *htest7, *htest8;

  public:
    static G4bool store_runID;
    static G4bool store_eventID;
    static G4bool store_weight;
    static G4bool store_timeToNextEvent;
    static G4bool store_BFieldAtDecay;
    static G4bool store_muIniTime;
    static G4bool store_muIniPosX;
    static G4bool store_muIniPosY;
    static G4bool store_muIniPosZ;
    static G4bool store_muIniMomX;
    static G4bool store_muIniMomY;
    static G4bool store_muIniMomZ;
    static G4bool store_muIniPolX;
    static G4bool store_muIniPolY;
    static G4bool store_muIniPolZ;
    static G4bool store_muDecayDetID;
    static G4bool store_muDecayPosX;
    static G4bool store_muDecayPosY;
    static G4bool store_muDecayPosZ;
    static G4bool store_muDecayTime;
    static G4bool store_muDecayPolX;
    static G4bool store_muDecayPolY;
    static G4bool store_muDecayPolZ;
    static G4bool store_muTargetTime;
    static G4bool store_muTargetPolX;
    static G4bool store_muTargetPolY;
    static G4bool store_muTargetPolZ;
    static G4bool store_muTargetMomX;
    static G4bool store_muTargetMomY;
    static G4bool store_muTargetMomZ;
    static G4bool store_muM0Time;
    static G4bool store_muM0PolX;
    static G4bool store_muM0PolY;
    static G4bool store_muM0PolZ;
    static G4bool store_muM1Time;
    static G4bool store_muM1PolX;
    static G4bool store_muM1PolY;
    static G4bool store_muM1PolZ;
    static G4bool store_muM2Time;
    static G4bool store_muM2PolX;
    static G4bool store_muM2PolY;
    static G4bool store_muM2PolZ;
    static G4bool store_posIniMomX;
    static G4bool store_posIniMomY;
    static G4bool store_posIniMomZ;
    static G4bool store_nOptPhot;
    static G4bool store_nOptPhotDet;
    static G4bool store_phot_time;
    static G4bool store_det_ID;
    static G4bool store_det_edep;
    static G4bool store_det_edep_el;
    static G4bool store_det_edep_pos;
    static G4bool store_det_edep_gam;
    static G4bool store_det_edep_mup;
    static G4bool store_det_nsteps;
    static G4bool store_det_length;
    static G4bool store_det_start;
    static G4bool store_det_end;
    static G4bool store_det_x;
    static G4bool store_det_y;
    static G4bool store_det_z;
    static G4bool store_det_kine;
    static G4bool store_det_VrtxKine;
    static G4bool store_det_VrtxX;
    static G4bool store_det_VrtxY;
    static G4bool store_det_VrtxZ;
    static G4bool store_det_VrtxVolID;
    static G4bool store_det_VrtxProcID;
    static G4bool store_det_VrtxTrackID;
    static G4bool store_det_VrtxParticleID;
    static G4bool store_det_VvvKine;
    static G4bool store_det_VvvX;
    static G4bool store_det_VvvY;
    static G4bool store_det_VvvZ;
    static G4bool store_det_VvvVolID;
    static G4bool store_det_VvvProcID;
    static G4bool store_det_VvvTrackID;
    static G4bool store_det_VvvParticleID;
    static G4bool store_fieldNomVal;
    static G4bool store_fieldIntegralBx;
    static G4bool store_fieldIntegralBy;
    static G4bool store_fieldIntegralBz;
    static G4bool store_fieldIntegralBz1;
    static G4bool store_fieldIntegralBz2;
    static G4bool store_fieldIntegralBz3;
    static G4bool store_odet_ID;
    static G4bool store_odet_nPhot;
    static G4bool store_odet_nPhotPrim;
    static G4bool store_odet_timeFirst;
    static G4bool store_odet_timeSecond;
    static G4bool store_odet_timeThird;
    static G4bool store_odet_timeA;
    static G4bool store_odet_timeB;
    static G4bool store_odet_timeC;
    static G4bool store_odet_timeD;
    static G4bool store_odet_timeMean;
    static G4bool store_odet_timeLast;
    static G4bool store_odet_timeCFD;
    static G4bool store_odet_amplCFD;
    static G4bool store_odet_timeCFDarray;
    static G4bool store_odet_timeC1;

    static G4int oldEventNumberInG4EqEMFieldWithSpinFunction;

  private:
    TFile* rootFile;
    TTree* rootTree;
    static musrRootOutput* pointerToRoot;
    static const Int_t maxNGeantParameters=30;
    char   rootOutputDirectoryName[1000];
    Double_t GeantParametersD[maxNGeantParameters];   // parameters transfered from GEANT to Root
          // 0 ... fieldOption:  0 ... no field, 1 ... uniform, 2 ... gaussian, 3 ... from table
          // 1 ... fieldValue:   intensity of the magnetic field         
          // 2 ... minimum of the generated decay time of the muon (in microsecond)
          // 3 ... maximum of the generated decay time of the muon (in microsecond)
          // 4 ... muon mean life time (in microsecond)
          // 5 ... nr. of the last generated event
          // 6 ... run number
          // 7 ... numberOfGeneratedEvents (i.e. number of the generated events; 
          //                    in case of Turtle nr. of events tried); 

  // Variables common to the whole event:
    Int_t runID;
    Int_t eventID;
    Double_t weight;
    Double_t timeToNextEvent;
    Double_t B_t[6];
    Double_t muIniTime;
    Double_t muIniPosX, muIniPosY, muIniPosZ;
    Double_t muIniMomX, muIniMomY, muIniMomZ;
    Double_t muIniPolX, muIniPolY, muIniPolZ;
    Int_t    muDecayDetID;
    Double_t muDecayPolX, muDecayPolY, muDecayPolZ;
    Double_t muTargetTime, muTargetPolX, muTargetPolY, muTargetPolZ;
    Double_t muTargetMomX, muTargetMomY, muTargetMomZ;
    Double_t muM0Time, muM0PolX, muM0PolY, muM0PolZ;
    Double_t muM1Time, muM1PolX, muM1PolY, muM1PolZ;
    Double_t muM2Time, muM2PolX, muM2PolY, muM2PolZ;
    Double_t muDecayPosX, muDecayPosY, muDecayPosZ;
    Double_t muDecayTime;
    Double_t posIniMomx, posIniMomy, posIniMomz;
    Int_t    nOptPhot, nOptPhotDet;
    static const Int_t maxNOptPhotDet=10000;
    Double_t phot_time[maxNOptPhotDet];

  public:
    static const Int_t maxNFieldnNominalValues=30;
  private:
    Int_t     nFieldNomVal;
    Double_t  fieldNomVal[maxNFieldnNominalValues];
    Double_t  BxIntegral, ByIntegral, BzIntegral;
    Double_t  BzIntegral1, BzIntegral2, BzIntegral3;

  // Variables for a particle in a given detector within the event
  public:
    static const Int_t maxNSubTracks=30;
  private:
  // Variables for the activity inside a given detector
  public:
    static const Int_t det_nMax=100;    // must be by 1 higher than the real number of detector "hits", because
                                                // else the detector nr. 0 is counted (0 is used if no
                                                // SensDetectorMapping correspond to a given logical volume).
  private:
    G4int     det_n;
    G4int     det_ID[det_nMax];
    G4double  det_edep[det_nMax];
    G4int     det_nsteps[det_nMax];
    G4double  det_length[det_nMax];
    G4double  det_edep_el[det_nMax];
    G4double  det_edep_pos[det_nMax];
    G4double  det_edep_gam[det_nMax];
    G4double  det_edep_mup[det_nMax];
    G4double  det_time_start[det_nMax];
    G4double  det_time_end[det_nMax];
    G4double  det_x[det_nMax];
    G4double  det_y[det_nMax];
    G4double  det_z[det_nMax];
    G4double  det_kine[det_nMax];
    G4double  det_VrtxKine[det_nMax];
    G4double  det_VrtxX[det_nMax];
    G4double  det_VrtxY[det_nMax];
    G4double  det_VrtxZ[det_nMax];
    G4int     det_VrtxVolID[det_nMax];
    G4int     det_VrtxProcID[det_nMax];
    G4int     det_VrtxTrackID[det_nMax];
    G4int     det_VrtxParticleID[det_nMax];
    G4double  det_VvvKine[det_nMax];
    G4double  det_VvvX[det_nMax];
    G4double  det_VvvY[det_nMax];
    G4double  det_VvvZ[det_nMax];
    G4int     det_VvvVolID[det_nMax];
    G4int     det_VvvProcID[det_nMax];
    G4int     det_VvvTrackID[det_nMax];
    G4int     det_VvvParticleID[det_nMax];

  public:
    static const Int_t odet_nMax=det_nMax;    
  
  private:
    G4int odet_n;
    G4int odet_ID[odet_nMax];
    G4int odet_nPhot[odet_nMax];
    G4int odet_nPhotPrim[odet_nMax];
    G4double odet_timeFirst[odet_nMax];
    G4double odet_timeSecond[odet_nMax];
    G4double odet_timeThird[odet_nMax];
    G4double odet_timeA[odet_nMax];
    G4double odet_timeB[odet_nMax];
    G4double odet_timeC[odet_nMax];
    G4double odet_timeD[odet_nMax];
    G4double odet_timeMean[odet_nMax];
    G4double odet_timeLast[odet_nMax];
    G4double odet_timeCFD[odet_nMax];
    G4double odet_amplCFD[odet_nMax];
    G4double odet_timeCFD100[odet_nMax];
    G4double odet_timeCFD101[odet_nMax];
    G4double odet_timeCFD102[odet_nMax];
    G4double odet_timeCFD103[odet_nMax];
    G4double odet_timeCFD104[odet_nMax];
    G4double odet_timeCFD105[odet_nMax];
    G4double odet_timeCFD106[odet_nMax];
    G4double odet_timeCFD107[odet_nMax];
    G4double odet_timeCFD108[odet_nMax];
    G4double odet_timeCFD109[odet_nMax];
    G4double odet_timeCFD110[odet_nMax];
    G4double odet_timeCFD111[odet_nMax];
    G4double odet_timeCFD112[odet_nMax];
    G4double odet_timeCFD200[odet_nMax];
    G4double odet_timeCFD201[odet_nMax];
    G4double odet_timeCFD202[odet_nMax];
    G4double odet_timeCFD203[odet_nMax];
    G4double odet_timeCFD204[odet_nMax];
    G4double odet_timeCFD205[odet_nMax];
    G4double odet_timeCFD206[odet_nMax];
    G4double odet_timeCFD207[odet_nMax];
    G4double odet_timeCFD208[odet_nMax];
    G4double odet_timeCFD209[odet_nMax];
    G4double odet_timeCFD210[odet_nMax];
    G4double odet_timeCFD211[odet_nMax];
    G4double odet_timeCFD212[odet_nMax];
    G4double odet_timeCFD300[odet_nMax];
    G4double odet_timeCFD301[odet_nMax];
    G4double odet_timeCFD302[odet_nMax];
    G4double odet_timeCFD303[odet_nMax];
    G4double odet_timeCFD304[odet_nMax];
    G4double odet_timeCFD305[odet_nMax];
    G4double odet_timeCFD306[odet_nMax];
    G4double odet_timeCFD307[odet_nMax];
    G4double odet_timeCFD308[odet_nMax];
    G4double odet_timeCFD309[odet_nMax];
    G4double odet_timeCFD310[odet_nMax];
    G4double odet_timeCFD311[odet_nMax];
    G4double odet_timeCFD312[odet_nMax];
    G4double odet_timeCFD400[odet_nMax];
    G4double odet_timeCFD401[odet_nMax];
    G4double odet_timeCFD402[odet_nMax];
    G4double odet_timeCFD403[odet_nMax];
    G4double odet_timeCFD404[odet_nMax];
    G4double odet_timeCFD405[odet_nMax];
    G4double odet_timeCFD406[odet_nMax];
    G4double odet_timeCFD407[odet_nMax];
    G4double odet_timeCFD408[odet_nMax];
    G4double odet_timeCFD409[odet_nMax];
    G4double odet_timeCFD410[odet_nMax];
    G4double odet_timeCFD411[odet_nMax];
    G4double odet_timeCFD412[odet_nMax];
    G4double odet_timeCFD500[odet_nMax];
    G4double odet_timeCFD501[odet_nMax];
    G4double odet_timeCFD502[odet_nMax];
    G4double odet_timeCFD503[odet_nMax];
    G4double odet_timeCFD504[odet_nMax];
    G4double odet_timeCFD505[odet_nMax];
    G4double odet_timeCFD506[odet_nMax];
    G4double odet_timeCFD507[odet_nMax];
    G4double odet_timeCFD508[odet_nMax];
    G4double odet_timeCFD509[odet_nMax];
    G4double odet_timeCFD510[odet_nMax];
    G4double odet_timeCFD511[odet_nMax];
    G4double odet_timeCFD512[odet_nMax];
    G4double odet_timeC1[odet_nMax];
    G4double odet_timeC2[odet_nMax];
    G4double odet_timeC3[odet_nMax];
    G4double odet_timeC4[odet_nMax];
    G4double odet_timeC5[odet_nMax];
    G4double odet_timeC6[odet_nMax];
    G4double odet_timeC7[odet_nMax];
    G4double odet_timeC8[odet_nMax];
    G4double odet_timeC9[odet_nMax];
    G4double odet_timeC10[odet_nMax];
    G4double odet_timeC11[odet_nMax];
    G4double odet_timeC12[odet_nMax];
    G4double odet_timeC13[odet_nMax];
    G4double odet_timeC14[odet_nMax];
    G4double odet_timeC15[odet_nMax];
    G4double odet_timeC16[odet_nMax];
    G4double odet_timeC17[odet_nMax];
    G4double odet_timeC18[odet_nMax];
    G4double odet_timeC19[odet_nMax];
    G4double odet_timeC20[odet_nMax];
    G4double odet_timeC21[odet_nMax];
    G4double odet_timeC22[odet_nMax];
    G4double odet_timeC23[odet_nMax];
    G4double odet_timeC24[odet_nMax];
    G4double odet_timeC25[odet_nMax];
    G4double odet_timeC26[odet_nMax];
    G4double odet_timeC27[odet_nMax];
    G4double odet_timeC28[odet_nMax];
    G4double odet_timeC29[odet_nMax];
    G4double odet_timeC30[odet_nMax];
    G4double odet_timeC31[odet_nMax];
    G4double odet_timeC32[odet_nMax];
    G4double odet_timeC33[odet_nMax];
    G4double odet_timeC34[odet_nMax];
    G4double odet_timeC35[odet_nMax];
    G4double odet_timeC36[odet_nMax];
    G4double odet_timeC37[odet_nMax];
    G4double odet_timeC38[odet_nMax];
    G4double odet_timeC39[odet_nMax];
    G4double odet_timeC40[odet_nMax];

  public:
    static const Int_t save_nMax=1000;

  private:
    G4int    save_n;
    G4int    save_detID[save_nMax];
    G4int    save_particleID[save_nMax];
    G4double save_ke[save_nMax];
    G4double save_time[save_nMax];
    G4double save_x[save_nMax];
    G4double save_y[save_nMax];
    G4double save_z[save_nMax];
    G4double save_px[save_nMax];
    G4double save_py[save_nMax];
    G4double save_pz[save_nMax];
    G4double save_polx[save_nMax];
    G4double save_poly[save_nMax];
    G4double save_polz[save_nMax];

    G4bool boolIsAnySpecialSaveVolumeDefined;

    std::map<std::string,int> SensDetectorMapping;
    std::map<std::string,int> ProcessIDMapping;
};

#endif

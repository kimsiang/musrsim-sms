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

#ifndef musrScintSD_h
#define musrScintSD_h 1

#include "G4VSensitiveDetector.hh"
#include "musrScintHit.hh"
#include "musrRootOutput.hh"

class G4Step;
class G4HCofThisEvent;

class APDidClass {
  public:
    APDidClass(G4int xAPDcellID, G4int xAPDcellNphot) {
      APDcellID    = xAPDcellID;
      APDcellNphot = xAPDcellNphot;
    }
    ~APDidClass() {}
    G4int GetAPDcellID() const {return APDcellID;}
    G4int GetAPDcellNphot() const {return APDcellNphot;}
    void  IncreaseAPDcellNphot() {APDcellNphot++;}

  private:    
    G4int APDcellID;
    G4int APDcellNphot;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class signalInfo {
  public:
    signalInfo(G4int id, G4int nP, G4int nPprim, G4double tFirst, G4double tSecond, G4double tThird, 
  	       G4double tA, G4double tB, G4double tC, G4double tD, G4double tE, G4double tLast,
               G4double tCFD, G4double aCFD, G4double tCFDarray[1000], G4double tOPSAC1array[41]) 
               {
                 detID=id; nPhot=nP; nPhotPrim=nPprim; timeFirst=tFirst; timeSecond=tSecond; timeThird=tThird; 
                 timeA=tA; timeB=tB; timeC=tC; timeD=tD; timeE=tE; timeLast=tLast;
                 timeCFD=tCFD; amplCFD=aCFD; 
                 if (musrRootOutput::store_odet_timeCFDarray) {for(int i=0;i<1000;i++) {timeCFDarray[i]=tCFDarray[i];}}
		 if (musrRootOutput::store_odet_timeC1)       {for(int i=1;i<41;i++)   {timeC1array[i] =tOPSAC1array[i];}}
               }
    ~signalInfo() {}
    void transferDataToRoot(musrRootOutput* myRootOut, G4int nn) {
               myRootOut->SetOPSAinfo(nn,detID,nPhot,nPhotPrim,timeFirst,timeSecond,timeThird,
				      timeA,timeB,timeC,timeD,timeE,timeLast,timeCFD,amplCFD);
	       if (musrRootOutput::store_odet_timeCFDarray) {
		 for (Int_t kk=0; kk<13; kk++) {
		   for (Int_t ll=0; ll<5; ll++) {
		     int index = (ll+1)*100+kk;
		     myRootOut -> SetCFDSpecialInfo(index,timeCFDarray[index]);
		   }
		 }
	       }

	       if (musrRootOutput::store_odet_timeC1) myRootOut -> SetTimeC1SpecialInfo(timeC1array);
    }

  private:
    G4int detID;
    G4int nPhot;
    G4int nPhotPrim;
    G4int nPhot_abs;
    G4int nPhot_refl;
    G4int nPhot_other;
    G4double timeFirst;
    G4double timeSecond;
    G4double timeThird;
    G4double timeA;
    G4double timeB;
    G4double timeC;
    G4double timeD;
    G4double timeE;
    G4double timeLast;
    G4double timeCFD;
    G4double amplCFD;
    G4double timeCFDarray[1000];
    G4double timeC1array[1000];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class musrScintSD : public G4VSensitiveDetector
{
  public:
    static musrScintSD* GetInstance();
    musrScintSD(G4String);
    ~musrScintSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);
    void EndOfRun();

    // Optical Photon Signal Analysis (OPSA)
    void Set_OPSA_minNrOfDetectedPhotons(G4int val) {OPSA_minNrOfDetectedPhotons=val;}
    void Set_OPSA_SignalSeparationTime(G4double val) {OPSA_signalSeparationTime=val;}
    void Set_OPSA_variables(G4double a, G4double b, G4double c, G4double d) 
             {OPSA_fracA=a; OPSA_fracB=b; OPSA_C_threshold=c; OPSA_D_threshold=d;}
    void Set_OPSA_CFD(G4double a1, G4double delay, G4double timeShiftOffset) 
        {OPSA_CFD_a1=a1; OPSA_CFD_delay=delay; OPSA_CFD_timeShiftOffset = timeShiftOffset;}
    void AddEventIDToMultimapOfEventIDsForOPSAhistos (G4int ev_ID, G4int detector_ID) {
         bool_multimapOfEventIDsForOPSAhistosEXISTS=true;
	 if (ev_ID==-1)  bool_StoreThisOPSAhistSUMMED = true;
	 if (ev_ID==-2)  bool_StoreThisOPSAhistALL = true;
	 multimapOfEventIDsForOPSAhistos.insert(std::pair<G4int,G4int>(ev_ID,detector_ID));
    }
    void SetOPSAhistoBinning(Int_t nBins, Double_t min, Double_t max) {
      OPSAhistoNbin=nBins; OPSAhistoMin=min; OPSAhistoMax=max;
      OPSAhistoBinWidth=(max-min)/nBins; OPSAhistoBinWidth1000=OPSAhistoBinWidth*1000;
    }
    void ProcessOpticalPhoton(G4Step* aStep);
    void EndOfEvent_OptiacalPhotons();
    void ReadInPulseShapeArray(const char* filename);
    void FindCFDtime(G4double& OPSA_CFD_time, G4double& OPSA_CFD_ampl, G4double timeOfFirstPhoton);
    G4int FindAPDcellID(G4Step* aStep);
    void SetAPDcellSizes(G4int nx, G4int ny, G4int nz, G4double half_ax, G4double half_ay, G4double half_az) {
      APDcellsEffectRequested = true;
      APDcell_nx = nx; APDcell_ny = ny; APDcell_nz = nz; 
      APDcell_ax = 2*half_ax/ nx; APDcell_ay = 2*half_ay/ ny; APDcell_az = 2*half_az/ nz; 
    }
    void SetAPDcellsTimeVariationSigma(G4double sigma) {
      if (sigma!=0) APDcellsTimeVariationRequested = true;
      APDcellsTimeVariationSigma=sigma;
    }
    void SetAPDcrossTalk(G4double crosstalkProb) {
      if (crosstalkProb>0) APDcrossTalkRequested = true;
      APDcrossTalkProb = crosstalkProb;
    }

  private:
    static musrScintSD* pointer;
    musrScintHitsCollection* scintCollection;
    G4bool   myStoreOnlyEventsWithHits;
    G4int    myStoreOnlyEventsWithHitInDetID;
    G4double mySignalSeparationTime;
    G4bool   myStoreOnlyTheFirstTimeHit;
    G4int    myStoreOnlyEventsWithMuonsDecayedInDetID;
    G4bool   boolIsVvvInfoRequested;
    musrRootOutput* myRootOutput;

    // for optical photon counting
  //    typedef std::multimap<G4double,Int_t> optHitDetectorMapType;
    typedef std::multimap<G4double,APDidClass> optHitDetectorMapType;
    typedef std::map<Int_t,optHitDetectorMapType*> optHitMapType;
    optHitMapType optHitMap;
    // for optical photon signal analysis (OPSA)
    G4int    OPSA_minNrOfDetectedPhotons;
    G4double OPSA_signalSeparationTime;
    G4double OPSA_fracA;
    G4double OPSA_fracB;
    G4double OPSA_C_threshold;
    G4double OPSA_D_threshold;
    typedef std::multimap<G4int,signalInfo*>  OPSA_signal_MapType;
    OPSA_signal_MapType OPSA_signal_Map;
  
    G4bool bool_multimapOfEventIDsForOPSAhistosEXISTS;
    G4bool bool_StoreThisOPSAhistSUMMED;
    G4bool bool_StoreThisOPSAhistALL;
    typedef std::multimap<G4int,G4int> multimapOfEventIDsForOPSAhistos_Type;
    multimapOfEventIDsForOPSAhistos_Type multimapOfEventIDsForOPSAhistos;
    TH1D* OPSAhisto;
    Int_t OPSAhistoNbin;
    Double_t OPSAhistoMin, OPSAhistoMax, OPSAhistoBinWidth, OPSAhistoBinWidth1000;
    typedef std::map<std::string,TH1D*>  mapOfOPSAsumHistograms_Type;
    mapOfOPSAsumHistograms_Type mapOfOPSAsumHistograms;
    mapOfOPSAsumHistograms_Type mapOfOPSAsum0Histograms;
    TH1D* OPSAhistoSUM;
    TH1D* OPSAhistoSUM0;
    TH1D* OPSAshape;
    G4bool bool_pulseShapeExists;
    G4int iPSmax;
    G4double pulseShapeArray[10001];
    TH1D* OPSA_CFD;
    G4double OPSA_CFD_a1,OPSA_CFD_delay, OPSA_CFD_timeShiftOffset;
    G4bool   APDcellsEffectRequested;              // simulate effects of finite cell number
    G4int    APDcell_nx, APDcell_ny, APDcell_nz;   // number of cells in APD along x, y and z directions
    G4double APDcell_ax, APDcell_ay, APDcell_az;   // dimensions of APD cells
    G4bool   APDcellsTimeVariationRequested;       // simmulate effect of detection time variations between differenct cells
    G4double APDcellsTimeVariationSigma;           // sigma of the detection time variations between differenct cells
    G4bool   APDcrossTalkRequested;
    G4double APDcrossTalkProb;

    void FireAPDcell(optHitDetectorMapType* optHitDetectorMap, G4int APDcellID, G4double time, G4int nTruePhe);

    static const G4double OPSA_C1_threshold;
    static const G4double OPSA_C2_threshold;
    static const G4double OPSA_C3_threshold;
    static const G4double OPSA_C4_threshold;
    static const G4double OPSA_C5_threshold;
    static const G4double OPSA_C6_threshold;
    static const G4double OPSA_C7_threshold;
    static const G4double OPSA_C8_threshold;
    static const G4double OPSA_C9_threshold;
    static const G4double OPSA_C10_threshold;
    static const G4double OPSA_C11_threshold;
    static const G4double OPSA_C12_threshold;
    static const G4double OPSA_C13_threshold;
    static const G4double OPSA_C14_threshold;
    static const G4double OPSA_C15_threshold;
    static const G4double OPSA_C16_threshold;
    static const G4double OPSA_C17_threshold;
    static const G4double OPSA_C18_threshold;
    static const G4double OPSA_C19_threshold;
    static const G4double OPSA_C20_threshold;
    static const G4double OPSA_C21_threshold;
    static const G4double OPSA_C22_threshold;
    static const G4double OPSA_C23_threshold;
    static const G4double OPSA_C24_threshold;
    static const G4double OPSA_C25_threshold;
    static const G4double OPSA_C26_threshold;
    static const G4double OPSA_C27_threshold;
    static const G4double OPSA_C28_threshold;
    static const G4double OPSA_C29_threshold;
    static const G4double OPSA_C30_threshold;
    static const G4double OPSA_C31_threshold;
    static const G4double OPSA_C32_threshold;
    static const G4double OPSA_C33_threshold;
    static const G4double OPSA_C34_threshold;
    static const G4double OPSA_C35_threshold;
    static const G4double OPSA_C36_threshold;
    static const G4double OPSA_C37_threshold;
    static const G4double OPSA_C38_threshold;
    static const G4double OPSA_C39_threshold;
    static const G4double OPSA_C40_threshold;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


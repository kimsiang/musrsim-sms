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

#ifndef musrPrimaryGeneratorAction_h
#define musrPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include <stdio.h>
#include <vector>

class G4GeneralParticleSource;
class G4ParticleGun;
class G4Event;
class musrDetectorConstruction;
class musrPrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class musrPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    musrPrimaryGeneratorAction(musrDetectorConstruction*);    
   ~musrPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);
    void SetRndmFlag(G4String val)          { rndmFlag = val;}
    void Setvertex(G4ThreeVector v)         {x0=v[0]; y0=v[1]; z0=v[2];}
    void SetvertexSigma(G4ThreeVector v)    {xSigma=v[0]; ySigma=v[1]; zSigma=v[2];}
    void SetvertexBoundary(G4ThreeVector v) {rMaxAllowed=v[0]; zMinAllowed=v[1]; zMaxAllowed=v[2];}
    void SetvertexRelativeR(G4double val)   {relativeRMaxAllowed=val;}
    void SetboxBoundary(G4ThreeVector v)    {xMaxSource=v[0]; yMaxSource=v[1]; zMaxSource=v[2];}          //P.B. 15 Dec 2009
    void SetboxBoundaryCentre(G4ThreeVector v) {xMaxSource0=v[0]; yMaxSource0=v[1]; zMaxSource0=v[2];}    //P.B. 15 Dec 2009
    void SetMeanArrivalTime(G4double val) {meanArrivalTime=val;}
    void SetMuonTime(G4double val)          {t0=val;}     //P.B. 13 May 2009
    void SetMuonTimeSigma(G4double val)     {tSigma=val;} //P.B. 13 May 2009
    void SetKEnergy(G4double val);
    void SetMomentum(G4double val)          {p0=val;}
    void SetMomentumSmearing(G4double val)  {pSigma=val;}
    void SetMomentumBoundary(G4ThreeVector v){pMinAllowed=v[0]; pMaxAllowed=v[1];}
    void SetTilt(G4ThreeVector v)           {xangle0=v[0]; yangle0=v[1];}
    void SetSigmaTilt(G4ThreeVector v)      {xangleSigma=v[0]; yangleSigma=v[1];zangleSigma=v[2];}
    void SetPitch(G4double val)             {pitch=val;}
    void SetBeamDirection(G4ThreeVector vIniDir);
    void SetInitialMuonPolariz(G4ThreeVector vIniPol);
    void SetInitialPolarizFraction(G4double val) {
      if ((val>1.)||(val<-1.)) {
	G4cout<<"musrPrimaryGeneratorAction.hh: SetInitialPolarizFraction(): polarisation fraction out of range ("<<val<<")"<<G4endl; 
	exit(1);
      }
      polarisFraction=val;
    }
    void SetMuonDecayTimeLimits(G4ThreeVector decayTimeLimits);
    void SetTurtleInput(G4String turtleFileName);
    void SetTurtleInputFileToEventNo(G4int lineNumberOfTurtleFile);
    void SetTurtleZ0(G4double val)          {z0_InitialTurtle=val;}
    void SetTurtleInterpretAxes(G4String interpretAxes){turtleInterpretAxes=interpretAxes;}
  //    void SetOrReadTheRandomNumberSeeds(G4int eventID);
    void SetOrReadTheRandomNumberSeeds(G4Event* anEvent);
    void SetTurtleMomentumBite (G4ThreeVector smearingParam)
                               {turtleMomentumBite=true; turtleMomentumP0=smearingParam[0]*CLHEP::MeV; turtleSmearingFactor=smearingParam[1]*0.01;}
    void SetTurtleMomentumScalingFactor(G4double momentumScaling) {turtleMomentumScalingFactor=momentumScaling;}
    void SetPrimaryParticule(G4String particleName);
    static G4String GetPrimaryName();                

  private:
    G4ParticleGun*                particleGun;	  // pointer a to G4 service class
    G4GeneralParticleSource*      particleSource; // pointer to the G4GeneralParticleSource, needed for radioactive samples
    musrDetectorConstruction*      musrDetector;  // pointer to the geometry
    
    musrPrimaryGeneratorMessenger* gunMessenger;  // messenger of this class
    G4String                      rndmFlag;	  // flag for a random impact point  
    G4String                      turtleInterpretAxes;  // specifies how to intrpret the TURTLE position axes x, y and z

  // cks delete    G4ParticleDefinition* muonMinusParticle;
    // cks Alpha and proton particles implemented for the simulation of Juan Pablo Urrego
    G4ParticleDefinition* alphaParticle;
    G4ParticleDefinition* protonParticle;
    // csk

    static G4String thePrimaryParticleName ;

    G4double x0, y0, z0, xSigma, ySigma, zSigma, rMaxAllowed, zMinAllowed, zMaxAllowed;
    G4double meanArrivalTime;
    G4double t0, tSigma;                                 //P.B. 13 May 2009
    G4double relativeRMaxAllowed;
    G4double xMaxSource0, yMaxSource0, zMaxSource0;      //P.B. 15 Dec 2009
    G4double xMaxSource, yMaxSource, zMaxSource;         //P.B. 15 Dec 2009
    G4double p0, pSigma, pMinAllowed, pMaxAllowed;
    G4double xangle0, yangle0, xangleSigma, yangleSigma, zangleSigma,pitch;
    G4bool   UnpolarisedMuonBeam, TransversalyUnpolarisedMuonBeam;
    G4double xPolarisIni, yPolarisIni, zPolarisIni;
    G4double xDirection, yDirection, zDirection;
    G4double polarisFraction;
    G4double muonDecayTimeMin;
    G4double muonDecayTimeMax;
    G4double muonMeanLife;
  // For the Turtle input:
    FILE*    fTurtleFile;
    G4bool   takeMuonsFromTurtleFile;
    G4double z0_InitialTurtle;                   // z0 at whith the turtle file was generated.
    G4int    numberOfGeneratedEvents;            // number of generated events at the input (e.g. including Turtle events out of the beampipe)
    G4bool   boolPrintInfoAboutGeneratedParticles;
    G4bool   turtleMomentumBite;
    G4double turtleMomentumP0;
    G4double turtleSmearingFactor;
    G4double turtleMomentumScalingFactor;
    void     swapTheAxisInTurtle(float& x_x, float& x_xprime, float& y_y, float& y_yprime);

public: 
    static G4bool setRandomNrSeedAccordingEventNr;
    static G4bool setRandomNrSeedFromFile;
    static G4bool setRandomNrSeedFromFile_RNDM;
    static G4int  nRndmEventToSaveSeeds;
    static std::vector<int> * GetPointerToSeedVector();
    static G4int lastEventID_in_pointerToSeedVector;
    G4double  decaytime;

private:
  static std::vector<int> * pointerToSeedVector;

};

#endif



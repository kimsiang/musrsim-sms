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

#ifndef musrScintHit_h
#define musrScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"
#include "globals.hh"
#include "G4ios.hh"
//  ROOT
#include "TFile.h"
#include "TTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class musrScintHit : public G4VHit
{
  public:

      musrScintHit();
     ~musrScintHit();
      musrScintHit(const musrScintHit&);
      const musrScintHit& operator=(const musrScintHit&);
      G4int operator==(const musrScintHit&) const;
  //      bool operator() (musrScintHit hit1, musrScintHit hit2) { return (hit1.globalTime<hit2.globalTime);}
  
      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:

      void SetParticleName (G4String name) {particleName = name; }
      void SetParticleID (G4int id) {particleID = id; }
      void SetTrackID  (G4int track)      { trackID = track; }
      void SetEdep     (G4double de)      { edep = de; }
      void SetPrePos      (G4ThreeVector xyz){ pre_pos = xyz; }
      void SetPostPos      (G4ThreeVector xyz){ post_pos = xyz; }
      void SetPol      (G4ThreeVector ijk){pol = ijk;}
      void SetLogVolName   (G4String logivol) {logicalVolume = logivol;}
      void SetGlobTime     (G4double gt)      { globalTime = gt;}
      void SetFirstStepInVolumeFlag (G4bool flag) { firstStepInVolume=flag;}
      void SetLastStepInVolumeFlag (G4bool flag)  { lastStepInVolume=flag;}
      void SetKineticEnergy     (G4double en)     { kineticEnergy = en;}
      void SetStepLength     (G4double length)    { stepLength = length;}
      void SetRunID(G4int i)                    {runID=i;}
      void SetEventID(G4int i)                  {eventID=i;}
      void SetVertexPosition(G4ThreeVector xyz) {vertexPosition = xyz; }
      void SetVertexKineticEnergy(G4double Ek)  {vertexKineticEnergy = Ek; }
      void SetLogicalVolumeAtVertex(G4String logivol) {logicalVolumeAtVertex = logivol; }
      void SetCreatorProcessName(G4String name) {creatorProcess = name; }
  //      void SetVerboseLevel   (G4int n) {  G4int musrScintHit::verboseLevel=n;};

      G4String GetParticleName() {return particleName; }
      G4int GetParticleID() {return particleID; }
      G4int GetTrackID()    { return trackID; }
      G4double GetEdep()    { return edep; }    
      G4ThreeVector GetPrePos(){ return pre_pos; }
      G4ThreeVector GetPostPos(){ return post_pos; }
      G4ThreeVector GetPol(){ return pol; }
      G4String GetLogVolName() { return logicalVolume; }
      G4double GetGlobalTime() { return globalTime; }
      G4bool GetFirstStepInVolumeFlag() {return firstStepInVolume;}
      G4bool GetLastStepInVolumeFlag() {return lastStepInVolume;}
      G4double GetKineticEnergy() { return kineticEnergy; }
      G4double GetStepLength() {return stepLength; }
      G4double* GetBField () {return BF;}
      G4int GetRunID() {return runID;}
      G4int GetEventID() {return eventID;}
      G4ThreeVector GetVertexPosition()   { return vertexPosition; }
      G4double GetVertexKineticEnergy()   { return vertexKineticEnergy; }
      G4String GetLogicalVolumeAtVertex() { return logicalVolumeAtVertex; }
      G4String GetCreatorProcessName()         { return creatorProcess; }
      G4double point[4];
      G4double B[6];
      const G4Field *mfield;

  private:

      G4String      particleName;
      G4int         particleID;
      G4int         trackID;
      G4double      edep;
      G4ThreeVector pre_pos;
      G4ThreeVector post_pos;
      G4ThreeVector pol;
      G4String      logicalVolume;
      G4double      globalTime;
      G4bool        firstStepInVolume;
      G4bool        lastStepInVolume;
      G4double      kineticEnergy;
      G4double      stepLength;
      G4int         eventID;
      G4int         runID;
      G4double      BF[6];
      G4ThreeVector vertexPosition;  // Position of the vertex at which the actual track was created
      G4double      vertexKineticEnergy;  // Kinetic energy of the track when it was created 
      G4String      logicalVolumeAtVertex;
      G4String      creatorProcess;

      static G4int ScintMultihit;
      static G4int runIDoldScint;
      static G4int eventIDoldScint;
      static G4int NIS;
      static G4int ScintChamberNbold;

      static G4int verboseLevel;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<musrScintHit> musrScintHitsCollection;

extern G4Allocator<musrScintHit> musrScintHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* musrScintHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) musrScintHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void musrScintHit::operator delete(void *aHit)
{
  musrScintHitAllocator.FreeSingle((musrScintHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

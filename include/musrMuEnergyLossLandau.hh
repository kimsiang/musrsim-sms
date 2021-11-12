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

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//  Muon energy loss in thin C-foil adding a Landau distribution to the energy 
//  loss
//  Id    : musrMuEnergyLossLandau, v 1.0
//  Author: Thomas Prokscha
//  Date  : 2016-08
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#ifndef   musrMuEnergyLossLandau_h
#define   musrMuEnergyLossLandau_h 1

#include "G4VDiscreteProcess.hh"
#include "G4ParticleTable.hh"

#include <TMath.h>
#include <TRandom.h>


/*! musrMuEnergyLossLandau class defines the energy loss given by the Landau distribution,
   in, for example, the thin Carbon foil of the LEM beam line 
  */

class musrMuEnergyLossLandau : public G4VDiscreteProcess
{
 public:
   
   musrMuEnergyLossLandau(const G4String& name = "MuEnergyLossLandau", // process description
		   G4ProcessType aType = fElectromagnetic);

  ~musrMuEnergyLossLandau();
  
  //! - Main method. muon energy loss process is executed at the END of a step. */
  G4VParticleChange* PostStepDoIt(
			     const G4Track&,
			     const G4Step&);
  
  G4double GetMeanFreePath(const G4Track& aTrack,
			   G4double previousStepSize,
			   G4ForceCondition* condition);


  //! Condition for process application (step Object).
  G4bool CheckCondition(const G4Step& aStep);
  
  //! Condition for process application (step Pointer).
  G4bool CheckCondition(const G4Step* aStep);
  
  
  G4String  p_name;
  G4bool condition;
  
  
  void GetFinalEnergy( const G4Step* aStep);
  TRandom *random;
  static double landauMPV;
  static double landauSigma;

  // model parameters
  G4ParticleTable* particleTable; 
  G4ParticleDefinition* particle;
  G4double rnd;
  G4DynamicParticle *DP;
  
  //! The particle change object.
  G4VParticleChange fParticleChange; 
 
  void  PrepareSecondary(const G4Track&);
  G4Track* aSecondary;

  void InitializeSecondaries(const G4Track&);
};

#endif

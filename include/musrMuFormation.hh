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
//  Muonium Formation according to yield.cc function (through GetYields method).
//  Id    : musrMuFormation.hh, v 1.4
//  Author: Taofiq PARAISO, T. Shiroka
//  Date  : 2007-12
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                              
#ifndef   musrMuFormation_h
#define   musrMuFormation_h 1

#include "G4VDiscreteProcess.hh"
#include "G4ParticleTable.hh"

#include "yields.hh"

/*! musrMuFormation class defines the muonium formation process in the Carbon foil 
 *  according to yields from Gonin's paper Sci. Rev. Instrum. 65(3), 648-652 (1994).
 * \image html yields3.gif The muonium formation yields.
 * The main parameters are the foil thickness and muon energy. For a given energy, 
 * a corresponding proportion of the muons will be converted into Muonium.
 * Concretely, the muon is eliminated and replaced by a Muonium with identical 
 * properties, including time, energy, momentum, position etc.
 *
 * The process is executed at the END of a step, i.e. the muon is converted into 
 * Muonium AFTER flying through the Carbon foil (see also yields.hh). */

class musrMuFormation : public G4VDiscreteProcess
{
 public:
   
   musrMuFormation(const G4String& name = "MuFormation", // process description
		   G4ProcessType aType = fElectromagnetic);

  ~musrMuFormation();
  
  //! - Main method. Muonium formation process is executed at the END of a step. */
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
  
  
  void GetDatas( const G4Step* aStep);
  // model parameters
  G4ParticleTable* particleTable; 
  G4ParticleDefinition* particle;
  Yields   Gonin;
  G4double yvector[3];
  G4double rnd;
  G4DynamicParticle *DP;
  
  //! The particle change object.
  G4VParticleChange fParticleChange; 
 
  void  PrepareSecondary(const G4Track&);
  G4Track* aSecondary;

  void InitializeSecondaries(const G4Track&);
};

#endif

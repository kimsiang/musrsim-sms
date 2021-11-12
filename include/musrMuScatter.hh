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

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//  Muonium "Scattering"
//  Id    : musrMuScatter.hh, v 1.4
//  Author: Taofiq PARAISO, T. Shiroka
//  Date  : 2007-12
//  Notes : Simplified model for Mu scattering. Spin effects have been excluded.
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

/*!
 * musrMuScatter class defines the Muonium scattering process. It implements a very 
 * basic model which assumes Muonium looses its electron as soon as it enters any 
 * material (except for vacuum and CFoil). The class has only a PostStepDoIt method.
 * The in-flight Muonium spin precession has been supressed. */

#ifndef musrMuScatter_h
#define musrMuScatter_h 1

#include "G4VDiscreteProcess.hh"

class musrMuScatter : public G4VDiscreteProcess
{
 public:
   musrMuScatter(const G4String& name="MuScatt", // process description
		G4ProcessType aType = fGeneral);
  
  ~musrMuScatter();

  //! \mm The actions are taken at the end of the step.
  G4VParticleChange* PostStepDoIt(const G4Track&,
  				  const G4Step&);

  G4double GetMeanFreePath(const G4Track& aTrack,
			   G4double previousStepSize,
			   G4ForceCondition* condition);
  
  //! The condition for applying the process.
  G4bool CheckCondition(const G4Step& aStep);

  
  G4bool   condition;
  G4double itime, gtime, ftime,deltatime;
  G4String p_name;
  G4DynamicParticle *DP;
  G4ParticleChange  fParticleChange; 
 
  void  PrepareSecondary(const G4Track&);
  G4Track* aSecondary;

  void InitializeSecondaries(const G4Track&);
};

#endif

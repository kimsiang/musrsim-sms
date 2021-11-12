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
//  Id    : musrMuScatter.cc, v 1.4
//  Author: Taofiq PARAISO, T. Shiroka
//  Date  : 2007-12
//  Notes : Simplified model for Mu scattering. Spin effects have been excluded.
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#include "musrMuScatter.hh"

using namespace std;

musrMuScatter::musrMuScatter(const G4String& name,
			     G4ProcessType  aType)
               : G4VDiscreteProcess(name, aType){}

musrMuScatter:: ~musrMuScatter(){}

// At the end of the step, the current volume is checked and if Muonium is in a solid 
// material (except for the carbon foil where it is generated), it is stopped immediately.
G4VParticleChange* musrMuScatter::PostStepDoIt(const G4Track& trackData,
                                               const G4Step& aStep)
{
  fParticleChange.Initialize(trackData);
  
  //! Tao - Get time information */
  itime = trackData.GetProperTime();
  gtime = trackData.GetGlobalTime();
  ftime = trackData.GetDynamicParticle()->GetPreAssignedDecayProperTime(); 
  
  deltatime = ftime - itime;
  fParticleChange.ProposeGlobalTime(deltatime + itime -gtime);
  
  /*! - Set position, momentum, energy and time of the particle change. */
  fParticleChange.ProposePosition(trackData.GetPosition());
  fParticleChange.ProposeMomentumDirection(trackData.GetMomentumDirection());
  fParticleChange.ProposeEnergy(trackData.GetKineticEnergy());
  fParticleChange.ProposeGlobalTime(gtime);
  fParticleChange.ProposeProperTime(itime);
  fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;
  
  /*! - Verify the condition of applying the process: if Mu is in a material 
        different than vacuum and carbon foil, then stop it directly. */
  if( CheckCondition(aStep))
    {
      fParticleChange.ProposePosition(trackData.GetStep()->GetPreStepPoint()->GetPosition());
      fParticleChange.ProposeTrackStatus(fStopButAlive) ;
    }
  
  /*! - Return the changed particle object. */
  return &fParticleChange;
}


/*! - Muonium will be stopped as soon as it enters a material different than vacuum or C foil. */
G4bool musrMuScatter::CheckCondition(const G4Step& aStep)
{
  G4bool condition = false;
  p_name = aStep.GetTrack()->GetDefinition()->GetParticleName(); // particle name  
  if(p_name == "Mu" && aStep.GetTrack()->GetVolume()->GetLogicalVolume()->GetName()!="log_CFoil" &&
  aStep.GetTrack()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName()!="G4_Galactic")
    {
      condition=true;
    }
  return condition;
}  


G4double musrMuScatter::GetMeanFreePath(const G4Track&,
					G4double,
					G4ForceCondition* condition)
{
  *condition = Forced;
   return DBL_MAX;
}


void musrMuScatter::PrepareSecondary(const G4Track& track)
{
  aSecondary = new G4Track(DP,track.GetDynamicParticle()->GetPreAssignedDecayProperTime(),track.GetPosition());
}

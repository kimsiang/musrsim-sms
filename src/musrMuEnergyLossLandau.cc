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

#include "musrMuEnergyLossLandau.hh"

using namespace std;

musrMuEnergyLossLandau::musrMuEnergyLossLandau(const G4String& name, G4ProcessType aType)
               : G4VDiscreteProcess(name, aType)
{
  random = new TRandom();
}

musrMuEnergyLossLandau::~musrMuEnergyLossLandau(){}

double musrMuEnergyLossLandau::landauMPV   = 0.1;
double musrMuEnergyLossLandau::landauSigma = 0.1;


G4VParticleChange* musrMuEnergyLossLandau::PostStepDoIt(const G4Track& trackData,
                                                 const G4Step&  aStep)
{ // Initialize ParticleChange  (by setting all its members equal to
  //                             the corresponding members in G4Track)
  fParticleChange.Initialize(trackData);

  G4Track theNewTrack;
  if(CheckCondition(aStep))
    {
      GetFinalEnergy(&aStep);
      G4Step theStep;
      PrepareSecondary( trackData);      
      fParticleChange.AddSecondary(aSecondary);
      fParticleChange.ProposeTrackStatus(fStopAndKill) ;
    }  
  else
    {
      fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;
    }
  return &fParticleChange;
}


G4bool musrMuEnergyLossLandau::CheckCondition(const G4Step& aStep)
{ // Decide when to call the MuEnergyLossLandau - i.e. for muons going through the C foil.
  condition=false;
  p_name = aStep.GetTrack()->GetDefinition()->GetParticleName(); // particle name  
  std::string logVolName = aStep.GetTrack()->GetVolume()->GetLogicalVolume()->GetName();
  if(p_name == "mu+" && ((logVolName=="log_coulombCFoil")||(logVolName=="log_CFoil")))
    {
      condition=true;
    }
  return condition;
}


G4double musrMuEnergyLossLandau::GetMeanFreePath(const G4Track&, 
                                                G4double,
                                                G4ForceCondition* condition1)
{
  *condition1 = Forced;
   return DBL_MAX;
}


void musrMuEnergyLossLandau::GetFinalEnergy(const G4Step* aStep)
{    // Particle generation
     particleTable=G4ParticleTable::GetParticleTable();     
     G4double Eloss, Efinal;
     G4double E = aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy()/CLHEP::keV;
     
     // Positive muon
     if(p_name=="mu+"){
      do{
        Eloss = random->Landau(landauMPV,landauSigma);
	if (Eloss > 0.) break;
      } while (1);     
     
     Efinal = E - Eloss; 
     if (Efinal < 0. ) Efinal = 0.; 
//      cout << "E, Eloss, Efinal = " << E << ", " << Eloss << ", "  << Efinal << endl;
 	     
     particle = particleTable->FindParticle(p_name) ;

     // Set the new dynamic particle DP
     DP = new G4DynamicParticle(particle, 
	         aStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection(),
	         Efinal/1000.);
	 
/*    IMPORTANT : COPY THOSE DATA TO GET THE SAME PARTICLE PROPERTIES!!!
      SHOULD BE KEPT WHEN BUILDING A PARTICLE CHANGE  */
      DP->SetProperTime(  aStep->GetTrack()->GetDynamicParticle()->GetProperTime());
      DP->SetPolarization(aStep->GetTrack()->GetDynamicParticle()->GetPolarization().x(),
	                     aStep->GetTrack()->GetDynamicParticle()->GetPolarization().y(),
			     aStep->GetTrack()->GetDynamicParticle()->GetPolarization().z());
      DP->SetPreAssignedDecayProperTime(
	 aStep->GetTrack()->GetDynamicParticle()->GetPreAssignedDecayProperTime());
     }
}


void musrMuEnergyLossLandau::PrepareSecondary(const G4Track& track)
{
  if(p_name=="mu+")
    {
     aSecondary = new G4Track(DP,track.GetGlobalTime(),track.GetPosition());
    }
}

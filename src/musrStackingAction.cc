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

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "musrStackingAction.hh"
#include "musrParameters.hh"
#include "musrErrorMessage.hh"
#include "musrRootOutput.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrStackingAction::musrStackingAction()
: opticalPhotonCounter(0)
{
  //std::cout<<" DEBUG KAMIL: StackingAction() called"<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

musrStackingAction::~musrStackingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack musrStackingAction::ClassifyNewTrack(const G4Track * aTrack) {
  if (aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
    opticalPhotonCounter++;
    return fWaiting;          // Optical photons will be treated after all other particles in order
                              // to have the reproducibility of the other particle tracks if a parameter
                              // for the optical photons is changed in the next run.
  }
  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrStackingAction::NewStage() {  
  //    G4cout << "Number of optical photons produces in this event : "
  //           << opticalPhotonCounter << G4endl;
  musrRootOutput::GetRootInstance()->SetNOptPhot(opticalPhotonCounter);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void musrStackingAction::PrepareNewEvent()
{ 
  opticalPhotonCounter = 0; 
}

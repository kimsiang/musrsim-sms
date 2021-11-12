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

#ifndef musrParameters_h
#define musrParameters_h 1

#include "globals.hh"

class musrParameters {
  public:
    musrParameters(G4String steeringFileName);
    ~musrParameters();

    static musrParameters* GetInstance();

    static G4String mySteeringFileName;      // name of the steering file (e.g. the *.mac file)
    static G4String myStopFileName;          // name of the stop file (e.g. the *.stop file), which will stop the run if it is created
    static G4String myRandomNumberFileName;  // name of the file with random numbers for RandomOption=3 in musrDetectorMessenger
    static G4bool storeOnlyEventsWithHits;   // variable specifying whether to store interesting 
                                             // or all events into the ROOT tree. (default = true)
    static G4int  storeOnlyEventsWithHitInDetID; // simillar to "storeOnlyEventsWithHits".  The event is stored
                                             // only and only if there was a hit in the detector with the ID 
                                             // equal to storeOnlyEventsWithHitInDetID.
                                             // The storeOnlyEventsWithHitInDetID has to be non zero.
    static G4double signalSeparationTime;    // minimim time separation between two subsequent signal 
    static G4bool storeOnlyTheFirstTimeHit;  // if true, only the hit that happened first will be
                                             // stored, anything else will be ignored 
                                             // (usefull for some special studies, not for a serious simulation)
    static G4int storeOnlyEventsWithMuonsDecayedInDetID;
    static G4bool killAllElectrons;          // if true, all electron tracks will be deleted (usefull for the studies of the muon beam)
    static G4bool killAllPositrons;          // if true, all positron tracks will be deleted (usefull for the studies of the muon beam)
    static G4bool killAllGammas;             // 
    static G4bool killAllNeutrinos;          //
  //cks    static G4bool includeMuoniumProcesses;   // If true, includes Muonium formation and all
  //cks                                             // other Mu-related processes in the simulation
    static G4bool boolG4GeneralParticleSource; // if true, G4GeneralParticleSource  will be initialised instead of G4ParticleGun
                                             //           - needed for the radioactive source
    static G4bool boolG4OpticalPhotons;      // if true, optical photons will be used (in the sensitive scintillators)
    static G4bool boolG4OpticalPhotonsUnprocess;  // if true, optical photons will not be processed - it might be
                                             //      usefull if the user wants to preselect some interesting events,
                                             //      and then to run the time-consuming processing of opt. photons
                                             //      again only for the interesting events.  This way the random number
                                             //      generator can generate the reproducible events.
                                             //      This option only works with   "boolG4OpticalPhotons".
    static G4bool field_DecayWithSpin;       // if true, then the routins for calculating the magnetic field will
                                             // use more precise argument.  This variable is set to "true" by
                                             // the SteppinAction and reset to "false" in the GetFieldValue.
                                             // It is being changed on step by step basis.
    static G4int nrOfEventsToBeGenerated;    // Nr of events to be simulated in this run (set by /run/beamOn command)
    static G4bool finiteRiseTimeInScintillator; // if true, rise time will be simulated in the scintillator.  For some strange
                                             // reason, Geant4 requires that this is specifically allowed.  We set it true
                                             // as soon as "FASTSCINTILLATIONRISETIME" or "SLOWSCINTILLATIONRISETIME" is set.
    static G4int maximumTimePerEvent;        // maximum time (in seconds) allowed for an event simulation - if exceeded, kill the event and
                                             // proceed to the next one.
    static G4int maximumNrOfStepsPerTrack;   // Maximum number of steps per track - if exceeded, kill the track and proceed with a next one.
  private:
    static musrParameters* pointerToParameters;
    G4bool boolG4RegionRequested;  // variable used internally just to check that no new volume is defined after
                                   // a G4Region has been requested - perhaps unnecessary check, but just to be on the safe side

};

#endif

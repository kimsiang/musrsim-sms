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

#ifndef musrDetectorMessenger_h
#define musrDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class musrDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class musrDetectorMessenger: public G4UImessenger
{
  public:
    musrDetectorMessenger(musrDetectorConstruction*);
   ~musrDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    musrDetectorConstruction* myDetector;
    
    G4UIdirectory*             musrDir;
    G4UIdirectory*             detDir;
    G4UIdirectory*             runDir;
    G4UIcmdWithAString*        Ignore1Cmd;
    G4UIcmdWithAString*        Ignore2Cmd;
    G4UIcmdWithAnInteger*      RunIDSetCmd;
    G4UIcmdWithAnInteger*      RandomOptionCmd;
    G4UIcmdWithAnInteger*      HowOftenToPrintEventCmd;
    G4UIcmdWithAnInteger*      RndmEventToSaveSeedsCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;

  public:

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


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

#ifndef musrEventAction_h
#define musrEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"
class G4Timer;
class G4Event;
class musrMagneticField;
class musrTabulatedField3D;
class musrTabulatedField2D;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class musrEventAction : public G4UserEventAction
{
  public:
    musrEventAction();
   ~musrEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    static musrEventAction* GetInstance();

  private:
  //  pointer to this class
    static musrEventAction* pointer;
  //
    G4Timer* timer;
  //    Variables for the time-dependent magnetic field
    G4bool   timeDependentField;
    time_t   timeOfRunStart;

  public:
    static G4int  nHowOftenToPrintEvent;
    static G4double maximumRunTimeAllowed;
};

#endif

    

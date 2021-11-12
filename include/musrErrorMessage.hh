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

#ifndef musrErrorMessage_h
#define musrErrorMessage_h 1

#include <map>
#include "globals.hh"
enum SEVERITY {INFO, WARNING, SERIOUS, FATAL};

typedef struct 
{
    SEVERITY mesSeverity;
    int nTimes;
} ErrorStruct;

class musrErrorMessage {
  public:
    musrErrorMessage();
   ~musrErrorMessage();
    static musrErrorMessage* GetInstance();
    void musrError(SEVERITY severity, G4String message, G4bool silent);
    void PrintErrorSummary();

  private:
    static musrErrorMessage* pointerToErrors;
    G4int nErrors;
  //    std::map<G4String,int> ErrorMapping;
    std::map<G4String,ErrorStruct> ErrorMapping;
    std::map<SEVERITY,G4String> severityWord;
};
#endif

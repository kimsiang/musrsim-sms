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

#include "musrErrorMessage.hh"

musrErrorMessage::musrErrorMessage():nErrors(1)
{
  pointerToErrors=this;
  severityWord[INFO]="INFO";
  severityWord[WARNING]="WARNING";
  severityWord[SERIOUS]="SERIOUS";
  severityWord[FATAL]="FATAL";
}

musrErrorMessage::~musrErrorMessage() {}

musrErrorMessage* musrErrorMessage::pointerToErrors=NULL;
musrErrorMessage* musrErrorMessage::GetInstance() {
  return pointerToErrors;
}

void musrErrorMessage::musrError(SEVERITY severity, G4String message, G4bool silent) {
  std::map<G4String,ErrorStruct>::iterator it;
  it = ErrorMapping.find(message);
  if (it == ErrorMapping.end()) {      // The error message is called for the first time
    ErrorStruct actualErrorMessage;
    actualErrorMessage.mesSeverity = severity;
    actualErrorMessage.nTimes = 1;
    ErrorMapping[message]=actualErrorMessage;
    it = ErrorMapping.find(message);
    G4cout<<"!!!"<<severityWord[severity]<<"!!! "<<message<<"  (First time occurrence)"<<G4endl;
    nErrors++;
  }
  else {                                // The error message is called for more than the first time
    (*it).second.nTimes++;
  }

  // Print out the error message if required
  if ((!silent)||(severity==FATAL)) { 
    if ((*it).second.nTimes>1) {
      G4cout<<"!!!"<<severityWord[severity]<<"!!! "<<message
	    <<" ("<<(*it).second.nTimes<<" occurences)"<<G4endl;
    }
  }
  
  if (severity==FATAL) {
    G4cout<<"S T O P     F O R C E D!"<<G4endl;
    exit(1);
  }
}



void musrErrorMessage::PrintErrorSummary() {
  std::map<G4String,ErrorStruct>::iterator it;
  G4cout<<"------   ERROR SUMMARY:  ----------------------------------------------------------------"<<G4endl;
  for (G4int i=0; i<4; i++) {
    for ( it=ErrorMapping.begin() ; it != ErrorMapping.end(); it++ ) {
      if ((*it).second.mesSeverity==i) {
	G4cout<<severityWord[(*it).second.mesSeverity]<<" ("<<(*it).second.nTimes<<" times):"
	      << (*it).first <<G4endl;
      }
    }
  }
  G4cout<<"-----------------------------------------------------------------------------------------"<<G4endl;
}

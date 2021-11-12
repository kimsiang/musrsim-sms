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

#include "musrQuadrupole.hh"
#include "musrParameters.hh"

musrQuadrupole::musrQuadrupole(G4double halflengthVal, G4double fieldRadiusVal, G4double gradientVal, G4double fringeFactorVal, G4LogicalVolume* logVolume, G4ThreeVector positionOfTheCenter) : F04ElementField(positionOfTheCenter, logVolume)
{    
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      Quadrupole field "
	 << "\n-----------------------------------------------------------"
         << G4endl;
  
  gradient = gradientVal;                  //*(tesla/m);
  fieldRadius = fieldRadiusVal;            //*mm;
  halflength = halflengthVal;              //*mm;
  enge = BLEngeFunction(ENGE_QUAD);

  G4double fringeFactor = fringeFactorVal;   // the default should be 1.
  G4bool fringe = true;  if (fringeFactor==0.) {fringe=false;}
  fringeDepth = fringeFactor * fieldRadius * 2.0;

  if (!fringe) {
    enge.set(0,0,0,0,0,0);
    fringeMaxZ = halflength;
  }
  else {
    for(int i=0; i<1000; ++i) {
      fringeMaxZ = i*fieldRadius/10.0 + halflength;
      if(enge((fringeMaxZ-halflength)/fringeDepth) < FRINGE_ACCURACY)   break;
    }
  }

  G4cout << " Field gradient set to "<< gradient/(CLHEP::tesla/CLHEP::m) << " T/m"<< G4endl;
  G4cout << " Field radius set to " << fieldRadius/CLHEP::mm << " mm"<< G4endl;
  G4cout << " Field halflength set to " << halflength/CLHEP::mm << " mm"<<G4endl;
  G4cout << " Field fringeDepth set to " << fringeDepth/CLHEP::mm << " mm"<<G4endl;
  G4cout << " Field fringeMaxZ set to " << fringeMaxZ/CLHEP::mm << " mm"<<G4endl;
  G4cout << "\n-----------------------------------------------------------" << G4endl;
}


void musrQuadrupole::addFieldValue(const G4double point[4], G4double *field ) const {
  G4ThreeVector global(point[0],point[1],point[2]);
  G4ThreeVector local;
  
  local = global2local.TransformPoint(global);
  G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
  if (r > fieldRadius || fabs(local[2]) > fringeMaxZ) { return;}
  
  // apply enge() to the scalar potential phi=-G0*x*y*enge(z);
  // B is minus its gradient. Handle both edges properly.
  G4double G0 = gradient;
  G4double f,fp;
  if (fringeDepth!=0) {
    double fringeZ = (fabs(local[2])-halflength)/fringeDepth; 
    f = enge(fringeZ);
    fp = enge.prime(fringeZ)/fringeDepth;
  }
  else {
    f = ( fabs(local[2]) > halflength) ? 0:1;
    fp = 0;
  }
  G4ThreeVector B(G0*f*local[1],G0*f*local[0],G0*fp*local[0]*local[1]);
  if(local[2] < 0.0) B[2] = -B[2];
  
  G4ThreeVector finalField(B[0],B[1],B[2]);
  finalField = global2local.Inverse().TransformAxis(finalField);
  
  field[0] += finalField.x();
  field[1] += finalField.y();
  field[2] += finalField.z();

  //  G4cout<<"musrQuadrupole.cc: field: ("<<field[0]/tesla<<","<<field[1]/tesla<<","<<field[2]/tesla<<")"<<G4endl;

}



G4double musrQuadrupole::GetNominalFieldValue() {
  return gradient;
}

void musrQuadrupole::SetNominalFieldValue(G4double newFieldValue) {
  //  // Rescale the magnetic field for a new value of the magnetic field
  gradient=newFieldValue;
  G4cout<<"musrQuadrupole.cc:   gradient changed to="<< gradient/(CLHEP::tesla/CLHEP::m)<<" T/m"<<G4endl;
}


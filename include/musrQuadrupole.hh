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

#ifndef musrQuadrupole_h
#define musrQuadrupole_h 1

#include "globals.hh"
#include "F04ElementField.hh"
#include "F04GlobalField.hh"
//#include "G4ios.hh"
#include "BLEngeFunction.hh"
//#include <fstream>
//#include <vector>
//#include <cmath>
const G4double FRINGE_ACCURACY=1.0e-4;

class musrQuadrupole : public F04ElementField
{
public:
  musrQuadrupole(G4double halflengthVal, G4double fieldRadiusVal, G4double gradientVal, G4double fringeFactorVal, G4LogicalVolume* logVolume, G4ThreeVector positionOfTheCenter);
  // musrQuadrupole is based on "BLCMDgenericquad" class of the G4beamline package.
  // 

  ///  Destructor.
  virtual ~musrQuadrupole() {}

  ///  addFieldValue() adds the field for this solenoid into field[].
  ///  point[] is in global coordinates.
  void  addFieldValue( const  G4double Point[4], G4double* field) const;

  G4double GetNominalFieldValue();
  void SetNominalFieldValue(G4double newFieldValue);
  
  //  getWidth(), getHeight(), getLength(),  return the dimensions of the field 
  // (used to define the boundary of the field)
  virtual G4double getWidth()  { return 2*fieldRadius; }   // x coordinate
  virtual G4double getHeight() { return 2*fieldRadius; }   // y coordinate
  virtual G4double getLength() { return 2*fringeMaxZ; }    //{ return 2*halflength; }   // z coordinate


private:
  G4double gradient;
  G4double fieldRadius;
  G4double halflength;

  G4double fringeMaxZ;
  G4double fringeDepth;

  BLEngeFunction enge;
};

#endif

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

#ifndef UNIFORM_BFIELD_HH
#define UNIFORM_BFIELD_HH
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "F04ElementField.hh"
#include "F04GlobalField.hh"

// UniformField implements a constant electromagnetic field in any direction. TS

class musrUniformField : public F04ElementField
{
  public:

    musrUniformField(G4double EMF[6], G4double half_X, G4double half_Y, G4double half_Z, G4LogicalVolume*, G4ThreeVector);

    virtual ~musrUniformField() {}

    // TS: Define the two newly added VIRTUAL functions of F04ElementField
    G4double GetNominalFieldValue();
    void SetNominalFieldValue(G4double newFieldValue);
    
    virtual G4double getLength() { return fieldLength; }
    virtual G4double getWidth()  { return fieldWidth; }
    virtual G4double getHeight() { return fieldHeight; }

    G4bool isOutside(G4ThreeVector& local) const;
    G4bool isWithin (G4ThreeVector& local) const;

    void addFieldValue(const G4double point[4], G4double field[6]) const;

  private:

    G4double EMfield[6];

    G4double fieldLength;
    G4double fieldWidth;
    G4double fieldHeight;
};

#endif

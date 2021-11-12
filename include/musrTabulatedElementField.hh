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

#ifndef musrTabulatedElementField_h
#define musrTabulatedElementField_h 1

#include "globals.hh"
#include "F04ElementField.hh"
#include "F04GlobalField.hh"
#include "G4ios.hh"
#include <fstream>
#include <vector>
#include <cmath>


class musrTabulatedElementField : public F04ElementField
{
public:
  musrTabulatedElementField(const char* filename, const char* fldTableType, G4double fieldValue, G4LogicalVolume* logVolume, G4ThreeVector positionOfTheCenter);
  //  "lenUnit" is the unit in which the grid coordinates are specified in the table
  //  "fieldNormalisation" is the normalisation that has to be applied on the field values in the table
  //                       such that the values correspond do 1T nominal value
  //  "fieldValue" is the field value (in T) that is required (i.e. values normalised to 1T will be
  //                       multiplied by this value).

  ///  Destructor.
  virtual ~musrTabulatedElementField() {}

  ///  addFieldValue() adds the field for this solenoid into field[].
  ///  point[] is in global coordinates.
  void  addFieldValue( const  G4double Point[4], G4double* field) const;
  void  addFieldValue2D( const  G4double Point[4], G4double* field) const;
  void  addFieldValue3D( const  G4double Point[4], G4double* field) const;

  G4double GetNominalFieldValue();
  void SetNominalFieldValue(G4double newFieldValue);

  //  getWidth(), getHeight(), getLength(),  return the dimensions of the field
  // (used to define the boundary of the field)
  virtual G4double getWidth()  {return maximumWidth;}    // x coordinate
  virtual G4double getHeight() {return maximumHeight;}   // y coordinate
  virtual G4double getLength() {return maximumLength;}   // z coordinate


private:
  // Storage space for the table
  std::vector< std::vector< std::vector< double > > > xField;
  std::vector< std::vector< std::vector< double > > > yField;
  std::vector< std::vector< std::vector< double > > > zField;
  std::vector< std::vector< double > > xField2D;
  std::vector< std::vector< double > > zField2D;
  // The dimensions of the table
  int nx,ny,nz; 
  // The units of the field
  char fieldTableType[100];
  G4String fUnit;
  double fieUnit;
  char   fldType;
  int  fldDim;
  // The physical limits of the defined region
  double minimumx, maximumx, minimumy, maximumy, minimumz, maximumz;
  // The physical extent of the defined region
  double dx, dy, dz;
  double ffieldValue;
  double maximumWidth, maximumHeight, maximumLength;
  int    symmetryType;  // this variable defines whether (and how) should be the field extended
                        // in the case it is defined just in one octant of the Cartesian coordinates.
  char   variableIncreasingOrder[100];
  int    jx, jy, jz;
  void   Invert(const char* indexToInvert);

};

#endif

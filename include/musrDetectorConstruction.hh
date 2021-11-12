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

#ifndef musrDetectorConstruction_h
#define musrDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4OpticalSurface.hh"
#include <map>

class G4Tubs;
class G4Box;
class G4Cons;
class G4Trd;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class musrDetectorMessenger;
class musrScintSD;
class musrMuEnergyLossLandau;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class musrDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  musrDetectorConstruction(G4String steeringFileName);
  ~musrDetectorConstruction();

public:
  
  G4VPhysicalVolume* Construct();


  void UpdateGeometry();
  void SetInputParameterFileName(G4String fileName) {parameterFileName=fileName;};
  void ReportGeometryProblem(char myString[501]);
  void ReportProblemInStearingFile(char* myString);
  G4Material* CharToMaterial(char myString[100]);
  G4LogicalVolume* FindLogicalVolume(G4String LogicalVolumeName);
  G4VPhysicalVolume* FindPhysicalVolume(G4String PhysicalVolumeName);
  void SetColourOfLogicalVolume(G4LogicalVolume* pLogVol,char* colour);

private:
  G4String           parameterFileName;  // name of the file with the geometry parameters
  G4bool             checkOverlap;  // parameter to check ovelaping volumes
  G4double           largestAcceptableStep;  // parameter defining largest step in the magnetic field

  musrScintSD* aScintSD;  
  musrDetectorMessenger* detectorMessenger;  // pointer to the Messenger

  std::map<std::string,G4RotationMatrix*> pointerToRotationMatrix;
  std::map<std::string,double*>           pointerToArray;
  std::map<std::string,double*>::iterator iterArray;
  std::map<std::string,G4FieldManager*> pointerToField;

  std::map<std::string,G4MaterialPropertiesTable*> materialPropertiesTableMap;
  std::map<std::string,G4MaterialPropertiesTable*>::iterator  itMPT;

  
private:
  void DefineMaterials(); 
};

#endif

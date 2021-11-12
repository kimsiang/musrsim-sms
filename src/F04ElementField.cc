//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//

#include "G4GeometryManager.hh"
#include "F04ElementField.hh"
#include "F04GlobalField.hh"
#include "musrParameters.hh"
#include "musrErrorMessage.hh"

G4Navigator* F04ElementField::aNavigator;

F04ElementField::F04ElementField(G4ThreeVector c, G4LogicalVolume* lv)
{ 
  elementFieldName="NAME_NOT_DEFINED";
  center = c;

  minX = minY = minZ = -DBL_MAX;
  maxX = maxY = maxZ =  DBL_MAX;

  //  G4cout<<"Kamil: F04GlobalField: addElementField() will be called: "<<G4endl;
  F04GlobalField::getObject()->addElementField(this);

  color = "1,1,1";

  userLimits = new G4UserLimits();

  lvolume = lv;
  lvolume->SetVisAttributes(getVisAttribute(color));

  maxStep = 1*CLHEP::mm;

  userLimits->SetMaxAllowedStep(maxStep);

  userLimits->SetUserMaxTrackLength(500.*CLHEP::m);
  userLimits->SetUserMaxTime(10*CLHEP::ms);
//  userLimits->SetUserMinEkine(0.1*MeV);
//  userLimits->SetUserMinRange(1*mm);

  lvolume->SetUserLimits(userLimits);
}

void F04ElementField::construct()
{
  G4Navigator* theNavigator = 
                    G4TransportationManager::GetTransportationManager()->
                                                 GetNavigatorForTracking();

  if (!aNavigator) {
     aNavigator = new G4Navigator();
     if ( theNavigator->GetWorldVolume() )
               aNavigator->SetWorldVolume(theNavigator->GetWorldVolume());
   }

  G4GeometryManager* geomManager = G4GeometryManager::GetInstance();

  if (!geomManager->IsGeometryClosed()) {
     geomManager->OpenGeometry();
     geomManager->CloseGeometry(true);
  }

  G4cout<<"F04ElementField: center="<<center.x()/CLHEP::mm<<", "<<center.y()/CLHEP::mm<<", "<<center.z()/CLHEP::mm<<" mm"<<G4endl;
  aNavigator->LocateGlobalPointAndSetup(center,0,false);

  G4TouchableHistoryHandle fTouchable = aNavigator->
                                         CreateTouchableHistoryHandle();

  G4int depth = fTouchable->GetHistoryDepth();
  for (G4int i = 0; i<depth; ++i) {
    if(fTouchable->GetVolume()->GetLogicalVolume() == lvolume)break;
    fTouchable->MoveUpHistory();
  }
  // Check that the point of origin of the field matches with the logical volume, to which it is assigned:
  G4String volumeName=lvolume->GetName();
  if (fTouchable->GetVolume()->GetLogicalVolume() != lvolume) {
    char eMessage[200];
    sprintf(eMessage,"F04ElementField.cc::construct(): Centre (point of origin) of the field outside the assigned logical volume \"%s\".",
	    volumeName.c_str());
    musrErrorMessage::GetInstance()->musrError(FATAL,eMessage,false);
  }

  //            G4cout<<"+!+!+! global2local VOLUME NAME: "<<fTouchable->GetVolume()->GetName()<<G4endl;
  // Set global2local transform.  The centre of the transformation is set to the centre of the volume, not
  // to the point "center", as one could naively expect.  This is corrected a few lines later.
  global2local = fTouchable->GetHistory()->GetTopTransform();

  // Print out the point of the origin of the field in the local coordinate system of the logical volume:
  G4ThreeVector local_center = global2local.TransformPoint(center);
  G4cout<<"\t==> "<<elementFieldName<<" in vol.=\""<<volumeName<<"\", center(local coord. system): "
	<<local_center.x()/CLHEP::mm<<", "<<local_center.y()/CLHEP::mm<<", "<<local_center.z()/CLHEP::mm<<" mm."<<G4endl;

  // Now move the centre of the transformation such that it coincides with the point "center":
  global2local*=G4AffineTransform(-local_center);

  //  After moving the centre of the transformation, the point "local_center" should be set to 0., 0., 0. in
  //  the following print-out message:
  //  local_center = global2local.TransformPoint(center);
  //  G4cout<<"\t==> "<<elementFieldName<<"  (volume=\""<<volumeName<<"\", center of the field in local coordinate syst 2: "
  // 	<<local_center.x()/mm<<", "<<local_center.y()/mm<<", "<<local_center.z()/mm<<" mm.)"<<G4endl;

  // set global bounding box
  G4double local[4], global[4];

  G4ThreeVector globalPosition;
  local[3] = 0.0;
  for (int i=0; i<2; ++i) {
      local[0] = (i==0 ? -1.0 : 1.0) * getWidth()/2.;
      for (int j=0; j<2; ++j) {
          local[1] = (j==0 ? -1.0 : 1.0) * getHeight()/2.;
          for (int k=0; k<2; ++k) {
              local[2] = (k==0 ? -1.0 : 1.0) * getLength()/2.;
              G4ThreeVector localPosition(local[0],local[1],local[2]);
              globalPosition = 
                         global2local.Inverse().TransformPoint(localPosition);
              global[0] = globalPosition.x();
              global[1] = globalPosition.y();
              global[2] = globalPosition.z();
              setGlobalPoint(global);
           }
      }
  }
}

G4VisAttributes* F04ElementField::getVisAttribute(G4String color)
{
   G4VisAttributes* p = NULL;
   if(color.size() > 0 &&
     (isdigit(color.c_str()[0]) || color.c_str()[0] == '.')) {
        G4double red=0.0, green=0.0, blue=0.0;
        if (sscanf(color.c_str(),"%lf,%lf,%lf",&red,&green,&blue) == 3) {
           p = new G4VisAttributes(true,G4Color(red,green,blue));
        } else {
           G4cout << " Invalid color " << color << G4endl;
        }
   }

   if (!p) p = new G4VisAttributes(G4VisAttributes::Invisible);
   p->SetDaughtersInvisible(false);

   return p;
}


void F04ElementField::SetEventNrDependentField(G4double initialField, G4double finalField, G4int nrOfSteps) {
  G4double eventNrStep        =  float(musrParameters::nrOfEventsToBeGenerated)/(nrOfSteps);
  G4double fieldStep          =  (finalField-initialField)/(nrOfSteps-1);
  //  G4cout<<"musrParameters::nrOfEventsToBeGenerated="<<musrParameters::nrOfEventsToBeGenerated<<G4endl;
  //  G4cout<<"fieldStep="<<fieldStep<<"     eventNrStep="<<eventNrStep<<G4endl;
  for (G4int i=1; i<=nrOfSteps; i++) {
    G4int    eventNumber = int(i*eventNrStep);
    G4double field       = initialField + i*fieldStep;
    changeFieldInStepsMap[eventNumber]=field;
  }
  
  G4cout << "Setting field in steps for field "<<elementFieldName<<G4endl;
  std::map<G4int,G4double>::iterator it;
  for ( it=changeFieldInStepsMap.begin() ; it != changeFieldInStepsMap.end(); it++ ) {
    G4cout << "Field will be changed at event "<< (*it).first << " to the value of " << (*it).second/CLHEP::tesla<<" T" << G4endl;
    //    G4double nominalFieldValue=it->second;
    //    it->SetNominalFieldValue(nominalFieldValue);
  }
}


void F04ElementField::SetElementFieldValueIfNeeded(G4int eventNr) {
  std::map<G4int,G4double>::iterator itr;
  itr = changeFieldInStepsMap.find(eventNr);
  if (itr==F04ElementField::changeFieldInStepsMap.end()) {
    // eventNr was not found in the map ==> field is not going to change at this eventNr
  } 
  else {
    G4double newFieldValue = itr->second;
    SetNominalFieldValue(newFieldValue);
    //	G4cout<<"Nominal Field changed for "<<G4endl;
  }
}

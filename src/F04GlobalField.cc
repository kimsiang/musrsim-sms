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

#include <time.h>

#include "Randomize.hh"
#include "G4TransportationManager.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"

#include "F04GlobalField.hh"
// #include "F04GlobalFieldOldPoints.hh"     //cks F04GlobalFieldOldPoints were introduced in order to keep the last few
//                                                 values of the field, so that they do not have to be recalculated, if the
//                                                 field at recent point is requested.  It turned out not to be useful
//                                                 (i.e. it did not speed-up the simulation)
#include "musrRootOutput.hh"

F04GlobalField* F04GlobalField::object = 0;

F04GlobalField::F04GlobalField() : G4ElectroMagneticField(),
                             minStep(0.01*CLHEP::mm), deltaChord(3.0*CLHEP::mm),
                             deltaOneStep(0.01*CLHEP::mm), deltaIntersection(0.1*CLHEP::mm),
                             epsMin(2.5e-7*CLHEP::mm), epsMax(0.05*CLHEP::mm),
                             fEquation(0), fFieldManager(0), 
                             fFieldPropagator(0), fStepper(0), fChordFinder(0)
// F04GlobalField::F04GlobalField() : G4MagneticField(),
//                             minStep(0.01*mm), deltaChord(3.0*mm),
//                             deltaOneStep(0.01*mm), deltaIntersection(0.1*mm),
//                             epsMin(2.5e-7*mm), epsMax(0.05*mm),
//                             fEquation(0), fFieldManager(0),
//  			     fFieldPropagator(0), fStepper(0), fChordFinder(0)
{

  fFieldMessenger = new F04FieldMessenger(this);

  fields = new FieldList();

  fStepperType = 4 ;       // ClassicalRK4 is default stepper

  //  set object

  object = this;

  updateField();

  //cks  myOldFieldPoints = F04GlobalFieldOldPoints();

}

F04GlobalField::~F04GlobalField()
{
  // cks 2009_05_14
  //      For some reason, the "clear()" on the next line was causing the crash in some cases
  //      (probably when more "Elementary fields" were created).  The problem appeared
  //      when running the 1050.mac file for the lem4.  Fixed by commenting out the clear():
  //  clear();
  delete fFieldMessenger;
  if (fEquation)        delete fEquation;
  //  if (fFieldManager)    delete fFieldManager;
  //  if (fFieldPropagator) delete fFieldPropagator;
  if (fStepper)         delete fStepper;
  if (fChordFinder)     delete fChordFinder;
}

void F04GlobalField::updateField()
{
  first = true;

  nfp = 0;
  fp = 0;

  clear();

  //  Construct equ. of motion of particles through B fields
  //  fEquation = new G4Mag_EqRhs(this);
  //  Construct equ. of motion of particles through e.m. fields
  //  fEquation = new G4EqMagElectricField(this);
  //  Construct equ. of motion of particles including spin through B fields
  //  fEquation = new G4Mag_SpinEqRhs(this);
  //  Construct equ. of motion of particles including spin through e.m. fields
  fEquation = new G4EqEMFieldWithSpin(this);

  //  Get transportation, field, and propagator managers
  G4TransportationManager* fTransportManager =
         G4TransportationManager::GetTransportationManager();

  fFieldManager = GetGlobalFieldManager();

  fFieldPropagator = fTransportManager->GetPropagatorInField();

  //  Need to SetFieldChangesEnergy to account for a time varying electric
  //  field (r.f. fields)
  fFieldManager->SetFieldChangesEnergy(true);

  //  Set the field
  fFieldManager->SetDetectorField(this);

  //  Choose a stepper for integration of the equation of motion
  SetStepper();

  //  Create a cord finder providing the (global field, min step length,
  //  a pointer to the stepper)
  fChordFinder = new G4ChordFinder((G4MagneticField*)this,minStep,fStepper);

  // Set accuracy parameters
  fChordFinder->SetDeltaChord( deltaChord );

  fFieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

  fFieldManager->SetDeltaIntersection(deltaIntersection);

  fFieldPropagator->SetMinimumEpsilonStep(epsMin);
  fFieldPropagator->SetMaximumEpsilonStep(epsMax);

  //  G4cout << "Accuracy Parameters:" <<
  //            " MinStep=" << minStep <<
  //            " DeltaChord=" << deltaChord <<
  //            " DeltaOneStep=" << deltaOneStep << G4endl;
  //  G4cout << "                    " <<
  //            " DeltaIntersection=" << deltaIntersection <<
  //            " EpsMin=" << epsMin <<
  //            " EpsMax=" << epsMax <<  G4endl;
  
  fFieldManager->SetChordFinder(fChordFinder);

}

F04GlobalField* F04GlobalField::getObject()
{
  if (!object) new F04GlobalField();
  return object;
}

void F04GlobalField::SetStepper()
{
  if(fStepper) delete fStepper;

  switch ( fStepperType )
  {
    case 0:
//      fStepper = new G4ExplicitEuler( fEquation, 8 ); // no spin tracking
      fStepper = new G4ExplicitEuler( fEquation, 12 ); // with spin tracking
      G4cout << "G4ExplicitEuler is called" << G4endl;
      break;
    case 1:
//      fStepper = new G4ImplicitEuler( fEquation, 8 ); // no spin tracking
      fStepper = new G4ImplicitEuler( fEquation, 12 ); // with spin tracking
      G4cout << "G4ImplicitEuler is called" << G4endl;
      break;
    case 2:
//      fStepper = new G4SimpleRunge( fEquation, 8 ); // no spin tracking
      fStepper = new G4SimpleRunge( fEquation, 12 ); // with spin tracking
      G4cout << "G4SimpleRunge is called" << G4endl;
      break;
    case 3:
//      fStepper = new G4SimpleHeum( fEquation, 8 ); // no spin tracking
      fStepper = new G4SimpleHeum( fEquation, 12 ); // with spin tracking
      G4cout << "G4SimpleHeum is called" << G4endl;
      break;
    case 4:
//      fStepper = new G4ClassicalRK4( fEquation, 8 ); // no spin tracking
      fStepper = new G4ClassicalRK4( fEquation, 12 ); // with spin tracking
      G4cout << "G4ClassicalRK4 (default) is called" << G4endl;
      break;
    case 5:
//      fStepper = new G4CashKarpRKF45( fEquation, 8 ); // no spin tracking
      fStepper = new G4CashKarpRKF45( fEquation, 12 ); // with spin tracking
      G4cout << "G4CashKarpRKF45 is called" << G4endl;
      break;
    default: fStepper = 0;
  }
}

G4FieldManager* F04GlobalField::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
                                ->GetFieldManager();
}

void F04GlobalField::GetFieldValue(const G4double* point, G4double* field) const
{
  // NOTE: this routine dominates the CPU time for tracking.
  // Using the simple array fp[] instead of fields[] 
  // directly sped it up

  //    G4cout<<"GlobalField: Field requested at point ("<<point[0]<<","<<point[1]<<","<<point[2]<<")"<<G4endl;
  //cks   Check whether the requested point is the same as in the previous call
  //    if (myOldFieldPoints.ThisPointWasCalculatedRecently(point,field)) {return;}


  field[0] = field[1] = field[2] = field[3] = field[4] = field[5] = 0.0;

  // protect against Geant4 bug that calls us with point[] NaN.
  if(point[0] != point[0]) return;

  // (can't use nfp or fp, as they may change)
  if (first) ((F04GlobalField*)this)->setupArray();   // (cast away const)

  for (int i=0; i<nfp; ++i) {
      const F04ElementField* p = fp[i];
      if (p->isInBoundingBox(point)) {
         p->addFieldValue(point,field);
      }
  }

  // cks  NOT SURE WHETHER THE FOLLOWING WAS STILL NEEDED, HOWEVER REMOVED NOT 
  //      TO DISTURB THE LOW ENERGY MUONS.
  //      Set some small field if field is almost zero (to avoid internal problems 
  //      of Geant at almost zero magnetic fields).
  //  if (sqrt(field[0]*field[0]+field[1]*field[1]+field[2]*field[2])<0.00001*tesla) {
  //    field[2] = 0.00001*tesla;
  //  }
  
  // cks  myOldFieldPoints.StoreTheFieldPointForFutureReuse(point,field);
}

void F04GlobalField::clear()
{
  if (fields) {
     if (fields->size()>0) {
        FieldList::iterator i;
	//cks 2009_05_14 : The following line seems to cause problems (sometimes)
	//                 See the comment in F04GlobalField::~F04GlobalField() for more details.
	for (i=fields->begin(); i!=fields->end(); ++i)  delete *i;
        fields->clear();
     }
  }

  if (fp) delete[] fp;
  first = true;

  nfp = 0;
  fp = NULL;
}

void F04GlobalField::setupArray()
{
  first = false;
  nfp = fields->size();
  fp = new const F04ElementField* [nfp+1]; // add 1 so it's never 0
  for (int i=0; i<nfp; ++i) {
    fp[i] = (*fields)[i];
    //
    // Find the event numbers, for which the field changes, for each Element Field.
    // Then mark these event numbers in the Global Field, such that it can be efficiently
    // find out during the run-time, at which event number the field may change.
    std::map<G4int,G4double> localChangeFieldInStepsMap = fp[i] ->GetEventNrDependentField();
    std::map<G4int,G4double>::iterator it;
    for ( it=localChangeFieldInStepsMap.begin() ; it != localChangeFieldInStepsMap.end(); it++ ) {
      G4int eventNr = it->first;
      G4cout<<"globalChangeFieldInStepsMap  set for event number "<<eventNr<<G4endl;
      globalChangeFieldInStepsMap[eventNr]=true;
    }
  }
  musrRootOutput::GetRootInstance()->SetNrFieldNomVal(nfp);
}


void F04GlobalField::CheckWhetherAnyNominalFieldValueNeedsToBeChanged(G4int eventNumber) {
  if (globalChangeFieldInStepsMap[eventNumber]) {
    //    G4cout<<"We should check each Element Field Object whether its field needs to be changed:"<<G4endl;
    G4int jjj=0;
    FieldList::iterator i;
    for (i=fields->begin(); i!=fields->end(); ++i) {

      // Set the nominal field value for the given field, if that has been requested for the given field
      (*i)->SetElementFieldValueIfNeeded(eventNumber);

      // Get the nominal field value for the given field and store it in the Root output 
      G4double nomFieldValue = (*i)->GetNominalFieldValue();
      musrRootOutput::GetRootInstance()->SetFieldNomVal(jjj,nomFieldValue);
      jjj++;
    }
  }
}



//  Print field value at all points requested by the user:
void F04GlobalField::PrintFieldAtRequestedPoints() const {
  G4ThreeVector p;
  G4double point[4];
  G4double delta=0.1*CLHEP::mm;
  G4double point2[4]; G4double point3[4]; G4double point4[4]; 
  G4double Bfi[6]={0,0,0,0,0,0};
  G4double BfiX[6]={0,0,0,0,0,0};
  G4double BfiY[6]={0,0,0,0,0,0};
  G4double BfiZ[6]={0,0,0,0,0,0};
  for (unsigned int i=0; i < pointsAtWhichUserWantsToPrintFieldValue.size(); i++) {
    p        = pointsAtWhichUserWantsToPrintFieldValue[i];
    point[0] = p.x();
    point[1] = p.y();
    point[2] = p.z();
    point[3] = 0.;
    object->GetFieldValue(point,Bfi);
    //    printf ("   Magnetic Field at %f, %f, %f mm  is   B= %10.10f, %10.10f, %10.10f tesla.\n", 
    //	    point[0]/mm,point[1]/mm,point[2]/mm,Bfi[0]/tesla,Bfi[1]/tesla,Bfi[2]/tesla);
        printf ("  Field at (%.2f, %.2f, %.2f) mm is:  B = (%0.10g, %0.10g, %0.10g) T,  E = (%0.10g, %0.10g, %0.10g) kV/mm\n",
        point[0]/CLHEP::mm, point[1]/CLHEP::mm, point[2]/CLHEP::mm,
        Bfi[0]/CLHEP::tesla,         Bfi[1]/CLHEP::tesla,         Bfi[2]/CLHEP::tesla,
        Bfi[3]/(CLHEP::kilovolt/CLHEP::mm), Bfi[4]/(CLHEP::kilovolt/CLHEP::mm), Bfi[5]/(CLHEP::kilovolt/CLHEP::mm));
  }

  if (pointsAtWhichUserWantsToPrintFieldDerivative.size()>0) {
    G4cout<<G4endl<<"       ===="<<G4endl<<"Derivatives of the magnetic field"<<G4endl;
    G4cout<<"   x    y    z (mm)          dBx/dx  dBx/dy  dBx/dz         dBy/dx  dBy/dy  dBy/dz        dBz/dx  dBz/dy  dBz/dz (G/mm)"<<G4endl; 
    //    G4cout<<"tesla = "<<tesla<<"  gauss="<<gauss<<"  ==>  tesla = "<<tesla/gauss<<"gauss"<<G4endl;
  }
  for (unsigned int i=0; i < pointsAtWhichUserWantsToPrintFieldDerivative.size(); i++) {
    p        = pointsAtWhichUserWantsToPrintFieldDerivative[i];
    point[0] = p.x();
    point[1] = p.y();
    point[2] = p.z();
    point[3] = 0.;
    point2[0] = point[0]+delta; point2[1] = point[1];       point2[2] = point[2];       point2[3] = point[3];
    point3[0] = point[0];       point3[1] = point[1]+delta; point3[2] = point[2];       point3[3] = point[3];
    point4[0] = point[0];       point4[1] = point[1];       point4[2] = point[2]+delta; point4[3] = point[3];

    object->GetFieldValue(point,Bfi);
    object->GetFieldValue(point2,BfiX);
    object->GetFieldValue(point3,BfiY);
    object->GetFieldValue(point4,BfiZ);
    //    printf ("   Magnetic Field at %f, %f, %f mm  is   B= %10.10f, %10.10f, %10.10f tesla.\n", 
    //	    point[0]/mm,point[1]/mm,point[2]/mm,Bfi[0]/tesla,Bfi[1]/tesla,Bfi[2]/tesla);
    printf ("  %.2f  %.2f  %.2f          %.2f  %.2f  %.2f             %.2f  %.2f  %.2f             %.2f  %.2f  %.2f\n",
	    point[0]/CLHEP::mm, point[1]/CLHEP::mm, point[2]/CLHEP::mm,
	    (BfiX[0]-Bfi[0])/CLHEP::gauss/delta*CLHEP::mm,  (BfiY[0]-Bfi[0])/CLHEP::gauss/delta*CLHEP::mm,  (BfiZ[0]-Bfi[0])/CLHEP::gauss/delta*CLHEP::mm,
	    (BfiX[1]-Bfi[1])/CLHEP::gauss/delta*CLHEP::mm,  (BfiY[1]-Bfi[1])/CLHEP::gauss/delta*CLHEP::mm,  (BfiZ[1]-Bfi[1])/CLHEP::gauss/delta*CLHEP::mm,
	    (BfiX[2]-Bfi[2])/CLHEP::gauss/delta*CLHEP::mm,  (BfiY[2]-Bfi[2])/CLHEP::gauss/delta*CLHEP::mm,  (BfiZ[2]-Bfi[2])/CLHEP::gauss/delta*CLHEP::mm );
  }
}



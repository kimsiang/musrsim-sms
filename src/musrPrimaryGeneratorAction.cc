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

#include "musrPrimaryGeneratorAction.hh"
#include "musrDetectorConstruction.hh"
#include "musrPrimaryGeneratorMessenger.hh"
#include "musrParameters.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"
#include "G4Gamma.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"
#include "time.h"
#include <iomanip>
#include "musrRootOutput.hh"   //cks for storing some info in the Root output file
#include "musrErrorMessage.hh"


G4bool musrPrimaryGeneratorAction::setRandomNrSeedAccordingEventNr=0;
G4bool musrPrimaryGeneratorAction::setRandomNrSeedFromFile=0;
G4bool musrPrimaryGeneratorAction::setRandomNrSeedFromFile_RNDM=0;
G4int  musrPrimaryGeneratorAction::nRndmEventToSaveSeeds=-2;
G4int  musrPrimaryGeneratorAction::lastEventID_in_pointerToSeedVector=0;

std::vector<int> * musrPrimaryGeneratorAction::pointerToSeedVector=NULL;
std::vector<int> * musrPrimaryGeneratorAction::GetPointerToSeedVector() {
  return pointerToSeedVector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
musrPrimaryGeneratorAction::musrPrimaryGeneratorAction(
                                            musrDetectorConstruction* musrDC)
  :musrDetector(musrDC), x0(0), y0(0), z0(-10*CLHEP::cm), xSigma(0), ySigma(0), zSigma(0), 
   rMaxAllowed(1e10*CLHEP::mm), zMinAllowed(-1e10*CLHEP::mm), zMaxAllowed(1e10*CLHEP::mm),
   meanArrivalTime(1./30000.*CLHEP::second), t0(0), tSigma(0),
   relativeRMaxAllowed(1e10*CLHEP::mm), 
   xMaxSource0(0), yMaxSource0(0), zMaxSource0(0),
   xMaxSource(1e10*CLHEP::mm), yMaxSource(1e10*CLHEP::mm), zMaxSource(1e10*CLHEP::mm),
   p0(0), pSigma(0), pMinAllowed(0), pMaxAllowed(1e10*CLHEP::mm),
   xangle0(0), yangle0(0), xangleSigma(0), yangleSigma(0), pitch(0),
   UnpolarisedMuonBeam(false), TransversalyUnpolarisedMuonBeam(false), xPolarisIni(1.), yPolarisIni(0.), zPolarisIni(0.),
   xDirection(0), yDirection(0), zDirection(1.),
   polarisFraction(1.),
   muonDecayTimeMin(-1), muonDecayTimeMax(-1), muonMeanLife(2197.03*CLHEP::ns),
   takeMuonsFromTurtleFile(false), z0_InitialTurtle(0),
   numberOfGeneratedEvents(0), 
   turtleMomentumBite(false), turtleMomentumP0(0.), turtleSmearingFactor(0.), turtleMomentumScalingFactor(1.)
						      //, firstCall(true)
{
  //create a messenger for this class
  gunMessenger = new musrPrimaryGeneratorMessenger(this);

  // create a vector for storing the event numbers as seeds for the random number generator
  pointerToSeedVector = new std::vector<int>;

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* muonParticle= particleTable->FindParticle("mu+");
  // cks delete muonMinusParticle= particleTable->FindParticle("mu-");
  // cks Implement also alpha and proton particles for the simulation of Juan Pablo Urrego
  alphaParticle= particleTable->FindParticle("alpha");
  protonParticle= particleTable->FindParticle("proton");
  turtleInterpretAxes="undefined";
  // csk

  G4int n_particle = 1;
  if (musrParameters::boolG4GeneralParticleSource) {
    G4cout<<"musrPrimaryGeneratorAction:  G4GeneralParticleSource is going to be initialised"<<G4endl;
    particleSource = new G4GeneralParticleSource ();
  }
  else {
    G4cout<<"musrPrimaryGeneratorAction:  G4ParticleGun is going to be initialised"<<G4endl;
    particleGun = new G4ParticleGun(n_particle);
    particleGun->SetParticleDefinition(muonParticle);
  }
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

musrPrimaryGeneratorAction::~musrPrimaryGeneratorAction()
{
  if (musrParameters::boolG4GeneralParticleSource) {delete  particleSource;}
  else {delete particleGun;}
  delete gunMessenger;
  if (takeMuonsFromTurtleFile) {fclose(fTurtleFile);}
  G4cout<<"musrPrimaryGeneratorAction:   Number of Generated Events = "<<numberOfGeneratedEvents<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void musrPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event.

  // Clear Root variables
  musrRootOutput* myRootOutput = musrRootOutput::GetRootInstance();
  myRootOutput->ClearAllRootVariables();  // Note that musrPrimaryGeneratorAction::GeneratePrimaries 
                                          // is called before the musrEventAction::BeginOfEventAction.
                                          // Therefore "ClearAllRootVariables" is called already here
                                          // (before the "SetInitialMuonParameters".

  // Set or read the seeds of random number generator 
  boolPrintInfoAboutGeneratedParticles=false;
  SetOrReadTheRandomNumberSeeds(anEvent);

  // If radioactive source is used, use  G4GeneralParticleSource :
  if (musrParameters::boolG4GeneralParticleSource) {
    particleSource->GeneratePrimaryVertex(anEvent);
    return;
  }

  G4double x, y, z;
  G4double p;
  G4double xangle, yangle;

  if (takeMuonsFromTurtleFile) {
    char  line[501];  
    G4int checkNrOfCounts=0;
    do {
      float xTmp, yTmp, xAngleTmp, yAngleTmp, pTmp;
      float dummy1, dummy2;
      int Ztmp=-1, Atmp=-1;
      fgets(line,500,fTurtleFile);
      if (feof(fTurtleFile)) {
	rewind(fTurtleFile);
	G4cout<<"End of TurtleFile, lets start from the beginning  (numberOfGeneratedEvents = "<<numberOfGeneratedEvents<<")"<<G4endl;
	fgets(line,500,fTurtleFile);
      }
      numberOfGeneratedEvents++;
      sscanf(&line[0],"%g %g %g %g %g %g %g %d %d",&xTmp,&xAngleTmp,&yTmp,&yAngleTmp,&pTmp,&dummy1,&dummy2,&Ztmp,&Atmp);
      if (turtleInterpretAxes!="undefined") swapTheAxisInTurtle(xTmp,xAngleTmp,yTmp,yAngleTmp);
      if (boolPrintInfoAboutGeneratedParticles) {
	G4cout<<"musrPrimaryGeneratorAction::GeneratePrimaries:  Turtle input for this event: "
	      <<xTmp<<", "<<xAngleTmp<<" "<<yTmp<<" "<<yAngleTmp<<" "<< pTmp<<G4endl;
      }
      //cks Implement also alpha and proton particles for the simulation of Juan Pablo Urrego
      if ((Ztmp==1)&&(Atmp==1)) {particleGun->SetParticleDefinition(protonParticle);}// G4cout<<"proton"<<G4endl;}
      else if ((Ztmp==2)&&(Atmp==4)) {particleGun->SetParticleDefinition(alphaParticle);}// G4cout<<"alpha particle"<<G4endl;}
      else if ((Ztmp==-1)&&(Atmp==-1)) {;}
      else {
	G4cout<<"musrPrimaryGeneratorAction: Unknown particle requested in the TURTLE input file: Z="
	      <<Ztmp<<", A="<<Atmp<<G4endl<<"S T O P     F O R C E D" << G4endl;
	G4cout<<xTmp<<", "<<xAngleTmp<<", "<<yTmp<<", "<<yAngleTmp<<", "<<pTmp<<", "<<dummy1<<", "<<dummy2<<", "<<Ztmp<<", "<<Atmp<<G4endl;
	exit(1);
      }
      //csk
      xangle = xAngleTmp*CLHEP::mrad;
      yangle = yAngleTmp*CLHEP::mrad;
      x      = xTmp*CLHEP::cm + (z0-z0_InitialTurtle)*tan(xangle) ;       // usually z0 is negative
      y      = yTmp*CLHEP::cm + (z0-z0_InitialTurtle)*tan(yangle) ;       // z0_InitialTurtle is the z0 at whith the turtle file was generated.
      p      = pTmp*CLHEP::GeV*turtleMomentumScalingFactor;
      // add some offset, if requested:
      x      = x + x0;
      y      = y + y0;
      // add some beam tilt, if requested:
      xangle = xangle + xangle0;
      yangle = yangle + yangle0;
      // add some beam pitch, if requested:
      if (pitch!=0) {
	xangle += - pitch * (x-x0);
	yangle += - pitch * (y-y0);
      }
      // add/remove some momentum smearing, if requested
      if (turtleMomentumBite) {
	p = turtleMomentumP0 - (turtleMomentumP0-p)*turtleSmearingFactor;
      }
      checkNrOfCounts++;
      if (checkNrOfCounts>1000) {
	G4cout<<"musrPrimaryGeneratorAction::GeneratePrimaries:  Too strict requirements on the r position!"<<G4endl;
      }
    } while( (x*x+y*y)>(rMaxAllowed*rMaxAllowed) );
    z=z0;
    //      G4cout<<"x,y,z=("<<x/mm<<","<<y/mm<<","<<z/mm<<"),  angles="<<xangle/mrad<<","<<yangle/mrad<<"  p="<<p/MeV<<G4endl;
  }

  else {          // Generate the starting position of the muon by random
    //  rMaxAllowed             ... maximal radius, within which the muon can be generated
    //  relativeRMaxAllowed     ... maximal radius, within which the muon can be generated, relative to x0 and y0 position
    //  xMaxSource0, yMaxSource0...
    //  xMaxSource, yMaxSource  ...
    //  x0, y0, z0              ... central point around which the muons are generated
    //  xSigma, ySigma, zSigma  ... sigma of the (gaussian) distributions of the beam
    //  x, y, z                 ... actual initial position of the generated muon
    
    G4int checkNrOfCounts=0;
    numberOfGeneratedEvents++;
    do {
      if (xSigma>0)      {x = G4RandGauss::shoot(x0,xSigma);}           //  Gaussian distribution
      else if (xSigma<0) {x = x0 + xSigma*(G4UniformRand()*2.-1.);}     //  Uniform step distribution
      else               { x = x0;}                                     //  Point-like

      if (ySigma>0)      {y = G4RandGauss::shoot(y0,ySigma);}
      else if (ySigma<0) {y = y0 + ySigma*(G4UniformRand()*2.-1.);}
      else               {y = y0;}

      if (zSigma>0)      {z = G4RandGauss::shoot(z0,zSigma);}
      else if (zSigma<0) {z = z0 + zSigma*(G4UniformRand()*2.-1.);}
      else               {z = z0;}

      checkNrOfCounts++;
      if (checkNrOfCounts>1000) {
	G4cout<<"musrPrimaryGeneratorAction::GeneratePrimaries:  Too strict requirements on the r or z position!"<<G4endl;
      }
    } while( 
	    ((x*x+y*y)>(rMaxAllowed*rMaxAllowed))||(z>zMaxAllowed)||(z<zMinAllowed) ||
	    (((x-x0)*(x-x0)+(y-y0)*(y-y0))>(relativeRMaxAllowed*relativeRMaxAllowed))||
	    (fabs(x-xMaxSource0)>xMaxSource)||(fabs(y-yMaxSource0)>yMaxSource)||(fabs(z-zMaxSource0)>zMaxSource)
	    ); // The generated muon has to stay within some well defined region,  e.g. within the beampipe
  
    // Now generate the momentum
    checkNrOfCounts=0;
    do {
      if (pSigma>0) {p = G4RandGauss::shoot(p0,pSigma);}
      else {p=p0;}
      checkNrOfCounts++;
      if (checkNrOfCounts>1000) {
	G4cout<<"musrPrimaryGeneratorAction::GeneratePrimaries:  Too strict requirements on the momentum!"<<G4endl;
      }
    } while ( (p>pMaxAllowed)||(p<pMinAllowed) );
    
    
    // Add some initial angle (px and py component of the momentum)
    if (xangleSigma>0) { xangle = G4RandGauss::shoot(xangle0,xangleSigma); }
    else { xangle = xangle0; }
    //  Add the beam tilt, which depends on the distance from the beam centre.
    //    if (xSigma>0) {xangle += - pitch * (x-x0)/xSigma; } 
    if (pitch!=0) {xangle += - pitch * (x-x0); } 

    if (yangleSigma>0) { yangle = G4RandGauss::shoot(yangle0,yangleSigma); }
    else { yangle = yangle0; }
    //  Add the beam tilt, which depends on the distance from the beam centre.
    //    if (ySigma>0) {yangle += - pitch * (y-y0)/ySigma; } 
    if (pitch!=0) {yangle += - pitch * (y-y0); } 
    
  }  // end of the part specific for the muons generated by random rather then from TURTLE
     
  //------- generate cos(Theta) distribution if zangleSigma < 0 and calculate final momentum
  G4double sinXangle, sinYangle, phi;
  G4double px, py, pz;
  if (zangleSigma<0) { 
    sinXangle=sqrt(G4UniformRand());
    phi=2*CLHEP::pi*G4UniformRand();
    px = p*sinXangle*cos(phi);
    py = p*sinXangle*sin(phi);
    pz = std::sqrt(p*p - px*px - py*py);
  }
  else {
    sinXangle=sin(xangle); 
    sinYangle=sin(yangle);       
    px = p*sinXangle;
    py = p*sinYangle;
    pz = std::sqrt(p*p - px*px - py*py);
    //  printf("px, py, pz = %f, %f, %f\n", px, py, pz);
  }
  
  // Calculate particle (muon) starting time
  G4double ParticleTime; //P.B. 13 May 2009
  if (tSigma>0)      {ParticleTime = G4RandGauss::shoot(t0,tSigma);}         //  Gaussian distribution       P.B. 13 May 2009
  else if (tSigma<0) {ParticleTime = t0 + tSigma*(G4UniformRand()*2.-1.);}   //  Uniform step distribution   P.B. 13 May 2009
  else               {ParticleTime = t0;}                                    //  Point-like                  P.B. 13 May 2009
  

  

  // Now rotate so that beam propagates along /gun/direction

  // Proper rotation of beam direction as follows:
  // hat{n} is beam direction unit vecrot and vec{p} is momentum vector
  // new momentum vector is vec{p'}=R vec{p}
  // where R is a rotation matrix around (hat{z} x hat{n}) - cross product
  // with an angle of acos(hat{z} dot hat{n}) - dot product
  // i.e. R is a rotation around (-n_y, n_x, 0) with angle acos(n_z). 

  if ((xDirection == 0) && (yDirection == 0)) {
    // Rotation does not work for beam direction along z.
    pz = zDirection * pz;
    // No change to the beam spot...
  } else {
    G4ThreeVector* PVec = new G4ThreeVector(px,py,pz);
    PVec->rotate(std::acos(zDirection),G4ThreeVector(-yDirection,xDirection,0));
    px = PVec->x();
    py = PVec->y();
    pz = PVec->z();

    // Rotate also beam spot
    G4ThreeVector* RVec = new G4ThreeVector(x-x0,y-y0,z-z0);
    RVec->rotate(std::acos(zDirection),G4ThreeVector(-yDirection,xDirection,0));
    x = x0+RVec->x();
    y = y0+RVec->y();
    z = z0+RVec->z();
  }

  // Assign spin 
  G4double xpolaris=0, ypolaris=0, zpolaris=0;
  if (UnpolarisedMuonBeam) {
    // for genarating random numbers on the sphere see  http://mathworld.wolfram.com/SpherePointPicking.html
    G4double thetaTMP=CLHEP::pi/2;
    if(!TransversalyUnpolarisedMuonBeam) thetaTMP = acos(2. * G4UniformRand()-1);
    G4double phiTMP = 2. * CLHEP::pi * G4UniformRand();

    xpolaris = std::sin(thetaTMP) * std::cos(phiTMP);;
    ypolaris = std::sin(thetaTMP) * std::sin(phiTMP);
    zpolaris = std::cos(thetaTMP);
  }
  else {
    if (G4UniformRand()>((1.-polarisFraction)/2.)) {
      xpolaris = xPolarisIni;  ypolaris = yPolarisIni;  zpolaris = zPolarisIni;
      //      G4cout<<"spin up"<<G4endl;
    }
    else {
      xpolaris = -xPolarisIni;  ypolaris = -yPolarisIni;  zpolaris = -zPolarisIni;
      //      G4cout<<"spin down"<<G4endl;
    }
  }
  
  particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
  G4double particle_mass = particleGun->GetParticleDefinition()->GetPDGMass();
  G4double particleEnergy = std::sqrt(p*p+particle_mass*particle_mass)-particle_mass;
  particleGun->SetParticleEnergy(particleEnergy);
  particleGun->SetParticleTime(ParticleTime);  //P.B. 13 May 2009
  particleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
  particleGun->SetParticlePolarization(G4ThreeVector(xpolaris,ypolaris,zpolaris));
  particleGun->GeneratePrimaryVertex(anEvent);

  //  G4cout<<"musrPrimaryGeneratorAction: Parameters:"<<G4endl;
  //  G4cout<<"     x0,y0,z0="<<x0/mm<<","<<y0/mm<<","<<z0/mm<<"   Sigma="<<xSigma/mm<<","<<ySigma/mm<<","<<zSigma/mm<<G4endl;
  //  G4cout<<"     rMaxAllowed="<<rMaxAllowed/mm<<"   zMaxAllowed="<<zMaxAllowed/mm<<"   zMinAllowed="<<zMinAllowed/mm<<G4endl;
  //  G4cout<<"     p0="<<p0/MeV<<"   pSigma="<<pSigma/MeV
  //	<<"   pMinAllowed="<<pMinAllowed/MeV<<"  pMaxAllowed=<<"<<pMaxAllowed/MeV<<G4endl;
  //  G4cout<<"     angle0="<<xangle0/deg<<","<<yangle0/deg<<",nic"
  //	<<"   Sigma="<<xangleSigma/deg<<","<<yangleSigma/deg<<",nic"<<G4endl;
  //  G4cout<<"     pitch="<<pitch/deg<<G4endl;
  //
  //  G4cout<<"musrPrimaryGeneratorAction: Generated muon:"<<G4endl;
  //  G4cout<<"     x,y,z="<<x/mm<<","<<y/mm<<","<<z/mm<<"      angle="<<xangle/deg<<","<< yangle/deg<<",nic"<<G4endl;
  //  G4cout<<"     p="<<px/MeV<<","<<py/MeV<<","<<pz/MeV<<"    E="<< (particleGun->GetParticleEnergy())/MeV<<G4endl;
  //  G4cout<<"     polarisation="<<xpolaris<<","<<ypolaris<<","<<zpolaris<<G4endl;

  


  // if requested by "/gun/decaytimelimits", set the decay time of the muon such that it is within
  // the required time window.  Otherwise the decay time is set internally by Geant.
  if (muonDecayTimeMax>0.) {
    //    G4cout<<"muonDecayTimeMin="<<muonDecayTimeMin/ns<<" ns ,  muonDecayTimeMax="<<muonDecayTimeMax/ns
    //	  <<" ns ,   muonMeanLife="<<muonMeanLife/ns<<" ns."<<G4endl;
    // find the primary muon
    G4PrimaryParticle* generatedMuon = anEvent->GetPrimaryVertex(0)->GetPrimary(0);
    //    G4double decayLowerLimit = 1-exp(-muonDecayTimeMin/muonMeanLife);
    //    G4double decayUpperLimit = 1-exp(-muonDecayTimeMax/muonMeanLife);
    //    G4double randomVal       = G4UniformRand()*(decayUpperLimit-decayLowerLimit) + decayLowerLimit;
    //    G4double decaytime       = -muonMeanLife*log(1-randomVal);
    //
    //  The following code is numerically more stable compared to the commented lines above:
    G4double decaytime;
    if (muonDecayTimeMin==muonDecayTimeMax) {decaytime=muonDecayTimeMin;}
    else {
      G4double expMin    = exp(-muonDecayTimeMin/muonMeanLife);
      G4double expMax    = exp(-muonDecayTimeMax/muonMeanLife);
      decaytime = -muonMeanLife * log(G4UniformRand()*(expMax-expMin)+expMin);
    }
    //    G4cout<<"decaytime="<<decaytime/ns<<"ns."<< G4endl;
    generatedMuon->SetProperTime(decaytime);
  }

  // Set the variable "timeToNextEvent", which is the time difference between this event and the next one
  // at a continuous muon beam.
  G4double timeIntervalBetweenTwoEvents = meanArrivalTime * CLHEP::RandExponential::shoot(1);

  // Save variables into ROOT output file:
  myRootOutput->SetInitialMuonParameters(x,y,z,px,py,pz,xpolaris,ypolaris,zpolaris,ParticleTime);
  myRootOutput->StoreGeantParameter(7,float(numberOfGeneratedEvents));
  myRootOutput->SetTimeToNextEvent(timeIntervalBetweenTwoEvents);
  if (boolPrintInfoAboutGeneratedParticles) {
    G4cout<<"musrPrimaryGeneratorAction::GeneratePrimaries:  x="<<x<<", y="<<y<<", z="<<z<<G4endl;
    G4cout<<"      px="<<px<<", py="<<py<<", pz="<<pz<<", xpolaris="<<xpolaris<<", ypolaris="<<ypolaris<<", zpolaris="<<zpolaris<<G4endl;
    G4cout<<"      numberOfGeneratedEvents="<<numberOfGeneratedEvents<<G4endl;
    G4cout<<"      ------------------------------------"<<G4endl;
  }
}



//===============================================================================
void musrPrimaryGeneratorAction::SetInitialMuonPolariz(G4ThreeVector vIniPol)
{
  G4double magnitude=vIniPol.mag();
  if(magnitude<0.00000001) {
    G4cout<< "Unpolarised initial muons"<<G4endl;
    UnpolarisedMuonBeam=true;
    if ((magnitude<0.0000000085)&&(magnitude>0.0000000075)) {
      G4cout<< "Transversaly unpolarised initial muons"<<G4endl;
      TransversalyUnpolarisedMuonBeam=true;
    }
  }
  else {
    xPolarisIni=vIniPol(0)/magnitude;
    yPolarisIni=vIniPol(1)/magnitude;
    zPolarisIni=vIniPol(2)/magnitude;
    G4cout<< "Initial Muon Polarisation set to ("<<xPolarisIni<<","<<yPolarisIni<<","<<zPolarisIni<<")"<<G4endl;
  }
}

//===============================================================================
void musrPrimaryGeneratorAction::SetBeamDirection(G4ThreeVector vIniDir)
{
  G4double magnitude=vIniDir.mag();
  xDirection=vIniDir(0)/magnitude;
  yDirection=vIniDir(1)/magnitude;
  zDirection=vIniDir(2)/magnitude;
  G4cout<< "Initial Beam Direction set to ("<<xDirection<<","<<yDirection<<","<<zDirection<<")"<<G4endl;
}

//===============================================================================
void musrPrimaryGeneratorAction::SetMuonDecayTimeLimits(G4ThreeVector decayTimeLimits) {
  muonDecayTimeMin = decayTimeLimits[0];
  muonDecayTimeMax = decayTimeLimits[1];
  muonMeanLife     = decayTimeLimits[2];
  // store the muon decay time parameters to the Root output
  musrRootOutput* myRootOutput = musrRootOutput::GetRootInstance();
  myRootOutput->StoreGeantParameter(2,muonDecayTimeMin/CLHEP::microsecond);
  myRootOutput->StoreGeantParameter(3,muonDecayTimeMax/CLHEP::microsecond); 
  myRootOutput->StoreGeantParameter(4,muonMeanLife/CLHEP::microsecond); 
}

//===============================================================================
void musrPrimaryGeneratorAction::SetTurtleInput(G4String turtleFileName) {
  takeMuonsFromTurtleFile = true;
  fTurtleFile = fopen(turtleFileName.c_str(),"r");
  if (fTurtleFile==NULL) {
    G4cout << "E R R O R :    Failed to open TURTLE input file \"" << turtleFileName 
	   <<"\"."<< G4endl;
    G4cout << "S T O P    F O R C E D" << G4endl;
    exit(1);
  }
  else {G4cout << "Turtle input file  \"" << turtleFileName <<"\" opened."<< G4endl;}
}

//===============================================================================
void musrPrimaryGeneratorAction::SetTurtleInputFileToEventNo(G4int lineNumberOfTurtleFile) {
  if (fTurtleFile==NULL) {
    G4cout << "musrPrimaryGeneratorAction::SetTurtleInputFileToEventNo:"
	   <<"  TURTLE input file not found - line number can not be set."<<G4endl;
  }
  else {
    char  line[501];
    for (Int_t i=0; i<lineNumberOfTurtleFile; i++) {
      if (feof(fTurtleFile)) rewind(fTurtleFile);
      fgets(line,500,fTurtleFile);
    }
    G4cout << "musrPrimaryGeneratorAction::SetTurtleInputFileToEventNo: Turtle input file will start at line no.:"
	   << lineNumberOfTurtleFile <<G4endl;
  }
}
//===============================================================================
void musrPrimaryGeneratorAction::SetOrReadTheRandomNumberSeeds(G4Event* anEvent) {
  G4int eventID = anEvent->GetEventID();
  if (eventID == nRndmEventToSaveSeeds) {
    G4cout<<"musrPrimaryGeneratorAction::SetOrReadTheRandomNumberSeeds:  S A V I N G   R A N D O M    N O.   S E E D S"<<G4endl;
    G4cout<<"                                                                     (for even nr. "<<eventID<<")"<<G4endl;
    G4RunManager::GetRunManager()->rndmSaveThisEvent();
    boolPrintInfoAboutGeneratedParticles = true;
  }
  if (eventID == 0) {
    if (setRandomNrSeedFromFile_RNDM) {
      G4cout<<"musrPrimaryGeneratorAction::SetOrReadTheRandomNumberSeeds: Restoring random number seeds from file kamil.rndm"<<G4endl;
      G4RunManager::GetRunManager()->RestoreRandomNumberStatus("kamil.rndm");
      boolPrintInfoAboutGeneratedParticles = true;
    }
  }  
  if (setRandomNrSeedFromFile) {
    //    //    G4cout<<"RandomNrInitialisers.size()="<<RandomNrInitialisers->size()<<G4endl;
    //    if (eventID < int(pointerToSeedVector->size())) {
    //    //    if (eventID < lastEventID_in_pointerToSeedVector)
    //      G4cout <<"musrEventAction.cc: seed will be set to="<< pointerToSeedVector->at(eventID)<<G4endl;
    //      anEvent -> SetEventID(pointerToSeedVector->at(eventID));
    //      //      CLHEP::HepRandom::setTheSeed(pointerToSeedVector->at(eventID));
    //      CLHEP::HepRandom::setTheSeed(eventID);
    //      CLHEP::RandGauss::setFlag(false);
    //      boolPrintInfoAboutGeneratedParticles = true;
    //    }
    if (!pointerToSeedVector->empty()) {
      eventID = pointerToSeedVector->back();
      pointerToSeedVector->pop_back();
    }
    else {
      eventID = ++lastEventID_in_pointerToSeedVector;
    }
    anEvent -> SetEventID(eventID);
    CLHEP::HepRandom::setTheSeed(eventID);
    CLHEP::RandGauss::setFlag(false);
    //    G4cout <<"musrPrimaryGeneratorAction.cc: seed will be set to="<< eventID<<G4endl;
  }
  else if (setRandomNrSeedAccordingEventNr) {
    //    long seeds[2];
    //    seeds[0] = (long) 234567890+thisEventNr*117;
    //    seeds[1] = (long) 333222111+thisEventNr*173;
    //
    //    //    seeds[1] = (long) (evt->GetEventID());
    //    //    seeds[0] = (long) 123456789;    // This leads to a gap in the decay time histogram fro N=100000 events
    //    //    seeds[1] = (long) 333222111+thisEventNr;     // ----------------------------||------------------------------------
    //    thisEventNr++;
    //    CLHEP::HepRandom::setTheSeeds(seeds);
    //    //    G4cout << "seed1: " << seeds[0] << "; seed2: " << seeds[1] << G4endl;
    //
    //    G4cout <<"      thisEventNr="<<thisEventNr;
    CLHEP::HepRandom::setTheSeed(eventID);
    //    G4cout <<"     getTheSeed="<<CLHEP::HepRandom::getTheSeed()<< G4endl; 
    CLHEP::RandGauss::setFlag(false);
  }
}

void musrPrimaryGeneratorAction::SetKEnergy(G4double val)  {
  //G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //G4double mu_mass = particleTable->FindParticle("mu+")->GetPDGMass();
  G4double particle_mass = particleGun->GetParticleDefinition()->GetPDGMass();
  p0=std::sqrt(val*val + 2*particle_mass*val);
/*   G4cout<<"musrPrimaryGeneratorAction::SetKEnergy: Particle mass " <<particle_mass<<" particle kinetic energy of "
  	<<val<<" MeV requested ==> initial particle momentum set to "<<p0<<" MeV/c"<<G4endl;*/
}

//===============================================================================
void musrPrimaryGeneratorAction::SetPrimaryParticule(G4String particleName) {
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particule= particleTable->FindParticle(particleName);
  if (particule!=NULL) {
    G4cout<<"   musrPrimaryGeneratorAction::SetPrimaryParticleMuMinus():  USING "
	  <<particleName<<" AS PRIMARY PARTICLE!"<<G4endl; 
    particleGun->SetParticleDefinition(particule);
  }
  else {
    G4cout<<"\n   musrPrimaryGeneratorAction::SetPrimaryParticleMuMinus():  Particle "
	  <<particleName<<" required as primary particle, but not found !!!";
    G4cout<<"S T O P     F O R C E D" << G4endl;
    exit(1);
  }
}

//===============================================================================
void musrPrimaryGeneratorAction::swapTheAxisInTurtle(float& x_x, float& x_xprime ,float& y_y, float& y_yprime) {
  //  G4cout<<"turtleInterpretAxes = "<<turtleInterpretAxes<<G4endl;
  if (turtleInterpretAxes=="-xy") {
    x_x = -x_x; x_xprime = -x_xprime;
  }
  else if (turtleInterpretAxes=="x-y") {
    y_y = -y_y; y_yprime = -y_yprime;
  }
  else if (turtleInterpretAxes=="-x-y") {
    x_x = -x_x; x_xprime = -x_xprime; y_y = -y_y; y_yprime = -y_yprime;
  }
  else if (turtleInterpretAxes=="yx") {
    float tmpX = x_x; float tmpXprime = x_xprime;
    x_x = y_y; x_xprime = y_yprime; y_y = tmpX; y_yprime = tmpXprime;
  }
  else if (turtleInterpretAxes=="-yx") {
    float tmpX = x_x; float tmpXprime = x_xprime;
    x_x = y_y; x_xprime = y_yprime; y_y = -tmpX; y_yprime = -tmpXprime;
  }
  else if (turtleInterpretAxes=="y-x") {
    float tmpX = x_x; float tmpXprime = x_xprime;
    x_x = -y_y; x_xprime = -y_yprime; y_y = tmpX; y_yprime = tmpXprime;
  }
  else if (turtleInterpretAxes=="-y-x") {
    float tmpX = x_x; float tmpXprime = x_xprime;
    x_x = -y_y; x_xprime = -y_yprime; y_y = -tmpX; y_yprime = -tmpXprime;
  }
  else {
    G4cout<<"musrPrimaryGeneratorAction::swapTheAxisInTurtle: Not known how to inpterpret turtleInterpretAxes="<<turtleInterpretAxes<<G4endl;
    G4cout<<"S T O P     F O R C E D" << G4endl;
    exit(1);
  }
}

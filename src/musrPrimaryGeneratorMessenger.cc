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

#include "musrPrimaryGeneratorMessenger.hh"
#include "musrPrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

musrPrimaryGeneratorMessenger::musrPrimaryGeneratorMessenger(musrPrimaryGeneratorAction* musrGun)
:musrAction(musrGun)
{ 
  setPrimaryParticleCmd = new G4UIcmdWithAString("/gun/primaryparticle",this);
  setPrimaryParticleCmd -> SetGuidance("Set the primary particle (default is mu+).");
  setPrimaryParticleCmd ->AvailableForStates(G4State_PreInit,G4State_Idle);

  setvertexCmd = new G4UIcmdWith3VectorAndUnit("/gun/vertex",this);
  setvertexCmd->SetGuidance(" Set x0, y0, z0 of the generated muons (with unit)");
  setvertexCmd->SetParameterName("mes_x0","mes_y0","mes_z0",true,true);
  setvertexCmd->SetDefaultUnit("mm");

  setvertexSigmaCmd = new G4UIcmdWith3VectorAndUnit("/gun/vertexsigma",this);
  setvertexSigmaCmd->SetGuidance(" Set xSigma, ySigma, ySigma of the generated muons (with unit)");
  setvertexSigmaCmd->SetParameterName("mes_xSigma","mes_ySigma","mes_zSigma",true,true);
  setvertexSigmaCmd->SetDefaultUnit("mm");

  setvertexBoundaryCmd = new G4UIcmdWith3VectorAndUnit("/gun/vertexboundary",this);
  setvertexBoundaryCmd->SetGuidance(" Set maximum allowed radius, zmin, zmax of the generated vertex (with unit)");
  setvertexBoundaryCmd->SetParameterName("mes_rMaxAllowed","mes_zMinAllowed","mes_zMaxAllowed",true,true);
  setvertexBoundaryCmd->SetDefaultUnit("mm");

  setvertexRelativeRCmd = new G4UIcmdWithADoubleAndUnit("/gun/vertexrelativer",this);
  setvertexRelativeRCmd->SetGuidance(" Set maximum allowed radius of the beam relative to x0 and y0,");
  setvertexRelativeRCmd->SetGuidance("    i.e. relative to the centre of the beam (with unit)");
  setvertexRelativeRCmd->SetParameterName("mes_relativeRMaxAllowed",true);
  setvertexRelativeRCmd->SetDefaultUnit("mm"); 

  setMeanArrivalTimeCmd = new G4UIcmdWithADoubleAndUnit("/gun/meanarrivaltime",this);
  setMeanArrivalTimeCmd->SetGuidance(" Set mean arrival time difference between two subsequent muons (at continuos beam)");
  setMeanArrivalTimeCmd->SetParameterName("mes_meanarrivaltime",true);
  setMeanArrivalTimeCmd->SetDefaultUnit("ns");

  setStarttimeCmd = new G4UIcmdWithADoubleAndUnit("/gun/starttime",this);                //P.B. 13 May 2009
  setStarttimeCmd->SetGuidance(" Set start time t of the generated muons (with unit)");  //P.B. 13 May 2009
  setStarttimeCmd->SetParameterName("mes_t0",true);                                      //P.B. 13 May 2009
  setStarttimeCmd->SetDefaultUnit("ns");                                                 //P.B. 13 May 2009

  setStarttimeSigmaCmd = new G4UIcmdWithADoubleAndUnit("/gun/starttimesigma",this);                      //P.B. 13 May 2009
  setStarttimeSigmaCmd->SetGuidance(" Set start time sigma tSigma of the generated muons (with unit)");  //P.B. 13 May 2009
  setStarttimeSigmaCmd->SetParameterName("mes_tSigma",true);                                             //P.B. 13 May 2009
  setStarttimeSigmaCmd->SetDefaultUnit("ns");                                                            //P.B. 13 May 2009

  setboxBoundaryCentreCmd = new G4UIcmdWith3VectorAndUnit("/gun/boxboundarycentre",this);                                   //P.B. 15 Dec 2009
  setboxBoundaryCentreCmd->SetGuidance(" Set centre point (x,y,z)MaxSource0 (with unit) of the restriction source box");    //P.B. 15 Dec 2009
  setboxBoundaryCentreCmd->SetParameterName("mes_xMaxSource0","mes_yMaxSource0","mes_zMaxSource0",true,true);               //P.B. 15 Dec 2009
  setboxBoundaryCentreCmd->SetDefaultUnit("mm");                                                                            //P.B. 15 Dec 2009

  setboxBoundaryCmd = new G4UIcmdWith3VectorAndUnit("/gun/boxboundary",this);                                               //P.B. 15 Dec 2009
  setboxBoundaryCmd->SetGuidance(" Set maximum deviation from (x,y,z)MaxSource0 (with unit)");                              //P.B. 15 Dec 2009
  setboxBoundaryCmd->SetParameterName("mes_xMaxSource","mes_yMaxSource","mes_zMaxSource",true,true);                        //P.B. 15 Dec 2009
  setboxBoundaryCmd->SetDefaultUnit("mm");                                                                                  //P.B. 15 Dec 2009

  setKEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gun/kenergy",this);
  setKEnergyCmd->SetGuidance(" Set kinetic energy of the generated muons (with unit)");
  setKEnergyCmd->SetParameterName("mes_E0",true);
  setKEnergyCmd->SetDefaultUnit("MeV");

  setMomentumCmd = new G4UIcmdWithADoubleAndUnit("/gun/momentum",this);
  setMomentumCmd->SetGuidance(" Set mean momentum of the generated muons (with unit)");
  setMomentumCmd->SetParameterName("mes_p0",true);
  setMomentumCmd->SetDefaultUnit("MeV");

  setMomentumSmearingCmd = new G4UIcmdWithADoubleAndUnit("/gun/momentumsmearing",this);
  setMomentumSmearingCmd->SetGuidance(" Set sigma of the momentum of the generated muons (with unit)");
  setMomentumSmearingCmd->SetParameterName("mes_pSigma",true);
  setMomentumSmearingCmd->SetDefaultUnit("MeV");

  setMomentumBoundaryCmd = new G4UIcmdWith3VectorAndUnit("/gun/momentumboundary",this);
  setMomentumBoundaryCmd->SetGuidance(" Set minimum and maximum momentum allowed (with unit, z component ignored)");
  setMomentumBoundaryCmd->SetParameterName("mes_pMinAllowed","mes_pMaxAllowed","mes_dummy",true,true);
  setMomentumBoundaryCmd->SetDefaultUnit("MeV");

  setTiltAngleCmd = new G4UIcmdWith3VectorAndUnit("/gun/tilt",this);
  setTiltAngleCmd->SetGuidance(" Set tilt angle of the generated muons (with unit, z component ignored)");
  setTiltAngleCmd->SetParameterName("mes_xangle","mes_yangle","dummy",true,true);
  setTiltAngleCmd->SetDefaultUnit("deg");

  setSigmaTiltAngleCmd = new G4UIcmdWith3VectorAndUnit("/gun/tiltsigma",this);
  setSigmaTiltAngleCmd->SetGuidance(" Set sigma of the tilt angle (with unit, z component ignored)");
  setSigmaTiltAngleCmd->SetParameterName("mes_xangleSigma","mes_yangleSigma","dummy",true,true);
  setSigmaTiltAngleCmd->SetDefaultUnit("deg");

  setPitchCmd = new G4UIcmdWithADoubleAndUnit("/gun/pitch",this);
  setPitchCmd->SetGuidance(" Set pitch angle of the generated muons (with unit)");
  setPitchCmd->SetParameterName("mes_pitch",true);
  setPitchCmd->SetDefaultUnit("deg");

  //  setMuonPolarizCmd = new G4UIcmdWithAnInteger("/gun/muonpolarization",this);
  //  setMuonPolarizCmd->SetGuidance(" Set initial mu polariz: 0=transverse, 1=longitudinal ");
  //  setMuonPolarizCmd->SetParameterName("IniPol",true);
  //  setMuonPolarizCmd->SetDefaultValue(0) ;

  setMuonPolarizCmd = new G4UIcmdWith3Vector("/gun/muonPolarizVector",this);
  setMuonPolarizCmd->SetGuidance("Set initial mu polarisation as a vector (without units)");
  setMuonPolarizCmd->SetGuidance("    The vector does not have to be normalised to 1");
  setMuonPolarizCmd->SetParameterName("mes_polarisX","mes_polarisY","mes_polarisZ",true,true);

  setDirectionCmd = new G4UIcmdWith3Vector("/gun/direction",this);
  setDirectionCmd->SetGuidance("Set initial mu beam direction as a vector (without units)");
  setDirectionCmd->SetGuidance("    The vector does not have to be normalised to 1");
  setDirectionCmd->SetParameterName("DirectionX","DirectionY","DirectionZ",true,true);

  setMuonPolarizFractionCmd = new G4UIcmdWithADouble("/gun/muonPolarizFraction",this);
  setMuonPolarizFractionCmd->SetGuidance(" Set the fraction of the muon polarisation (in the range of -1 to 1),");
  setMuonPolarizFractionCmd->SetGuidance(" where fraction = (N_up_spin - N_down_spin) / (N_up_spin + N_down_spin)");
  setMuonPolarizFractionCmd->SetParameterName("mes_polarisFraction",true);

  setMuonDecayTimeCmd = new G4UIcmdWith3VectorAndUnit("/gun/decaytimelimits",this);
  setMuonDecayTimeCmd->SetGuidance(" Set minimum and maximum decay time and the muon mean life ");
  setMuonDecayTimeCmd->SetParameterName("decayMin","decayMax","decayTime",true,true);
  setMuonDecayTimeCmd->SetDefaultUnit("ns");

  setTurtleCmd = new G4UIcmdWithAString("/gun/turtlefilename",this);
  setTurtleCmd->SetGuidance("Set the filename of the TURTLE input file.");
  setTurtleCmd->SetGuidance("If this varialble is set, TURTLE input will be used for initial muons");
  setTurtleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  setTurtleZ0Cmd = new G4UIcmdWithADoubleAndUnit("/gun/turtleZ0position",this);
  setTurtleZ0Cmd->SetGuidance("Set the z0, with which the TURTLE input file has been generated.");
  setTurtleZ0Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  setTurtleZ0Cmd->SetParameterName("mes_z0Turtle",true);
  setTurtleZ0Cmd->SetDefaultUnit("mm");

  setTurtleInterpretAxesCmd = new G4UIcmdWithAString("/gun/turtleInterpretAxes",this);
  setTurtleInterpretAxesCmd->SetGuidance("Specify whether the x and y axes of the position in TURTLE should be interpretted differently.");
  setTurtleInterpretAxesCmd->SetGuidance("The following options are supported:  x-y, -xy, -x-y, yx, y-x, -yx, -y-x .");
  setTurtleInterpretAxesCmd->SetGuidance("Example: the option y-x means that first four coordinates in the TURTLE input file");
  setTurtleInterpretAxesCmd->SetGuidance("         are interpreded as  y, yprime, -x, -xprime.");
  setTurtleInterpretAxesCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  setTurtleMomentumBite = new G4UIcmdWith3Vector("/gun/turtleMomentumBite",this);
  setTurtleMomentumBite->SetGuidance(" Modify smearing of the turtle momentum bite.  The first value is the mean momentum in MeV/c,");
  setTurtleMomentumBite->SetGuidance("   the second value is the smearing factor in per cent, by which the momentum bite,");
  setTurtleMomentumBite->SetGuidance("   will be increased/decreased around the mean momemtum.  100 per cent correspond to no");
  setTurtleMomentumBite->SetGuidance("   change, 0 per cent will create monoenergetic beam, 200 per cent will create a two times");
  setTurtleMomentumBite->SetGuidance("   broader beam.  The third parameter is dummy.");
  setTurtleMomentumBite->SetParameterName("mes_turtleMomentumP0","mes_turtleSmearingFactor","mes_turtleSmearingDummy",true,true);

  setTurtleMomentumScalingFactor = new G4UIcmdWithADouble("/gun/turtleMomentumScalingFactor",this);
  setTurtleMomentumScalingFactor->SetGuidance(" Modify the turtle momentum by this scaling factor - e.g. if the user wants to shift");
  setTurtleMomentumScalingFactor->SetGuidance("   the mean momentum from 28 MeV/c to 14 MeV/c, this scaling factor has to be set to 0.5 .");
  setTurtleMomentumScalingFactor->SetGuidance("   Note that also the momentum bite will be reduced by the same scaling factor.");
  setTurtleMomentumScalingFactor->SetParameterName("mes_turtleMomentumScalingFactor",true);


  setTurtleEventNrCmd = new G4UIcmdWithAnInteger("/gun/turtleFirstEventNr",this);
  setTurtleEventNrCmd->SetGuidance("Set the line number that should be taken as the first event from the turtle input file.");
  setTurtleEventNrCmd->SetParameterName("mes_turtleFirstEvent",true);
  setTurtleEventNrCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  setTurtleEventNrCmd->SetDefaultValue(0) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

musrPrimaryGeneratorMessenger::~musrPrimaryGeneratorMessenger()
{
  delete setPrimaryParticleCmd;
  delete setvertexCmd;
  delete setvertexSigmaCmd;
  delete setvertexBoundaryCmd;
  delete setvertexRelativeRCmd;
  delete setMeanArrivalTimeCmd;
  delete setStarttimeCmd;         //P.B. 13 May 2009
  delete setStarttimeSigmaCmd;    //P.B. 13 May 2009
  delete setboxBoundaryCentreCmd; //P.B. 15 Dec 2009
  delete setboxBoundaryCmd;       //P.B. 15 Dec 2009
  delete setKEnergyCmd;
  delete setMomentumCmd;
  delete setMomentumSmearingCmd;
  delete setMomentumBoundaryCmd;
  delete setTiltAngleCmd;
  delete setSigmaTiltAngleCmd;
  delete setPitchCmd;
  delete setMuonPolarizCmd;
  delete setDirectionCmd;
  delete setMuonPolarizFractionCmd;
  delete setMuonDecayTimeCmd;
  delete setTurtleCmd;
  delete setTurtleZ0Cmd;
  delete setTurtleInterpretAxesCmd;
  delete setTurtleMomentumBite;
  delete setTurtleMomentumScalingFactor;
  delete setTurtleEventNrCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void musrPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == setPrimaryParticleCmd)
    { musrAction->SetPrimaryParticule(newValue); }
  if( command == setvertexCmd)
    { musrAction->Setvertex(setvertexCmd->GetNew3VectorValue(newValue));}
  if( command == setvertexSigmaCmd)
    { musrAction->SetvertexSigma(setvertexSigmaCmd->GetNew3VectorValue(newValue));}
  if( command == setvertexBoundaryCmd)
    { musrAction->SetvertexBoundary(setvertexBoundaryCmd->GetNew3VectorValue(newValue));}
  if( command == setvertexRelativeRCmd)
    { musrAction->SetvertexRelativeR(setvertexRelativeRCmd->GetNewDoubleValue(newValue));}
  if( command == setMeanArrivalTimeCmd)
    { musrAction->SetMeanArrivalTime(setMeanArrivalTimeCmd->GetNewDoubleValue(newValue));}
  if( command == setStarttimeCmd)                                                       //P.B. 13 May 2009
    { musrAction->SetMuonTime(setStarttimeCmd->GetNewDoubleValue(newValue));}           //P.B. 13 May 2009
  if( command == setStarttimeSigmaCmd)                                                  //P.B. 13 May 2009
    { musrAction->SetMuonTimeSigma(setStarttimeSigmaCmd->GetNewDoubleValue(newValue));} //P.B. 13 May 2009
  if( command == setboxBoundaryCentreCmd)                                               //P.B. 15 Dec 2009
    { musrAction->SetboxBoundaryCentre(setboxBoundaryCentreCmd->GetNew3VectorValue(newValue));} //P.B. 15 Dec 2009
  if( command == setboxBoundaryCmd)                                                      //P.B. 15 Dec 2009
    { musrAction->SetboxBoundary(setboxBoundaryCmd->GetNew3VectorValue(newValue));}      //P.B. 15 Dec 2009
  if( command == setKEnergyCmd)
   { musrAction->SetKEnergy(setKEnergyCmd->GetNewDoubleValue(newValue));}
  if( command == setMomentumCmd)
   { musrAction->SetMomentum(setMomentumCmd->GetNewDoubleValue(newValue));}
  if( command == setMomentumSmearingCmd)
    { musrAction->SetMomentumSmearing(setMomentumSmearingCmd->GetNewDoubleValue(newValue));}
  if( command == setMomentumBoundaryCmd)
   { musrAction->SetMomentumBoundary(setMomentumBoundaryCmd->GetNew3VectorValue(newValue));}
  if( command == setTiltAngleCmd)
   { musrAction->SetTilt(setTiltAngleCmd->GetNew3VectorValue(newValue));}
  if( command == setSigmaTiltAngleCmd)
   { musrAction->SetSigmaTilt(setSigmaTiltAngleCmd->GetNew3VectorValue(newValue));}
  if( command == setPitchCmd)
    { musrAction->SetPitch(setPitchCmd->GetNewDoubleValue(newValue));}
  if( command == setMuonPolarizCmd)
    { musrAction->SetInitialMuonPolariz(setMuonPolarizCmd->GetNew3VectorValue(newValue));}
  if( command == setDirectionCmd)
    { musrAction->SetBeamDirection(setDirectionCmd->GetNew3VectorValue(newValue));}
  if( command == setMuonPolarizFractionCmd)
    { musrAction->SetInitialPolarizFraction(setMuonPolarizFractionCmd->GetNewDoubleValue(newValue));}
  if( command == setMuonDecayTimeCmd)
    { musrAction->SetMuonDecayTimeLimits(setMuonDecayTimeCmd->GetNew3VectorValue(newValue)); }
  if( command == setTurtleCmd)
    { musrAction->SetTurtleInput(newValue); }
  if( command == setTurtleZ0Cmd)
    { musrAction->SetTurtleZ0(setTurtleZ0Cmd->GetNewDoubleValue(newValue)); }
  if( command == setTurtleInterpretAxesCmd)
    { musrAction->SetTurtleInterpretAxes(newValue);}
  if( command == setTurtleMomentumBite)
    { musrAction-> SetTurtleMomentumBite(setTurtleMomentumBite->GetNew3VectorValue(newValue)); }
  if( command == setTurtleMomentumScalingFactor)
    { musrAction-> SetTurtleMomentumScalingFactor(setTurtleMomentumScalingFactor->GetNewDoubleValue(newValue));}
  if( command == setTurtleEventNrCmd)
    { musrAction->SetTurtleInputFileToEventNo(setTurtleEventNrCmd->GetNewIntValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


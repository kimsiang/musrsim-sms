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

#include "musrMuonium.hh"
#include "G4ParticleTable.hh"

#include "MuDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                          MUONIUM                               ###
// ######################################################################
musrMuonium* musrMuonium::theInstance = 0;

musrMuonium* musrMuonium::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "Mu";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance ==0)
  {
  // create particle
  //
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding
  anInstance = new G4ParticleDefinition(
                 name,   0.1056584*CLHEP::GeV, 2.99591e-16*CLHEP::MeV,   0.*CLHEP::eplus, 
		    1,               0,             0,          
		    0,               0,             0,             
	     "lepton",              -1,             0,         -1313,
		false,      2197.03*CLHEP::ns,          NULL,
             false,           "mu"
              );
   // Bohr magnetron of Muonium - T. Shiroka
   // The magnetic moment of Mu is the sum of those of mu+ and e- with
   // the respective gyromagnetic ratio anomalies as coefficients
   
   G4double muBmu =  0.5*CLHEP::eplus*CLHEP::hbar_Planck/(0.10565840*CLHEP::GeV/CLHEP::c_squared);
   G4double muBel = -0.5*CLHEP::eplus*CLHEP::hbar_Planck/(0.51099906*CLHEP::MeV/CLHEP::c_squared);
   G4double muB   =  1.0011659208*muBmu + 1.0011596521859*muBel;
   
   anInstance->SetPDGMagneticMoment( muB );

  //create Decay Table 
  G4DecayTable* table = new G4DecayTable();
  // create a decay channel
  G4VDecayChannel* mode = new MuDecayChannel("Mu",1.00);
  table->Insert(mode);
  anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<musrMuonium*>(anInstance);
  return theInstance;
}

musrMuonium*  musrMuonium::MuoniumDefinition()
{
  return Definition();
}

musrMuonium*  musrMuonium::Muonium()
{
  return Definition();
}


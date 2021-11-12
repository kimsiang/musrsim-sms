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
// $Id: G4MuonDecayChannel.hh,v 1.6 2006/06/29 19:23:35 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      30 May 1997 H.Kurashige
//
//      2008-05 Toni Shiroka: Extention of the code for the muonium decay
// ------------------------------------------------------------
#ifndef MuDecayChannel_h
#define MuDecayChannel_h 1

#include <CLHEP/Units/PhysicalConstants.h>
#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"

class MuDecayChannel :public G4VDecayChannel
{
  // Class Decription
  //  This class describes muon decay kinemtics.
  //  This version neglects muon polarization  
  //              assumes the pure V-A coupling
  //              gives incorrect energy spectrum for neutrinos
  //

  public:  // With Description
    //Constructors 
      MuDecayChannel(const G4String& theParentName,
			   G4double  theBR);
    //  Destructor
      virtual ~MuDecayChannel();

  public:  // With Description
     virtual G4DecayProducts *DecayIt(G4double);     
  
};


#endif

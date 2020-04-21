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
// $Id$
//
/// \file MonoChromEventAction.cc
/// \brief Implementation of the MonoChromEventAction class

#include "MonoChromEventAction.hh"
#include "MonoChromRunAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonoChromEventAction::MonoChromEventAction(MonoChromRunAction* runAction)
 : G4UserEventAction(), 
   fRunAction(runAction),
   fCollectionIDphotodiode(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonoChromEventAction::~MonoChromEventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MonoChromEventAction::BeginOfEventAction(const G4Event* /*evt*/)
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MonoChromEventAction::EndOfEventAction(const G4Event* evt )
{
  //Hits collections
  //  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  
  // Get hits collections IDs
  if (fCollectionIDphotodiode < 0) {
    G4SDManager* SDMan = G4SDManager::GetSDMpointer();  
    fCollectionIDphotodiode = SDMan->GetCollectionID("photodiode/edep");
  }
  
  //Energy in crystals : identify 'good events'
  //
  const G4double eThreshold = 0.001*eV;
  G4int nbOfFired = 0;
  
  G4THitsMap<G4double>* evtMap = 
    (G4THitsMap<G4double>*)(HCE->GetHC(fCollectionIDphotodiode));
  
  G4double edep = 0.0;
  std::map<G4int,G4double*>::iterator itr;
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    //G4int copyNb  = (itr->first);
    edep = *(itr->second);
    
    if (edep > eThreshold) nbOfFired++;
    
    //G4cout << "\n  photodiode " << edep/eV << " eV " << G4endl;
  }  
  
  if(edep/eV > 0.1){
    fRunAction->AddEnergy(edep);
    fRunAction->CountEvent();
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

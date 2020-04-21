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
/// \file MonoChrom/include/MonoChromDetectorConstruction.hh
/// \brief Definition of the MonoChromDetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MonoChromDetectorConstruction_h
#define MonoChromDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MonoChromDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  MonoChromDetectorConstruction();
  virtual ~MonoChromDetectorConstruction();
  
public:
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
private:

  void SetMaterials();
  G4Material* air = nullptr;
  G4Material* polyEth = nullptr;
  G4Material* tyvek = nullptr;

  static const G4int nIndices = 2;
  G4double energies[nIndices];

  G4double air_RIs[nIndices];  
  G4double pe_RIs[nIndices];  

  void SetVolumes();
  G4VPhysicalVolume* expHall_phys = nullptr;
  G4VPhysicalVolume* darkBox_phys = nullptr;
  G4VPhysicalVolume* reflectSheet_phys = nullptr;


  void SetDimensions();
  G4double fExpHall_x;
  G4double fExpHall_y;
  G4double fExpHall_z;
  
  G4double fDarkBox_x;
  G4double fDarkBox_y;
  G4double fDarkBox_z;
  
  G4double fReflectSheet_x;
  G4double fReflectSheet_y;
  G4double fReflectSheet_z;

  G4double fDetRMin;
  G4double fDetRMax;
  G4double fDetDz;
  G4double fDetSPhi;
  G4double fDetDPhi;

  void SetSurfaces();
  
  


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*MonoChromDetectorConstruction_h*/

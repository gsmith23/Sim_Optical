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
/// \file MonoChrom/src/MonoChromDetectorConstruction.cc
/// \brief Implementation of the MonoChromDetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MonoChromDetectorConstruction.hh"

#include "G4NistManager.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonoChromDetectorConstruction::MonoChromDetectorConstruction()
 : G4VUserDetectorConstruction()
{
  // half lengths for boxes
  fExpHall_x = fExpHall_y = fExpHall_z = 0.55*m;
  fDarkBox_x = fDarkBox_y = fDarkBox_z = 0.5*m;

  fTyvek_x   = fTyvek_y   = 0.25*m;
  // http://protectiontechnologies.dupont.com/tyvek-graphics-na-product-selector
  // 1082D 257 um thick, 105 gsm  
  fTyvek_z   = 0.01375*cm; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonoChromDetectorConstruction::~MonoChromDetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* MonoChromDetectorConstruction::Construct()
{

  // --------------------------------------
  // ------------- G4Material -------------
  
  G4NistManager* nist = G4NistManager::Instance();
  
//   G4double a, z, density;
//   G4int nelements;
  
  //-------------
  // Air
  //
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

  //  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  //  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  //  G4Material* air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  //  air->AddElement(N, 70.*perCent);
  //  air->AddElement(O, 30.*perCent);
  //-------------

  //-------------
  // Polyethelene
  //  /material/nist/listMaterials
  
  // 1082D 257 um thick, 105 gsm
  // 275 / 105 = 0.382
  
  // http://protectiontechnologies.dupont.com/tyvek-graphics-na-product-selector
  //  G4double tyvek_thickness = 0.275 * m;
  //  G4double tyvek_gsm       = 105.0 * g / m2;
  G4double tyvek_density   = 0.382 * g / cm3; 
  
  //G4Material* polyEth = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  
  G4Material* tyvek = 
    nist->BuildMaterialWithNewDensity("TYVEK",
				      "G4_POLYETHYLENE",
				      tyvek_density);
  //-------------
  
  // array elements for optics
  static const G4int nElements = 2;
  G4double photon_Energies[] = {1.7712027,6.1992093}; //[700,200] nm
  
  
  //------------------------------------------------------
  //---------- Volumes G4MaterialPropertiesTables --------
  // 
  /* Water
   
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);
  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);
  
  //
  // ------------ Generate & Add Material Properties Table ------------
  //
  
  // [610,300]
  G4double photon_Energies[] =
    { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
      2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
      2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
      2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
      2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
      3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
      3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
      3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
  
  const G4int nElements = sizeof(photon_Energies)/sizeof(G4double);

  //
  // Water
  //
  G4double refractiveIndex1[] =
    { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
      1.346,  1.3465, 1.347,  1.3475, 1.348,
      1.3485, 1.3492, 1.35,   1.3505, 1.351,
      1.3518, 1.3522, 1.3530, 1.3535, 1.354,
      1.3545, 1.355,  1.3555, 1.356,  1.3568,
      1.3572, 1.358,  1.3585, 1.359,  1.3595,
      1.36,   1.3608};
  
  assert(sizeof(refractiveIndex1) == sizeof(photon_Energies));
  
  G4double absorption[] =
    {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
     15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
     45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
     52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
     30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
     17.500*m, 14.500*m };
  
  assert(sizeof(absorption) == sizeof(photon_Energies));
  
  G4double scintilFast[] =
    { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00 };
  
  assert(sizeof(scintilFast) == sizeof(photon_Energies));
  
  G4double scintilSlow[] =
    { 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
      7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00,
      3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00,
      4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
      7.00, 6.00, 5.00, 4.00 };
  
  assert(sizeof(scintilSlow) == sizeof(photon_Energies));
  
  G4MaterialPropertiesTable* water_Table = new G4MaterialPropertiesTable();
  
  water_Table->AddProperty("RINDEX",       photon_Energies, refractiveIndex1,nElements)
    ->SetSpline(true);
  water_Table->AddProperty("ABSLENGTH",    photon_Energies, absorption,     nElements)
    ->SetSpline(true);
  water_Table->AddProperty("FASTCOMPONENT",photon_Energies, scintilFast,     nElements)
    ->SetSpline(true);
  water_Table->AddProperty("SLOWCOMPONENT",photon_Energies, scintilSlow,     nElements)
    ->SetSpline(true);
  
  water_Table->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  water_Table->AddConstProperty("RESOLUTIONSCALE",1.0);
  water_Table->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  water_Table->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  water_Table->AddConstProperty("YIELDRATIO",0.8);
  
  // [790,200] nm
  G4double energy_water[] = {
    1.56962*eV, 1.58974*eV, 1.61039*eV, 1.63157*eV,
    1.65333*eV, 1.67567*eV, 1.69863*eV, 1.72222*eV,
    1.74647*eV, 1.77142*eV, 1.7971 *eV, 1.82352*eV,
    1.85074*eV, 1.87878*eV, 1.90769*eV, 1.93749*eV,
    1.96825*eV, 1.99999*eV, 2.03278*eV, 2.06666*eV,
    2.10169*eV, 2.13793*eV, 2.17543*eV, 2.21428*eV,
    2.25454*eV, 2.29629*eV, 2.33962*eV, 2.38461*eV,
    2.43137*eV, 2.47999*eV, 2.53061*eV, 2.58333*eV,
    2.63829*eV, 2.69565*eV, 2.75555*eV, 2.81817*eV,
    2.88371*eV, 2.95237*eV, 3.02438*eV, 3.09999*eV,
    3.17948*eV, 3.26315*eV, 3.35134*eV, 3.44444*eV,
    3.54285*eV, 3.64705*eV, 3.75757*eV, 3.87499*eV,
    3.99999*eV, 4.13332*eV, 4.27585*eV, 4.42856*eV,
    4.59258*eV, 4.76922*eV, 4.95999*eV, 5.16665*eV,
    5.39129*eV, 5.63635*eV, 5.90475*eV, 6.19998*eV
  };
  
  const G4int numentries_water = sizeof(energy_water)/sizeof(G4double);
  
  //assume 100 times larger than the rayleigh scattering for now.
  G4double mie_water[] = {
    167024.4*m, 158726.7*m, 150742  *m,
    143062.5*m, 135680.2*m, 128587.4*m,
    121776.3*m, 115239.5*m, 108969.5*m,
    102958.8*m, 97200.35*m, 91686.86*m,
    86411.33*m, 81366.79*m, 76546.42*m,
    71943.46*m, 67551.29*m, 63363.36*m,
    59373.25*m, 55574.61*m, 51961.24*m,
    48527.00*m, 45265.87*m, 42171.94*m,
    39239.39*m, 36462.50*m, 33835.68*m,
    31353.41*m, 29010.30*m, 26801.03*m,
    24720.42*m, 22763.36*m, 20924.88*m,
    19200.07*m, 17584.16*m, 16072.45*m,
    14660.38*m, 13343.46*m, 12117.33*m,
    10977.70*m, 9920.416*m, 8941.407*m,
    8036.711*m, 7202.470*m, 6434.927*m,
    5730.429*m, 5085.425*m, 4496.467*m,
    3960.210*m, 3473.413*m, 3032.937*m,
    2635.746*m, 2278.907*m, 1959.588*m,
    1675.064*m, 1422.710*m, 1200.004*m,
    1004.528*m, 833.9666*m, 686.1063*m
  };
  
  assert(sizeof(mie_water) == sizeof(energy_water));
  
  // gforward, gbackward, forward backward ratio
  G4double mie_water_const[3]={0.99,0.99,0.8};
  
  water_Table->AddProperty("MIEHG",energy_water,mie_water,numentries_water)
    ->SetSpline(true);
  water_Table->AddConstProperty("MIEHG_FORWARD",mie_water_const[0]);
  water_Table->AddConstProperty("MIEHG_BACKWARD",mie_water_const[1]);
  water_Table->AddConstProperty("MIEHG_FORWARD_RATIO",mie_water_const[2]);
  
  G4cout << "Water G4MaterialPropertiesTable" << G4endl;
  water_Table->DumpTable();
  
  water->SetMaterialPropertiesTable(water_Table);
  
  // Set the Birks Constant for the Water scintillator
  
  water->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
  
  */
  //----------------------------------------------------
  
  //
  // Air
  //
  //   G4double refractiveIndex2[] =
  //     { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
  //       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
  //       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
  //       1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
  //       1.00, 1.00, 1.00, 1.00 };
  
  G4double air_R_Indices[] = {1.0,1.0};
  G4double pe_R_Indices[]  = {1.5,1.5};
  
  G4MaterialPropertiesTable* air_Table = new G4MaterialPropertiesTable();
  air_Table->AddProperty("RINDEX", photon_Energies, air_R_Indices, nElements);
  
  G4MaterialPropertiesTable* tyvek_Table = new G4MaterialPropertiesTable();
  tyvek_Table->AddProperty("RINDEX", photon_Energies, pe_R_Indices, nElements);
  
  G4cout << "\n Air G4MaterialPropertiesTable" << G4endl;
  air_Table->DumpTable();
  
  G4cout << "\n Tyvek G4MaterialPropertiesTable" << G4endl;
  tyvek_Table->DumpTable();
  
  air->SetMaterialPropertiesTable(air_Table);
  tyvek->SetMaterialPropertiesTable(tyvek_Table);
  

  // --------------------------------------
  // ------------- Volumes ----------------
  
  // The Experimental Hall
  //
  G4Box* expHall_box = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);
  
  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,air,"World",0,0,0);
  
  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);
  
  // The Dark Box
  //
  G4Box* darkBox_box = new G4Box("DarkBox",fDarkBox_x,fDarkBox_y,fDarkBox_z);
  
  G4LogicalVolume* darkBox_log
    = new G4LogicalVolume(darkBox_box,air,"DarkBox",0,0,0);
  
  G4VPhysicalVolume* darkBox_phys
    = new G4PVPlacement(0,G4ThreeVector(),darkBox_log,"DarkBox",
                        expHall_log,false,0);
  
  // The Tyvek Panel
  //
  G4Box* tyvekSheet_box = new G4Box("Tyvek",fTyvek_x,fTyvek_y,fTyvek_z);
  
  G4LogicalVolume* tyvekSheet_log
    = new G4LogicalVolume(tyvekSheet_box,
			  tyvek,
			  "Tyvek",
			  0,0,0);
  
  G4RotationMatrix* yRot = new G4RotationMatrix;
  yRot->rotateY(-M_PI/4.*rad);
  
  G4VPhysicalVolume* tyvekSheet_phys =
    new G4PVPlacement(yRot,G4ThreeVector(0,0,0),tyvekSheet_log,"Tyvek",
		      darkBox_log,false,0);
  
  // ------------- Surfaces --------------
  //
  // DarkBox
  //
  G4OpticalSurface* opDarkBoxSurface = new G4OpticalSurface("DarkBoxSurface");
  //  opDarkBoxSurface->SetType(dielectric_dielectric);
  //  opDarkBoxSurface->SetFinish(ground);
  //  opDarkBoxSurface->SetModel(unified);
  
  // G4cout << G4endl;
//   opDarkBoxSurface->SetType(dielectric_LUTDAVIS);
//   opDarkBoxSurface->SetFinish(Rough_LUT);
//   opDarkBoxSurface->SetModel(DAVIS);

  opDarkBoxSurface->SetType(dielectric_metal);
  opDarkBoxSurface->SetFinish(polished);
  opDarkBoxSurface->SetModel(glisur);
  opDarkBoxSurface->SetPolish(0.0); // 0 is max roughness
  
  //
  // Tyvek
  //
  G4OpticalSurface* opTyvekSurface = new G4OpticalSurface("TyvekSurface");
  
  // // below are the defaults anyway
  //   opTyvekSurface->SetType(dielectric_dielectric);
  //   //opTyvekSurface->SetType(dielectric_metal);
  //   opTyvekSurface->SetFinish(polished);
  //   opTyvekSurface->SetModel(glisur);
  //   // 0.0 max roughness - Lambertian
  //   // 1.0     smooth - Snell's 
  //   opTyvekSurface->SetPolish(0.0); 
  
  // UNIFIED MODEL
  //http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/BackupVersions/V10.4/html/_images/UNIFIED_model_diagram.png
     
  opTyvekSurface->SetModel(unified);
  
  
  /*
    dielectric_metal
     polished             - no refraction, if not absorbed then specular spike reflection 
     ground               - no refraction, if reflected then coefficients
  */
  
  // Apply coefficients but no transmission
  // dielectric_metal gives no transmission ie no refraction
  opTyvekSurface->SetType(dielectric_metal);
  opTyvekSurface->SetFinish(ground);
  
  /* 
     dielectric_dielectric
     polished             - if not absorbed then apply Snell's law - ie refract or reflect 
     ground               - as above but if reflected use coefficients
     polishedfrontpainted - no refraction, if not absorbed then specular spike reflection (yes)
     groundfrontpainted   - lambertian reflection or absorption only  
     polishedbackpainted  - out-in: absorption or specular spike, in-out: absorption or coefficients
     groundbackpainted    - out-in: absorption or lambertian, in-out: absorption or coefficients    (ID option)
  */

  // Apply coefficients
  //   opTyvekSurface->SetType(dielectric_dielectric);
  //   opTyvekSurface->SetFinish(polishedbackpainted);
  
  opTyvekSurface->SetSigmaAlpha(0.23); // 0.0 rad means no lobe
  
  G4LogicalBorderSurface* darkBoxSurface =
    new G4LogicalBorderSurface("DarkBoxSurface",
			       darkBox_phys,
			       expHall_phys,
			       opDarkBoxSurface); // opTyvekSurface
 
//   G4LogicalSkinSurface* darkBoxSurface =
//     new G4LogicalSkinSurface("DarkBoxSurface",
// 			     darkBox_phys,
// 			     expHall_phys,
// 			     opDarkBoxSurface); 
  // print info
  G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
    (darkBoxSurface->GetSurface(darkBox_phys,expHall_phys)->
     GetSurfaceProperty());
  
  if (opticalSurface) {
    G4cout << "\n Dark Box surface properties " << G4endl;
    opticalSurface->DumpInfo();
  }


  //   G4LogicalSkinSurface* tyvekSurface =
  //     new G4LogicalSkinSurface("TyvekSurface",
  // 			     tyvekSheet_log,
  // 			     opTyvekSurface);
  

  G4LogicalBorderSurface* tyvekSurface =
    new G4LogicalBorderSurface("TyvekSurface",
			       darkBox_phys,
			       tyvekSheet_phys,
			       opTyvekSurface);
  
  
  // print info
  opticalSurface = dynamic_cast <G4OpticalSurface*>
    (tyvekSurface->GetSurface(darkBox_phys,tyvekSheet_phys)->
     GetSurfaceProperty());
  //(tyvekSurface->GetSurface(tyvekSheet_log)->GetSurfaceProperty());
  
  if (opticalSurface) {
    G4cout << "\n Tyvek Sheet surface properties " << G4endl;
    opticalSurface->DumpInfo();
  }
  
  //------------------------------------------------------
  //---------- Surfaces G4MaterialPropertiesTables -------
  //
  // Generate & Add Material Properties Table attached to the optical surfaces
  //
  
  //--------------------------
  //  Dark Box Surface (border)
  G4double refractiveIndex[nElements] = {1.5, 1.5};
  // G4double specularLobe[nElements]    = {0.0, 0.0};
  // G4double specularSpike[nElements]   = {0.0, 0.0};
  // G4double backScatter[nElements]     = {0.0, 0.0};  

  G4double darkBox_reflectivity[nElements]  = {0.0, 0.0};
  G4double darkBox_efficiency[nElements]    = {0.0, 0.0}; 
  
  G4MaterialPropertiesTable* darkBoxSurf_Table = new G4MaterialPropertiesTable();

  darkBoxSurf_Table->AddProperty("REFLECTIVITY"         ,photon_Energies,
				 darkBox_reflectivity ,nElements);
  darkBoxSurf_Table->AddProperty("EFFICIENCY"           ,photon_Energies,
				 darkBox_efficiency   ,nElements);
  darkBoxSurf_Table->AddProperty("RINDEX",                photon_Energies, refractiveIndex, nElements);

  // darkBoxSurf_Table->AddProperty("SPECULARLOBECONSTANT",  photon_Energies, specularLobe,    nElements);
  // darkBoxSurf_Table->AddProperty("SPECULARSPIKECONSTANT", photon_Energies, specularSpike,   nElements);
  // darkBoxSurf_Table->AddProperty("BACKSCATTERCONSTANT",   photon_Energies, backScatter,     nElements);
  
  G4cout << "\n Dark Box Surface G4MaterialPropertiesTable" << G4endl;
  darkBoxSurf_Table->DumpTable();
  G4cout << G4endl;
  
  opDarkBoxSurface->SetMaterialPropertiesTable(darkBoxSurf_Table);
  
  //--------------------------
  //  Tyvek Surface 
  
  G4double tyvek_R[nElements]   = {0.9, 0.9};
  // no transmission for dielectric_metal
  //  G4double tyvek_Eff[nElements] = {1.0, 1.0}; 
  
  // diffuse Lobe = 1 - (Csl + Css + cBS)
  G4double tyvek_Csl[nElements] = {0.2, 0.2}; // spec lobe
  G4double tyvek_Css[nElements] = {0.0, 0.0}; // spec spike
  G4double tyvek_Cbs[nElements] = {0.0, 0.0}; // back spike
  
  G4MaterialPropertiesTable *tyvekSurf_Table = new G4MaterialPropertiesTable();
  
  tyvekSurf_Table->AddProperty("REFLECTIVITY"         ,photon_Energies,
			       tyvek_R ,nElements);
//   tyvekSurf_Table->AddProperty("EFFICIENCY"           ,photon_Energies,
// 			       tyvek_Eff   ,nElements);
  tyvekSurf_Table->AddProperty("RINDEX"               ,photon_Energies,
			       pe_R_Indices ,nElements);
  tyvekSurf_Table->AddProperty("SPECULARLOBECONSTANT" ,photon_Energies,
			       tyvek_Csl ,nElements);
  tyvekSurf_Table->AddProperty("SPECULARSPIKECONSTANT",photon_Energies,
			       tyvek_Css,nElements);
  tyvekSurf_Table->AddProperty("BACKSCATTERCONSTANT"  ,photon_Energies,
			       tyvek_Cbs  ,nElements);
  
  G4cout << "\n Tyvek Surface G4MaterialPropertiesTable" << G4endl;
  tyvekSurf_Table->DumpTable();
  G4cout << G4endl;
  
  opTyvekSurface->SetMaterialPropertiesTable(tyvekSurf_Table);

  // always return the physical World
  return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

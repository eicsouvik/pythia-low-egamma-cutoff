/**!
 * Example macro for running an event generator in standalone mode.
 *
 * Usage:
 *
 * root4star
 * .L standalone.pythia6.C
 * int nevents=100;
 * standalone( nevents )
 */

#include "TLorentzVector.h"

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary = 0;

class StarPythia6;
StarPythia6 *pythia6 = 0;

class StarPythia8;
StarPythia8 *pythia8 = 0;

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void simrun15mb( Int_t nevents=1e2,
		 UInt_t rngSeed = 12345,
		 TString simulator = "py6:all",
		 //TString simulator = "py8:inel",
		 Double_t cms = 40,
		 TString outfile = "test.root",
		 TString coll = "ep")
{ 

  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = "tables nodefault";
    bfc(0, simple );
  }

  gSystem->Load( "libVMC.so");
  gSystem->Load( "St_g2t.so" );
  gSystem->Load( "St_geant_Maker.so" );
 
  gSystem->Load( "StarGeneratorUtil.so" );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so" );

  gSystem->Load( "libMathMore.so"   );  
  //
  // Create the primary event generator and insert it
  // before the geant maker
  //
  //  StarPrimaryMaker *
  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName( Form("%s.%1.0fGeV.%s.%s",simulator.Data(),cms,coll.Data(),outfile.Data()));
    //  chain -> AddBefore( "geant", _primary );
  }

  //
  // Initialize random number generator
  //
  StarRandom &random = StarRandom::Instance();
  random.capture(); // maps all ROOT TRandoms to StarRandom
  random.seed( rngSeed );

  //
  // Setup cuts on which particles get passed to geant for
  //   simulation.  (To run generator in standalone mode,
  //   set ptmin=1.0E9.)
  //                    ptmin  ptmax
  _primary->SetPtRange  (1.0E9,  -1.0);         // GeV
  //                    etamin etamax
  _primary->SetEtaRange ( -100.0, +100.0 );
  //                    phimin phimax
  _primary->SetPhiRange ( 0., TMath::TwoPi() );
  
  
  // 
  // Setup a realistic z-vertex distribution:
  //   x = 0 gauss width = 1mm
  //   y = 0 gauss width = 1mm
  //   z = 0 gauss width = 30cm
  // 
  _primary->SetVertex( 0., 0., 0. );
  _primary->SetSigma( 0., 0., 0. );

  //
  // Setup an event generator
  //
  if(simulator.Contains("py6",TString::kIgnoreCase))
    {
      Pythia6( simulator, cms);
    }
  else if(simulator.Contains("py8",TString::kIgnoreCase))
    {
      Pythia8( simulator, cms, coll);
    }

  //
  // Initialize primary event generator and all sub makers
  //
  _primary -> Init();

  if(simulator.Contains("py6",TString::kIgnoreCase))
    {
	(StarPythia6::pydat3()).mdcy(102,1)=0; // PI0 111
	(StarPythia6::pydat3()).mdcy(106,1)=0; // PI+ 211
	(StarPythia6::pydat3()).mdcy(109,1)=0; // ETA 221
	(StarPythia6::pydat3()).mdcy(116,1)=0; // K+ 321
	(StarPythia6::pydat3()).mdcy(112,1)=1; // K_SHORT 310
	(StarPythia6::pydat3()).mdcy(105,1)=1; // K_LONG 130
	(StarPythia6::pydat3()).mdcy(164,1)=1; // LAMBDA0 3122
	(StarPythia6::pydat3()).mdcy(167,1)=1; // SIGMA0 3212
	(StarPythia6::pydat3()).mdcy(162,1)=1; // SIGMA- 3112
	(StarPythia6::pydat3()).mdcy(169,1)=1; // SIGMA+ 3222
	(StarPythia6::pydat3()).mdcy(172,1)=1; // Xi- 3312
	(StarPythia6::pydat3()).mdcy(174,1)=1; // Xi0 3322
	(StarPythia6::pydat3()).mdcy(176,1)=1; // OMEGA- 3334
    }
  else if(simulator.Contains("py8",TString::kIgnoreCase))
    {
      pythia8->Set("111:onMode=0"); // pi0 stable to permit mother/daughter in star record
      pythia8->Set("211:onMode=0"); // pi+/- stable
      pythia8->Set("221:onMode=0"); // eta stable
      pythia8->Set("321:onMode=0"); // K+/- stable
      pythia8->Set("310:onMode=1"); // K short
      pythia8->Set("130:onMode=1"); // K long
      pythia8->Set("3122:onMode=1"); // Lambda 0 
      pythia8->Set("3112:onMode=1"); // Sigma -
      pythia8->Set("3222:onMode=1"); // Sigma +
      pythia8->Set("3212:onMode=1"); // Sigma 0
      pythia8->Set("3312:onMode=1"); // Xi -
      pythia8->Set("3322:onMode=1"); // Xi 0
      pythia8->Set("3334:onMode=1"); // Omega -
    }
  //
  // Trigger on nevents
  //
  trig( nevents );

}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void trig( Int_t n=1 )
{
  for ( Int_t i=0; i<n; i++ ) {
    chain->Clear();
    chain->Make();
  }
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void Pythia6( TString mode="py6:minbias", double cms = 200)
{
  gSystem->Load("libLHAPDF.so");
  gSystem->Load("libPythia6_4_28.so");

  pythia6 = new StarPythia6("pythia6");
  pythia6->SetFrame("CMS", cms );
  pythia6->SetBlue("proton");
  pythia6->SetYell("proton");
  PySubs_t &pysubs = pythia6->pysubs();
  PyPars_t &pypars = pythia6->pypars();
  PyDat3_t &pydat3 = pythia6->pydat3();

  if ( mode == "py6.mb" )
    {
      pysubs.msel = 1; 
    }

  if ( mode == "py6.all" )
    {
      pysubs.msel = 2;
      pysubs.msub(91) = 0;
      pysubs.msub(92) = 0;
      pysubs.msub(93) = 0;
    }

  if ( mode == "py6.mbP0" )
    {
      pysubs.msel = 1; 
      pythia6->PyTune( 320 );
    }

  if ( mode == "py6.allP0" )
    {
      pysubs.msel = 2; 
      pythia6->PyTune( 320 );
    }

  if ( mode.Contains("spin" ) )
    {
      if( mode == "py6.spin" ) 
	{
	  pysubs.msel = 1;
	}
      if ( mode == "py6.spinall" )
	{
	  pysubs.msel = 2;
	}
     if(mode == "py6:spinsd")
	{
	  pysubs.msel = 0;
	  pysubs.msub(92) = 1;
	  pysubs.msub(93) = 1;
	}
      if(mode == "py6:spindd")
	{
	  pysubs.msel = 0;
	  pysubs.msub(94) = 1;
	}

      pypars.mstp(5) = 370;     // Perugia 2012 tune Pythia 6.4.28  
      pypars.mstp(51) = 10042;  // CTEQ6L1 PDF 10042 MSTW 21000

      pythia6->Init();
      pypars.mstp(5) = 0;

      //See page 59 of https://drupal.star.bnl.gov/STAR/system/files/AnalysisNoteV1.pdf
      pypars.parp(90) = 0.213; //exponent for pT_0 which controls underlying events
    }

  if ( mode.Contains("P12" ) )
    {
      if ( mode == "py6.mbP12" )
	{
	  pysubs.msel = 1; 
	}
      if ( mode == "py6.allP12" )
	{
	  pysubs.msel = 2;
	}
      pypars.mstp(5) = 370;     // Perugia 2012 tune Pythia 6.4.28  
      pypars.mstp(51) = 10042;  // CTEQ6L1 PDF 10042 MSTW 21000
      pythia6->Init();
    }

  _primary->AddGenerator(pythia6);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
//
//
void Pythia8( TString mode="py8:mb", double cms = 200, TString coll = "ep")
{
  gSystem->Load( "libLHAPDF.so"  ); // LHAPDF needs to be called before PYTHIA8

  if( coll == "pn" ) gSystem->Load( "libPythia8_1_86.so");
  else    gSystem->Load( "libPythia8_3_03.so");

    Double_t Pe[4];
    Double_t Pp[4];

    TLorentzVector P_temp1;
    TLorentzVector P_temp2;

    const double Ee = 10.0, Ep = 100.0;
    const double me = 0.000511, mp = 0.938272;

    P_temp1.SetPxPyPzE(0, 0, -sqrt(Ee*Ee - me*me), Ee);
    P_temp2.SetPxPyPzE(0, 0, +sqrt(Ep*Ep - mp*mp), Ep);

    Pe[0] = P_temp1.Px(); // Px
    Pe[1] = P_temp1.Py(); // Py
    Pe[2] = P_temp1.Pz(); // Pz
    Pe[3] = P_temp1.E();  // Energy
    std::cout << "Pe: Px = " << Pe[0] << ", Py = " << Pe[1] << ", Pz = " << Pe[2] << ", E = " << Pe[3] << std::endl;

    Pp[0] = P_temp2.Px(); // Px
    Pp[1] = P_temp2.Py(); // Py
    Pp[2] = P_temp2.Pz(); // Pz
    Pp[3] = P_temp2.E();  // Energy
    std::cout << "Pp: Px = " << Pp[0] << ", Py = " << Pp[1] << ", Pz = " << Pp[2] << ", E = " << Pp[3] << std::endl;


// Calculate sqrt(s) using s = (P1 + P2)^2
     Double_t s1 = (P_temp1 + P_temp2).M2(); // Square of invariant mass, which is s

         std::cout << "Method 1 sqrt(s): " << std::sqrt(s1) << std::endl;


// Calculate sqrt(s)
     Double_t s2 = std::sqrt(std::pow(Pe[3] + Pp[3], 2) - (std::pow(Pe[0] + Pp[0], 2) + std::pow(Pe[1] + Pp[1], 2) + std::pow(Pe[2] + Pp[2], 2)));

         std::cout << "Method 2 sqrt(s): " << s2 << std::endl;


 pythia8 = new StarPythia8();    

  pythia8->SetFrame("4MOM", Pe,Pp);

  pythia8->SetBlue("electron");
  pythia8->SetYell("proton");
  
  const double beamRapidity = P_temp2.Rapidity();
   
  std::cout<<"Proton rapidity = " << beamRapidity << std::endl;

  pythia8->Set("PDF:lepton2gamma = on");
  pythia8->Set("WeakBosonExchange:ff2ff(t:gmZ) = off");
  pythia8->Set("Photon:Q2max = 0.01");
  pythia8->Set("Photon:Wmin = 5.0");
  pythia8->Set("Photon:Wmax = -1.0");

  pythia8->Set("Photon:ProcessType = 1");

  pythia8->Set("SoftQCD:all = off");
  pythia8->Set("HardQCD:all = off");
  pythia8->Set("PhotonParton:all = off");

  pythia8->Set("SoftQCD:nonDiffractive = on");

  pythia8->Set("PartonLevel:MPI = on");
  pythia8->Set("PartonLevel:ISR = on");
  pythia8->Set("PartonLevel:FSR = on");
  pythia8->Set("HadronLevel:all = on");

  pythia8->Set("TimeShower:QEDshowerByL = off");

  pythia8->Set("ColourReconnection:reconnect = on");

    if (mode.Contains("CR2", TString::kIgnoreCase)) {
      pythia8->Set("StringPT:sigma = 0.335");
      pythia8->Set("StringZ:aLund = 0.36");
      pythia8->Set("StringFlav:probQQtoQ = 0.078");
      pythia8->Set("StringFlav:probStoUD = 0.2");
      pythia8->Set(
          "StringFlav:probQQ1toQQ0join = "
          "0.0275, 0.0275, 0.0275, 0.0275");

      pythia8->Set("MultipartonInteractions:pT0Ref = 2.15");
      pythia8->Set("BeamRemnants:remnantMode = 1");
      pythia8->Set("BeamRemnants:saturation = 5");

      pythia8->Set("ColourReconnection:mode = 1");
      pythia8->Set("ColourReconnection:allowDoubleJunRem = off");
      pythia8->Set("ColourReconnection:m0 = 0.3");
      pythia8->Set("ColourReconnection:allowJunctions = on");
      pythia8->Set("ColourReconnection:junctionCorrection = 1.20");
      pythia8->Set("ColourReconnection:timeDilationMode = 2");
      pythia8->Set("ColourReconnection:timeDilationPar = 0.18");
    }
    else {
      pythia8->Set("ColourReconnection:mode = 0");
    }

 // Helpful during setup:
 pythia8->Set("Init:showProcesses = on");
 pythia8->Set("Init:showChangedSettings = on");

//

  _primary->AddGenerator(pythia8);
}
// ----------------------------------------------------------------------------

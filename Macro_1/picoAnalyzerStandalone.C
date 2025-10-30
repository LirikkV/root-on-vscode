/**
 * \brief Example of how to read a file (list of files) using StPicoEvent class in a standalone mode
 *
 * picoAnalyzer.C is an example of reading STAR picoDst format in a standalone mode
 * on your laptop or local computer farm.
 * Prerequisites:
 * - StPicoEvent directory with the library (libStPicoDst.so) compiled
 * - CERN ROOT package
 * - g++ >= 4.8
 * - Makefile
 *
 * First, the program must be compiled with the Makefile, with simple command in the bash shell:
 * make
 *
 * Then the executable file picoAnalyzerStandalone will be created. The current version of the program
 * expects 3 arguments: ./picoAnalyzerStandalone inputFile outputFile
 * The first one is the program name, the second one is the name of the inputfile that 
 * maybe either the picoDst file itself, in a format dummyname.picoDst.root or a list of
 * such files called dummyname.list or dummyname.lis. The outputFile assumes the some_output_name.root.
 *
 * \author Grigory Nigmatkulov
 * \date August 6, 2019
 * \email nigmatkulov@gmail.com
 */

// C++ headers
#include <iostream>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

// PicoDst headers
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"
#include "StPicoTrack.h"
#include "StPicoBTofHit.h"
#include "StPicoBTowHit.h"
#include "StPicoEmcTrigger.h"
#include "StPicoBTofPidTraits.h"
#include "StPicoTrackCovMatrix.h"
#include "StPicoFmsHit.h"
#include "StPicoETofHit.h"
#include "StPicoEpdHit.h"

//_________________
int main(int argc, char* argv[]) {

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  R__LOAD_LIBRARY(libStPicoDst);
#else
  gSystem->Load("../libs/libStPicoDst.so");
#endif

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  const char* fileName;
  const char* oFileName;

  switch (argc) {
  case 3:
    fileName = argv[1];
    oFileName = argv[2];
    break;
  default:
    std::cout << "Usage: picoAnalyzerStandalone inputFileName outputFileName.root" << std::endl;
    return -1;
  }
  std::cout << " inputFileName : " << fileName << std::endl;
  std::cout << " outputFileName: " << oFileName << std::endl;
  
  StPicoDstReader* picoReader = new StPicoDstReader(fileName);
  picoReader->Init();

  // This is a way if you want to spead up I/O
  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event*", 1);
  picoReader->SetStatus("Track*", 1);
  picoReader->SetStatus("BTofHit*", 1);
  picoReader->SetStatus("BTofPidTraits*", 1);
  picoReader->SetStatus("BTowHit*", 1);
  picoReader->SetStatus("ETofHit*", 1);
  picoReader->SetStatus("EpdHit*", 1);
  //picoReader->SetStatus("EmcTrigger",0);
  //picoReader->SetStatus("TrackCovMatrix",1);
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  if( !picoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;
  Long64_t events2read = picoReader->chain()->GetEntries();

  std::cout << "Number of events to read: " << events2read
	    << std::endl;


  TFile *oFile = new TFile(oFileName, "recreate");
  
  // Histogramming
  // Event
  TH1F *hRefMult = new TH1F("hRefMult",
			    "Reference multiplicity;refMult",
			    500, -0.5, 499.5);
  TH2F *hVtxXvsY = new TH2F("hVtxXvsY",
			    "hVtxXvsY",
			    200,-10.,10.,200,-10.,10.);
  TH1F *hVtxZ = new TH1F("hVtxZ","hVtxZ",
			 140, -70., 70.);
  //after event cut
  TH2F *hVtxXvsY_cut = new TH2F("hVtxXvsY_cut",
			    "hVtxXvsY_cut",
			    200,-10.,10.,200,-10.,10.);
  TH1F *hVtxZ_cut = new TH1F("hVtxZ_cut","hVtxZ_cut",
			 140, -70., 70.);

  // Track
  TH1F *hGlobalPtot = new TH1F("hGlobalPtot",
			       "Global track momentum;p (GeV/c)",
			       100, 0., 1. );
  TH1F *hGlobalPtotCut = new TH1F("hGlobalPtotCut",
				  "Global track momentum after cut;p (GeV/c)",
				  100, 0., 1. );
  TH1F *hPrimaryPtot = new TH1F("hPrimaryPtot",
				"Primary track momentum;p (GeV/c)",
			       100, 0., 2. );
  //after track cut:
  TH1F *hPrimaryPtot_cut = new TH1F("hPrimaryPtot_cut",
				   "Primary track momentum after cut;p (GeV/c)",
				  100, 0., 2. );
  TH1F *hTransvMomentum = new TH1F("hTransvMomentum",
				   "Track transverse momentum;p_{T} (GeV/c)",
				   200, 0., 2.);
  TH2F *hGlobalPhiVsPt[2];
  for(int i=0; i<2; i++) {
    hGlobalPhiVsPt[i] = new TH2F(Form("hGlobalPhiVsPt_%d",i),
				 Form("#phi vs. p_{T} for charge: %d;p_{T} (GeV/c);#phi (rad)", (i==0) ? 1 : -1),
				 300, 0., 3.,
				 630, -3.15, 3.15);
  }
  TH1F *hNSigmaPion = new TH1F("hNSigmaPion",
			       "n#sigma(#pi);n#sigma(#pi)",
			       400, -10., 10.);
  TH1F *hNSigmaElectron = new TH1F("hNSigmaElectron",
				   "n#sigma(e);n#sigma(e)",
				   400,-10.,10.);
  TH1F *hNSigmaKaon = new TH1F("hNSigmaKaon",
			       "n#sigma(K);n#sigma(K)",
			       400, -10., 10.);
  TH1F *hNSigmaProton = new TH1F("hNSigmaProton",
				 "n#sigma(p);n#sigma(p)",
				 400, -10., 10.);       
  TH1I *hNFitHits = new TH1I("hNFitHits", "Number of hits in TPC for track fit",
          61, -0.5,60.5);
  TH1I *hNFitHits_cut = new TH1I("hNFitHits_cut", "Number of hits in TPC for track fit after cut",
          61, -0.5,60.5);
  TH1F *hDCA = new TH1F("hDCA", "DCA to ptimary vertex", 
          100, 0.0,5.0);
  TH1F *hDCA_cut = new TH1F("hDCA_cut", "DCA to ptimary vertex after cut", 
          100, 0.0,5.0);
  TH1F *hPrimaryPtrans = new TH1F("hPrimaryPtrans",
				"Primary track trancverse momentum;p (GeV/c)",
			       100, 0., 2. );
  TH1F *hPrimaryPtrans_cut = new TH1F("hPrimaryPttrans_cut",
				   "Primary track transverce momentum after cut;p (GeV/c)",
				  100, 0., 2. );
  TH1F *hPrimaryPseudorap = new TH1F("hPrimaryPseudorap",
				   "Primary track pseudorapidity",
				  100, -2., 2. );
  TH1F *hPrimaryPseudorap_cut = new TH1F("hPrimaryPseudorap_cut",
				   "Primary track pseudorapidity after cut",
				  100, -2., 2. );
  
  TH2F *hNSigmPion_vs_pPrimTotDevQ = new TH2F("hNSigmPion_vs_pPrimTotDevQ",
			    "nSigma(pion) vs P_prim_tot/q;;nSigma",
			    400,-2.,2.,400,-60.,60.);
  TH2F *hNSigmPion_vs_pPrimTotDevQ_cut_PID = new TH2F("hNSigmPion_vs_pPrimTotDevQ_cut_PID",
			    "nSigma(pion) vs P_prim_tot/q after cut;;nSigma",
			    400,-2.,2.,400,-60.,60.);

  TH2F *h1_OverBeta_vs_pPrimTotDevQ = new TH2F("1_OverBeta_vs_pPrimTotDevQ",
			    "1/beta vs P_prim_tot/q;;1/beta",
			    200,-2.,2.,200,-1.,10.);
  TH2F *h1_OverBeta_vs_pPrimTotDevQ_cut_PID = new TH2F("1_OverBeta_vs_pPrimTotDevQ_cut_PID",
			    "1/beta vs P_prim_tot/q after PID;;1/beta",
			    200,-2.,2.,200,-1.,10.);
  
  TH2F *hm2_vs_pPrimTotDevQ = new TH2F("hm2_vs_pPrimTotDevQ",
			    "m^2 vs P_prim_tot/q;;m^2(Gev/c)",
			    200,-2.,2.,200,-1.,10.);
  TH2F *hm2_vs_pPrimTotDevQ_cut_PID = new TH2F("hm2_vs_pPrimTotDevQ_cut_PID",
			    "m^2 vs P_prim_tot/q after PID;;m^2(Gev/c)",
			    200,-2.,2.,200,-1.,10.);
  
  TH2F *h1_OverBetaDelta_vs_pPrimTotDevQ = new TH2F("1_OverBetaDelta_vs_pPrimTotDevQ",
			    "|1/beta - 1/beta_exp| vs P_prim_tot/q;;|1/beta - 1/beta_exp|",
			    200,-2.,2.,200,-1.,10.);
  TH2F *h1_OverBetaDelta_vs_pPrimTotDevQ_cut_PID = new TH2F("1_OverBetaDelta_vs_pPrimTotDevQ_cut_PID",
			    "|1/beta - 1/beta_exp| vs P_prim_tot/q after PID;;|1/beta - 1/beta_exp|",
			    200,-2.,2.,200,-2.,2.);


  // BTof pid traits
  TH1F *hTofBeta = new TH1F("hTofBeta", "BTofPidTraits #beta;#beta",
			    2000, 0., 2.);

  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    std::cout << "Working on event #[" << (iEvent+1)
	      << "/" << events2read << "]" << std::endl;

    Bool_t readEvent = picoReader->readPicoEvent(iEvent);
    if( !readEvent ) {
      std::cout << "Something went wrong, Master! Nothing to analyze..."
		<< std::endl;
      break;
    }

    // Retrieve picoDst
    StPicoDst *dst = picoReader->picoDst();

    // Retrieve event information
    StPicoEvent *event = dst->event();
    if( !event ) {
      std::cout << "Something went wrong, Master! Event is hiding from me..."
		<< std::endl;
      break;
    }
    hRefMult->Fill( event->refMult() );

    TVector3 pVtx = event->primaryVertex();
    hVtxXvsY->Fill( event->primaryVertex().X(), event->primaryVertex().Y() );
    hVtxZ->Fill( event->primaryVertex().Z() );

    // Lirikk's cuts:
    // variables for event cut:
    double_t Vtx_r_Max = 2.0;  // cm
    double_t Vtx_z_Max = 40.0; // cm

    // Event selection:
    Bool_t is_Vtx_r_cut = event->primaryVertex().Perp() < Vtx_r_Max;
    Bool_t is_abs_Vtx_z_cut = fabs(event->primaryVertex().Z()) < Vtx_z_Max;
    if (is_Vtx_r_cut && is_abs_Vtx_z_cut)
    {
      //filling QA hists after cut:
      hVtxXvsY_cut->Fill(event->primaryVertex().X(), event->primaryVertex().Y());
      hVtxZ_cut->Fill(event->primaryVertex().Z());
    

    // Track analysis
    Int_t nTracks = dst->numberOfTracks();
    Int_t nMatrices = dst->numberOfTrackCovMatrices();
    if(nTracks != nMatrices) {
      //std::cout << "Number of tracks and matrices do not match!" << std::endl;
    }
    //std::cout << "Number of tracks in event: " << nTracks << std::endl;
    
    // Track loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

      // Retrieve i-th pico track
      StPicoTrack *picoTrack = dst->track(iTrk);
      
      if(!picoTrack) continue;
      //std::cout << "Track #[" << (iTrk+1) << "/" << nTracks << "]"  << std::endl;

      hGlobalPtot->Fill( picoTrack->gMom().Mag() );

      // Simple single-track cut
      if( picoTrack->gMom().Mag() < 0.1 ||
	  picoTrack->gDCA(pVtx).Mag()>50. ) {
	continue;
      } 

      hGlobalPtotCut->Fill( picoTrack->gMom().Mag() );
      


      if( picoTrack->charge() > 0 ) {
	hGlobalPhiVsPt[0]->Fill( picoTrack->gMom().Pt(),
				 picoTrack->gMom().Phi() );
      }
      else {
	hGlobalPhiVsPt[1]->Fill( picoTrack->gMom().Pt(),
				 picoTrack->gMom().Phi() );	
      }
      hNSigmaElectron->Fill( picoTrack->nSigmaElectron() );
      hNSigmaPion->Fill( picoTrack->nSigmaPion() );
      hNSigmaKaon->Fill( picoTrack->nSigmaKaon() );
      hNSigmaProton->Fill( picoTrack->nSigmaProton() );
      
      hTransvMomentum->Fill( picoTrack->gMom().Pt() );

      // Check if track has TOF signal
      if( picoTrack->isTofTrack() ) {
	// Retrieve corresponding trait
	StPicoBTofPidTraits *trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
	if( !trait ) {
	  std::cout << "O-oh... No BTofPidTrait # " << picoTrack->bTofPidTraitsIndex()
		    << " for track # " << iTrk << std::endl;
	  std::cout << "Check that you turned on the branch!" << std::endl;
	  continue;
	}
	// Fill beta
	hTofBeta->Fill( trait->btofBeta() );
      } //if( isTofTrack() )
      
    //Lirikk's QA hists before track selection:
    hNFitHits->Fill(picoTrack->nHitsFit());
    if (picoTrack->isPrimary())
    {
      hPrimaryPtot->Fill(picoTrack->pMom().Mag());
      hPrimaryPtrans->Fill(picoTrack->pMom().Pt());
      hPrimaryPseudorap->Fill(picoTrack->pMom().Eta());
    }
    hDCA->Fill(picoTrack->gDCA(pVtx).Mag());
    //Track selection:
    // variables for track cut:
    Int_t N_TPC_fit_hits_min = 15;
    Double_t DCA_max = 3.0;//cm
    Double_t p_tot_prim_min = 0.15;//Gev/c
    Double_t p_tot_prim_max = 1.5;//Gev/c
    Double_t p_trans_prim_min = 0.15;
    Double_t p_trans_prim_max = 1.5;
    Double_t pseudo_rap_prim_max = 1.0;
    //track selection:
    Bool_t is_N_TPC_fit_hits_cut = picoTrack->nHitsFit()>=N_TPC_fit_hits_min;
    Bool_t is_p_tot_prim_cut = picoTrack->isPrimary() && 
                            p_tot_prim_min < picoTrack->pMom().Mag() &&
                            picoTrack->pMom().Mag()<p_tot_prim_max;
    Bool_t is_DCA_abs_cut = picoTrack->gDCA(pVtx).Mag()<DCA_max;
    Bool_t is_p_trans_prim_cut = picoTrack->isPrimary() && 
                            p_trans_prim_min<picoTrack->pMom().Pt() &&
                            picoTrack->pMom().Pt()<p_trans_prim_max;
    Bool_t is_pseudo_prm_cut = picoTrack->isPrimary() &&
                            fabs(picoTrack->pMom().Eta())<pseudo_rap_prim_max;
    if(is_N_TPC_fit_hits_cut && is_p_tot_prim_cut && 
      is_DCA_abs_cut && is_p_trans_prim_cut && 
      is_pseudo_prm_cut) 
    {
      hNFitHits_cut->Fill(picoTrack->nHitsFit());
	    hPrimaryPtot_cut->Fill( picoTrack->pMom().Mag() );
      hDCA_cut->Fill(picoTrack->gDCA(pVtx).Mag());
      hPrimaryPtrans_cut->Fill(picoTrack->pMom().Pt());
      hPrimaryPseudorap_cut->Fill(picoTrack->pMom().Eta());

      
      
      //start of PID:
      //variables for PID:
      Double_t PtotPrimQ = (picoTrack->pMom().Mag())/(picoTrack->charge());

      //constants for PID:
      Double_t p_tot_prim_mid_PID = 0.55;//Gev/c
      Double_t one_over_beta_delta_max = 0.015;
      Float_t nSigmaPion_max = 3.0;

      //Lirikk's QA before PID TOF check:
      hNSigmPion_vs_pPrimTotDevQ->Fill(PtotPrimQ, picoTrack->nSigmaPion());
      
      //TPC only:

      //TPC+TOF:
      if( picoTrack->isTofTrack() ) {
      StPicoBTofPidTraits *trait = dst->btofPidTraits(picoTrack->bTofPidTraitsIndex()); // Retrieve corresponding trait

      //variables for QA:
      Double_t m_square = picoTrack->pMom().Mag2()*(1./((trait->btofBeta())*(trait->btofBeta()))-1.);
      Double_t m_Pion = 0.13957039;
      Double_t one_beta_expect = sqrt(m_Pion*m_Pion+(picoTrack->pMom().Mag2()))/(picoTrack->pMom().Mag());

      //QA histograms after TOF check before PID
      h1_OverBeta_vs_pPrimTotDevQ->Fill(PtotPrimQ,1./(trait->btofBeta()));
      hm2_vs_pPrimTotDevQ->Fill(PtotPrimQ,m_square);
      h1_OverBetaDelta_vs_pPrimTotDevQ->Fill(PtotPrimQ, fabs(1./(trait->btofBeta()) - one_beta_expect));

      Bool_t is_p_tot_tof = p_tot_prim_mid_PID<picoTrack->pMom().Mag(); //p_max already set by track choise
      Bool_t is_nSigma_Pion = fabs(picoTrack->nSigmaPion())<nSigmaPion_max;
      Bool_t is_1_beta_delta = fabs(1./(trait->btofBeta())-one_beta_expect)<one_over_beta_delta_max;
      if(is_p_tot_tof && is_nSigma_Pion 
          && is_1_beta_delta)
      {
        //Fill hists after PID cut:
        hNSigmPion_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaPion());
        h1_OverBeta_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ,1./(trait->btofBeta()));
        hm2_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ,m_square);
        h1_OverBetaDelta_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, fabs(1./(trait->btofBeta()) - one_beta_expect));
      }//end of PID

    }//end of TOF

    }//end of track selection
    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    
    }//end of event selection
  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  picoReader->Finish();
  oFile->Write();
  oFile->Close();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	    << std::endl;
}

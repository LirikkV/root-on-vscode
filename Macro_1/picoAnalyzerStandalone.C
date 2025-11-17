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
#include <vector>
#include <deque>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include <TLorentzVector.h>

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
void comparePionsFillHistB(const std::vector<TLorentzVector>& new_Pions_Arr,
                           const std::deque<std::vector<TLorentzVector>>& event_Queue,
                           TH1D* hist_B)
{
  //loop over queue vectors of pions:
  for(size_t i=0;i<event_Queue.size();i++)
  {
    //loop over pions from new event
    for(size_t j=0;j<new_Pions_Arr.size();j++)
    {
      //loop over pions from selected queue event:
      for(size_t k=0;k<event_Queue[i].size();k++)
      {
      TLorentzVector delta_4_momenta = new_Pions_Arr[j] - event_Queue[i][k];
      double_t q_inv = sqrt(-delta_4_momenta.Mag2());

      hist_B->Fill(q_inv);
      }
    }
  }
}

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
			    500, -0.5, 600.);
  TH2F *hVtxXvsY = new TH2F("hVtxXvsY",
			    "Vtx XvsY",
			    200,-10.,10.,200,-10.,10.);
  TH2F *hVtxXvsY_cut = new TH2F("hVtxXvsY_cut",
			    "Vtx XvsY after event cut",
			    200,-10.,10.,200,-10.,10.);

  TH1F *hVtxZ = new TH1F("hVtxZ","Vtx Z",
			 140, -70., 70.);
  TH1F *hVtxZ_cut = new TH1F("hVtxZ_cut","Vtx Z after event cut",
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
  TH1F *hPrimaryPtot_cut = new TH1F("hPrimaryPtot_cut",
				   "Primary track momentum after track cut;p (GeV/c)",
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
			    200,-2.,2.,200,-60.,60.);
  TH2F *hNSigmPion_vs_pPrimTotDevQ_cut_PID = new TH2F("hNSigmPion_vs_pPrimTotDevQ_cut_PID",
			    "nSigma(pion) vs P_prim_tot/q after cut;;nSigma",
			    200,-2.,2.,200,-5.,5.);

  TH2F *hNSigmKaon_vs_pPrimTotDevQ = new TH2F("hNSigmKaon_vs_pPrimTotDevQ",
			    "nSigma(Kaon) vs P_prim_tot/q;;nSigma",
			    200,-2.,2.,200,-60.,60.);
  TH2F *hNSigmKaon_vs_pPrimTotDevQ_cut_PID = new TH2F("hNSigmKaon_vs_pPrimTotDevQ_cut_PID",
			    "nSigma(Kaon) vs P_prim_tot/q after cut;;nSigma",
			    200,-2.,2.,200,-30.,5.);
  
  TH2F *hNSigmProton_vs_pPrimTotDevQ = new TH2F("hNSigmProton_vs_pPrimTotDevQ",
			    "nSigma(Proton) vs P_prim_tot/q;;nSigma",
			    200,-2.,2.,200,-60.,60.);
  TH2F *hNSigmProton_vs_pPrimTotDevQ_cut_PID = new TH2F("hNSigmProton_vs_pPrimTotDevQ_cut_PID",
			    "nSigma(Proton) vs P_prim_tot/q after cut;;nSigma",
			    200,-2.,2.,200,-35.,5.);
  
  TH2F *hNSigmElectron_vs_pPrimTotDevQ = new TH2F("hNSigmElectron_vs_pPrimTotDevQ",
			    "nSigma(Electron) vs P_prim_tot/q;;nSigma",
			    200,-2.,2.,200,-60.,60.);
  TH2F *hNSigmElectron_vs_pPrimTotDevQ_cut_PID = new TH2F("hNSigmElectron_vs_pPrimTotDevQ_cut_PID",
			    "nSigma(Electron) vs P_prim_tot/q after cut;;nSigma",
			    200,-2.,2.,200,-10.,5.);

  TH2F *h1_OverBeta_vs_pPrimTotDevQ = new TH2F("1_OverBeta_vs_pPrimTotDevQ",
			    "1/beta vs P_prim_tot/q;;1/beta",
			    200,-2.,2.,200,-1.,10.);
  TH2F *h1_OverBeta_vs_pPrimTotDevQ_cut_PID = new TH2F("1_OverBeta_vs_pPrimTotDevQ_cut_PID",
			    "1/beta vs P_prim_tot/q after PID;;1/beta",
			    400,-2.,2.,200,0.95,1.1);
  
  TH2F *hm2_vs_pPrimTotDevQ = new TH2F("hm2_vs_pPrimTotDevQ",
			    "m^2 vs P_prim_tot/q;;m^2(Gev/c)",
			    200,-2.,2.,200,-1.,10.);
  TH2F *hm2_vs_pPrimTotDevQ_cut_PID = new TH2F("hm2_vs_pPrimTotDevQ_cut_PID",
			    "m^2 vs P_prim_tot/q after PID;;m^2(Gev/c)",
			    400,-2.,2.,200,-0.1,0.1);
  
  TH2F *h1_OverBetaDelta_vs_pPrimTotDevQ = new TH2F("1_OverBetaDelta_vs_pPrimTotDevQ",
			    "1/beta - 1/beta_exp vs P_prim_tot/q;;1/beta - 1/beta_exp",
			    200,-2.,2.,200,-1.,10.);
  TH2F *h1_OverBetaDelta_vs_pPrimTotDevQ_cut_PID = new TH2F("1_OverBetaDelta_vs_pPrimTotDevQ_cut_PID",
			    "1/beta - 1/beta_exp vs P_prim_tot/q after PID;;1/beta - 1/beta_exp",
			    400,-2.,2.,200,-0.02,0.04);

  TH2F *hdEdx_vs_pPrimTotDevQ = new TH2F("hdEdx_vs_pPrimTotDevQ",
			    "dE/dx vs P_prim_tot/q;;dE/dx(GeV/cm)",
			    200,-1.6,1.6,200,-0.2,10.);
  TH2F *hdEdx_vs_pPrimTotDevQ_cut_PID = new TH2F("hdEdx_vs_pPrimTotDevQ_cut_PID",
			    "dE/dx vs P_prim_tot/q after PID;;dE/dx(GeV/cm)",
			    200,-1.6,1.6,200,-0.2,10. );
  
  TH2F *h2DpPrimTr_vs_etaPtim = new TH2F("h2DpPrimTr_vs_etaPtim",
			    "p_prim_T vs eta_prim;p_T;eta",
			    200,-0.5,2.,200,-4.,4.);
  TH2F *h2DpPrimTr_vs_etaPtim_cut = new TH2F("h2DpPrimTr_vs_etaPtim_cut",
			    "p_prim_T vs eta_prim after track cut;p_T;eta",
			    200,-0.5,2.,200,-4.,4.);
        
  TH2F *h2DpPrimTr_vs_etaPtim_equal_P = new TH2F("h2DpPrimTr_vs_etaPtim_equal_P",
			    "p_prim_T vs eta_prim after track cut with equal P;p_T;eta",
			    200,-0.5,2.,200,-4.,4.);

  // BTof pid traits
  TH1F *hTofBeta = new TH1F("hTofBeta", "BTofPidTraits #beta;#beta",
			    2000, 0., 2.);

  //Correlation function:

  TH1D *hA_q_inv_ALL = new TH1D("hA_q_inv_ALL",
				   "Numerator of Corr.Funct with both TPC & TPC+TOF methods",
				  300, 0., 3.0 );
  TH1D *hB_q_inv_ALL = new TH1D("hB_q_inv_ALL",
				   "Denumerator of Corr.Funct with both TPC & TPC+TOF methods",
				  300, 0., 3.0 );
  
  //cuts by Vz: 4 cats; Vz from -40 to 40
  //cuts by refMult: 10 cuts; RefMult from 0 to 600
  const Int_t nVzCuts = 4;
  const Int_t nRefMultCuts = 10;
  const Double_t VzBins[nVzCuts+1] = {-40., -20., 0., 20., 40.};
  const Double_t RefMultBins[nRefMultCuts+1] = {0.,60.,120.,180.,240.,300.,360.,420.,480.,540.,600};
  /*
  TH1D* hB_hists[nVzCuts][nRefMultCuts]; //array for 40 hists, first index - Vz, second index - RefMult

  for(int iVz=0;iVz<nVzCuts;iVz++)
  {
    for(int iRefM=0;iRefM<nRefMultCuts;iRefM++)
    {
      TString B_histName = Form("hist_%d_%d", iVz, iRefM);
      TString B_histTitle = Form("Denum. of CF Vz [%.0f,%.0f], RefMult [%.0f,%.0f]",VzBins[iVz],VzBins[iVz+1], 
                                                                        RefMultBins[iRefM],RefMultBins[iRefM + 1]);
      hB_hists[iVz][iRefM] = new TH1D(B_histName,B_histTitle,100,-0.1,3.0);
    }
  }
  */

  //for mixing events:
  const Int_t BUFFER_SIZE = 5;
  std::deque<std::vector<TLorentzVector>> Pions_mix_queue_Arr_4_mom; //this is queue from events; just queue from vectors from 4-vectors of particle
  std::vector<std::vector<std::deque<std::vector<TLorentzVector>>>> Pions_Buffer(
                                          nVzCuts, 
                                          std::vector<std::deque<std::vector<TLorentzVector>>>(nRefMultCuts));//this is 4*10 queues from events; just 2D vector of previous queues for each cut

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

    //let's create a c++ vector with 4-momenta of pions in one event:
    std::vector<TLorentzVector> Pions_4_momenta_Arr_TPC_ONLY;
    std::vector<TLorentzVector> Pions_4_momenta_Arr_TOF_TPC;
    std::vector<TLorentzVector> Pions_4_momenta_Arr_ALL;


    // Retrieve picoDst
    StPicoDst *dst = picoReader->picoDst();

    // Retrieve event information
    StPicoEvent *event = dst->event();
    if( !event ) {
      std::cout << "Something went wrong, Master! Event is hiding from me..."
		<< std::endl;
      break;
    }
    

    TVector3 pVtx = event->primaryVertex();
    //QA hists before event selection:
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
      //QA hists after event cut:
      hVtxXvsY_cut->Fill(event->primaryVertex().X(), event->primaryVertex().Y());
      hVtxZ_cut->Fill(event->primaryVertex().Z());
      hRefMult->Fill( event->refMult() );

    // Track analysis
    Int_t nTracks = dst->numberOfTracks();
    // Track loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

      // Retrieve i-th pico track
      StPicoTrack *picoTrack = dst->track(iTrk);
      
      if(!picoTrack) continue;
      
    //QA hists before track selection:
    hNFitHits->Fill(picoTrack->nHitsFit());
    hDCA->Fill(picoTrack->gDCA(pVtx).Mag());  
    if (picoTrack->isPrimary())
    {
      hPrimaryPtot->Fill(picoTrack->pMom().Mag());
      hPrimaryPtrans->Fill(picoTrack->pMom().Pt());
      hPrimaryPseudorap->Fill(picoTrack->pMom().Eta());
      h2DpPrimTr_vs_etaPtim->Fill(picoTrack->pMom().Pt(),picoTrack->pMom().Eta());
    }
    
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
      //QA after track selection:
      hNFitHits_cut->Fill(picoTrack->nHitsFit());
	    hPrimaryPtot_cut->Fill( picoTrack->pMom().Mag() );
      hDCA_cut->Fill(picoTrack->gDCA(pVtx).Mag());
      hPrimaryPtrans_cut->Fill(picoTrack->pMom().Pt());
      hPrimaryPseudorap_cut->Fill(picoTrack->pMom().Eta());
      h2DpPrimTr_vs_etaPtim_cut->Fill(picoTrack->pMom().Pt(),picoTrack->pMom().Eta());

      //Lines of equal P at 2D hist P_prim_T VS Pseudorap
      for(int i=1;i<=21;i+=2)
      {
        Bool_t is_P_const = (p_tot_prim_max-p_tot_prim_min)/(20.)*(i)<=picoTrack->pMom().Mag() &&
                            picoTrack->pMom().Mag()<= (p_tot_prim_max-p_tot_prim_min)/(20.)*(i+1);
        if(is_P_const)
        {
        h2DpPrimTr_vs_etaPtim_equal_P->Fill(picoTrack->pMom().Pt(),picoTrack->pMom().Eta());
        }
      }

      //start of PID:
      //variables for PID:
      Double_t PtotPrimQ = (picoTrack->pMom().Mag())/(picoTrack->charge());

      //constants for PID:
      Double_t p_tot_prim_mid_PID = 0.55;//Gev/c
      Double_t one_over_beta_delta_max = 0.015;
      Float_t nSigmaPion_max_TOF = 3.0;
      Double_t m2_min = -0.05;//Gev^2/c^4
      Double_t m2_max = 0.08; //Gev^2/c^4
      Float_t nSigmaPion_max_TPC = 2.0;
      Float_t nSigmaElectron_min_TPC = 2.0;
      Float_t nSigmaKaon_min_TPC = 2.0;
      Float_t nSigmaProton_min_TPC = 2.0;


      //QA before PID TPC or TPC+TOF check:
      hNSigmPion_vs_pPrimTotDevQ->Fill(PtotPrimQ, picoTrack->nSigmaPion());
      hNSigmKaon_vs_pPrimTotDevQ->Fill(PtotPrimQ, picoTrack->nSigmaKaon());
      hNSigmProton_vs_pPrimTotDevQ->Fill(PtotPrimQ, picoTrack->nSigmaProton());
      hNSigmElectron_vs_pPrimTotDevQ->Fill(PtotPrimQ, picoTrack->nSigmaElectron());
      hdEdx_vs_pPrimTotDevQ->Fill(PtotPrimQ, picoTrack->dEdx());
      
      Bool_t is_TPC_momenta_range = picoTrack->pMom().Mag()<p_tot_prim_mid_PID;
      Bool_t is_TPC_TOF_momenta_range = picoTrack->isTofTrack() && p_tot_prim_mid_PID<=picoTrack->pMom().Mag();
      //TPC only PID:
      if(is_TPC_momenta_range) //p_min already set by track choise
      {

        Bool_t is_nSigma_Pion_TPC = fabs(picoTrack->nSigmaPion())<nSigmaPion_max_TPC;
        Bool_t is_nSigma_Kaon_TPC = fabs(picoTrack->nSigmaKaon())>nSigmaKaon_min_TPC;
        Bool_t is_nSigma_Proton_TPC = fabs(picoTrack->nSigmaProton())>nSigmaProton_min_TPC;
        Bool_t is_nSigma_Electron_TPC = fabs(picoTrack->nSigmaElectron())>nSigmaElectron_min_TPC;
        if(is_nSigma_Pion_TPC && is_nSigma_Electron_TPC 
            && is_nSigma_Kaon_TPC && is_nSigma_Proton_TPC)
        {
          //QA histis filling after PID but after TPC only:
          hNSigmPion_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaPion());
          hNSigmKaon_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaKaon());
          hNSigmProton_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaProton());
          hNSigmElectron_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaElectron());
          hdEdx_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ,picoTrack->dEdx());
          //let's fill c++ vector of Pions after TPC only:
          Double_t m_Pion = 0.13957039;//GeV
          Double_t temp_pion_Energy_TPC_ONLY = sqrt(picoTrack->pMom().Mag2()+m_Pion*m_Pion);
          TLorentzVector temp_four_vector(picoTrack->pMom(),temp_pion_Energy_TPC_ONLY);

          Pions_4_momenta_Arr_ALL.push_back(temp_four_vector);
        }

      }//end of TPC only
      //TPC+TOF PID:
      else if( is_TPC_TOF_momenta_range ) //p_max already set by track choise
      {
      StPicoBTofPidTraits *trait = dst->btofPidTraits(picoTrack->bTofPidTraitsIndex()); // Retrieve corresponding trait
      if (!trait)
      {
        std::cout << "O-oh... No BTofPidTrait # " << picoTrack->bTofPidTraitsIndex()
                  << " for track # " << iTrk << std::endl;
        std::cout << "Check that you turned on the branch!" << std::endl;
        continue;
      }
      
      //variables for QA:
      Double_t m_square = picoTrack->pMom().Mag2()*(1./((trait->btofBeta())*(trait->btofBeta()))-1.);
      Double_t m_Pion = 0.13957039;//GeV
      Double_t one_beta_expect = sqrt(m_Pion*m_Pion+(picoTrack->pMom().Mag2()))/(picoTrack->pMom().Mag());

      //QA hists before PID but after TOF check (all that contains trait):
      h1_OverBeta_vs_pPrimTotDevQ->Fill(PtotPrimQ,1./(trait->btofBeta()));
      hm2_vs_pPrimTotDevQ->Fill(PtotPrimQ,m_square);
      h1_OverBetaDelta_vs_pPrimTotDevQ->Fill(PtotPrimQ, 1./(trait->btofBeta()) - one_beta_expect);
      
      Bool_t is_nSigma_Pion = fabs(picoTrack->nSigmaPion())<nSigmaPion_max_TOF;
      Bool_t is_1_beta_delta = fabs(1./(trait->btofBeta())-one_beta_expect)<one_over_beta_delta_max;
      Bool_t is_m2 = m2_min<m_square && m_square<m2_max;
      if(is_nSigma_Pion && is_1_beta_delta 
          && is_m2)
      {
        //QA hists filling after PID TPC+TOF cut:
        hNSigmPion_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaPion());
        hNSigmKaon_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaKaon());
        hNSigmProton_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaProton());
        hNSigmElectron_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaElectron());
        hdEdx_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ,picoTrack->dEdx());

        h1_OverBeta_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ,1./(trait->btofBeta()));
        h1_OverBetaDelta_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, 1./(trait->btofBeta()) - one_beta_expect);
        hm2_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ,m_square);
        //let's fill c++ vector of Pions after TPC & TOF:
        Double_t m_Pion = 0.13957039;//GeV
        Double_t temp_pion_Energy_TOF_TPC = sqrt(picoTrack->pMom().Mag2()+m_Pion*m_Pion);
        TLorentzVector temp_four_vector(picoTrack->pMom(),temp_pion_Energy_TOF_TPC);

        Pions_4_momenta_Arr_ALL.push_back(temp_four_vector);
      }//end of PID

    }//end of TOF + TPC

    }//end of track selection
    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    //now let's build A(q_inv) - Numerator of correlation function (Pions from one event):

    if(!Pions_4_momenta_Arr_ALL.empty())
    {
      Int_t N_of_Pions = Pions_4_momenta_Arr_ALL.size();
      for (Int_t i = 0; i < N_of_Pions; i++)
      {
        for (Int_t j = i+1; j < N_of_Pions; j++)
        {
          TLorentzVector delta_4_momenta = Pions_4_momenta_Arr_ALL[i]-Pions_4_momenta_Arr_ALL[j];
          double_t q_inv = sqrt(-delta_4_momenta.Mag2());

          hA_q_inv_ALL->Fill(q_inv);
        }
      }
    }

    //let's mix events:
    //queue structure: NEW element goes to BACK --- OLD elements pops out of FRONT

    // if(event->primaryVertex().Z()<20. && -20.<event->primaryVertex().Z() && event->refMult()>0. && event->refMult()<60.)
    // {
    //   //compare new vector of pions and all vectors of pions in qeue:
    //   comparePionsFillHistB(Pions_4_momenta_Arr_ALL,Pions_mix_queue_Arr_4_mom,hB_q_inv_ALL);
    //   Pions_mix_queue_Arr_4_mom.push_back(Pions_4_momenta_Arr_ALL);
    //   // clear buffer:
    //   if (Pions_mix_queue_Arr_4_mom.size() > BUFFER_SIZE)
    //   {
    //     Pions_mix_queue_Arr_4_mom.pop_front(); // delete the oldest element
    //   }
    // }

    //cycle for each cut by Vz & RefMult:
    for(int iVz=0;iVz<nVzCuts;iVz++)
    {
      for(int iRefM=0;iRefM<nRefMultCuts;iRefM++)
      {
        Bool_t is_prim_Vz_cut = VzBins[iVz]<event->primaryVertex().Z() && event->primaryVertex().Z()<VzBins[iVz+1];
        Bool_t is_Ref_Mult_cut = RefMultBins[iRefM]<event->refMult() && event->refMult()<RefMultBins[iRefM+1];
        if(is_prim_Vz_cut && is_Ref_Mult_cut)
        {
          //compare new vector of pions and all vectors of pions in qeue:
          comparePionsFillHistB(Pions_4_momenta_Arr_ALL,Pions_Buffer[iVz][iRefM],hB_q_inv_ALL);
          Pions_Buffer[iVz][iRefM].push_back(Pions_4_momenta_Arr_ALL);
          // clear buffer:
          if(Pions_Buffer[iVz][iRefM].size() > BUFFER_SIZE)
          {
            Pions_Buffer[iVz][iRefM].pop_front(); // delete the oldest element
          }

        }
      }
    }

    }//end of event selection
  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  picoReader->Finish();
  oFile->Write();
  oFile->Close();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	    << std::endl;
}

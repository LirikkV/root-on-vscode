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

#include <bitset>

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
#include "TStopwatch.h"

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

//It's better use ROOT::Math::LorentzVector instead of TLorenzVector
const Double_t m_Pion = 0.13957039;//GeV
using My_LorenzVector = ROOT::Math::PxPyPzEVector;

struct My_ParticleTrackInfo{
  My_LorenzVector p4;
  UInt_t topologyMap0;
  UInt_t topologyMap1;
  ULong64_t iTpcTopologyMap;
  Int_t Nhits;
  StPicoPhysicalHelix helix;
};
//Spliting level: 
const double maxSplitLevel = 0.7;
//+1 if 01 or 10 at the same bit positions of two track hit info
//-1 if 11
// 0 if 00
// __builtin_popcount(a) gives number of 1 in a in binary representation
// 7 in binary is 111 =>__builtin_popcount(7)=3
//31 = 1111 => __builtin_popcount(31)=5
//if we compare 100 and 110, we get 0=(-1)+(+1)+(0) because it's pairs 11 01 00
//so if we have: int a b
//then xor(a,b) return +1 if 10 or 01 
//at c++ xor(a,b) is a^b
//__builtin_popcount(a^b) returnes number of 10 or 01 pairs
//__builtin_popcount(a&b) returnes number of 11 pairs
//00 pairs don't effect on sum

//we need to start this sum from 8-th bit
//0xFFFFFF00 = 11111111111111111111111100000000 - first bits from 0 to 7 are zeroes
//0x1FFFFF =              111111111111111111111 - only first bits from 0 to 20 are not zeroes 
//a&0xFFFFFF00 sets first 8 bits of a to 0
double getSplitLevel(const My_ParticleTrackInfo& tr_1,const My_ParticleTrackInfo& tr_2)
{
  int minusOnes = __builtin_popcount( (tr_1.topologyMap0&tr_2.topologyMap0) & 0xFFFFFF00)+
                     __builtin_popcount( (tr_1.topologyMap1&tr_2.topologyMap1) & 0x1FFFFF);
  int sumNhits = tr_1.Nhits+tr_2.Nhits;
  return(static_cast<double>(sumNhits - 3 * minusOnes) / sumNhits);
}

//_________________
//functions for fraction of merged rows:
const double FMR_max = 0.7;
int TpcLocalTransform( TVector3& aPoint, int& aSector, int& aRow, float& aU, double& aPhi);
//_________________
void calculateTpcExitAndEntrancePoints(StPicoPhysicalHelix tHelix,TVector3 PrimVert,TVector3 SecVert,TVector3 tmpTpcEntrancePoint,TVector3 tmpTpcExitPoint,TVector3* tmpPosSample,float* tmpZ,float* tmpU,int* tmpSect);
//_________________
double fractionOfMergedRow(StPicoPhysicalHelix helixTrk1, StPicoPhysicalHelix helixTrk2);
//_________________

void fill_A_qinv(const std::vector<My_ParticleTrackInfo>& Pions_4_momenta_hits_Arr, TH1D* hist_A, TH2F* hSL_A)
{
  if(!Pions_4_momenta_hits_Arr.empty())
    {
      Int_t N_of_Pions = Pions_4_momenta_hits_Arr.size();

      for (Int_t i = 0; i < N_of_Pions; i++)
      {
          const double px =  Pions_4_momenta_hits_Arr[i].p4.Px();
          const double py = Pions_4_momenta_hits_Arr[i].p4.Py();
          const double pz = Pions_4_momenta_hits_Arr[i].p4.Pz();
          const double e = Pions_4_momenta_hits_Arr[i].p4.E();
        for (Int_t j = i+1; j < N_of_Pions; j++)
        {
          double SplitLevel = getSplitLevel(Pions_4_momenta_hits_Arr[j],Pions_4_momenta_hits_Arr[i]);
          if(SplitLevel>maxSplitLevel)continue; //if SplitLevel is wery high - we stop iteration and get to next step over j
          // double FMR = fractionOfMergedRow(Pions_4_momenta_hits_Arr[j].helix,Pions_4_momenta_hits_Arr[i].helix);
          // if(FMR>FMR_max) continue;

          const double dpx = px - Pions_4_momenta_hits_Arr[j].p4.Px();
          const double dpy = py - Pions_4_momenta_hits_Arr[j].p4.Py();
          const double dpz = pz - Pions_4_momenta_hits_Arr[j].p4.Pz();
          const double de = e - Pions_4_momenta_hits_Arr[j].p4.E();

          // //std::bitset<32>(int_number).to_string() - gets string in binary form
          // std::cout<<Pions_4_momenta_hits_Arr[j].trHits_1 <<"  "
          //          <<Pions_4_momenta_hits_Arr[j].trHits_2<<std::endl;
          // std::cout<<Pions_4_momenta_hits_Arr[i].trHits_1<<"  "
          //          <<Pions_4_momenta_hits_Arr[i].trHits_2<<std::endl;
          // std::cout<<std::bitset<32>(Pions_4_momenta_hits_Arr[j].trHits_1).to_string() <<"  "
          //          <<std::bitset<32>(Pions_4_momenta_hits_Arr[j].trHits_2).to_string()
          //          <<std::endl;
          // std::cout<<std::bitset<32>(Pions_4_momenta_hits_Arr[i].trHits_1).to_string() <<"  "
          //          <<std::bitset<32>(Pions_4_momenta_hits_Arr[i].trHits_2).to_string()
          //          <<std::endl;
          // std::cout<<getSplitLevel(Pions_4_momenta_hits_Arr[j],Pions_4_momenta_hits_Arr[i])<<std::endl;
          // std::cout<<Pions_4_momenta_hits_Arr[j].Nhits<<" "<< Pions_4_momenta_hits_Arr[i].Nhits<<std::endl;

          double q_inv_2 = dpx*dpx+dpy*dpy+dpz*dpz-de*de;
          if(q_inv_2 > 0.)
          {
            double q_inv = sqrt(q_inv_2);
            hSL_A->Fill(q_inv, SplitLevel);
            hist_A->Fill(q_inv);
          }
        }
      }
    }
}
//this function compare pions from "new" event with events from buffer and fills Numerator of CF
void comparePionsFillHistB(const std::vector<My_ParticleTrackInfo>& new_Pions_Arr,
                           const std::deque<std::vector<My_ParticleTrackInfo>>& event_Queue,
                           TH1D* hist_B, TH2F* hSL_B)
{
  const size_t nNewPions = new_Pions_Arr.size();
  const size_t nEventsInQueue = event_Queue.size();

  //loop over queue vectors of pions:
  for(size_t i=0;i<nEventsInQueue;i++)
  {
    const std::vector<My_ParticleTrackInfo>& queue_Pions_Arr = event_Queue[i];
    const size_t nQueuePions = queue_Pions_Arr.size();

    //loop over pions from new event
    for(size_t j=0;j<nNewPions;j++)
    {
      const double px1 = new_Pions_Arr[j].p4.Px();
      const double py1 = new_Pions_Arr[j].p4.Py();
      const double pz1 = new_Pions_Arr[j].p4.Pz();
      const double e1 = new_Pions_Arr[j].p4.E();

      //loop over pions from selected queue event:
      for(size_t k=0;k<nQueuePions;k++)
      {
        double SplitLevel = getSplitLevel(new_Pions_Arr[j],queue_Pions_Arr[k]);
        if(SplitLevel>maxSplitLevel)continue;
        // double FMR = fractionOfMergedRow(new_Pions_Arr[j].helix, queue_Pions_Arr[k].helix);
        // if(FMR>FMR_max) continue;

        const double dpx = px1-queue_Pions_Arr[k].p4.Px();
        const double dpy = py1-queue_Pions_Arr[k].p4.Py();
        const double dpz = pz1-queue_Pions_Arr[k].p4.Pz();
        const double de = e1-queue_Pions_Arr[k].p4.E();

        const double q_inv_2 = dpx*dpx+dpy*dpy+dpz*dpz-de*de;
        if(q_inv_2>0.)
        {
        double_t q_inv = sqrt(q_inv_2);
        hSL_B->Fill(q_inv, SplitLevel);
        hist_B->Fill(q_inv);
        }

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

  //let's see how much time program get:
  TStopwatch timer;

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
			    200,-3.5,3.5,200,-3.5,3.5);
  TH2F *hVtxXvsY_cut = new TH2F("hVtxXvsY_cut",
			    "Vtx XvsY after event cut",
			    200,-3.5,3.5,200,-3.5,3.5);

  TH1F *hVtxZ = new TH1F("hVtxZ","Vtx Z",
			 140, -70., 70.);
  TH1F *hVtxZ_cut = new TH1F("hVtxZ_cut","Vtx Z after event cut",
			 140, -70., 70.);

  // Track
  // TH1F *hGlobalPtot = new TH1F("hGlobalPtot",
	// 		       "Global track momentum;p (GeV/c)",
	// 		       100, 0., 1. );
  // TH1F *hGlobalPtotCut = new TH1F("hGlobalPtotCut",
	// 			  "Global track momentum after cut;p (GeV/c)",
	// 			  100, 0., 1. );
  TH1F *hPrimaryPtot = new TH1F("hPrimaryPtot",
				"Primary track momentum;p (GeV/c)",
			       100, 0., 2. );
  TH1F *hPrimaryPtot_cut = new TH1F("hPrimaryPtot_cut",
				   "Primary track momentum after track cut;p (GeV/c)",
				  100, 0., 2. );

  // TH1F *hTransvMomentum = new TH1F("hTransvMomentum",
	// 			   "Track transverse momentum;p_{T} (GeV/c)",
	// 			   200, 0., 2.);

  // TH2F *hGlobalPhiVsPt[2];
  // for(int i=0; i<2; i++) {
  //   hGlobalPhiVsPt[i] = new TH2F(Form("hGlobalPhiVsPt_%d",i),
	// 			 Form("#phi vs. p_{T} for charge: %d;p_{T} (GeV/c);#phi (rad)", (i==0) ? 1 : -1),
	// 			 300, 0., 3.,
	// 			 630, -3.15, 3.15);
  // }     
      
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
			    200,-2.,2.,200,-40.,40.);
  TH2F *hNSigmPion_vs_pPrimTotDevQ_cut_PID = new TH2F("hNSigmPion_vs_pPrimTotDevQ_cut_PID",
			    "nSigma(pion) vs P_prim_tot/q after cut;;nSigma",
			    200,-2.,2.,200,-40.,40.);

  TH2F *hNSigmKaon_vs_pPrimTotDevQ = new TH2F("hNSigmKaon_vs_pPrimTotDevQ",
			    "nSigma(Kaon) vs P_prim_tot/q;;nSigma",
			    200,-2.,2.,200,-40.,40.);
  TH2F *hNSigmKaon_vs_pPrimTotDevQ_cut_PID = new TH2F("hNSigmKaon_vs_pPrimTotDevQ_cut_PID",
			    "nSigma(Kaon) vs P_prim_tot/q after cut;;nSigma",
			    200,-2.,2.,200,-40.,40.);
  
  TH2F *hNSigmProton_vs_pPrimTotDevQ = new TH2F("hNSigmProton_vs_pPrimTotDevQ",
			    "nSigma(Proton) vs P_prim_tot/q;;nSigma",
			    200,-2.,2.,200,-40.,40.);
  TH2F *hNSigmProton_vs_pPrimTotDevQ_cut_PID = new TH2F("hNSigmProton_vs_pPrimTotDevQ_cut_PID",
			    "nSigma(Proton) vs P_prim_tot/q after cut;;nSigma",
			    200,-2.,2.,200,-40.,40.);
  
  TH2F *hNSigmElectron_vs_pPrimTotDevQ = new TH2F("hNSigmElectron_vs_pPrimTotDevQ",
			    "nSigma(Electron) vs P_prim_tot/q;;nSigma",
			    200,-2.,2.,200,-40.,40.);
  TH2F *hNSigmElectron_vs_pPrimTotDevQ_cut_PID = new TH2F("hNSigmElectron_vs_pPrimTotDevQ_cut_PID",
			    "nSigma(Electron) vs P_prim_tot/q after cut;;nSigma",
			    200,-2.,2.,200,-40.,40.);

  TH2F *h1_OverBeta_vs_pPrimTotDevQ = new TH2F("1_OverBeta_vs_pPrimTotDevQ",
			    "1/beta vs P_prim_tot/q;;1/beta",
			    200,-2.,2.,200,0.5,2.0);
  TH2F *h1_OverBeta_vs_pPrimTotDevQ_cut_PID = new TH2F("1_OverBeta_vs_pPrimTotDevQ_cut_PID",
			    "1/beta vs P_prim_tot/q after PID;;1/beta",
			    400,-2.,2.,200,0.5,2.0);
  
  TH2F *hm2_vs_pPrimTotDevQ = new TH2F("hm2_vs_pPrimTotDevQ",
			    "m^2 vs P_prim_tot/q;;m^2(Gev/c)",
			    400,-2.,2.,200,-0.1,0.1);
  TH2F *hm2_vs_pPrimTotDevQ_cut_PID = new TH2F("hm2_vs_pPrimTotDevQ_cut_PID",
			    "m^2 vs P_prim_tot/q after PID;;m^2(Gev/c)",
			    400,-2.,2.,200,-0.1,0.1);

  TH2F *hTEST_P2_pPrimTotDevQ = new TH2F("hTEST_P2_pPrimTotDevQ",
			    "P_tot^2 vs P_prim_tot/q;;m^2(Gev/c)",
			    400,-2.,2.,200, 0.0,2.5);
  TH2F *hTEST_Beta_pPrimTotDevQ = new TH2F("hTEST_sqrt_Beta_pPrimTotDevQ",
			    "(1/Beta^2-1) vs P_prim_tot/q;;m^2(Gev/c)",
			    400,-2.,2.,200,-0.1,0.1);
  
  TH2F *h1_OverBetaDelta_vs_pPrimTotDevQ = new TH2F("1_OverBetaDelta_vs_pPrimTotDevQ",
			    "1/beta - 1/beta_exp vs P_prim_tot/q;;1/beta - 1/beta_exp",
			    200,-2.,2.,200,-0.2,1.2);
  TH2F *h1_OverBetaDelta_vs_pPrimTotDevQ_cut_PID = new TH2F("1_OverBetaDelta_vs_pPrimTotDevQ_cut_PID",
			    "1/beta - 1/beta_exp vs P_prim_tot/q after PID;;1/beta - 1/beta_exp",
			    400,-2.,2.,200,-0.2,1.2);

  TH2F *hdEdx_vs_pPrimTotDevQ = new TH2F("hdEdx_vs_pPrimTotDevQ",
			    "dE/dx vs P_prim_tot/q;;dE/dx(keV/cm)",
			    200,-1.6,1.6,200,-0.2,10.);
  TH2F *hdEdx_vs_pPrimTotDevQ_cut_PID = new TH2F("hdEdx_vs_pPrimTotDevQ_cut_PID",
			    "dE/dx vs P_prim_tot/q after PID;;dE/dx(keV/cm)",
			    200,-1.6,1.6,200,-0.2,10. );
  
  TH2F *h2DpPrimTr_vs_etaPtim = new TH2F("h2DpPrimTr_vs_etaPtim",
			    "p_prim_T vs eta_prim;p_T;eta",
			    200,-0.5,2.,200,-4.,4.);
  TH2F *h2DpPrimTr_vs_etaPtim_cut = new TH2F("h2DpPrimTr_vs_etaPtim_cut",
			    "p_prim_T vs eta_prim after track cut;p_T;eta",
			    200,-0.5,2.,200,-4.,4.);
        
  TH2F *hTEST2DpPrimTr_vs_etaPtim_equal_P = new TH2F("hTEST2DpPrimTr_vs_etaPtim_equal_P",
			    "p_prim_T vs eta_prim after track cut with equal P;p_T;eta",
			    200,-0.5,2.,200,-4.,4.);

  // BTof pid traits
  TH1F *hTofBeta = new TH1F("hTofBeta", "BTofPidTraits #beta;#beta",
			    2000, 0., 2.);

  //Correlation function:

  TH1D *hA_Pi_Plus_q_inv_ALL = new TH1D("hA_Pi_Plus_q_inv_ALL",
				   "Numerator of Corr.Funct Pi+ Pi+ with both TPC & TPC+TOF methods;q_inv;A",
				  375, 0., 3.0 );
  TH1D *hB_Pi_Plus_q_inv_ALL = new TH1D("hB_Pi_Plus_q_inv_ALL",
				   "Denumerator of Corr.Funct Pi+ Pi+ with both TPC & TPC+TOF methods;q_inv;B",
				  375, 0., 3.0 );

  TH1D *hA_Pi_Minus_q_inv_ALL = new TH1D("hA_Pi_Minus_q_inv_ALL",
				   "Numerator of Corr.Funct Pi- Pi- with both TPC & TPC+TOF methods;q_inv;A",
				  375, 0., 3.0 );
  TH1D *hB_Pi_Minus_q_inv_ALL = new TH1D("hB_Pi_Minus_q_inv_ALL",
				   "Denumerator of Corr.Funct Pi- Pi- with both TPC & TPC+TOF methods;q_inv;B",
				  375, 0., 3.0 );
  //Splittiong level:
  TH2F *hSL_A = new TH2F("SL_vs_qinv_A",
			    "Split level vs q_inv A;q_inv;SL",
			    400,-0.,0.2,200,-0.5,1.);
  TH2F *hSL_B = new TH2F("SL_vs_qinv_B",
			    "Split level vs q_inv B;q_inv;SL",
			    400,-0.,0.2,200,-0.5,1.);
  //cuts by Vz: 4 cats; Vz from -40 to 40
  //cuts by refMult: 10 cuts; RefMult from 0 to 600
  const Int_t nVzCuts = 4;
  const Int_t nRefMultCuts = 10;
  const Double_t VzBins[nVzCuts+1] = {-40., -20., 0., 20., 40.};
  const Double_t RefMultBins[nRefMultCuts+1] = {0.,60.,120.,180.,240.,300.,360.,420.,480.,540.,600};

  //for mixing events:
  const Int_t BUFFER_SIZE = 5;
  //this is queue from events; just queue from vectors from 4-vectors and hits information of particle's track
  std::vector<std::vector<std::deque<std::vector<My_ParticleTrackInfo>>>> Pions_Plus_Buffer(
                                          nVzCuts, 
                                          std::vector<std::deque<std::vector<My_ParticleTrackInfo>>>(nRefMultCuts));
                                          //this is 4*10 queues; each queue consists from events; just 2D vector of previous queues for each cut
  std::vector<std::vector<std::deque<std::vector<My_ParticleTrackInfo>>>> Pions_Minus_Buffer(
                                          nVzCuts, 
                                          std::vector<std::deque<std::vector<My_ParticleTrackInfo>>>(nRefMultCuts));
  

  //let's create a c++ vector with 4-momenta of pions in one event:
  std::vector<My_ParticleTrackInfo> Pions_Plus_4_momenta_hits_Arr_ALL;
  std::vector<My_ParticleTrackInfo> Pions_Minus_4_momenta_hits_Arr_ALL;

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

    //at each step we will clear vectors with 4-momenta of pions:
    Pions_Plus_4_momenta_hits_Arr_ALL.clear();
    Pions_Minus_4_momenta_hits_Arr_ALL.clear();

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
    hVtxXvsY->Fill( pVtx.X(), pVtx.Y() );
    hVtxZ->Fill( pVtx.Z() );

    // Lirikk's cuts:
    // variables for event cut:
    double_t Vtx_r_Max = 2.0;  // cm
    double_t Vtx_z_Max = 40.0; // cm

    // Event selection:
    Bool_t is_Vtx_r_cut = pVtx.Perp() < Vtx_r_Max;
    Bool_t is_abs_Vtx_z_cut = fabs(pVtx.Z()) < Vtx_z_Max;
    if (is_Vtx_r_cut && is_abs_Vtx_z_cut)
    {
      //QA hists after event cut:
      hVtxXvsY_cut->Fill(pVtx.X(), pVtx.Y());
      hVtxZ_cut->Fill(pVtx.Z());
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
    if(!is_N_TPC_fit_hits_cut) continue;
    Bool_t is_DCA_abs_cut = picoTrack->gDCA(pVtx).Mag()<DCA_max;
    if(!is_DCA_abs_cut) continue;
    Bool_t is_pseudo_prm_cut = picoTrack->isPrimary() &&
                            fabs(picoTrack->pMom().Eta())<pseudo_rap_prim_max;
    if(!is_pseudo_prm_cut) continue;
    Bool_t is_p_tot_prim_cut = picoTrack->isPrimary() && 
                            p_tot_prim_min < picoTrack->pMom().Mag() &&
                            picoTrack->pMom().Mag()<p_tot_prim_max;
    if(!is_p_tot_prim_cut) continue;
    Bool_t is_p_trans_prim_cut = picoTrack->isPrimary() && 
                            p_trans_prim_min<picoTrack->pMom().Pt() &&
                            picoTrack->pMom().Pt()<p_trans_prim_max;
    if(is_p_trans_prim_cut) 
    {
      //QA after track selection:
      hNFitHits_cut->Fill(picoTrack->nHitsFit());
	    hPrimaryPtot_cut->Fill( picoTrack->pMom().Mag() );
      hDCA_cut->Fill(picoTrack->gDCA(pVtx).Mag());
      hPrimaryPtrans_cut->Fill(picoTrack->pMom().Pt());
      hPrimaryPseudorap_cut->Fill(picoTrack->pMom().Eta());
      h2DpPrimTr_vs_etaPtim_cut->Fill(picoTrack->pMom().Pt(),picoTrack->pMom().Eta());

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
        if(!is_nSigma_Pion_TPC) continue;
        Bool_t is_nSigma_Kaon_TPC = fabs(picoTrack->nSigmaKaon())>nSigmaKaon_min_TPC;
        if(!is_nSigma_Kaon_TPC) continue;
        Bool_t is_nSigma_Proton_TPC = fabs(picoTrack->nSigmaProton())>nSigmaProton_min_TPC;
        if(!is_nSigma_Proton_TPC) continue;
        Bool_t is_nSigma_Electron_TPC = fabs(picoTrack->nSigmaElectron())>nSigmaElectron_min_TPC;
        if(is_nSigma_Electron_TPC)
        {
          //QA histis filling after PID but after TPC only:
          hNSigmPion_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaPion());
          hNSigmKaon_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaKaon());
          hNSigmProton_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaProton());
          hNSigmElectron_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ, picoTrack->nSigmaElectron());
          hdEdx_vs_pPrimTotDevQ_cut_PID->Fill(PtotPrimQ,picoTrack->dEdx());
          //let's fill c++ vector of Pions after TPC only:
          Double_t temp_pion_Energy_TPC_ONLY = sqrt(picoTrack->pMom().Mag2()+m_Pion*m_Pion);
          My_LorenzVector temp_four_vector(picoTrack->pMom().Px(), picoTrack->pMom().Py(), 
                              picoTrack->pMom().Pz(), temp_pion_Energy_TPC_ONLY);

          My_ParticleTrackInfo temp_MyParticle;
          temp_MyParticle.p4 = temp_four_vector;
          temp_MyParticle.topologyMap0 = picoTrack->topologyMap(0);
          temp_MyParticle.topologyMap1 = picoTrack->topologyMap(1);
          temp_MyParticle.iTpcTopologyMap = picoTrack->iTpcTopologyMap();
          temp_MyParticle.Nhits = picoTrack->nHits();
          temp_MyParticle.helix = picoTrack->helix( event->bField() );

          //separation to Pi+Pi+ & Pi-Pi- pairs
          if(picoTrack->charge()>0.)
          {
            //let's fill
            Pions_Plus_4_momenta_hits_Arr_ALL.push_back(temp_MyParticle);
          }
          else if(picoTrack->charge()<0.)
          {
            Pions_Minus_4_momenta_hits_Arr_ALL.push_back(temp_MyParticle);
          }
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
      
      hTEST_Beta_pPrimTotDevQ->Fill(PtotPrimQ,(1./((trait->btofBeta())*(trait->btofBeta()))-1.));
      hTEST_P2_pPrimTotDevQ->Fill(PtotPrimQ, picoTrack->pMom().Mag2());

      Bool_t is_nSigma_Pion = fabs(picoTrack->nSigmaPion())<nSigmaPion_max_TOF;
      if(!is_nSigma_Pion) continue;
      Bool_t is_1_beta_delta = fabs(1./(trait->btofBeta())-one_beta_expect)<one_over_beta_delta_max;
      if(!is_1_beta_delta) continue;
      Bool_t is_m2 = m2_min<m_square && m_square<m2_max;
      if(is_m2)
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
        
        //let's fill c++ vector of Pions after TPC & TOF PID:
        Double_t temp_pion_Energy_TOF_TPC = sqrt(picoTrack->pMom().Mag2()+m_Pion*m_Pion);
        My_LorenzVector temp_four_vector(picoTrack->pMom().Px(), picoTrack->pMom().Py(),
                              picoTrack->pMom().Pz(), temp_pion_Energy_TOF_TPC);

          My_ParticleTrackInfo temp_MyParticle;
          temp_MyParticle.p4 = temp_four_vector;
          temp_MyParticle.topologyMap0 = picoTrack->topologyMap(0);
          temp_MyParticle.topologyMap1 = picoTrack->topologyMap(1);
          temp_MyParticle.iTpcTopologyMap = picoTrack->iTpcTopologyMap();
          temp_MyParticle.Nhits = picoTrack->nHits();
          temp_MyParticle.helix = picoTrack->helix( event->bField() );

        //separation to Pi+Pi+ & Pi-Pi- pairs
        if (picoTrack->charge() > 0.)
        {
          Pions_Plus_4_momenta_hits_Arr_ALL.push_back(temp_MyParticle);
        }
        else if (picoTrack->charge() < 0.)
        {
          Pions_Minus_4_momenta_hits_Arr_ALL.push_back(temp_MyParticle);
        }
      }//end of PID

    }//end of TOF + TPC

    }//end of track selection
    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    //now let's build A(q_inv) - Numerator of correlation function (Pions from one event):
    //for Pi+Pi+ & Pi-Pi- pairs:
    fill_A_qinv(Pions_Plus_4_momenta_hits_Arr_ALL,hA_Pi_Plus_q_inv_ALL,hSL_A);
    fill_A_qinv(Pions_Minus_4_momenta_hits_Arr_ALL,hA_Pi_Minus_q_inv_ALL,hSL_A);


    //let's mix events:
    //queue structure: NEW element goes to BACK --- OLD elements pops out of FRONT

    
    //instead of cycles better use BinarySearch to search for interval our event belongs to

    Int_t iVz = TMath::BinarySearch(nVzCuts + 1, VzBins, (Double_t)pVtx.Z());
    Int_t iRefM = TMath::BinarySearch(nRefMultCuts +1, RefMultBins, (Double_t)event->refMult());

    //check indexies for preventing segmentation fault:
    if(iVz>=0 && iVz<nVzCuts && iRefM>=0 && iRefM<nRefMultCuts)
    {

      // compare new vector of pions and all vectors of pions in qeue:
      // for Pi+Pi+ & Pi-Pi- pions:
      comparePionsFillHistB(Pions_Plus_4_momenta_hits_Arr_ALL, Pions_Plus_Buffer[iVz][iRefM], hB_Pi_Plus_q_inv_ALL,hSL_B);
      comparePionsFillHistB(Pions_Minus_4_momenta_hits_Arr_ALL, Pions_Minus_Buffer[iVz][iRefM], hB_Pi_Minus_q_inv_ALL,hSL_B);

      Pions_Plus_Buffer[iVz][iRefM].push_back(Pions_Plus_4_momenta_hits_Arr_ALL);
      Pions_Minus_Buffer[iVz][iRefM].push_back(Pions_Minus_4_momenta_hits_Arr_ALL);

      // clear buffer:
      if (Pions_Plus_Buffer[iVz][iRefM].size() > BUFFER_SIZE)
      {
        Pions_Plus_Buffer[iVz][iRefM].pop_front(); // delete the oldest element
      }
      if (Pions_Minus_Buffer[iVz][iRefM].size() > BUFFER_SIZE)
      {
        Pions_Minus_Buffer[iVz][iRefM].pop_front(); // delete the oldest element
      }
    }

    }//end of event selection
  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)


  //Some beautify of histograms:
  Double_t max_hNFitHits = hNFitHits->GetMaximum();
  hNFitHits_cut->SetMaximum(max_hNFitHits);

  picoReader->Finish();
  oFile->Write();
  oFile->Close();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	    << std::endl;
  
  timer.Stop();
  std::cout <<"Real time:"<< timer.RealTime()<<std::endl;
  std::cout <<"CPU time:"<< timer.CpuTime()<<std::endl;
}



//here TPC topology and some useful functions:
//_________________
float mInnerTpcRadius  = 50.f;   //[cm]
float mOuterTpcRadius = 200.f;   //[cm]
float mTpcHalfLength  = 200.f;   //[cm]

/// Number of points where perform calculations
static const unsigned short mNumberOfPoints = 11;
static const unsigned short mNumberOfPadrows = 45;
static const int mNumberOfPadrowsInner = 13;
static const int mNumberOfSectors = 24;
/// Radius at each padrow
static const float tRowRadius[mNumberOfPadrows] = { 60, 64.8,69.6,74.4,79.2,84,88.8,93.6,98.8,
																								      104,109.2,114.4,119.6,127.195,129.195,131.195,
																								      133.195,135.195,137.195,139.195,141.195,
																								      143.195,145.195,147.195,149.195,151.195,
																								      153.195,155.195,157.195,159.195,161.195,
																								      163.195,165.195,167.195,169.195,171.195,
																								      173.195,175.195,177.195,179.195,181.195,
																								      183.195,185.195,187.195,189.195 };

//_________________
int TpcLocalTransform( TVector3& aPoint, int& aSector, int& aRow,
		       						 float& aU, double& aPhi) {

  // Number of pad in a row for each padrow
  static int tNPadAtRow[mNumberOfPadrows] = {
    88,96,104,112,118,126,134,142,150,158,166,174,182,
    98,100,102,104,106,106,108,110,112,112,114,116,118,120,122,122,
    124,126,128,128,130,132,134,136,138,138,140,142,144,144,144,144
  };
  // Phi angle of each of 24 sectors
  static double tSectToPhi[mNumberOfSectors] = {
    2., 1., 0., 11., 10., 9., 8., 7., 6., 5., 4., 3.,
    4., 5., 6., 7., 8., 9., 10., 11., 0., 1., 2., 3.
  };

  // Pad size for innner and outer sector
  static double tPadWidthInner = 0.335;
  static double tPadWidthOuter = 0.67;

  static double tPi = TMath::Pi();

  // Find sector number
  aPhi = aPoint.Phi();
  if( aPhi<0. ) {
    aPhi+=(2*tPi);
  }
  aPhi += tPi/12.;

  if( aPhi>2*tPi ) {
    aPhi-=2*tPi;
  }

  int tiPhi = (int) (aPhi/tPi*6.);
  if( aPoint.Z() < 0 ) {
    aSector = (tiPhi<3) ? 3-tiPhi : 15-tiPhi;
  }
  else{
    aSector = (tiPhi<4) ? 21+tiPhi : 9+tiPhi;
  }
  aPhi = tSectToPhi[aSector-1] * tPi/6.;

  //if((fabs(aPhi-aPoint.phi())>(tPi/12)){
  //cout << "Sector missmatch " << aPhi << " " << aPoint.phi() << " "
  // << aSector << endl;
  //}

  // Calculate local coordinate
  float tR = aPoint.X() * TMath::Cos(aPhi) + aPoint.Y() * TMath::Sin(aPhi);
  aU =      -aPoint.X() * TMath::Sin(aPhi) + aPoint.Y() * TMath::Cos(aPhi);

  // Find pad row
  if(tR<57.6) {
    aRow = 0;
    return 1;
  }
  float radmax = 62.4;
  float spacing= 4.8;
  aRow=1;
  while ( ( tR > radmax) && ( aRow < 46 ) ) {
    aRow++;
    if ( aRow == 8 ) {
      radmax = 96.2;
      spacing = 5.2;
    }
    else {
      // Row 13 is the last one of the inner sector
      if ( aRow == mNumberOfPadrowsInner ) {
				radmax = 126.195;
				spacing = 2.0;
      }
      else {
				radmax += spacing;
      }
    } //else
  } //while ( ( tR > radmax) && ( aRow < 46 ) )
  if ( aRow > mNumberOfPadrows ) {
    //cout << "No pad row " << tR << endl;
    return 2;
  }

  // Check if u (=aU) inbound
  double tPadWidth = (aRow < 14 ) ? tPadWidthInner : tPadWidthOuter;
  if( TMath::Abs(aU) > ( tNPadAtRow[aRow-1] * tPadWidth/2. ) ) {
    return 3;
  }

  return 0;
}



//_________________
void calculateTpcExitAndEntrancePoints(StPicoPhysicalHelix tHelix,
                                       TVector3 PrimVert,
                                       TVector3 SecVert,
                                       TVector3 tmpTpcEntrancePoint,
                                       TVector3 tmpTpcExitPoint,
                                       TVector3 *tmpPosSample,
                                       float* tmpZ,
                                       float* tmpU,
                                       int* tmpSect) {

  // This calculates the exit point of a secondary track,
  // either through the endcap or through the Outer Field Cage
  // We assume the track to start at tHelix.origin-PrimaryVertex
  // it also calculates the entrance point of the secondary track,
  // which is the point at which it crosses the inner field cage.
  TVector3 ZeroVec(0.,0.,0.);

  ZeroVec.SetX(SecVert.x()-PrimVert.x());
  ZeroVec.SetY(SecVert.y()-PrimVert.y());
  ZeroVec.SetZ(SecVert.z()-PrimVert.z());

  double dip, curv, phase;
  int h;
  curv = tHelix.curvature();
  dip  = tHelix.dipAngle();
  phase= tHelix.phase();
  h    = tHelix.h();

  StPicoHelix hel(curv, dip, phase, ZeroVec, h);

  std::pair< double, double > candidates;
  // This is how much length to go to leave through sides of TPC
  double sideLength;
  // This is how much length to go to leave through endcap of TPC
  double endLength;

  // Figure out how far to go to leave through side...
  candidates = hel.pathLength( mTpcHalfLength );
  sideLength = (candidates.first > 0) ? candidates.first : candidates.second;

  static TVector3 WestEnd( 0., 0., mTpcHalfLength  );
  static TVector3 EastEnd( 0., 0., -mTpcHalfLength );
  static TVector3 EndCapNormal( 0., 0., 1.0);

  endLength = hel.pathLength( WestEnd, EndCapNormal);
  if (endLength < 0.0) {
    endLength = hel.pathLength( EastEnd, EndCapNormal);
  }
  if (endLength < 0.0) {
    std::cout << "void StHbtParticle::calculateTpcExitAndEntrancePoints : "
	      << "Hey -- I cannot find an exit point out endcaps" << std::endl;
  }

  // OK, firstExitLength will be the shortest way out of the detector...
  double firstExitLength = (endLength < sideLength) ? endLength : sideLength;

  // Now then, let's return the POSITION at which particle leaves TPC...
  tmpTpcExitPoint = hel.at( firstExitLength );

  // Finally, calculate the position at which the track crosses the inner field cage
  candidates = hel.pathLength( mInnerTpcRadius );

  sideLength = (candidates.first > 0) ? candidates.first : candidates.second;

  tmpTpcEntrancePoint = hel.at(sideLength);

  // Check that the entrance point exists and was found
  if (std::isnan( tmpTpcEntrancePoint.X() ) ||
      std::isnan( tmpTpcEntrancePoint.Y() ) ||
      std::isnan( tmpTpcEntrancePoint.Z() ) ) {

    std::cout << "void StHbtParticle::calculateTpcExitAndEntrancePoints : NAN" << std::endl;
    std::cout << "tmpNominalTpcEntrancePoint = ( "
	      << tmpTpcEntrancePoint.X() << " , "
	      << tmpTpcEntrancePoint.Y() << " , "
	      << tmpTpcEntrancePoint.Z() << " ) " << std::endl;
    tmpTpcEntrancePoint.SetX(-9999.);
    tmpTpcEntrancePoint.SetY(-9999.);
    tmpTpcEntrancePoint.SetZ(-9999.);
  }

  // Check that the exit point exists and was found
  if ( std::isnan( tmpTpcExitPoint.X() ) ||
       std::isnan( tmpTpcExitPoint.Y() ) ||
       std::isnan( tmpTpcExitPoint.Z() ) ) {

    std::cout << "void StHbtParticle::calculateTpcExitAndEntrancePoints : NAN" << std::endl;
    std::cout << "tmpNominalTpcExitPoint = ( "
	      << tmpTpcExitPoint.X() << " , "
	      << tmpTpcExitPoint.Y() << " , "
	      << tmpTpcExitPoint.Z() << " ) " << std::endl;
    tmpTpcExitPoint.SetX(-9999.);
    tmpTpcExitPoint.SetY(-9999.);
    tmpTpcExitPoint.SetZ(-9999.);
  }

  // Mike Lisa: OK, let's try something a little more along the lines
  // of NA49 and E895 strategy. Calculate the "nominal" position at N
  // radii (say N=11) within the TPC, and for a pair cut use the average
  // separation of these N
  // Grigory Nigmatkulov: For the future measurements N was changed
  // to mNumberOfPoints and the *magic numbers* were changed to the
  // mInnerTpcRadius and mOuterTpcRadius
  int irad = 0;
  candidates = hel.pathLength( mInnerTpcRadius );
  float step  = (mOuterTpcRadius - mInnerTpcRadius) / ( mNumberOfPoints - 1 ) ;
  sideLength = (candidates.first > 0) ? candidates.first : candidates.second;

  // Declare and initialize variable outside the loop
  float radius = mInnerTpcRadius;

  // Loop over radii
  while( irad<mNumberOfPoints && !std::isnan( sideLength ) ) {

    radius = mInnerTpcRadius + irad * step;
    candidates = hel.pathLength( radius );
    sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
    tmpPosSample[irad] = hel.at(sideLength);

    if(std::isnan(tmpPosSample[irad].x()) ||
       std::isnan(tmpPosSample[irad].y()) ||
       std::isnan(tmpPosSample[irad].z()) ) {

      std::cout << "tmpPosSample for radius = " << radius << " NAN"<< std::endl;
      std::cout << "tmpPosSample = ( "
		<< tmpPosSample[irad].X() << " , "
		<< tmpPosSample[irad].Y() << " , "
		<< tmpPosSample[irad].Z() << " ) " << std::endl;
      tmpPosSample[irad] =  TVector3(-9999.,-9999.,-9999);
    }

    // Do not forget to increment radii
    irad++;

    if ( irad<mNumberOfPoints ) {
      float radius = mInnerTpcRadius + irad*step;
      candidates = hel.pathLength(radius);
      sideLength = (candidates.first > 0) ? candidates.first : candidates.second;
    }
  } //while( irad<11 && !std::isnan(sideLength) )

  // In case, when track left TPC, the rest of postions will
  // be set to unphysical values
  for (int i=irad; i<11; i++) {
    tmpPosSample[i] = TVector3(-9999.,-9999.,-9999);
  }
  /*
  // Lets see if it will work out and move it to the protected area
  static float tRowRadius[mNumberOfPadrows] = {60,64.8,69.6,74.4,79.2,84,88.8,93.6,98.8,
					       104,109.2,114.4,119.6,127.195,129.195,131.195,
					       133.195,135.195,137.195,139.195,141.195,
					       143.195,10.45195,147.195,149.195,151.195,
					       153.195,155.195,157.195,159.195,161.195,
					       163.195,165.195,167.195,169.195,171.195,
					       173.195,175.195,177.195,179.195,181.195,
					       183.195,185.195,187.195,189.195};
  */

  int tRow,tSect,tOutOfBound;
  double tLength,tPhi;
  float tU;
  TVector3 tPoint(0,0,0);
  TVector3 tn(0,0,0);
  TVector3 tr(0,0,0);
  int ti = 0;

  // Test to enter the loop
  candidates =  hel.pathLength(tRowRadius[ti]);
  tLength = (candidates.first > 0) ? candidates.first : candidates.second;

  if ( std::isnan(tLength) ) {

    std::cout << "tLength Init tmp NAN" << std::endl;
    std::cout << "padrow number = " << ti << " not reached " << std::endl;
    std::cout << "*** DO NOT ENTER THE LOOP***" << std::endl;

    // Sector
    tmpSect[ti] = -1;
  } //if ( std::isnan(tLength) )

  // Start iteration over all padrows
  while( ti<mNumberOfPadrows && !std::isnan(tLength) ) {

    candidates =  hel.pathLength(tRowRadius[ti]);
    tLength = (candidates.first > 0) ? candidates.first : candidates.second;

    if ( std::isnan(tLength) ) {

      std::cout << "tLength loop 1st NAN" << std::endl;
      std::cout << "padrow number =  " << ti << " not reached" << std::endl;
      std::cout << "*** THIS IS AN ERROR SHOULDN'T  LOOP ***" << std::endl;
      // Sector
      tmpSect[ti] = -1;
    } // if ( std::isnan(tLength) )

    tPoint = hel.at(tLength);
    // Find which sector it is on
    TpcLocalTransform( tPoint, tmpSect[ti], tRow, tU, tPhi );
    // int tmpSectti = -1;
    // TpcLocalTransform( tPoint, tmpSectti, tRow, tU, tPhi );
    // tmpSect[ti] = tmpSectti;

    if ( std::isnan(tmpSect[ti]) ) {
      std::cout << "***ERROR tmpSect" << std::endl;
    }
    if ( std::isnan(tRow) ) {
      std::cout << "***ERROR tRow" << std::endl;
    }
    if ( std::isnan(tU) ) {
      std::cout << "***ERROR tU" << std::endl;
    }
    if ( std::isnan(tPhi) ) {
      std::cout << "***ERROR tPhi" << std::endl;
    }

    // calculate crossing plane
    tn.SetX( TMath::Cos(tPhi) );
    tn.SetY( TMath::Sin(tPhi) );
    tr.SetX( tRowRadius[ti] * TMath::Cos(tPhi) );
    tr.SetY( tRowRadius[ti] * TMath::Sin(tPhi) );

    // Find crossing point
    tLength = hel.pathLength(tr,tn);
    if ( std::isnan(tLength) ) {

      std::cout <<"tLength loop 2nd NAN" << std::endl;
      std::cout <<"padrow number = " << ti << " not reached" << std::endl;
      // Sector
      tmpSect[ti] = -2;
    } //if ( std::isnan(tLength) )

    tPoint = hel.at(tLength);
    tmpZ[ti] = tPoint.z();
    tOutOfBound = TpcLocalTransform( tPoint, tSect, tRow, tmpU[ti],  tPhi );
    if ( std::isnan(tSect) ) {
      std::cout << "***ERROR tSect 2" << std::endl;
    }
    if ( std::isnan(tRow) ) {
      std::cout << "***ERROR tRow 2" << std::endl;
    }
    if ( std::isnan(tmpU[ti]) ) {
      std::cout << "***ERROR tmpU[ti] 2" << std::endl;
    }
    if ( std::isnan(tPhi) ) {
      std::cout << "***ERROR tPhi 2 " << std::endl;
    }

    if( tOutOfBound ||
				( ( tmpSect[ti] == tSect ) && ( tRow != (ti+1) ) ) ) {
      // Sector
      tmpSect[ti]=-2;
    }
    else {
      if( tmpSect[ti] != tSect ) {
        // Try again on the other sector
        tn.SetX( TMath::Cos(tPhi) );
        tn.SetY( TMath::Sin(tPhi) );
        tr.SetX( tRowRadius[ti] * TMath::Cos(tPhi) );
        tr.SetY( tRowRadius[ti] * TMath::Sin(tPhi) );

        // Find crossing point
        tLength = hel.pathLength(tr,tn);
        tPoint = hel.at(tLength);
        if ( std::isnan(tLength) ) {
         std::cout << "tLength loop 3rd NAN" << std::endl;
         std::cout << "padrow number = " << ti << " not reached" << std::endl;
         // Sector
         tmpSect[ti] = -1;
        } //if ( std::isnan(tLength) )
        tmpZ[ti] = tPoint.z();
        tmpSect[ti] = tSect;
        tOutOfBound = TpcLocalTransform( tPoint, tSect, tRow, tmpU[ti], tPhi);

        if ( std::isnan(tSect) ) {
        std::cout << "***ERROR tSect 3" << std::endl;
        }
        if ( std::isnan(tRow) ) {
        std::cout << "***ERROR tRow 3" << std::endl;
        }
        if ( std::isnan(tmpU[ti]) ) {
        std::cout << "***ERROR tmpU[ti] 3" << std::endl;
        }
        if ( std::isnan(tPhi) ) {
        std::cout << "***ERROR tPhi 3 " << std::endl;
        }

        if( tOutOfBound || ( tSect != tmpSect[ti] ) || ( tRow!=(ti+1) ) ) {
        // Sector
        tmpSect[ti] = -1;
        } //if( tOutOfBound || ( tSect != tmpSect[ti] ) || ( tRow!=(ti+1) ) )
      } //if(tmpSect[ti] != tSect)
    } //else

    if (std::isnan(tmpSect[ti])){
      std::cout << "*******************ERROR***************************" << std::endl;
      std::cout <<"StHbtParticle--Fctn tmpSect=" << tmpSect[ti] << std::endl;
      std::cout << "*******************ERROR***************************" << std::endl;
    }
    if (std::isnan(tmpU[ti])){
      std::cout << "*******************ERROR***************************" << std::endl;
      std::cout <<"StHbtParticle--Fctn tmpU=" << tmpU[ti] << std::endl;
      std::cout << "*******************ERROR***************************" << std::endl;
    }
    if (std::isnan(tmpZ[ti])){
      std::cout << "*******************ERROR***************************" << std::endl;
      std::cout <<"StHbtParticle--Fctn tmpZ=" << tmpZ[ti] << std::endl;
      std::cout << "*******************ERROR***************************" << std::endl;
    }

    // If padrow ti not reached all other beyond are not reached
    // in this case set sector to -1
    if ( tmpSect[ti] == -1 ) {
      for (int tj=ti; tj<mNumberOfPadrows; tj++) {
				tmpSect[tj] = -1;
				ti=mNumberOfPadrows;
      }
    } //if ( tmpSect[ti] == -1 )

    // Increment padrow
    ti++;

    if ( ti<mNumberOfPadrows ) {
      candidates =  hel.pathLength(tRowRadius[ti]);
      tLength = (candidates.first > 0) ? candidates.first : candidates.second;
    } //if ( ti<mNumberOfPadrows )
  } //while( ti<mNumberOfPadrows && !std::isnan(tLength) )

}

                                    
const float mMaxDuInner = 0.8f;
const float mMaxDzInner = 3.0f;
const float mMaxDuOuter = 1.4f;
const float mMaxDzOuter = 3.2f;

//_________________
double fractionOfMergedRow(StPicoPhysicalHelix helixTrk1, StPicoPhysicalHelix helixTrk2) {

  // Calculate merging factor for the pair in STAR TPC
  double tDu, tDz;
  int tN = 0;
  double mFracOfMergedRow = 0.;
  double tDist;
  double tDistMax = 200.;
  
  // Primary and secondary vertex positions INTENTIOANLLY set to (0,0,0)
  // in order to make all estimations for future pair cuts,
  // i.e. in order to remove merged tracks.
  // This also implies, that all tracks originate from (0,0,0)

  TVector3 pVtx( 0., 0., 0. );
  TVector3 sVtx( 0., 0., 0. );
  TVector3 entrancePoint(0, 0, 0);
  TVector3 exitPoint(0, 0, 0);
  
  /// Information about hit position in TPC local coordinate system
  TVector3 posSampleTrk1[mNumberOfPoints];
  float mZTrk1[mNumberOfPadrows];
  float mUTrk1[mNumberOfPadrows];
  int mSectTrk1[mNumberOfPadrows];
  TVector3 posSampleTrk2[mNumberOfPoints];
  float mZTrk2[mNumberOfPadrows];
  float mUTrk2[mNumberOfPadrows];
  int mSectTrk2[mNumberOfPadrows];

  calculateTpcExitAndEntrancePoints(helixTrk1, pVtx, sVtx, entrancePoint, exitPoint,
                                    &posSampleTrk1[0], &mZTrk1[0], &mUTrk1[0], &mSectTrk1[0] );
  calculateTpcExitAndEntrancePoints(helixTrk2, pVtx, sVtx, entrancePoint, exitPoint,
                                    &posSampleTrk2[0], &mZTrk2[0], &mUTrk2[0], &mSectTrk2[0] );


  for ( int ti=0; ti < mNumberOfPadrows ; ti++ ) {

    if ( ( mSectTrk1[ti] == mSectTrk2[ti] ) &&
	       mSectTrk1[ti] != -1 ) {
      tDu = TMath::Abs( mUTrk1[ti] - mUTrk2[ti] );
      tDz = TMath::Abs( mZTrk1[ti] - mZTrk2[ti] );
      tN++;

      if ( ti<13 ) {
      	mFracOfMergedRow += ( tDu<mMaxDuInner && tDz<mMaxDzInner );
      	tDist = TMath::Sqrt( tDu * tDu / mMaxDuInner / mMaxDuInner +
      		tDz * tDz / mMaxDzInner / mMaxDzInner );
	      //mFracOfMergedRow += (tDu<mMaxDuInner && tDz<mMaxDzInner);
      }
      else {
      	mFracOfMergedRow += ( tDu<mMaxDuOuter && tDz<mMaxDzOuter );
      	tDist = TMath::Sqrt( tDu * tDu / mMaxDuOuter / mMaxDuOuter +
			    tDz * tDz / mMaxDzOuter / mMaxDzOuter );
	      //mFracOfMergedRow += (tDu<mMaxDuOuter && tDz<mMaxDzOuter);
      }

      if ( tDist<tDistMax ) {
        tDistMax = tDist;
      }
    }
  } // for ( int ti=0; ti < mTrack1->mNumberOfPadrows ; ti++ )

  if ( tN>0 ) {
    mFracOfMergedRow /= tN;
  }
  else {
    mFracOfMergedRow = -1.;
  }
  return mFracOfMergedRow;
}
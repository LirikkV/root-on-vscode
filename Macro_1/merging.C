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

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVector3.h"

// PicoDst headers
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"
#include "StPicoTrack.h"
#include "StPicoPhysicalHelix.h"
#include "StPicoBTofHit.h"
#include "StPicoBTowHit.h"
#include "StPicoEmcTrigger.h"
#include "StPicoBTofPidTraits.h"
#include "StPicoTrackCovMatrix.h"
#include "StPicoFmsHit.h"
#include "StPicoETofHit.h"
#include "StPicoEpdHit.h"

#include <algorithm>
//_________________
int TpcLocalTransform( TVector3& aPoint, int& aSector, int& aRow, float& aU, double& aPhi);
//_________________
void calculateTpcExitAndEntrancePoints(StPicoPhysicalHelix tHelix,TVector3 PrimVert,TVector3 SecVert,TVector3 tmpTpcEntrancePoint,TVector3 tmpTpcExitPoint,TVector3* tmpPosSample,float* tmpZ,float* tmpU,int* tmpSect);
//_________________
double fractionOfMergedRow(StPicoPhysicalHelix helixTrk1, StPicoPhysicalHelix helixTrk2);

const float fPionMassSq = 0.019480064;

struct MyParticle {
	TVector3 P_tot;
	UInt_t topologyMap0;
	UInt_t topologyMap1;
	ULong64_t iTpcTopologyMap;
	Int_t nHits;
  StPicoPhysicalHelix helix;
};


// ./picoAnalyzerStandalone /home/ubuntu/Data/st_physics_21029052_raw_6500014.picoDst.root oTest.root

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
  TH2F* hFMRvsQinv = new TH2F("hFMRvsQinv",";q_{inv} [GeV/c];Fraction of Merged Rows",100,0.,1.,220,-1.1,1.1);
  StPicoDstReader* picoReader = new StPicoDstReader(fileName);
  picoReader->Init();

  // This is a way if you want to spead up I/O
  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event*", 1);
  picoReader->SetStatus("Track*", 1);
  picoReader->SetStatus("BTofHit*", 1);
  picoReader->SetStatus("BTofPidTraits*", 1);
  // picoReader->SetStatus("BTowHit*", 1);
  // picoReader->SetStatus("ETofHit*", 1);
  // picoReader->SetStatus("EpdHit*", 1);
  // picoReader->SetStatus("EmcTrigger",0);
  // picoReader->SetStatus("TrackCovMatrix",1);
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


    TVector3 pVtx = event->primaryVertex();

    // Cut in Vertex Position
    float Vx = pVtx.X();
    float Vy = pVtx.Y();
    float Vz = pVtx.Z();
    // collider mode
    // if( TMath::Abs(Vz) > 145 ) continue; 
    // if( TMath::Sqrt( Vx*Vx + Vy*Vy ) > 2. ) continue;
    // FXT mode
    if( TMath::Abs(Vz-200) > 2. ) continue; 
    if( TMath::Sqrt( Vx*Vx + (Vy+2.)*(Vy+2) ) > 2. ) continue;

	


    // Track analysis
    Int_t nTracks = dst->numberOfTracks();
    Int_t nMatrices = dst->numberOfTrackCovMatrices();
    if(nTracks != nMatrices) {
      std::cout << "Number of tracks and matrices do not match!" << std::endl;
    }
    std::vector< MyParticle > vPions[2];
    // Track loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
      // Retrieve i-th pico track
      StPicoTrack *picoTrack = dst->track(iTrk);
      if(!picoTrack) continue;

      // cut for the primary track
      if( !picoTrack->isPrimary() ) {
              continue;
      }
      if( picoTrack->nHitsFit() < 15 ) continue;
      // cut on total momentum
      if( 0.15 > picoTrack->pMom().Mag() ||  picoTrack->pMom().Mag() > 1.5 ) continue;
      Double_t DCA = picoTrack->gDCA(pVtx).Mag();
      if( TMath::Abs(DCA) > 3) continue;
      Double_t Pt = picoTrack->pMom().Pt();
      Double_t eta = picoTrack->pMom().Eta();
      if( 0.15 >  Pt || Pt >  1.5 ) continue;
	    if(TMath::Abs(eta) > 1.5) continue;
      Int_t iCh = (picoTrack->charge() > 0) ? 0 : 1;
      MyParticle myPart;
      myPart.P_tot = picoTrack->pMom();
      myPart.topologyMap0 = picoTrack->topologyMap(0);
      myPart.topologyMap1 = picoTrack->topologyMap(1);
      myPart.iTpcTopologyMap = picoTrack->iTpcTopologyMap();
      myPart.nHits = picoTrack->nHits();
      myPart.helix = picoTrack->helix( event->bField() );



      // Cut for TPC-only identification
	    if( 0.15 < picoTrack->pMom().Mag() &&  picoTrack->pMom().Mag() < 0.55 ) {
      	if( TMath::Abs(picoTrack->nSigmaPion()) > 2 ) continue;
        if( TMath::Abs(picoTrack->nSigmaElectron()) < 2  ) continue;
        if( TMath::Abs(picoTrack->nSigmaKaon()) < 2 ) continue;   
        if( TMath::Abs(picoTrack->nSigmaProton()) < 2 ) continue;
        vPions[iCh].push_back( myPart );
	    } 
      // Check if track has TOF signal
      if( picoTrack->isTofTrack() ) {
	      // Retrieve corresponding trait
	      StPicoBTofPidTraits *trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
	      Double_t beta =  trait->btofBeta();
	      if( 0.55 < picoTrack->pMom().Mag() &&  picoTrack->pMom().Mag() < 1.5 ) {
          if( TMath::Abs(picoTrack->nSigmaPion()) > 3 ) continue;
          Double_t betaExptdPion = TMath::Sqrt(picoTrack->pMom().Mag()*picoTrack->pMom().Mag()/(picoTrack->pMom().Mag()*picoTrack->pMom().Mag()+0.14*0.14));
          Double_t m2 = Pt*Pt*((1/(beta*beta)) - 1);   // m^2=(1/beta^2-1)*p^2
          if(TMath::Abs(1/beta - 1/betaExptdPion) > 0.015) continue;
          if( m2 < -0.05 || m2 > 0.08 ) continue;
          vPions[iCh].push_back( myPart );
	      }
      } // if( picoTrack->isTofTrack() )
    } // for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    /////////////////////////////////////////
    ////// Fraction of Merged Rows //////////
    /////////////////////////////////////////
    for (int iCh=0; iCh<2; iCh++) {
      for( int i = 0; i < int(vPions[iCh].size()); i++  ) {
        for( int j = i+1; j< int(vPions[iCh].size()); j++) {
          // fPionMass
          double dQinv = TMath::Sqrt((vPions[iCh][i].P_tot - vPions[iCh][j].P_tot).Mag2() -
                                     TMath::Power(
                                      TMath::Sqrt(vPions[iCh][i].P_tot.Mag2() + fPionMassSq) - 
                                      TMath::Sqrt(vPions[iCh][j].P_tot.Mag2() + fPionMassSq), 2.) );
          double dFractionOfMergedRow = fractionOfMergedRow(vPions[iCh][i].helix, vPions[iCh][j].helix);
          hFMRvsQinv->Fill(dQinv,dFractionOfMergedRow);
		    } // for( int j = i+1; j< int(vPions[iCh].size()); j++)
	    } // for( int i = 0; i < int(vPions[iCh].size()); i++  )
    } // for (int iCh=0; iCh<2; iCh++)
    for(int i = 0; i < 2; i++) vPions[i].clear();

  } // for event loop (Long64_t iEvent=0; iEvent<events2read; iEvent++)



  TFile* fo = new TFile(oFileName,"recreate");
  hFMRvsQinv->Write();
  fo->Close();
  picoReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
  << std::endl;
}
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


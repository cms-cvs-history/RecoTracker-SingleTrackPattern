#ifndef RecoTracker_SingleTrackPattern_AnalyzeHitEff_h
#define RecoTracker_SingleTrackPattern_AnalyzeHitEff_h

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoTracker/SingleTrackPattern/interface/CosmicTrajectoryBuilder.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementVector.h"
#include "RecoTracker/SingleTrackPattern/test/TrackLocalAngle.h"

#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <vector>
#include <TTree.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <Riostream.h>
#include "TRandom2.h"

class AnalyzeHitEff : public edm::EDAnalyzer
{
 public:
  
  explicit AnalyzeHitEff(const edm::ParameterSet& conf);
  
  virtual ~AnalyzeHitEff();
  virtual void beginJob(const edm::EventSetup& c);
  virtual void endJob(); 
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);

 private:

  
  edm::ParameterSet conf_;
  TFile* hFile;
  TTree* n;
  int events,EventSeedCKF,EventTrackCKF,EventTriggCKF;


  ofstream out1,out2;
  TrackLocalAngle* theAngleFinder;
  uint layers;
  vector<uint> ActiveLayStereo, ActiveLayMono;
  std::vector<const SiStripRecHit2D*> vRPhi_SiStripRecHit2D;
  std::vector<const SiStripRecHit2D*> vSte_SiStripRecHit2D;
  TRandom2 RanGen2;

  TH1F *AngleLayTIB;
  TH1F *CKFLayEff, *InvCKFLayEff;
  TH1F *TIFInvModEffStereo, *TIFInvModEffStereo_all, *TIFModEffStereo;
  TH1F *TIFInvModEffRPhi, *TIFInvModEffRPhi_all, *TIFModEffRPhi;
  TH1F *ResidualXValidSte, *ResidualXValidSte_2;
  TH1F *ResidualXValidRPhi, *ResidualXValidRPhi_2;

  // Tree declarations
  // All RecHit
  int RHNum;
  float	RHgpx[200];
  float	RHgpy[200];
  float	RHgpz[200];
  int RHMod[200];
  int RphiSte[200];

  // Track RecHit
  int TRHNum;
  float	TRHgpx[200];
  float	TRHgpy[200];
  float	TRHgpz[200];
  int TRHLay[200];
  uint TRHMod[200];

};


#endif

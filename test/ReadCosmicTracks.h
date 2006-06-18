#ifndef RecoTracker_SingleTrackPattern_ReadCosmicTracks_h
#define RecoTracker_SingleTrackPattern_ReadCosmicTracks_h

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Vector/interface/GlobalPoint.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>

class CompareTOF {
 public:
  CompareTOF(){};
   bool operator()( PSimHit ps1,
		    PSimHit ps2){

     return ps1.tof()<ps2.tof();};
 };
class ReadCosmicTracks : public edm::EDAnalyzer
{
 typedef TrajectoryStateOnSurface     TSOS;
 public:
  
  explicit ReadCosmicTracks(const edm::ParameterSet& conf);
  
  virtual ~ReadCosmicTracks();
  virtual void beginJob(const edm::EventSetup& c);
  virtual void endJob(); 
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
 
  void makeResiduals(const TrajectorySeed& seed,
		     const TrackingRecHitCollection & hits,
		     const edm::Event& e, 
		     const edm::EventSetup& es);
  TSOS startingTSOS(const TrajectorySeed& seed)const;
 private:
  edm::ParameterSet conf_;
  std::vector<PSimHit> theStripHits;
  edm::ESHandle<TrackerGeometry> tracker;
  edm::ESHandle<MagneticField> magfield;
  TFile* hFile;
  TH1F  *hptres;
  TH1F *hptres1,*hptres2,*hptres3,*hptres4,*hptres5,*hptres6,*hptres7,*hptres8,*hptres9,*hptres10;
 TH1F  *hchiSq;
  TH1F *hchiSq1,*hchiSq2,*hchiSq3,*hchiSq4,*hchiSq5,*hchiSq6,*hchiSq7,*hchiSq8,*hchiSq9,*hchiSq10;
  TH1F  *heffpt,*hrespt,*hchipt,*hchiR;
  TH1F  *heffhit,*hcharge;
  TH1F  *hresTIB,*hresTOB,*hresTID,*hresTEC;
  uint inum[10];
  uint iden[10];
  uint inum2[9];
  uint iden2[9];
  float ichiR[10];
  uint iden3[10];
  bool seed_plus;
  PropagatorWithMaterial  *thePropagator;
  PropagatorWithMaterial  *thePropagatorOp;
  KFUpdator *theUpdator;
  Chi2MeasurementEstimator *theEstimator;
  const TransientTrackingRecHitBuilder *RHBuilder;
  const KFTrajectorySmoother * theSmoother;
  const KFTrajectoryFitter * theFitter;
  TrajectoryStateTransform tsTransform;
};


#endif

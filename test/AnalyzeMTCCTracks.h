#ifndef RecoTracker_SingleTrackPattern_AnalyzeMTCCTracks_h
#define RecoTracker_SingleTrackPattern_AnalyzeMTCCTracks_h

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Vector/interface/GlobalPoint.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoTracker/SingleTrackPattern/interface/CosmicTrajectoryBuilder.h"
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>

class AnalyzeMTCCTracks : public edm::EDAnalyzer
{
typedef TrajectoryStateOnSurface     TSOS;
 public:
  
  explicit AnalyzeMTCCTracks(const edm::ParameterSet& conf);
  
  virtual ~AnalyzeMTCCTracks();
  virtual void beginJob(const edm::EventSetup& c);
  virtual void endJob(); 
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  void makeResiduals(const TrajectorySeed& seed,
		     const TrackingRecHitCollection &hits,
		     const edm::Event& e, 
		     const edm::EventSetup& es);

  std::vector<TrajectoryMeasurement> seedMeasurements(const TrajectorySeed& seed) const;
  TrajectoryStateOnSurface startingTSOS(const TrajectorySeed& seed)const;
  Trajectory createStartingTrajectory( const TrajectorySeed& seed) const;
 private:
  edm::ParameterSet conf_;

  TFile* hFile;
  TH1F  *hphi, *hnhit,*hchi,*hresTOB,*hresTIB,*heta;
  bool seed_plus;
  PropagatorWithMaterial  *thePropagator;
  PropagatorWithMaterial  *thePropagatorOp;
  KFUpdator *theUpdator;
  Chi2MeasurementEstimator *theEstimator;
  const TransientTrackingRecHitBuilder *RHBuilder;
  const KFTrajectorySmoother * theSmoother;
  const KFTrajectoryFitter * theFitter;
  TrajectoryStateTransform tsTransform;
  edm::ESHandle<TrackerGeometry> tracker;
  edm::ESHandle<MagneticField> magfield;
  const  CosmicTrajectoryBuilder *pippo;
};


#endif

#ifndef RecoTracker_SingleTrackPattern_AnalyzeMTCCTracks_h
#define RecoTracker_SingleTrackPattern_AnalyzeMTCCTracks_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Vector/interface/GlobalPoint.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>

class AnalyzeMTCCTracks : public edm::EDAnalyzer
{
 public:
  
  explicit AnalyzeMTCCTracks(const edm::ParameterSet& conf);
  
  virtual ~AnalyzeMTCCTracks();
  virtual void beginJob(const edm::EventSetup& c);
  virtual void endJob(); 
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  
 private:
  edm::ParameterSet conf_;
  std::vector<PSimHit> theStripHits;
  TFile* hFile;
  TH1F  *hphi, *hnhit;
};


#endif

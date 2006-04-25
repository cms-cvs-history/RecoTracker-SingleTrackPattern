#ifndef RecoTracker_SingleTrackPattern_ReadCosmicTracks_h
#define RecoTracker_SingleTrackPattern_ReadCosmicTracks_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>


class ReadCosmicTracks : public edm::EDAnalyzer
{
 public:
  
  explicit ReadCosmicTracks(const edm::ParameterSet& conf);
  
  virtual ~ReadCosmicTracks();
  
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  
 private:
  edm::ParameterSet conf_;
  std::vector<PSimHit> theStripHits;
};


#endif

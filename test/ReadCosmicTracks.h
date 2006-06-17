#ifndef RecoTracker_SingleTrackPattern_ReadCosmicTracks_h
#define RecoTracker_SingleTrackPattern_ReadCosmicTracks_h

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

class CompareTOF {
 public:
  CompareTOF(){};
   bool operator()( PSimHit ps1,
		    PSimHit ps2){

     return ps1.tof()<ps2.tof();};
 };
class ReadCosmicTracks : public edm::EDAnalyzer
{
 public:
  
  explicit ReadCosmicTracks(const edm::ParameterSet& conf);
  
  virtual ~ReadCosmicTracks();
  virtual void beginJob(const edm::EventSetup& c);
  virtual void endJob(); 
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  
 private:
  edm::ParameterSet conf_;
  std::vector<PSimHit> theStripHits;
  TFile* hFile;
  TH1F  *hptres;
  TH1F *hptres1,*hptres2,*hptres3,*hptres4,*hptres5,*hptres6,*hptres7,*hptres8,*hptres9,*hptres10;
 TH1F  *hchiSq;
  TH1F *hchiSq1,*hchiSq2,*hchiSq3,*hchiSq4,*hchiSq5,*hchiSq6,*hchiSq7,*hchiSq8,*hchiSq9,*hchiSq10;
  TH1F  *heffpt,*hrespt,*hchipt,*hchiR;
  TH1F  *heffhit;
  uint inum[10];
  uint iden[10];
  uint inum2[9];
  uint iden2[9];
  float ichiR[10];
  uint iden3[10];
};


#endif

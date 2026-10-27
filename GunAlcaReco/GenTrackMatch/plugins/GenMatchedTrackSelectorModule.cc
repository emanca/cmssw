// GenMatchedTrackSelectorModule
//
// Stream selector module for gen-matched track selection in particle-gun
// ALCARECO production. Modeled exactly on AlignmentTrackSelectorModule:
// the selection logic lives in the GenMatchedTrackSelector functor
// (GunAlcaReco/GenTrackMatch/interface), and ObjectSelectorStream supplies
// the complete stock persistence family for the selection:
//
//   vector<reco::Track>                  (matched clones)
//   vector<reco::TrackExtra>             (cloned extras, refs re-pointed)
//   TrackingRecHits owned               (cloned hits of the clones)
//   edmNew::DetSetVector<SiPixelCluster> (cloned clusters used by the hits)
//   edmNew::DetSetVector<SiStripCluster> (cloned clusters used by the hits)
//
// With filter = True the event passes iff at least one track is gen-matched,
// so pileup-only events are dropped. This is the stock-ALCARECO-shaped
// replacement for the GenMatchedTrackFilter EDFilter (whose clones' extras
// dangled, see CAMPAIGN_NOTES v2.0 entry).

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelectorStream.h"

// Crucial: pulls the TrackCollection StoreManager specialization that clones the
// full sibling family (extras + owns hits + clusters). Without it the module
// compiles fine but emits reco::Tracks only (comment verbatim from the TkAl
// selector module pattern this follows).
#include "CommonTools/RecoAlgos/interface/TrackSelector.h"

#include "GunAlcaReco/GenTrackMatch/interface/GenMatchedTrackSelector.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

struct GenMatchedTrackSelectorConfig {
  typedef GenMatchedTrackSelector::container container;
  typedef container::const_iterator const_iterator;
  typedef reco::TrackCollection collection;

  GenMatchedTrackSelectorConfig(const edm::ParameterSet& cfg, edm::ConsumesCollector&& iC)
      : theSelector(cfg, iC) {}

  const_iterator begin() const { return theSelectedTracks.begin(); }
  const_iterator end() const { return theSelectedTracks.end(); }
  size_t size() const { return theSelectedTracks.size(); }

  void select(const edm::Handle<reco::TrackCollection>& c,
              const edm::Event& evt,
              const edm::EventSetup& eSetup) {
    theSelectedTracks.clear();
    for (reco::TrackCollection::const_iterator i = c.product()->begin(); i != c.product()->end(); ++i) {
      theSelectedTracks.push_back(&*i);
    }
    theSelectedTracks = theSelector.select(theSelectedTracks, evt, eSetup);
  }

private:
  container theSelectedTracks;

  GenMatchedTrackSelector theSelector;
};

typedef ObjectSelectorStream<GenMatchedTrackSelectorConfig> GenMatchedTrackSelectorModule;

DEFINE_FWK_MODULE(GenMatchedTrackSelectorModule);

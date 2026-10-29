/// \class AlCaRecoMuonRelinker
///
/// Writes a reco::Muon collection for an ALCARECO stream whose track
/// references all resolve inside the ALCARECO file:
///
///   - innerTrack is re-pointed to the matching entry of the stream's
///     selected (cloned) track collection,
///   - the standalone, global and TeV refit (TPFMS, Picky, DYT) tracks are
///     cloned into small track-only collections (TrackExtraRef nulled; the
///     track parameters, covariance, hit pattern and chi2 are kept) and the
///     muon references are re-pointed to the clones.
///
/// A selected track is matched to its source track (e.g. generalTracks) by
/// exact equality of charge, parameters and number of valid hits: the
/// ALCARECO track selectors copy the source tracks without modification, so
/// this works for any selector and needs no index map. A muon is linked if
/// its innerTrack is a source track that was selected.
///
/// Also writes an edm::Association<reco::MuonCollection> keyed on the
/// selected tracks (instance "trackToMuon") pointing into the output muons.
///
/// With onlyLinked = True, muons without a selected inner track are dropped;
/// otherwise they are kept with a null innerTrack (and a muonBestTrack that
/// is null if the best track type is InnerTrack).

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/OrphanHandle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <array>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
  // muon track types cloned into their own collections, with the instance labels
  const std::array<std::pair<reco::Muon::MuonTrackType, const char*>, 5> kClonedTypes{
      {{reco::Muon::OuterTrack, "standAlone"},
       {reco::Muon::CombinedTrack, "global"},
       {reco::Muon::TPFMS, "tpfms"},
       {reco::Muon::Picky, "picky"},
       {reco::Muon::DYT, "dyt"}}};

  bool sameTrack(const reco::Track& a, const reco::Track& b) {
    return a.charge() == b.charge() && a.parameters() == b.parameters() &&
           a.numberOfValidHits() == b.numberOfValidHits();
  }
}  // namespace

class AlCaRecoMuonRelinker : public edm::stream::EDProducer<> {
public:
  explicit AlCaRecoMuonRelinker(const edm::ParameterSet& cfg)
      : muonsToken_(consumes<reco::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
        selectedTracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("selectedTracks"))),
        sourceTracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("sourceTracks"))),
        onlyLinked_(cfg.getParameter<bool>("onlyLinked")) {
    produces<reco::MuonCollection>();
    for (const auto& ct : kClonedTypes) {
      produces<reco::TrackCollection>(ct.second);
    }
    produces<edm::Association<reco::MuonCollection>>("trackToMuon");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
    desc.add<edm::InputTag>("selectedTracks", edm::InputTag(""));
    desc.add<edm::InputTag>("sourceTracks", edm::InputTag("generalTracks"));
    desc.add<bool>("onlyLinked", true);
    descriptions.add("alcaRecoMuonRelinker", desc);
  }

  void produce(edm::Event& evt, const edm::EventSetup&) override {
    edm::Handle<reco::MuonCollection> muonsH;
    evt.getByToken(muonsToken_, muonsH);
    edm::Handle<reco::TrackCollection> selH;
    evt.getByToken(selectedTracksToken_, selH);
    edm::Handle<reco::TrackCollection> srcH;
    evt.getByToken(sourceTracksToken_, srcH);

    // source track index -> selected track index
    std::unordered_map<unsigned int, unsigned int> srcToSel;
    for (unsigned int j = 0; j < selH->size(); ++j) {
      for (unsigned int k = 0; k < srcH->size(); ++k) {
        if (sameTrack((*selH)[j], (*srcH)[k])) {
          srcToSel.emplace(k, j);
          break;
        }
      }
    }

    std::map<std::string, std::unique_ptr<reco::TrackCollection>> clones;
    std::map<std::string, reco::TrackRefProd> cloneRefProds;
    std::map<std::string, std::map<std::pair<edm::ProductID, unsigned int>, unsigned int>> cloneIndex;
    for (const auto& ct : kClonedTypes) {
      clones[ct.second] = std::make_unique<reco::TrackCollection>();
      cloneRefProds[ct.second] = evt.getRefBeforePut<reco::TrackCollection>(ct.second);
    }

    auto outMuons = std::make_unique<reco::MuonCollection>();
    std::vector<int> selToMuon(selH->size(), -1);

    for (const auto& mu : *muonsH) {
      const reco::TrackRef& inner = mu.innerTrack();
      int selIdx = -1;
      if (inner.isNonnull() && inner.id() == srcH.id()) {
        const auto it = srcToSel.find(inner.key());
        if (it != srcToSel.end())
          selIdx = it->second;
      }
      if (selIdx < 0 && onlyLinked_)
        continue;

      reco::Muon outMu(mu);
      outMu.setInnerTrack(selIdx >= 0 ? reco::TrackRef(selH, selIdx) : reco::TrackRef());

      for (const auto& ct : kClonedTypes) {
        const reco::TrackRef ref = mu.muonTrack(ct.first);
        reco::TrackRef newRef;
        if (ref.isNonnull() && ref.isAvailable()) {
          auto& index = cloneIndex[ct.second];
          const auto key = std::make_pair(ref.id(), static_cast<unsigned int>(ref.key()));
          auto it = index.find(key);
          if (it == index.end()) {
            reco::Track clone(*ref);
            clone.setExtra(reco::TrackExtraRef());
            clones[ct.second]->push_back(clone);
            it = index.emplace(key, clones[ct.second]->size() - 1).first;
          }
          newRef = reco::TrackRef(cloneRefProds[ct.second], it->second);
        }
        outMu.setMuonTrack(ct.first, newRef);
      }

      if (selIdx >= 0)
        selToMuon[selIdx] = outMuons->size();
      outMuons->push_back(outMu);
    }

    for (const auto& ct : kClonedTypes) {
      evt.put(std::move(clones[ct.second]), ct.second);
    }
    const edm::OrphanHandle<reco::MuonCollection> muonsOut = evt.put(std::move(outMuons));

    auto assoc = std::make_unique<edm::Association<reco::MuonCollection>>(muonsOut);
    edm::Association<reco::MuonCollection>::Filler filler(*assoc);
    filler.insert(selH, selToMuon.begin(), selToMuon.end());
    filler.fill();
    evt.put(std::move(assoc), "trackToMuon");
  }

private:
  const edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  const edm::EDGetTokenT<reco::TrackCollection> selectedTracksToken_;
  const edm::EDGetTokenT<reco::TrackCollection> sourceTracksToken_;
  const bool onlyLinked_;
};

DEFINE_FWK_MODULE(AlCaRecoMuonRelinker);

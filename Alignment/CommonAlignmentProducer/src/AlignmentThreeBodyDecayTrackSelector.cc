#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <vector>

#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"

#include "Alignment/CommonAlignmentProducer/interface/AlignmentThreeBodyDecayTrackSelector.h"

using namespace edm;
using namespace std;

namespace {
  struct ThreeBodyCandidate {
    double pt;
    const reco::Track* kaon;
    const reco::Track* pion;
    const reco::Track* softPion;
  };
}

AlignmentThreeBodyDecayTrackSelector::AlignmentThreeBodyDecayTrackSelector(const edm::ParameterSet& cfg,
                                                                           edm::ConsumesCollector& iC) {
  (void)iC;
  LogDebug("Alignment") << "> applying three body decay Trackfilter ...";

  theMassrangeSwitch = cfg.getParameter<bool>("applyMassrangeFilter");
  theIntermediateMassSwitch = cfg.getParameter<bool>("applyIntermediateMassrangeFilter");
  theMassDifferenceSwitch = cfg.getParameter<bool>("applyMassDifferenceFilter");
  theChargeSwitch = cfg.getParameter<bool>("applyChargeFilter");

  theMinMass = cfg.getParameter<double>("minXMass");
  theMaxMass = cfg.getParameter<double>("maxXMass");
  theMinIntermediateMass = cfg.getParameter<double>("minIntermediateMass");
  theMaxIntermediateMass = cfg.getParameter<double>("maxIntermediateMass");
  theMinMassDifference = cfg.getParameter<double>("minMassDifference");
  theMaxMassDifference = cfg.getParameter<double>("maxMassDifference");
  theFirstDaughterMass = cfg.getParameter<double>("firstDaughterMass");
  theSecondDaughterMass = cfg.getParameter<double>("secondDaughterMass");
  theThirdDaughterMass = cfg.getParameter<double>("thirdDaughterMass");
  theCandNumber = cfg.getParameter<unsigned int>("numberOfCandidates");

  theCharge = cfg.getParameter<int>("charge");
  theUnsignedSwitch = cfg.getParameter<bool>("useUnsignedCharge");
}

AlignmentThreeBodyDecayTrackSelector::~AlignmentThreeBodyDecayTrackSelector() {}

bool AlignmentThreeBodyDecayTrackSelector::useThisFilter() {
  return theMassrangeSwitch || theIntermediateMassSwitch || theMassDifferenceSwitch || theChargeSwitch;
}

AlignmentThreeBodyDecayTrackSelector::Tracks AlignmentThreeBodyDecayTrackSelector::select(const Tracks& tracks,
                                                                                          const edm::Event& iEvent,
                                                                                          const edm::EventSetup& iSetup) {
  (void)iEvent;
  (void)iSetup;
  Tracks result = tracks;

  if (useThisFilter())
    result = checkMass(result);

  LogDebug("Alignment") << ">  ThreeBodyDecay tracks all,kept: " << tracks.size() << "," << result.size();
  return result;
}

AlignmentThreeBodyDecayTrackSelector::Tracks AlignmentThreeBodyDecayTrackSelector::checkMass(const Tracks& cands) const {
  Tracks result;

  if (cands.size() < 3)
    return result;

  vector<ThreeBodyCandidate> candCollection;

  for (unsigned int iK = 0; iK < cands.size(); ++iK) {
    const reco::Track* kaon = cands[iK];
    const int kaonCharge = kaon->charge();

    TLorentzVector kaonP4;
    kaonP4.SetXYZT(kaon->px(),
                   kaon->py(),
                   kaon->pz(),
                   std::sqrt(kaon->p() * kaon->p() + theFirstDaughterMass * theFirstDaughterMass));

    for (unsigned int iPi = 0; iPi < cands.size(); ++iPi) {
      if (iPi == iK)
        continue;

      const reco::Track* pion = cands[iPi];
      if (pion->charge() != -kaonCharge)
        continue;

      TLorentzVector pionP4;
      pionP4.SetXYZT(pion->px(),
                     pion->py(),
                     pion->pz(),
                     std::sqrt(pion->p() * pion->p() + theSecondDaughterMass * theSecondDaughterMass));

      const TLorentzVector intermediate = kaonP4 + pionP4;
      if (theIntermediateMassSwitch &&
          !(intermediate.M() > theMinIntermediateMass && intermediate.M() < theMaxIntermediateMass))
        continue;

      bool acceptedCandidate = false;
      ThreeBodyCandidate bestCandidate{0., nullptr, nullptr, nullptr};

      for (unsigned int iPis = 0; iPis < cands.size(); ++iPis) {
        if (iPis == iK || iPis == iPi)
          continue;

        const reco::Track* softPion = cands[iPis];
        if (softPion->charge() != pion->charge())
          continue;

        if (theChargeSwitch && !checkCharge(kaon, pion, softPion))
          continue;

        TLorentzVector softPionP4;
        softPionP4.SetXYZT(softPion->px(),
                           softPion->py(),
                           softPion->pz(),
                           std::sqrt(softPion->p() * softPion->p() + theThirdDaughterMass * theThirdDaughterMass));

        const TLorentzVector mother = intermediate + softPionP4;
        const double massDifference = mother.M() - intermediate.M();

        if (theMassrangeSwitch && !(mother.M() > theMinMass && mother.M() < theMaxMass))
          continue;

        if (theMassDifferenceSwitch &&
            !(massDifference > theMinMassDifference && massDifference < theMaxMassDifference))
          continue;

        if (!acceptedCandidate || mother.Pt() > bestCandidate.pt) {
          bestCandidate = {mother.Pt(), kaon, pion, softPion};
          acceptedCandidate = true;
        }
      }

      if (acceptedCandidate)
        candCollection.push_back(bestCandidate);
    }
  }

  if (candCollection.empty())
    return result;

  sort(candCollection.begin(), candCollection.end(), [](const ThreeBodyCandidate& a, const ThreeBodyCandidate& b) {
    return a.pt > b.pt;
  });

  map<const reco::Track*, unsigned int> uniqueTrackIndex;
  for (unsigned int i = 0; i < candCollection.size() && i < theCandNumber; ++i) {
    const ThreeBodyCandidate& candidate = candCollection[i];
    if (uniqueTrackIndex.find(candidate.kaon) == uniqueTrackIndex.end()) {
      result.push_back(candidate.kaon);
      uniqueTrackIndex[candidate.kaon] = i;
    }
    if (uniqueTrackIndex.find(candidate.pion) == uniqueTrackIndex.end()) {
      result.push_back(candidate.pion);
      uniqueTrackIndex[candidate.pion] = i;
    }
    if (uniqueTrackIndex.find(candidate.softPion) == uniqueTrackIndex.end()) {
      result.push_back(candidate.softPion);
      uniqueTrackIndex[candidate.softPion] = i;
    }
  }

  return result;
}

bool AlignmentThreeBodyDecayTrackSelector::checkCharge(const reco::Track* kaon,
                                                       const reco::Track* pion,
                                                       const reco::Track* softPion) const {
  const bool validD0Charges = (kaon->charge() + pion->charge() == 0);
  const bool validSoftCharge = (pion->charge() == softPion->charge());
  if (!validD0Charges || !validSoftCharge)
    return false;

  int totalCharge = kaon->charge() + pion->charge() + softPion->charge();
  if (theUnsignedSwitch)
    totalCharge = std::abs(totalCharge);

  return totalCharge == theCharge;
}

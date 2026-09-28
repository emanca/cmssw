#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace {
constexpr double kKaonMass = 0.493677;
constexpr double kPionMass = 0.139570;
constexpr double kD0Mass = 1.86483;
constexpr double kDstMass = 2.01026;
constexpr double kNominalDeltaMass = kDstMass - kD0Mass;
constexpr double kMassError = 1.e-6;
constexpr float kFloatSentinel = -99.f;
constexpr int kIntSentinel = -1;

struct GenChain {
  const reco::Candidate* dst = nullptr;
  const reco::Candidate* d0 = nullptr;
  const reco::Candidate* kaon = nullptr;
  const reco::Candidate* pion = nullptr;
  const reco::Candidate* softPion = nullptr;
  int dstMotherPdgId = 0;
  bool prompt = false;
  bool fromB = false;
  bool fromW = false;
  bool fromCharm = false;
};

struct TruthParticleRecord {
  int pdgId = 0;
  int status = 0;
  int charge = 0;
  float pt = kFloatSentinel;
  float eta = kFloatSentinel;
  float phi = kFloatSentinel;
  float mass = kFloatSentinel;
  float vx = kFloatSentinel;
  float vy = kFloatSentinel;
  float vz = kFloatSentinel;
};

struct RecoTrackRecord {
  int index = kIntSentinel;
  int charge = 0;
  float deltaR = kFloatSentinel;
  float pt = kFloatSentinel;
  float eta = kFloatSentinel;
  float phi = kFloatSentinel;
  float mass = kFloatSentinel;
  float px = kFloatSentinel;
  float py = kFloatSentinel;
  float pz = kFloatSentinel;
  float vx = kFloatSentinel;
  float vy = kFloatSentinel;
  float vz = kFloatSentinel;
  float dxyBeamSpot = kFloatSentinel;
  float dzBeamSpot = kFloatSentinel;
  float dxyPrimaryVertex = kFloatSentinel;
  float dzPrimaryVertex = kFloatSentinel;
  float normalizedChi2 = kFloatSentinel;
  int nValidHits = kIntSentinel;
};

struct FitRecord {
  bool valid = false;
  float chi2 = kFloatSentinel;
  float dof = kFloatSentinel;
  float pval = kFloatSentinel;
  float mass = kFloatSentinel;
};

struct MatchResult {
  bool valid = false;
  std::size_t chainIndex = 0;
  std::array<unsigned int, 3> trackIndices = {{0, 0, 0}};
  std::array<double, 3> deltaRs = {{0., 0., 0.}};
  double score = std::numeric_limits<double>::max();
};

struct TrackMatchOption {
  unsigned int index = 0;
  double deltaR = 0.;
};

int candidateCharge(const reco::Candidate& candidate) {
  return static_cast<int>(std::lround(candidate.charge()));
}

double deltaPhi(double first, double second) {
  double result = first - second;
  while (result > M_PI)
    result -= 2. * M_PI;
  while (result <= -M_PI)
    result += 2. * M_PI;
  return result;
}

double deltaR(const reco::Candidate& candidate, const reco::Track& track) {
  const double dEta = candidate.eta() - track.eta();
  const double dPhi = deltaPhi(candidate.phi(), track.phi());
  return std::sqrt(dEta * dEta + dPhi * dPhi);
}

TLorentzVector trackP4(const reco::Track& track, double mass) {
  TLorentzVector result;
  result.SetXYZM(track.px(), track.py(), track.pz(), mass);
  return result;
}

const reco::Candidate* lastCopy(const reco::Candidate* particle) {
  if (particle == nullptr)
    return nullptr;

  std::set<const reco::Candidate*> visited;
  const reco::Candidate* current = particle;
  while (visited.insert(current).second) {
    const reco::Candidate* next = nullptr;
    for (std::size_t i = 0; i < current->numberOfDaughters(); ++i) {
      const reco::Candidate* daughter = current->daughter(i);
      if (daughter != nullptr && daughter->pdgId() == current->pdgId()) {
        next = daughter;
        break;
      }
    }
    if (next == nullptr)
      break;
    current = next;
  }
  return current;
}

std::vector<const reco::Candidate*> collapsedDaughters(const reco::Candidate& parent) {
  std::vector<const reco::Candidate*> daughters;
  std::set<const reco::Candidate*> unique;
  for (std::size_t i = 0; i < parent.numberOfDaughters(); ++i) {
    const reco::Candidate* daughter = lastCopy(parent.daughter(i));
    if (daughter != nullptr && unique.insert(daughter).second)
      daughters.push_back(daughter);
  }
  return daughters;
}

bool hasFlavor(int pdgId, int flavor) {
  int id = std::abs(pdgId);
  if (id == flavor)
    return true;
  if (id < 100 || id >= 1000000000)
    return false;

  id /= 10;  // Remove the spin digit before inspecting hadron constituents.
  for (unsigned int digit = 0; digit < 3 && id > 0; ++digit) {
    if (id % 10 == flavor)
      return true;
    id /= 10;
  }
  return false;
}

void collectAncestors(const reco::Candidate* particle,
                      std::set<const reco::Candidate*>& visited,
                      std::vector<const reco::Candidate*>& ancestors) {
  if (particle == nullptr || !visited.insert(particle).second)
    return;

  for (std::size_t i = 0; i < particle->numberOfMothers(); ++i) {
    const reco::Candidate* mother = particle->mother(i);
    if (mother == nullptr || visited.count(mother) != 0)
      continue;
    ancestors.push_back(mother);
    collectAncestors(mother, visited, ancestors);
  }
}

int firstDistinctMotherPdgId(const reco::Candidate& particle) {
  const reco::Candidate* current = &particle;
  std::set<const reco::Candidate*> visited;
  while (current != nullptr && visited.insert(current).second) {
    const reco::Candidate* samePdgMother = nullptr;
    const reco::Candidate* distinctMother = nullptr;
    for (std::size_t i = 0; i < current->numberOfMothers(); ++i) {
      const reco::Candidate* mother = current->mother(i);
      if (mother == nullptr)
        continue;
      if (mother->pdgId() == particle.pdgId())
        samePdgMother = mother;
      else if (distinctMother == nullptr)
        distinctMother = mother;
    }
    if (distinctMother != nullptr)
      return distinctMother->pdgId();
    current = samePdgMother;
  }
  return 0;
}

void classifyOrigin(GenChain& chain) {
  std::set<const reco::Candidate*> visited;
  std::vector<const reco::Candidate*> ancestors;
  collectAncestors(chain.dst, visited, ancestors);
  for (const reco::Candidate* ancestor : ancestors) {
    if (ancestor->pdgId() == chain.dst->pdgId())
      continue;
    const int id = std::abs(ancestor->pdgId());
    chain.fromB = chain.fromB || hasFlavor(id, 5);
    chain.fromW = chain.fromW || id == 24;
    chain.fromCharm = chain.fromCharm || (id > 100 && hasFlavor(id, 4));
  }
  // Prompt charm includes direct production and feed-down from excited charm
  // states.  Beauty feed-down is the non-prompt category; the other flags are
  // retained independently to characterize the production chain.
  chain.prompt = !chain.fromB;
}

bool buildStrictChain(const reco::Candidate& inputDst, GenChain& chain) {
  const reco::Candidate* dst = lastCopy(&inputDst);
  if (dst == nullptr || std::abs(dst->pdgId()) != 413)
    return false;

  const int expectedD0Id = dst->pdgId() > 0 ? 421 : -421;
  const int expectedSoftPionId = dst->pdgId() > 0 ? 211 : -211;
  const reco::Candidate* d0 = nullptr;
  const reco::Candidate* softPion = nullptr;
  unsigned int chargedDstDaughters = 0;
  for (const reco::Candidate* daughter : collapsedDaughters(*dst)) {
    if (candidateCharge(*daughter) != 0)
      ++chargedDstDaughters;
    if (daughter->pdgId() == expectedD0Id) {
      if (d0 != nullptr)
        return false;
      d0 = daughter;
    } else if (daughter->pdgId() == expectedSoftPionId) {
      if (softPion != nullptr)
        return false;
      softPion = daughter;
    } else if (candidateCharge(*daughter) != 0) {
      return false;
    }
  }
  if (d0 == nullptr || softPion == nullptr || chargedDstDaughters != 1)
    return false;

  d0 = lastCopy(d0);
  const int expectedKaonId = d0->pdgId() > 0 ? -321 : 321;
  const int expectedPionId = d0->pdgId() > 0 ? 211 : -211;
  const reco::Candidate* kaon = nullptr;
  const reco::Candidate* pion = nullptr;
  unsigned int chargedD0Daughters = 0;
  for (const reco::Candidate* daughter : collapsedDaughters(*d0)) {
    if (candidateCharge(*daughter) != 0)
      ++chargedD0Daughters;
    if (daughter->pdgId() == expectedKaonId) {
      if (kaon != nullptr)
        return false;
      kaon = daughter;
    } else if (daughter->pdgId() == expectedPionId) {
      if (pion != nullptr)
        return false;
      pion = daughter;
    } else if (candidateCharge(*daughter) != 0) {
      return false;
    }
  }
  if (kaon == nullptr || pion == nullptr || chargedD0Daughters != 2)
    return false;

  chain.dst = dst;
  chain.d0 = d0;
  chain.kaon = lastCopy(kaon);
  chain.pion = lastCopy(pion);
  chain.softPion = lastCopy(softPion);
  chain.dstMotherPdgId = firstDistinctMotherPdgId(*dst);
  classifyOrigin(chain);
  return true;
}

std::vector<TrackMatchOption> compatibleTracks(const reco::Candidate& particle,
                                               const reco::TrackCollection& tracks,
                                               double maxDeltaR) {
  std::vector<TrackMatchOption> options;
  const int charge = candidateCharge(particle);
  for (unsigned int index = 0; index < tracks.size(); ++index) {
    if (tracks[index].charge() != charge)
      continue;
    const double dr = deltaR(particle, tracks[index]);
    if (dr <= maxDeltaR)
      options.push_back({index, dr});
  }
  return options;
}

MatchResult bestMatchForChain(const GenChain& chain,
                              std::size_t chainIndex,
                              const reco::TrackCollection& tracks,
                              double maxDeltaR) {
  const std::array<std::vector<TrackMatchOption>, 3> options = {
      {compatibleTracks(*chain.kaon, tracks, maxDeltaR),
       compatibleTracks(*chain.pion, tracks, maxDeltaR),
       compatibleTracks(*chain.softPion, tracks, maxDeltaR)}};

  MatchResult best;
  best.chainIndex = chainIndex;
  for (const auto& kaon : options[0]) {
    for (const auto& pion : options[1]) {
      if (pion.index == kaon.index)
        continue;
      for (const auto& softPion : options[2]) {
        if (softPion.index == kaon.index || softPion.index == pion.index)
          continue;
        const double score = kaon.deltaR * kaon.deltaR + pion.deltaR * pion.deltaR +
                             softPion.deltaR * softPion.deltaR;
        if (!best.valid || score < best.score) {
          best.valid = true;
          best.trackIndices = {{kaon.index, pion.index, softPion.index}};
          best.deltaRs = {{kaon.deltaR, pion.deltaR, softPion.deltaR}};
          best.score = score;
        }
      }
    }
  }
  return best;
}

TruthParticleRecord makeTruthRecord(const reco::Candidate& particle) {
  TruthParticleRecord record;
  record.pdgId = particle.pdgId();
  record.status = particle.status();
  record.charge = candidateCharge(particle);
  record.pt = particle.pt();
  record.eta = particle.eta();
  record.phi = particle.phi();
  record.mass = particle.mass();
  record.vx = particle.vx();
  record.vy = particle.vy();
  record.vz = particle.vz();
  return record;
}

RecoTrackRecord makeRecoRecord(const reco::Track& track,
                               unsigned int index,
                               double dr,
                               double mass,
                               const reco::Vertex* primaryVertex,
                               const reco::BeamSpot* beamSpot) {
  RecoTrackRecord record;
  record.index = static_cast<int>(index);
  record.charge = track.charge();
  record.deltaR = dr;
  record.pt = track.pt();
  record.eta = track.eta();
  record.phi = track.phi();
  record.mass = mass;
  record.px = track.px();
  record.py = track.py();
  record.pz = track.pz();
  record.vx = track.vx();
  record.vy = track.vy();
  record.vz = track.vz();
  record.normalizedChi2 = track.normalizedChi2();
  record.nValidHits = track.numberOfValidHits();
  if (beamSpot != nullptr) {
    record.dxyBeamSpot = track.dxy(beamSpot->position());
    record.dzBeamSpot = track.dz(beamSpot->position());
  }
  if (primaryVertex != nullptr) {
    record.dxyPrimaryVertex = track.dxy(primaryVertex->position());
    record.dzPrimaryVertex = track.dz(primaryVertex->position());
  }
  return record;
}

FitRecord fitRecord(RefCountedKinematicTree tree) {
  FitRecord record;
  if (!tree || !tree->isValid())
    return record;

  tree->movePointerToTheTop();
  const RefCountedKinematicVertex vertex = tree->currentDecayVertex();
  const RefCountedKinematicParticle particle = tree->currentParticle();
  if (!vertex || !vertex->vertexIsValid() || !particle)
    return record;

  record.valid = true;
  record.chi2 = vertex->chiSquared();
  record.dof = vertex->degreesOfFreedom();
  record.mass = particle->currentState().mass();
  if (record.chi2 >= 0.f && record.dof > 0.f)
    record.pval = TMath::Prob(record.chi2, static_cast<int>(std::lround(record.dof)));
  return record;
}
}  // namespace

class DstToD0PiMCTruthProducer : public edm::EDProducer {
public:
  explicit DstToD0PiMCTruthProducer(const edm::ParameterSet&);
  ~DstToD0PiMCTruthProducer() override;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  RefCountedKinematicTree fitD0(const edm::EventSetup&, const reco::Track&, const reco::Track&) const;
  RefCountedKinematicTree fitDst(const edm::EventSetup&,
                                 const reco::Track&,
                                 const RefCountedKinematicParticle&) const;
  void bookEventTree();
  void bookMatchedCandidateTree();
  void bookTruthParticle(TTree*, const std::string&, TruthParticleRecord&);
  void bookRecoTrack(const std::string&, RecoTrackRecord&);
  void resetEventRecord(const edm::Event&);
  void resetCandidateRecord();

  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> genParticleToken_;

  double maxDeltaR_ = 0.03;
  double hadPtMin_ = 0.6;
  // A negative value preserves the historical common threshold for the soft pion.
  double softPionPtMin_ = -1.0;
  double pvalChi2VtxMin_ = 0.005;
  double dmD0Max_ = 0.05;
  double dmDstMax_ = 0.15;
  double dmDstMinusD0Max_ = 0.003;
  double fiducialEtaMax_ = 2.5;

  TTree* eventTree_ = nullptr;
  TTree* matchedCandidateTree_ = nullptr;

  unsigned int eventRun_ = 0;
  unsigned int eventLumi_ = 0;
  unsigned long long eventNumber_ = 0;
  int eventNInputTracks_ = 0;
  int eventNGenChains_ = 0;
  int eventNFiducialGenChains_ = 0;
  int eventNMatchedChains_ = 0;
  bool eventHasGenParticles_ = false;
  bool eventHasGenChain_ = false;
  bool eventHasFiducialGenChain_ = false;
  bool eventHasUniqueMatch_ = false;
  bool eventD0FitValid_ = false;
  bool eventDstFitValid_ = false;
  bool eventPassDataSelection_ = false;

  unsigned int candidateRun_ = 0;
  unsigned int candidateLumi_ = 0;
  unsigned long long candidateEvent_ = 0;
  int candidateNGenChains_ = 0;
  int candidateNMatchedChains_ = 0;
  float candidateMatchScore_ = kFloatSentinel;
  int candidateDstMotherPdgId_ = 0;
  bool candidatePrompt_ = false;
  bool candidateFromB_ = false;
  bool candidateFromW_ = false;
  bool candidateFromCharm_ = false;
  TruthParticleRecord truthDst_;
  TruthParticleRecord truthD0_;
  TruthParticleRecord truthKaon_;
  TruthParticleRecord truthPion_;
  TruthParticleRecord truthSoftPion_;
  RecoTrackRecord recoKaon_;
  RecoTrackRecord recoPion_;
  RecoTrackRecord recoSoftPion_;
  float rawD0Mass_ = kFloatSentinel;
  float rawDstMass_ = kFloatSentinel;
  float rawDeltaMass_ = kFloatSentinel;
  FitRecord d0Fit_;
  bool dstFitAvailable_ = false;
  FitRecord dstFit_;
  bool passTrackPt_ = false;
  bool passD0FitValid_ = false;
  bool passD0Pval_ = false;
  bool passD0Mass_ = false;
  bool passDstFitValid_ = false;
  bool passDstPval_ = false;
  bool passDstMass_ = false;
  bool passDeltaMass_ = false;
  bool passDataSelection_ = false;

  unsigned long long counterEvents_ = 0;
  unsigned long long counterGenChains_ = 0;
  unsigned long long counterMatchedChains_ = 0;
  unsigned long long counterMatchedEvents_ = 0;
  unsigned long long counterD0FitValid_ = 0;
  unsigned long long counterDstFitValid_ = 0;
  unsigned long long counterFullSelection_ = 0;
};

DstToD0PiMCTruthProducer::DstToD0PiMCTruthProducer(const edm::ParameterSet& config)
    : trackToken_(consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("src"))),
      vertexToken_(consumes<reco::VertexCollection>(config.getParameter<edm::InputTag>("vertices"))),
      beamSpotToken_(consumes<reco::BeamSpot>(config.getParameter<edm::InputTag>("beamSpot"))),
      genParticleToken_(
          consumes<edm::View<reco::Candidate>>(config.getParameter<edm::InputTag>("genParticles"))) {
  if (config.existsAs<double>("maxDeltaR"))
    maxDeltaR_ = config.getParameter<double>("maxDeltaR");
  if (config.existsAs<double>("hadPtMin"))
    hadPtMin_ = config.getParameter<double>("hadPtMin");
  if (config.existsAs<double>("softPionPtMin"))
    softPionPtMin_ = config.getParameter<double>("softPionPtMin");
  if (config.existsAs<double>("pvalChi2VtxMin"))
    pvalChi2VtxMin_ = config.getParameter<double>("pvalChi2VtxMin");
  if (config.existsAs<double>("dmD0Max"))
    dmD0Max_ = config.getParameter<double>("dmD0Max");
  if (config.existsAs<double>("dmDstMax"))
    dmDstMax_ = config.getParameter<double>("dmDstMax");
  if (config.existsAs<double>("dmDstMinusD0Max"))
    dmDstMinusD0Max_ = config.getParameter<double>("dmDstMinusD0Max");
  if (config.existsAs<double>("fiducialEtaMax"))
    fiducialEtaMax_ = config.getParameter<double>("fiducialEtaMax");

  produces<reco::TrackCollection>("truthMatchedDstarTracks");

  edm::Service<TFileService> fileService;
  eventTree_ = fileService->make<TTree>("Events", "MC D* truth matching event summary");
  matchedCandidateTree_ =
      fileService->make<TTree>("MatchedCandidates", "MC truth-matched D* candidates; filled downstream");
  bookEventTree();
  bookMatchedCandidateTree();
}

DstToD0PiMCTruthProducer::~DstToD0PiMCTruthProducer() {
  std::cout << "DstToD0PiMCTruthProducer counters:"
            << " events=" << counterEvents_ << " genChains=" << counterGenChains_
            << " matchedChains=" << counterMatchedChains_ << " matchedEvents=" << counterMatchedEvents_
            << " d0FitValid=" << counterD0FitValid_ << " dstFitValid=" << counterDstFitValid_
            << " fullSelection=" << counterFullSelection_ << std::endl;
}

void DstToD0PiMCTruthProducer::bookEventTree() {
  eventTree_->Branch("MC_run", &eventRun_);
  eventTree_->Branch("MC_lumi", &eventLumi_);
  eventTree_->Branch("MC_event", &eventNumber_);
  eventTree_->Branch("MC_nInputTracks", &eventNInputTracks_);
  eventTree_->Branch("truth_nGenChains", &eventNGenChains_);
  eventTree_->Branch("truth_nFiducialGenChains", &eventNFiducialGenChains_);
  eventTree_->Branch("MC_nMatchedChains", &eventNMatchedChains_);
  eventTree_->Branch("truth_hasGenParticles", &eventHasGenParticles_);
  eventTree_->Branch("truth_hasGenChain", &eventHasGenChain_);
  eventTree_->Branch("truth_hasFiducialGenChain", &eventHasFiducialGenChain_);
  eventTree_->Branch("MC_hasUniqueMatch", &eventHasUniqueMatch_);
  eventTree_->Branch("MC_D0_fit_valid", &eventD0FitValid_);
  eventTree_->Branch("MC_Dst_fit_valid", &eventDstFitValid_);
  eventTree_->Branch("MC_passDataSelection", &eventPassDataSelection_);
  eventTree_->Branch("truth_Dst_motherPdgId", &candidateDstMotherPdgId_);
  eventTree_->Branch("truth_Dst_prompt", &candidatePrompt_);
  eventTree_->Branch("truth_Dst_fromB", &candidateFromB_);
  eventTree_->Branch("truth_Dst_fromW", &candidateFromW_);
  eventTree_->Branch("truth_Dst_fromCharm", &candidateFromCharm_);
  bookTruthParticle(eventTree_, "truth_Dst_", truthDst_);
  bookTruthParticle(eventTree_, "truth_D0_", truthD0_);
  bookTruthParticle(eventTree_, "truth_K_", truthKaon_);
  bookTruthParticle(eventTree_, "truth_pi_", truthPion_);
  bookTruthParticle(eventTree_, "truth_pi_s_", truthSoftPion_);
}

void DstToD0PiMCTruthProducer::bookTruthParticle(TTree* tree,
                                                 const std::string& prefix,
                                                 TruthParticleRecord& record) {
  tree->Branch((prefix + "pdgId").c_str(), &record.pdgId);
  tree->Branch((prefix + "status").c_str(), &record.status);
  tree->Branch((prefix + "charge").c_str(), &record.charge);
  tree->Branch((prefix + "pt").c_str(), &record.pt);
  tree->Branch((prefix + "eta").c_str(), &record.eta);
  tree->Branch((prefix + "phi").c_str(), &record.phi);
  tree->Branch((prefix + "mass").c_str(), &record.mass);
  tree->Branch((prefix + "vx").c_str(), &record.vx);
  tree->Branch((prefix + "vy").c_str(), &record.vy);
  tree->Branch((prefix + "vz").c_str(), &record.vz);
}

void DstToD0PiMCTruthProducer::bookRecoTrack(const std::string& prefix, RecoTrackRecord& record) {
  matchedCandidateTree_->Branch((prefix + "index").c_str(), &record.index);
  matchedCandidateTree_->Branch((prefix + "charge").c_str(), &record.charge);
  matchedCandidateTree_->Branch((prefix + "deltaR").c_str(), &record.deltaR);
  matchedCandidateTree_->Branch((prefix + "pt").c_str(), &record.pt);
  matchedCandidateTree_->Branch((prefix + "eta").c_str(), &record.eta);
  matchedCandidateTree_->Branch((prefix + "phi").c_str(), &record.phi);
  matchedCandidateTree_->Branch((prefix + "mass").c_str(), &record.mass);
  matchedCandidateTree_->Branch((prefix + "px").c_str(), &record.px);
  matchedCandidateTree_->Branch((prefix + "py").c_str(), &record.py);
  matchedCandidateTree_->Branch((prefix + "pz").c_str(), &record.pz);
  matchedCandidateTree_->Branch((prefix + "vx").c_str(), &record.vx);
  matchedCandidateTree_->Branch((prefix + "vy").c_str(), &record.vy);
  matchedCandidateTree_->Branch((prefix + "vz").c_str(), &record.vz);
  matchedCandidateTree_->Branch((prefix + "dxyBeamSpot").c_str(), &record.dxyBeamSpot);
  matchedCandidateTree_->Branch((prefix + "dzBeamSpot").c_str(), &record.dzBeamSpot);
  matchedCandidateTree_->Branch((prefix + "dxyPrimaryVertex").c_str(), &record.dxyPrimaryVertex);
  matchedCandidateTree_->Branch((prefix + "dzPrimaryVertex").c_str(), &record.dzPrimaryVertex);
  matchedCandidateTree_->Branch((prefix + "normalizedChi2").c_str(), &record.normalizedChi2);
  matchedCandidateTree_->Branch((prefix + "nValidHits").c_str(), &record.nValidHits);
}

void DstToD0PiMCTruthProducer::bookMatchedCandidateTree() {
  matchedCandidateTree_->Branch("MC_run", &candidateRun_);
  matchedCandidateTree_->Branch("MC_lumi", &candidateLumi_);
  matchedCandidateTree_->Branch("MC_event", &candidateEvent_);
  matchedCandidateTree_->Branch("truth_nGenChains", &candidateNGenChains_);
  matchedCandidateTree_->Branch("MC_nMatchedChains", &candidateNMatchedChains_);
  matchedCandidateTree_->Branch("MC_matchScoreDeltaR2", &candidateMatchScore_);
  matchedCandidateTree_->Branch("truth_Dst_motherPdgId", &candidateDstMotherPdgId_);
  matchedCandidateTree_->Branch("truth_Dst_prompt", &candidatePrompt_);
  matchedCandidateTree_->Branch("truth_Dst_fromB", &candidateFromB_);
  matchedCandidateTree_->Branch("truth_Dst_fromW", &candidateFromW_);
  matchedCandidateTree_->Branch("truth_Dst_fromCharm", &candidateFromCharm_);

  bookTruthParticle(matchedCandidateTree_, "truth_Dst_", truthDst_);
  bookTruthParticle(matchedCandidateTree_, "truth_D0_", truthD0_);
  bookTruthParticle(matchedCandidateTree_, "truth_K_", truthKaon_);
  bookTruthParticle(matchedCandidateTree_, "truth_pi_", truthPion_);
  bookTruthParticle(matchedCandidateTree_, "truth_pi_s_", truthSoftPion_);
  bookRecoTrack("MC_reco_K_", recoKaon_);
  bookRecoTrack("MC_reco_pi_", recoPion_);
  bookRecoTrack("MC_reco_pi_s_", recoSoftPion_);

  matchedCandidateTree_->Branch("MC_raw_D0_mass", &rawD0Mass_);
  matchedCandidateTree_->Branch("MC_raw_Dst_mass", &rawDstMass_);
  matchedCandidateTree_->Branch("MC_raw_DeltaM", &rawDeltaMass_);
  matchedCandidateTree_->Branch("MC_D0_fit_valid", &d0Fit_.valid);
  matchedCandidateTree_->Branch("MC_D0_fit_chi2", &d0Fit_.chi2);
  matchedCandidateTree_->Branch("MC_D0_fit_dof", &d0Fit_.dof);
  matchedCandidateTree_->Branch("MC_D0_fit_pval", &d0Fit_.pval);
  matchedCandidateTree_->Branch("MC_D0_fit_mass", &d0Fit_.mass);
  matchedCandidateTree_->Branch("MC_Dst_fit_available", &dstFitAvailable_);
  matchedCandidateTree_->Branch("MC_Dst_fit_valid", &dstFit_.valid);
  matchedCandidateTree_->Branch("MC_Dst_fit_chi2", &dstFit_.chi2);
  matchedCandidateTree_->Branch("MC_Dst_fit_dof", &dstFit_.dof);
  matchedCandidateTree_->Branch("MC_Dst_fit_pval", &dstFit_.pval);
  matchedCandidateTree_->Branch("MC_Dst_fit_mass", &dstFit_.mass);
  matchedCandidateTree_->Branch("MC_pass_trackPt", &passTrackPt_);
  matchedCandidateTree_->Branch("MC_pass_D0_fitValid", &passD0FitValid_);
  matchedCandidateTree_->Branch("MC_pass_D0_pval", &passD0Pval_);
  matchedCandidateTree_->Branch("MC_pass_D0_mass", &passD0Mass_);
  matchedCandidateTree_->Branch("MC_pass_Dst_fitValid", &passDstFitValid_);
  matchedCandidateTree_->Branch("MC_pass_Dst_pval", &passDstPval_);
  matchedCandidateTree_->Branch("MC_pass_Dst_mass", &passDstMass_);
  matchedCandidateTree_->Branch("MC_pass_DstMinusD0", &passDeltaMass_);
  matchedCandidateTree_->Branch("MC_passDataSelection", &passDataSelection_);
}

void DstToD0PiMCTruthProducer::resetEventRecord(const edm::Event& event) {
  eventRun_ = event.run();
  eventLumi_ = event.luminosityBlock();
  eventNumber_ = event.id().event();
  eventNInputTracks_ = 0;
  eventNGenChains_ = 0;
  eventNFiducialGenChains_ = 0;
  eventNMatchedChains_ = 0;
  eventHasGenParticles_ = false;
  eventHasGenChain_ = false;
  eventHasFiducialGenChain_ = false;
  eventHasUniqueMatch_ = false;
  eventD0FitValid_ = false;
  eventDstFitValid_ = false;
  eventPassDataSelection_ = false;
}

void DstToD0PiMCTruthProducer::resetCandidateRecord() {
  candidateRun_ = eventRun_;
  candidateLumi_ = eventLumi_;
  candidateEvent_ = eventNumber_;
  candidateNGenChains_ = 0;
  candidateNMatchedChains_ = 0;
  candidateMatchScore_ = kFloatSentinel;
  candidateDstMotherPdgId_ = 0;
  candidatePrompt_ = false;
  candidateFromB_ = false;
  candidateFromW_ = false;
  candidateFromCharm_ = false;
  truthDst_ = TruthParticleRecord();
  truthD0_ = TruthParticleRecord();
  truthKaon_ = TruthParticleRecord();
  truthPion_ = TruthParticleRecord();
  truthSoftPion_ = TruthParticleRecord();
  recoKaon_ = RecoTrackRecord();
  recoPion_ = RecoTrackRecord();
  recoSoftPion_ = RecoTrackRecord();
  rawD0Mass_ = kFloatSentinel;
  rawDstMass_ = kFloatSentinel;
  rawDeltaMass_ = kFloatSentinel;
  d0Fit_ = FitRecord();
  dstFitAvailable_ = false;
  dstFit_ = FitRecord();
  passTrackPt_ = false;
  passD0FitValid_ = false;
  passD0Pval_ = false;
  passD0Mass_ = false;
  passDstFitValid_ = false;
  passDstPval_ = false;
  passDstMass_ = false;
  passDeltaMass_ = false;
  passDataSelection_ = false;
}

RefCountedKinematicTree DstToD0PiMCTruthProducer::fitD0(const edm::EventSetup& setup,
                                                        const reco::Track& kaon,
                                                        const reco::Track& pion) const {
  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

  KinematicParticleFactoryFromTransientTrack particleFactory;
  std::vector<RefCountedKinematicParticle> particles;
  float chi2 = 0.f;
  float dof = 0.f;
  float kaonMassError = kMassError;
  float pionMassError = kMassError;
  particles.push_back(particleFactory.particle(trackBuilder->build(kaon), kKaonMass, chi2, dof, kaonMassError));
  particles.push_back(particleFactory.particle(trackBuilder->build(pion), kPionMass, chi2, dof, pionMassError));

  KinematicParticleVertexFitter fitter;
  return fitter.fit(particles);
}

RefCountedKinematicTree DstToD0PiMCTruthProducer::fitDst(
    const edm::EventSetup& setup,
    const reco::Track& softPion,
    const RefCountedKinematicParticle& d0) const {
  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

  KinematicParticleFactoryFromTransientTrack particleFactory;
  std::vector<RefCountedKinematicParticle> particles;
  float chi2 = 0.f;
  float dof = 0.f;
  float pionMassError = kMassError;
  float d0Mass = d0->currentState().mass();
  float d0MassError = std::sqrt(std::abs(d0->currentState().kinematicParametersError().matrix()(6, 6)));
  particles.push_back(
      particleFactory.particle(trackBuilder->build(softPion), kPionMass, chi2, dof, pionMassError));
  particles.push_back(
      particleFactory.particle(d0->refittedTransientTrack(), d0Mass, chi2, dof, d0MassError));

  KinematicParticleVertexFitter fitter;
  return fitter.fit(particles);
}

void DstToD0PiMCTruthProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  ++counterEvents_;
  resetEventRecord(event);
  resetCandidateRecord();
  auto outputTracks = std::make_unique<reco::TrackCollection>();

  edm::Handle<reco::TrackCollection> tracks;
  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<reco::BeamSpot> beamSpot;
  edm::Handle<edm::View<reco::Candidate>> genParticles;
  event.getByToken(trackToken_, tracks);
  event.getByToken(vertexToken_, vertices);
  event.getByToken(beamSpotToken_, beamSpot);
  event.getByToken(genParticleToken_, genParticles);

  if (tracks.isValid())
    eventNInputTracks_ = static_cast<int>(tracks->size());
  eventHasGenParticles_ = genParticles.isValid();

  std::vector<GenChain> chains;
  if (genParticles.isValid()) {
    std::set<const reco::Candidate*> acceptedLastCopies;
    for (const reco::Candidate& particle : *genParticles) {
      if (std::abs(particle.pdgId()) != 413)
        continue;
      const reco::Candidate* dst = lastCopy(&particle);
      if (dst == nullptr || acceptedLastCopies.count(dst) != 0)
        continue;
      GenChain chain;
      if (buildStrictChain(particle, chain)) {
        acceptedLastCopies.insert(chain.dst);
        chains.push_back(chain);
      }
    }
  }

  eventNGenChains_ = static_cast<int>(chains.size());
  eventHasGenChain_ = !chains.empty();
  candidateNGenChains_ = eventNGenChains_;
  counterGenChains_ += chains.size();
  const double effectiveSoftPionPtMin = softPionPtMin_ >= 0. ? softPionPtMin_ : hadPtMin_;
  for (const GenChain& chain : chains) {
    const bool fiducial = chain.kaon->pt() >= hadPtMin_ && chain.pion->pt() >= hadPtMin_ &&
                          chain.softPion->pt() >= effectiveSoftPionPtMin &&
                          std::abs(chain.kaon->eta()) <= fiducialEtaMax_ &&
                          std::abs(chain.pion->eta()) <= fiducialEtaMax_ &&
                          std::abs(chain.softPion->eta()) <= fiducialEtaMax_;
    if (fiducial)
      ++eventNFiducialGenChains_;
  }
  eventHasFiducialGenChain_ = eventNFiducialGenChains_ > 0;

  // Preserve a representative generated chain even when none of its daughters
  // can be matched to the ALCARECO track pool.  If a match is found below,
  // these records are replaced by the matched chain before either tree fills.
  if (!chains.empty()) {
    const auto representative = std::max_element(
        chains.begin(), chains.end(), [](const GenChain& first, const GenChain& second) {
          return first.dst->pt() < second.dst->pt();
        });
    candidateDstMotherPdgId_ = representative->dstMotherPdgId;
    candidatePrompt_ = representative->prompt;
    candidateFromB_ = representative->fromB;
    candidateFromW_ = representative->fromW;
    candidateFromCharm_ = representative->fromCharm;
    truthDst_ = makeTruthRecord(*representative->dst);
    truthD0_ = makeTruthRecord(*representative->d0);
    truthKaon_ = makeTruthRecord(*representative->kaon);
    truthPion_ = makeTruthRecord(*representative->pion);
    truthSoftPion_ = makeTruthRecord(*representative->softPion);
  }

  MatchResult globalBestMatch;
  if (tracks.isValid()) {
    for (std::size_t chainIndex = 0; chainIndex < chains.size(); ++chainIndex) {
      const MatchResult chainMatch = bestMatchForChain(chains[chainIndex], chainIndex, *tracks, maxDeltaR_);
      if (!chainMatch.valid)
        continue;
      ++eventNMatchedChains_;
      if (!globalBestMatch.valid || chainMatch.score < globalBestMatch.score)
        globalBestMatch = chainMatch;
    }
  }

  candidateNMatchedChains_ = eventNMatchedChains_;
  counterMatchedChains_ += static_cast<unsigned long long>(eventNMatchedChains_);
  eventHasUniqueMatch_ = globalBestMatch.valid;
  if (globalBestMatch.valid) {
    ++counterMatchedEvents_;
    const GenChain& chain = chains[globalBestMatch.chainIndex];
    const reco::Vertex* primaryVertex = nullptr;
    if (vertices.isValid() && !vertices->empty()) {
      const auto firstGoodVertex = std::find_if(
          vertices->begin(), vertices->end(), [](const reco::Vertex& vertex) { return !vertex.isFake(); });
      primaryVertex = firstGoodVertex != vertices->end() ? &*firstGoodVertex : &vertices->front();
    }
    const reco::BeamSpot* validBeamSpot = beamSpot.isValid() ? &*beamSpot : nullptr;

    const reco::Track& kaon = tracks->at(globalBestMatch.trackIndices[0]);
    const reco::Track& pion = tracks->at(globalBestMatch.trackIndices[1]);
    const reco::Track& softPion = tracks->at(globalBestMatch.trackIndices[2]);
    outputTracks->push_back(kaon);
    outputTracks->push_back(pion);
    outputTracks->push_back(softPion);

    candidateMatchScore_ = globalBestMatch.score;
    candidateDstMotherPdgId_ = chain.dstMotherPdgId;
    candidatePrompt_ = chain.prompt;
    candidateFromB_ = chain.fromB;
    candidateFromW_ = chain.fromW;
    candidateFromCharm_ = chain.fromCharm;
    truthDst_ = makeTruthRecord(*chain.dst);
    truthD0_ = makeTruthRecord(*chain.d0);
    truthKaon_ = makeTruthRecord(*chain.kaon);
    truthPion_ = makeTruthRecord(*chain.pion);
    truthSoftPion_ = makeTruthRecord(*chain.softPion);
    recoKaon_ = makeRecoRecord(
        kaon, globalBestMatch.trackIndices[0], globalBestMatch.deltaRs[0], kKaonMass, primaryVertex, validBeamSpot);
    recoPion_ = makeRecoRecord(
        pion, globalBestMatch.trackIndices[1], globalBestMatch.deltaRs[1], kPionMass, primaryVertex, validBeamSpot);
    recoSoftPion_ = makeRecoRecord(softPion,
                                          globalBestMatch.trackIndices[2],
                                          globalBestMatch.deltaRs[2],
                                          kPionMass,
                                          primaryVertex,
                                          validBeamSpot);

    const TLorentzVector rawD0 = trackP4(kaon, kKaonMass) + trackP4(pion, kPionMass);
    const TLorentzVector rawDst = rawD0 + trackP4(softPion, kPionMass);
    rawD0Mass_ = rawD0.M();
    rawDstMass_ = rawDst.M();
    rawDeltaMass_ = rawDstMass_ - rawD0Mass_;

    RefCountedKinematicTree d0Tree;
    try {
      d0Tree = fitD0(setup, kaon, pion);
    } catch (const std::exception&) {
      d0Tree = RefCountedKinematicTree();
    }
    d0Fit_ = fitRecord(d0Tree);
    if (d0Fit_.valid) {
      ++counterD0FitValid_;
      d0Tree->movePointerToTheTop();
      dstFitAvailable_ = true;
      RefCountedKinematicTree dstTree;
      try {
        dstTree = fitDst(setup, softPion, d0Tree->currentParticle());
      } catch (const std::exception&) {
        dstTree = RefCountedKinematicTree();
      }
      dstFit_ = fitRecord(dstTree);
      if (dstFit_.valid)
        ++counterDstFitValid_;
    }

    passTrackPt_ = kaon.pt() >= hadPtMin_ && pion.pt() >= hadPtMin_ &&
                   softPion.pt() >= effectiveSoftPionPtMin;
    passD0FitValid_ = d0Fit_.valid;
    passD0Pval_ = d0Fit_.valid && d0Fit_.pval > pvalChi2VtxMin_;
    passD0Mass_ = d0Fit_.valid && std::abs(d0Fit_.mass - kD0Mass) <= dmD0Max_;
    passDstFitValid_ = dstFit_.valid;
    passDstPval_ = dstFit_.valid && dstFit_.pval > pvalChi2VtxMin_;
    passDstMass_ = dstFit_.valid && std::abs(dstFit_.mass - kDstMass) <= dmDstMax_;
    passDeltaMass_ = d0Fit_.valid && dstFit_.valid &&
                     std::abs((dstFit_.mass - d0Fit_.mass) - kNominalDeltaMass) <= dmDstMinusD0Max_;
    passDataSelection_ = passTrackPt_ && passD0FitValid_ && passD0Pval_ && passD0Mass_ &&
                         passDstFitValid_ && passDstPval_ && passDstMass_ && passDeltaMass_;
    if (passDataSelection_)
      ++counterFullSelection_;
  }

  eventD0FitValid_ = d0Fit_.valid;
  eventDstFitValid_ = dstFit_.valid;
  eventPassDataSelection_ = passDataSelection_;
  eventTree_->Fill();
  event.put(std::move(outputTracks), "truthMatchedDstarTracks");
}

DEFINE_FWK_MODULE(DstToD0PiMCTruthProducer);

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <unordered_set>
#include <vector>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

namespace {
struct FitQuality {
  bool isValid = false;
  double chi2 = -1.;
  double dof = -1.;
  double pval = -1.;
};

struct CheapLine {
  double x = 0.;
  double y = 0.;
  double z = 0.;
  double ux = 0.;
  double uy = 0.;
  double uz = 0.;
};

struct CheapCAResult {
  double dxy = 999.;
  double dz = 999.;
  double dist3d = 999.;
  double midX = 0.;
  double midY = 0.;
  double midZ = 0.;
};

struct DistanceWithError {
  double value = -1.;
  double error = -1.;
  double significance = -1.;
};

CheapLine makeCheapLine(const reco::Track& track) {
  const double norm = std::max(track.p(), 1.e-12);
  return {track.vx(), track.vy(), track.vz(), track.px() / norm, track.py() / norm, track.pz() / norm};
}

CheapLine makeCheapLine(double x, double y, double z, double px, double py, double pz) {
  const double norm = std::max(std::sqrt(px * px + py * py + pz * pz), 1.e-12);
  return {x, y, z, px / norm, py / norm, pz / norm};
}

CheapCAResult cheapClosestApproachXY(const CheapLine& a, const CheapLine& b) {
  const double dx = b.x - a.x;
  const double dy = b.y - a.y;
  const double aa = a.ux * a.ux + a.uy * a.uy;
  const double bb = b.ux * b.ux + b.uy * b.uy;
  const double ab = a.ux * b.ux + a.uy * b.uy;
  const double ad = a.ux * dx + a.uy * dy;
  const double bd = b.ux * dx + b.uy * dy;
  const double det = aa * bb - ab * ab;

  double sa = 0.;
  double sb = 0.;
  if (std::abs(det) > 1.e-12) {
    sa = (ad * bb - bd * ab) / det;
    sb = (ad * ab - bd * aa) / det;
  }

  const double ax = a.x + sa * a.ux;
  const double ay = a.y + sa * a.uy;
  const double az = a.z + sa * a.uz;
  const double bx = b.x + sb * b.ux;
  const double by = b.y + sb * b.uy;
  const double bz = b.z + sb * b.uz;
  const double dxy = std::sqrt((ax - bx) * (ax - bx) + (ay - by) * (ay - by));
  const double dz = az - bz;
  return {dxy, dz, std::sqrt(dxy * dxy + dz * dz), 0.5 * (ax + bx), 0.5 * (ay + by), 0.5 * (az + bz)};
}

double deltaPhi(double phi1, double phi2) {
  double dphi = phi1 - phi2;
  while (dphi > M_PI)
    dphi -= 2. * M_PI;
  while (dphi <= -M_PI)
    dphi += 2. * M_PI;
  return dphi;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
  const double dphi = deltaPhi(phi1, phi2);
  const double deta = eta1 - eta2;
  return std::sqrt(deta * deta + dphi * dphi);
}

double safeSignificance(double value, double error) {
  return error > 0. ? std::abs(value) / error : -1.;
}

double safeTrackDzSignificance(const reco::Track& track, const reco::Vertex& pv) {
  return safeSignificance(track.dz(pv.position()), track.dzError());
}

double safeTrackDxySignificance(const reco::Track& track, const reco::Vertex& pv) {
  return safeSignificance(track.dxy(pv.position()), track.dxyError());
}

double pointingCos3D(const reco::Vertex& pv, const RefCountedKinematicVertex& vertex, const TLorentzVector& p4) {
  const double dx = vertex->position().x() - pv.x();
  const double dy = vertex->position().y() - pv.y();
  const double dz = vertex->position().z() - pv.z();
  const double l = std::sqrt(dx * dx + dy * dy + dz * dz);
  const double p = p4.P();
  if (l <= 0. || p <= 0.)
    return -99.;
  return (dx * p4.Px() + dy * p4.Py() + dz * p4.Pz()) / (l * p);
}

double pointingCos2D(const reco::Vertex& pv, const RefCountedKinematicVertex& vertex, const TLorentzVector& p4) {
  const double dx = vertex->position().x() - pv.x();
  const double dy = vertex->position().y() - pv.y();
  const double l = std::sqrt(dx * dx + dy * dy);
  const double p = p4.Pt();
  if (l <= 0. || p <= 0.)
    return -99.;
  return (dx * p4.Px() + dy * p4.Py()) / (l * p);
}

DistanceWithError vertexDistance3D(const reco::Vertex& pv, const RefCountedKinematicVertex& vertex) {
  DistanceWithError result;
  const double dx = vertex->position().x() - pv.x();
  const double dy = vertex->position().y() - pv.y();
  const double dz = vertex->position().z() - pv.z();
  result.value = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (result.value <= 0.)
    return result;

  const double direction[3] = {dx / result.value, dy / result.value, dz / result.value};
  const auto pvCov = pv.covariance();
  const auto vtxCov = vertex->error().matrix();
  double error2 = 0.;
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      error2 += direction[i] * (pvCov.At(i, j) + vtxCov.At(i, j)) * direction[j];
    }
  }
  result.error = error2 > 0. ? std::sqrt(error2) : -1.;
  result.significance = safeSignificance(result.value, result.error);
  return result;
}

double vertexTransverseDistance(const reco::Vertex& pv, const RefCountedKinematicVertex& vertex) {
  const double dx = vertex->position().x() - pv.x();
  const double dy = vertex->position().y() - pv.y();
  return std::sqrt(dx * dx + dy * dy);
}

double vertexTransverseDistanceFromBeamSpot(const reco::BeamSpot& beamSpot, const RefCountedKinematicVertex& vertex) {
  const double z = vertex->position().z();
  const double bx = beamSpot.x(z);
  const double by = beamSpot.y(z);
  const double dx = vertex->position().x() - bx;
  const double dy = vertex->position().y() - by;
  return std::sqrt(dx * dx + dy * dy);
}

TLorentzVector p4(const reco::Track& track, double mass) {
  TLorentzVector v;
  v.SetXYZM(track.px(), track.py(), track.pz(), mass);
  return v;
}

constexpr double kNominalDeltaM = 0.145425;

unsigned long long pairKey(unsigned int first, unsigned int second) {
  if (second < first)
    std::swap(first, second);
  return (static_cast<unsigned long long>(first) << 32) | second;
}

struct SelectedCandidate {
  unsigned int kaonIndex = 0;
  unsigned int pionIndex = 0;
  unsigned int softPionIndex = 0;
  reco::Track kaon;
  reco::Track pion;
  float rawD0Mass = 0.f;
  float rawD0MassKK = 0.f;
  float rawD0MassKpiConj = 0.f;
  float rawDstMass = 0.f;
  float rawDeltaM = 0.f;
  float dstPt = 0.f;
  float dstConePtSum = 0.f;
  float dstIso = 0.f;
  float pisDRD0Raw = 0.f;
  float cheapCAScore = 0.f;
  float deltaDxyBSKpi = 0.f;
  float deltaDzBSKpi = 0.f;
  int isSameChargeD0 = 0;
  int isOppositeChargeD0 = 0;
  float cheapCAKpiDxy = 0.f;
  float cheapCAKpiDz = 0.f;
  float cheapCAKpi3D = 0.f;
  float cheapCAD0PisDxy = 0.f;
  float cheapCAD0PisDz = 0.f;
  float cheapCAD0Pis3D = 0.f;
  int KCharge = 0;
  int piCharge = 0;
  int pisCharge = 0;
  float KPt = 0.f;
  float KEta = 0.f;
  float KPhi = 0.f;
  float KE = 0.f;
  int KPdgId = 0;
  float KNormChi2 = 0.f;
  float KValidHits = 0.f;
  float KPixelHits = 0.f;
  float KLostInnerHits = 0.f;
  float KDxyBS = 0.f;
  float KDxyErrBS = 0.f;
  float KDxySigBS = 0.f;
  float KDxyErrPV = 0.f;
  float KDxyPV = 0.f;
  float KDzPV = 0.f;
  float KDxySigPV = 0.f;
  float KDzSigPV = 0.f;
  float piPt = 0.f;
  float piEta = 0.f;
  float piPhi = 0.f;
  float piE = 0.f;
  int piPdgId = 0;
  float piNormChi2 = 0.f;
  float piValidHits = 0.f;
  float piPixelHits = 0.f;
  float piLostInnerHits = 0.f;
  float piDxyBS = 0.f;
  float piDxyErrBS = 0.f;
  float piDxySigBS = 0.f;
  float piDxyErrPV = 0.f;
  float piDxyPV = 0.f;
  float piDzPV = 0.f;
  float piDxySigPV = 0.f;
  float piDzSigPV = 0.f;
  float pisPt = 0.f;
  float pisEta = 0.f;
  float pisPhi = 0.f;
  float pisE = 0.f;
  int pisPdgId = 0;
  float pisNormChi2 = 0.f;
  float pisValidHits = 0.f;
  float pisPixelHits = 0.f;
  float pisLostInnerHits = 0.f;
  float pisDxyBS = 0.f;
  float pisDxyErrBS = 0.f;
  float pisDxySigBS = 0.f;
  float pisDxyErrPV = 0.f;
  float pisDxyPV = 0.f;
  float pisDzPV = 0.f;
  float pisDxySigPV = 0.f;
  float pisDzSigPV = 0.f;
  float d0Pt = 0.f;
  float d0Eta = 0.f;
  float d0Phi = 0.f;
  float d0FitVertexProb = 0.f;
  float d0FitCos3D = 0.f;
  float d0FitCos2D = 0.f;
  float d0FitDecayLength3D = 0.f;
  float d0FitDecayLength3DErr = 0.f;
  float d0FitDecayLength3DSig = 0.f;
  float d0FitMass = 0.f;
  float d0FitChi2 = 0.f;
  float d0FitDof = 0.f;
  float d0FitPval = 0.f;
  float d0FitDxyBS = -99.f;
  float d0FitDxySigBS = -99.f;
  float d0FitTransverseDecayLength = -99.f;
  float d0FitTransverseDecayLengthSig = -99.f;
  float dstFitMass = 0.f;
  float dstFitDeltaM = 0.f;
  float dstFitChi2 = 0.f;
  float dstFitDof = 0.f;
  float dstFitPval = 0.f;
  float dstFitCos3D = 0.f;
  float dstFitCos2D = 0.f;
  float dstFitDecayLength3D = 0.f;
  float dstFitDecayLength3DErr = 0.f;
  float dstFitDecayLength3DSig = 0.f;
  float dstFitTransverseDecayLength = 0.f;
  float dstFitTransverseDecayLengthSig = -99.f;
};
}  // namespace

class DstToD0PiCascadeRefitProducer : public edm::one::EDProducer<edm::one::SharedResources> {
public:
  explicit DstToD0PiCascadeRefitProducer(const edm::ParameterSet&);
  ~DstToD0PiCascadeRefitProducer() override;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  bool checkCharge(const reco::Track& kaon, const reco::Track& pion, const reco::Track& softPion) const;
  FitQuality fitQuality(RefCountedKinematicTree tree) const;
  RefCountedKinematicTree fitD0(const edm::EventSetup&, const reco::Track&, const reco::Track&) const;
  RefCountedKinematicTree fitDst(const edm::EventSetup&, const reco::Track&, const RefCountedKinematicParticle&) const;
  void clearVectors();
  void bookTree();

  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

  double kaonMass_;
  double pionMass_;
  double softPionMass_;
  double kaonMassErr_;
  double pionMassErr_;
  double hadPtMin_;
  double softPionPtMin_;
  double maxTrackDzPV_;
  double maxKpiDeltaR_;
  double minD0Mass_;
  double maxD0Mass_;
  double minDstMass_;
  double maxDstMass_;
  double minDeltaMass_;
  double maxDeltaMass_;
  double minD0FitPval_;
  double minD0FitMass_;
  double maxD0FitMass_;
  double minDstFitPval_;
  double minDstFitMass_;
  double maxDstFitMass_;
  double minDstFitDeltaMass_;
  double maxDstFitDeltaMass_;
  double maxCheapCAKpi3D_;
  double maxCheapCAD0Pis3D_;
  double minDstIso_;
  double minDstPt_;
  double maxSoftPiDRD0_;
  double dstIsoCone_;
  bool applyChargeFilter_;
  bool useUnsignedCharge_;
  bool applyUniqueSelection_;
  bool writeTree_;
  int targetCharge_;

  TTree* tree_ = nullptr;
  unsigned int run_ = 0;
  unsigned int lumi_ = 0;
  unsigned long long event_ = 0;
  int nInputTracks_ = 0;
  int nRawCandidates_ = 0;
  int nAfterKpiCA_ = 0;
  int nAfterD0pisCA_ = 0;
  int nAfterDstIso_ = 0;
  int nAfterDstPt_ = 0;
  int nAfterSoftPiDR_ = 0;
  int nAfterD0Fit_ = 0;
  int nAfterDstFit_ = 0;
  int nAfterUnique_ = 0;
  int nOutputPairs_ = 0;

  std::vector<float> rawD0Mass_;
  std::vector<float> rawD0MassKK_;
  std::vector<float> rawD0MassKpiConj_;
  std::vector<float> rawDstMass_;
  std::vector<float> rawDeltaM_;
  std::vector<float> dstPt_;
  std::vector<float> dstConePtSum_;
  std::vector<float> dstIso_;
  std::vector<float> pisDRD0Raw_;
  std::vector<float> cheapCAScore_;
  std::vector<float> deltaDxyBSKpi_;
  std::vector<float> deltaDzBSKpi_;
  std::vector<int> isSameChargeD0_;
  std::vector<int> isOppositeChargeD0_;
  std::vector<float> cheapCAKpiDxy_;
  std::vector<float> cheapCAKpiDz_;
  std::vector<float> cheapCAKpi3D_;
  std::vector<float> cheapCAD0PisDxy_;
  std::vector<float> cheapCAD0PisDz_;
  std::vector<float> cheapCAD0Pis3D_;
  std::vector<int> KCharge_;
  std::vector<int> piCharge_;
  std::vector<int> pisCharge_;
  std::vector<float> KPt_;
  std::vector<float> KEta_;
  std::vector<float> KPhi_;
  std::vector<float> KE_;
  std::vector<int> KPdgId_;
  std::vector<float> KNormChi2_;
  std::vector<float> KValidHits_;
  std::vector<float> KPixelHits_;
  std::vector<float> KLostInnerHits_;
  std::vector<float> KDxyBS_;
  std::vector<float> KDxyErrBS_;
  std::vector<float> KDxySigBS_;
  std::vector<float> KDxyErrPV_;
  std::vector<float> KDxyPV_;
  std::vector<float> KDzPV_;
  std::vector<float> KDxySigPV_;
  std::vector<float> KDzSigPV_;
  std::vector<float> piPt_;
  std::vector<float> piEta_;
  std::vector<float> piPhi_;
  std::vector<float> piE_;
  std::vector<int> piPdgId_;
  std::vector<float> piNormChi2_;
  std::vector<float> piValidHits_;
  std::vector<float> piPixelHits_;
  std::vector<float> piLostInnerHits_;
  std::vector<float> piDxyBS_;
  std::vector<float> piDxyErrBS_;
  std::vector<float> piDxySigBS_;
  std::vector<float> piDxyErrPV_;
  std::vector<float> piDxyPV_;
  std::vector<float> piDzPV_;
  std::vector<float> piDxySigPV_;
  std::vector<float> piDzSigPV_;
  std::vector<float> pisPt_;
  std::vector<float> pisEta_;
  std::vector<float> pisPhi_;
  std::vector<float> pisE_;
  std::vector<int> pisPdgId_;
  std::vector<float> pisNormChi2_;
  std::vector<float> pisValidHits_;
  std::vector<float> pisPixelHits_;
  std::vector<float> pisLostInnerHits_;
  std::vector<float> pisDxyBS_;
  std::vector<float> pisDxyErrBS_;
  std::vector<float> pisDxySigBS_;
  std::vector<float> pisDxyErrPV_;
  std::vector<float> pisDxyPV_;
  std::vector<float> pisDzPV_;
  std::vector<float> pisDxySigPV_;
  std::vector<float> pisDzSigPV_;
  std::vector<float> d0Pt_;
  std::vector<float> d0Eta_;
  std::vector<float> d0Phi_;
  std::vector<float> d0FitVertexProb_;
  std::vector<float> d0FitCos3D_;
  std::vector<float> d0FitCos2D_;
  std::vector<float> d0FitDecayLength3D_;
  std::vector<float> d0FitDecayLength3DErr_;
  std::vector<float> d0FitDecayLength3DSig_;
  std::vector<float> d0FitMass_;
  std::vector<float> d0FitChi2_;
  std::vector<float> d0FitDof_;
  std::vector<float> d0FitPval_;
  std::vector<float> d0FitDxyBS_;
  std::vector<float> d0FitDxySigBS_;
  std::vector<float> d0FitTransverseDecayLength_;
  std::vector<float> d0FitTransverseDecayLengthSig_;
  std::vector<float> dstFitMass_;
  std::vector<float> dstFitDeltaM_;
  std::vector<float> dstFitChi2_;
  std::vector<float> dstFitDof_;
  std::vector<float> dstFitPval_;
  std::vector<float> dstFitCos3D_;
  std::vector<float> dstFitCos2D_;
  std::vector<float> dstFitDecayLength3D_;
  std::vector<float> dstFitDecayLength3DErr_;
  std::vector<float> dstFitDecayLength3DSig_;
  std::vector<float> dstFitTransverseDecayLength_;
  std::vector<float> dstFitTransverseDecayLengthSig_;

  unsigned long long eventsProcessed_ = 0;
  unsigned long long eventsWithRawCandidates_ = 0;
  unsigned long long eventsWithCascadeCandidates_ = 0;
  unsigned long long eventsWithOutputPairs_ = 0;
  unsigned long long totalRawCandidates_ = 0;
  unsigned long long totalOutputPairs_ = 0;
};

DstToD0PiCascadeRefitProducer::DstToD0PiCascadeRefitProducer(const edm::ParameterSet& iConfig)
    : trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      beamSpotToken_(consumes<reco::BeamSpot>(
          iConfig.existsAs<edm::InputTag>("beamSpot") ? iConfig.getParameter<edm::InputTag>("beamSpot")
                                                       : edm::InputTag("offlineBeamSpot"))),
      kaonMass_(iConfig.getParameter<double>("kaonMass")),
      pionMass_(iConfig.getParameter<double>("pionMass")),
      softPionMass_(iConfig.getParameter<double>("softPionMass")),
      kaonMassErr_(iConfig.existsAs<double>("kaonMassErr") ? iConfig.getParameter<double>("kaonMassErr") : 1.e-6),
      pionMassErr_(iConfig.existsAs<double>("pionMassErr") ? iConfig.getParameter<double>("pionMassErr") : 1.e-6),
      hadPtMin_(iConfig.getParameter<double>("hadPtMin")),
      softPionPtMin_(iConfig.getParameter<double>("softPionPtMin")),
      maxTrackDzPV_(iConfig.existsAs<double>("maxTrackDzPV") ? iConfig.getParameter<double>("maxTrackDzPV") : -1.),
      maxKpiDeltaR_(iConfig.existsAs<double>("maxKpiDeltaR") ? iConfig.getParameter<double>("maxKpiDeltaR") : -1.),
      minD0Mass_(iConfig.getParameter<double>("minD0Mass")),
      maxD0Mass_(iConfig.getParameter<double>("maxD0Mass")),
      minDstMass_(iConfig.getParameter<double>("minDstMass")),
      maxDstMass_(iConfig.getParameter<double>("maxDstMass")),
      minDeltaMass_(iConfig.getParameter<double>("minDeltaMass")),
      maxDeltaMass_(iConfig.getParameter<double>("maxDeltaMass")),
      minD0FitPval_(iConfig.existsAs<double>("minD0FitPval") ? iConfig.getParameter<double>("minD0FitPval") : -1.),
      minD0FitMass_(iConfig.existsAs<double>("minD0FitMass") ? iConfig.getParameter<double>("minD0FitMass") : -1.),
      maxD0FitMass_(iConfig.existsAs<double>("maxD0FitMass") ? iConfig.getParameter<double>("maxD0FitMass") : -1.),
      minDstFitPval_(iConfig.existsAs<double>("minDstFitPval") ? iConfig.getParameter<double>("minDstFitPval") : -1.),
      minDstFitMass_(iConfig.existsAs<double>("minDstFitMass") ? iConfig.getParameter<double>("minDstFitMass") : -1.),
      maxDstFitMass_(iConfig.existsAs<double>("maxDstFitMass") ? iConfig.getParameter<double>("maxDstFitMass") : -1.),
      minDstFitDeltaMass_(
          iConfig.existsAs<double>("minDstFitDeltaMass") ? iConfig.getParameter<double>("minDstFitDeltaMass") : -1.),
      maxDstFitDeltaMass_(
          iConfig.existsAs<double>("maxDstFitDeltaMass") ? iConfig.getParameter<double>("maxDstFitDeltaMass") : -1.),
      maxCheapCAKpi3D_(iConfig.getParameter<double>("maxCheapCAKpi3D")),
      maxCheapCAD0Pis3D_(iConfig.getParameter<double>("maxCheapCAD0Pis3D")),
      minDstIso_(iConfig.getParameter<double>("minDstIso")),
      minDstPt_(iConfig.getParameter<double>("minDstPt")),
      maxSoftPiDRD0_(iConfig.getParameter<double>("maxSoftPiDRD0")),
      dstIsoCone_(iConfig.getParameter<double>("dstIsoCone")),
      applyChargeFilter_(iConfig.getParameter<bool>("applyChargeFilter")),
      useUnsignedCharge_(iConfig.getParameter<bool>("useUnsignedCharge")),
      applyUniqueSelection_(
          iConfig.existsAs<bool>("applyUniqueSelection") ? iConfig.getParameter<bool>("applyUniqueSelection") : true),
      writeTree_(iConfig.existsAs<bool>("writeTree") ? iConfig.getParameter<bool>("writeTree") : true),
      targetCharge_(iConfig.getParameter<int>("charge")) {
  usesResource("TFileService");
  produces<reco::TrackCollection>();
  produces<reco::TrackCollection>("dstarTracks");
  if (writeTree_) {
    bookTree();
  }
}

DstToD0PiCascadeRefitProducer::~DstToD0PiCascadeRefitProducer() {
  std::cout << "DstToD0PiCascadeRefitProducer summary"
            << " eventsProcessed=" << eventsProcessed_
            << " eventsWithRawCandidates=" << eventsWithRawCandidates_
            << " eventsWithCascadeCandidates=" << eventsWithCascadeCandidates_
            << " eventsWithOutputPairs=" << eventsWithOutputPairs_
            << " totalRawCandidates=" << totalRawCandidates_
            << " totalOutputPairs=" << totalOutputPairs_ << std::endl;
}

void DstToD0PiCascadeRefitProducer::bookTree() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("DstCascade", "Dst cascade and transient-track fit diagnostics");
  tree_->Branch("run", &run_);
  tree_->Branch("lumi", &lumi_);
  tree_->Branch("event", &event_);
  tree_->Branch("nInputTracks", &nInputTracks_);
  tree_->Branch("nRawCandidates", &nRawCandidates_);
  tree_->Branch("nAfterKpiCA", &nAfterKpiCA_);
  tree_->Branch("nAfterD0pisCA", &nAfterD0pisCA_);
  tree_->Branch("nAfterDstIso", &nAfterDstIso_);
  tree_->Branch("nAfterDstPt", &nAfterDstPt_);
  tree_->Branch("nAfterSoftPiDR", &nAfterSoftPiDR_);
  tree_->Branch("nAfterD0Fit", &nAfterD0Fit_);
  tree_->Branch("nAfterDstFit", &nAfterDstFit_);
  tree_->Branch("nAfterUnique", &nAfterUnique_);
  tree_->Branch("nOutputPairs", &nOutputPairs_);
  tree_->Branch("massb_piK", &rawD0Mass_);
  tree_->Branch("massb_KK", &rawD0MassKK_);
  tree_->Branch("massb_Kpi_conj", &rawD0MassKpiConj_);
  tree_->Branch("raw_D0_mass", &rawD0Mass_);
  tree_->Branch("raw_Dst_mass", &rawDstMass_);
  tree_->Branch("raw_DeltaM", &rawDeltaM_);
  tree_->Branch("Dst_pt", &dstPt_);
  tree_->Branch("Dst_sumPtCone04", &dstConePtSum_);
  tree_->Branch("Dst_iso", &dstIso_);
  tree_->Branch("preIso_Dst_sumPtCone04", &dstConePtSum_);
  tree_->Branch("preIso_Dst_iso", &dstIso_);
  tree_->Branch("pis_dR_D0raw", &pisDRD0Raw_);
  tree_->Branch("pis_dR_D0", &pisDRD0Raw_);
  tree_->Branch("isSameChargeD0", &isSameChargeD0_);
  tree_->Branch("isOppositeChargeD0", &isOppositeChargeD0_);
  tree_->Branch("prefit_isSameChargeD0", &isSameChargeD0_);
  tree_->Branch("prefit_isOppositeChargeD0", &isOppositeChargeD0_);
  tree_->Branch("prefit_raw_D0_mass", &rawD0Mass_);
  tree_->Branch("prefit_raw_Dst_mass", &rawDstMass_);
  tree_->Branch("prefit_raw_DeltaM", &rawDeltaM_);
  tree_->Branch("prefit_Dst_pt", &dstPt_);
  tree_->Branch("prefit_pis_dR_D0raw", &pisDRD0Raw_);
  tree_->Branch("prefit_deltaDxyBS_Kpi", &deltaDxyBSKpi_);
  tree_->Branch("prefit_deltaDzBS_Kpi", &deltaDzBSKpi_);
  tree_->Branch("prefit_Dst_sumPtCone04", &dstConePtSum_);
  tree_->Branch("prefit_Dst_iso", &dstIso_);
  tree_->Branch("cheapCA_Kpi_dxy", &cheapCAKpiDxy_);
  tree_->Branch("cheapCA_Kpi_zAtCA", &cheapCAKpiDz_);
  tree_->Branch("cheapCA_Kpi_3d", &cheapCAKpi3D_);
  tree_->Branch("cheapCA_D0pis_dxy", &cheapCAD0PisDxy_);
  tree_->Branch("cheapCA_D0pis_zAtCA", &cheapCAD0PisDz_);
  tree_->Branch("cheapCA_D0pis_3d", &cheapCAD0Pis3D_);
  tree_->Branch("cheapCA_score", &cheapCAScore_);
  tree_->Branch("prefit_cheapCA_Kpi_dxy", &cheapCAKpiDxy_);
  tree_->Branch("prefit_cheapCA_Kpi_zAtCA", &cheapCAKpiDz_);
  tree_->Branch("prefit_cheapCA_Kpi_3d", &cheapCAKpi3D_);
  tree_->Branch("prefit_cheapCA_D0pis_dxy", &cheapCAD0PisDxy_);
  tree_->Branch("prefit_cheapCA_D0pis_zAtCA", &cheapCAD0PisDz_);
  tree_->Branch("prefit_cheapCA_D0pis_3d", &cheapCAD0Pis3D_);
  tree_->Branch("prefit_cheapCA_score", &cheapCAScore_);
  tree_->Branch("K_charge", &KCharge_);
  tree_->Branch("pi_charge", &piCharge_);
  tree_->Branch("pis_charge", &pisCharge_);
  tree_->Branch("K_pt", &KPt_);
  tree_->Branch("K_E", &KE_);
  tree_->Branch("K_eta", &KEta_);
  tree_->Branch("K_phi", &KPhi_);
  tree_->Branch("K_pdgId", &KPdgId_);
  tree_->Branch("K_norm_chi2", &KNormChi2_);
  tree_->Branch("K_normChi2", &KNormChi2_);
  tree_->Branch("K_N_valid_hits", &KValidHits_);
  tree_->Branch("K_validHits", &KValidHits_);
  tree_->Branch("K_pixelHits", &KPixelHits_);
  tree_->Branch("K_lostInnerHits", &KLostInnerHits_);
  tree_->Branch("K_dxy_BS", &KDxyBS_);
  tree_->Branch("K_dxyErr_BS", &KDxyErrBS_);
  tree_->Branch("K_sigdxy_BS", &KDxySigBS_);
  tree_->Branch("K_dxy_PV", &KDxyPV_);
  tree_->Branch("K_dxyErr_PV", &KDxyErrPV_);
  tree_->Branch("K_sigdxy_PV", &KDxySigPV_);
  tree_->Branch("K_dxyPV", &KDxyPV_);
  tree_->Branch("K_dzPV", &KDzPV_);
  tree_->Branch("K_dxySigPV", &KDxySigPV_);
  tree_->Branch("K_dzSigPV", &KDzSigPV_);
  tree_->Branch("pi_pt", &piPt_);
  tree_->Branch("pi_E", &piE_);
  tree_->Branch("pi_eta", &piEta_);
  tree_->Branch("pi_phi", &piPhi_);
  tree_->Branch("pi_pdgId", &piPdgId_);
  tree_->Branch("pi_norm_chi2", &piNormChi2_);
  tree_->Branch("pi_normChi2", &piNormChi2_);
  tree_->Branch("pi_N_valid_hits", &piValidHits_);
  tree_->Branch("pi_validHits", &piValidHits_);
  tree_->Branch("pi_pixelHits", &piPixelHits_);
  tree_->Branch("pi_lostInnerHits", &piLostInnerHits_);
  tree_->Branch("pi_dxy_BS", &piDxyBS_);
  tree_->Branch("pi_dxyErr_BS", &piDxyErrBS_);
  tree_->Branch("pi_sigdxy_BS", &piDxySigBS_);
  tree_->Branch("pi_dxy_PV", &piDxyPV_);
  tree_->Branch("pi_dxyErr_PV", &piDxyErrPV_);
  tree_->Branch("pi_sigdxy_PV", &piDxySigPV_);
  tree_->Branch("pi_dxyPV", &piDxyPV_);
  tree_->Branch("pi_dzPV", &piDzPV_);
  tree_->Branch("pi_dxySigPV", &piDxySigPV_);
  tree_->Branch("pi_dzSigPV", &piDzSigPV_);
  tree_->Branch("pis_pt", &pisPt_);
  tree_->Branch("pis_E", &pisE_);
  tree_->Branch("pis_eta", &pisEta_);
  tree_->Branch("pis_phi", &pisPhi_);
  tree_->Branch("pis_pdgId", &pisPdgId_);
  tree_->Branch("pis_norm_chi2", &pisNormChi2_);
  tree_->Branch("pis_normChi2", &pisNormChi2_);
  tree_->Branch("pis_N_valid_hits", &pisValidHits_);
  tree_->Branch("pis_validHits", &pisValidHits_);
  tree_->Branch("pis_pixelHits", &pisPixelHits_);
  tree_->Branch("pis_lostInnerHits", &pisLostInnerHits_);
  tree_->Branch("pis_dxy_BS", &pisDxyBS_);
  tree_->Branch("pis_dxyErr_BS", &pisDxyErrBS_);
  tree_->Branch("pis_sigdxy_BS", &pisDxySigBS_);
  tree_->Branch("pis_dxy_PV", &pisDxyPV_);
  tree_->Branch("pis_dxyErr_PV", &pisDxyErrPV_);
  tree_->Branch("pis_sigdxy_PV", &pisDxySigPV_);
  tree_->Branch("pis_dxyPV", &pisDxyPV_);
  tree_->Branch("pis_dzPV", &pisDzPV_);
  tree_->Branch("pis_dxySigPV", &pisDxySigPV_);
  tree_->Branch("pis_dzSigPV", &pisDzSigPV_);
  tree_->Branch("D0_pt", &d0Pt_);
  tree_->Branch("D0_eta", &d0Eta_);
  tree_->Branch("D0_phi", &d0Phi_);
  tree_->Branch("D0_fit_vertexProb", &d0FitVertexProb_);
  tree_->Branch("D0_fit_cos3D", &d0FitCos3D_);
  tree_->Branch("D0_fit_cos2D", &d0FitCos2D_);
  tree_->Branch("D0_fit_decayLength3D", &d0FitDecayLength3D_);
  tree_->Branch("D0_fit_decayLength3DErr", &d0FitDecayLength3DErr_);
  tree_->Branch("D0_fit_decayLength3DSig", &d0FitDecayLength3DSig_);
  tree_->Branch("cos_D0_PV", &d0FitCos3D_);
  tree_->Branch("cosT_D0_PV", &d0FitCos2D_);
  tree_->Branch("d_vtxD0_PV", &d0FitDecayLength3D_);
  tree_->Branch("sigd_vtxD0_PV", &d0FitDecayLength3DSig_);
  tree_->Branch("dxy_vtxD0_BS", &d0FitDxyBS_);
  tree_->Branch("sigdxy_vtxD0_BS", &d0FitDxySigBS_);
  tree_->Branch("chi2_piK", &d0FitChi2_);
  tree_->Branch("dof_piK", &d0FitDof_);
  tree_->Branch("pval_piK", &d0FitPval_);
  tree_->Branch("isValid_piK", &isOppositeChargeD0_);
  tree_->Branch("isGood_piK", &isOppositeChargeD0_);
  tree_->Branch("mass_piK", &d0FitMass_);
  tree_->Branch("D0_fit_mass", &d0FitMass_);
  tree_->Branch("D0_fit_chi2", &d0FitChi2_);
  tree_->Branch("D0_fit_dof", &d0FitDof_);
  tree_->Branch("D0_fit_pval", &d0FitPval_);
  tree_->Branch("Dst_fit_mass", &dstFitMass_);
  tree_->Branch("mass_D0pis", &dstFitMass_);
  tree_->Branch("Dst_fit_deltaM", &dstFitDeltaM_);
  tree_->Branch("deltaM_D0pis_piK", &dstFitDeltaM_);
  tree_->Branch("chi2_D0pis", &dstFitChi2_);
  tree_->Branch("dof_D0pis", &dstFitDof_);
  tree_->Branch("pval_D0pis", &dstFitPval_);
  tree_->Branch("isValid_D0pis", &isOppositeChargeD0_);
  tree_->Branch("isGood_D0pis", &isOppositeChargeD0_);
  tree_->Branch("cos_Dst_PV", &dstFitCos3D_);
  tree_->Branch("cosT_Dst_PV", &dstFitCos2D_);
  tree_->Branch("d_vtxDst_PV", &dstFitDecayLength3D_);
  tree_->Branch("sigd_vtxDst_PV", &dstFitDecayLength3DSig_);
  tree_->Branch("dxy_vtxDst_PV", &dstFitTransverseDecayLength_);
  tree_->Branch("sigdxy_vtxDst_PV", &dstFitTransverseDecayLengthSig_);
  tree_->Branch("Dst_fit_chi2", &dstFitChi2_);
  tree_->Branch("Dst_fit_dof", &dstFitDof_);
  tree_->Branch("Dst_fit_pval", &dstFitPval_);
}

void DstToD0PiCascadeRefitProducer::clearVectors() {
  rawD0Mass_.clear();
  rawD0MassKK_.clear();
  rawD0MassKpiConj_.clear();
  rawDstMass_.clear();
  rawDeltaM_.clear();
  dstPt_.clear();
  dstConePtSum_.clear();
  dstIso_.clear();
  pisDRD0Raw_.clear();
  cheapCAScore_.clear();
  deltaDxyBSKpi_.clear();
  deltaDzBSKpi_.clear();
  isSameChargeD0_.clear();
  isOppositeChargeD0_.clear();
  cheapCAKpiDxy_.clear();
  cheapCAKpiDz_.clear();
  cheapCAKpi3D_.clear();
  cheapCAD0PisDxy_.clear();
  cheapCAD0PisDz_.clear();
  cheapCAD0Pis3D_.clear();
  KCharge_.clear();
  piCharge_.clear();
  pisCharge_.clear();
  KPt_.clear();
  KEta_.clear();
  KPhi_.clear();
  KE_.clear();
  KPdgId_.clear();
  KNormChi2_.clear();
  KValidHits_.clear();
  KPixelHits_.clear();
  KLostInnerHits_.clear();
  KDxyBS_.clear();
  KDxyErrBS_.clear();
  KDxySigBS_.clear();
  KDxyErrPV_.clear();
  KDxyPV_.clear();
  KDzPV_.clear();
  KDxySigPV_.clear();
  KDzSigPV_.clear();
  piPt_.clear();
  piEta_.clear();
  piPhi_.clear();
  piE_.clear();
  piPdgId_.clear();
  piNormChi2_.clear();
  piValidHits_.clear();
  piPixelHits_.clear();
  piLostInnerHits_.clear();
  piDxyBS_.clear();
  piDxyErrBS_.clear();
  piDxySigBS_.clear();
  piDxyErrPV_.clear();
  piDxyPV_.clear();
  piDzPV_.clear();
  piDxySigPV_.clear();
  piDzSigPV_.clear();
  pisPt_.clear();
  pisEta_.clear();
  pisPhi_.clear();
  pisE_.clear();
  pisPdgId_.clear();
  pisNormChi2_.clear();
  pisValidHits_.clear();
  pisPixelHits_.clear();
  pisLostInnerHits_.clear();
  pisDxyBS_.clear();
  pisDxyErrBS_.clear();
  pisDxySigBS_.clear();
  pisDxyErrPV_.clear();
  pisDxyPV_.clear();
  pisDzPV_.clear();
  pisDxySigPV_.clear();
  pisDzSigPV_.clear();
  d0Pt_.clear();
  d0Eta_.clear();
  d0Phi_.clear();
  d0FitVertexProb_.clear();
  d0FitCos3D_.clear();
  d0FitCos2D_.clear();
  d0FitDecayLength3D_.clear();
  d0FitDecayLength3DErr_.clear();
  d0FitDecayLength3DSig_.clear();
  d0FitMass_.clear();
  d0FitChi2_.clear();
  d0FitDof_.clear();
  d0FitPval_.clear();
  d0FitDxyBS_.clear();
  d0FitDxySigBS_.clear();
  d0FitTransverseDecayLength_.clear();
  d0FitTransverseDecayLengthSig_.clear();
  dstFitMass_.clear();
  dstFitDeltaM_.clear();
  dstFitChi2_.clear();
  dstFitDof_.clear();
  dstFitPval_.clear();
  dstFitCos3D_.clear();
  dstFitCos2D_.clear();
  dstFitDecayLength3D_.clear();
  dstFitDecayLength3DErr_.clear();
  dstFitDecayLength3DSig_.clear();
  dstFitTransverseDecayLength_.clear();
  dstFitTransverseDecayLengthSig_.clear();
}

bool DstToD0PiCascadeRefitProducer::checkCharge(const reco::Track& kaon,
                                                const reco::Track& pion,
                                                const reco::Track& softPion) const {
  if (!applyChargeFilter_)
    return true;

  if (kaon.charge() + pion.charge() != 0)
    return false;
  if (softPion.charge() != pion.charge())
    return false;

  int totalCharge = kaon.charge() + pion.charge() + softPion.charge();
  if (useUnsignedCharge_)
    totalCharge = std::abs(totalCharge);
  return totalCharge == targetCharge_;
}

FitQuality DstToD0PiCascadeRefitProducer::fitQuality(RefCountedKinematicTree tree) const {
  FitQuality quality;
  if (!tree || !tree->isValid())
    return quality;
  tree->movePointerToTheTop();
  RefCountedKinematicVertex vertex = tree->currentDecayVertex();
  if (!vertex || !vertex->vertexIsValid())
    return quality;
  quality.isValid = true;
  quality.chi2 = vertex->chiSquared();
  quality.dof = vertex->degreesOfFreedom();
  quality.pval = TMath::Prob(quality.chi2, quality.dof);
  return quality;
}

RefCountedKinematicTree DstToD0PiCascadeRefitProducer::fitD0(const edm::EventSetup& iSetup,
                                                             const reco::Track& kaon,
                                                             const reco::Track& pion) const {
  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);

  KinematicParticleFactoryFromTransientTrack particleFactory;
  std::vector<RefCountedKinematicParticle> parts;
  float chi = 0.f;
  float ndf = 0.f;
  float kaonMassErr = kaonMassErr_;
  float pionMassErr = pionMassErr_;
  reco::TransientTrack kaonTT = ttBuilder->build(kaon);
  reco::TransientTrack pionTT = ttBuilder->build(pion);
  parts.push_back(particleFactory.particle(kaonTT, kaonMass_, chi, ndf, kaonMassErr));
  parts.push_back(particleFactory.particle(pionTT, pionMass_, chi, ndf, pionMassErr));

  KinematicParticleVertexFitter fitter;
  return fitter.fit(parts);
}

RefCountedKinematicTree DstToD0PiCascadeRefitProducer::fitDst(const edm::EventSetup& iSetup,
                                                              const reco::Track& softPion,
                                                              const RefCountedKinematicParticle& d0) const {
  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);

  KinematicParticleFactoryFromTransientTrack particleFactory;
  std::vector<RefCountedKinematicParticle> parts;
  float chi = 0.f;
  float ndf = 0.f;
  float softPionMassErr = pionMassErr_;
  float d0Mass = d0->currentState().mass();
  float d0MassErr = std::sqrt(std::abs(d0->currentState().kinematicParametersError().matrix()(6, 6)));
  reco::TransientTrack softPionTT = ttBuilder->build(softPion);
  reco::TransientTrack d0TT = d0->refittedTransientTrack();
  parts.push_back(particleFactory.particle(softPionTT, softPionMass_, chi, ndf, softPionMassErr));
  parts.push_back(particleFactory.particle(d0TT, d0Mass, chi, ndf, d0MassErr));

  KinematicParticleVertexFitter fitter;
  return fitter.fit(parts);
}

void DstToD0PiCascadeRefitProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  ++eventsProcessed_;
  clearVectors();
  run_ = iEvent.run();
  lumi_ = iEvent.luminosityBlock();
  event_ = iEvent.id().event();
  nInputTracks_ = 0;
  nRawCandidates_ = 0;
  nAfterKpiCA_ = 0;
  nAfterD0pisCA_ = 0;
  nAfterDstIso_ = 0;
  nAfterDstPt_ = 0;
  nAfterSoftPiDR_ = 0;
  nAfterD0Fit_ = 0;
  nAfterDstFit_ = 0;
  nAfterUnique_ = 0;
  nOutputPairs_ = 0;

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackToken_, tracks);
  edm::Handle<reco::VertexCollection> vertices;
  try {
    iEvent.getByToken(vertexToken_, vertices);
  } catch (const cms::Exception&) {
  }
  edm::Handle<reco::BeamSpot> beamSpot;
  try {
    iEvent.getByToken(beamSpotToken_, beamSpot);
  } catch (const cms::Exception&) {
  }
  const reco::Vertex* primaryVertex = nullptr;
  if (vertices.isValid()) {
    for (const auto& vertex : *vertices) {
      if (!vertex.isFake() && vertex.isValid() && vertex.ndof() > 4.) {
        primaryVertex = &vertex;
        break;
      }
    }
  }

  auto selectedPairs = std::make_unique<reco::TrackCollection>();
  auto selectedTriplets = std::make_unique<reco::TrackCollection>();
  std::vector<SelectedCandidate> selectedCandidates;

  if (tracks.isValid()) {
    nInputTracks_ = tracks->size();
    for (unsigned int iK = 0; iK < tracks->size(); ++iK) {
      const reco::Track& kaon = tracks->at(iK);
      if (kaon.pt() < hadPtMin_ || kaon.isLooper())
        continue;
      if (maxTrackDzPV_ >= 0. && primaryVertex != nullptr &&
          std::abs(kaon.dz(primaryVertex->position())) > maxTrackDzPV_)
        continue;

      for (unsigned int iPi = 0; iPi < tracks->size(); ++iPi) {
        if (iPi == iK)
          continue;
        const reco::Track& pion = tracks->at(iPi);
        if (pion.pt() < hadPtMin_ || pion.isLooper())
          continue;
        if (maxTrackDzPV_ >= 0. && primaryVertex != nullptr &&
            std::abs(pion.dz(primaryVertex->position())) > maxTrackDzPV_)
          continue;
        if (maxKpiDeltaR_ >= 0. && deltaR(kaon.eta(), kaon.phi(), pion.eta(), pion.phi()) > maxKpiDeltaR_)
          continue;
        if (pion.charge() != -kaon.charge())
          continue;

        const TLorentzVector rawD0 = p4(kaon, kaonMass_) + p4(pion, pionMass_);
        if (rawD0.M() < minD0Mass_ || rawD0.M() > maxD0Mass_)
          continue;
        const TLorentzVector rawD0KK = p4(kaon, kaonMass_) + p4(pion, kaonMass_);
        const TLorentzVector rawD0KpiConj = p4(kaon, pionMass_) + p4(pion, kaonMass_);

        const CheapCAResult kpiCA = cheapClosestApproachXY(makeCheapLine(kaon), makeCheapLine(pion));
        for (unsigned int iPis = 0; iPis < tracks->size(); ++iPis) {
          if (iPis == iK || iPis == iPi)
            continue;
          const reco::Track& softPion = tracks->at(iPis);
          if (softPion.pt() < softPionPtMin_ || softPion.isLooper())
            continue;
          if (maxTrackDzPV_ >= 0. && primaryVertex != nullptr &&
              std::abs(softPion.dz(primaryVertex->position())) > maxTrackDzPV_)
            continue;
          if (!checkCharge(kaon, pion, softPion))
            continue;

          const TLorentzVector rawDst = rawD0 + p4(softPion, softPionMass_);
          const double rawDeltaM = rawDst.M() - rawD0.M();
          if (rawDst.M() < minDstMass_ || rawDst.M() > maxDstMass_)
            continue;
          if (rawDeltaM < minDeltaMass_ || rawDeltaM > maxDeltaMass_)
            continue;

          ++nRawCandidates_;
          const CheapLine rawD0Line = makeCheapLine(kpiCA.midX, kpiCA.midY, kpiCA.midZ, rawD0.Px(), rawD0.Py(), rawD0.Pz());
          const CheapCAResult d0pisCA = cheapClosestApproachXY(rawD0Line, makeCheapLine(softPion));
          const double cheapCAScore = kpiCA.dist3d / 0.05 + d0pisCA.dist3d / 0.10 +
                                      std::abs(rawDeltaM - kNominalDeltaM) / 0.006 +
                                      std::abs(rawD0.M() - 1.86483) / 0.20;
          if (maxCheapCAKpi3D_ >= 0. && kpiCA.dist3d > maxCheapCAKpi3D_)
            continue;
          ++nAfterKpiCA_;
          if (maxCheapCAD0Pis3D_ >= 0. && d0pisCA.dist3d > maxCheapCAD0Pis3D_)
            continue;
          ++nAfterD0pisCA_;

          double dstConePtSum = 0.;
          for (const auto& isoTrack : *tracks) {
            if (deltaR(isoTrack.eta(), isoTrack.phi(), rawDst.Eta(), rawDst.Phi()) < dstIsoCone_)
              dstConePtSum += isoTrack.pt();
          }
          const double dstIso = dstConePtSum > 0. ? rawDst.Pt() / dstConePtSum : -1.;
          if (minDstIso_ >= 0. && dstIso < minDstIso_)
            continue;
          ++nAfterDstIso_;
          if (minDstPt_ >= 0. && rawDst.Pt() < minDstPt_)
            continue;
          ++nAfterDstPt_;
          const double pisDRD0Raw = deltaR(softPion.eta(), softPion.phi(), rawD0.Eta(), rawD0.Phi());
          if (maxSoftPiDRD0_ >= 0. && pisDRD0Raw > maxSoftPiDRD0_)
            continue;
          ++nAfterSoftPiDR_;

          RefCountedKinematicTree d0Tree;
          try {
            d0Tree = fitD0(iSetup, kaon, pion);
          } catch (const std::exception&) {
            continue;
          }
          const FitQuality d0Quality = fitQuality(d0Tree);
          if (!d0Quality.isValid)
            continue;
          ++nAfterD0Fit_;
          d0Tree->movePointerToTheTop();
          RefCountedKinematicVertex d0Vertex = d0Tree->currentDecayVertex();
          RefCountedKinematicParticle d0 = d0Tree->currentParticle();
          const double d0FitMass = d0->currentState().mass();
          if (minD0FitPval_ >= 0. && d0Quality.pval < minD0FitPval_)
            continue;
          if (minD0FitMass_ >= 0. && d0FitMass < minD0FitMass_)
            continue;
          if (maxD0FitMass_ >= 0. && d0FitMass > maxD0FitMass_)
            continue;

          RefCountedKinematicTree dstTree;
          try {
            dstTree = fitDst(iSetup, softPion, d0);
          } catch (const std::exception&) {
            continue;
          }
          const FitQuality dstQuality = fitQuality(dstTree);
          if (!dstQuality.isValid)
            continue;
          ++nAfterDstFit_;
          dstTree->movePointerToTheTop();
          RefCountedKinematicVertex dstVertex = dstTree->currentDecayVertex();
          RefCountedKinematicParticle dst = dstTree->currentParticle();
          const double dstFitMass = dst->currentState().mass();
          const double dstFitDeltaM = dstFitMass - d0FitMass;
          if (minDstFitPval_ >= 0. && dstQuality.pval < minDstFitPval_)
            continue;
          if (minDstFitMass_ >= 0. && dstFitMass < minDstFitMass_)
            continue;
          if (maxDstFitMass_ >= 0. && dstFitMass > maxDstFitMass_)
            continue;
          if (minDstFitDeltaMass_ >= 0. && dstFitDeltaM < minDstFitDeltaMass_)
            continue;
          if (maxDstFitDeltaMass_ >= 0. && dstFitDeltaM > maxDstFitDeltaMass_)
            continue;

          SelectedCandidate candidate;
          candidate.kaonIndex = iK;
          candidate.pionIndex = iPi;
          candidate.softPionIndex = iPis;
          candidate.kaon = kaon;
          candidate.pion = pion;
          candidate.rawD0Mass = rawD0.M();
          candidate.rawD0MassKK = rawD0KK.M();
          candidate.rawD0MassKpiConj = rawD0KpiConj.M();
          candidate.rawDstMass = rawDst.M();
          candidate.rawDeltaM = rawDeltaM;
          candidate.dstPt = rawDst.Pt();
          candidate.dstConePtSum = dstConePtSum;
          candidate.dstIso = dstIso;
          candidate.pisDRD0Raw = pisDRD0Raw;
          candidate.cheapCAScore = cheapCAScore;
          candidate.deltaDxyBSKpi =
              beamSpot.isValid() ? std::abs(kaon.dxy(beamSpot->position()) - pion.dxy(beamSpot->position())) : -99.;
          candidate.deltaDzBSKpi =
              beamSpot.isValid() ? std::abs(kaon.dz(beamSpot->position()) - pion.dz(beamSpot->position())) : -99.;
          candidate.isSameChargeD0 = (kaon.charge() == pion.charge());
          candidate.isOppositeChargeD0 = (kaon.charge() == -pion.charge());
          candidate.cheapCAKpiDxy = kpiCA.dxy;
          candidate.cheapCAKpiDz = kpiCA.dz;
          candidate.cheapCAKpi3D = kpiCA.dist3d;
          candidate.cheapCAD0PisDxy = d0pisCA.dxy;
          candidate.cheapCAD0PisDz = d0pisCA.dz;
          candidate.cheapCAD0Pis3D = d0pisCA.dist3d;
          candidate.KCharge = kaon.charge();
          candidate.piCharge = pion.charge();
          candidate.pisCharge = softPion.charge();
          candidate.KPt = kaon.pt();
          candidate.KEta = kaon.eta();
          candidate.KPhi = kaon.phi();
          candidate.KE = p4(kaon, kaonMass_).E();
          candidate.KPdgId = kaon.charge() > 0 ? 321 : -321;
          candidate.KNormChi2 = kaon.normalizedChi2();
          candidate.KValidHits = kaon.numberOfValidHits();
          candidate.KPixelHits = kaon.hitPattern().numberOfValidPixelHits();
          candidate.KLostInnerHits = kaon.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
          candidate.KDxyBS = beamSpot.isValid() ? kaon.dxy(beamSpot->position()) : -99.;
          candidate.KDxyErrBS = kaon.dxyError();
          candidate.KDxySigBS = safeSignificance(candidate.KDxyBS, candidate.KDxyErrBS);
          candidate.KDxyPV = primaryVertex ? kaon.dxy(primaryVertex->position()) : -99.;
          candidate.KDxyErrPV = kaon.dxyError();
          candidate.KDzPV = primaryVertex ? kaon.dz(primaryVertex->position()) : -99.;
          candidate.KDxySigPV = primaryVertex ? safeTrackDxySignificance(kaon, *primaryVertex) : -1.;
          candidate.KDzSigPV = primaryVertex ? safeTrackDzSignificance(kaon, *primaryVertex) : -1.;
          candidate.piPt = pion.pt();
          candidate.piEta = pion.eta();
          candidate.piPhi = pion.phi();
          candidate.piE = p4(pion, pionMass_).E();
          candidate.piPdgId = pion.charge() > 0 ? 211 : -211;
          candidate.piNormChi2 = pion.normalizedChi2();
          candidate.piValidHits = pion.numberOfValidHits();
          candidate.piPixelHits = pion.hitPattern().numberOfValidPixelHits();
          candidate.piLostInnerHits = pion.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
          candidate.piDxyBS = beamSpot.isValid() ? pion.dxy(beamSpot->position()) : -99.;
          candidate.piDxyErrBS = pion.dxyError();
          candidate.piDxySigBS = safeSignificance(candidate.piDxyBS, candidate.piDxyErrBS);
          candidate.piDxyPV = primaryVertex ? pion.dxy(primaryVertex->position()) : -99.;
          candidate.piDxyErrPV = pion.dxyError();
          candidate.piDzPV = primaryVertex ? pion.dz(primaryVertex->position()) : -99.;
          candidate.piDxySigPV = primaryVertex ? safeTrackDxySignificance(pion, *primaryVertex) : -1.;
          candidate.piDzSigPV = primaryVertex ? safeTrackDzSignificance(pion, *primaryVertex) : -1.;
          candidate.pisPt = softPion.pt();
          candidate.pisEta = softPion.eta();
          candidate.pisPhi = softPion.phi();
          candidate.pisE = p4(softPion, softPionMass_).E();
          candidate.pisPdgId = softPion.charge() > 0 ? 211 : -211;
          candidate.pisNormChi2 = softPion.normalizedChi2();
          candidate.pisValidHits = softPion.numberOfValidHits();
          candidate.pisPixelHits = softPion.hitPattern().numberOfValidPixelHits();
          candidate.pisLostInnerHits = softPion.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
          candidate.pisDxyBS = beamSpot.isValid() ? softPion.dxy(beamSpot->position()) : -99.;
          candidate.pisDxyErrBS = softPion.dxyError();
          candidate.pisDxySigBS = safeSignificance(candidate.pisDxyBS, candidate.pisDxyErrBS);
          candidate.pisDxyPV = primaryVertex ? softPion.dxy(primaryVertex->position()) : -99.;
          candidate.pisDxyErrPV = softPion.dxyError();
          candidate.pisDzPV = primaryVertex ? softPion.dz(primaryVertex->position()) : -99.;
          candidate.pisDxySigPV = primaryVertex ? safeTrackDxySignificance(softPion, *primaryVertex) : -1.;
          candidate.pisDzSigPV = primaryVertex ? safeTrackDzSignificance(softPion, *primaryVertex) : -1.;
          candidate.d0Pt = rawD0.Pt();
          candidate.d0Eta = rawD0.Eta();
          candidate.d0Phi = rawD0.Phi();
          candidate.d0FitVertexProb = d0Quality.pval;
          candidate.d0FitCos3D = primaryVertex ? pointingCos3D(*primaryVertex, d0Vertex, rawD0) : -99.;
          candidate.d0FitCos2D = primaryVertex ? pointingCos2D(*primaryVertex, d0Vertex, rawD0) : -99.;
          const DistanceWithError d0Distance = primaryVertex ? vertexDistance3D(*primaryVertex, d0Vertex) : DistanceWithError();
          candidate.d0FitDecayLength3D = d0Distance.value;
          candidate.d0FitDecayLength3DErr = d0Distance.error;
          candidate.d0FitDecayLength3DSig = d0Distance.significance;
          candidate.d0FitMass = d0FitMass;
          candidate.d0FitChi2 = d0Quality.chi2;
          candidate.d0FitDof = d0Quality.dof;
          candidate.d0FitPval = d0Quality.pval;
          candidate.d0FitDxyBS = beamSpot.isValid() ? vertexTransverseDistanceFromBeamSpot(*beamSpot, d0Vertex) : -99.;
          candidate.d0FitDxySigBS = -99.;
          candidate.d0FitTransverseDecayLength = primaryVertex ? vertexTransverseDistance(*primaryVertex, d0Vertex) : -99.;
          candidate.d0FitTransverseDecayLengthSig = -99.;
          candidate.dstFitMass = dstFitMass;
          candidate.dstFitDeltaM = dstFitDeltaM;
          candidate.dstFitChi2 = dstQuality.chi2;
          candidate.dstFitDof = dstQuality.dof;
          candidate.dstFitPval = dstQuality.pval;
          candidate.dstFitCos3D = primaryVertex ? pointingCos3D(*primaryVertex, dstVertex, rawDst) : -99.;
          candidate.dstFitCos2D = primaryVertex ? pointingCos2D(*primaryVertex, dstVertex, rawDst) : -99.;
          const DistanceWithError dstDistance = primaryVertex ? vertexDistance3D(*primaryVertex, dstVertex) : DistanceWithError();
          candidate.dstFitDecayLength3D = dstDistance.value;
          candidate.dstFitDecayLength3DErr = dstDistance.error;
          candidate.dstFitDecayLength3DSig = dstDistance.significance;
          candidate.dstFitTransverseDecayLength = primaryVertex ? vertexTransverseDistance(*primaryVertex, dstVertex) : -99.;
          candidate.dstFitTransverseDecayLengthSig = -99.;
          selectedCandidates.push_back(candidate);
        }
      }
    }
  }

  std::sort(selectedCandidates.begin(),
            selectedCandidates.end(),
            [](const SelectedCandidate& lhs, const SelectedCandidate& rhs) {
              const double lhsDelta = std::abs(lhs.rawDeltaM - kNominalDeltaM);
              const double rhsDelta = std::abs(rhs.rawDeltaM - kNominalDeltaM);
              if (lhsDelta != rhsDelta)
                return lhsDelta < rhsDelta;
              return lhs.dstPt > rhs.dstPt;
            });

  std::unordered_set<unsigned int> usedTracks;
  std::unordered_set<unsigned long long> usedKpiPairs;
  for (const auto& candidate : selectedCandidates) {
    if (applyUniqueSelection_) {
      const unsigned long long kpiKey = pairKey(candidate.kaonIndex, candidate.pionIndex);
      if (usedKpiPairs.count(kpiKey) != 0)
        continue;
      if (usedTracks.count(candidate.kaonIndex) != 0 || usedTracks.count(candidate.pionIndex) != 0 ||
          usedTracks.count(candidate.softPionIndex) != 0)
        continue;

      usedKpiPairs.insert(kpiKey);
      usedTracks.insert(candidate.kaonIndex);
      usedTracks.insert(candidate.pionIndex);
      usedTracks.insert(candidate.softPionIndex);
    }
    ++nAfterUnique_;

    rawD0Mass_.push_back(candidate.rawD0Mass);
    rawD0MassKK_.push_back(candidate.rawD0MassKK);
    rawD0MassKpiConj_.push_back(candidate.rawD0MassKpiConj);
    rawDstMass_.push_back(candidate.rawDstMass);
    rawDeltaM_.push_back(candidate.rawDeltaM);
    dstPt_.push_back(candidate.dstPt);
    dstConePtSum_.push_back(candidate.dstConePtSum);
    dstIso_.push_back(candidate.dstIso);
    pisDRD0Raw_.push_back(candidate.pisDRD0Raw);
    cheapCAScore_.push_back(candidate.cheapCAScore);
    deltaDxyBSKpi_.push_back(candidate.deltaDxyBSKpi);
    deltaDzBSKpi_.push_back(candidate.deltaDzBSKpi);
    isSameChargeD0_.push_back(candidate.isSameChargeD0);
    isOppositeChargeD0_.push_back(candidate.isOppositeChargeD0);
    cheapCAKpiDxy_.push_back(candidate.cheapCAKpiDxy);
    cheapCAKpiDz_.push_back(candidate.cheapCAKpiDz);
    cheapCAKpi3D_.push_back(candidate.cheapCAKpi3D);
    cheapCAD0PisDxy_.push_back(candidate.cheapCAD0PisDxy);
    cheapCAD0PisDz_.push_back(candidate.cheapCAD0PisDz);
    cheapCAD0Pis3D_.push_back(candidate.cheapCAD0Pis3D);
    KCharge_.push_back(candidate.KCharge);
    piCharge_.push_back(candidate.piCharge);
    pisCharge_.push_back(candidate.pisCharge);
    KPt_.push_back(candidate.KPt);
    KEta_.push_back(candidate.KEta);
    KPhi_.push_back(candidate.KPhi);
    KE_.push_back(candidate.KE);
    KPdgId_.push_back(candidate.KPdgId);
    KNormChi2_.push_back(candidate.KNormChi2);
    KValidHits_.push_back(candidate.KValidHits);
    KPixelHits_.push_back(candidate.KPixelHits);
    KLostInnerHits_.push_back(candidate.KLostInnerHits);
    KDxyBS_.push_back(candidate.KDxyBS);
    KDxyErrBS_.push_back(candidate.KDxyErrBS);
    KDxySigBS_.push_back(candidate.KDxySigBS);
    KDxyErrPV_.push_back(candidate.KDxyErrPV);
    KDxyPV_.push_back(candidate.KDxyPV);
    KDzPV_.push_back(candidate.KDzPV);
    KDxySigPV_.push_back(candidate.KDxySigPV);
    KDzSigPV_.push_back(candidate.KDzSigPV);
    piPt_.push_back(candidate.piPt);
    piEta_.push_back(candidate.piEta);
    piPhi_.push_back(candidate.piPhi);
    piE_.push_back(candidate.piE);
    piPdgId_.push_back(candidate.piPdgId);
    piNormChi2_.push_back(candidate.piNormChi2);
    piValidHits_.push_back(candidate.piValidHits);
    piPixelHits_.push_back(candidate.piPixelHits);
    piLostInnerHits_.push_back(candidate.piLostInnerHits);
    piDxyBS_.push_back(candidate.piDxyBS);
    piDxyErrBS_.push_back(candidate.piDxyErrBS);
    piDxySigBS_.push_back(candidate.piDxySigBS);
    piDxyErrPV_.push_back(candidate.piDxyErrPV);
    piDxyPV_.push_back(candidate.piDxyPV);
    piDzPV_.push_back(candidate.piDzPV);
    piDxySigPV_.push_back(candidate.piDxySigPV);
    piDzSigPV_.push_back(candidate.piDzSigPV);
    pisPt_.push_back(candidate.pisPt);
    pisEta_.push_back(candidate.pisEta);
    pisPhi_.push_back(candidate.pisPhi);
    pisE_.push_back(candidate.pisE);
    pisPdgId_.push_back(candidate.pisPdgId);
    pisNormChi2_.push_back(candidate.pisNormChi2);
    pisValidHits_.push_back(candidate.pisValidHits);
    pisPixelHits_.push_back(candidate.pisPixelHits);
    pisLostInnerHits_.push_back(candidate.pisLostInnerHits);
    pisDxyBS_.push_back(candidate.pisDxyBS);
    pisDxyErrBS_.push_back(candidate.pisDxyErrBS);
    pisDxySigBS_.push_back(candidate.pisDxySigBS);
    pisDxyErrPV_.push_back(candidate.pisDxyErrPV);
    pisDxyPV_.push_back(candidate.pisDxyPV);
    pisDzPV_.push_back(candidate.pisDzPV);
    pisDxySigPV_.push_back(candidate.pisDxySigPV);
    pisDzSigPV_.push_back(candidate.pisDzSigPV);
    d0Pt_.push_back(candidate.d0Pt);
    d0Eta_.push_back(candidate.d0Eta);
    d0Phi_.push_back(candidate.d0Phi);
    d0FitVertexProb_.push_back(candidate.d0FitVertexProb);
    d0FitCos3D_.push_back(candidate.d0FitCos3D);
    d0FitCos2D_.push_back(candidate.d0FitCos2D);
    d0FitDecayLength3D_.push_back(candidate.d0FitDecayLength3D);
    d0FitDecayLength3DErr_.push_back(candidate.d0FitDecayLength3DErr);
    d0FitDecayLength3DSig_.push_back(candidate.d0FitDecayLength3DSig);
    d0FitMass_.push_back(candidate.d0FitMass);
    d0FitChi2_.push_back(candidate.d0FitChi2);
    d0FitDof_.push_back(candidate.d0FitDof);
    d0FitPval_.push_back(candidate.d0FitPval);
    d0FitDxyBS_.push_back(candidate.d0FitDxyBS);
    d0FitDxySigBS_.push_back(candidate.d0FitDxySigBS);
    d0FitTransverseDecayLength_.push_back(candidate.d0FitTransverseDecayLength);
    d0FitTransverseDecayLengthSig_.push_back(candidate.d0FitTransverseDecayLengthSig);
    dstFitMass_.push_back(candidate.dstFitMass);
    dstFitDeltaM_.push_back(candidate.dstFitDeltaM);
    dstFitChi2_.push_back(candidate.dstFitChi2);
    dstFitDof_.push_back(candidate.dstFitDof);
    dstFitPval_.push_back(candidate.dstFitPval);
    dstFitCos3D_.push_back(candidate.dstFitCos3D);
    dstFitCos2D_.push_back(candidate.dstFitCos2D);
    dstFitDecayLength3D_.push_back(candidate.dstFitDecayLength3D);
    dstFitDecayLength3DErr_.push_back(candidate.dstFitDecayLength3DErr);
    dstFitDecayLength3DSig_.push_back(candidate.dstFitDecayLength3DSig);
    dstFitTransverseDecayLength_.push_back(candidate.dstFitTransverseDecayLength);
    dstFitTransverseDecayLengthSig_.push_back(candidate.dstFitTransverseDecayLengthSig);

    selectedPairs->push_back(candidate.kaon);
    selectedPairs->push_back(candidate.pion);
    selectedTriplets->push_back(candidate.kaon);
    selectedTriplets->push_back(candidate.pion);
    selectedTriplets->push_back(tracks->at(candidate.softPionIndex));
  }

  nOutputPairs_ = selectedPairs->size() / 2;
  totalRawCandidates_ += nRawCandidates_;
  totalOutputPairs_ += nOutputPairs_;
  if (nRawCandidates_ > 0)
    ++eventsWithRawCandidates_;
  if (nAfterSoftPiDR_ > 0)
    ++eventsWithCascadeCandidates_;
  if (nOutputPairs_ > 0)
    ++eventsWithOutputPairs_;

  if (writeTree_) {
    tree_->Fill();
  }
  iEvent.put(std::move(selectedPairs));
  iEvent.put(std::move(selectedTriplets), "dstarTracks");
}

DEFINE_FWK_MODULE(DstToD0PiCascadeRefitProducer);

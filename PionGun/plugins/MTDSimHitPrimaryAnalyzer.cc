#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "SimG4CMS/Forward/interface/MtdHitCategory.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include "TH1D.h"

#include <algorithm>
#include <unordered_map>
#include <vector>
#include <tuple>
#include <functional>

class MTDSimHitPrimaryAnalyzer : public edm::one::EDAnalyzer<> {
public:
  explicit MTDSimHitPrimaryAnalyzer(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

private:
  struct TrackAccum {
    unsigned nETL = 0;   // counts only PDG-matched primary-associated ETL hits
    double tMin = 1e99;
    double tMax = -1e99;
    bool hasHit = false;
  };

  struct DiskKey {
    int zside = 0;  // ±1
    int disc  = 0;  // ETL disk index
    bool operator==(const DiskKey& other) const { return zside == other.zside && disc == other.disc; }
  };

  struct DiskKeyHash {
    std::size_t operator()(const DiskKey& k) const {
      // simple hash combining zside/disc
      return std::hash<int>{}(k.zside * 100 + k.disc);
    }
  };

  struct FaceSeen {
    bool seenPosFace = false;  // local entry z > 0
    bool seenNegFace = false;  // local entry z < 0
    unsigned nHits = 0;
    unsigned nFromBack = 0;
  };

  edm::EDGetTokenT<std::vector<PSimHit>> tokETL_;
  edm::EDGetTokenT<std::vector<PSimHit>> tokBTL_;  // kept for interface stability
  edm::EDGetTokenT<edm::SimTrackContainer> tokSimTrk_;
  edm::EDGetTokenT<edm::SimVertexContainer> tokSimVtx_;

  TH1D* h_nETLPrimaryPdgMatched_;
  TH1D* h_timeSpreadETLPrimaryPdgMatched_;
  TH1D* h_backFracPerEventPrimaryETL_PdgMatched_;

  // New histograms for "both faces of same disk by a single track"
  TH1D* h_nDiskBothFacesPerTrack_;
  TH1D* h_hasAnyDiskBothFacesPerEvent_;

  // Global counters across all events (primary-associated + PDG exact match only)
  unsigned long long totalETLHitsPrimaryPdgMatchedAll_ = 0;
  unsigned long long totalETLHitsPrimaryPdgMatchedFromBack_ = 0;
  unsigned long long totalEvents_ = 0;

  // Global counters for both-face check
  unsigned long long totalTrackDiskPairsChecked_ = 0;
  unsigned long long totalTrackDiskPairsBothFaces_ = 0;
};

MTDSimHitPrimaryAnalyzer::MTDSimHitPrimaryAnalyzer(const edm::ParameterSet& iConfig) {
  tokETL_ = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("etlSimHits"));
  tokBTL_ = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("btlSimHits"));
  tokSimTrk_ = consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("simTracks"));
  tokSimVtx_ = consumes<edm::SimVertexContainer>(iConfig.getParameter<edm::InputTag>("simVertices"));

  edm::Service<TFileService> fs;
  h_nETLPrimaryPdgMatched_ = fs->make<TH1D>(
      "h_nETLPrimaryPdgMatched",
      "Primary-associated ETL PSimHit multiplicity per primary SimTrack (PDG exact match);N_{ETL}^{primary-associated,PDG};Tracks",
      20, -0.5, 19.5);

  h_timeSpreadETLPrimaryPdgMatched_ = fs->make<TH1D>(
      "h_timeSpreadETLPrimaryPdgMatched",
      "Primary-associated ETL time spread (PDG exact match);#Delta t=t_{max}-t_{min} [ns];Tracks with N_{ETL}#geq1",
      20, 0.0, 0.2);

  h_backFracPerEventPrimaryETL_PdgMatched_ = fs->make<TH1D>(
      "h_backFracPerEventPrimaryETL_PdgMatched",
      "Event-level ETL from-back fraction (primary-associated, PDG exact match);N_{ETL,fromBack}^{assoc,PDG}/N_{ETL,all}^{assoc,PDG};Events",
      110, 0.0, 1.1);

  h_nDiskBothFacesPerTrack_ = fs->make<TH1D>(
      "h_nDiskBothFacesPerTrack",
      "Number of ETL disks with hits on both faces for a single primary track (PDG exact match);N_{disk}^{both faces};Tracks",
      10, -0.5, 9.5);

  h_hasAnyDiskBothFacesPerEvent_ = fs->make<TH1D>(
      "h_hasAnyDiskBothFacesPerEvent",
      "Whether event contains any (track,disk) with hits on both faces (PDG exact match);Flag (0/1);Events",
      2, -0.5, 1.5);
}

void MTDSimHitPrimaryAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  ++totalEvents_;

  edm::Handle<edm::SimTrackContainer> hTrk;
  edm::Handle<std::vector<PSimHit>> hETL;
  edm::Handle<std::vector<PSimHit>> hBTL;

  iEvent.getByToken(tokSimTrk_, hTrk);
  iEvent.getByToken(tokETL_, hETL);
  iEvent.getByToken(tokBTL_, hBTL);  // intentionally unused
  (void)hBTL;

  // SimTrack ID -> (PDG, isPrimaryApprox)
  std::unordered_map<unsigned, std::pair<int, bool>> trkInfo;
  trkInfo.reserve(hTrk->size());

  // primary SimTrack ID -> accumulators (PDG-matched ETL hits only)
  std::unordered_map<unsigned, TrackAccum> acc;
  acc.reserve(hTrk->size());

  for (auto const& trk : *hTrk) {
    const bool isPrimary = (trk.genpartIndex() >= 0);  // particle gun 场景通常够用
    const unsigned tid = trk.trackId();                // SimTrack ID
    trkInfo.emplace(tid, std::make_pair(trk.type(), isPrimary));
    if (isPrimary) {
      acc.emplace(tid, TrackAccum{});  // include nETL=0 tracks
    }
  }

  // Per primary track: per (zside,disc) face occupancy using local entry z sign
  std::unordered_map<unsigned, std::unordered_map<DiskKey, FaceSeen, DiskKeyHash>> trackDiskFaces;
  trackDiskFaces.reserve(hTrk->size());

  // Event-level counters (ONLY PDG exact-match sample)
  unsigned long long nETLHitsPrimaryPdgMatchedAll = 0;
  unsigned long long nETLHitsPrimaryPdgMatchedFromBack = 0;

  // Fill ETL hit info (primary-associated + PDG exact-match only)
  for (auto const& hit : *hETL) {
    const unsigned tidOrig = hit.originalTrackId();  // decoded base ID for matching SimTrack
    const unsigned tidOff  = hit.offsetTrackId();    // offset category (e.g. k_idETLfromBack)

    auto it = trkInfo.find(tidOrig);
    if (it == trkInfo.end()) continue;

    const int primaryPdg = it->second.first;
    const bool isPrimary = it->second.second;
    if (!isPrimary) continue;  // primary-associated only

    const int hitPdg = hit.particleType();
    if (hitPdg != primaryPdg) continue;  // only keep PDG exact match version

    auto ait = acc.find(tidOrig);
    if (ait == acc.end()) continue;  // defensive
    auto& a = ait->second;

    // Track-level accumulators (PDG-matched sample only)
    a.nETL++;
    a.hasHit = true;

    const double t = hit.timeOfFlight();
    a.tMin = std::min(a.tMin, t);
    a.tMax = std::max(a.tMax, t);

    // Event-level counters
    ++nETLHitsPrimaryPdgMatchedAll;
    if (tidOff == MtdHitCategory::k_idETLfromBack) {
      ++nETLHitsPrimaryPdgMatchedFromBack;
    }

    // ---- NEW: identify disk and which face the hit entered from ----
    const ETLDetId etlId(hit.detUnitId());
    DiskKey dk;
    dk.zside = etlId.zside();
    dk.disc  = etlId.nDisc();  // if compile error: try etlId.disc()

    // Use local entry z sign as face label.
    // For grazing cases entry z can be ~0; fall back to exit z if needed.
    const float zin  = hit.entryPoint().z();
    const float zout = hit.exitPoint().z();

    auto& fs = trackDiskFaces[tidOrig][dk];
    if (zin > 0.f) {
      fs.seenPosFace = true;
    } else if (zin < 0.f) {
      fs.seenNegFace = true;
    } else {
      if (zout > 0.f) {
        fs.seenPosFace = true;
      } else if (zout < 0.f) {
        fs.seenNegFace = true;
      }
    }

    fs.nHits++;
    if (tidOff == MtdHitCategory::k_idETLfromBack) fs.nFromBack++;
  }

  // Accumulate global counters
  totalETLHitsPrimaryPdgMatchedAll_ += nETLHitsPrimaryPdgMatchedAll;
  totalETLHitsPrimaryPdgMatchedFromBack_ += nETLHitsPrimaryPdgMatchedFromBack;

  edm::LogPrint("MTDSimHitPrimaryAnalyzer")
      << "Run " << iEvent.id().run() << " Event " << iEvent.id().event()
      << " : primary-associated ETL simhit multiplicity + time spread (PDG exact match only)"
      << " (trackId: pdgId, nETL, dT[ns])";

  if (acc.empty()) {
    edm::LogPrint("MTDSimHitPrimaryAnalyzer")
        << "  (no primary SimTracks found by current definition: genpartIndex()>=0)";
  } else {
    // First, both-faces summary per primary track
    bool eventHasAnyTrackDiskBothFaces = false;

    for (auto const& kv : acc) {
      const unsigned tid = kv.first;
      unsigned nDiskBothFacesThisTrack = 0;

      auto itTD = trackDiskFaces.find(tid);
      if (itTD != trackDiskFaces.end()) {
        for (auto const& kv2 : itTD->second) {
          const DiskKey& dk = kv2.first;
          const FaceSeen& fs = kv2.second;

          ++totalTrackDiskPairsChecked_;
          const bool bothFaces = (fs.seenPosFace && fs.seenNegFace);

          if (bothFaces) {
            ++nDiskBothFacesThisTrack;
            ++totalTrackDiskPairsBothFaces_;
            eventHasAnyTrackDiskBothFaces = true;

            edm::LogPrint("MTDSimHitPrimaryAnalyzer")
                << "  [BOTH_FACES] trackId=" << tid
                << " zside=" << dk.zside
                << " disc=" << dk.disc
                << " nHits=" << fs.nHits
                << " nFromBack=" << fs.nFromBack;
          }
        }
      }

      h_nDiskBothFacesPerTrack_->Fill(nDiskBothFacesThisTrack);
    }

    h_hasAnyDiskBothFacesPerEvent_->Fill(eventHasAnyTrackDiskBothFaces ? 1 : 0);

    // Original per-track printout and hist filling
    for (auto const& kv : acc) {
      const unsigned tid = kv.first;  // SimTrack ID
      const auto& a = kv.second;
      const int pdg = trkInfo[tid].first;

      // Includes nETL=0 primary tracks in the PDG-matched sample
      h_nETLPrimaryPdgMatched_->Fill(a.nETL);

      if (a.hasHit) {
        const double dT = a.tMax - a.tMin;
        h_timeSpreadETLPrimaryPdgMatched_->Fill(dT);

        edm::LogPrint("MTDSimHitPrimaryAnalyzer")
            << "  trackId=" << tid << " pdgId=" << pdg
            << " nETL=" << a.nETL
            << " dT=" << dT;
      } else {
        edm::LogPrint("MTDSimHitPrimaryAnalyzer")
            << "  trackId=" << tid << " pdgId=" << pdg
            << " nETL=0 (no ETL simhits in PDG-matched sample)";
      }
    }
  }

  // Event-level ETL from-back fraction (PDG exact-match sample only)
  if (nETLHitsPrimaryPdgMatchedAll > 0) {
    const double fracBackPdg =
        static_cast<double>(nETLHitsPrimaryPdgMatchedFromBack) /
        static_cast<double>(nETLHitsPrimaryPdgMatchedAll);

    h_backFracPerEventPrimaryETL_PdgMatched_->Fill(fracBackPdg);

    edm::LogPrint("MTDSimHitPrimaryAnalyzer")
        << "  [ETL primary-associated + PDG exact-match hits] total="
        << nETLHitsPrimaryPdgMatchedAll
        << " fromBack(offset=" << MtdHitCategory::k_idETLfromBack << ")="
        << nETLHitsPrimaryPdgMatchedFromBack
        << " fraction=" << fracBackPdg;
  } else {
    edm::LogPrint("MTDSimHitPrimaryAnalyzer")
        << "  [ETL primary-associated + PDG exact-match hits] total=0 fromBack=0 fraction=nan";
  }
}

void MTDSimHitPrimaryAnalyzer::endJob() {
  edm::LogPrint("MTDSimHitPrimaryAnalyzer")
      << "================ Final summary (all events) ================";

  edm::LogPrint("MTDSimHitPrimaryAnalyzer")
      << "Processed events = " << totalEvents_;

  edm::LogPrint("MTDSimHitPrimaryAnalyzer")
      << "[Primary-associated + PDG exact-match] total ETL hits = "
      << totalETLHitsPrimaryPdgMatchedAll_;

  edm::LogPrint("MTDSimHitPrimaryAnalyzer")
      << "[Primary-associated + PDG exact-match] ETL from-back hits (offset="
      << MtdHitCategory::k_idETLfromBack << ") = "
      << totalETLHitsPrimaryPdgMatchedFromBack_;

  if (totalETLHitsPrimaryPdgMatchedAll_ > 0) {
    const double frac =
        static_cast<double>(totalETLHitsPrimaryPdgMatchedFromBack_) /
        static_cast<double>(totalETLHitsPrimaryPdgMatchedAll_);
    edm::LogPrint("MTDSimHitPrimaryAnalyzer")
        << "[Primary-associated + PDG exact-match] global ETL from-back fraction = " << frac;
  } else {
    edm::LogPrint("MTDSimHitPrimaryAnalyzer")
        << "[Primary-associated + PDG exact-match] global ETL from-back fraction = nan (no ETL hits)";
  }

  edm::LogPrint("MTDSimHitPrimaryAnalyzer")
      << "[Both-faces check] checked (track,disk) pairs = " << totalTrackDiskPairsChecked_;

  edm::LogPrint("MTDSimHitPrimaryAnalyzer")
      << "[Both-faces check] (track,disk) pairs with hits on both faces = "
      << totalTrackDiskPairsBothFaces_;

  if (totalTrackDiskPairsChecked_ > 0) {
    const double fracBoth =
        static_cast<double>(totalTrackDiskPairsBothFaces_) /
        static_cast<double>(totalTrackDiskPairsChecked_);
    edm::LogPrint("MTDSimHitPrimaryAnalyzer")
        << "[Both-faces check] fraction of (track,disk) pairs with both-face hits = " << fracBoth;
  } else {
    edm::LogPrint("MTDSimHitPrimaryAnalyzer")
        << "[Both-faces check] fraction of (track,disk) pairs with both-face hits = nan";
  }

  edm::LogPrint("MTDSimHitPrimaryAnalyzer")
      << "============================================================";
}

DEFINE_FWK_MODULE(MTDSimHitPrimaryAnalyzer);
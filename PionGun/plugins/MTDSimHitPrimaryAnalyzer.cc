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

#include "TH1D.h"

#include <algorithm>
#include <unordered_map>
#include <vector>

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

  edm::EDGetTokenT<std::vector<PSimHit>> tokETL_;
  edm::EDGetTokenT<std::vector<PSimHit>> tokBTL_;  // kept for interface stability
  edm::EDGetTokenT<edm::SimTrackContainer> tokSimTrk_;
  edm::EDGetTokenT<edm::SimVertexContainer> tokSimVtx_;

  TH1D* h_nETLPrimaryPdgMatched_;
  TH1D* h_timeSpreadETLPrimaryPdgMatched_;
  TH1D* h_backFracPerEventPrimaryETL_PdgMatched_;

  // Global counters across all events (primary-associated + PDG exact match only)
  unsigned long long totalETLHitsPrimaryPdgMatchedAll_ = 0;
  unsigned long long totalETLHitsPrimaryPdgMatchedFromBack_ = 0;
  unsigned long long totalEvents_ = 0;
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
      200, 0.0, 10.0);

  h_backFracPerEventPrimaryETL_PdgMatched_ = fs->make<TH1D>(
      "h_backFracPerEventPrimaryETL_PdgMatched",
      "Event-level ETL from-back fraction (primary-associated, PDG exact match);N_{ETL,fromBack}^{assoc,PDG}/N_{ETL,all}^{assoc,PDG};Events",
      110, 0.0, 1.1);
}

void MTDSimHitPrimaryAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  ++totalEvents_;

  edm::Handle<edm::SimTrackContainer> hTrk;
  edm::Handle<std::vector<PSimHit>> hETL;
  edm::Handle<std::vector<PSimHit>> hBTL;

  iEvent.getByToken(tokSimTrk_, hTrk);
  iEvent.getByToken(tokETL_, hETL);
  iEvent.getByToken(tokBTL_, hBTL);  // intentionally unused

  // SimTrack ID -> (PDG, isPrimaryApprox)
  std::unordered_map<unsigned, std::pair<int, bool>> trkInfo;
  trkInfo.reserve(hTrk->size());

  // primary SimTrack ID -> accumulators (PDG-matched ETL hits only)
  std::unordered_map<unsigned, TrackAccum> acc;
  acc.reserve(hTrk->size());

  for (auto const& trk : *hTrk) {
    const bool isPrimary = (trk.genpartIndex() >= 0);  // particle gun场景通常够用
    const unsigned tid = trk.trackId();                // SimTrack ID
    trkInfo.emplace(tid, std::make_pair(trk.type(), isPrimary));
    if (isPrimary) {
      acc.emplace(tid, TrackAccum{});  // include nETL=0 tracks
    }
  }

  // Event-level counters (ONLY PDG exact-match sample)
  unsigned long long nETLHitsPrimaryPdgMatchedAll = 0;
  unsigned long long nETLHitsPrimaryPdgMatchedFromBack = 0;

  // Fill ETL hit info (primary-associated + PDG exact-match only)
  for (auto const& hit : *hETL) {
    const unsigned tidOrig = hit.originalTrackId();  // decoded base ID for matching SimTrack
    const unsigned tidOff  = hit.offsetTrackId();    // offset category (e.g. 4 = ETL from-back)

    auto it = trkInfo.find(tidOrig);
    if (it == trkInfo.end()) continue;

    const int primaryPdg = it->second.first;
    const bool isPrimary = it->second.second;
    if (!isPrimary) continue;  // primary-associated only

    const int hitPdg = hit.particleType();
    if (hitPdg != primaryPdg) continue;  // <-- only keep PDG exact match version

    auto ait = acc.find(tidOrig);
    if (ait == acc.end()) continue;  // defensive
    auto& a = ait->second;

    // Track-level accumulators (now for PDG-matched sample only)
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
      << "============================================================";
}

DEFINE_FWK_MODULE(MTDSimHitPrimaryAnalyzer);
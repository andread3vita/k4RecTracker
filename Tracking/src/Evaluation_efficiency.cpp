#include "Evaluation_efficiency.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <vector>

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHit3D.h"
#include "extension/MutableTrackerHit3D.h"

double pt_thr   = 0.100;
double cos_thr  = 0.99;
double hits_thr = 4;

double pur_thr = 0.66;
double eff_thr = 0.05;

std::vector<double> computeCosThetaAndPt(double px, double py, double pz) {


  double P        = std::sqrt(px * px + py * py + pz * pz);
  double cosTheta = std::abs(pz / P);

  double pt = std::sqrt(px * px + py * py);
  
  return {cosTheta, pt};

}


DECLARE_COMPONENT(Evaluation_efficiency)

Evaluation_efficiency::Evaluation_efficiency(const std::string& aName, ISvcLocator* aSvcLoc): Gaudi::Algorithm(aName, aSvcLoc) {

    declareProperty("inputTracks", m_input_tracks, "Input track collection name");
    declareProperty("inputMCparticles", m_input_MCparticles, "Input MC particles");

    declareProperty("inputHits_CDC_sim", m_input_hits_CDC_sim, "Input CDC sim hit collection name");
    declareProperty("inputHits_VTXIB_sim", m_input_hits_VTXIB_sim, "Input VTXIB sim hit collection name");
    declareProperty("inputHits_VTXD_sim", m_input_hits_VTXD_sim, "Input VTXD sim hit collection name");
    declareProperty("inputHits_VTXOB_sim", m_input_hits_VTXOB_sim, "Input VTXOB sim hit collection name");

    // declareProperty("outputTracksEff", m_output_trackEff, "output track efficiencies");
    // declareProperty("outputMCPur", m_output_MCpur, "output MC purities");

    declareProperty("outputTrackingEff", m_output_tracking, "output tracking efficiency");
}

Evaluation_efficiency::~Evaluation_efficiency() {}

StatusCode Evaluation_efficiency::initialize() {
    if (Gaudi::Algorithm::initialize().isFailure()) {
        return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
}

StatusCode Evaluation_efficiency::execute(const EventContext&) const {

    /// Get tracks and MC_particles
    const edm4hep::MCParticleCollection* MC_particles = m_input_MCparticles.get();
    const extension::TrackCollection*          tracks       = m_input_tracks.get();

    const edm4hep::SimTrackerHitCollection* inputHits_VTXD_sim  = m_input_hits_VTXD_sim.get();
    const edm4hep::SimTrackerHitCollection* inputHits_VTXIB_sim = m_input_hits_VTXOB_sim.get();
    const edm4hep::SimTrackerHitCollection* inputHits_VTXOB_sim = m_input_hits_VTXOB_sim.get();
    const edm4hep::SimTrackerHitCollection* inputHits_CDC_sim   = m_input_hits_CDC_sim.get();

    std::vector<double> costheta_mc = {};
    std::vector<double> pt_mc = {};
    std::vector<int>  ListMC;
    for (const auto& MC_par : *MC_particles) {
      // define a lorentz vector and compute theta and pt
      double         px     = MC_par.getMomentum().x;
      double         py     = MC_par.getMomentum().y;
      double         pz     = MC_par.getMomentum().z;
      std::vector<double> info = computeCosThetaAndPt(px, py, pz);

      double pt       = info[1];
      double costheta = info[0];

      // extract the object ID
      auto index_MC = static_cast<int>(MC_par.getObjectID().index);

      // fill all the vectors with those information
      pt_mc.push_back(pt);
      costheta_mc.push_back(costheta);
      ListMC.push_back(index_MC);
    }

    /// compute the number of hits for each particle
    std::vector<int> particle_hits(ListMC.size(), 0);
    for (const auto& hit : *inputHits_VTXD_sim) {

        int part_index = hit.getParticle().getObjectID().index;
        particle_hits[part_index] +=1;
    }

    for (const auto& hit : *inputHits_VTXIB_sim) {
      int part_index = hit.getParticle().getObjectID().index;
      particle_hits[part_index] += 1;
    }

    for (const auto& hit : *inputHits_VTXOB_sim) {
      int part_index = hit.getParticle().getObjectID().index;
      particle_hits[part_index] += 1;
    }

    for (const auto& hit : *inputHits_CDC_sim) {
      int part_index = hit.getParticle().getObjectID().index;
      particle_hits[part_index] += 1;
    }

    size_t numPart = particle_hits.size();

    std::cout << "MC particles properties:" << std::endl;
    for (size_t i = 0; i < numPart; i++) {
      std::cout << "Particle:" << i << "\tcostheta:" << costheta_mc[i] << "\tpt:" << pt_mc[i] << "\tnumHits:" << particle_hits[i] << std::endl;
    }

    /// check which particle is reconstrutable
    std::vector<bool> isReco;
    for (size_t i = 0; i < numPart; i++) {

        bool isRec = (pt_mc[i] > pt_thr) && (costheta_mc[i] < cos_thr) && (particle_hits[i] > hits_thr);
        isReco.push_back(isRec);
    
    }

    /// Efficiency
    int total_hits = 0;
    std::vector<std::map<int, double>> eff_mc;
    std::vector<std::map<int, int>>    hits_Tracks;
    // auto output_eff = m_output_trackEff.createAndPut();

    for (const auto& track : *tracks) {

        auto hits_in_track = track.getTrackerHits();
        int  numHits       = hits_in_track.size() > 0 ? hits_in_track.size() : 1;

        std::vector<int> hit_idx;
        hit_idx.reserve(hits_in_track.size());
        for (int k = 0; k < numHits; ++k) {
          int idx = static_cast<int>(hits_in_track[k].getEDep());
          hit_idx.push_back(idx);
          total_hits += 1;
        }

        std::map<int, int> counter;
        for (int idx : hit_idx) {
            counter[idx]++;
        }

        std::map<int, double> EFF;
        std::map<int, int>    hitsPart;
        for (size_t value = 0; value < numPart; ++value) {
            
            EFF[value]      = counter[value] / static_cast<double>(numHits);
            hitsPart[value] = counter[value];
        
        }

        eff_mc.push_back(EFF);
        hits_Tracks.push_back(hitsPart);

        // output_eff->push_back(EFF);

    }

    /// Purity
    // auto output_pur = m_output_MCpur.createAndPut();
    std::vector<std::map<int, double>> pur_tracks;
    int temp_count = 0;
    for (const auto& part : *MC_particles) {
        int numHitsTracker = 0;

        for (const auto& hits_map : hits_Tracks) {
          auto it = hits_map.find(temp_count);
          if (it != hits_map.end()) {
            numHitsTracker += it->second;
          }
        }

        numHitsTracker = (numHitsTracker > 0) ? numHitsTracker : 1;

        std::map<int, double> PUR;
        for (size_t value = 0; value < hits_Tracks.size(); ++value) {
            
            auto it       = hits_Tracks[value].find(temp_count);
            int  hitCount = (it != hits_Tracks[value].end()) ? it->second : 0;
            PUR[value]    = static_cast<double>(hitCount) / numHitsTracker;
        }

        temp_count += 1;
        pur_tracks.push_back(PUR);
        // output_pur->push_back(PUR);
    }

   
    /// compute tracking efficiency
    std::cout << "\nSummary:" << std::endl;
    int numberOfMatched = 0;
    for (size_t i = 0; i < pur_tracks.size(); ++i) {
        for (size_t s = 0; s < eff_mc.size(); ++s) {
            double pur = 0.0;
            double eff = 0.0;

            // Find the purity
            auto pur_it = pur_tracks[i].find(s);
            if (pur_it != pur_tracks[i].end()) {
                pur = pur_it->second;
            }

            // Find the efficiency
            auto eff_it = eff_mc[s].find(i);
            if (eff_it != eff_mc[s].end()) {
                eff = eff_it->second;
            }

            std::cout << "Particle: " << i << " Track: " << s << "\t\tPur: " << pur << "\t\tEff: " << eff << std::endl;

            // Check if purity and efficiency are above the thresholds
            if (pur > pur_thr && eff > eff_thr && isReco[i] == true) {
                numberOfMatched++;
            }

      }
    }

    int tracking_eff_denom                     = std::accumulate(isReco.begin(), isReco.end(), 0);
    int tracking_eff_num                       = numberOfMatched;

    auto output_tracking = m_output_tracking.createAndPut();
    output_tracking->push_back(tracking_eff_num);
    output_tracking->push_back(tracking_eff_denom);

    return StatusCode::SUCCESS;
}

StatusCode Evaluation_efficiency::finalize() { return StatusCode::SUCCESS; }

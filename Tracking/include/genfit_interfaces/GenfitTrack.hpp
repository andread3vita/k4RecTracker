/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GENFIT_TRACK_H
#define GENFIT_TRACK_H

// Standard Library
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include<ranges>

// ROOT
#include <TVector3.h>

// GenFit
#include <AbsBField.h>
#include <AbsTrackRep.h>
#include <DAF.h>
#include <FieldManager.h>
#include <RKTrackRep.h>
#include <Track.h>
#include <MaterialEffects.h>

#include "GenfitPlanarMeasurement.hpp"
#include "GenfitWireMeasurement.hpp"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Fields.h"
#include "DDRec/Vector3D.h"
#include "DDSegmentation/BitFieldCoder.h"

// EDM4hep
#include "edm4hep/MutableTrack.h"
#include "edm4hep/TrackCollection.h"

// Local
#include "utils.hpp"

/** @class GenfitTrack
 *
 *  Internal helper class that bridges an EDM4hep track with the GENFIT track representation.
 *  This class is responsible for preparing, initializing, and managing the GENFIT track object,
 *  starting from an EDM4hep `edm4hep::Track`, and performing the track fit using the GENFIT library.
 *
 *  The `GenfitTrack` encapsulates the logic required to:
 *    - extract and convert the track parameters (position, momentum, charge)
 *    - create and configure the appropriate `genfit::AbsTrackRep` based on the particle hypothesis
 *    - construct the `genfit::Track` object including measurement hits
 *    - execute the GENFIT fitting procedure
 *
 *  The class maintains access to both the EDM4hep and GENFIT representations of the track, and stores
 *  relevant geometry information such as the drift chamber (DCH) description and segmentation decoder.
 *
 *
 *  Author: Andrea De Vita  
 *  Date  : 2025-11
 *
 */
namespace GenfitInterface {

    class GenfitTrack {
    public:

        GenfitTrack(const edm4hep::Track& track, const dd4hep::rec::DCH_info* dch_info, const dd4hep::DDSegmentation::BitFieldCoder* decoder, const int particle_hypothesis, std::optional<int> maxHitForLoopers = std::nullopt, std::optional<int> maxHitForLoopers_back = std::nullopt, std::optional<TVector3> initial_position = std::nullopt, std::optional<TVector3> initial_momentum = std::nullopt, std::optional<TVector3> initial_position_back = std::nullopt, std::optional<TVector3> initial_momentum_back = std::nullopt);
        ~GenfitTrack();

        void checkInitialization();
        void init(const edm4hep::Track& track_init, std::optional<int> maxHitForLoopers = std::nullopt, std::optional<int> maxHitForLoopers_back = std::nullopt, std::optional<TVector3> initial_position = std::nullopt, std::optional<TVector3> initial_momentum = std::nullopt, std::optional<TVector3> initial_position_back = std::nullopt, std::optional<TVector3> initial_momentum_back = std::nullopt);
        
        void createGenFitTrack(int dir, int debug_lvl, std::optional<double> pz_initial = std::nullopt, std::optional<double> pT_initial = std::nullopt);
        bool fit(int dir, double Beta_init, double Beta_final, double Beta_steps, double Bz, int debug_lvl);

        genfit::Track* getTrack_genfit() { return m_genfitTrack; }
        genfit::AbsTrackRep* getRep_genfit() { return m_genfitTrackRep; }

        edm4hep::MutableTrack& getTrack_edm4hep(int dir = 1) { 
            
            if (dir > 0) return m_edm4hepTrack;
            else  return m_edm4hepTrack_back;
 }

        void getTrack_init(int dir = 1) { 


            if (dir > 0) {
                std::cout << "GENFIT Initial position: (" << m_posInit.X() << ", " << m_posInit.Y() << ", " << m_posInit.Z() << ")" << std::endl;
                std::cout << "GENFIT Initial momentum: (" << m_momInit.X() << ", " << m_momInit.Y() << ", " << m_momInit.Z() << ")" << std::endl;
            }
            else
            {
                std::cout << "GENFIT Initial position (backward): (" << m_posInit_back.X() << ", " << m_posInit_back.Y() << ", " << m_posInit_back.Z() << ")" << std::endl;
                std::cout << "GENFIT Initial momentum (backward): (" << m_momInit_back.X() << ", " << m_momInit_back.Y() << ", " << m_momInit_back.Z() << ")" << std::endl;
            }
        }

    private:

        int m_particle_hypothesis;
        TVector3 m_posInit{0., 0., 0.};
        TVector3 m_momInit{0., 0., 0.};

        genfit::AbsTrackRep* m_genfitTrackRep;    
        genfit::Track* m_genfitTrack;
        
        edm4hep::MutableTrack m_edm4hepTrack;  

        const dd4hep::rec::DCH_info* m_dch_info;
        const dd4hep::DDSegmentation::BitFieldCoder* m_dc_decoder;

        TVector3 m_IP_referencePoint{0., 0., 0.};
        TVector3 m_firstHit_referencePoint;
        TVector3 m_lastHit_referencePoint;


        TVector3 m_posInit_back{0., 0., 0.};
        TVector3 m_momInit_back{0., 0., 0.};

        genfit::AbsTrackRep* m_genfitTrackRep_back;    
        genfit::Track* m_genfitTrack_back;
        
        edm4hep::MutableTrack m_edm4hepTrack_back;  

        const dd4hep::rec::DCH_info* m_dch_info_back;
        const dd4hep::DDSegmentation::BitFieldCoder* m_dc_decoder_back;

        TVector3 m_IP_referencePoint_back{0., 0., 0.};
        TVector3 m_firstHit_referencePoint_back;
        TVector3 m_lastHit_referencePoint_back;

    };

}  // namespace GenfitInterface

#endif // GENFIT_TRACK_H
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

#include "Planar_measurement.hpp"
#include "Wire_measurement.hpp"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Fields.h"
#include "DDRec/Vector3D.h"
#include "DDSegmentation/BitFieldCoder.h"

// EDM4hep
#include "extension/MutableTrack.h"
#include "extension/TrackCollection.h"

// Local
#include "utils.hpp"

/** @class GenfitTrack
 *
 *  Internal helper class that bridges an EDM4hep track with the GENFIT track representation.
 *  This class is responsible for preparing, initializing, and managing the GENFIT track object,
 *  starting from an EDM4hep `extension::Track`, and performing the track fit using the GENFIT library.
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
 *  Date  : 2025-06
 *
 */


namespace GenfitInterface {

    class GenfitTrack {
    public:
        GenfitTrack(const extension::Track& track, const dd4hep::rec::DCH_info* dch_info,const dd4hep::DDSegmentation::BitFieldCoder* decoder, const int particle_hypothesis);
        ~GenfitTrack();

        void checkInitialization();
        void init(const extension::Track& track_init);

        void createGenFitTrack(int debug_lvl);
        bool fit(double Beta_init, double Beta_final, double Beta_steps);

        genfit::Track* getTrack_genfit() { return genfitTrack_; }
        genfit::AbsTrackRep* getRep_genfit() { return genfitTrackRep_; }
        extension::MutableTrack& getTrack_edm4hep() { return edm4hepTrack_; }

    private:

        int _particle_hypothesis;
        TVector3 _posInit{0., 0., 0.};
        TVector3 _momInit{0., 0., 0.};

        genfit::AbsTrackRep* genfitTrackRep_;    
        genfit::AbsTrackRep* genfitTrackRepBack_;    

        genfit::Track* genfitTrack_;
        genfit::Track* backwardGenfitTrack_;
        
        extension::MutableTrack edm4hepTrack_;  

        const dd4hep::rec::DCH_info* _dch_info;
        const dd4hep::DDSegmentation::BitFieldCoder* _dc_decoder;

        TVector3 IP_referencePoint{0., 0., 0.};
        TVector3 firstHit_referencePoint;
        TVector3 lastHit_referencePoint;

        
    };

}  // namespace GenfitInterface

#endif // GENFIT_TRACK_H
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

#ifndef WIRE_MEASUREMENT_H
#define WIRE_MEASUREMENT_H

// Standard
#include <iostream>

// ROOT
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>

// GenFit
#include <TrackPoint.h>
#include <WirePointMeasurement.h>

// EDM4hep
#include "extension/SenseWireHitCollection.h"
#include "DDSegmentation/BitFieldCoder.h"
#include "DDRec/DCH_info.h"

/** @class Wire_measurement
 *
 *  Interface class that wraps a SenseWireHit object into a GENFIT-compatible WirePointMeasurement.
 *  
 *  This class converts a hit in a drift chamber wire (represented by `extension::SenseWireHit`) into a 
 *  `genfit::WirePointMeasurement`, which can be used in GENFIT's track fitting. The conversion requires 
 *  additional detector information such as geometry and segmentation, provided through DD4hep.
 *
 *  The constructor initializes the GENFIT measurement using drift radius and wire geometry, and assigns 
 *  detector and hit indices for identification and debugging purposes.
 *
 *
 *  Author: Andrea De Vita  
 *  Date  : 2025-06
 *
 */

namespace GenfitInterface {

    class Wire_measurement {
    public:

        // Constructor: takes a SenseWireHit and creates a WirePointMeasurement
        Wire_measurement(const extension::SenseWireHit& hit, const dd4hep::rec::DCH_info* dch_info, const dd4hep::DDSegmentation::BitFieldCoder* decoder, const int det_idx, const int hit_idx,const int debug_lvl);

        // Getter for genfit::WirePointMeasurement
        genfit::WirePointMeasurement* getGenFit() const;

    private:
        genfit::WirePointMeasurement* genfitHit_; 
    };

} // namespace GenfitInterface

#endif 
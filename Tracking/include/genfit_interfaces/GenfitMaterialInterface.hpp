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

#ifndef GENFIT_MATERIALINTERFACE_H
#define GENFIT_MATERIALINTERFACE_H

// GenFit
#include "AbsMaterialInterface.h"
#include "Exception.h"
#include "MaterialEffects.h"

// ROOT
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"

// Gaudi
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/ISvcLocator.h"

// DD4hep
#include "DD4hep/Detector.h"

namespace GenfitInterface{

        class GenfitMaterialInterface : public genfit::AbsMaterialInterface{
        public:

                GenfitMaterialInterface(const dd4hep::Detector* dd4hepGeo);
                virtual ~GenfitMaterialInterface(){};
                static GenfitMaterialInterface* getInstance(const dd4hep::Detector* dd4hepGeo);
                void destruct();

                void setDetector(dd4hep::Detector*);
                bool initTrack(double posX, double posY, double posZ, double dirX, double dirY, double dirZ) override;

                genfit::Material getMaterialParameters() override;

                double findNextBoundary(const genfit::RKTrackRep* rep,
                                        const genfit::M1x7& state7,
                                        double sMax,
                                        bool varField = true) override;

                void setMinSafetyDistanceCut(double safeDistCut=1e-7) {m_safeDistCut=safeDistCut;}
                virtual void setDebugLvl(unsigned int lvl = 1) {debugLvl_ = lvl;}

        private:

                static GenfitMaterialInterface* m_instance;
                TGeoManager* m_geoManager;
                double m_safeDistCut;

                TGeoManager* getGeoManager();
                double getSafeDistance();
                double getStep();
                TGeoNode* findNextBoundary(double stepmax, const char* path="", bool frombdr=false);
                bool isSameLocation(double posX, double posY, double posZ, bool change=false);
                void setCurrentDirection(double nx, double ny, double nz);
               

        };

}

#endif // GenfitMaterialInterface

/**************************************************************************
  * basf2 (Belle II Analysis Software Framework)                           *
  * Author: The Belle II Collaboration                                     *
  *                                                                        *
  * See git log for contributors and copyright holders.                    *
  * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
  **************************************************************************/
  
 #pragma once
  
 #include <cdc/dataobjects/CDCHit.h>
 #include <cdc/translators/ADCCountTranslatorBase.h>
 #include <cdc/translators/CDCGeometryTranslatorBase.h>
 #include <cdc/translators/TDCCountTranslatorBase.h>
  
 #include <genfit/AbsMeasurement.h>
 #include <genfit/MeasurementOnPlane.h>
 #include <genfit/TrackCandHit.h>
 #include <genfit/HMatrixU.h>
  
 #include <framework/geometry/B2Vector3.h>
  
 #include <memory>
  
  
 namespace Belle2 {
   class CDCRecoHit : public genfit::AbsMeasurement  {
  
   public:
     CDCRecoHit();
  
     CDCRecoHit(const CDCHit* cdcHit, const genfit::TrackCandHit* trackCandHit);
  
     ~CDCRecoHit() {}
  
     CDCRecoHit* clone() const override;
  
     WireID getWireID() const
     {
       return m_wireID;
     }
  
     static void setTranslators(CDC::ADCCountTranslatorBase*    const adcCountTranslator,
                                CDC::CDCGeometryTranslatorBase* const cdcGeometryTranslator,
                                CDC::TDCCountTranslatorBase*    const tdcCountTranslator,
                                //temp4cosmics                               bool useTrackTime = false);
                                bool useTrackTime = false, bool cosmics = false);
  
     genfit::SharedPlanePtr constructPlane(const genfit::StateOnPlane& state) const override;
  
     std::vector<genfit::MeasurementOnPlane*> constructMeasurementsOnPlane(const genfit::StateOnPlane& state) const override;
  
     virtual const genfit::HMatrixU* constructHMatrix(const genfit::AbsTrackRep*) const override;
  
     std::vector<double> timeDerivativesMeasurementsOnPlane(const genfit::StateOnPlane& state) const;
  
     bool getFlyByDistanceVector(B2Vector3D& pointingVector, B2Vector3D& trackDir,
                                 const genfit::AbsTrackRep* rep = nullptr,
                                 bool usePlaneFromFit = false);
  
     void setLeftRightResolution(int lr) { m_leftRight = lr; }
  
     bool isLeftRightMeasurement() const override { return true; }
  
     int getLeftRightResolution() const override { return m_leftRight; }
  
  
     const CDCHit* getCDCHit() const
     {
       return m_cdcHit;
     }
  
   protected:
 #ifndef __CINT__ // rootcint doesn't know smart pointers
     static std::unique_ptr<CDC::ADCCountTranslatorBase>     s_adcCountTranslator;
  
     static std::unique_ptr<CDC::CDCGeometryTranslatorBase>  s_cdcGeometryTranslator;
  
     static std::unique_ptr<CDC::TDCCountTranslatorBase>    s_tdcCountTranslator;
  
     static bool s_useTrackTime;
     static bool s_cosmics;
  
 #endif
  
     unsigned short m_tdcCount;
  
     unsigned short m_adcCount;
  
     WireID m_wireID;
  
     const CDCHit* m_cdcHit;  
  
     signed char m_leftRight;
  
     ClassDefOverride(CDCRecoHit, 10);
     // Version history:
     // ver 10: ICalibrationParametersDerivatives interface moved to derived class.
     //         ClassDef -> ClassDefOverride + consistent override keyword usage.
     //         Private members -> protected for access from derived class.
     // ver 9: Derives from ICalibrationParametersDerivatives to expose
     //        alignment/calibration interface
     // ver 8: Rewrite to deal with realistic translators.  No longer
     //        derives from genfit::WireMeasurement.
   };
 }
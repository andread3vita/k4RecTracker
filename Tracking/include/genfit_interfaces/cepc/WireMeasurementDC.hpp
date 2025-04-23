#ifndef WireMeasurementDC_h
#define WireMeasurementDC_h

#include "AbsMeasurement.h"
#include "TrackPoint.h"
#include "SharedPlanePtr.h"
#include "MeasurementOnPlane.h"
#include "StateOnPlane.h"
#include "AbsHMatrix.h"
// #include "GenfitHit.h"


#include "extension/TrackerHit.h"
namespace extension{
    class TrackerHit;
}


using namespace genfit;
class WireMeasurementDC : public genfit::AbsMeasurement{

 public:
  WireMeasurementDC();
  WireMeasurementDC(double driftDistance, double driftDistanceError, const TVector3& endPoint1, const TVector3& endPoint2, int detId, int hitId, TrackPoint* trackPoint);
  // WireMeasurementDC(const GenfitHit* genfitHit,int iHit);

  virtual ~WireMeasurementDC() {;}

  virtual WireMeasurementDC* clone() const override {return new WireMeasurementDC(*this);}

  virtual SharedPlanePtr constructPlane(const StateOnPlane& state) const override;

  /**  Hits with a small drift distance get a higher weight, whereas hits with
    * big drift distances become weighted down.
    * When these initial weights are used by the DAF, the smoothed track will be closer to the real
    * trajectory than if both sides are weighted with 0.5 regardless of the drift distance.
    * This helps a lot when resolving l/r ambiguities with the DAF.
    * The idea is that for the first iteration of the DAF, the wire positions are taken.
    * For small drift radii, the wire position does not bend the fit away from the
    * trajectory, whereas the wire position for hits with large drift radii is further away
    * from the trajectory and will therefore bias the fit if not weighted down.
    */
  virtual std::vector<MeasurementOnPlane*> constructMeasurementsOnPlane(const StateOnPlane& state) const override;

  virtual const AbsHMatrix* constructHMatrix(const AbsTrackRep*) const override;

  /** Reset the wire end points.
   */
  void setWireEndPoints(const TVector3& endPoint1, const TVector3& endPoint2);

  /** Set maximum drift distance. This is used to calculate the start weights of the two
   * measurementsOnPlane.
   */
  void setMaxDistance(double d){maxDistance_ = d;}
  /**
   * select how to resolve the left/right ambiguity:
   * -1: negative (left) side on vector (wire direction) x (track direction)
   * 0: mirrors enter with same weight, DAF will decide.
   * 1: positive (right) side on vector (wire direction) x (track direction)
   * where the wire direction is pointing from endPoint1 to endPoint2
   */
  void setLeftRightResolution(int lr);

  virtual bool isLeftRigthMeasurement() const {return true;}
  double getMaxDistance(){return maxDistance_;}
  int getLeftRightResolution() const override {return leftRight_;}
  int getWireMeasurementDCLayer(){return layer_;}
  void setTrackerHit(const extension::TrackerHit& trackerHit,int l,int c,double t){
      trackerHit_=&trackerHit;
      layer_=l;
      cell_=c;
      time_=t;
  }
  void setTrackerHitObj(const extension::TrackerHit& trackerHit,int l,int c,double t){
      trackerHit_= nullptr;
      trackerHitObj_=trackerHit;
      layer_=l;
      cell_=c;
      time_=t;
  }
  // void setSimTrackerHit(const edm4hep::SimTrackerHit& simTrackerHit){
  //     simTrackerHit_=&simTrackerHit;
  // }
  void print();

  const extension::TrackerHit* getTrackerHit(){return trackerHit_;}
  extension::TrackerHit getTrackerHitObj(){return trackerHitObj_;}
  // const GenfitHit* getGenfitHit(){return genfitHit_;}

 protected:

  double wireEndPoint1_[3]; //! Wire end point 1 (X, Y, Z)
  double wireEndPoint2_[3]; //! Wire end point 2 (X, Y, Z)
  double maxDistance_;
  double leftRight_;

  int layer_;
  int cell_;
  double time_;
  // const edm4hep::SimTrackerHit* simTrackerHit_;
  const extension::TrackerHit* trackerHit_;
  extension::TrackerHit trackerHitObj_;
  // const GenfitHit* genfitHit_;

 public:

  //ClassDefOverride(WireMeasurementDC, 1)

};

/** @} */

#endif // WireMeasurementDC_h
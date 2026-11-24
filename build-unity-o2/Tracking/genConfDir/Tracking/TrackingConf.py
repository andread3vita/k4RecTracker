#Mon Nov 16 11:16:21 2026"""Automatically generated. DO NOT EDIT please"""
from GaudiKernel.GaudiHandles import *
from GaudiKernel.Proxy.Configurable import *

class GGTFTrackFinder( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : {  },
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputPlanarHitCollections' : [ 'InputPlanarHitCollections' ],
    'InputWireHitCollections' : [ 'InputWireHitCollections' ],
    'OutputTracksGGTF' : [ 'OutputTracksGGTF' ],
    'ModelPath' : '',
    'EmbeddingCoordinates' : 3,
    'InferenceThreads' : 4,
    'MaxHits' : 20000,
    'SessionRecreateEvents' : 0,
    'Tbeta' : 0.6,
    'Td' : 0.3,
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputPlanarHitCollections' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection> (std::vector<edm4hep::TrackerHitPlaneCollection const*,std::allocator<edm4hep::TrackerHitPlaneCollection const*> > const&,std::vector<edm4hep::SenseWireHitCollection const*,std::allocator<edm4hep::SenseWireHitCollection const*> > const&),Gaudi::Functional::Traits::use_<> >] """,
    'InputWireHitCollections' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection> (std::vector<edm4hep::TrackerHitPlaneCollection const*,std::allocator<edm4hep::TrackerHitPlaneCollection const*> > const&,std::vector<edm4hep::SenseWireHitCollection const*,std::allocator<edm4hep::SenseWireHitCollection const*> > const&),Gaudi::Functional::Traits::use_<> >] """,
    'OutputTracksGGTF' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection> (std::vector<edm4hep::TrackerHitPlaneCollection const*,std::allocator<edm4hep::TrackerHitPlaneCollection const*> > const&,std::vector<edm4hep::SenseWireHitCollection const*,std::allocator<edm4hep::SenseWireHitCollection const*> > const&),Gaudi::Functional::Traits::use_<> >] """,
    'ModelPath' : """ Path to the ONNX model file [GGTFTrackFinder] """,
    'EmbeddingCoordinates' : """ Number of embedding-space coordinates per hit in the ONNX model output; the model must produce one additional trailing value per hit, the beta logit (the sigmoid of the logit is the cluster potential compared against Tbeta) [GGTFTrackFinder] """,
    'InferenceThreads' : """ Number of intra-operation threads used by the ONNX Runtime session; inference time plateaus beyond about 4-8 threads for the current GATr model, so increase further only if the job has many free cores [GGTFTrackFinder] """,
    'MaxHits' : """ Maximum number of hits accepted per event [GGTFTrackFinder] """,
    'SessionRecreateEvents' : """ Recreate the ONNX session every N events (0 disables); bounds the per-input-shape allocation state ONNX Runtime retains for ever-changing hit counts, at the price of one session reload per interval [GGTFTrackFinder] """,
    'Tbeta' : """ Threshold on sigmoid(beta logit) used to identify cluster core points [GGTFTrackFinder] """,
    'Td' : """ Radius used to assign nearby hits to a cluster core [GGTFTrackFinder] """,
  }
  __declaration_location__ = 'GGTFTrackFinder.cpp:654'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(GGTFTrackFinder, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'Tracking'
  def getType( self ):
      return 'GGTFTrackFinder'
  pass # class GGTFTrackFinder

class GenfitTrackFitter( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : {  },
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputTracks' : [ 'InputTracks' ],
    'OutputFittedTracks' : [ 'Fitted_tracks' ],
    'OutputFittedTracksWithFilteredHits' : [ 'Fitted_tracks_with_filtered_hits' ],
    'OutputFittedHits' : [ 'Fitted_hits' ],
    'GeoSvc' : ServiceHandle('GeoSvc'),
    'WireTrackerName' : 'DCH_v2',
    'UseBrems' : False,
    'FitterType' : 'DAF',
    'BetaInit' : 100.0,
    'BetaFinal' : 0.05,
    'BetaSteps' : 15,
    'ParticleHypothesisList' : [ 211 ],
    'UseFirstHitAsReference' : False,
    'InitializationType' : 1,
    'Sigma_d0' : 0.05,
    'Sigma_phi' : 0.1,
    'OmegaFactor' : 0.5,
    'Z0Factor' : 0.1,
    'Sigma_tanLambda' : 0.1,
    'TrackStateLocation' : 2,
    'InitPosition' : [ 0 , 0 , 0 ],
    'InitMomentum' : [ 0 , 0 , 0 ],
    'Epsilon' : 1e-04,
    'SmoothWindow' : 5,
    'SkipTrackOrdering' : False,
    'ListOfTypesToSkip' : [  ],
    'RunSingleEvaluation' : False,
    'FilterTrackHits' : True,
    'RunCalorimeterExtrapolation' : True,
    'RadialThresholdPromptTrack' : 100.0,
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputTracks' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection,edm4hep::TrackCollection,edm4hep::TrackerHitPlaneCollection> (edm4hep::TrackCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'OutputFittedTracks' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection,edm4hep::TrackCollection,edm4hep::TrackerHitPlaneCollection> (edm4hep::TrackCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'OutputFittedTracksWithFilteredHits' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection,edm4hep::TrackCollection,edm4hep::TrackerHitPlaneCollection> (edm4hep::TrackCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'OutputFittedHits' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection,edm4hep::TrackCollection,edm4hep::TrackerHitPlaneCollection> (edm4hep::TrackCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'GeoSvc' : """ Detector geometry service [GenfitTrackFitter] """,
    'WireTrackerName' : """ WireTrackerName in the detector description (used to retrieve Wire Tracker geometry and material information) [GenfitTrackFitter] """,
    'UseBrems' : """ Include Bremsstrahlung energy loss and noise in the fit [GenfitTrackFitter] """,
    'FitterType' : """ Fitter type to use [https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf]: 'DAF', 'KALMAN', or 'KALMAN_REF' [GenfitTrackFitter] """,
    'BetaInit' : """ Initial annealing parameter (if m_Fitter_type == 'DAF') [https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf] [GenfitTrackFitter] """,
    'BetaFinal' : """ Final annealing parameter (if m_Fitter_type == 'DAF') [https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf] [GenfitTrackFitter] """,
    'BetaSteps' : """ Number of steps in the annealing schedule (if m_Fitter_type == 'DAF') [https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf] [GenfitTrackFitter] """,
    'ParticleHypothesisList' : """ List of particle hypotheses to consider: {11,13,211,321,2212} -> e, mu, pi, K, p [GenfitTrackFitter] """,
    'UseFirstHitAsReference' : """ If true, the first hit is used as the reference point for track initialization [GenfitTrackFitter] """,
    'InitializationType' : """ Method for initializing track parameters before fitting: 0: first two hits (position from first hit, momentum from direction between first and second hit); 1: refined initialization using ComputeInitialParameters(Bz); 2: initialize from a track state (position from reference point, momentum from helix parameters); 3: use user-provided InitPosition and InitMomentum [GenfitTrackFitter] """,
    'Sigma_d0' : """ Initial uncertainty on d0 (mm) for the initial covariance matrix when InitializationType = 0,1,3. [GenfitTrackFitter] """,
    'Sigma_phi' : """ Initial uncertainty on phi (rad) for the initial covariance matrix when InitializationType = 0,1,3. [GenfitTrackFitter] """,
    'OmegaFactor' : """ Scaling factor for omega uncertainty for the initial covariance matrix when InitializationType = 0,1,3.The actual sigma_omega is computed as OmegaFactor * |omega| [GenfitTrackFitter] """,
    'Z0Factor' : """ Scaling factor for z0 uncertainty for the initial covariance matrix when InitializationType = 0,1,3.The actual sigma_z0 is computed as Z0Factor * |z0 (in mm)| [GenfitTrackFitter] """,
    'Sigma_tanLambda' : """ Initial uncertainty on tanLambda for the initial covariance matrix when InitializationType = 0,1,3. [GenfitTrackFitter] """,
    'TrackStateLocation' : """ TrackState location used for initialization when InitializationType = 2. Defines where the reference point and helix parameters are taken from [GenfitTrackFitter] """,
    'InitPosition' : """ User-defined initial position [in mm] for InitializationType = 3 [GenfitTrackFitter] """,
    'InitMomentum' : """ User-defined initial momentum [in GeV/c] for InitializationType = 3 [GenfitTrackFitter] """,
    'Epsilon' : """ Threshold to detect sign changes in the first derivative during track initialization. Used to identify hits from the first round of the looper [GenfitTrackFitter] """,
    'SmoothWindow' : """ Number of hits used for smoothing before computing the first derivative in track initialization [GenfitTrackFitter] """,
    'SkipTrackOrdering' : """ Skip hit ordering before fitting [GenfitTrackFitter] """,
    'ListOfTypesToSkip' : """ List of track types to skip during fitting (e.g. unmatched tracks with type = 0) [GenfitTrackFitter] """,
    'RunSingleEvaluation' : """ If true, only the first particle hypothesis is evaluated. If false, all hypotheses are scanned and the best fit is chosen based on chi2/ndf [GenfitTrackFitter] """,
    'FilterTrackHits' : """ Filter track hits after fitting based on their quality (weights assigned by the fitter). This also resolves the left/right ambiguity in drift chamber hits: only the hit with the higher weight is considered. [GenfitTrackFitter] """,
    'RunCalorimeterExtrapolation' : """ If true, the track is extrapolated to the calorimeter surfaces (barrel and endcap)and the corresponding track states are stored in the output track.This is applied only if a non-zero magnetic field (Bz > 0) is found at the last hit position. [GenfitTrackFitter] """,
    'RadialThresholdPromptTrack' : """ Radius [mm] defining a spherical region centered at {0,0,0}. Tracks produced within this radius are classified as prompt, while those outside are considered displaced. This parameter is used to select different track initialization strategies: prompt tracks are initialized using the PCA of the hits with respect to a reference point (IP or VP), whereas displaced tracks are initialized using the PCA of the first hit. [GenfitTrackFitter] """,
  }
  __declaration_location__ = 'GenfitTrackFitter.cpp:834'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(GenfitTrackFitter, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'Tracking'
  def getType( self ):
      return 'GenfitTrackFitter'
  pass # class GenfitTrackFitter

class PerfectTrackFinder( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : {  },
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputPlanarHitCollections' : [ 'InputPlanarHitCollections' ],
    'InputWireHitCollections' : [ 'InputWireHitCollections' ],
    'InputMCParticles' : [ 'InputMCParticles' ],
    'OutputPerfectTracks' : [ 'OutputPerfectTracks' ],
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputPlanarHitCollections' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection> (std::vector<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*,std::allocator<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*> > const&,std::vector<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*,std::allocator<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*> > const&,edm4hep::MCParticleCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'InputWireHitCollections' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection> (std::vector<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*,std::allocator<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*> > const&,std::vector<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*,std::allocator<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*> > const&,edm4hep::MCParticleCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'InputMCParticles' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection> (std::vector<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*,std::allocator<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*> > const&,std::vector<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*,std::allocator<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*> > const&,edm4hep::MCParticleCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'OutputPerfectTracks' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection> (std::vector<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*,std::allocator<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*> > const&,std::vector<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*,std::allocator<podio::LinkCollection<edm4hep::TrackerHit,edm4hep::SimTrackerHit> const*> > const&,edm4hep::MCParticleCollection const&),Gaudi::Functional::Traits::use_<> >] """,
  }
  __declaration_location__ = 'PerfectTrackFinder.cpp:164'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(PerfectTrackFinder, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'Tracking'
  def getType( self ):
      return 'PerfectTrackFinder'
  pass # class PerfectTrackFinder

class PlotTrackHitDistances( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : {  },
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputSimTrackerHits' : [ 'DCHCollection' ],
    'InputTracksFromGenParticlesAssociation' : [ 'TracksFromGenParticlesAssociation' ],
    'Bz' : 2.0,
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputSimTrackerHits' : """  [k4FWCore::details::Consumer<void (edm4hep::SimTrackerHitCollection const&,podio::LinkCollection<edm4hep::Track,edm4hep::MCParticle> const&),Gaudi::Functional::Traits::use_<> >] """,
    'InputTracksFromGenParticlesAssociation' : """  [k4FWCore::details::Consumer<void (edm4hep::SimTrackerHitCollection const&,podio::LinkCollection<edm4hep::Track,edm4hep::MCParticle> const&),Gaudi::Functional::Traits::use_<> >] """,
    'Bz' : """ Z component of the (assumed constant) magnetic field in Tesla. [PlotTrackHitDistances] """,
  }
  __declaration_location__ = 'PlotTrackHitResiduals.cpp:83'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(PlotTrackHitDistances, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'Tracking'
  def getType( self ):
      return 'PlotTrackHitDistances'
  pass # class PlotTrackHitDistances

class TrackMerger( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : {  },
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputInnerTracks' : 'InnerTracks',
    'InputOuterTracks' : 'OuterTracks',
    'OutTracks' : 'MyCandidateMergedTracks',
    'Greedy' : True,
    'UseSignificance' : True,
    'D0Tolerance' : 0.5,
    'Z0Tolerance' : 2.5,
    'PhiTolerance' : -1.0,
    'OmegaTolerance' : -1.0,
    'TanLambdaTolerance' : -1.0,
    'D0Significance' : 3.0,
    'Z0Significance' : 3.0,
    'PhiSignificance' : 3.0,
    'PtSignificance' : 2.0,
    'TanLambdaSignificance' : 3.0,
    'GeoSvc' : ServiceHandle('GeoSvc'),
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputInnerTracks' : """  [k4FWCore::details::Transformer<edm4hep::TrackCollection (edm4hep::TrackCollection const&,edm4hep::TrackCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'InputOuterTracks' : """  [k4FWCore::details::Transformer<edm4hep::TrackCollection (edm4hep::TrackCollection const&,edm4hep::TrackCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'OutTracks' : """  [k4FWCore::details::Transformer<edm4hep::TrackCollection (edm4hep::TrackCollection const&,edm4hep::TrackCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'Greedy' : """ If true, each track is used only once. [TrackMerger] """,
    'UseSignificance' : """ If true (default), match on the *Significance cuts, i.e. |diff| / sqrt(var(inner) + var(outer)) using the track states' covariance matrices. If false, match on the *Tolerance cuts (absolute differences). [TrackMerger] """,
    'D0Tolerance' : """ Absolute mode: maximum allowed |D0(inner) - D0(outer)| [mm] for a match. Negative disables this criterion. [TrackMerger] """,
    'Z0Tolerance' : """ Absolute mode: maximum allowed |Z0(inner) - Z0(outer)| [mm] for a match. Negative disables this criterion. [TrackMerger] """,
    'PhiTolerance' : """ Absolute mode: maximum allowed |phi(inner) - phi(outer)| [rad], wrapped into [-pi, pi], for a match. Negative disables this criterion. [TrackMerger] """,
    'OmegaTolerance' : """ Absolute mode: maximum allowed |omega(inner) - omega(outer)| [1/mm] for a match. Negative disables this criterion. [TrackMerger] """,
    'TanLambdaTolerance' : """ Absolute mode: maximum allowed |tanLambda(inner) - tanLambda(outer)| for a match. Negative disables this criterion. [TrackMerger] """,
    'D0Significance' : """ Significance mode: maximum allowed D0 significance for a match. Negative disables this criterion. [TrackMerger] """,
    'Z0Significance' : """ Significance mode: maximum allowed Z0 significance for a match. Negative disables this criterion. [TrackMerger] """,
    'PhiSignificance' : """ Significance mode: maximum allowed phi significance (difference wrapped into [-pi, pi]) for a match. Negative disables this criterion. [TrackMerger] """,
    'PtSignificance' : """ Significance mode: maximum allowed significance of the signed pT (charge sign times pT, derived from omega and the magnetic field) for a match. Requires a geometry with a field. Negative disables this criterion. [TrackMerger] """,
    'TanLambdaSignificance' : """ Significance mode: maximum allowed tanLambda significance for a match. Negative disables this criterion. [TrackMerger] """,
    'GeoSvc' : """ Detector geometry service [TrackMerger] """,
  }
  __declaration_location__ = 'TrackMerger.cpp:289'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(TrackMerger, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'Tracking'
  def getType( self ):
      return 'TrackMerger'
  pass # class TrackMerger

class TrackdNdxDelphesBased( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : {  },
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputLinkCollection' : [ 'TrackMCParticleLinks' ],
    'HeaderName' : [ 'EventHeader' ],
    'OutputCollection' : [ 'RecDqdxCollection' ],
    'ZmaxParameterName' : 'DCH_gas_Lhalf',
    'ZminParameterName' : 'DCH_gas_Lhalf',
    'RminParameterName' : 'DCH_gas_inner_cyl_R',
    'RmaxParameterName' : 'DCH_gas_outer_cyl_R',
    'GasSel' : 0,
    'FillFactor' : 1.0,
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputLinkCollection' : """  [k4FWCore::details::Transformer<edm4hep::RecDqdxCollection (podio::LinkCollection<edm4hep::Track,edm4hep::MCParticle> const&,edm4hep::EventHeaderCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'HeaderName' : """  [k4FWCore::details::Transformer<edm4hep::RecDqdxCollection (podio::LinkCollection<edm4hep::Track,edm4hep::MCParticle> const&,edm4hep::EventHeaderCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'OutputCollection' : """  [k4FWCore::details::Transformer<edm4hep::RecDqdxCollection (podio::LinkCollection<edm4hep::Track,edm4hep::MCParticle> const&,edm4hep::EventHeaderCollection const&),Gaudi::Functional::Traits::use_<> >] """,
    'ZmaxParameterName' : """ Name of XML file parameter describing the z extent of the active volume (max value in +z direction) [TrackdNdxDelphesBased] """,
    'ZminParameterName' : """ Name of XML file parameter describing the z extent of the active volume (min value in -z direction). For forward-backward symmetric detectors, use same name as for ZmaxParameterName. This value is then automatically converted into the negative value. [TrackdNdxDelphesBased] """,
    'RminParameterName' : """ Name of XML file parameter describing the inner radius of the active volume. [TrackdNdxDelphesBased] """,
    'RmaxParameterName' : """ Name of XML file parameter describing the outer radius of the active volume. [TrackdNdxDelphesBased] """,
    'GasSel' : """ Gas selection: 0: He(90%)-Isobutane(10%), 1: pure He, 2: Ar(50%)-Ethane(50%), 3: pure Ar. [TrackdNdxDelphesBased] """,
    'FillFactor' : """ Factor (between 0 and 1) describing the fraction of the detector volume that is active (e.g., for Straw Tube Tracker the factor is significantly below 1 due to gaps between tubes). The factor is used to scale the calculated track length and thus the corresponding number of clusters. [TrackdNdxDelphesBased] """,
  }
  __declaration_location__ = 'TrackdNdxDelphesBased.cpp:21'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(TrackdNdxDelphesBased, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'Tracking'
  def getType( self ):
      return 'TrackdNdxDelphesBased'
  pass # class TrackdNdxDelphesBased

class TracksFromGenParticles( ConfigurableAlgorithm ) :
  __slots__ = { 
    'ExtraInputs' : set(),
    'ExtraOutputs' : set(),
    'OutputLevel' : 0,
    'Enable' : True,
    'AuditInitialize' : False,
    'AuditReinitialize' : False,
    'AuditRestart' : False,
    'AuditExecute' : False,
    'AuditFinalize' : False,
    'AuditStart' : False,
    'AuditStop' : False,
    'Timeline' : True,
    'MonitorService' : 'MonitorSvc',
    'RegisterForContextService' : False,
    'Cardinality' : 0,
    'NeededResources' : {  },
    'Asynchronous' : False,
    'FilterCircularDependencies' : True,
    'InputGenParticles' : [ 'MCParticles' ],
    'InputSimTrackerHits' : [ 'SimTrackerHits' ],
    'OutputTracks' : [ 'TracksFromGenParticles' ],
    'OutputMCRecoTrackParticleAssociation' : [ 'TracksFromGenParticlesAssociation' ],
    'ExtrapolateToECal' : False,
    'KeepOnlyBestExtrapolation' : True,
    'MinimumParticleMomentum' : 0.01,
    'OnlyCaloReachingParticles' : True,
    'TrackerIDs' : [  ],
    'SystemEncoding' : 'system:5',
  }
  _propertyDocDct = { 
    'ExtraInputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'ExtraOutputs' : """  [DataHandleHolderBase<PropertyHolder<CommonMessaging<implements<IAlgorithm,IDataHandleHolder,IProperty,IStateful> > > >] """,
    'OutputLevel' : """ output level [Gaudi::Algorithm] """,
    'Enable' : """ should the algorithm be executed or not [Gaudi::Algorithm] """,
    'AuditInitialize' : """ trigger auditor on initialize() [Gaudi::Algorithm] """,
    'AuditReinitialize' : """ trigger auditor on reinitialize() [Gaudi::Algorithm] """,
    'AuditRestart' : """ trigger auditor on restart() [Gaudi::Algorithm] """,
    'AuditExecute' : """ trigger auditor on execute() [Gaudi::Algorithm] """,
    'AuditFinalize' : """ trigger auditor on finalize() [Gaudi::Algorithm] """,
    'AuditStart' : """ trigger auditor on start() [Gaudi::Algorithm] """,
    'AuditStop' : """ trigger auditor on stop() [Gaudi::Algorithm] """,
    'Timeline' : """ send events to TimelineSvc [Gaudi::Algorithm] """,
    'MonitorService' : """ name to use for Monitor Service [Gaudi::Algorithm] """,
    'RegisterForContextService' : """ flag to enforce the registration for Algorithm Context Service [Gaudi::Algorithm] """,
    'Cardinality' : """ how many clones to create - 0 means algo is reentrant [Gaudi::Algorithm] """,
    'NeededResources' : """ named resources needed during event looping [Gaudi::Algorithm] """,
    'Asynchronous' : """ whether algorithm is asynchronous and uses Boost Fiber to suspend while offloaded code is running. [Gaudi::Algorithm] """,
    'FilterCircularDependencies' : """ filter out circular data dependencies [Gaudi::Algorithm] """,
    'InputGenParticles' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection,podio::LinkCollection<edm4hep::Track,edm4hep::MCParticle> > (edm4hep::MCParticleCollection const&,std::vector<edm4hep::SimTrackerHitCollection const*,std::allocator<edm4hep::SimTrackerHitCollection const*> > const&),Gaudi::Functional::Traits::use_<> >] """,
    'InputSimTrackerHits' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection,podio::LinkCollection<edm4hep::Track,edm4hep::MCParticle> > (edm4hep::MCParticleCollection const&,std::vector<edm4hep::SimTrackerHitCollection const*,std::allocator<edm4hep::SimTrackerHitCollection const*> > const&),Gaudi::Functional::Traits::use_<> >] """,
    'OutputTracks' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection,podio::LinkCollection<edm4hep::Track,edm4hep::MCParticle> > (edm4hep::MCParticleCollection const&,std::vector<edm4hep::SimTrackerHitCollection const*,std::allocator<edm4hep::SimTrackerHitCollection const*> > const&),Gaudi::Functional::Traits::use_<> >] """,
    'OutputMCRecoTrackParticleAssociation' : """  [k4FWCore::details::MultiTransformer<std::tuple<edm4hep::TrackCollection,podio::LinkCollection<edm4hep::Track,edm4hep::MCParticle> > (edm4hep::MCParticleCollection const&,std::vector<edm4hep::SimTrackerHitCollection const*,std::allocator<edm4hep::SimTrackerHitCollection const*> > const&),Gaudi::Functional::Traits::use_<> >] """,
    'ExtrapolateToECal' : """ Calculate track state at ECal inner face or not [TracksFromGenParticles] """,
    'KeepOnlyBestExtrapolation' : """ Keep only extrapolation with shortest time or not [TracksFromGenParticles] """,
    'MinimumParticleMomentum' : """ Keep only particles with momentum (in GeV) greater than MinimumParticleMomentum [TracksFromGenParticles] """,
    'OnlyCaloReachingParticles' : """ Build tracks only for charged particles reaching the calorimeter (not isDecayedInTracker) [TracksFromGenParticles] """,
    'TrackerIDs' : """ System IDs of tracking subdetectors [TracksFromGenParticles] """,
    'SystemEncoding' : """ System encoding string [TracksFromGenParticles] """,
  }
  __declaration_location__ = 'TracksFromGenParticles.cpp:468'
  def __init__(self, name = Configurable.DefaultName, **kwargs):
      super(TracksFromGenParticles, self).__init__(name)
      for n,v in kwargs.items():
         setattr(self, n, v)
  def getDlls( self ):
      return 'Tracking'
  def getType( self ):
      return 'TracksFromGenParticles'
  pass # class TracksFromGenParticles

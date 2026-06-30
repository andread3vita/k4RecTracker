# Example steering for the LeastSquaresVertexFitter Gaudi transformer.
#
# It reads a podio/EDM4hep file that already contains a reconstructed
# edm4hep::TrackCollection and fits those tracks to a common vertex, writing an
# edm4hep::VertexCollection to the output file.
#
# Adapt "Input", the track collection name and (optionally) the beam-spot
# settings to your sample.

import os
import math
from k4FWCore import ApplicationMgr, IOSvc
from Gaudi.Configuration import INFO, DEBUG
from Configurables import EventDataSvc, UniqueIDGenSvc, RndmGenSvc
from Configurables import GeoSvc
from k4FWCore.parseArgs import parser


################## Parser
parser.add_argument("--inputFile", help="InputFile")
parser.add_argument("--outputFile", help="OutputFile")
args = parser.parse_args()

# ################## InputOutput
svc = IOSvc("IOSvc")
svc.Input = args.inputFile
svc.Output = args.outputFile

################ Detector geometry
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
detectors_to_use = ["FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"]
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

from Configurables import DDPlanarDigi
from Configurables import DCHdigi_v02

############### Vertex Digitizer

innerVertexResolution_x = 0.003  # [mm], assume 3 µm resolution for ARCADIA sensor
innerVertexResolution_y = 0.003  # [mm], assume 3 µm resolution for ARCADIA sensor
innerVertexResolution_t = 1000  # [ns]
outerVertexResolution_x = 0.050 / math.sqrt(12)  # [mm], assume ATLASPix3 sensor with 50 µm pitch
outerVertexResolution_y = 0.150 / math.sqrt(12)  # [mm], assume ATLASPix3 sensor with 150 µm pitch
outerVertexResolution_t = 1000  # [ns]

vtxb_digitizer = DDPlanarDigi("VTXBdigitizer")
vtxb_digitizer.SubDetectorName = "Vertex"
vtxb_digitizer.IsStrip = False
vtxb_digitizer.ResolutionU = [
    innerVertexResolution_x,
    innerVertexResolution_x,
    innerVertexResolution_x,
    outerVertexResolution_x,
    outerVertexResolution_x,
]
vtxb_digitizer.ResolutionV = [
    innerVertexResolution_y,
    innerVertexResolution_y,
    innerVertexResolution_y,
    outerVertexResolution_y,
    outerVertexResolution_y,
]
vtxb_digitizer.ResolutionT = [
    innerVertexResolution_t,
    innerVertexResolution_t,
    innerVertexResolution_t,
    outerVertexResolution_t,
    outerVertexResolution_t,
]
vtxb_digitizer.SimTrackHitCollectionName = ["VertexBarrelCollection"]
vtxb_digitizer.SimTrkHitRelCollection = ["VTXBSimDigiLinks"]
vtxb_digitizer.TrackerHitCollectionName = ["VTXBDigis"]
vtxb_digitizer.ForceHitsOntoSurface = True

vtxd_digitizer = DDPlanarDigi("VTXDdigitizer")
vtxd_digitizer.SubDetectorName = "Vertex"
vtxd_digitizer.IsStrip = False
vtxd_digitizer.ResolutionU = [
    outerVertexResolution_x,
    outerVertexResolution_x,
    outerVertexResolution_x,
]
vtxd_digitizer.ResolutionV = [
    outerVertexResolution_y,
    outerVertexResolution_y,
    outerVertexResolution_y,
]
vtxd_digitizer.ResolutionT = [
    outerVertexResolution_t,
    outerVertexResolution_t,
    outerVertexResolution_t,
]
vtxd_digitizer.SimTrackHitCollectionName = ["VertexEndcapCollection"]
vtxd_digitizer.SimTrkHitRelCollection = ["VTXDSimDigiLinks"]
vtxd_digitizer.TrackerHitCollectionName = ["VTXDDigis"]
vtxd_digitizer.ForceHitsOntoSurface = True

############### Wrapper Digitizer
siWrapperResolution_x = 0.050 / math.sqrt(12)  # [mm]
siWrapperResolution_y = 1.0 / math.sqrt(12)  # [mm]
siWrapperResolution_t = 0.040  # [ns], assume 40 ps timing resolution for a single layer -> Should lead to <30 ps resolution when >1 hit

siwrb_digitizer = DDPlanarDigi("SiWrBdigitizer")
siwrb_digitizer.SubDetectorName = "SiWrB"
siwrb_digitizer.IsStrip = False
siwrb_digitizer.ResolutionU = [siWrapperResolution_x, siWrapperResolution_x]
siwrb_digitizer.ResolutionV = [siWrapperResolution_y, siWrapperResolution_y]
siwrb_digitizer.ResolutionT = [siWrapperResolution_t, siWrapperResolution_t]
siwrb_digitizer.SimTrackHitCollectionName = ["SiWrBCollection"]
siwrb_digitizer.SimTrkHitRelCollection = ["SiWrBSimDigiLinks"]
siwrb_digitizer.TrackerHitCollectionName = ["SiWrBDigis"]
siwrb_digitizer.ForceHitsOntoSurface = True

siwrd_digitizer = DDPlanarDigi("SiWrDdigitizer")
siwrd_digitizer.SubDetectorName = "SiWrD"
siwrd_digitizer.IsStrip = False
siwrd_digitizer.ResolutionU = [siWrapperResolution_x, siWrapperResolution_x]
siwrd_digitizer.ResolutionV = [siWrapperResolution_y, siWrapperResolution_y]
siwrd_digitizer.ResolutionT = [siWrapperResolution_t, siWrapperResolution_t]
siwrd_digitizer.SimTrackHitCollectionName = ["SiWrDCollection"]
siwrd_digitizer.SimTrkHitRelCollection = ["SiWrDSimDigiLinks"]
siwrd_digitizer.TrackerHitCollectionName = ["SiWrDDigis"]
siwrd_digitizer.ForceHitsOntoSurface = True


############### DCH Digitizer

dch_digitizer = DCHdigi_v02(
    "DCHdigi2",
    InputSimHitCollection=["DCHCollection"],
    OutputDigihitCollection=["DCH_DigiCollection"],
    OutputLinkCollection=["DCH_DigiSimAssociationCollection"],
    DCH_name="DCH_v2",
    zResolution_mm=30.0,  # in mm
    xyResolution_mm=0.1,  # in mm
    Deadtime_ns=400.0,  # in ns
    GasType=0,  # 0: He(90%)-Isobutane(10%), 1: pure He, 2: Ar(50%)-Ethane(50%), 3: pure Ar
    ReadoutWindowStartTime_ns=1.0,  # in ns (taking into account time of flight, drift, and signal travel)
    ReadoutWindowDuration_ns=450.0,  # in ns
    DriftVelocity_um_per_ns=-1.0,  # in um/ns, if negative, automatically chosen based on GasType
    SignalVelocity_mm_per_ns=200.0,  # in mm/ns (Default: 2/3 of the speed of light)
    OutputLevel=INFO,
)

############### Perfect Track Finder
from Configurables import PerfectTrackFinder

perfectFinder = PerfectTrackFinder(
    "PerfectTracking",
    InputPlanarHitCollections=[
        "VTXDSimDigiLinks",
        "VTXBSimDigiLinks",
        "SiWrBSimDigiLinks",
        "SiWrDSimDigiLinks",
    ],
    InputWireHitCollections=["DCH_DigiSimAssociationCollection"],
    InputMCParticles=["MCParticles"],
    OutputPerfectTracks=["PerfectTracks"],
    OutputLevel=INFO,
)

############### Genfit Track Fitter
from Configurables import GenfitTrackFitter

trackFitter = GenfitTrackFitter(
    "GenfitTrackFitter",
    InputTracks=["PerfectTracks"],
    OutputFittedTracks=["FittedTracks"],
    OutputFittedTracksWithFilteredHits=["FittedTracksWithFilteredHits"],
    OutputFittedHits=["FittedHits"],
    OutputLevel=DEBUG,
)

trackFitter.RunSingleEvaluation = True
trackFitter.UseBrems = False

trackFitter.BetaInit = 100
trackFitter.BetaFinal = 0.5
trackFitter.BetaSteps = 15

trackFitter.InitializationType = 1

trackFitter.SkipTrackOrdering = True
trackFitter.ListOfTypesToSkip = [0]

trackFitter.FilterTrackHits = True
trackFitter.FitterType = "DAF"

trackFitter.UseFirstHitAsReference = False

trackFitter.ParticleHypothesisList = [13]

from Configurables import LeastSquaresVertexFitter

vertexFitter = LeastSquaresVertexFitter(
    "LeastSquaresVertexFitter",
    InputTracks=["FittedTracks"],
    OutputVertices=["PrimaryVertices"],
    TrackStateLocation=1,
    MinimumNumberOfTracks=2,
    # fit control (defaults shown explicitly for clarity)
    MaxIterations=100,
    ConvergenceThreshold=1.0e-12,
    SeedStartRadius=-1.0,  # negative => start the seed at the perigee
    # Beam-spot (prior-vertex) constraint.
    #   - Enable it for PROMPT vertices: for back-to-back / collinear tracks
    #     (e.g. Z->mumu at rest) the vertex is otherwise unconstrained along the
    #     dimuon axis and scatters by tens of mm; the constraint pins it and the
    #     resolution collapses to a few um.
    #   - Keep it OFF for DISPLACED / secondary vertices (e.g. Bs->mumu), where
    #     it would wrongly pull the vertex toward the origin.
    UseBeamSpotConstraint=False,
    BeamSpotPosition=[0.0, 0.0, 0.0],  # mm
    BeamSpotSize=[0.010, 0.010, 4.0],  # mm (sigma_x, sigma_y, sigma_z), illustrative FCC-ee-like
    MarkAsPrimaryVertex=True,
    OutputLevel=DEBUG,
)

from Configurables import EventDataSvc

############### Application Manager
import subprocess

ifilename = "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/IDEA/DataAlgFORGEANT.root"
subprocess.run(["wget", "--no-clobber", ifilename])

ApplicationMgr(
    TopAlg=[
        dch_digitizer,
        vtxb_digitizer,
        vtxd_digitizer,
        siwrb_digitizer,
        siwrd_digitizer,
        perfectFinder,
        trackFitter,
        vertexFitter,
    ],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice, EventDataSvc("EventDataSvc"), UniqueIDGenSvc("uidSvc"), RndmGenSvc()],
    StopOnSignal=True,
)

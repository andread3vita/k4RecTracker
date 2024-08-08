from Gaudi.Configuration import INFO, WARNING
from k4FWCore import ApplicationMgr
from Configurables import k4DataSvc, PodioOutput, PodioInput
from Configurables import AuditorSvc, ChronoAuditor
from k4FWCore.parseArgs import parser
from Configurables import HiveSlimEventLoopMgr, HiveWhiteBoard, AvalancheSchedulerSvc

import os

evtslots = 10
threads = 12

whiteboard = HiveWhiteBoard(
    "EventDataSvc",
    EventSlots=evtslots,
    ForceLeaves=True,
)

slimeventloopmgr = HiveSlimEventLoopMgr(
    "HiveSlimEventLoopMgr", SchedulerName="AvalancheSchedulerSvc", OutputLevel=WARNING
)

scheduler = AvalancheSchedulerSvc(ThreadPoolSize=threads, ShowDataFlow=True, OutputLevel=WARNING)

################ parser

parser.add_argument("--inputFile", default="/afs/cern.ch/user/a/adevita/public/workDir/test/dataset/muons_eta0017_mom10/out_sim_edm4hep/out_sim_edm4hep_gun_10.root", help="InputFile")
parser.add_argument("--outputFile", default="output_eff.root", help="OutputFile")
args = parser.parse_args()

################ input & output

from k4FWCore import IOSvc
io_svc = IOSvc("IOSvc")
io_svc.input = args.inputFile
io_svc.output = args.outputFile


################# Simulation setup
# Detector geometry
from Configurables import GeoSvc

geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
print(path_to_detector)
detectors_to_use = ["/afs/cern.ch/user/a/adevita/public/workDir/k4geo/FCCee/IDEA/compact/IDEA_o1_v02/IDEA_o1_v02.xml"]
# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

# digitize vertex hits
from Configurables import VTXdigitizer
import math

innerVertexResolution_x = 0.003  # [mm], assume 5 µm resolution for ARCADIA sensor
innerVertexResolution_y = 0.003  # [mm], assume 5 µm resolution for ARCADIA sensor
innerVertexResolution_t = 1000  # [ns]
outerVertexResolution_x = 0.050 / math.sqrt(
    12
)  # [mm], assume ATLASPix3 sensor with 50 µm pitch
outerVertexResolution_y = 0.150 / math.sqrt(
    12
)  # [mm], assume ATLASPix3 sensor with 150 µm pitch
outerVertexResolution_t = 1000  # [ns]

vtxib_digitizer = VTXdigitizer(
    "VTXIBdigitizer",
    inputSimHits="VTXIBCollection",
    outputDigiHits="VTXIBDigis",
    detectorName="Vertex",
    readoutName="VTXIBCollection",
    xResolution=innerVertexResolution_x,  # mm, r-phi direction
    yResolution=innerVertexResolution_y,  # mm, z direction
    tResolution=innerVertexResolution_t,
    forceHitsOntoSurface=False,
    OutputLevel=INFO,
)

vtxob_digitizer = VTXdigitizer(
    "VTXOBdigitizer",
    inputSimHits="VTXOBCollection",
    outputDigiHits="VTXOBDigis",
    detectorName="Vertex",
    readoutName="VTXOBCollection",
    xResolution=outerVertexResolution_x,  # mm, r-phi direction
    yResolution=outerVertexResolution_y,  # mm, z direction
    tResolution=outerVertexResolution_t,  # ns
    forceHitsOntoSurface=False,
    OutputLevel=INFO,
)

vtxd_digitizer = VTXdigitizer(
    "VTXDdigitizer",
    inputSimHits="VTXDCollection",
    outputDigiHits="VTXDDigis",
    detectorName="Vertex",
    readoutName="VTXDCollection",
    xResolution=outerVertexResolution_x,  # mm, r direction
    yResolution=outerVertexResolution_y,  # mm, phi direction
    tResolution=outerVertexResolution_t,  # ns
    forceHitsOntoSurface=False,
    OutputLevel=INFO,
)

# digitize drift chamber hits
from Configurables import DCHsimpleDigitizerExtendedEdm

dch_digitizer = DCHsimpleDigitizerExtendedEdm(
    "DCHsimpleDigitizerExtendedEdm",
    inputSimHits="CDCHHits",
    outputDigiHits="CDCHDigis",
    outputSimDigiAssociation="CDCHDigisAssociation",
    readoutName="CDCHHits",
    xyResolution=0.1,  # mm
    zResolution=1,  # mm
    OutputLevel=INFO,
)

################ tracker
from Configurables import tracking_func

dch_tracking = tracking_func(
    "GenFitter",
    inputHits_CDC=["CDCHDigis"],
    inputHits_VTXD=["VTXDDigis"],
    inputHits_VTXIB=["VTXIBDigis"],
    inputHits_VTXOB=["VTXOBDigis"],
    inputHits_CDC_sim=["CDCHHits"],
    inputHits_VTXD_sim=["VTXDCollection"],
    inputHits_VTXIB_sim=["VTXIBCollection"],
    inputHits_VTXOB_sim=["VTXOBCollection"],
    outputTracks=["CDCHTracks"],
    outputHits=["Hits"],
    clustering_space=["clustering_space"],
    clustering_space_tracks=["clustering_space_tracks"],
    OutputLevel=INFO,
)

################ efficiency
from Configurables import Efficiency_calc

trackin_eff = Efficiency_calc(
    "Evaluation_efficiency",
    InputCollectionTracks=["CDCHTracks"],
    InputCollectionParticles=["MCParticles"],
    inputHits_CDC_sim=["CDCHHits"],
    inputHits_VTXD_sim=["VTXDCollection"],
    inputHits_VTXIB_sim=["VTXIBCollection"],
    inputHits_VTXOB_sim=["VTXOBCollection"],
    outputTrackingEff_primary=["outputTrackingEff_primary"],
    outputTrackingEff_secondary=["outputTrackingEff_secondary"],
    OutputLevel=INFO,
)

################ Application
# from Configurables import EventDataSvc

# chra = ChronoAuditor()
# audsvc = AuditorSvc()
# audsvc.Auditors = [chra]


# ApplicationMgr(
#     TopAlg=[
#         vtxib_digitizer,
#         vtxob_digitizer,
#         vtxd_digitizer,
#         dch_digitizer,
#         dch_tracking,
#         trackin_eff
#     ],
#     EvtSel="NONE",
#     EvtMax=100,
#     ExtSvc=[geoservice, EventDataSvc("EventDataSvc"), audsvc],
#     StopOnSignal=True,
#     OutputLevel=INFO,
# )

ApplicationMgr(
    TopAlg=[
        vtxib_digitizer,
        vtxob_digitizer,
        vtxd_digitizer,
        dch_digitizer,
        dch_tracking,
        trackin_eff
    ],
    EvtSel="NONE",
    StopOnSignal=True,
    EvtMax=-1,
    ExtSvc=[whiteboard],
    EventLoop=slimeventloopmgr,
    MessageSvcType="InertMessageSvc",
    OutputLevel=INFO,
)
import os
import math 

from Gaudi.Configuration import INFO,DEBUG

from Configurables import EventDataSvc,UniqueIDGenSvc
from Configurables import RndmGenSvc
from Configurables import SimG4SaveTrackerHits
from Configurables import GeoSvc

from k4FWCore import IOSvc,ApplicationMgr
from k4FWCore.parseArgs import parser


################## Parser
parser.add_argument("--inputFile", default="out_sim_edm4hep.root", help="InputFile")
parser.add_argument("--outputFile", default="output_digi_v3_ddplanar.root", help="OutputFile")
parser.add_argument("--detector", default="IDEA_v3_o1", help="IDEA version")
args = parser.parse_args()

# ################## InputOutput
svc = IOSvc("IOSvc")
svc.Input = args.inputFile
svc.Output = args.outputFile

version = args.detector
if version == "IDEA_v3_o1":
    
    ################ Detector geometry
    geoservice = GeoSvc("GeoSvc")
    path_to_detector = os.environ.get("K4GEO", "")
    print(path_to_detector)
    detectors_to_use=['FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml']
    geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
    geoservice.OutputLevel = INFO
    

    from Configurables import DCHdigi_v01
    from Configurables import DDPlanarDigi
    
    idea_vtxb_digitizer_digi = DDPlanarDigi("VTXB_Digi")
    idea_vtxb_digitizer_digi.SubDetectorName = "Vertex"
    idea_vtxb_digitizer_digi.IsStrip = False
    idea_vtxb_digitizer_digi.ResolutionU = [0.003, 0.003, 0.003, 0.003, 0.003, 0.003]
    idea_vtxb_digitizer_digi.ResolutionV = [0.003, 0.003, 0.003, 0.003, 0.003, 0.003]
    idea_vtxb_digitizer_digi.ResolutionT = [1000]
    idea_vtxb_digitizer_digi.SimTrackerHitCollectionName = ["VertexBarrelCollection"]
    idea_vtxb_digitizer_digi.SimTrkHitRelCollection = ["VertexBarrel_Association"]
    idea_vtxb_digitizer_digi.TrackerHitCollectionName = ["VertexBarrelCollection_digi"]
    
    
    idea_vtxd_digitizer = DDPlanarDigi("VTXD_Digi")
    idea_vtxd_digitizer.SubDetectorName = "Vertex"
    idea_vtxd_digitizer.IsStrip = False
    idea_vtxd_digitizer.ResolutionU = [0.003, 0.003, 0.003, 0.003, 0.003, 0.003]
    idea_vtxd_digitizer.ResolutionV = [0.003, 0.003, 0.003, 0.003, 0.003, 0.003]
    idea_vtxd_digitizer.ResolutionT = [1000]
    idea_vtxd_digitizer.SimTrackerHitCollectionName = ["VertexEndcapCollection"]
    idea_vtxd_digitizer.SimTrkHitRelCollection = ["VertexEndcap_Association"]
    idea_vtxd_digitizer.TrackerHitCollectionName = ["VertexEndcapCollection_digi"]
    
    ################ DC digitizer
    DCHdigi = DCHdigi_v01("DCHdigi")
    DCHdigi.DCH_simhits=["DCHCollection"]
    DCHdigi.DCH_name="DCH_v2"
    # DCHdigi.fileDataAlg="DataAlgFORGEANT.root"
    DCHdigi.calculate_dndx=True
    DCHdigi.create_debug_histograms=False
    DCHdigi.zResolution_mm=30
    DCHdigi.xyResolution_mm=0.1
    DCHdigi.OutputLevel=INFO
    
    
    mgr = ApplicationMgr(TopAlg=[idea_vtxb_digitizer_digi,idea_vtxd_digitizer,DCHdigi],
        EvtSel="NONE",
        EvtMax=-1,
        ExtSvc=[geoservice,EventDataSvc("EventDataSvc"),UniqueIDGenSvc("uidSvc"),RndmGenSvc()],
        OutputLevel=INFO,
        )
    
if version == "IDEA_v2_o1":

    # Detector geometry
    from Configurables import GeoSvc
    geoservice = GeoSvc("GeoSvc")
    path_to_detector = os.environ.get("K4GEO", "")
    print(path_to_detector)
    detectors_to_use = ["/afs/cern.ch/user/a/adevita/public/workDir/k4geo/FCCee/IDEA/compact/IDEA_o1_v02/IDEA_o1_v02.xml"]
    geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
    geoservice.OutputLevel = INFO

    # digitize vertex hits
    from Configurables import VTXdigitizer_IDEAv2
    innerVertexResolution_x = 0.003 # [mm], assume 5 µm resolution for ARCADIA sensor
    innerVertexResolution_y = 0.003 # [mm], assume 5 µm resolution for ARCADIA sensor
    innerVertexResolution_t = 1000 # [ns]
    outerVertexResolution_x = 0.050/math.sqrt(12) # [mm], assume ATLASPix3 sensor with 50 µm pitch
    outerVertexResolution_y = 0.150/math.sqrt(12) # [mm], assume ATLASPix3 sensor with 150 µm pitch
    outerVertexResolution_t = 1000 # [ns]

    vtxib_digitizer = VTXdigitizer_IDEAv2("VTXIBdigitizer",
        inputSimHits="VTXIBCollection",
        outputDigiHits="VTXIBDigis",
        outputSimDigiAssociation = "VTXIB_links",
        detectorName = "Vertex",
        readoutName = "VTXIBCollection",
        xResolution = innerVertexResolution_x, # mm, r-phi direction
        yResolution = innerVertexResolution_y, # mm, z direction
        tResolution = innerVertexResolution_t,
        forceHitsOntoSurface = False,
        OutputLevel = INFO
    )

    vtxob_digitizer = VTXdigitizer_IDEAv2("VTXOBdigitizer",
        inputSimHits="VTXOBCollection",
        outputDigiHits="VTXOBDigis",
        outputSimDigiAssociation = "VTXOB_links",
        detectorName = "Vertex",
        readoutName = "VTXOBCollection",
        xResolution = outerVertexResolution_x, # mm, r-phi direction
        yResolution = outerVertexResolution_y, # mm, z direction
        tResolution = outerVertexResolution_t, # ns
        forceHitsOntoSurface = False,
        OutputLevel = INFO
    )

    vtxd_digitizer  = VTXdigitizer_IDEAv2("VTXDdigitizer",
        inputSimHits="VTXDCollection",
        outputDigiHits="VTXDDigis",
        outputSimDigiAssociation = "VTXD_links",
        detectorName = "Vertex",
        readoutName = "VTXDCollection",
        xResolution = outerVertexResolution_x, # mm, r direction
        yResolution = outerVertexResolution_y, # mm, phi direction
        tResolution = outerVertexResolution_t, # ns
        forceHitsOntoSurface = False,
        OutputLevel = INFO
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
    
    mgr = ApplicationMgr(TopAlg=[vtxib_digitizer,vtxob_digitizer,vtxd_digitizer,dch_digitizer],
        EvtSel="NONE",
        EvtMax=-1,
        ExtSvc=[geoservice,EventDataSvc("EventDataSvc"),UniqueIDGenSvc("uidSvc"),RndmGenSvc()],
        OutputLevel=INFO,
        )
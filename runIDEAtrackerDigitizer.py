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
parser.add_argument("--outputFile", default="output_digi.root", help="OutputFile")
parser.add_argument("--detector", default="IDEA_v3_o1", help="IDEA version")
args = parser.parse_args()

# ################## InputOutput
svc = IOSvc("IOSvc")
svc.input = args.inputFile
svc.output = args.outputFile

version = args.detector
if version == "IDEA_v3_o1":
    
    from Configurables import VTXdigitizer
    from Configurables import DCHdigi_v01
    
    ################## Vertex sensor resolutions
    innerVertexResolution_x = 0.003 # [mm], assume 3 µm resolution for ARCADIA sensor
    innerVertexResolution_y = 0.003 # [mm], assume 3 µm resolution for ARCADIA sensor
    innerVertexResolution_t = 1000 # [ns]
    outerVertexResolution_x = 0.050/math.sqrt(12) # [mm], assume ATLASPix3 sensor with 50 µm pitch
    outerVertexResolution_y = 0.150/math.sqrt(12) # [mm], assume ATLASPix3 sensor with 150 µm pitch
    outerVertexResolution_t = 1000 # [ns]
    
    ################ Detector geometry
    geoservice = GeoSvc("GeoSvc")
    path_to_detector = os.environ.get("K4GEO", "")
    print(path_to_detector)
    detectors_to_use=['FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml']
    geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
    geoservice.OutputLevel = INFO
    
    
    ############### Vertex Digitizer
    idea_vtxb_digitizer = VTXdigitizer("VTXBdigitizer",
        inputSimHits = "VertexBarrelCollection",
        outputDigiHits = "VertexBarrelCollection_digi",
        outputSimDigiAssociation = "VertexBarrel_Association",
        detectorName = "Vertex",
        readoutName = "VertexBarrelCollection",
        xResolution = [innerVertexResolution_x, innerVertexResolution_x, innerVertexResolution_x, outerVertexResolution_x, outerVertexResolution_x], # mm, r-phi direction
        yResolution = [innerVertexResolution_y, innerVertexResolution_y, innerVertexResolution_y, outerVertexResolution_y, outerVertexResolution_y], # mm, z direction
        tResolution = [innerVertexResolution_t, innerVertexResolution_t, innerVertexResolution_t, outerVertexResolution_t, outerVertexResolution_t], # ns
        forceHitsOntoSurface = False,
        OutputLevel = INFO
    )

    idea_vtxd_digitizer  = VTXdigitizer("VTXDdigitizer",
        inputSimHits = "VertexEndcapCollection",
        outputDigiHits = "VertexEndcapCollection_digi",
        outputSimDigiAssociation = "VertexEndcap_Association",
        detectorName = "Vertex",
        readoutName = "VertexEndcapCollection",
        xResolution = [outerVertexResolution_x, outerVertexResolution_x, outerVertexResolution_x], # mm, r direction
        yResolution = [outerVertexResolution_y, outerVertexResolution_y, outerVertexResolution_y], # mm, phi direction
        tResolution = [outerVertexResolution_t, outerVertexResolution_t, outerVertexResolution_t], # ns
        forceHitsOntoSurface = False,
        OutputLevel = INFO
    )
    
    ################ DC digitizer
    DCHdigi = DCHdigi_v01("DCHdigi")
    DCHdigi.DCH_simhits=["DCHCollection"]
    DCHdigi.DCH_name="DCH_v2"
    # DCHdigi.fileDataAlg="DataAlgFORGEANT.root"
    DCHdigi.calculate_dndx=True
    DCHdigi.create_debug_histograms=False
    DCHdigi.zResolution_mm=1
    DCHdigi.xyResolution_mm=0.1
    DCHdigi.OutputLevel=INFO
    
    rndm_gen_svc = RndmGenSvc()
    
    mgr = ApplicationMgr(TopAlg=[DCHdigi,idea_vtxb_digitizer,idea_vtxd_digitizer],
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
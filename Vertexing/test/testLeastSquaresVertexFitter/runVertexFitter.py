from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser
from Configurables import EventDataSvc, LeastSquaresVertexFitter
from Gaudi.Configuration import INFO


parser.add_argument("--inputFile", help="Input file produced by the vertex finder")
parser.add_argument("--outputFile", help="Output file containing fitted vertices")
args = parser.parse_args()

io_svc = IOSvc("IOSvc")
io_svc.Input = args.inputFile
io_svc.Output = args.outputFile

vertex_fitter = LeastSquaresVertexFitter(
    "LeastSquaresVertexFitter",
    InputVerticesCandidates=["VertexCandidates"],
    OutputFittedVertices=["FittedVertices"],
    UseBeamSpotConstraint=True,
    BeamSpotPosition=[0.0, 0.0, 0.0],
    BeamSpotSize=[0.01, 0.01, 0.1],
    OutputLevel=INFO,
)

ApplicationMgr(
    TopAlg=[vertex_fitter],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    StopOnSignal=True,
)

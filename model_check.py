import onnx
from onnx import checker

# Carica il modello
model = onnx.load("/afs/cern.ch/user/a/adevita/public/workDir/k4RecTracker/Tracking/model_multivector_1_input.onnx")

# Verifica il modello
checker.check_model(model)
print("Modello ONNX valido")
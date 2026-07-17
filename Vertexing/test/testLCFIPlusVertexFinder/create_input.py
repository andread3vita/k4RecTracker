import edm4hep
import podio
from podio import root_io


def add_track(collection, z0, phi, tan_lambda, omega):
    track = collection.create()
    state = edm4hep.TrackState()
    state.location = edm4hep.TrackState.AtIP
    state.D0 = 0.0
    state.Z0 = z0
    state.phi = phi
    state.tanLambda = tan_lambda
    state.omega = omega
    covariance = [0.0] * 21
    covariance[0] = 1.0e-4
    covariance[2] = 1.0e-6
    covariance[5] = 1.0e-10
    covariance[9] = 1.0e-4
    covariance[14] = 1.0e-6
    state.covMatrix = covariance
    track.addToTrackStates(state)


tracks = edm4hep.TrackCollection()
add_track(tracks, 0.0, -2.0, 0.2, 0.0010)
add_track(tracks, 0.0, 0.2, -0.1, 0.0011)
add_track(tracks, 0.0, 2.2, 0.3, 0.0009)
add_track(tracks, 10.0, -1.0, 0.8, 0.0012)
add_track(tracks, 10.0, 0.7, 0.9, 0.0010)
add_track(tracks, 10.0, 2.5, 0.7, 0.0008)

frame = podio.Frame()
frame.put(tracks, "InputTracks")
writer = root_io.Writer("input.root")
writer.write_frame(frame, "events")
writer.finish()

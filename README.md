# k4RecTracker

This repository hosts Gaudi components related to vertex and tracker reconstruction as well as tracking.

## Dependencies

* ROOT
* PODIO
* EDM4HEP
* Gaudi
* k4FWCore
* DD4HEP

## Installation
```
#go somewhere public (to help support team helping you to debug your code)
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
git clone git@github.com:key4hep/k4RecTracker.git
cd k4RecTracker
git checkout -b trackingPipeline origin/trackingPipeline
k4_local_repo
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install
ctest -R test_TrackFinder -V
ctest -R test_TrackFitter -V
```

## Repository content

* `DCHdigi`: drift chamber digitization (for now, this step produces 'reco' collection)
* `ARCdigi`: ARC digitization (for now, this step produces 'reco' collection)
* `VTXdigi`: vertex detector digitization (for now, this step produces 'reco' collection)
* `Tracking`: tracking algorithms orchestrating [GenFit](https://github.com/GenFit/GenFit)

## References:

These could perhaps be useful for newcomers:
1. [lhcb-98-064 COMP](https://cds.cern.ch/record/691746/files/lhcb-98-064.pdf)
2. [Hello World in the Gaudi Framework](https://lhcb.github.io/DevelopKit/02a-gaudi-helloworld)

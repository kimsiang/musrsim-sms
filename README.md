# musrsim-sms
Geant4 package for musrSim (shanghai muon source dedicated)

# Environment

### Dependencies
Geant4 version `10.7.2`

### Setup environment on INPAC-cluster

```
source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh
```

# Updates

### 2021-11-17 (CC)
Enabled the customizing of crosssection factors on mac steering file, the default value is set to 1.0 if not specified

Example:
```
# set gmumu xsection factor to 1000.0
/musr/command G4EmExtraPhysics SetCrossSecFactor gmumuFactor 1000.0
```

Now can specify a name when launch the job:
```
../musrSim 1003.mac name
```
The output file will be `musr_1003_name.root`. 
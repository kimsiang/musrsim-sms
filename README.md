
# musrsim-sms
Geant4 package for musrSim (shanghai muon source dedicated)

# Tutorial

### Setup environment
This package based on Geant4 version `10.7.2`.

On `INPAC-cluster`, you can setup Geant4 enviroment with:

```
source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh
```
### Download and compile the package

```
git clone https://github.com/kimsiang/musrsim-sms.git
cd musrsim-sms
mkdir build
cd build
make -j4
```
### Create working directory

```
mkdir run
cd run
cp ../../run/1000_Laser.mac ../../run/visVRML.mac .
```
### Start simulation

```
../musrsim 1000_Laser.mac test_run
```
### Run on condor
First copy condor scripts
```
cp ../../run/submit.condor ../../run/run.sh .
```
Then submit 
```
condor_submit submit.condor
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

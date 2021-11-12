musr_X.root:
X=

==============================================================================================
Microbeam - Alex & Elvezio, April 2012

260034 ... ideal parallel beam,  middle collimator 1mm, sample collimator 0.1mm
260033 ... ideal parallel beam,  middle collimator 10mm, sample collimator 1mm
260032 ... ideal parallel beam,  middle collimator 1mm
260031 ... ideal parallel beam,  middle collimator 10mm
260024 ... FOR075_microbeam.DAT, middle collimator 1mm, sample collimator 0.1mm
260023 ... FOR075_microbeam.DAT, middle collimator 10mm, sample collimator 1mm
260022 ... FOR075_microbeam.DAT, middle collimator 1mm
260021 ... FOR075_microbeam.DAT, middle collimator 10mm
260008-9 ... files for tuning the quadrupole triplets using the "optimize.py" python script
260001 ... first trial

==============================================================================================
NEUTRON DETECTOR - Alexey, November 2011:
250001 - first tests for a neutron detector

==============================================================================================
GPS

5018? ... like 5013? but with the new (pyramidal) forward veto - a short one
5017? ... like 5013? but with the new (pyramidal) forward veto - a middle one
5016? ... like 5013? but with the new (pyramidal) forward veto - a long one
5013? ... like 50111 (but larger statistics)
5012? ... like 50110 (but larger statistics)
50111 ... like 50081 but after geometry bug corrected (new geometry)
50110 ... like 50089 but after geometry bug corrected (old geometry)
          CORRECTED BUG IN THE GEOMETRY (FOUND BY ALEX AMATO) 1.2.2012
50108 ... like 50102, old geometry
50107 ... like 50101, old geometry
50102 ... like 50101, 150G
50101 ... based on 50081 but with the longitudinal polarisatoin, 50G
50091-93 ... like 50071-73 but with Geant4.9.5 (13.1.2012)
50090 ... like 50081 but with the old geometry
50089 ... like 50081 but with the old geometry (repetition of the 50021 geometry)
50081-83 ... like 50061-63 but with Geant4.9.5 (13.1.2012)
50080    ... like 50081 but at 0 magnetic field
50071-73 ... like 50051-50053 but with the new GPS design from 9.12.2011
50061-63 ... like 50021-50023 but with the new GPS design from 9.12.2011
50051 ... like 50021 but more focused beam;  some detectors have different numbers!!!
50023 ... like 50021 but B=0.5 T
50022 ... like 50021 but B=0.05 T
50021 ... old GPS geometry, wide beam, B=0.005 T
50001 ... GPS geometry

==============================================================================================
GPD

40171 ... like 40111, but with larger angular coverage by another 2 "backward" detectors.
40161,2  ... like 40121, but the collimator is thicker only in the insert region of 2x4cm.
40152 ... like 40111 but with collimator 2 shifted: -25cm -> -45cm
40151 ... like 40111 but with collimator 2 shifted: -25cm -> -35cm
40141-3 ... like 40111-3 but with both lead collimators thicker
40131-3 ... like 40111-3 but with collimator at 25cm made thicker
40127   ... like 40123 but with m-counter embedded inside another veto detector
40126   ... like 40123 but with a narrow active collimator (to apply it in anticoincidence with backw. counters)
40125   ... like 40123 but with smaller opening in active collimator
40124   ... like 40123 but with larger m-counter
40121-3 ... like 40111-3 but with collimator at 10cm made thicker (3cm -> 5cm)
40113 ... like 40111 but GPDcollimator removed and added scintillators 21 and 22
                         Al plate made of scintillator (veto counter)
          ... i.e. like 40106 (40093)
40112 ... like 40111 but with GPDcollimator made of Pb (instead of Cu)
40111 ... like 40063 but with Geant4.9.4 (9.2.2011), 8000000 events

40106 ... like 40093 but with Geant4.9.4 (31.1.2011) less events than 40105
40105 ... like 40093 but with Geant4.9.4 (28.1.2011)
40104 ... like 40103 but GPDcollimator made of Pb (instead of Cu)
40103 ... like 40101 but scintillators 21 and 22 removed and GPDcollimator returned back
               -> i.e. this run should be very close to 40063 but with thicker GPDcollimator2 collim
40102 ... like 40101 but scintillators 21 and 22 removed (probably not very interesting run)
40101 ... like 40093 but with GPDcollimator2 made thicker (3cm -> 5cm)
40093 ... like 40063 but GPDcollimator removed and added scintillators 21 and 22
                         Al plate made of scintillator (veto counter)
40082  ... like 40063 but GPDcollimator made of Pb
40081  ... like 40063 but thicker GPDcollimator
40061-2 ... SetLargestAcceptableStep=200, beam according to measurement
            done by Alexander with the camera
40051-2 ... test runs without detector to see where the muons stop in sample plane
40041-5 ... significantly improved geometry
40033 ... 0.03T, vertexsigma 5mm, p=100.0 +- 3 MeV
40032 ... 0.03T, vertexsigma 5mm, p=100.0 +- 0 MeV
40031 ... OT,    vertexsigma 5mm, p=100.0 +- 0 MeV
40030 ... OT,    vertexsigma 0mm, p=100.0 +- 0 MeV
40001 ... initial GPD design
===============================================================================================

Dec 2011 tests of Hight Field

31044 ... similar to 31035 but at B= 50 gauss
31043 ... similar to 31035 but at B=500 gauss
31042 ... similar to 31035 but at B= 50 gauss and without the second collimator "collimatorA2"
31041 ... similar to 31035 but at B=500 gauss and without the second collimator "collimatorA2"

May 2011 tests of High Field

31036 ... like 31035 but B=8.0T.
31035 ... 10M, B=0.5T, improved geometry (Ti windows, Kapton window, sample holder)
31006 ... large stat. simulation, B=8.0T
31005 ... large stat. simulation, B=0.5T
31004 ... large stat. simulation, B=50G
31003 ... large stat. simulation, B=0T.
31002 ... modified beam pipe; field rising from 0 to 10 T
31001 ... first tests to see the dependence of muon rate on field


===============================================================================================
Dec 2010 test of High Field


30711-36 ... large statistics simulation like 30605 and other fields
306??  ...  using 2 adjacent detector segments (more statistics than 304?? and a small bug in *.mac corrected)
30583 ... like 30582, but with the goulombRegion commented out
30582 ... like 30581, but with the test size of positron counters
30581 ... like 30580, but with the ALC fieldmap
30580 ... like 8805
30512-3 ... using all 8 segments of detectors
304?? ...   using 2 adjacent detector segments
30220-5 ... like 30200-5 but with longitudinal muon polarisation
30200-5 ... like 30100-5 but with TURTLE beam
30100-5 ... December 2009 test for the High Field instrument done in ALC area - point-like beam ; transverse spin polaris.
==================================================================================================
SPIN ROTATOR STUDY

20810 ... like 20661 but after the bug correction
20809 ... like 20641 but after the bug correction
20808 ... like 20611 but after the bug correction
20807 ... like 20612 but after the bug correction
20806 ... like 20501 but after the bug correction
20805 ... like 20502 but after the bug correction
20804 ... like 20503 but after the bug correction
20803 ... like 20504 but after the bug correction
20802 ... like 20505 but after the bug correction
20801 ... like 20511 but after the bug correction
6.9.2010 - bug found and corrected in the field-map extrapolation of the electric and magnetic fields
           for the spin-rotator (3DEOpera and 3D and symmetryType=1 or 2).

20709 ... like 20708, but with R1IN and R1OU shifted, and beam starting just at spinrot entrance.
20708 ... like 20702, but without the windows
20707 ... like 20702, p=28.0 MeV/c.
20706 ... like 20702, p=27.8 MeV/c.
20705 ... like 20702, p=27.6 MeV/c.
20704 ... like 20702, p=27.2 MeV/c.
20703 ... like 20702, p=27.0 MeV/c.
20702 ... like 20701 but with parallel and monoenergetic beam (27.4 MeV/c)
20701 ... repetition of 20661 but with Geant 4.9.4.p02 (24.8.2011)

20391-4...like 20361-4 but without muon decays
20381 ... like 20351   but without muon decays
20371 ... like 20321   but without muon decays
20361 ... like 20351 but with monoenergetic beam
20351 ... like 20321 but with the matched B and E fields according to the calculation (calculations for 27.4 MeV/c muons)
20331 ... like 20321 but with uniform B and E fields
20321 ... repetition of 20313
20311-6 ... based on 20303 and 20211 using the newest TURTLE (from February 2010)
20303 ... like 20181 but using Geant4.9.3 with appropriately modified steering file
20302 ... repeating 20181 but (i.e. using Geant4.9.p02)
20301 ... like 20181 but using Geant4.9.3
2023? ... like 20181 but using NEW VERSION OF TURTLE
2022? ... like 2021? but with the turtle file generated using NEW VERSION OF TURTLE
2021? ... like 2020? but with larger gap in between the magnets of the spin rotators
2020? ... like 20181 but with correction of "/gun/vertex 0 0 -12032 mm" and different "/gun/turtleInterpretAxes"
2019? ... like 2018? but vith v04 instead of b05c
2018? ... like 2013? but with a little bit larger e01
2016? ... like 20131 but with varying QSE4
2015? ... like 2013?  but with   /gun/turtleInterpretAxes x-y
2014? ... electrodes = e01,  magnetic field = v04,  Turtle = reggiani_Jan2010_z-12032,         varying e01
2013? ... electrodes = e01,  magnetic field = b05c,  Turtle = reggiani_Jan2010_z-12032,        varying e01
2012? ... electrodes = e01,  magnetic field = b05c,  Turtle = pie3_HiField_d05_x30_z-15111
2011? ... electrodes = e01,  magnetic field = v04,  Turtle = pie3_HiField_d05_x30_z-15111
20024 ... like 20023 but with y coordinate interpreted as -y
20023 ... like 20021 but with E and B fields negative
20022 ... like 20021 but with y coordinate interpreted as -y
20021 ... Spin Rotator study - using output from TURTLE generated just in front of the spin rotator

==================================================================================================
HIGH FIELD
14602-12 ... based on 14601 (generated by perl script), but varying the sample (and also cryostat) position.
	   All at 1 T field.
14601 ... based on 14600, but with wide TURTLE beam, and with some modifications in the detector numbering
          (to be closer to my earlier numbering).  Field of 1T.
14600 ... file from J. Rodrigez (ALC-modified flow cryostat)
14593-4 ... like 14591-2 but with point-like beam.
14592 ... like 14591 but 10M events and only interesting events saved 
14591 ... like 145?? but with longitudinal polarisation, field only 50 Gaus (transvers.)
          --> Robert wants to know the long. ampl. of asymm.
145?? ... like 144??, but all events saved and statistics only 1M.
1448? ... like 14452 but very large statistics
144?? ... based on 14060; updated geometry according to the latest knowledge; 8+1 validation counters
143?? ... 28.3.2012 quick (dirthy) test with flow cryostat requested by Robert - using the 0.4 mm Cu shield  
142?? ... 28.3.2012 quick (dirthy) test with flow cryostat requested by Robert - using the 1.0 mm Al shield  
141?? ... like 140?? but with TURTLE with zero momentum bite
14020 ... like 14010 but smaller statistics and storing all events (not just the ones with hit).
14001 ... 15.3.2012; Geant 4.9.5;  based on 12700 and 31043
-----
13850-61 ... like 13750 but with the M-counter as in 12700
13750-61 ... like 13550 but with the old beamline+old collimators (12250) shifted 9mm upstream.
13650-61 ... like 12250 (repetition to see whether I can reproduce 12250).
13550-61 ... like 13150 but with the modified M-counter region and the old beamline+old collimators (12250).
13450-61 ... like 12260 but with the beamline and collimators as in 12950.
13350-61 ... like 13250 but with both collimators as in 12260 (just on different place).
13250-61 ... like 13150 but with the old (first) collimator (as in 12260).
13150-61 ... like 12950 but using 400 mum M-counter (300 mum in previous runs)
13050-61 ... like 12950 but using 400 mum M-counter (300 mum in previous runs)  - bug (M0 overlap)
12950-61 ... like 12700 but using smaller sample holder (diameter=3mm instead of 8mm)
12800-61 ... like 12700 but using larger sample (8mm diameter instead of 5mm diameter)
12700-61 ... NEW DR geometry (based on status of 26.8.2011) (geometry started from runs 12360+30711)
12400-61 ... like 12300 but with TURTLE with zero momentum bite
12300-61 ... based on 8800, but with the different magnet and Turtle files (flow cryostat)
12200-61 ... like 12100 but with Geant4.9.4.p02
12199 ... like 12160 but with Geant4.9.4 (28.1.2011)
12149 ... like 12110 but with Geant4.9.4 (28.1.2011)
12000-11 ... like 10900-10911 but with Geant 4.9.3 and with a new Turtle file (FOR077_reggiani_Feb2010.dat)
11650-9 ... Geant 4.9.3 using Coulomb scattering for muons (16.3.2010)
11631-5 ... different parameters of the processes in Geant 4.9.3.  Problem with generated nr of muons
            due to a bad compilation on merlin00 (effects triggered/generated rate, should
            not cause other problems) (15.3.2010)
116?? ... like 109?? but with Geant 4.9.3, however not modifying the processes (i.e. processes as in 4.9.2.p02 *.mac files)
115?? ... high statistics, old TURTLE Geant 4.9.3
114?? ... like 109?? but with February TURTLE (new turtle executable) and Geant 4.9.3 and 12M of events
113?? ... like 109?? but with January TURTLE (new turtle executable) and Geant 4.9.3
11210 ... like 11110 but with Geant 4.9.3
111?? ... like 109?? but with the NEW TURTLE file from Davide using the NEW TURTLE executable
110?? ... like 109?? but with the newest Cryogenics magnet (from 30.11.2009)
109?? ... like 106?? but with the newest OI magnet (from 3.11.2009)
108?? ... like 106?? but with the longitudinal field (also OK for the pile-up background study)
107?? ... like 106?? but with the field rotation
106?? ... like 102??, but without any field rotation and using new TURTLE file from Davide (pie3_HiField_d05_x30.dat)
105?? ... like 104?? ... (generated beam) but without the rotation of the magnet
10420 ... like 10410 but with magnetic field made positive instead of negative.  10420 and 10410 are expected to be
          identical, if difference found, would be very strange! 
104?? ... like 102?? but with the generated beam (instead of Turtle) with tiltsigma and pitch.
10330-41  ... like 10310-11 but only 5mrad rotation
10310-11  ... like 10210 but with the different directions of the field map rotation (to find out "the worst case scenario")
              (5mrad rotation)
102?? ... like 98?? but with the negative field map (to have the same field map sign as in the old OI map)
101?? ... repetition of 91??
10017-10021 ... like 10007 - 10011 but with the field map rotated in oposit direction (should not change anything)
100?? ... like 97?? but with musrSim and field of the solenoid tilted by 5 mrad
99?? ... like 97?? but with musrSim and field of the solenoid tilted by 20 mrad
98?? ... like 97?? but with musrSim and field of the solenoid tilted by 10 mrad
97?? ... like 85?? but with the new magnet OI_082009
96?? ... like 87?? but with musrSim and field of the solenoid tilted by 20 mrad
95?? ... like 87?? but with musrSim and field of the solenoid tilted by 10 mrad
94?? ... like 85?? but with musrSim and field of the solenoid tilted by 20 mrad
93?? ... like 85?? but with musrSim and field of the solenoid tilted by 40 mrad
92?? ... like 85?? but with musrSim and field of the solenoid tilted by 100 mrad
91?? ... like 85?? but with musrSim and field of the solenoid tilted by 10 mrad
90?? ... like 88?? but using the merged code (musrSim)
8900-11,50-61 ... design of Alexey's dilution cryostat
8800-11,50-61 ... like 7600 but with massiver collimators (small flow cryostat) +  sample holder
8700-11,50-61 ... like 8000 but with massiver collimators (56 mm cryostat)
8600-11,50-61 ... like 84?? but with massiver collimators
8500-11,50-61 ... like 8333 = new suggestion from Robert + more massive collimators
8450-61 ... pileup for 8400_11
8400-11 ... small dilution cryostat (new suggestion from Robert) but using the veto + veto plate a'la Elvezio
8350-61 ... pileup for 8300_11
8300-11 ... small dilution cryostat (new suggestion from Robert)
8250-61 ... pileup for 8200_11
8200-11 ... like 8100-10 but using Elvezio's idea of veto cylinder (no veto behind the sample!)
8150-61 ... pileup for 8100_11
8100-11 ... like 8000-10 but for 48 mm dilution cryostat
8050-61 ... pileup for 8000_11
8000-11 ... 56 mm dilution cryostat + small active collimator size + sample holder rod + 4 cryostat walls
7700-11 ... like 7600-11 but for pileup
7600-11 ... like 7500-11 but with small active collimator size
7590-96 ... different outer radius of the active collimator
7520-40 ... like 7500-11 but with fine step in magnetic field
7500-11 ... like 7400-11 but for the pile-up
7417 ... like 7405 but with the passive collimator made active (in order to investigate why so many muons stop
            in the M0)
7400-11 ... similar to 7180-91 (old version of GEANT = 4.9.1.p03, but not with the changes by Gumplinger) with:
            1) Veto cylinder 2mm thick, with 0.5mm gap around wals, and veto + cylinder combined together
            2) 1 mm Al plate around M0
            3) Larger windows in cryostat flanges, one of them changed from Ti to Mylar
7300-11 ... like 7180-91 but with new GEANT version 4.9.2
7250-61 ... like 7221 (detector radius 50mm) but for fields between 0 and 10 T
7230-41 ... like 7222 (detector radius 60mm) but for fields between 0 and 10 T.
7220-02 ... like 7190 (at 9.99 T) but with larger inner detector radius (40, 50 and 60 mm) (+ modified Active Collim. to fit in pos. counters)
7200-11 ... like 7180-91 but with smaller m0 counter (3 mm instead of 6mm in diameter)
7198    ... Robert's independent test using unpolarised muons (to see that the asymmetry amplitude is zero)
7180-91 ... corrected 7160-71 but with an 8 mm gap between active collimator and the positron detector,
            and only a half veto cylinder (forward half removed), and thinner veto (2mm instead of 2.4mm), 
            and larger m0 counter (6 mm instead of 3 mm in diameter)
7160-71 ... like 7140-50 but with the short cryostat, but with an conflict in geometry (veto cylinder goes through 
            cryostat window).
7140-50 ... like 6720-30 but with the flange on the cryostat 1 outside the cylinder
7120-30 ... like 6720-30 but with a box-like 10T field in the volume of 15x10cm around the sample (z=15cm, R=10cm).
7100-10 ... like 6720-30 but with JASTEC field map
70?? ... like 69?? but for the pileup background (no event reweighting, decay times -1 -1)
6900-90 ... like 6720-30 but for different magnets:  2420d, 2420e, alc_kamil_2D, worsened_A_OI_10T_grid_xz, ISIS
68?? ... like 67?? but for the pileup background (no event reweighting, decay times -1 -1)
6760-70 ... like 6720-30 but with FOR077_pie3_HiField_d02_x15.dat and correct distnace (1.3 meters) from the sample!!!
6740-50 ... like 6720-30 but with FOR077_pie3_HiField_d01_x15.dat and correct distnace (2.5 meters) from the sample!!!
6720-30 ... repeating 6020-30 with musr_Gumplinger_and_Apos executable

6670,6690 ... like 6630 and 6650  but with no limit on muon decay time and all events saved (for the pile-up investigations)
6630,6650 ... like 6530 and 6550 but with new (longer) turtle files from Davide and without any event reweighting
6560, 61, 70 ... like 6540,41,50 but with larger TURTLE input file and musr_Gumplinger_and_Apos executable
6540-50 ... like 6020-30 but for the new turtle input file (focused  beam) from Davide FOR077_pie3_HiField_d02_x15.dat
6520-30 ... like 6020-30 but for the new turtle input file (parallel beam) from Davide FOR077_pie3_HiField_d01_x15.dat 
63?? ... for the pile-up background (shound not be done with event reweighting next time!)
...
...
6020-30 ... (2008.11.  ) ... like 6000-10 but with additional cylindrical veto counter.
6000-10 ... (2008.11.  ) ... cryostat with 30mm outer diameter, 3 walls
5900-10 ... (2008.11.05) ... like 5500-10 but with 60mm diameter detector (dilution-like cryostat).
5800-10 ... (2008.10.24) ... like 5500-10 but with the ISIS magnet
5700-10 ... (2008.10.24) ... like 5500-10 but with the ALC magnet
==== 2008.10.24 ======
FOR070_2008_09_29_XXII.DAT:  RMS of x_ini(z=-1.0m) = 35-40 mm; RMS of y_ini(z=-1.0m) = 27-30 mm;
                             boundary R=70mm;
                             p_ini = 27.74 +/- 0.72 MeV/c;  p_min_allowed=25.2, p_max_allowed=29.8.
==== 2008.10.21 ==== change on Merlin00 - for some reason musr executable crashes now
5600-10 ... like 5500-10 but with the Cryogenics magnet
5540-50 ... (2008.11.05) like 5500-10 but with larger detector and cryostat (inner detector diameter = 40mm).
               The size of veto counter is the same as in 5520-30, larger than in 5500-5510. 
5520-30 ... like 5500-10 but with larger detector, cryostat and veto counter (inner detector diameter = 30mm)
5500-10 ... like 5370-80 but with FOR070_2008_09_29_XXII.DAT turtle file (broader momentum bite)
==== some changes in the code
5370-80 ... like 5350-60 but with FOR070_2008_09_29_XII.DAT turtle file
5350-60 ... 2 milion of events, 0-10T field, FOR070_2008_09_29_XI.DAT turgle file, OI magnet
5220-30 
5200-10
5182 ... 3 milion of events with APDs at 0 T field
5181+3 ... 2 and 1 milion of events with APDs (for radiation studies) at 10 T
5170-80
5161-2 ... for Alexey - test of the solid angle of the detector and APD's in 10T and 0T.
5150-60 ... like 5144 but with FOR070_2008_09_29_I.DAT turtle input file (focused to z=0 at B=0) and B from 0 to 10 T
5145 ... like 5181 but with 50cm long narrow-part beam-pipe and newly positioned conical collimators
5143 ... like 5181 but with 50cm long narrow-part beam-pipe
5142 ... like 5181 (should be exactly the same as 5181)
5141 ... like 5181 but with ALC-like beam pipe
====2008.09.29 ==== 1) implementing G4cons (for conical collimators)
                    2) numberOfGeneratedEvents in G4PrimaryGeneratorAction (i.e. number of the generated events; 
                                  in case of Turtle nr. of events tried); 
====2008.09.25 ==== just a minor change in program - typeOfProcess variable removed from musrParameters.cc and other routines
                    + steering file restructured
5139 ... ISIS high field magnet;
5132-6 ... like 5081 but with field homogeneity worsened 10 times around 10 Tesla; decay times aroud 0.4, 5.4, 10.4, 20.4 and 30.4 microseconds.
5131 ... like 5081 but with 16 detectors instead of 8.
5121-4  ... like 5081 but with sigma(p) of 0, 1, 2 and 4 %.       (5081 = 3%) 
5110-20 ... like 5083 (i.e. with the Robert's TURTLE input file) and field from 0 to 10 T in steps of 1T.
5081, 5091-6 ... repeated once again for the updated version of the program
====2008.09.18.===
             1) Adding the possibility to store the muon polarisation and time
                when first hit occurred in the detector called M0, M1 and M2.
                It is done in a way very similar to "Target" detector.
	     2) Nr. of really generated events and the run number is now
                stored into the "geantParameters" vector.
             3) Definitively removed variables "globalTime" and "fieldValue"
                from the root tree. 
5100-7 ... like 5081, 5091-7 but with 2T field in the sample region instead of 10T (but 10T field around the sample)
replaced 5091-7 ... like 5081 but with the M0-sample distance of 40, 60, 80, 100, 120, 150, 200 mm (instead of 22 mm)
5085 ... like 5081 but with the 5 T field inside the sample (10 T outside) and M0-sample distance only 7mm (instead of 22mm).
5084 ... like 5081 but with the 5 T field inside the sample (10 T outside)
5083 ... like 5081 but with the Robert's TURTLE
5082 ... like 5081 but with my TURTLE
replaced 5081 ... similar to 5002 but with 2.4 mm thick veto detector and broather initial beam
=====2008.09.10==== Some further small improvements
5067 ... like 5061 but using Coulomb scattering in a G4Region for e+ and e-          8.54   11.51
5066 ... like 5061 but with Vvv info commented out.                                  5.04    7.30
5065 ... like 5061 but with goulombRegion commented out.                             8.94   11.60
5064 ... like 5062 but with the decay times around 100400 ns.                        9.25   12.05
5063 ... like 5061 but seting production cuts in G4Region                            5.04    7.39
5062 ... like 5061 but using Coulomb scattering in a G4Region for muons              9.27   11.82
5061 ... should be very close to 5041.mac (after the modifications)                  8.93   11.59
=====2008.08.28==== Processes have to be defined in the steering file now.
                    G4Region implemented ==> Should allow to define coulomb scattering within G4Region
                                         ==> Should allow to set the production cuts in a G4Region
5042-52 ... like 5031 but with the decay times around 20400 40400 60400 80400 100400 150400 200400 250400 300400 400400 500400 ns.
5041 ... like 5031 after the change of code (no difference between 5031 expected)
=====2008.08.20==== 1) Change in event reweighting ===> The criterium is where the muon decays (and not 
                          through which volumes it penetrates)
                    2) Cleaning of the musrStepingAction.cc code
5032-4 ... like 5031 but with the decay times around 1400, 10400 and 100400 ns.
5031 ... like 5002 but after the bug in G4DecayWithSpin was removed  (2008.08.19)
=====2008.08.19==== Bug removed in G4DecayWithSpin causing that there was amplitude dumping as a function of time
5022 ... like 5002 but with 2.4 mm thick veto detector

5020 ... 100k events, decay time fixed to 100405 ns, muons generated at rest inside sample.
5019 ... like 5009 but with the decay times around 100000400 ns
5018 ... like 5009 but with the decay times around 10000400 ns
5017 ... like 5016 but with no muon beam momentum smearing.
5016 ... test run like 5015 but with decay time fixed to 100405 ns.
5015 ... like 5009 but just for testing
5014 ... like 5002 but with M0 counter 200mu thick (instead of 400)
5013 ... like 5002 but with momentum bite 2 times smaller than in 5002
5012 ... like 5002 but with the decay times around 1000400 ns
5010-11 ... like 5002 but with the distance between sample and M0 counter 7 and 12 mm (instead of 22 mm)
5007-9 ... like 5002, but with the decay times around 1400, 10400 and 100400 ns.
5006 ... like 5002, momentum bite of sigma 2 times larger than in 5002
5005 ... like 5002, d=30mm
5004 ... like 5002, d=24mm
5003 ... like 5002, d=14mm
5002 ... like 5001 but high stat.  d = 20mm
5001 ... based on 1351 but with thiner scintillator boxes (low statistics)




3200 ... (27.4.2010) like 3048 but with GEANT 4.9.3 (visualisation for the new publication)
3121 ... (5.10.2009) like 3120 but with only 25 mum thick sampleCellWindow (instead of 50 mum).
3120 ... (5.10.2009) like 3107 but with:  1) Reduced step lenght of muons to 0.01 mm
                                          2) Added M0 counter (300 mum of scintillator + 25 mum of Mylar) 
3112 ... like 3111 but with parallel beam
3111 ... like 3106 but with additional Cu collimator (close to sample, opening of 15mm in diameter)
3108 ... like 3106 but without the cryostat (Alexey was interested what is the effect of cryostat on asymmetry)
3107 ... like 3100 with the Cu sample cell filled with water with TURTLE from Davide
3106 ... like 3100 with the Al sample cell filled with water with TURTLE from Davide
3105 ... like 3100 with the Cu sample plate with TURTLE from Davide
3103 ... like 3100 with the Cu sample cell filled with water
3102 ... like 3100 with the Al sample cell filled with water
3101 ... like 3100 with the Cu sample plate
3100 ... like 3027 but with the new ALC cryostat - drawings from Alexey.  Test setup. (27.8.2009)
---
3068 ... like 3048 but with the new ALC field map form Oleg Kiselev ALC_2D_analytical_oleg.dat
3067 ... like 3048 but with some variables commented out and with the mean momentum of 26MeV
3066 ... like 3048 but with some variables commented out and with the mean momentum of 27MeV
3065 ... like 3048 but with some variables commented out and with the momentum smearing of 3MeV
3060-2 ... like 3051 but 3.33333333xstatistics in each run
3055-9 ... like 3050 but 2xstatistics in each run
3052-4 ... like 3049 but 3.33333333xstatistics in each run
3049-51 ... like 3048 but with 3 different TURTLE input files
3048 ... (7.10.2009) like 3043, small step size also for kapton (not only for target)
3047 ... (9.9.2009) like 3044 but after implementing muTargetMomZ (and smaller statistics)
3046 ... like 3027 but with muon momentum of 26 MeV
3045 ... like 3027 but with muon momentum of 27 MeV
3044 ... like 3027 but with muon momentum smearing of 3 MeV
3043 ... like 3027 but with no shield-like volumes to allow for positron annihilation
---
3042 ... like 3021 but with shifted sample by -2mm
3041 ... like 3021 but with shifted sample by +2mm
3040 ... like 3021 but with shifted ALC field by 5cm (centre of the field is at -5.02 cm)
3039 ... like 3021 but with sample thickness of 1.2 mm
3038 ... like 3021 but with shifted ALC field by 5cm (centre of the field is at +5.02 cm)
3037 ... should be (almost) identical 3021 
3036 ... like 3021 but with tilted ALC field by 0.3 deg
3035 ... like 3021 but with tilted ALC field by 1 deg
3034 ... like 2021 but with sample thickness only 0.8 mm
3033 ... like 2021 but with tilted ALC field by 4 deg
3032 ... like 2021 but with tilted sample by 4 deg
3031 ... like 2021 but with tilted sample by 1 deg
3028 ... like 3027 but with step limit in the sample holder of 0.01 mm (for muons and e+, e-)
3027 ... like 3021 but with step limit in the sample of 0.01 mm (for muons and e+, e-).
3026 ... like 3023 but with G4MuPairProduction & G4MuBremsstrahlung & G4MultipleScattering off
3025 ... like 3023 but with G4MuPairProduction & G4MuBremsstrahlung off
3024 ... like 3023 but with G4MuPairProduction off
3023 ... like 3021 but with monoenergetic muon beam
3022 ... like 3021 but with negative muons (mu-) - the difference in low energy treatment of e+ and e- !!!
3021 ... like 3001 but with sensitive target & all Root varibles stored &
                   field volume has ID number 98 (not 99 as target) & new version of musrSim (after implementing
                   negative muons).
3006 ... like 3002 but with negative muons (mu-) - the difference in low energy treatment of e+ and e- !!!
3005 ... like 3002 but with G4MuPairProduction & G4MuBremsstrahlung & G4MultipleScattering off
3004 ... like 3002 but with G4MuPairProduction & G4MuBremsstrahlung off
3003 ... like 3002 but with G4MuPairProduction off
3002 ... like 3001 but with fewer setps in field (less events), sensitive target & all Root varibles stored &
                   field volume has ID number 98 (not 99 as target)
3001 ... like 2511 but with musrSim (2009-05-27)
2921 ... (2008.10.28) like 2201 but after many many modifications in the program 
=== 2008.10.28 ========================
2901-05 ... like 2511 but with coulomb scattering in target (for both e+ and mu+) and 1/5 statistics

2621-05 ...  like 2601-05 but with TURTLE file ...*KAMIL*...
2601-05 ... mass production: like 2511 but with the TURTLE file   FOR061_2008_04_22.DAT
2591-95 ... mass production: like 2511 but with the point-like beam (/gun/vertexsigma 0 0 0 mm)
2581-85 ... mass production: like 2511 but with the homogeneous beam (/gun/vertexsigma -20 -20 0 mm)
2571-75 ... mass production: like 2511 but with the 28+-3 MeV muon beam
2561-65 ... mass production: like 2511 but with the 26 MeV muon beam
2551-55 ... mass production: like 2511 but with the 27 MeV muon beam
2541-45 ... mass production: like 2511 but with the 3mm sample and 1mm Pb shield
2531-35 ... mass production: like 2511 but with the 3mm sample
2521-25 ... mass production: like 2511 but with the 1mm Pb shield
2511-15 ... mass production: basic distribution, similar to 2501 but with a tiny change in the sample holder
2504 ... like 2501 but with the 3mm sample and 1mm Pb shield
2503 ... like 2501 but with the 3mm sample
2502 ... like 2501 but with the 1mm Pb shield
2501 ... new default; like 2500 but with less variables stored
2500 ... like 2178 (but with the corrected GlobalField implementation):
	    p=28MeV; Kapton; det_particleID added;
            new field from Toni; polarisation =95%; all variables stored

2425 ... like 2422 but with the det_Vvv info added
2422 ... like 2421 but with the det_VrtxtrackID negative if the signal consist of more than just one track 
2421 ... like 2178 but with some new variables of det_??? added

2417 ... like 2181 but with 0.5 deg tilt in y
2416 ... like 2181 but with Al sample holder 4mm thick (instead of 2mm)

2301-2305 ... like 2201-2205 but for the simple tube geometry mimiking the old ALC.
2221-2226 ... like 2211-2216 but with the Cu cryostat added in the gap between the sample and rod.
2216 ... like 2211 but moved by 120 mm outside
2215 ... like 2211 but moved by 100 mm outside
2214 ... like 2211 but moved by 80 mm outside
2213 ... like 2211 but moved by 60 mm outside
2212 ... like 2211 but moved by 40 mm outside
2211 ... like 2202 but with the Cu rod added (simulates the holder of the cryostat)
2205 ... like 2201; Sample = 7mm Cu
2204 ... like 2201; Sample = 5mm Cu
2203 ... like 2201; Sample = 3mm Cu
2202 ... like 2201; Sample = 2mm Cu
2201 ... simple tube geometry mimiking the new ALC, Alexey wanted for the paper; Sample = 1mm Cu


2195-7 ... like 2181 but with the initial beam shifted by 2, 5 and 10 mm in y coordinate
2194 ... like 2181 but with oposite polarisation (100%)
2193 ... like 2181 but with polarisation = 100%
2192 ... like 2181 but with the 3mm sample and 1mm Pb shield
2191 ... like 2181 but with the 3mm sample
2190 ... like 2181 but with the 1mm Pb shield
2189 ... like 2181 but with the sample moved by 1mm forwards
2188 ... like 2181 but with the sample moved by 1mm backwards
2187 ... like 2181 but with the TURTLE beam
2186 ... like 2181 but with the point-like beam
2185 ... like 2181 but with homogeneous beam
2184 ... like 2181 but with momentum width 3MeV
2183 ... like 2181 but with p=26MeV
2182 ... like 2181 but with p=27MeV
2181 ... new default; like 2178 but with less variables stored
2178 ... new default: like 2113 but p=28MeV; Kapton; det_particleID added;
         new field from Toni; polarisation =95%; all variables stored
2177 ... like 2113 but with the random polarisation after removing bug from 2176
2176 ... !! buggy !! like 2113 but with the random polarisation
2175 ... like 2113 but with the polarisation factor -100%
2174 ... like 2113 but with the oposit polarisation
2173 ... like 2113 but with the polarisation factor 0%
2172 ... like 2113 but with the polarisation factor 80%
2171 ... like 2113 but with the polarisation factor 95%
2162-8 ... like 2113 but with the muons generated just in front of the sample 
             at y = 0, 2.5, 5, 7.5, 10, 12.5 and 15 mm
2161 ... like 2113 but with the reuse of the calculated field values
2160 ... special file to calculate the mass of the ALC detector
2159 ... special file to study the effect of the kapton and Ti foils on the muon momentum
2158 ... like 2113 but with new ALC soleniod field map (calculated by Femlab (Tony)) with momentum of 26.00 MeV
2157 ... like 2113 but with new ALC soleniod field map (calculated by Femlab (Tony)) with momentum of 28.00 MeV
2156 ... like 2113 but with new ALC soleniod field map (calculated by Femlab (Tony)) with momentum of 29.79 MeV
2155 ... like 2113 but with quads (as found with Alexey) and parallel beam before quads.
2154 ... like 2113 but with TURTLE (by Kamil  90 cm before solenoid with momentum of 29.79MeV)
2153 ... like 2113 but with TURTLE (by Kamil  90 cm before solenoid)
2152 ... like 2113 but with TURTLE (by Robert 90 cm before solenoid)
2151 ... like 2113 but with the muon mean momentum of 28 MeV/c
2150 ... repete of 2113 (to test how long it runs)
2140-2146 like 2130-2136 but with oposite polarisation
2135 ... like 2130 but with cryostat of length=68mm.
2134 ... like 2130 but with cryostat of length=48mm.
2133 ... like 2130 but with cryostat of length=28mm.
2132 ... like 2130 but with cryostat of length=8mm.
2131 ... like 2130 but with cryostat of length=0mm.
2130 ... based on 2113 but sample holder removed (new ALC basic geometry to be compared with 2131-2135)
2128 ... like 2121 but with flat distributed beam in 20mm radius
2127 ... like 2121 but with just point-like beam
2126 ... like 2113 but with some more variables stored and with the beam offset and tilt that should compensate
2125 ... like 2113 but with some more variables stored and with the beam tilt
2124 ... like 2113 but with some more variables stored and with the beam offset
2123 ... like 2122 but with the oposite numbering of the scintillator detectors and reversed B field.
2122 ... like 2121 but with the oposite numbering of the scintillator detectors
2121 ... like 2116 but with some more variables stored
2111-2117 ... like 2101-2107 but with GEANT 4.9.1 and changes as in 2110
2110 ... like 2100 but with muon decay into two channels (also the radiative one)
2100 ... like 2051 but with GEANT 4.9.1
2101-2107 ... like 2091-2097, but after removing a small bug from ElementField
2097 ... new 3D field map ALC, globalfield, field in log_fieldvolume at 0.5, rotation 2 degree
2096 ... new 3D field map ALC, globalfield, field in log_fieldvolume at 0.5, rotation 1 degree
2095 ... new 3D field map ALC, globalfield, field in log_fieldvolume at 0.5, rotation 0.5 degree
2094 ... new 3D field map ALC, globalfield, field in log_fieldvolume at 0.5, rotation 0.25 degree
2093 ... new 3D field map ALC, globalfield, field in log_fieldvolume at 0.5
2092 ... new 3D field map ALC, globalfield, field in log_target at 0.5
2091 ... new 3D field map ALC, globalfield, field in log_target at 0.01
2090 ... test with the new 3D field map ALC
2084 ... like 2051 but with the signalSeparationTime = 200ns
2083 ... like 2051 but with the signalSeparationTime = 100ns
2082 ... like 2051 but with the signalSeparationTime = 25ns
2081 ... like 2051 but with the signalSeparationTime = 10ns
2078 ... like 2051 but muons generated just in front of the target at point at y=12
2076 ... like 2051 but muons generated just in front of the target at point at y=10
2075 ... like 2051 but muons generated just in front of the target at point at y=8
2074 ... like 2051 but muons generated just in front of the target at point at y=6
2073 ... like 2051 but muons generated just in front of the target at point at y=4
2072 ... like 2051 but muons generated just in front of the target at point at y=2
2071 ... like 2051 but muons generated just in front of the target at point at y=0
2070 ... like 2051 but muons generated just in front of the target (gaussian profile in x&y)
2065 ... Cryogenics, globalfield, field in log_fieldvolume at 0.5, rotation 2 degree
2064 ... Cryogenics, globalfield, field in log_fieldvolume at 0.5, rotation 1 degree
2063 ... Cryogenics, globalfield, field in log_fieldvolume at 0.5, rotation 0.5 degree
2062 ... Cryogenics, globalfield, field in log_fieldvolume at 0.5, rotation 0.25 degree
2061 ... Cryogenics, globalfield, field in log_fieldvolume at 0.5
2060 ... Cryogenics, globalfield, field in log_target at 0.5
2059 ... Cryogenics, globalfield, field in log_target at 0.01
2058 ... link to 2051
2057 ... repete of 2041  (just to be sure that I can reproduce the results of 2041)
2056 ... like 2021 but with the Cryogenics magnetic field
2055 ... like 2051 but with tilt in y of 0.5 degree
2054 ... like 2051 but with tilt in y of 0.25 degree
2053 ... like 2051 but with offset in y of 10 mm
2052 ... like 2051 but with offset in y of 5 mm
2051 ... like 2021 but with some unused variables commented out (to have smaller root output file)
2049 ... like 2021 but with the pitch 0.05
2048 ... like 2046 but with the pitch 0.05 (should slightly overfocuse focuse at 0T)
2047 ... like 2021 but with the tilt and offset (offset should add to the tilt at 0T)
2046 ... like 2021 but with the tilt and offset (offset should compensate the tilt at 0T)
2045 ... like 2021 but with the setlargestacceptablestep to 5 mm (instead of 200)
2044 ... like 2041 but with the setlargestacceptablestep to 5 mm (instead of 200)
2043 ... like 2041 but using the global field instead of magnetic field applied to log_World
2042 ... like 2041 but using the global field instead of magnetic field applied to log_target
2041 ... like 2021 but with the OI magnetic field
2034-6 ... like 2026 but with beam pitch
2031-3 ... like 2021 but with beam pitch
2029 ... like 2021 but with the turtle input (new turtle generated for high field)
2028 ... like 2027 but with the reversed initial polarisation
2027 ... like 2021 but with the muons generated just in front of the sample
2026 ... like 2021 but with the M counter
2025 ... like 2021 but with no sample (investigation of background)
2024 ... like 2023 + 1mm shield of Pb in the forward detector.
2023 ... like 2021 but with 3mm Cu Sample
2022 ... like 2021 + 1mm shield of Pb in the forward detector.
2021 ... 1mm Cu Sample; similar to 2003, but collimators at other position, modified sample holder
2015 ... like 2003 but with the 5mm sample
2014 ... like 2003 but with the far-away colimator off-axis by 10 mm
2013 ... like 2003 but with the near colimator off-axis by 10 mm
2012 ... like 2003 but with the beam tilt of 1.5 deg (no pitch)
2011 ... like 2003 but with the beam tilt of 0.5 deg (no pitch)
2008-10 ... different offsets of the magnetic field: 1cm in x, +/-3 cm in z
2007 ... like 2004 but with beam pitch and tilt
2006 ... like 2004 but with beam pitch and tilt
2005 ... like 2004 but with beam ofset by 10 mm.
2004 ... like 2002 but with thin tube detector holder  (should correspond to alc_07_0523.dat)
2003 ... like 2002 but with large tube detector holder
2002 ... like 2001 but with larger target
2001 ... first test steering for the ALC test in December 2007  (the same steering as in 605.mac)

1400-1415 ... like 1371 but with large beampipe (larger than magnet allows, but the same as was in 1351)
1391 ... like 1371 (to have more statistics)   ... 1371_1391 = 18*300000 events
1381 ... like 1371 but at B=0T
1371 ... like 1351 but with the split-coil Cryogenics field map (and narrower beampipe)
1361 ... like 1351 but at B=0T
1351 ... like 1339 but after some bugs removed; higher statistics;
1341 ... like 1338 but after some bugs removed
1339 ... like 1338 but with reweighting of events with muons in collimators
1338 ... active collimator with the segmented beam-pipe and with the TURTLE input file
1337 ... like 1336 but this time with event reweighting according to the collimator1 (reweighting by factor of 5)
1336 ... similer to 1335 but with broader beam and after modifying the simulation code to allow 
         us the event reweighting (however not used in  this run)
1335 ... like 1332 but with the new version of the active collimator
1331-4 ... like 1321-4 but ???
1321-4 ... like 1311-4 but after migrating from magnetic to electromagnetic field in the F04GlobalField class
1314 ... based on 1313 but with the latest TURTLE input file 
1313 ... based on 1312 but without the old "active collimator"
1312 ... like 1311 but with the G4Beamline-like magnetic field
1311 ... like 1174 (repete on 4 March 2008, after long time)
( old: 1200-13... like 1100-16 but with Cryogenics magnet instead of the OI. )
119? ... like 118? but with narrow beam
118? ... like 117? but with all events saved and Cu used in the active collimator instead of scintillator
117? ... With Coulomb scattering, Ti cryostat, active collimator and the following options
                  M0 counter       M1 counter       Cryostat Cu-wall
     1172           200 mu           -                      -
     1173           200 mu           -                      yes
     1174           400 mu           -                      yes
     1175           200 mu          200 mu                  yes
     1176           400 mu          400 mu                  yes
116? ... like 114? but after removing the bug
115? ... without Coulomb scatt. to compare with 114?
114? ... with Coulomb scattering (but buggy)
1131-4 ... like 1100 but with different l (1131 without M0 counter)
1120-21... just for testing with fine field map inside the sample. 
1100-16... like 1050-1059 but with geant version 4.9.0 and larger target, sample and wider beam
1050-9...like 1035, but with 2000000 of events and different decaytimelimits
1029 ... l=20mm; d=20mm;   scint= 3.63x3.63x3.63; decaytimelimits = 400-420ns, wide beam
1028 ... l=20mm; d=15mm;   scint= 2.9x2.9x2.9; decaytimelimits = 400-420ns, wide beam
1025 ... l=20mm; d=12.2mm; scint= 2.5x2.5x2.5; decaytimelimits = 400-420ns, wide beam
1022 ... l=20mm; d=12.2mm; scint= 2.5x2.5x2.5; decaytimelimits = 10400-10420ns, wide beam
1021 ... l=20mm; d=12.2mm; scint= 2.5x2.5x2.5; decaytimelimits = 100400-100420ns, wide beam
1003-9 ... like 183-9 but with decaytimelimits = 10400-10420ns
1002 ... l=5mm; d=12.2mm,decaytimelimits = 100400-100420ns, wide beam
1001 ... l=5mm; d=12.2mm,decaytimelimits = 10400-10420ns, wide beam
1000 ... like 1002 but with the homogenious beam (not the Gaussian one)
193-9 ... like 183-9 but with focused beam
183-9 ... like 173-9 but with SetLargestAcceptableStep=4 and decaytime between 400 and 420 ns
179 ... like 173, with:  l=20mm; d=20mm; scint= 3.63x3.63x3.63
178 ... like 173, with:  l=20mm; d=15mm; scint= 2.9x2.9x2.9
177 ... like 173, with:  l=20mm; d=20mm; scint= 2.5x2.5x2.5
176 ... like 173, with:  l=20mm; d=15mm; scint= 2.5x2.5x2.5
175 ... like 173, with:  l=20mm; d=12.2mm
174 ... like 173, with:  l=10mm; d=12.2mm
173 ... New simple HiFi geometry suggested on the meeting on 13th June (Alex, Robert, Dirk)
        l=5mm (distance from centre to cryostat backward wall); 
        d=12.2mm (inner distance of opposite scintilators)
172 ... like 170 but with old 3D field map.
171 ... like 170 but all events stored
170 ... like 165 but after changing the way how initial polarisation is defined.
168,9 ... like 165,6 but with old 3D field map.
167 ... like 166 but smaller beam momentum smearing.
166 ... like 165 but at 10 T
165 ... like 156 but with O.I. field map (now field map is read in in polar coordinates)
161 ... like 156 but no magnetic field
160 ... like 156 but no beam momentum and generated position smearing
159 ... like 156 but no beam momentum smearing
158 ... like 156 but without cryostat
156,157 ... like 154 and 155 but with GEANT 4.8.3; corrected titanium windows in the cryostat;  SetLargestAcceptableStep=20;
155  ... like 154 but with decaytime between 1400 and 1440 ns
154  ... tubsbox; otherwise like 152 
153  ... like 141 but with time info of the hits (signalSeparationTime=1ns)
152  ... like 141 but with time info of the hits (signalSeparationTime=0.1ns)
146  ... like 142 but now using G4CoulombScattering for muons
145  ... like 144 but with restricted beam width (sigma=2.5mm)
144  ... like 143 but without the magnetic field
143  ... like 142 but without the cryostats
142  ... like 141 but without beam smearing
141  ... HiFi with tubs; now with smaller range cut and beam size restricted to the beam 
136  ... HiFi;  TF=3T;  (now with smaller range cut (was 1cm before, now 0.1mm and in target 0.01mm)

701 ... like 314 (i.e. the test ALC geometry to our best knowledge) for the updated code/steering.
        (to reproduce 314, sign of the polarisation had to be changed from 1 to -1)   (planned for October 2007 seminar).
615-8 ... reversally polarised muons to 605-8
605-8 ... like 602-4 but with all events (not only those with hits)  (and plus 1.5 mm thickness) ;
602-4 ... new ALC with Cu sample of 1, 2 and 3 mm.
601 ... new ALC design with many details (test run)
561-6 ... like 543,544 but with changing the radius of the iron inner tube for the backward detector (not the forward!)
          (R_Fe_back=30,40,45 mm)   (+34mm for the runs 543,544).
551-8 ... like 543,544 but with changing the thickness of the iron inner tube for the backward detector (not the forward!)
          (Fe_back_thickness=1,3,4,5 mm)  (+2mm for the runs 543,544).
541-6 ... like 525&526 with the correct (inner) R_det=47 & 56 mm & different iron inner tubes of wall
          thickness of 1,2 and 3 mm.
BUG IN SOME PREVIOUS RUNS:  R_det was in fact not the inner radius but the radius in the middle of the scint.
531-538 ... like 521-530 but with l_target=3mm and l_det_backward=80,100,140,160mm (+runs 525&526 for 120mm).
521-530 ... segmented det; l_sep=2*2mm; l_det=120mm; Cu target; l_target= 1,2,3,4,5mm; 
           R_forw=47 mm; R_back=56 mm; different scintilator widths; added Fe tube 2mm;
511-520 ... like 501-510 but with the segmented detectors
501-510 ... l_sep=2*2mm; l_det=120mm; Cu target; l_target= 1,2,3,4,5mm; R_forw=42 mm; R_back=51 mm; added Fe tube 2mm;
491-500 ... l_sep=0; l_det=120mm; Cu target; l_target= 1,2,3,4,5mm; R_forw=45 mm; R_back=54 mm; added Fe tube 2mm;
481-484 ... like 445-6 but with added steel (2mm) or aluminium (4mm) tube 
471-478 ... l_sep=0; l_det=150, 120, 100, 70 mm; Cu target; l_target= 3mm; R_det=47.5 mm  + runs 449-450, 445-446
461-464 ... l_sep=0; l_det=200mm; Cu target; l_target= 3mm; R_back= 47.5 mm; R_forw=42.5, 37.5 mm
451-458 ... l_sep=0; l_det=200mm; Cu target; l_target= 3mm; R_forw= 47.5 mm; R_back=52.5, 57.5, 62.5, 67.5 mm
449-450 ... l_sep=0; l_det=170mm; Cu target; l_target= 3mm; R_det=47.5 mm
441-446 ... l_sep=0; l_det=200mm; Cu target; l_target= 3mm; R_det=27.5, 37.5 & 47.5 mm
431-436 ... l_sep=2*20mm; R_det=27.5mm; l_det=50, 100, 150 mm; Cu target; l_target= 3mm; (+ RUN 425 & 426 for l_det=200mm)
421-428 ... l_sep=2*20mm; R_det=27.5mm; l_det=200mm; Cu target; l_target=0.5, 1, 3, 5 mm;
409-14 ... like 389 and 408 but with the separation of 2*10, 2*20 and 2*30mm.
408   ... like 389 but with reversed polarisation
405-6 ... like previous but with the separation of 2*30mm
403-4 ... like previous but with the separation of 2*20mm
401-2 ... like 382 and 393 but with the separation between detectors of 2*10mm
393 ... like 382 but with reversed polarisation (0,0,1) instead of (0,0,-1)
389-392 ... like 382,386-8 with adding 3mm Cu after the target
386-8 ... like 382 but with increasing R_of_the_backward_detector=27.5, 30, 32.5 and 35
381-5 ... like 375-9 but NO COLLIMATOR !
376-9 ... like 375, but with R_beam_max = 20, 22.5, 25 and 27mm;
375 ... sigma_beam=11mm; R_beam_max=15mm; collimators on
370-1 ... like 364-5 but with the point-like shifted by 27 mm on the x axis
368-9 ... like 364-5 but with the point-like shifted by 18 mm on the x axis
366-7 ... like 364-5 but with the point-like shifted by 9  mm on the x axis
364-5 ... big target, but just a point-like beam
362-3 ... like 358-9 but with beam limited to just r=14mm (like in 354-7) 
360-1 ... like 358-9 but with only water target without brass
358-9 ... like 356-7 but with wider target and beam
354-7 ... like 350-3, but with 11mm sigma of the beam width
350-3 ... like 346-9, but with 5mm sigma of the beam width
348-9 ... like 346-7, but with 200mm long detectors 
346-7 ... water target, 100mm long detectors, point-like beam;
344-5 ... like 332-3 but after the implementation of 2D field possibility (just to check all remains OK)
342-3 ... like 332-3 but with 2mm target with centre at 0mm
340-1 ... like 332-3 but with 2mm target with centre at 1mm
338-9 ... like 332-3 but with beam of spacial spread of generated muons with sigma = 11mm
330-7 ... like 322-9 but with corrected detector diameter to 55mm, no separation gap, point-like beam
          and thinner Cu target.
322-9 ... Simple ALC-like design (just tubes) with the length of scintilators = 50, 100, 200 and 300 mm 
          (study for Alexey and Robert).
320-1 ... like 314-315 but with the target shifted by 10 mm towards the back detector
318-9 ... like 314-315 but with the forward detector shifted by 5 mm in x direction
314-7 ... like 310-3 but with time info of the hits
310-3 ... (NEW DEFINITION!!!) like 290-4 but with 2mm PVC added behind the target.  
312-3 ... (REWRITTEN!!!) like 292-3 but with target thickness of 1 mm instead of 0.8
310-1 ... (REWRITTEN!!!) like 290-1 but with Geant 4.8.3  - no change observed
304-9 ... like 290-1 but with G4CoulombScattering and just for 5 B-field values and 10000 events/fieldValue
300-3 ... like 290-3 but with G4CoulombScattering instead of G4MultipleScattering
292-3 ... Al target
290-1 ... like 260-1 but with solid aluminium tube in forward region added
278-83... like 268-73 but without target
274-5 ... repeated 260-1
272-3 ... like 268-9 but with beam tilt of 0.2 and 0.5 degree
270-1 ... like 268-9 but with vertex shifted by 4 and 9 mm
268-9 ... like 260-1 but with TURTLE BEAM PROFILE
262-7 ... different runs (all?) without target
260-1 ... NEW DEFAULT like 234-5 but with much smaller range cut (was 1cm before, now 0.1mm and in target 0.01mm)
254-5 ... like 234-5 but with penelope processes
252-3 ... like 234-5 but with lowenergy processes
250-1 ... like 234-5 but with forward dectector shifted by 10 mm
240-1 ... like 234-5 but with very simplified geometry
238-9 ... like 234-5 but with random tilt with sigma of 0.5 degree
236-7 ... like 234-5 but with tilt in y of 0.5 degree
234-5 ... ALC; target = 2mm Cu; changed beam properties 
232-3 ... like 222-3 but with the beamwidth of 1mm
230-1 ... like 222-3 but with the momentum smearing of 10%
228-9 ... like 222-3 but with target shifted by 2mm in the forward direction
226-7 ... like 222-3 but without target
224-5 ... like 222-3 but with Cu 2mm target;
222-3 ... ALC; target = 0.8mm Al; detailed geometry description
220, 221 like 210&211 but muon source shifted to ~-1m
218, 219 like 210&211 but without target
216, 217 like 210&211 but with target shifted by 2mm in forward direction
214, 215 like 210&211 but with Cu 2mm target
212, 213 like 210&211 but with 1mm beam width
211 ... ALC; target = Al; B=2.3-5T; only three modules
210 ... ALC; target = Al; B=0-2T; only three modules
208 ... like 205 but with target = Al
207 ... like 204 but with target = Al
206 ... like 204 but with almost no beam width
205 ... ALC; target=Cu; B=2.3-5T
204 ... ALC; target=Cu; B=0-2T
201 ... ALC; design

116 ... like 110 but with the tube target (was box overlaping with detector region in 110)
115 ... 100k like 116 but with TF=3 Tesla;
114 ... 1M; tubsbox;   TF=3 Tesla   , emu20.dat, step=100mm; P=29MeV; cryostat included   31.049206h, 30.420839h
113 ... 1M; tubsbox;   TF=3000 Gauss, emu20.dat, step=100mm; P=29MeV; cryostat included   2.5324222h, 2.4793444h
112 ... 1M; tubsbox;   TF=3000 Gauss, uniform,   step=200mm; P=29MeV; cryostat included   65.873625h, 64.539083h
111 ... 1M; tubsbox;    TF=300 Gauss, uniform,   step=200mm; P=29MeV; cryostat included;  0.97465h,   0.95433889h
110 ... 1M; small tube; TF=300 Gauss, uniform,   step=200mm; P=29MeV; cryostat included;  1.1270083h, 1.0706889h
106 ... like 102 but putting in many more volumes
105 ... 1M; small tube; TF=3 Tesla from emu20.dat; P=29MeV;
104 ... 1M; small tube; TF=3000 Gauss from emu20.dat; P=29MeV;
103 ... 1M; small tube; TF=3000 Gauss; step=200mm; P=29MeV;
102 ... 1M; small tube; TF=300 Gauss; step=200mm; P=29MeV;

1 ... Transverse field 30 Gauss;  strange ratio of Forward/Backward - why 
        not 1?
2 ... Like 1 but two more variable added to the Root Tree.
3 ... Transverse field 10e-6 T.
16... 1M; No field; 
17... 1M; TF=30 Gauss
18... 1M; TF=300 Gauss
19... 1M; TF=30 Gauss; Delta_P/P = 0.05; 
20... 1M; LF=30 Gauss; Delta_P/P = 0.
21... 1M; LF=30 Gauss; P=1MeV and Delta_P/P=0.3
22... 1M; TF=30 Gauss; P=1MeV and Delta_P/P=0.3
23... like 17 (to retest version 17, because 22 seems to be buggy)
24... 0.1M; TF=30 Gauss; P=3MeV and Delta_P/P=0.05
25... 1M; TF=300 Gauss; P=3MeV and Delta_P/P=0.05
26...  debugging
27... 1M; TF=300 Gauss; P=29MeV; Gaussian magnetic field with B_RMS=5mm
33... 0.1M; TF=300 Gaus; P=29MeV; Delta_P/P=0.5 (with down limit)
34... like 33 but changed way of reading the posIniMomX Y Z in the program 
36... 1M; TF=30 Gauss; P=29MeV;         finer range cut = 1mm (was 1cm before)
37... 1M; TF=300 Gauss; P=29MeV;
38... 1M; TF=300 Gauss; Gaussian field (sigma=10mm)
39... 1M; TF=300 Gauss; Gaussian field; P=29MeV; momentum smearing
40... 10M; like 37
41... like 37 but with LF
42... like 38 but gaussian sigma = 5mm
43... like 38 but gaussian sigma = 8mm
44... 1M; P=29MeV; zero field
45... 1M; TF=300 Gauss; P=29MeV
46... 1M; TF=300 Gauss; P=29MeV; momentum smearing
48... 1M; TF=300 Gauss; P=3MeV; momentum smearing with sigma=1.5MeV, cutted
49... 100k like 45, however with new version of musr with the Geant4 muon decay with spin 
51,53,54,55,58 ... like 4? but with the new version of musr with the Geant4 muon decay
70... 1M; TF=300 Gauss; P=29MeV
71... 1M; TF=300 Gauss; P=29MeV; tube geometry
72... 1M; TF=3000 Gauss; P=29MeV; delta P=0.05; tube geom; 
73... 1M; TF=3000 Gauss; P=29MeV; delta P=0.5; tube geom;
74... like 71 but with separation gap of 10degrees between the scintilators.
75... 5M; like 73
76... 5M; like 72,  random seed initialisation = 2
77... 1M; TF=0; P=29MeV; tube geom; random seed init=2
78... like 71 but with GEANT 4.8.2.p01  (and random seed initialsation =2)
79... like 78 but with the removed classes "musrParticleGun.cc" and "musrParticleGunMessenger.cc"
      (was maybe using old geant version ?)
80... 100k; TF=300 Gauss; P=29MeV; tube geometry; GEANT 4.8.2.p01; random seed initialisation = 2
81... 100k; TF=300 Gauss; P=29MeV; tube geometry; GEANT 4.8.2.p01; random seed initialisation not set
82... like 80 but after change of the way how the random seed initialisation is done.
83... 10k; TF=3000 Gauss; P=29MeV; delta P=0.5; tube geom; GEANT 4.8.2.p01 
84... like 80 but with another change how the random seed is initialised
85... 1M; TF=300 Gauss; P=29MeV; tube geometry; GEANT 4.8.2.p01; randomOption = 2 (back to code of 82)
87... like 71 but with the new way how the geometry is read in; GEANT 4.8.2.p01
88... test
90... like 71 but with the geometry read in from the steering file
91... like 90 but using "mm" instead of "cm" for the detector design (was wrong previously!!!)
      - wrong thickness!!!  Was 0.5mm, should be 5mm
92... correceted 91 (thickness of the scintilator 5mm)
98... like 71 but after merging the geometry_XX.txt file with the run_XX.mac file
99... like 92 ----------------------||-------------------------------------------
100... like 99 but with the magnetic field read in from the field map


130001 ... like 102202 (repetition)

120072 ... like 120062 but without G4MultipleScattering, G4MuBremsstrahlung and G4MuPairProduction
12006? ... like 12004? but without SetUserLimits log_??? 0.01
12005[2-4] ... (1.10.2009)  like 12004[2-4] but without /musr/ignore SetUserLimits log_plate  0.01
120049 ... (1.10.2009)  like 120014 but with SetUserLimits log_??? 0.01
12004[2-4] ... (30.9.2009)  like 12003[2-4] but with momentum spread of the initial muons
120038 ... like 120032 but without SetUserLimits log_??? 0.01
12003[2-4] ... like 12002[2-4] but with no initial momentum spred of the muons and with 
           /musr/ignore SetUserLimits log_plate  0.01, and with start of the Cu material at z=0mm.
120024 ... (24.9.2009) like 120022 but Cu sample divided into two pieces (the first is 0.04 mm thick)
120023 ... (24.9.2009) like 120022 but Cu sample divided into two pieces (the first is 0.02 mm thick)
120022 ... (24.9.2009) like 120016 but Cu sample a little bit thicker
120020 ... (18.9.2009) like 120016, but with kapton tape moved just in front of the Cu sample
120019 ... (18.9.2009) like 120018, but with Cu moved just in front of the Cu sample 
120018 ... (18.9.2009) like 120016, but with Cu tape instaed of kapton
120017 ... (10.9.2009) like 120014, but no air and no kapton (i.e. only Cu sample plate)
120016 ... (10.9.2009) like 120014, but no air
120015 ... (10.9.2009) like 120014, but no kapton
120014 ... (10.9.2009) similar to previous, but with kapton tape and air (like in the real ALC experiment)
120010 ... (9.9.2009) test of the penetration depth in the ALC experiment (Cu sample plate) vs. SRIM
120000 ... test with the 76MeV muons stopping in GPD experiment (done for Rustem)

======================================================================================
10433? ... like 10432? but with 10% (?) momentum smearing of initial muons
10432?     like 10431? but without the collimators
104315 ... like 104313 but with mu-
104311-4 ... experiment for David Taqqu with the preassure of 5, 20, 40, 60 mbar
10430? ... modified experiments for David Taqqu
10424? ... large stat, 500 V electric field (oposit to 10423?)
10423? ... large stat, 500 V electric field
104220 ... large statistics
104212 ... based on 104202 but with Geant4.9.5, larger stat. and an additional collimator
104211 ... based on 104201 but with Geant4.9.5, larger stat. and an additional collimator
104208 ... like 104202 but with Geant4.9.5
104207 ... like 104201 but with Geant4.9.5

104202 ... like 104201 but with muons generated at rest in the volume of interest - to be
           used for the detector efficiency calculations.
104201 ... real geometry as obtained from Aldo on 13.12.2011

104061 | 15MeV/c| 61 mum | 
104062 | 15MeV/c| 62 mum | 
104063 | 15MeV/c| 63 mum | 
104064 | 15MeV/c| 64 mum | 
104065 | 15MeV/c| 65 mum | 
104066 | 15MeV/c| 66 mum | 
104067 | 15MeV/c| 67 mum | 

104051 | 12MeV/c| 44 mum | 
104052 | 12MeV/c| 45 mum | 
104053 | 12MeV/c| 46 mum | 
104054 | 12MeV/c| 47 mum | 
104055 | 12MeV/c| 48 mum | 
104056 | 12MeV/c| 49 mum | 

104041 | 10MeV/c| 19 mum | 
104042 | 10MeV/c| 20 mum | 
104043 | 10MeV/c| 21 mum | 
104044 | 10MeV/c| 22 mum | 

104020 | No collim   | 1M events
104021 | 70 mm brass | 1M events
104022 | 70 mm PVC   | 1M events
104023 | 20 mm brass | 1M events

104001  Experiment of Aldo (piE1 in October-November 2011):

======================================================================================
BEAMLINE OSCILLATIONS for Aldo:
103000: for Aldo   ;  beampipe for ALDO
103001: ALC magnet ;  #piE3 TURTLE p=27.5#pm0.52MeV/c;  measured (?) ALC field map
103002: ALC magnet ;  #piE3 TURTLE p=10.0#pm0.19MeV/c;  measured (?) ALC field map
103003: ALC magnet ;  parallel beam p=10.0#pm0.30MeV/c;  measured (?) ALC field map
103004: ALC magnet ;  parallel beam p=10.0#pm0.30MeV/c;  calculated ALC field map
103005: ALC magnet ;  #piE3 TURTLE p=10.0#pm0.19MeV/c shifted;  measured (?) ALC field map
103006: ALC magnet ;  ideal focused beam, p=10.0#pm0.30MeV/c;   measured (?) ALC field map
103007: ALC magnet ;  prallel beam, p=10.0#pm0.30MeV/c; 25#mum scint.;  measured (?) ALC field map
103008: ALC magnet ;  prallel beam, p=10.0#pm0.30MeV/c; 29#mum scint.;  measured (?) ALC field map



======================================================================================

102303 ... like 102302 but with a newer Turtle file from Davide Reggiani
102302 ... like 102301 but with larger statistics
102301 ... (21.3.2012)  based on 102202 but updated to the latest design.
----
102202 ... like 102002 but using FOR070_2008_10_17_XXII.DAT turtle input file.
102102 ... (2008.10.20.) similar to 102002 but using FOR070_2008_10_17_XXII.DAT and including 
		new method to calculate the field integral along the muon path.
1020?? ... beam envelope using  TURTLE input file "FOR070_2008_09_29_XII.DAT" 
                ... like FOR070_2008_09_29_II.DAT but with delta p/p ~ 3% instead of ~1.5%
1020?? ... beam envelope using  TURTLE input file "FOR070_2008_09_29_XI.DAT" 
                ... like FOR070_2008_09_29_I.DAT but with delta p/p ~ 3% instead of ~1.5%
102002 ... beam envelope using  TURTLE input file "FOR070_2008_09_29_II.DAT"
                ... based on the initial Robert's turtle file, but focused
                    better to z=-20cm at zero field (to have a different
                    TURTLE file to see a possible influence of the Turtle files)
102001 ... beam envelope using  TURTLE input file "FOR070_2008_09_29_I.DAT"  
                ... very similar to the initial Robert's turtle file, but focused
                    better to z=0 at zero field



101061-2 ... like 101051-2 but for the 2m long beampipe (i.e. for the ALC design)
101053 ... like 101052 but with the better focuse to the centre of the solenoid
101052 ... like 101051 but 10 times more events
101051 ... like 101021 but with some more save volumes also before the quadrupole triplet
101033 ... like 101032 but with the r=15mm collimator at z=-10cm
101032 ... like 101031 but with the r=20mm collimator at z=-80cm
101031 ... real quads focused to -20cm, TURTLE 60
101023 ... like 101022 but with the r=15mm collimator at z=-10cm
101022 ... like 101021 but with the r=20mm collimator at z=-80cm
101021 ... real quads focused to +10cm, TURTLE 60
101019 ... like 101013 but with the restricted beam radius at the quadrupoles (should be directly comparable to 101014)
101018 ... real quads (as found out with Alexey) with field gradients adujsted to focuse parallel beam to the centre
101014 ... Turtle61
101013 ... quadrupoles as in the TURTLE input file from Robert,  the TURTLE output generated before the triplet (Turtle60)
101012 ... like 101011 but with the TURTLE beam before the quadrupole
101011 ... ALC with parallel beam before the quadrupoles, quadrupoles as measured with Alexey
101002 ... like 100392, but with the global field (implemented from F04 example from Gumplinger)
100411-24 ... like 100311-24 but for Cryogenics
100392 ... only "save" volumes, new TURTLE from Robert as is
100391 ... like 100387 but with less steps in the magnetic field
100390 ... only "save" volumes, gaussian parallel beam, tilt of the beam in y (corresponds to 5mm in the zero field)
100389 ... only "save" volumes, gaussian parallel beam, offset of the beam in y (5mm)
100388 ... only "save" volumes, gaussian parallel beam, smaller sigma in y than in x
100387 ... only "save" volumes, gaussian parallel beam
100386 ... only "save" volumes, focused beam
100385 ... only "save" volumes, TURTLE focused more or less to the centre.
100384 ... only "save" volumes, TURTLE as is
100381-3 ... using the new possibility of the "save" volumes
100331
100321-9 ... like 100311-4 but with small radius of M0 counter and of the Sample
100319 ... OI magnet; FWHM(p)=3%; pitch such that beam focused 10cm after the sample at 0 field
100318 ... OI magnet; FWHM(p)=3%; pitch such that beam focused 50cm after the sample at 0 field
100317 ... OI magnet; FWHM(p)=3%; TURTLE input file; pitch
100316 ... OI magnet; FWHM(p)=3%; TURTLE input file; no pitch
100315 ... OI magnet; FWHM(p)=3%; pitch such that beam focused 20cm before the sample at 0 field
100314 ... OI magnet; FWHM(p)=3%; pitch such that beam focused 5cm before the sample at 0 field
100313 ... OI magnet; FWHM(p)=3%; pitch such that beam focused to the sample at 0 field
100312 ... OI magnet; FWHM(p)=8%; no pitch
100311 ... OI magnet; FWHM(p)=3%; no pitch
100301 ... like 100161 but with the beam focused at 0 Field to the centre of the sample, i.e. to (0,0,0).
100221-4 ... like 100121-4 but with momentum spread FWHM=8%
100171-4 ... M0 counter = 0.4 mm , momentum spread FWHM=3%; R=5, 10, 15 & 20 mm;  l=20mm with Cryogenics magnet
100161   ... M0 counter = 0.4 mm , momentum spread FWHM=3%; Gaussian beam;  l=20mm, 
100151-4 ... M0 counter = 0.4 mm , momentum spread FWHM=3%; R=5, 10, 15 & 20 mm;  l=5mm
100141-4 ... M0 counter = 0.4 mm , momentum spread FWHM=3%; R=5, 10, 15 & 20 mm;  l=20mm
100131-4 ... M0 counter = 0.4 mm , no momentum spread, R=5, 10, 15 & 20 mm; l=20mm
100121-4 ... no M0 counter, momentum spread FWHM=3%; R=5, 10, 15 & 20 mm; l=20mm
100111-4 ... no M0 counter, no momentum spread, R=5, 10, 15 & 20 mm; l=20mm
100101-2 ... 0.4 mm M0 counter, momentum spread, R=10&20mm; l=20mm
100002 ... like 100001 but muons start at -50cm.
100001 ... muon beam properties in the target: OI magnet, muons start at z =~ -130cm.

-------------------------------------------------------------------------------------
24632? ... like 24628 but with 100% diffusion scattering.
24631? ... like 24630? but air properties commented out - 10 times higher number of photons!
24630? ... to investigate why 24628? failed (World replaced by "vzduch") - no difference found
24629? ... like 24624? but with "polished unified" and diffusive scat. near edges.
24628? ... like 24621? but with 10% diffusive scat.
24627? ... like 24622? but 5x5 mm.
24626? ... like 24621? but 5x5 mm.
24625? ... like 24620? but 5x5 mm.
24624? ... like 24620? but with the "groundair LUT" near edges (0.1 mm)
24623? ... like 24620? but with the bevelled edges
24622? ... like 24620? but with "polishedair LUT" surface
24621? ... like 24620? but with "polished unified" surface
24620? ... B418, SPTR=230, ideal box, ideal surface  (20.7.2012)

24615? ... like 24604? but with polishedteflonair LUT
24614? ... like 24604? but with the beam off-center (Alexey wanted to see if it plays a role)
24613? ... like 24605 but for thinner scintillator bar
24610-2? ... like 24607?-9? but for the scintillator BC418
24607-9? ... like 24604?-6? but for attenuation length of 80 mm
24606? ... like 24603?, SPTR=230; tau1=0.35; tau2=1.6
24605? ... like 24603?, SPTR=120; tau1=0.35; tau2=1.6
24604? ... like 24603?, SPTR=  0; tau1=0.35; tau2=1.6
24603? ... study on the long scintillator bar: SPTR=0; tau1=0; tau2=0

2455?? ... like 2448?? but with 60x40 APD cells
2454?? ... like 2448?? but with just 1 APD cell
2453?? ... like 2448?? but with 6x4 APD cells
2452?? ... like 2447?? but with just 1 APD cell
2451?? ... like 2447?? but with 6x4 APD cells
2450?? ... like 2443?? but with even small detector (1000 x smaller than normal size).
2448?? ... like 2424?? but with 100% reflectivity
2447?? ... like 2437?? but with 100% reflectivity
2446?? ... like 2436?? but with 100% reflectivity
2445?? ... like 2435?? but with 100% reflectivity
2444?? ... like 2434?? but with 100% reflectivity
2443?? ... like 2442?? but with extremely small detector
2442?? ... like 2434?? but with rise times of 0
2441?? ... based on 2434?? for rise times of 0, 0.063, 0.35, 1.0
                          and decay times of 1, 1.32, 1.6, 5., 15 and 43.8
2440?? ... like 2422 but with the random time jitter of 120 ps
243941  ... like 243841 but for LaBr3
243841  ... like 243810 but with the RESOLUTIONSCALE=2.5  (to have amplitude smearing 
            for the comparison with Seifert)
243[8-9] ... like 2427?? and 2428?? but random time jitter of 120 ps
243[4-7] ... like 2420 2423 2424 2425 but random time jitter of 120 ps
243[0-3] ... like 2420 2423 2424 2425 but with rise/decay times of 0.063 and 1.32
2428?? ... like 2424??, but for LaBr3 scintillator (Seifert)
2427?? ... like 2424??, but for LYSO scintillator (Seifert)
2426?? ... like 2422 but with the shorter rise and decay times                                 
2425?? ... 2424?? + 2400 APD cells
2424?? ... 2423?? + random time jitter of 230 ps (i.e. should replace previous 224950)
2423?? ... 2420?? + cross talk of 0.11           (i.e. should replace previous 224900)
2422?? ... based on 2218?? - i.e. experimental setup modified according to 2420??, higher light yield
2421?? ... based on 2218?? - i.e. experimental setup modified according to 2420??
2420?? ... based on 2248?? but with GEANT 4.9.5 (10.1.2012) run on Merlin4,
           rise/decay times of 0.35 and 1.6 respectively, and CFD delay of 1.5 ns.
-------------------------------------------------------------------------------------
230001-5 ... test for Alexey: muons passing through scintillator of different thicknesses
-------------------------------------------------------------------------------------

2256?? ...  like 2248?? but run on Merlin4 (9.1.2012)
2251??-  ... like 2249?? + finite nr. of cells effect (60x60 cells)
225001-  ... like 224801 + finite nr. of cells effect (6x6 cells)
224901-  ... like 224801 + cross talk
224801-  ... repetition of 223301 after modifications of the code
 --------- end of October 2011: Changes in the programm (cross-talk added, cell jitter changed)

224201-62 ... like 2239?? but with the finite nr. of cells effect (1x1 cells)
224101-62 ... like 2239?? but with the finite nr. of cells effect (6x6 cells)
224001-62 ... like 2239?? but with the finite nr. of cells effect (60x60 cells)
223901-62 ... wide electron beam covering the full scintillator area
223301-862 ... like the previous runs, but after correcting the "dielectric_LUT" bug.
           ... statistics increased.
MISTAKE in the previous runs - surface of the scintillator was wrongly defined
        ("dielectric_dielectric" instead of "dielectric_LUT")
222801-62 ... like 222801-62 but with statistical fluctuations for photon creation
222701-62 ... like 222301-62 but with statistical fluctuations for ionis. losses
222601-62 ... like 222301-62 but rise time=0.35, decay time=1.32.
222501-62 ... like 222301-62 but rise time=0.,   decay time=1.32.
222401-62 ... like 222301-62 but rise time=0.35, decay time=1.6.
222351-62 ... like 222301-12 but sigma(APD)=100ps.
222301-12 ... light transport in a single counter; ionisation losses and photon creation
              without statistical fluctuations.
              Rise time=0.063, decay time=1.32, sigma(APD)=0ps.

222021-2
222001 ... simplified simulation of 221131, where there is just one counter, and
           2 MeV electron is generated at its entrance.

-------------------------------------------------------------------------------------

221951-67 ... like 221921-37 but with a different APD responce function
              (APDpulseShapeFile_2011_Fitted.txt)
221921-37 ... like 221801 but with CFD information
221801 ...  4.11.2011 - based on 225101 and 221701
221401 ....13.10.2011 - based on 223901 and 221131

221301 ... like 221131 but smaller statistic and without plastic holders
221211 ... like 221131 but smaller statistic and scint rise time = 1.0ns
221201 ... like 221131 but smaller statistic and scint rise time = 0.75ns
221191 ... like 221131 but smaller statistic and scint rise time = 0.5ns
221181 ... like 221131 but smaller statistic and scint rise time = 0.35ns
221171 ... like 221131 but smaller statistic and scint rise time = 0.2ns
221161 ... like 221131 but smaller statistic and scint rise time = 0.1ns
221151 ... like 221131 but smaller statistic and scint rise time = 0.03ns
221141 ... like 221131 but smaller statistic and scint rise time = 0.0ns

221134 ... like 221131 but using 0.1 ps of random time scattering (instead of 100ps)
221133 ... like 221131 but using APDpulseShapeFile_2011_1e_H.txt
221132 ... like 221131 but delay time back to 1.5 ns
221131 ... like 221101_8 but after implementing CFD array & delay time=0.5 ns

221123 ... like 221101-8 but using CFD delay of 0.8 ns
221122 ... like 221101-8 but using CFD delay of 0.7 ns
221121 ... like 221101-8 but using CFD delay of 0.6 ns
221120 ... like 221101-8 but using CFD delay of 0.4 ns
221119 ... like 221101-8 but using CFD delay of 0.3 ns

221118 ... like 221101-8 but using CFD delay of 4.0 ns
221117 ... like 221101-8 but using CFD delay of 3.0 ns
221116 ... like 221101-8 but using CFD delay of 2.0 ns
221115 ... like 221101-8 but using CFD delay of 1.0 ns
221114 ... like 221101-8 but using CFD delay of 0.5 ns
221113 ... like 221101-8 but using CFD delay of 0.2 ns 
221112 ... like 221101-8 but using APDpulseShapeFile_2011_1e_H.txt (steeper pulse rise)
221111 ... like 221101-8 but using 0.1 ps of random time scattering (instead of 100ps)
221101-8 ... new big simulation using APDpulseShapeFile_2011_1e_Z.txt
221086  ... like  221084 but after a change of stacking order for optical photons.
221083,84 ... with a newer APD response function
221061-4 ... like 221051-4 but after changing the CFD signal shaping and introducing some other
         ... definitions of odet_timeC and odet_timeD.
221051-4 ... implementing also rise time of the photons in the scintillator, but there were some bugs.
221041-4 ... like 221032 but after correcting bug in mean time + different random option in some runs
221033 ... like 221032 but APD readout from both sides
221032 ... like 221031 but /musr/command OPSA CFD -0.2 1.5 -1.8
221031 ... like 221011 but /musr/command OPSA CFD -0.3 2. -1.8 + plastic holders, huge statistics

221014 ... like 221011 but /musr/command OPSA CFD -0.3 0.5 -1.8
221013 ... like 221011 but /musr/command OPSA CFD -0.1 2. -1.8
221012 ... like 221011 but /musr/command OPSA CFD -0.2 2. -1.8
221011 ... like 104 (i.e. Sr source, C1 and C2 positron counters) including optical photons 
       ... /musr/command OPSA CFD -0.3 2. -1.8
221000 ... again some Sr Tests, after long time. Now with Geant 4.9.4. (Feb 2011).

220051-9 ... different fields
220001 ... muSR time asymmetry measurement in the ALC magent - test for Alexey

210013 ... TOF measurement of the G-APDs in the ALC magnet using C1 and C2 positron counters;  electrons as beam; B=4.8T - Coulomb sc.
210012 ... TOF measurement of the G-APDs in the ALC magnet using C1 and C2 positron counters;  electrons as beam; B=0T
210011 ... TOF measurement of the G-APDs in the ALC magnet using C1 and C2 positron counters;  electrons as beam; B=4.8T
210002 ... TOF measurement of the G-APDs in the ALC magnet using C1 and C2     muon counters;      muons as beam; B=0T
210001 ... TOF measurement of the G-APDs in the ALC magnet using C1 and C2     muon counters;      muons as beam; B=4.8T

200037 ... Sr test with just C1 and C2 of the size 3x3x2mm
200031-35 ... Sr test with just C1 and C2 of the size 4x4.6x2 mm
200021 ... muons stopped in the centre of the sphere; C1=5mm, C3=20mm
200019 ... like 200013 but with only C1 and C2 (to test new idea of Alexey)
200018 ... Sr experiment with C1 and C2 moved very close to C3
200013-7... Sr experiment in the set-up as in the case of the Alexey's experiment (C1=2mm, C2=2mm, C3=60mm)
200011 ... Strontium source emiting beta electrons towards the sphere detectors (C1=2mm, C2=2mm, C3=20mm)
200008 ... like 200007, but with 5 MeV positrons as a source
200007 ... like 200004, but with 30 MeV positrons as a source
200006 ... Strontium source emiting beta electrons towards the sphere detectors (the first one is 0.2 mm thick + 0.3mm Cu)
200005 ... Strontium source emiting beta electrons towards the sphere detectors (the first one is 0.2 mm thick)
200004 ... Strontium source emiting beta electrons towards the sphere detectors (the first one is 5 mm thick)
200003 ... Strontium source emiting beta electrons towards the sphere detectors (the first one is 2 mm thick)


PICTURES:

HiFi_01.eps:   - longitudinal cross-section through the HiFi detector.
  visual.mac  with  
/vis/viewer/zoom                   20
/vis/viewer/set/sectionPlane on 0 0 0 cm -1 0 0
  and setting z coordinate to 15 (mm?)


Field Integral for OI for Robert (mean value =1.45 Tm):
root [1] TFile* f=new TFile("data/musr_102151.root")
root [2] t1->Draw("BzIntegral","(BzIntegral>-2)&&(muDecayDetID==1201)")


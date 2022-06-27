#include "musrDetectorConstruction.hh"
#include "musrPhysicsList.hh"
#include "FTFP_BERT.hh"
#include "musrPrimaryGeneratorAction.hh"
#include "musrRunAction.hh"
#include "musrEventAction.hh"
#include "musrStackingAction.hh"
#include "musrSteppingAction.hh"
#include "musrSteppingVerbose.hh"
#include <string>

#include <fstream>
#include <sstream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#ifdef G4UI_USE_TCSH
#include "G4UItcsh.hh"  //JSL: should be commented on windows ?
#endif

#ifndef G4VIS_USE
#define G4VIS_USE
#endif

//#include <TApplication.h>
//#include <TSystem.h>


// The following two lines are needed to cope with the problem of 
// "Error in <TPluginManager::FindHandler>: Cannot find plugin handler for TVirtualStreamerInfo! 
//               Does $ROOTSYS/etc/plugins/TVirtualStreamerInfo exist?"
#include "TROOT.h"
#include "TPluginManager.h"

#include "Randomize.hh"

// #include <X11/Xlib.h>  //JSL

#ifdef G4VIS_USE
// #include "musrVisManager.hh"
#include "G4VisExecutive.hh"
#include "G4TrajectoryDrawByCharge.hh"    // TS Trajectory drawing by ID or charge
#endif

#include "musrRootOutput.hh"
#include "musrParameters.hh"
#include "musrErrorMessage.hh"
//#include "F04GlobalField.hh"

int main(int argc,char** argv) {
    //  XInitThreads();
    G4cout<<"\n\n*************************************************************"<<G4endl;
    G4cout<<" musrSim version 1.0.5 for Geant4.10.3, released on 20 Mar 2017"<<G4endl;
    G4cout<<"      WWW:  https://www.psi.ch/lmu/geant4-simulations"<<G4endl;
    // choose the Random engine
    //  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);   // the /musr/run/randomOption 2 does not work with RanecuEngine
    CLHEP::HepRandom::setTheEngine(new CLHEP::HepJamesRandom);

    //my Verbose output class
    G4VSteppingVerbose::SetInstance(new musrSteppingVerbose);

    // Run manager
    G4RunManager * runManager = new G4RunManager;

    // Create class "myParameters", which is a collection of many different parameters
    G4String steeringFileName=argv[1];
    musrParameters* myParameters = new musrParameters(steeringFileName);

    // Create class "musrErrorMessage"
    musrErrorMessage* myErrorMessage = new musrErrorMessage();

    //  TApplication* myapp=new TApplication("myapp",0,0);

    // Create Root class for storing the oMuMu Filterutput of the Geant simulation
    std::string name = "";
    int random_seed_offset = 0;

    if(argc > 2) {
        name = std::string(argv[2]);
        std::cout << "\nmusrSim.cc: Set output ROOT file: " << "data/musrSim_" << atoi(argv[1]) << "_" << name << ".root\n" << std::endl;
    }
    if(argc > 3) random_seed_offset = atoi(argv[3]);
    // Read from macro to customize output ROOT file name
    // Only if the name not specified in argv[2]
    std::ifstream fin(steeringFileName.c_str());
    std::string line;
    while (!(argc>2) && std::getline(fin, line)){
        std::istringstream ss(line);
        std::string cmd1, cmd2, cmd3;
        if (ss >> cmd1 >> cmd2 >> cmd3){
            if (!cmd1.compare("/musr/command") && !cmd2.compare("SetOutputFileName")){
                if (!cmd3.compare("DEFAULT")) {
                    std::cout << "\nmusrSim.cc: Set default output ROOT file: " <<"data/musrSim_" << atoi(argv[1]) << ".root\n" << std::endl;
                }
                else {
                    name = cmd3;
                    std::cout << "\nmusrSim.cc: Set output ROOT file: " << "data/musrSim_" << atoi(argv[1]) << "_" << cmd3 << ".root\n" << std::endl;
                }
            }
            // Set random seed in macro file
            if (!(argc >3) && !cmd1.compare("/musr/command") && !cmd2.compare("SetRndSeedOffset")){
                random_seed_offset = atoi(cmd3.c_str());
                std::cout << "musrSim.cc: Set random seed offset: " << random_seed_offset << std::endl;
            }
        }
    }

    musrRootOutput* myRootOutput = new musrRootOutput(name);

// The following command is needed to cope with the problem of 
// "Error in <TPluginManager::FindHandler>: Cannot find plugin handler for TVirtualStreamerInfo! 
//               Does $ROOTSYS/etc/plugins/TVirtualStreamerInfo exist?"
//  /* magic line from Rene - for future reference! */
    gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                          "*",
                                          "TStreamerInfo",
                                          "RIO",
                                          "TStreamerInfo()");


    // UserInitialization classes (mandatory)
    //  musrDetectorConstruction* musrdetector = new musrDetectorConstruction;
    if (argc>1) {
        G4int myRunNr=atoi(argv[1]);       // Get the run number from the name of the
        // parameter file, if it starts with a number.
        if (myRunNr>0)  {runManager->SetRunIDCounter(myRunNr);}
        //    musrdetector->SetInputParameterFileName(argv[1]);
    }
    musrDetectorConstruction* musrdetector = new musrDetectorConstruction(steeringFileName,random_seed_offset);
    runManager->SetUserInitialization(musrdetector);
    //runManager->SetUserInitialization(new musrPhysicsList);
    runManager->SetUserInitialization(new FTFP_BERT(1,steeringFileName));

#ifdef G4VIS_USE
    // Visualization, if you choose to have it!
    //  G4VisManager* visManager = new musrVisManager;
    G4VisManager* visManager = new G4VisExecutive;    // TS Trajectory drawing by ID or charge
    visManager->Initialize();
#endif

    // UserAction classes
    runManager->SetUserAction(new musrPrimaryGeneratorAction(musrdetector));
    runManager->SetUserAction(new musrRunAction);
    runManager->SetUserAction(new musrEventAction);
    // Initiate musrStackingAction only if optical photons are required (otherwise not needed)
    if (musrParameters::boolG4OpticalPhotons) runManager->SetUserAction(new musrStackingAction);
    runManager->SetUserAction(new musrSteppingAction);

    //Initialize G4 kernel
    runManager->Initialize();

    //get the pointer to the User Interface manager
    G4UImanager * UI = G4UImanager::GetUIpointer();

    if(argc==1)
        // Define (G)UI terminal for interactive mode
    {
        // G4UIterminal is a (dumb) terminal.
        G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
        session = new G4UIterminal(new G4UItcsh);
#else
        session = new G4UIterminal();
#endif

        UI->ApplyCommand("/control/execute vis.mac");
        session->SessionStart();
        delete session;
    }
    else
        // Batch mode
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UI->ApplyCommand(command+fileName);
        if (argc>2) {
            G4String SecondArgument = argv[2];
            if (SecondArgument=="idle") {
                G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
                session = new G4UIterminal(new G4UItcsh);
#else
                session = new G4UIterminal();
#endif
                G4cout<<"Go to idle state now:"<<G4endl;
                session->SessionStart();
                delete session;
            }
        }
    }

    //  myapp->Run(kTRUE);

#ifdef G4VIS_USE
    delete visManager;
#endif
    delete myRootOutput;
    delete myErrorMessage;
    delete myParameters;
    // cks  runManager->SetVerboseLevel(2);   // This line can help debug crashes during the runManager delete
    delete runManager;
    //  F04GlobalField* myGlobalField = F04GlobalField::getObject();
    //  if (myGlobalField!=NULL)  {delete myGlobalField;}

    return 0;
}


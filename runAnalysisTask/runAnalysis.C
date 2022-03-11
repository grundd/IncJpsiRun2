// runAnalysis.C 
// can be used to run the AnalysisTaskCentralJpsi_DG task
// in ~/alice/AliPhysics/PWGUD/UPC/

// options to set:
// iDataset:
// 0 => LHC18qr data
// 1 => kCohJpsiToEl
// 2 => kCohJpsiToMu
// 3 => kCohPsi2sToElPi
// 4 => kCohPsi2sToMuPi
// 5 => kIncohJpsiToEl
// 6 => kIncohJpsiToMu
// 7 => kIncohPsi2sToElPi
// 8 => kIncohPsi2sToMuPi
// 9 => kTwoGammaToElMedium
// 10 => kTwoGammaToMuMedium

TString strMC[10] = {"kCohJpsiToEl", "kCohJpsiToMu",
                     "kCohPsi2sToElPi", "kCohPsi2sToMuPi",
                     "kIncohJpsiToEl", "kIncohJpsiToMu",
                     "kIncohPsi2sToElPi", "kIncohPsi2sToMuPi",
                     "kTwoGammaToElMedium", "kTwoGammaToMuMedium"};

void runAnalysis(Int_t iDataset, Bool_t isNeutral = kFALSE, const char* suffix = "")
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t local = kTRUE;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t gridTest = kTRUE;

    // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif
     
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("myTask");
    // ESD input handler
    AliESDInputHandler *esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
    // MC handler (is iDataset > 0 => MC)
    AliMCEventHandler *MCHandler = NULL;
    if(iDataset > 0){
        MCHandler = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(MCHandler);
        if (!MCHandler) Printf("Could not retrieve MC event handler.");
    }

    // compile the class and load the add task macro
#if !defined (__CINT__) || defined (__CLING__)
    // ROOT6
    gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    char txt_cmd[120];
    sprintf(txt_cmd,"$ALICE_PHYSICS/PWGUD/UPC/AddTaskCentralJpsi_DG.C(%d,\"%s\")", isNeutral, suffix);
    AliAnalysisTaskCentralJpsi_DG *task = reinterpret_cast<AliAnalysisTaskCentralJpsi_DG*>(gInterpreter->ExecuteMacro(txt_cmd));
#else
    // ROOT5
    // ...
#endif

    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("esdTree");
        // add a few files to the chain (change this so that your local files are added)
        if(iDataset == 0) chain->Add("/home/david/alice/ALICEInputData_LocalTests/LHC18qr_pass1_000296623/101_AliESDs.root");
        if(iDataset > 0){
            TString path = "/home/david/alice/ALICEInputData_LocalTests/MC_LHC20g1_295585_001_" + strMC[iDataset-1] + "/AliESDs.root";
            chain->Add(path.Data());
        }
        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);

    // ================================================================
    // ======================== Running on GRID =======================
    // ================================================================

    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        /*
        alienHandler->SetAdditionalLibs("AliAnalysisTaskCentralJpsi_DG.cxx AliAnalysisTaskCentralJpsi_DG.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskCentralJpsi_DG.cxx");
        */
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20220301_ROOT6-1"); // Guillermo: new version might solve a problem with duplicity of some runs
        // Guillermo: vAN-20200221_ROOT6-1
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        
        // select the input files
        // if data
        if(iDataset == 0){
            alienHandler->SetGridDataDir("/alice/data/2018/LHC18q");   
            alienHandler->SetDataPattern("/*pass3/*AliESDs.root");
            alienHandler->SetRunPrefix("000");
            alienHandler->AddRunNumber(295937); 
            // working directory
            alienHandler->SetGridWorkingDir("try_data");
            alienHandler->SetExecutable("try_data.sh");
            alienHandler->SetJDLName("try_data.jdl");
        }
        // if MC
        if(iDataset > 0){
            alienHandler->SetGridDataDir(Form("/alice/sim/2020/LHC20g1/%s", strMC[iDataset-1].Data()));
            alienHandler->SetDataPattern("/*AliESDs.root");
            // no run prefix for MC!
            alienHandler->AddRunNumber(295585); 
            // working directory
            alienHandler->SetGridWorkingDir(Form("try_MC_%s", strMC[iDataset-1].Data()));
            alienHandler->SetExecutable(Form("try_MC_%s.sh", strMC[iDataset-1].Data()));
            alienHandler->SetJDLName(Form("try_MC_%s.jdl", strMC[iDataset-1].Data()));
        }

        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(40);
        // specify how many seconds your job may take
        alienHandler->SetTTL(10000);
        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate") 
        // (see below) mode, set SetMergeViaJDL(kFALSE) 
        // to collect final results
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(kTRUE); // set to kFALSE to download merged results
        // define the output folders
        alienHandler->SetGridOutputDir("myOutputDir");
        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);

        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test"); // "test" or "terminate"
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            alienHandler->SetRunMode("full"); // "full" or "terminate"
            mgr->StartAnalysis("grid");
        }
    }
    
    // if local or gridTest, create subfolder and store the results
    if(local == kTRUE || (local == kFALSE && gridTest == kTRUE)){
        TString str_folder = "";
        // if local or grid
        if(local) str_folder += "local_";
        else      str_folder += "grid_";
        // if LHC18qr data or MC
        if(iDataset == 0) str_folder += "data";
        else {
            str_folder += "MC_" + strMC[iDataset-1];
            // if feed-down datasets
            if(iDataset == 3 || iDataset == 4 || iDataset == 7 || iDataset == 8){
                if(isNeutral) str_folder += "_neutral";
                else          str_folder += "_charged";
            }
        }             
        // move the results to the new folder
        gSystem->Exec("rm -r " + str_folder);
        gSystem->Exec("mkdir -p " + str_folder);
        gSystem->Exec("mv AnalysisResults.root track_cuts.root " + str_folder);
    }

    return;
}

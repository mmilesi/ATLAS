// EDM include(s):
#include "xAODRootAccess/Init.h"

// SH include(s):
#include "SampleHandler/Sample.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "SampleHandler/DiskListLocal.h"
#include "SampleHandler/DiskListEOS.h"

// EL include(s):
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "EventLoop/LSFDriver.h"

// package include(s):
#include "HTopMultilepAnalysis/HTopMultilepAnalysis.h"
#include "HTopMultilepAnalysis/HTopMultilepEventSelector.h"
#include "HTopMultilepAnalysis/TruthMatchAlgo.h"
#include "HTopMultilepAnalysis/HTopMultilepTreeAlgo.h"
#include "xAODAnaHelpers/BasicEventSelection.h"
#include "xAODAnaHelpers/JetCalibrator.h"
#include "xAODAnaHelpers/JetSelector.h"
#include "xAODAnaHelpers/BJetEfficiencyCorrector.h"
#include "xAODAnaHelpers/JetHistsAlgo.h"
#include "xAODAnaHelpers/MuonCalibrator.h"
#include "xAODAnaHelpers/MuonEfficiencyCorrector.h"
#include "xAODAnaHelpers/MuonSelector.h"
#include "xAODAnaHelpers/ElectronCalibrator.h"
#include "xAODAnaHelpers/ElectronEfficiencyCorrector.h"
#include "xAODAnaHelpers/ElectronSelector.h"
#include "xAODAnaHelpers/OverlapRemover.h"

// c++ include(s):
#include <string>
#include <sstream>
#include <fstream>
#include <time.h>
#include <stdio.h>

// ROOT include(s):
#include "TSystem.h"
#include "TEnv.h"


using namespace std;

bool cmdline( int argc, char** argv, map<string,string>& opts );
void usage();
string currentDateTime(); 

int main ( int argc, char **argv ) {
    map<string, string> opts;
    if ( !cmdline(argc, argv, opts) ) { return 0; }

    bool grid(false), lxbatch(false), inclusive(false), smearing(false), weights(false), make_histos(false), debug(false); 
    int maxEvents(-1), skipTillEvent(-1);
    stringstream ss_grid;    ss_grid    << opts["grid"]          << " "; ss_grid    >> grid;
    stringstream ss_lxbatch; ss_lxbatch << opts["lxbatch"]       << " "; ss_lxbatch >> lxbatch;    
    stringstream ss_inc;     ss_inc     << opts["inclusive"]     << " "; ss_inc     >> inclusive;
    stringstream ss_sm;      ss_sm      << opts["smearing"]      << " "; ss_sm      >> smearing;
    stringstream ss_w;       ss_w       << opts["weights"]       << " "; ss_w       >> weights;
    stringstream ss_mh;      ss_mh      << opts["makeHistos"]    << " "; ss_mh      >> make_histos;    
    stringstream ss_d;       ss_d       << opts["debug"]         << " "; ss_d       >> debug;
    stringstream ss_maxevts; ss_maxevts << opts["maxEvents"]     << " "; ss_maxevts >> maxEvents;
    stringstream ss_skip;    ss_skip    << opts["skipTillEvent"] << " "; ss_skip    >> skipTillEvent;

    string inGridDSName = opts["inGridDSName"];
    string gridUser     = opts["gridUser"];    
    string inGridDSList = opts["inGridDSList"];
    string inDir        = opts["inDir"];
    string inDSName     = opts["inDSName"];
    string inDSType     = opts["inDSType"];    
    string inFileName   = opts["inFileName"];
    string outDir_in    = opts["outDir"];
    string outDSID      = opts["outDSID"];
    
    if ( outDir_in.empty() ) {
        cout << "Please enter output directory!" << endl;
        return 0;
    }
    string outDir = outDir_in;

    xAOD::Init().ignore();

    // instantiate SampleHandler object
    SH::SampleHandler sh;
   
    if ( grid ){
        if ( inGridDSName.empty() && inGridDSList.empty() ) {
            cout << "Name of dataset on grid required, or provide dataset list!" << endl;
            return 0;
        }
        if ( gridUser.empty() ) {
            cout << "Name of grid User required!" << endl;
            return 0;
        }	
        if ( outDSID.empty() ) {
            cout << "Must provide an ID string for the output dataset!" << endl;
            return 0;
        }
	if ( !inGridDSList.empty() ){
	  // run jobs on the datasets passed via input text file (regardless of "inGridDSName" being empty or not)
      	  char token[500];
      	  FILE *f = fopen(inGridDSList.c_str(),"r");
      	  if ( !f ) {
      	    cout << "ERROR! Cannot open " << inGridDSList.c_str() << " !" << endl;
      	    return 0;
      	  }
	  unsigned int iLine(0);
      	  while ( !feof(f) ) {
      	    if ( !fgets(token,1000,f) ) { break; }
      	    stringstream ss;
      	    string DSName;
      	    ss << token;  ss >> DSName;
	    
	    string::iterator beg_it = DSName.begin();
	    string::iterator end_it = DSName.end();
	    
	    if ( *beg_it == '#' ) { continue; }
	    
	    // NB: need to compare *iterator to a character, not a string literal. Character constants use single quotes, e.g., '/'.
	    if ( *(end_it-1) != '/' ) {
	      cout << "WARNING! Dataset " << DSName.c_str() << "\n should end with '/'. Will not submit it " << endl;
	      continue;
	    }
	    
      	    SH::scanDQ2 (sh, DSName);
      	    cout << "Resolving grid dataset: "<< DSName <<endl;
	    
	    ++iLine;
      	  }
	} else {
	  // run a job on grid on a single file passed by command line
          SH::scanDQ2(sh, inGridDSName);  
      	  cout << "Resolving grid dataset: "<< inGridDSName << endl;	  
	}
	
    } else {
        
	if ( inDir.empty() ) {
            cout << "Name of local input folder required!" << endl;
            return 0;
        }
	
	size_t is_eos = inDir.find("/eos/");
	if ( is_eos != string::npos ) {

	  string inDir_eos = "root://eosatlas/" + inDir;
	  SH::DiskListEOS list_eos(inDir, inDir_eos);
	  
	  if ( !inDSName.empty() ) {
	    cout << "Processing all files in dataset: " << inDSName << endl;
	    SH::scanDir (sh, list_eos, "*.root*", inDSName); // scans all files only in the dataset passed
	  } else if ( !inFileName.empty() ){
	    cout << "Processing only file: " << inFileName << endl;
            SH::scanDir(sh, list_eos, inFileName); // scans only the filename passed	  
	  } else {
	    cout << "Processing all files in all datasets" << endl;
	    SH::scanDir (sh, list_eos);  // scans all the files in all the dataset subfolders
	  }
	  
	} else {
	
	  SH::DiskListLocal list(inDir);
	  
	  if ( !inDSName.empty() ) {
	    cout << "Processing all files in dataset: " << inDSName << endl;
	    SH::scanDir (sh, list, "*.root*", inDSName); // scans all files only in the datset passed
	  } else if ( !inFileName.empty() ){
	    cout << "Processing only file: " << inFileName << endl;
            SH::scanDir(sh, list, inFileName); // scans only the filename passed	  
	  } else {
	    cout << "Processing all files in all datasets" << endl;
	    SH::scanDir (sh, inDir);  // scans all the files in all the dataset subfolders
	  }
	  
	}
    }

    sh.setMetaString("nc_tree", "CollectionTree");

    // set to process only root files in a dataset also containing log files
    sh.setMetaString ("nc_grid_filter", "*.root*");

    sh.print();

    // access a single sample
    // Sample *sample = sh.get ("mc14_13TeV. blahblahblah");
    // sample->setMetaString ("SimulationFlavour", "AFII");

    // ************************************************************
    // **** Instantiating algorithms   

    EL::Job job;
    job.sampleHandler(sh);

    std::string localDataDir = "$ROOTCOREBIN/data/HTopMultilepAnalysis/";

    //
    // 1. xAODAnaHelpers algorithms
    //    
    // basic event selection : this is mandatory!
    BasicEventSelection* baseEventSel             = new BasicEventSelection();
    baseEventSel->setName("baseEventSel")->setConfig(localDataDir+"Event/"+"baseEvent_HTopMultilep"+"_"+inDSType+".config");
    JetCalibrator* jetCalib                       = new JetCalibrator();
    jetCalib->setName("jetCalib_AntiKt4EMTopo")->setConfig(localDataDir+"Jets/"+"jetCalib_AntiKt4TopoEMCalib"+"_"+inDSType+".config");
    MuonCalibrator* muonCalib                     = new MuonCalibrator();
    muonCalib->setName("muonCalib")->setConfig(localDataDir+"Muons/"+"muonCalib"+"_"+inDSType+".config");
    ElectronCalibrator* electronCalib             = new ElectronCalibrator();
    electronCalib->setName("electronCalib")->setConfig(localDataDir+"Electrons/"+"electronCalib"+"_"+inDSType+".config");
    MuonEfficiencyCorrector*      muonEffCorr     = new MuonEfficiencyCorrector();
    muonEffCorr->setName("muonEfficiencyCorrector")->setConfig(localDataDir+"Muons/"+"muonEffCorr"+"_"+inDSType+".config");
    ElectronEfficiencyCorrector*  electronEffCorr = new ElectronEfficiencyCorrector();
    electronEffCorr->setName("electronEfficiencyCorrector")->setConfig(localDataDir+"Electrons/"+"electronEffCorr"+"_"+inDSType+".config");
    JetSelector* jetSelect_selection              = new JetSelector();
    jetSelect_selection->setName("jetSelect_selection")->setConfig(localDataDir+"Jets/"+"jetSelect_HTopMultilep"+"_"+inDSType+".config");
    MuonSelector* muonSelect_preselection         = new MuonSelector();
    muonSelect_preselection->setName("muonSelect_preselection")->setConfig(localDataDir+"Muons/"+"muonPreSelect_HTopMultilep"+"_"+inDSType+".config");
    ElectronSelector* electronSelect_preselection   = new ElectronSelector();
    electronSelect_preselection->setName("electronSelect_preselection")->setConfig(localDataDir+"Electrons/"+"electronPreSelect_HTopMultilep"+"_"+inDSType+".config");
    BJetEfficiencyCorrector* bjetEffCorr_btag     = new BJetEfficiencyCorrector();
    bjetEffCorr_btag->setName("bjetEffCor_btag")->setConfig(localDataDir+"Jets/"+"bjetEffCorr"+"_"+inDSType+".config");
    MuonSelector* muonSelect_selection            = new MuonSelector();
    muonSelect_selection->setName("muonSelect_selection")->setConfig(localDataDir+"Muons/"+"muonSelect_HTopMultilep"+"_"+inDSType+".config");
    ElectronSelector* electronSelect_selection   = new ElectronSelector();
    electronSelect_selection->setName("electronSelect_selection")->setConfig(localDataDir+"Electrons/"+"electronSelect_HTopMultilep"+"_"+inDSType+".config");
    OverlapRemover* overlapRemoval                = new OverlapRemover();
    overlapRemoval->setName("overlap_removal")->setConfig(localDataDir+"OverlapRemoval/"+"overlapRemoval_HTopMultilep"+"_"+inDSType+".config");

    //
    // 2. HTopMultilepAnalysis algorithms
    // 
    HTopMultilepEventSelector* eventSelect       = new HTopMultilepEventSelector();
    eventSelect->setName("eventSelect_skim")->setConfig(localDataDir+"Event/"+"eventSelect_HTopMultilep"+"_"+inDSType+".config");
    TruthMatchAlgo* truthMatching                = new TruthMatchAlgo();
    truthMatching->setName("truthMatching")->setConfig(localDataDir+"Analysis/"+"analysis_TruthMatchingLeptons"+"_"+inDSType+".config");
    HTopMultilepAnalysis* analysis 		 = new HTopMultilepAnalysis();
    analysis->setName("multilep_analysis")->setConfig(localDataDir+"Analysis/"+"analysis_HTopMultilep"+"_"+inDSType+".config");
    HTopMultilepTreeAlgo* out_tree 		 = new  HTopMultilepTreeAlgo();
    out_tree->setName("physics")->setConfig(localDataDir+"Tree/"+"tree_HTopMultilep"+"_"+inDSType+".config");

    // dump plots - jets
    JetHistsAlgo* jetHistsAlgo_signal(nullptr);
    JetHistsAlgo* jetHistsAlgo_all(nullptr); 
    if ( make_histos ) {
      JetHistsAlgo* jk_AntiKt4EM                   = new JetHistsAlgo();
      jk_AntiKt4EM->setName("AntiKt4EM/")->setConfig(localDataDir+"Jets/"+"jetHistsAlgo_signal.config");
    }  
          
    size_t isDxAOD = inDSType.find("DxAOD");

    // Add all the algorithms to the EL::job - Here order matters!
    job.algsAdd( baseEventSel );
    job.algsAdd( jetCalib );
    job.algsAdd( muonCalib ); 
    job.algsAdd( muonSelect_preselection );
    job.algsAdd( electronCalib );
    job.algsAdd( electronSelect_preselection );
    job.algsAdd( jetSelect_selection );
    job.algsAdd( bjetEffCorr_btag );
    job.algsAdd( overlapRemoval );
    job.algsAdd( muonSelect_selection );
    job.algsAdd( muonEffCorr );
    job.algsAdd( electronSelect_selection );
    job.algsAdd( electronEffCorr );
    job.algsAdd( eventSelect );
    job.algsAdd( truthMatching );
    job.algsAdd( analysis );
    if ( make_histos ) {
      job.algsAdd( jetHistsAlgo_all );
      job.algsAdd( jetHistsAlgo_signal );
    }
    job.algsAdd( out_tree );
    
    // ************************************************************
    // **** submit the job to the driver
    
    outDir = outDir + "_" + currentDateTime();
    
    // set maximum number of events to run on (not available for grid driver)
    if ( !grid && maxEvents > 0 ) { 
      job.options()->setDouble (EL::Job::optMaxEvents, maxEvents);
      cout << "Running only on " <<  maxEvents << " events" << endl;
    }
    // skip events till...
    if ( !grid && skipTillEvent > 0 ) {
      job.options()->setDouble (EL::Job::optSkipEvents, skipTillEvent);
      cout << "Skipping events up to " << skipTillEvent << endl;      
    }

    // TTreeCache - specify the size of the cache you would like for your job before submitting it (in this case 10MB)
    job.options()->setDouble (EL::Job::optCacheSize, 10*1024*1024);
    // assume data pattern will not change after n files processed 
    job.options()->setDouble (EL::Job::optCacheLearnEntries, 20);
 
    if ( lxbatch ) {
      
      EL::LSFDriver driver;
      driver.options()->setString (EL::Job::optSubmitFlags, "-L /bin/bash"); // or whatever shell you are using
      driver.shellInit = "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase && source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh && eval rcSetup Base,2.1.30"; // or whatever version of the Analysis Release you are using    
      driver.submitOnly(job, outDir);
    }
    else if ( grid ) {
        
	EL::PrunDriver driver;

	// force user to check output dataset ID string - recommended is 'HTopMultilep.some_number'
	cout << "Are you sure to apply the following ID tag to the grid output dataset name?\t" << outDSID << "\n" 
	     << "recommended name is HTopMultilep.test#SOME_NUMBER# "          << "\n"
	     << "(please answer (yes/y/Y) to confirm, any key to reject" << "\n"
	     << endl;
	bool check = false;
	while ( !check ) {
	  string decision; cin >> decision;
	  if ( decision == "yes" || decision == "y" || decision ==  "Y" )  { check = true; }
	  else                                                             { cout << "Please choose the ID tag again" << endl; continue; }
	}
	
	// *************************** //
	// Grid job submission options //
	// *************************** //	
	
	driver.options()->setString("nc_outputSampleName", "user."+gridUser+"."+outDSID+".%in:name[2]%"+".%in:name[3]%"); // numbers in "[]" refer to the part of the string separated by "." in the inDS name
        
	// this does not work!
	// driver.options()->setString(EL::Job::optGridNFilesPerJob, "MAX"); //By default, split in as few jobs as possible 	
	
	// other useful EL grid driver options
        driver.options()->setDouble(EL::Job::optGridMergeOutput, 1.0);  //run merging jobs for all samples before downloading (recommended) 
	driver.options()->setDouble(EL::Job::optGridNFilesPerJob, 1.0); // default: 50
	//driver.options()->setDouble(EL::Job::optGridMemory, 10240); // send jobs on nodes w/ high RAM (10 GB) - useful if you have memory leaks :)
 	// set 1 file per job only for data (speeds up a lot!)
        for ( SH::SampleHandler::iterator sample_it = sh.begin(); sample_it != sh.end(); ++sample_it ) {
          SH::Sample *sample = *sample_it;
          std::string  sample_name = sample->name();
          std::size_t found_data_in_name = sample_name.find("data");
          if ( found_data_in_name == std::string::npos)  { continue; }
          sh.get(sample_name)->setMetaDouble(EL::Job::optGridNFilesPerJob, 1.0);
        }
	// driver.options()->setString(EL::Job::optGridExcludedSite, "ANALY_IN2P3-CC-T2,ANALY_IN2P3-CC"); // can list them just by separating grid site names w/ comma 
	driver.options()->setString(EL::Job::optGridSite, "UKI-NORTHGRID-MAN-HEP_LOCALGROUPDISK"); // run on this specific grid site 
	
	cout << "outDS: " << driver.options()->castString("nc_outputSampleName") << endl;
	
	if ( !debug ) {
	  cout << "Submitting job to grid driver! " << endl;
	  // submit the job, returning control immediately if the driver supports it.
	  driver.submitOnly(job, outDir);
	}
	
    } else {
        EL::DirectDriver driver;
        driver.submit(job, outDir);
    }

    return 0;
}

bool cmdline( int argc, char** argv, map<string,string>& opts ) {
    opts.clear();

    // defaults
    opts["outDir"] = "";
    opts["grid"] = "0";
    opts["lxbatch"] = "0";    
    opts["gridUser"] = "";    
    opts["inDir"] = "";
    opts["inDSName"] = "";   
    opts["inDSType"] = "xAOD";      
    opts["inGridDSList"] = ""; 
    opts["outDSID"] = "";
    opts["inFileName"] = "";
    opts["inGridDSName"] = "";
    opts["inclusive"] = "0";
    opts["smearing"] = "0";
    opts["weights"] = "0";
    opts["makeHistos"] = "0";    
    opts["debug"] = "0";
    opts["maxEvents"] = "-1";
    opts["skipTillEvent"] = "-1";
    
    for ( int i=1 ;i<argc ; i++) {

        string opt=argv[i];

        if ( opt=="--help" ) { usage(); return false; }

        if ( 0 != opt.find("--") ) {
            cout<<"ERROR: options must start with '--'!"<<endl;
            return false;
        }
        opt.erase(0,2);
        if ( opts.find(opt) == opts.end() ) {
            cout<<"ERROR: invalid option '"<<opt<<"'!"<<endl;
            return false;
        }
        string nxtopt=argv[i+1];
        if ( 0 == nxtopt.find("--") || i+1 >= argc ) {
            cout<<"ERROR: option '"<<opt<<"' requires value!"<<endl;
            return false;
        }

        opts[opt] = nxtopt;
        i++;
    }

    return true;
}

void usage()
{
    cout<<"USAGE: run [-option value]\n\n"
        <<"options \t [default]:\n\n"
        <<"-outDir \t (required!)\n"
        <<"-grid \t [0] (set it to 1 to run on grid)\n"
        <<"-gridUser \t (required to run on grid)\n"	
        <<"-lxbatch \t [0] (set it to 1 to run on lxplus batch system)\n"	
        <<"-inDir \t (required for local run)\n"
        <<"-inDSName \t (only for local run)\n"
        <<"-inDSType \t (can be 'xAOD', 'DxAOD-19.1.4.7' etc.)\n"	
        <<"-inGridDSList \t (.txt file conatining list of datasets)\n"	
        <<"-outDSID \t (only for grid run)\n"
        <<"-inFileName \t (only for local run)\n"
        <<"-inGridDSName \t (only for grid run - wildcards accepted!)\n"
        <<"-inclusive \t [0]\n"
        <<"-smearing \t [0]\n"
        <<"-weight \t [0]\n"
        <<"-makeHistos \t [0] (set it to 1 if you want to dump histograms)\n"	
        <<"-debug \t [0]\n"
	<<"-maxEvents \t [ALL] \t (set number of events to run on - available only for local run)\n"
	<<"-skipTillEvent \t [NONE] \t (available only for local run)\n"	
        <<endl;

    return;
}


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

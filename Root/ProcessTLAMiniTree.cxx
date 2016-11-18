#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include "EventLoop/OutputStream.h"

#include <TLAAlgos/ProcessTLAMiniTree.h>
#include <xAODAnaHelpers/HelperFunctions.h>
#include <xAODAnaHelpers/tools/ReturnCheck.h>

#include "TFile.h"
#include "TKey.h"
#include "TLorentzVector.h"
#include "TEnv.h"
#include "TSystem.h"

#include <utility>      
#include <iostream>
#include <fstream>

using namespace std;

// this is needed to distribute the algorithm to the workers
ClassImp(ProcessTLAMiniTree)

ProcessTLAMiniTree :: ProcessTLAMiniTree () :
  m_applyGRL(false),
  // m_GRLxml("$ROOTCOREBIN/data/xAODAnaHelpers/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml"),
  m_GRLxml("$ROOTCOREBIN/data/TLAAlgos/data16_13TeV.periodAllYear_DetStatus-v82-pro20-13_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml"),
  m_applyTLALArEventVetoData(false),
  m_TLALArEventVetoFiles("$ROOTCOREBIN/data/TLAEventCleaning/event-veto-data/"),
  m_debug(false),
  m_doTrigger(false),
  m_isDijetNtupleTruth(false),
  m_isTLANtupleTruth(false),
  m_isTLANtupleTrig(false),
  m_isTLANtupleOffline(false),
  m_doTruthOnly(false),

  m_doSecondaryJets(false),
  m_secondaryJetName(""),
  m_doPileupFromMap(false),
  m_pileupMap(""),
  m_avgIntPerX(-1),

  m_YStarCut(99),
  m_YBoostCut(99),
  m_doBlind(false),
  m_doData(false),
  m_useWeighted(false),
  m_doCleaning(false),
  m_etaCut(1.8),
  m_leadJetPtCut(200.),
  m_subleadJetPtCut(200.),
  m_lumi(0),
  m_sampleEvents(0),
  m_grl(nullptr),
  m_jet_pt(0),
  m_jet_eta(0),
  m_jet_phi(0),
  m_jet_E(0),
  m_jet_muonSegments(0),
  m_jet_EMFrac(0),
  m_jet_HECFrac(0),
  m_LArError(false),
  m_jet_timing(0),
  m_jet_negativeE(0),
  m_MHT(0),
  m_mjj(0),
  m_jet_clean_passLooseBad(0),
  m_passedTriggers(nullptr),
  m_triggerPrescales(nullptr),
  /*m_applySF(false),
  m_hcalibration(nullptr),
  m_scaleFactorLocation(""),
  m_scaleFactorHistoName(""),
  m_pt_freeze(2500),
  m_eta_freeze(4.5),*/
  hIncl(nullptr),
  hCentral(nullptr),
  hCrack(nullptr),
  hEndcap(nullptr),
  hIncl_mjjWindow(nullptr),
  hCentral_mjjWindow(nullptr),
  hCrack_mjjWindow(nullptr),
  hEndcap_mjjWindow(nullptr),

  hSecIncl(nullptr),
  hSecCentral(nullptr),
  hSecCrack(nullptr),
  hSecEndcap(nullptr),
  hSecIncl_mjjWindow(nullptr),
  hSecCentral_mjjWindow(nullptr),
  hSecCrack_mjjWindow(nullptr),
  hSecEndcap_mjjWindow(nullptr),


  m_secJet_pt(0),
  m_secJet_eta(0),
  m_secJet_phi(0),
  m_secJet_E(0),
  m_secJet_muonSegments(0),
  m_secJet_EMFrac(0),
  m_secJet_HECFrac(0),
  m_secJet_timing(0),
  m_secJet_negativeE(0)

//  hOffline(nullptr),
//  hTrigger(nullptr),
//  hTruth(nullptr),
{
  Info("ProcessTLAMiniTree()", "Calling constructor");
    
  /*m_applyGRL      = true;
  m_GRLxml   = "$ROOTCOREBIN/data/xAODAnaHelpers/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml";  //https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/GoodRunListsForAnalysis
  m_GRLxml   = "$ROOTCOREBIN/data/xAODAnaHelpers/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml";  //https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/GoodRunListsForAnalysis*/
}


EL::StatusCode  ProcessTLAMiniTree :: configure ()
{

//  //
//  // Read Input from .config file
//  //
//  m_debug                    = config->GetValue("Debug" ,          m_debug );
//  m_doTrigger                = config->GetValue("DoTrigger" ,      m_doTrigger );
//  m_YStarCut                 = config->GetValue("YStarCut" ,       m_YStarCut);
//  m_YBoostCut                = config->GetValue("YBoostCut" ,      m_YBoostCut);
//  m_doBlind                  = config->GetValue("DoBlind" ,        m_doBlind );
//  m_doData                   = config->GetValue("DoData" ,         m_doData );
//  m_useWeighted              = config->GetValue("UseWeighted" ,    m_useWeighted );
//  m_doCleaning               = config->GetValue("DoCleaning" ,     m_doCleaning );
//  m_etaCut                   = config->GetValue("AbsEtaCut" ,      m_etaCut );
//  m_leadJetPtCut             = config->GetValue("LeadJetPtCut" ,   m_leadJetPtCut );
//  m_subleadJetPtCut          = config->GetValue("subLeadJetPtCut" ,   m_leadJetPtCut );
//  m_lumi                     = config->GetValue("Lumi",            m_lumi );
//  m_applyGRL                 = config->GetValue("ApplyGRL",        m_applyGRL);
//  m_GRLxml                   = config->GetValue("GRL",             m_GRLxml.c_str());

  //histograms that are always there, defaulting to what is in the 'jet' branch
  //this is the distribution we cut on
    
  if (m_isTLANtupleTrig) {
    hIncl = new eventHists("TriggerJets"  ,       wk());
    hCentral = new eventHists("TriggerJets_0-12"  ,       wk());
    hCrack = new eventHists("TriggerJets_12-16"  ,       wk());
    hEndcap = new eventHists("TriggerJets_16-28"  ,       wk());
      
    hIncl_mjjWindow = new eventHists("TriggerJets_mjjWindow"  ,       wk());
    hCentral_mjjWindow = new eventHists("TriggerJets_0-12_mjjWindow"  ,       wk());
    hCrack_mjjWindow = new eventHists("TriggerJets_12-16_mjjWindow"  ,       wk());
    hEndcap_mjjWindow = new eventHists("TriggerJets_16-28_mjjWindow"  ,       wk());
      
  }
  else if (m_isTLANtupleOffline) {
   hIncl = new eventHists("OfflineJets"  ,       wk());
   hCentral = new eventHists("OfflineJets_0-12"  ,       wk());
   hCrack = new eventHists("OfflineJets_12-16"  ,       wk());
   hEndcap = new eventHists("OfflineJets_16-28"  ,       wk());
      
   hIncl_mjjWindow = new eventHists("OfflineJets_mjjWindow"  ,       wk());
   hCentral_mjjWindow = new eventHists("OfflineJets_0-12_mjjWindow"  ,       wk());
   hCrack_mjjWindow = new eventHists("OfflineJets_12-16_mjjWindow"  ,       wk());
   hEndcap_mjjWindow = new eventHists("OfflineJets_16-28_mjjWindow"  ,       wk());

  }
  else if (m_isDijetNtuple) {
   hIncl    = new eventHists("TriggerJets"        ,       wk());
   hCentral = new eventHists("TriggerJets_0-12"   ,       wk());
   hCrack   = new eventHists("TriggerJets_12-16"  ,       wk());
   hEndcap  = new eventHists("TriggerJets_16-28"  ,       wk());
      
   hIncl_mjjWindow    = new eventHists("TriggerJets_mjjWindow"        ,       wk());
   hCentral_mjjWindow = new eventHists("TriggerJets_0-12_mjjWindow"   ,       wk());
   hCrack_mjjWindow   = new eventHists("TriggerJets_12-16_mjjWindow"  ,       wk());
   hEndcap_mjjWindow  = new eventHists("TriggerJets_16-28_mjjWindow"  ,       wk());

   if(m_doSecondaryJets) {
     hSecIncl    = new eventHists("TriggerJetsSec"        ,       wk());
     hSecCentral = new eventHists("TriggerJetsSec_0-12"   ,       wk());
     hSecCrack   = new eventHists("TriggerJetsSec_12-16"  ,       wk());
     hSecEndcap  = new eventHists("TriggerJetsSec_16-28"  ,       wk());
     
     hSecIncl_mjjWindow    = new eventHists("TriggerJetsSec_mjjWindow"        ,       wk());
     hSecCentral_mjjWindow = new eventHists("TriggerJetsSec_0-12_mjjWindow"   ,       wk());
     hSecCrack_mjjWindow   = new eventHists("TriggerJetsSec_12-16_mjjWindow"  ,       wk());
     hSecEndcap_mjjWindow  = new eventHists("TriggerJetsSec_16-28_mjjWindow"  ,       wk());
   }
  }

  else hIncl = new eventHists("Incl"  ,       wk());
  //histograms that are there for comparisons...for later
  /*if (m_isTLANtupleOffline) hOffline = new eventHists("Offline"  ,       wk());
  if (m_isTLANtupleTrig) hTrigger = new eventHists("Trigger"  ,       wk());
  if (m_isTLANtupleTruth || m_isDijetNtupleTruth) hTruth = new eventHists("Truth"  ,       wk());*/
  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ProcessTLAMiniTree :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessTLAMiniTree :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  Info("histInitialize()", "Calling histInitialize \n");

    //make the test histogram
    
  m_h2_LArError = new TH2D("h2_LArError", "h2_LArError", 2, 0, 2, 2, 0, 2);
  wk()->addOutput(m_h2_LArError);
    
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessTLAMiniTree :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  Info("fileExecute()", "Calling fileExecute"); // CWK debugging  

  TFile* inputFile = wk()->inputFile();
  cout << inputFile << ", " << inputFile->GetName() << endl;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessTLAMiniTree :: changeInput (bool firstFile)
{

  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  Info("changeInput()", "calling changeInput() \n"); // CWK debugging  

  if(firstFile){
    if ( this->configure() == EL::StatusCode::FAILURE ) {
      Error("initialize()", "Failed to properly configure. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    TFile* inputFile = wk()->inputFile();

    //do not always use cutflow to normalize, in case of power-law
    if (m_useCutflow) {
        
      TIter next(inputFile->GetListOfKeys());
      TKey *key;
      while ((key = (TKey*)next())) {
          
        std::string keyName = key->GetName();
          
        std::size_t found = keyName.find("cutflow");
        bool foundCutFlow = (found!=std::string::npos);

        found = keyName.find("weighted");
        bool foundWeighted = (found!=std::string::npos);

        if(foundCutFlow){
	
	      if(m_useWeighted && foundWeighted){
	        cout << "Getting NSample events from " << keyName << endl;
	        m_sampleEvents = ((TH1F*)key->ReadObj())->GetBinContent(1);
	        cout << "Setting Sample events to:" << m_sampleEvents << endl;
	        cout << "Setting Lumi to: " << m_lumi << endl;
	  
	        std::string inFileName = inputFile->GetName();
	        //std::size_t foundJZ5Pos = inFileName.find("JZ5");
            //bool foundJZ5 = (foundJZ5Pos!=std::string::npos);
	        //if(foundJZ5){
	        //  m_sampleEvents = (m_sampleEvents - 352219559.60316455) ;
	        //  cout << "Correcting Sample events to: " << m_sampleEvents << endl;
	        //}
	  
	      }//if using weighted cutflow
        }//Found CF
      }//over Keys
    }//if using cutflow
      
    else {
        
      cout << "Not using cutflow, assuming normalization will happen later" << m_sampleEvents << endl;
      cout << "Setting Lumi to: " << m_lumi << endl;
      cout << "Setting Sample events to 1.0" << endl;
      m_sampleEvents = 1.0;

    }
      
  }// do only on first file

  TTree *tree = wk()->tree();
  tree->SetBranchStatus ("*", 0);

  //event-level variables
    
  tree->SetBranchStatus  ("runNumber",    1);
  tree->SetBranchAddress ("runNumber",    &m_runNumber);

  tree->SetBranchStatus  ("eventNumber",    1);
  tree->SetBranchAddress ("eventNumber",    &m_eventNumber);


  if(!m_doTruthOnly){
    tree->SetBranchStatus  ("lumiBlock",    1);
    tree->SetBranchAddress ("lumiBlock",    &m_lumiBlock);

    tree->SetBranchStatus  ("LArError",    1);
    tree->SetBranchAddress ("LArError",    &m_LArError);
      
    tree->SetBranchStatus  ("timeStampNSOffset",    1);
    tree->SetBranchAddress ("timeStampNSOffset",    &m_timeStampNSOffset);
      
    tree->SetBranchStatus  ("timeStamp",    1);
    tree->SetBranchAddress ("timeStamp",    &m_timeStamp);
      
    tree->SetBranchStatus  ("NPV",    1);
    tree->SetBranchAddress ("NPV",    &m_NPV);

    if(m_doTrigger){
      tree->SetBranchStatus  ("passedTriggers", 1);
      tree->SetBranchAddress ("passedTriggers", &m_passedTriggers);
      
      tree->SetBranchStatus  ("triggerPrescales", 1);
      tree->SetBranchAddress ("triggerPrescales", &m_triggerPrescales);
    }
  }

  tree->SetBranchStatus  ("weight", 1);
  tree->SetBranchAddress ("weight", &m_weight);

  tree->SetBranchStatus  ("weight_xs", 1);
  tree->SetBranchAddress ("weight_xs", &m_weight_xs);

  //truth-level variables
  if (m_isTLANtupleTruth) {
    tree->SetBranchStatus  ("truthJet_pt", 1);
    tree->SetBranchAddress ("truthJet_pt", &m_jet_pt);

    tree->SetBranchStatus  ("truthJet_eta", 1);
    tree->SetBranchAddress ("truthJet_eta", &m_jet_eta);

    tree->SetBranchStatus  ("truthJet_phi", 1);
    tree->SetBranchAddress ("truthJet_phi", &m_jet_phi);

    tree->SetBranchStatus  ("truthJet_E", 1);
    tree->SetBranchAddress ("truthJet_E", &m_jet_E);
  }

  else if (m_isDijetNtupleTruth) {
    tree->SetBranchStatus  ("jet_pt", 1);
    tree->SetBranchAddress ("jet_pt", &m_jet_pt);

    tree->SetBranchStatus  ("jet_eta", 1);
    tree->SetBranchAddress ("jet_eta", &m_jet_eta);

    tree->SetBranchStatus  ("jet_phi", 1);
    tree->SetBranchAddress ("jet_phi", &m_jet_phi);

    tree->SetBranchStatus  ("jet_E", 1);
    tree->SetBranchAddress ("jet_E", &m_jet_E);
  }

  //trigger-level variables
  else if (m_isTLANtupleTrig) {
      tree->SetBranchStatus  ("trigJet_pt", 1);
      tree->SetBranchAddress ("trigJet_pt", &m_jet_pt);
      
      tree->SetBranchStatus  ("trigJet_eta", 1);
      tree->SetBranchAddress ("trigJet_eta", &m_jet_eta);
      
      tree->SetBranchStatus  ("trigJet_phi", 1);
      tree->SetBranchAddress ("trigJet_phi", &m_jet_phi);
      
      tree->SetBranchStatus  ("trigJet_E", 1);
      tree->SetBranchAddress ("trigJet_E", &m_jet_E);

      tree->SetBranchStatus  ("trigJet_GhostMuonSegmentCount", 1);
      tree->SetBranchAddress ("trigJet_GhostMuonSegmentCount", &m_jet_muonSegments);

      tree->SetBranchStatus  ("trigJet_EMFrac", 1);
      tree->SetBranchAddress ("trigJet_EMFrac", &m_jet_EMFrac);
      
      tree->SetBranchStatus  ("trigJet_HECFrac", 1);
      tree->SetBranchAddress ("trigJet_HECFrac", &m_jet_HECFrac);

      tree->SetBranchStatus  ("trigJet_MHT", 1);
      tree->SetBranchAddress ("trigJet_MHT", &m_MHT);
      
      tree->SetBranchStatus  ("trigJet_mjj", 1);
      tree->SetBranchAddress ("trigJet_mjj", &m_mjj);

      if(!m_doTruthOnly){
          tree->SetBranchStatus  ("trigJet_clean_passLooseBad", 1);
          tree->SetBranchAddress ("trigJet_clean_passLooseBad", &m_jet_clean_passLooseBad);
      }

  
  }
    
  //reco-level variables
  else if (m_isTLANtupleOffline) {
      
      tree->SetBranchStatus  ("jet_pt", 1);
      tree->SetBranchAddress ("jet_pt", &m_jet_pt);
      
      tree->SetBranchStatus  ("jet_eta", 1);
      tree->SetBranchAddress ("jet_eta", &m_jet_eta);
      
      tree->SetBranchStatus  ("jet_phi", 1);
      tree->SetBranchAddress ("jet_phi", &m_jet_phi);
      
      tree->SetBranchStatus  ("jet_E", 1);
      tree->SetBranchAddress ("jet_E", &m_jet_E);
      
      tree->SetBranchStatus  ("jet_GhostMuonSegmentCount", 1);
      tree->SetBranchAddress ("jet_GhostMuonSegmentCount", &m_jet_muonSegments);
      
      tree->SetBranchStatus  ("jet_EMFrac", 1);
      tree->SetBranchAddress ("jet_EMFrac", &m_jet_EMFrac);
      
      tree->SetBranchStatus  ("jet_HECFrac", 1);
      tree->SetBranchAddress ("jet_HECFrac", &m_jet_HECFrac);
      
      tree->SetBranchStatus  ("jet_MHT", 1);
      tree->SetBranchAddress ("jet_MHT", &m_MHT);

      tree->SetBranchStatus  ("jet_mjj", 1);
      tree->SetBranchAddress ("jet_mjj", &m_mjj);
      
      if(!m_doTruthOnly){
          tree->SetBranchStatus  ("jet_clean_passLooseBad", 1);
          tree->SetBranchAddress ("jet_clean_passLooseBad", &m_jet_clean_passLooseBad);
      }

  }

  //truth-level variables
  else if (m_isDijetNtuple) {
    
      tree->SetBranchStatus  ("jet_pt", 1);
      tree->SetBranchAddress ("jet_pt", &m_jet_pt);
      
      tree->SetBranchStatus  ("jet_eta", 1);
      tree->SetBranchAddress ("jet_eta", &m_jet_eta);
      
      tree->SetBranchStatus  ("jet_phi", 1);
      tree->SetBranchAddress ("jet_phi", &m_jet_phi);
      
      tree->SetBranchStatus  ("jet_E", 1);
      tree->SetBranchAddress ("jet_E", &m_jet_E);

      if(!m_doTruthOnly) {
	// old NTUPs- prior to Oct 2016
	// tree->SetBranchStatus  ("jet__EMFrac", 1);
	// tree->SetBranchAddress ("jet__EMFrac", &m_jet_EMFrac);
	
	// tree->SetBranchStatus  ("jet__HECFrac", 1);
	// tree->SetBranchAddress ("jet__HECFrac", &m_jet_HECFrac);
	
	// tree->SetBranchStatus  ("jet__GhostMuonSegmentCount", 1);
	// tree->SetBranchAddress ("jet__GhostMuonSegmentCount", &m_jet_muonSegments);

	// new NTUPs - from 15.10.16
	tree->SetBranchStatus  ("jet_EMFrac", 1);
	tree->SetBranchAddress ("jet_EMFrac", &m_jet_EMFrac);
	
	tree->SetBranchStatus  ("jet_HECFrac", 1);
	tree->SetBranchAddress ("jet_HECFrac", &m_jet_HECFrac);
	
	tree->SetBranchStatus  ("jet_GhostMuonSegmentCount", 1);
	tree->SetBranchAddress ("jet_GhostMuonSegmentCount", &m_jet_muonSegments);

	tree->SetBranchStatus  ("jet_clean_passLooseBad", 1);
	tree->SetBranchAddress ("jet_clean_passLooseBad", &m_jet_clean_passLooseBad);	

	tree->SetBranchStatus  ("jet_Timing", 1);
	tree->SetBranchAddress ("jet_Timing", &m_jet_timing);

	tree->SetBranchStatus  ("jet_NegativeE", 1);
	tree->SetBranchAddress ("jet_NegativeE", &m_jet_negativeE);

      }

      if(m_doSecondaryJets) {
	tree->SetBranchStatus  ((m_secondaryJetName+"_pt").c_str(), 1);
	tree->SetBranchAddress ((m_secondaryJetName+"_pt").c_str(), &m_secJet_pt);
	
	tree->SetBranchStatus  ((m_secondaryJetName+"_eta").c_str(), 1);
	tree->SetBranchAddress ((m_secondaryJetName+"_eta").c_str(), &m_secJet_eta);
	
	tree->SetBranchStatus  ((m_secondaryJetName+"_phi").c_str(), 1);
	tree->SetBranchAddress ((m_secondaryJetName+"_phi").c_str(), &m_secJet_phi);
	
	tree->SetBranchStatus  ((m_secondaryJetName+"_E").c_str(), 1);
	tree->SetBranchAddress ((m_secondaryJetName+"_E").c_str(), &m_secJet_E);
	
	if(!m_doTruthOnly) {
	  // new NTUPs - from 15.10.16
	  tree->SetBranchStatus  ((m_secondaryJetName+"_EMFrac").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetName+"_EMFrac").c_str(), &m_secJet_EMFrac);
	  
	  tree->SetBranchStatus  ((m_secondaryJetName+"_HECFrac").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetName+"_HECFrac").c_str(), &m_secJet_HECFrac);
	  
	  tree->SetBranchStatus  ((m_secondaryJetName+"_GhostMuonSegmentCount").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetName+"_GhostMuonSegmentCount").c_str(), &m_secJet_muonSegments);
	  // cleaning done from primary
	  tree->SetBranchStatus  ((m_secondaryJetName+"_Timing").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetName+"_Timing").c_str(), &m_secJet_timing);
	  
	  tree->SetBranchStatus  ((m_secondaryJetName+"_NegativeE").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetName+"_NegativeE").c_str(), &m_secJet_negativeE);

	}
      }

  }



  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessTLAMiniTree :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
  m_eventCounter = -1;
  Info("initialize()", "Calling initialize \n"); // CWK debugging  
  if(m_applyGRL){
    std::cout<<"I am about to initialise with m_GRLxml = "<<m_GRLxml<<std::endl; // CWK debugging
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    std::vector<std::string> vecStringGRL;
    m_GRLxml = gSystem->ExpandPathName( m_GRLxml.c_str() );
    vecStringGRL.push_back(m_GRLxml);
    RETURN_CHECK("BasicEventSelection::initialize()", m_grl->setProperty( "GoodRunsListVec", vecStringGRL), "");
    RETURN_CHECK("BasicEventSelection::initialize()", m_grl->setProperty("PassThrough", false), "");
    RETURN_CHECK("BasicEventSelection::initialize()", m_grl->initialize(), "");
  }
  
  if(m_applyTLALArEventVetoData){

    m_dataForLArEventVeto = new TLALArEventVetoData();
    std::string LArEventVetoExpandedPath = gSystem->ExpandPathName( m_TLALArEventVetoFiles.c_str() );
    // replace this path with the path to your event veto data directory
    m_dataForLArEventVeto->loadFromDirectory(LArEventVetoExpandedPath);
  }
    

  if(m_doPileupFromMap) {
    m_pileupMap = gSystem->ExpandPathName( m_pileupMap.c_str() ); // in terms of $ROOTCOREBIN
    std::cout << "I am getting the pileup map from " << m_pileupMap << std::endl;
    TFile* pileupFile = new TFile(m_pileupMap.c_str());
    m_h2_pileupMap = (TH2F*)pileupFile->Get("h_runNumber_LB_pileup");
    m_h2_pileupMap->SetDirectory(0);
    pileupFile->Close();
  }

  Info("initialize()", "Succesfully initialized! \n");


  /*
  if (m_applySF) {
      
    if (m_scaleFactorLocation.EqualTo("") || m_scaleFactorHistoName.EqualTo("")) {
        cout<< "No scale factor found, will not apply" <<endl;
        m_applySF = false;
    }
    else {
  	    TFile* caliFile = TFile::Open(m_scaleFactorLocation);
	    m_hcalibration = (TH2D*) caliFile->Get(m_scaleFactorHistoName);
	    m_hcalibration->SetDirectory(0);
	    caliFile->Close();
	    delete caliFile;

	    m_pt_freeze = m_hcalibration->GetXaxis()->GetXmax();
	    m_eta_freeze = m_hcalibration->GetYaxis()->GetXmax();
	    if(m_debug) cout<< "Calibration edge in pT:" << m_hcalibration->GetXaxis()->GetXmax()<<endl;
	    if(m_debug) cout<< "Calibration edge in eta:"<< m_hcalibration->GetYaxis()->GetXmax()<<endl;
    }
  }
  */

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ProcessTLAMiniTree :: execute ()
{
  // Here you do everything that needs to be done on every single
  // event, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  //if(m_debug) Info("execute()", "Processing Event");
  ++m_eventCounter;

  wk()->tree()->GetEntry (wk()->treeEntry());
  unsigned njets       = m_jet_pt->size();
  unsigned nsecJets    = 0;
  if(m_doSecondaryJets) nsecJets = m_secJet_pt->size();
  float prescaleWeight = 1.0;
    
  if(m_eventCounter % 10000 == 0) cout << "executing event " << m_eventCounter << endl;
  else if(m_debug) cout << "executing event " << m_eventCounter << endl;

  if(m_doData){
   
    if ( m_applyGRL ) {
      if ( !m_grl->passRunLB( m_runNumber, m_lumiBlock ) ) {
	if(m_debug) cout << "GRL:: Fail Event " << endl;
        return EL::StatusCode::SUCCESS; // go to next event
      }else{
	//if(m_debug) cout << "GRL:: Pass Event " << endl;
      }
    }

    if(m_doPileupFromMap) {
      // std::cout << "RunNum: " << m_runNumber << ", LB: " << m_lumiBlock << std::endl;
      int runNumberBin = m_h2_pileupMap->GetYaxis()->FindBin(to_string(m_runNumber).c_str());
      // std::cout << "RN bin: " << runNumberBin << std::endl;
      m_avgIntPerX = m_h2_pileupMap->GetBinContent(m_lumiBlock, runNumberBin);
      // std::cout << "m_avgIntPerX: " << m_avgIntPerX << std::endl;
    }

  }//end of do data

  //
  // Minimum selection is a dijet
  //
    
  if(njets < 2 && nsecJets < 2){
    if(m_debug) cout << " Fail NJets " << endl;
    return EL::StatusCode::SUCCESS;
  }
    
  //minimal NPV selection
  if(m_NPV < 1 && m_isTLANtupleOffline){
      if(m_debug) cout << " Fail NPV " << endl;
      return EL::StatusCode::SUCCESS;
  }
  bool failLArError = false;
  bool failToolError = false;
    
  //minimal LAr/event cleaning selection
  if (m_LArError && (m_isTLANtupleOffline || m_isTLANtupleTrig)) {
      //if(m_debug) cout << " Fail LArError " << endl;
      //cout << "eventNumber: " << m_eventNumber << " Fail LArError " << endl;
      //return EL::StatusCode::SUCCESS;
      failLArError=true;
  }


  if (m_applyTLALArEventVetoData && !m_isTLANtupleTruth && !m_isDijetNtupleTruth) {
        //if(m_debug) cout << " Fail LArError " << endl;
        // run, lbn, timestamp (seconds), timestamp_ns_offset
        bool veto = m_dataForLArEventVeto->shouldVeto(m_runNumber,m_lumiBlock,m_timeStamp,m_timeStampNSOffset);
        failToolError=veto;
        //if (veto) cout << "eventNumber: " << m_eventNumber << " Fail LArError, event veto " << endl;
        //return EL::StatusCode::SUCCESS;
    }
    
  if (m_applyTLALArEventVetoData) {
    if (failToolError && failLArError) m_h2_LArError->Fill(1.5,1.5);
    if (!failToolError && failLArError) m_h2_LArError->Fill(1.5,0.5);
    if (failToolError && !failLArError) m_h2_LArError->Fill(0.5,1.5);
    if (!failToolError && !failLArError) m_h2_LArError->Fill(0.5,0.5);
  }
    
  if (m_isTLANtupleOffline)
    if (failLArError) return EL::StatusCode::SUCCESS;

  if (m_isTLANtupleTrig && m_applyTLALArEventVetoData)
    if (failToolError) return EL::StatusCode::SUCCESS;

  /*if (m_applySF) {


      TLorentzVector jet1_t = TLorentzVector();
      jet1_t.SetPtEtaPhiE(m_jet_pt->at(0),m_jet_eta->at(0),m_jet_phi->at(0),m_jet_E->at(0));
      TLorentzVector jet2_t = TLorentzVector();
      jet2_t.SetPtEtaPhiE(m_jet_pt->at(1),m_jet_eta->at(1),m_jet_phi->at(1),m_jet_E->at(1));

      int EtaBin = m_hcalibration->GetYaxis()->FindBin(jet1_t.Eta());
      int PTBin = m_hcalibration->GetXaxis()->FindBin(jet1_t.Pt());

      if (jet1_t.Pt() >m_pt_freeze)
        PTBin = m_hcalibration->GetXaxis()->FindBin(m_pt_freeze - 1e-5);
      
      if (fabs(jet1_t.Eta()) >m_eta_freeze){
	    if(jet1_t.Eta() < 0)
		EtaBin = m_hcalibration->GetYaxis()->FindBin(-m_eta_freeze +1e-5);    
     	if(jet1_t.Eta() > 0)
		EtaBin = m_hcalibration->GetYaxis()->FindBin(m_eta_freeze - 1e-5);    
     }

      double scale = m_hcalibration->GetBinContent( PTBin, EtaBin);
      jet1_t = jet1_t * scale;

      EtaBin = m_hcalibration->GetYaxis()->FindBin(jet2_t.Eta());
      PTBin = m_hcalibration->GetXaxis()->FindBin(jet2_t.Pt());
      if (jet2_t.Pt() >m_pt_freeze)
	  PTBin = m_hcalibration->GetXaxis()->FindBin(m_pt_freeze - 1e-5);
      if (fabs(jet2_t.Eta()) > m_eta_freeze){
        if(jet2_t.Eta() < 0)
		  EtaBin = m_hcalibration->GetYaxis()->FindBin(-m_eta_freeze+1e-5);
        if(jet2_t.Eta() > 0)
		  EtaBin = m_hcalibration->GetYaxis()->FindBin(m_eta_freeze-1.e-5);
      }
      
      scale = m_hcalibration->GetBinContent( PTBin, EtaBin);

      jet2_t = jet2_t * scale;

      m_jet_pt->at(0) = jet1_t.Pt();
      m_jet_pt->at(1) = jet2_t.Pt();

      m_jet_eta->at(0) = jet1_t.Eta();
      m_jet_eta->at(1) = jet2_t.Eta();

      m_jet_phi->at(0) = jet1_t.Phi();
      m_jet_phi->at(1) = jet2_t.Phi();

      m_jet_E->at(0) = jet1_t.E();
      m_jet_E->at(1) = jet2_t.E();

  }*/


  // now I'm at the jet selection. Want to split into secondary and non-secondary
  // change to a loop, then return -> skipEvent=true; break;
  
  std::vector<bool> isSecondaryVec;
  isSecondaryVec.clear();
  isSecondaryVec.push_back(false);
  if(m_doSecondaryJets) isSecondaryVec.push_back(true);
  
  
  bool skipEvent = false;
  for(bool isSecondary: isSecondaryVec) {
    
    // veto njets < 2 - earlier I've done an or
    if(isSecondary && nsecJets < 2){
      if(m_debug) cout << " Fail NSecJets " << endl;
      continue;
    }
    else if(!isSecondary && njets < 2){
      if(m_debug) cout << " Fail NJets " << endl;
      continue;
    }

    // set jet_pt etc accordingly
    if(m_debug) {
      cout << "setting up jet collection: ";
      if(isSecondary) cout << "secondary" << endl;
      else cout << "primary" <<endl;
    }
    
    vector<float> *jet_pt = m_jet_pt;
    vector<float> *jet_eta = m_jet_eta;
    vector<int> *jet_clean_passLooseBad = m_jet_clean_passLooseBad;
    if(isSecondary) {
      jet_pt = m_secJet_pt;
      jet_eta = m_secJet_eta;
      // jet_clean_passLooseBad = m_secJet_clean_passLooseBad; // not calculated for secondary jets, but 1:1 correspondance. Might cause crashes in future if not 1:1...
    }
    
    if(m_debug) cout << jet_pt->size() << " jets in collection" << endl;

    if(jet_pt->at(0) < m_leadJetPtCut) {
      if(m_debug) cout << " Fail LeadJetPt " << endl;
      continue;
    }

    if(m_doTrigger){
      if(m_debug) Info("execute()", "Doing Trigger ");

      bool m_dumpTrig = false ;
      if(m_dumpTrig){
	cout << " --------" << endl;
	for(std::string& thisTrig: *m_passedTriggers)
	  cout << thisTrig << endl;
	}//end dump trigger

      std::string trig;
      if (jet_pt->at(0) < 85) continue;
      if (jet_pt->at(0) < 116) { trig = "HLT_j60"; }
      else if (jet_pt->at(0) < 172) { trig = "HLT_j85"; }
      else if (jet_pt->at(0) < 240) { trig = "HLT_j110"; }
      else if (jet_pt->at(0) < 318) { trig = "HLT_j200"; }
      else if (jet_pt->at(0) < 350) { trig = "HLT_j260"; }
      else if (jet_pt->at(0) < 410) { trig = "HLT_j320"; }
      else trig = "HLT_j360";
	std::vector<string>::iterator trigIt = std::find(m_passedTriggers->begin(), m_passedTriggers->end(), trig);
      if (trigIt == m_passedTriggers->end()) continue;
      else {
	prescaleWeight = m_triggerPrescales->at(std::distance(m_passedTriggers->begin(), trigIt));
	//std::cout << "trig: " << trig << ", distance from front of vector" << (std::distance(m_passedTriggers->begin(), trigIt)) << std::endl;
	//std::cout << "prescale: " << prescaleWeight << std::endl;
      }//end do trigger with prescales

      //bool passHLT_j360 = (find(m_passedTriggers->begin(), m_passedTriggers->end(), "HLT_j360" ) != m_passedTriggers->end());

    }// end of do trigger

    if(jet_pt->at(1) < m_subleadJetPtCut) {
      if(m_debug) cout << " Fail subLeadJetPt " << endl;
      continue;
    }

    if(fabs(jet_eta->at(0)) > m_etaCut) {
      if(m_debug) cout << " Fail LeadJetEta " << endl;
      continue;
    }


    if(fabs(jet_eta->at(1)) > m_etaCut) {
      if(m_debug) cout << " Fail subLeadJetEta " << endl;
      continue;
    }

    TLorentzVector jet1 = TLorentzVector();
    if(isSecondary) jet1.SetPtEtaPhiE(m_secJet_pt->at(0),m_secJet_eta->at(0),m_secJet_phi->at(0),m_secJet_E->at(0));
    else            jet1.SetPtEtaPhiE(m_jet_pt->at(0),m_jet_eta->at(0),m_jet_phi->at(0),m_jet_E->at(0));

    TLorentzVector jet2 = TLorentzVector();
    if(isSecondary) jet2.SetPtEtaPhiE(m_secJet_pt->at(1),m_secJet_eta->at(1),m_secJet_phi->at(1),m_secJet_E->at(1));
    else            jet2.SetPtEtaPhiE(m_jet_pt->at(1),m_jet_eta->at(1),m_jet_phi->at(1),m_jet_E->at(1));

    float yStar = ( jet1.Rapidity() - jet2.Rapidity() ) / 2.0;

    if(fabs(yStar) > m_YStarCut){
      if(m_debug) cout << " Fail Ystar " << endl;
      continue;
    }

    float yBoost = ( jet1.Rapidity() + jet2.Rapidity() ) / 2.0;
    if(fabs(yBoost) > m_YBoostCut){
      if(m_debug) cout << " Fail YBoost " << endl;
      continue;
    }

    //
    // doing cleaning on all jets above 50 GeV
    //
    bool  passCleaning   = true;
    if(!m_doTruthOnly){
      for(unsigned int i = 0;  i< njets; ++i){
	if(jet_pt->at(i) > 50){
	  if(fabs(jet_eta->at(i)) < m_etaCut){
	    if(!m_doCleaning && !jet_clean_passLooseBad->at(i)){
	      cout << "Skipping jet " << endl;
	      continue;
	    }

	    if(!jet_clean_passLooseBad->at(i)) passCleaning = false;
	  }
	}else{
	  continue;
	}
      }
    }

    if(m_debug) cout << " Pass All Cut " << endl;

    // get event weight 
    float eventWeight = m_weight;
    if(!m_doData) eventWeight = m_weight * m_lumi/m_sampleEvents;
    if(m_debug) cout << " lumi: " << m_lumi << endl;
    if(m_debug) cout << " weight: " << m_weight << endl;
    if(m_debug) cout << " sampleEvents: " << m_sampleEvents << endl;

    if(m_doData) eventWeight = 1.0;

    //if(m_useWeighted && (eventWeight > 1000)){
    //  cout << "skipping event with Weight: " << eventWeight << " " << m_lumi << " " << m_weight << " " << m_weight_xs << endl;
    //  continue;    
    //}

    //
    //  Jet Cleaning 
    //
    if(m_doCleaning && !passCleaning){
      if (m_debug) cout << "Fail Cleaning " << endl;
      continue;
    }

    if(m_debug) cout << " Weight: " << eventWeight << endl;


    // make event data
    if(m_debug) cout << " Make EventData " << endl;
    eventData thisEvent = eventData(m_runNumber,
				    m_eventNumber,
				    m_jet_pt,
				    m_jet_eta,
				    m_jet_phi,
				    m_jet_E,
				    m_jet_muonSegments,
				    m_jet_EMFrac,
				    m_jet_HECFrac,
				    m_jet_timing,
				    m_jet_negativeE,
				    m_MHT,
				    eventWeight,
				    prescaleWeight,
				    m_avgIntPerX);
    
    if(isSecondary) {
      thisEvent = eventData(m_runNumber,
			    m_eventNumber,
			    m_secJet_pt,
			    m_secJet_eta,
			    m_secJet_phi,
			    m_secJet_E,
			    m_secJet_muonSegments,
			    m_secJet_EMFrac,
			    m_secJet_HECFrac,
			    m_secJet_timing,
			    m_secJet_negativeE,
			    m_MHT, // this isn't filled anyway, currently
			    eventWeight,
			    prescaleWeight,
			    m_avgIntPerX);
    }

  
    // fill mjj-agnostic histograms
    if(isSecondary) {
      hSecIncl->Fill(thisEvent);
      //central: leading jet within 1.2
      if (fabs(m_secJet_eta->at(0))<1.2) hSecCentral->Fill(thisEvent);
      //crack: leading jet within 1.2-1.6
      else if (fabs(m_secJet_eta->at(0))>=1.2 && fabs(m_secJet_eta->at(0))<1.6) hSecCrack->Fill(thisEvent);
      //endcap: leading jet within 1.6-2.8
      else if (fabs(m_secJet_eta->at(0))>=1.6 && fabs(m_secJet_eta->at(0))<2.8) hSecEndcap->Fill(thisEvent);
      if(m_debug) cout << "filled hSecIncl, hSecCentral, hSecCrack and hSecEndcap" << endl;
    }
    else {
      hIncl->Fill(thisEvent);
      if (fabs(m_jet_eta->at(0))<1.2) hCentral->Fill(thisEvent);
      else if (fabs(m_jet_eta->at(0))>=1.2 && fabs(m_jet_eta->at(0))<1.6) hCrack->Fill(thisEvent);
      else if (fabs(m_jet_eta->at(0))>=1.6 && fabs(m_jet_eta->at(0))<2.8) hEndcap->Fill(thisEvent);
      if(m_debug) cout << "filled hIncl, hCentral, hCrack and hEndcap" << endl;
    }




    // mjj isn't in the new ntuples, calculate it here
    float mjj = 0;
    if(m_isDijetNtuple) {
      if(m_debug) cout << "calculating mjj" << endl;
      mjj = (thisEvent.jets.at(0).vec() + thisEvent.jets.at(1).vec()).M();
    }
    else {
      mjj = m_mjj;
    }

    // fill mjj-aware histograms
    if(isSecondary) {
      if ( mjj>394 && mjj<1236 ) {
	hSecIncl_mjjWindow->Fill(thisEvent);
	//central: leading jet within 1.2
	if (fabs(m_secJet_eta->at(0))<1.2) hSecCentral_mjjWindow->Fill(thisEvent);
	//crack: leading jet within 1.2-1.6
	else if (fabs(m_secJet_eta->at(0))>=1.2 && fabs(m_secJet_eta->at(0))<1.6) hSecCrack_mjjWindow->Fill(thisEvent);
	//endcap: leading jet within 1.6-2.8
	else if (fabs(m_secJet_eta->at(0))>=1.6 && fabs(m_secJet_eta->at(0))<2.8) hSecEndcap_mjjWindow->Fill(thisEvent);
      }
    }
    else {
      if ( mjj>394 && mjj<1236 ) {
	hIncl_mjjWindow->Fill(thisEvent);
	if (fabs(m_jet_eta->at(0))<1.2) hCentral_mjjWindow->Fill(thisEvent);
	else if (fabs(m_jet_eta->at(0))>=1.2 && fabs(m_jet_eta->at(0))<1.6) hCrack_mjjWindow->Fill(thisEvent);
	else if (fabs(m_jet_eta->at(0))>=1.6 && fabs(m_jet_eta->at(0))<2.8) hEndcap_mjjWindow->Fill(thisEvent);
      }
    }

    //hOffline->Fill(thisEvent_offline);
    //hTrigger->Fill(thisEvent_trigger);

  }

  return EL::StatusCode::SUCCESS;
}







EL::StatusCode ProcessTLAMiniTree :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessTLAMiniTree :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  //this makes ROOT angry
  //delete m_dataForLArEventVeto;
    
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessTLAMiniTree :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}


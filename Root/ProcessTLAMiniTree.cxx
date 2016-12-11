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

  m_applyLArEventCleaning(true),
  m_invertLArEventCleaning(false),
  m_applyTLALArEventVetoData(false),
  m_TLALArEventVetoFiles("$ROOTCOREBIN/data/TLAEventCleaning/event-veto-data/"),
  m_doCleaning(false),
  m_recalculateJetCleaning(false),
  m_invertJetCleaning(false),

  m_debug(false),
  m_doTrigger(false),
  m_doTrigger_j110(false),
  m_isDijetNtupleOffline(false),
  m_isDijetNtupleTrig(false),
  m_isDijetNtupleTruth(false),
  m_isTLANtupleTruth(false),
  m_isTLANtupleTrig(false),
  m_isTLANtupleOffline(false),
  m_doTruthOnly(false),

  m_primaryJetInName("jet"),
  m_primaryJetOutName("OfflineJets"),
  m_doSecondaryJets(false),
  m_secondaryJetInName("trigJet"),
  m_secondaryJetOutName("TriggerJets"),
  m_doPileupFromMap(false),
  m_pileupMap(""),
  m_avgIntPerX(-1),
  m_avgIntPerX_fromMap(-1),
  m_avgIntPerX_fromAOD(-1),

  m_YStarCut(99),
  m_YBoostCut(99),
  m_doData(false),
  m_useWeighted(false),
  m_etaCut(1.8),
  m_leadJetPtCut(200.),
  m_subleadJetPtCut(200.),
  m_lumi(0),
  m_sampleEvents(0),
  m_grl(nullptr),

  m_plotPtSlices(false),
  m_plotEtaSlices(false),
  m_plotMjjWindow(false),

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
  m_jet_clean_passLooseBad(0),
  m_jet_clean_passLooseBad_recalc(0),

  m_jet_LArQuality(0),
  m_jet_AverageLArQF(0),
  m_jet_HECQuality(0),
  m_jet_FracSamplingMax(0),
  m_jet_FracSamplingMaxIndex(0),
  m_jet_LeadingClusterPt(0),
  m_jet_LeadingClusterSecondLambda(0),
  m_jet_LeadingClusterCenterLambda(0),
  m_jet_LeadingClusterSecondR(0),

  m_MHT(0),
  m_mjj(0),
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
  hPt200(nullptr),
  hPt210(nullptr),
  hPt220(nullptr),
  hPt230(nullptr),
  hPt240(nullptr),
  hPt250(nullptr),

  hSecIncl(nullptr),
  hSecCentral(nullptr),
  hSecCrack(nullptr),
  hSecEndcap(nullptr),
  hSecIncl_mjjWindow(nullptr),
  hSecCentral_mjjWindow(nullptr),
  hSecCrack_mjjWindow(nullptr),
  hSecEndcap_mjjWindow(nullptr),
  hSecPt200(nullptr),
  hSecPt210(nullptr),
  hSecPt220(nullptr),
  hSecPt230(nullptr),
  hSecPt240(nullptr),
  hSecPt250(nullptr),

  m_secJet_pt(0),
  m_secJet_eta(0),
  m_secJet_phi(0),
  m_secJet_E(0),
  m_secJet_muonSegments(0),
  m_secJet_EMFrac(0),
  m_secJet_HECFrac(0),
  m_secJet_timing(0),
  m_secJet_negativeE(0),
  m_secJet_clean_passLooseBad(0),
  m_secJet_clean_passLooseBad_recalc(0),
  m_secJet_LArQuality(0),
  m_secJet_AverageLArQF(0),
  m_secJet_HECQuality(0),
  m_secJet_FracSamplingMax(0),
  m_secJet_FracSamplingMaxIndex(0),
  m_secJet_LeadingClusterPt(0),
  m_secJet_LeadingClusterSecondLambda(0),
  m_secJet_LeadingClusterCenterLambda(0),
  m_secJet_LeadingClusterSecondR(0)


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
  else if (m_isDijetNtupleTrig || m_isDijetNtupleOffline) {
   hIncl    = new eventHists((m_primaryJetOutName+"").c_str()        ,       wk());
   if (m_plotEtaSlices) {
     hCentral = new eventHists((m_primaryJetOutName+"_0-12").c_str()   ,       wk());
     hCrack   = new eventHists((m_primaryJetOutName+"_12-16").c_str()  ,       wk());
     hEndcap  = new eventHists((m_primaryJetOutName+"_16-28").c_str()  ,       wk());
   }
   if(m_plotMjjWindow) {
     hIncl_mjjWindow    = new eventHists((m_primaryJetOutName+"_mjjWindow").c_str()        ,       wk());
     if (m_plotEtaSlices) {
       hCentral_mjjWindow = new eventHists((m_primaryJetOutName+"_0-12_mjjWindow").c_str()   ,       wk());
       hCrack_mjjWindow   = new eventHists((m_primaryJetOutName+"_12-16_mjjWindow").c_str()  ,       wk());
       hEndcap_mjjWindow  = new eventHists((m_primaryJetOutName+"_16-28_mjjWindow").c_str()  ,       wk());
     }
   }
   if(m_plotPtSlices) {
     hPt200 = new eventHists((m_primaryJetOutName+"_j0Pt200").c_str()  ,       wk());
     hPt210 = new eventHists((m_primaryJetOutName+"_j0Pt210").c_str()  ,       wk());
     hPt220 = new eventHists((m_primaryJetOutName+"_j0Pt220").c_str()  ,       wk());
     hPt230 = new eventHists((m_primaryJetOutName+"_j0Pt230").c_str()  ,       wk());
     hPt240 = new eventHists((m_primaryJetOutName+"_j0Pt240").c_str()  ,       wk());
     hPt250 = new eventHists((m_primaryJetOutName+"_j0Pt250").c_str()  ,       wk());
   }

   if(m_doSecondaryJets) {
     hSecIncl    = new eventHists((m_secondaryJetOutName+"").c_str()        ,       wk());
     if (m_plotEtaSlices) {
       hSecCentral = new eventHists((m_secondaryJetOutName+"_0-12").c_str()   ,       wk());
       hSecCrack   = new eventHists((m_secondaryJetOutName+"_12-16").c_str()  ,       wk());
       hSecEndcap  = new eventHists((m_secondaryJetOutName+"_16-28").c_str()  ,       wk());
     }
     if (m_plotMjjWindow) {
       hSecIncl_mjjWindow    = new eventHists((m_secondaryJetOutName+"_mjjWindow").c_str()        ,       wk());
       if (m_plotEtaSlices) {
         hSecCentral_mjjWindow = new eventHists((m_secondaryJetOutName+"_0-12_mjjWindow").c_str()   ,       wk());
         hSecCrack_mjjWindow   = new eventHists((m_secondaryJetOutName+"_12-16_mjjWindow").c_str()  ,       wk());
         hSecEndcap_mjjWindow  = new eventHists((m_secondaryJetOutName+"_16-28_mjjWindow").c_str()  ,       wk());
       }
     }
     if (m_plotPtSlices) {
       hSecPt200 = new eventHists((m_secondaryJetOutName+"_j0Pt200").c_str()  ,       wk());
       hSecPt210 = new eventHists((m_secondaryJetOutName+"_j0Pt210").c_str()  ,       wk());
       hSecPt220 = new eventHists((m_secondaryJetOutName+"_j0Pt220").c_str()  ,       wk());
       hSecPt230 = new eventHists((m_secondaryJetOutName+"_j0Pt230").c_str()  ,       wk());
       hSecPt240 = new eventHists((m_secondaryJetOutName+"_j0Pt240").c_str()  ,       wk());
       hSecPt250 = new eventHists((m_secondaryJetOutName+"_j0Pt250").c_str()  ,       wk());
     }
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

  // validation histograms
  m_h2_LArError = new TH2D("h2_LArError", "h2_LArError;Offline: isLArError;Tool: isLArError", 2, 0, 2, 4, 0, 4);
  m_h2_LArError->GetYaxis()->SetBinLabel(1,"None");
  m_h2_LArError->GetYaxis()->SetBinLabel(2,"NoiseBurst");
  m_h2_LArError->GetYaxis()->SetBinLabel(3,"MiniNoiseBurst");
  m_h2_LArError->GetYaxis()->SetBinLabel(4,"Other/Both");
  wk()->addOutput(m_h2_LArError);

  m_h2_LArError_postSelection = new TH2D("h2_LArError_postSelection", "h2_LArError_postSelection;Offline: isLArError;Tool: isLArError", 2, 0, 2, 4, 0, 4);
  m_h2_LArError_postSelection->GetYaxis()->SetBinLabel(1,"None");
  m_h2_LArError_postSelection->GetYaxis()->SetBinLabel(2,"NoiseBurst");
  m_h2_LArError_postSelection->GetYaxis()->SetBinLabel(3,"MiniNoiseBurst");
  m_h2_LArError_postSelection->GetYaxis()->SetBinLabel(4,"Other/Both");
  wk()->addOutput(m_h2_LArError_postSelection);

  m_h2_LArError_postSelection_w = new TH2D("h2_LArError_postSelection_w", "h2_LArError_postSelection_weighted;Offline: isLArError;Tool: isLArError", 2, 0, 2, 4, 0, 4);
  m_h2_LArError_postSelection_w->GetYaxis()->SetBinLabel(1,"None");
  m_h2_LArError_postSelection_w->GetYaxis()->SetBinLabel(2,"NoiseBurst");
  m_h2_LArError_postSelection_w->GetYaxis()->SetBinLabel(3,"MiniNoiseBurst");
  m_h2_LArError_postSelection_w->GetYaxis()->SetBinLabel(4,"Other/Both");
  wk()->addOutput(m_h2_LArError_postSelection_w);


  m_h2_avgIntPerX_map_AOD = new TH2D("h2_avgIntPerX_map_AOD", "h2_avgIntPerX_map_AOD;from map;from AOD", 101,-1,100,101,-1,100);
  wk()->addOutput(m_h2_avgIntPerX_map_AOD);
    
  m_h2_jetCleaning = new TH2D("h2_jetCleaning", "h2_jetCleaning;jet pass LooseBad;manual calc", 2, 0, 2, 2, 0, 2);
  wk()->addOutput(m_h2_jetCleaning);

  // cutflow histograms
  m_h_cutflow_primary = new TH1D("h_cutflow_primary", "h_cutflow_primary", 1, 0, 1);
  m_h_cutflow_primary->SetBit(TH1::kCanRebin); // since I want to extend range with more entries in cutflow
  wk()->addOutput(m_h_cutflow_primary);

  m_h_cutflow_secondary = new TH1D("h_cutflow_secondary", "h_cutflow_secondary", 1, 0, 1);
  m_h_cutflow_secondary->SetBit(TH1::kCanRebin);
  wk()->addOutput(m_h_cutflow_secondary);

  m_h_cutflow_primary_w = new TH1D("h_cutflow_primary_w", "h_cutflow_primary_weighted", 1, 0, 1);
  m_h_cutflow_primary_w->SetBit(TH1::kCanRebin);
  wk()->addOutput(m_h_cutflow_primary_w);

  m_h_cutflow_secondary_w = new TH1D("h_cutflow_secondary_w", "h_cutflow_secondary_weighted", 1, 0, 1);
  m_h_cutflow_secondary_w->SetBit(TH1::kCanRebin);
  wk()->addOutput(m_h_cutflow_secondary_w);


  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessTLAMiniTree :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  Info("fileExecute()", "Calling fileExecute"); // CWK debugging  

  TFile* inputFile = wk()->inputFile();
  cout << inputFile << ", " << inputFile->GetName() << endl;

  // get nEvents values from NTUP and set bin errors appropriately
  TH1D* NTUP_MetaData = (TH1D*)inputFile->Get("MetaData_EventCount");
  m_h_cutflow_primary->Fill("NTUP initial", NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("nEvents initial")));
  m_h_cutflow_primary->Fill("NTUP selected", NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("nEvents selected")));
  m_h_cutflow_primary->SetBinError(1, sqrt(NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("nEvents initial"))) );
  m_h_cutflow_primary->SetBinError(2, sqrt(NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("nEvents selected"))) );
  // weighted
  m_h_cutflow_primary_w->Fill("NTUP initial", NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("sumOfWeights initial")));
  m_h_cutflow_primary_w->Fill("NTUP selected", NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("sumOfWeights selected")));
  m_h_cutflow_primary->SetBinError(1, sqrt(NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("sumOfWeightsSquared initial"))) );
  m_h_cutflow_primary->SetBinError(2, sqrt(NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("sumOfWeightsSquared selected"))) );

  // secondary hists
  m_h_cutflow_secondary->Fill("NTUP initial", NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("nEvents initial")));
  m_h_cutflow_secondary->Fill("NTUP selected", NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("nEvents selected")));
  m_h_cutflow_secondary->SetBinError(1, sqrt(NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("nEvents initial"))) );
  m_h_cutflow_secondary->SetBinError(2, sqrt(NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("nEvents selected"))) );
  // weighted
  m_h_cutflow_secondary_w->Fill("NTUP initial", NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("sumOfWeights initial")));
  m_h_cutflow_secondary_w->Fill("NTUP selected", NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("sumOfWeights selected")));
  m_h_cutflow_secondary->SetBinError(1, sqrt(NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("sumOfWeightsSquared initial"))) );
  m_h_cutflow_secondary->SetBinError(2, sqrt(NTUP_MetaData->GetBinContent(NTUP_MetaData->GetXaxis()->FindBin("sumOfWeightsSquared selected"))) );


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

    if(m_doTrigger || m_doTrigger_j110){
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
  else if (m_isDijetNtupleTrig || m_isDijetNtupleOffline) {
    
      tree->SetBranchStatus  ((m_primaryJetInName+"_pt").c_str(), 1);
      tree->SetBranchAddress ((m_primaryJetInName+"_pt").c_str(), &m_jet_pt);
      
      tree->SetBranchStatus  ((m_primaryJetInName+"_eta").c_str(), 1);
      tree->SetBranchAddress ((m_primaryJetInName+"_eta").c_str(), &m_jet_eta);
      
      tree->SetBranchStatus  ((m_primaryJetInName+"_phi").c_str(), 1);
      tree->SetBranchAddress ((m_primaryJetInName+"_phi").c_str(), &m_jet_phi);
      
      tree->SetBranchStatus  ((m_primaryJetInName+"_E").c_str(), 1);
      tree->SetBranchAddress ((m_primaryJetInName+"_E").c_str(), &m_jet_E);

      if(!m_doTruthOnly) {
	// old NTUPs- prior to Oct 2016
	// tree->SetBranchStatus  ("jet__EMFrac", 1);
	// tree->SetBranchAddress ("jet__EMFrac", &m_jet_EMFrac);
	
	// tree->SetBranchStatus  ("jet__HECFrac", 1);
	// tree->SetBranchAddress ("jet__HECFrac", &m_jet_HECFrac);
	
	// tree->SetBranchStatus  ("jet__GhostMuonSegmentCount", 1);
	// tree->SetBranchAddress ("jet__GhostMuonSegmentCount", &m_jet_muonSegments);

	// new NTUPs - from 15.10.16
	tree->SetBranchStatus  ((m_primaryJetInName+"_EMFrac").c_str(), 1);
	tree->SetBranchAddress ((m_primaryJetInName+"_EMFrac").c_str(), &m_jet_EMFrac);
	
	tree->SetBranchStatus  ((m_primaryJetInName+"_HECFrac").c_str(), 1);
	tree->SetBranchAddress ((m_primaryJetInName+"_HECFrac").c_str(), &m_jet_HECFrac);
	
	tree->SetBranchStatus  ((m_primaryJetInName+"_GhostMuonSegmentCount").c_str(), 1);
	tree->SetBranchAddress ((m_primaryJetInName+"_GhostMuonSegmentCount").c_str(), &m_jet_muonSegments);

	tree->SetBranchStatus  ((m_primaryJetInName+"_clean_passLooseBad").c_str(), 1);
	tree->SetBranchAddress ((m_primaryJetInName+"_clean_passLooseBad").c_str(), &m_jet_clean_passLooseBad);	

	tree->SetBranchStatus  ((m_primaryJetInName+"_Timing").c_str(), 1);
	tree->SetBranchAddress ((m_primaryJetInName+"_Timing").c_str(), &m_jet_timing);

	tree->SetBranchStatus  ((m_primaryJetInName+"_NegativeE").c_str(), 1);
	tree->SetBranchAddress ((m_primaryJetInName+"_NegativeE").c_str(), &m_jet_negativeE);

        tree->SetBranchStatus  ((m_primaryJetInName+"_LArQuality").c_str(), 1);
        tree->SetBranchAddress ((m_primaryJetInName+"_LArQuality").c_str(), &m_jet_LArQuality);
        tree->SetBranchStatus  ((m_primaryJetInName+"_AverageLArQF").c_str(), 1);
        tree->SetBranchAddress ((m_primaryJetInName+"_AverageLArQF").c_str(), &m_jet_AverageLArQF);
        tree->SetBranchStatus  ((m_primaryJetInName+"_HECQuality").c_str(), 1);
        tree->SetBranchAddress ((m_primaryJetInName+"_HECQuality").c_str(), &m_jet_HECQuality);
        tree->SetBranchStatus  ((m_primaryJetInName+"_FracSamplingMax").c_str(), 1);
        tree->SetBranchAddress ((m_primaryJetInName+"_FracSamplingMax").c_str(), &m_jet_FracSamplingMax);
        tree->SetBranchStatus  ((m_primaryJetInName+"_FracSamplingMaxIndex").c_str(), 1);
        tree->SetBranchAddress ((m_primaryJetInName+"_FracSamplingMaxIndex").c_str(), &m_jet_FracSamplingMaxIndex);
        tree->SetBranchStatus  ((m_primaryJetInName+"_LeadingClusterPt").c_str(), 1);
        tree->SetBranchAddress ((m_primaryJetInName+"_LeadingClusterPt").c_str(), &m_jet_LeadingClusterPt);
        tree->SetBranchStatus  ((m_primaryJetInName+"_LeadingClusterSecondLambda").c_str(), 1);
	tree->SetBranchAddress ((m_primaryJetInName+"_LeadingClusterSecondLambda").c_str(), &m_jet_LeadingClusterSecondLambda);
	tree->SetBranchStatus  ((m_primaryJetInName+"_LeadingClusterCenterLambda").c_str(), 1);
	tree->SetBranchAddress ((m_primaryJetInName+"_LeadingClusterCenterLambda").c_str(), &m_jet_LeadingClusterCenterLambda);
	tree->SetBranchStatus  ((m_primaryJetInName+"_LeadingClusterSecondR").c_str(), 1);
	tree->SetBranchAddress ((m_primaryJetInName+"_LeadingClusterSecondR").c_str(), &m_jet_LeadingClusterSecondR);

      }

      if(m_isDijetNtupleOffline) {
	tree->SetBranchStatus  ("averageInteractionsPerCrossing", 1);
	tree->SetBranchAddress ("averageInteractionsPerCrossing", &m_avgIntPerX_fromAOD);
      }


      if(m_doSecondaryJets) {
	tree->SetBranchStatus  ((m_secondaryJetInName+"_pt").c_str(), 1);
	tree->SetBranchAddress ((m_secondaryJetInName+"_pt").c_str(), &m_secJet_pt);
	
	tree->SetBranchStatus  ((m_secondaryJetInName+"_eta").c_str(), 1);
	tree->SetBranchAddress ((m_secondaryJetInName+"_eta").c_str(), &m_secJet_eta);
	
	tree->SetBranchStatus  ((m_secondaryJetInName+"_phi").c_str(), 1);
	tree->SetBranchAddress ((m_secondaryJetInName+"_phi").c_str(), &m_secJet_phi);
	
	tree->SetBranchStatus  ((m_secondaryJetInName+"_E").c_str(), 1);
	tree->SetBranchAddress ((m_secondaryJetInName+"_E").c_str(), &m_secJet_E);
	
	if(!m_doTruthOnly) {
	  // new NTUPs - from 15.10.16
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_EMFrac").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_EMFrac").c_str(), &m_secJet_EMFrac);
	  
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_HECFrac").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_HECFrac").c_str(), &m_secJet_HECFrac);
	  
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_GhostMuonSegmentCount").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_GhostMuonSegmentCount").c_str(), &m_secJet_muonSegments);
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_Timing").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_Timing").c_str(), &m_secJet_timing);
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_NegativeE").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_NegativeE").c_str(), &m_secJet_negativeE);
	  
	  // below this they are dodgy (irrespective of whe secondaryJetInName is)
          tree->SetBranchStatus  ((m_secondaryJetInName+"_clean_passLooseBad").c_str(), 1);
          tree->SetBranchAddress ((m_secondaryJetInName+"_clean_passLooseBad").c_str(), &m_secJet_clean_passLooseBad);
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_LArQuality").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_LArQuality").c_str(), &m_secJet_LArQuality);
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_AverageLArQF").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_AverageLArQF").c_str(), &m_secJet_AverageLArQF);

	  tree->SetBranchStatus  ((m_secondaryJetInName+"_HECQuality").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_HECQuality").c_str(), &m_secJet_HECQuality);
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_FracSamplingMax").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_FracSamplingMax").c_str(), &m_secJet_FracSamplingMax);
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_FracSamplingMaxIndex").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_FracSamplingMaxIndex").c_str(), &m_secJet_FracSamplingMaxIndex);
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_LeadingClusterPt").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_LeadingClusterPt").c_str(), &m_secJet_LeadingClusterPt);
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_LeadingClusterSecondLambda").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_LeadingClusterSecondLambda").c_str(), &m_secJet_LeadingClusterSecondLambda);
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_LeadingClusterCenterLambda").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_LeadingClusterCenterLambda").c_str(), &m_secJet_LeadingClusterCenterLambda);
	  tree->SetBranchStatus  ((m_secondaryJetInName+"_LeadingClusterSecondR").c_str(), 1);
	  tree->SetBranchAddress ((m_secondaryJetInName+"_LeadingClusterSecondR").c_str(), &m_secJet_LeadingClusterSecondR);

	}
	
      }

  }

  if(m_debug) {
    cout << "\nlist of activated branches:" << endl;
    for(int ob=0; ob < tree->GetListOfBranches()->GetEntries(); ob++ ) {
      
      TBranch* br = (TBranch*)tree->GetListOfBranches()->At(ob);
      TString brname = br->GetName();
      if ( tree->GetBranchStatus(brname) == 1 )
	cout << brname << endl;
    }
    cout << endl;
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
    cout << "I am loading LArEventVeto files from " << LArEventVetoExpandedPath << endl;
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

  m_avgIntPerX = -1;

  string cutname = "executed";
  m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
  m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(m_weight));
  m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
  m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(m_weight));

 
  if(m_doData){
   
    if ( m_applyGRL ) {
      if ( !m_grl->passRunLB( m_runNumber, m_lumiBlock ) ) {
	if(m_debug) cout << "GRL:: Fail Event " << endl;
        return EL::StatusCode::SUCCESS; // go to next event
      }else{
	if(m_debug) cout << "GRL:: Pass Event: runnum, LB = " << m_runNumber << ", " << m_lumiBlock << endl;
      }
    }

    if(m_doPileupFromMap) {
      // std::cout << "RunNum: " << m_runNumber << ", LB: " << m_lumiBlock << std::endl;
      int runNumberBin = m_h2_pileupMap->GetYaxis()->FindBin(to_string(m_runNumber).c_str());
      int lumiBlockBin = m_h2_pileupMap->GetXaxis()->FindBin(m_lumiBlock);
      // std::cout << "RN bin: " << runNumberBin << ", LB bin: " << lumiBlockBin << std::endl;
      m_avgIntPerX_fromMap = m_h2_pileupMap->GetBinContent(lumiBlockBin, runNumberBin);
      // std::cout << "m_avgIntPerX: " << m_avgIntPerX << std::endl;

      m_h2_avgIntPerX_map_AOD->Fill(m_avgIntPerX_fromMap, m_avgIntPerX_fromAOD);
      m_avgIntPerX = m_avgIntPerX_fromMap;
    }

    cutname = "pass GRL";
    m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
    m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(m_weight));
    m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
    m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(m_weight));

  }//end of m_doData

  // take avgIntPerX from the AOD if it exists, otherwise go for the 'from map' version (which will be -1 if unset)
  if(m_avgIntPerX_fromAOD != -1) {
    m_avgIntPerX = m_avgIntPerX_fromAOD;
  }



  //////////////////////////////
  // LAr event cleaning start //
  //////////////////////////////

  int LArError_offline = 0;
  int LArError_tool = 0;
  bool failLArError = false;
  bool failToolError = false;

  //minimal LAr/event cleaning selection
  if(!m_isTLANtupleTruth && !m_isDijetNtupleTruth) {
    if (m_LArError) {
      LArError_offline = 1;
      failLArError = true;
    }

    if (m_applyTLALArEventVetoData) {
        //if(m_debug) cout << " Fail LArError " << endl;
        // run, lbn, timestamp (seconds), timestamp_ns_offset
      bool veto = m_dataForLArEventVeto->shouldVeto(m_runNumber,m_lumiBlock,m_timeStamp,m_timeStampNSOffset);
      if(veto) {
	failToolError = true;
	LArError_tool = 1;
	string LArErrorType = m_dataForLArEventVeto->vetoType(m_runNumber,m_lumiBlock,m_timeStamp,m_timeStampNSOffset);
	if (LArErrorType=="NoiseBurst")          LArError_tool = 1;
	else if (LArErrorType=="MiniNoiseBurst") LArError_tool = 2;
	else                                     LArError_tool = 3;
      }
      m_h2_LArError->Fill(LArError_offline, LArError_tool);
    }
    
  }

  // apply decision of tool, if I'm applying it
  if (m_applyLArEventCleaning) {
    if (m_invertLArEventCleaning) {
      if( ( (m_isTLANtupleOffline || m_isDijetNtupleOffline) && !failLArError) || 
	  ( (m_isTLANtupleTrig || m_isDijetNtupleTrig) && !failToolError) ) {
	if(m_debug) cout << " Pass LAr event cleaning" << endl;
	return EL::StatusCode::SUCCESS;
      }
    }
    else {
      if( ( (m_isTLANtupleOffline || m_isDijetNtupleOffline) && failLArError) || 
	  ( (m_isTLANtupleTrig || m_isDijetNtupleTrig) && failToolError) ) {
	if(m_debug) cout << " Fail LAr event cleaning" << endl;
	return EL::StatusCode::SUCCESS;
      }
    }
  }

  cutname = "LAr event cleaning";
  m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
  m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(m_weight));
  m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
  m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(m_weight));

  //////////////////////////////
  // LAr event cleaning end   //
  //////////////////////////////



  //minimal NPV selection
  if(m_NPV < 1 && (m_isTLANtupleOffline || m_isDijetNtupleOffline)){
      if(m_debug) cout << " Fail NPV " << endl;
      return EL::StatusCode::SUCCESS;
  }

  cutname = "NPV";
  m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
  m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(m_weight));
  m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
  m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(m_weight));
    


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
    
    // veto njets < 2
    if(isSecondary && nsecJets < 2){
      if(m_debug) cout << " Fail NSecJets " << endl;
      continue;
    }
    else if(!isSecondary && njets < 2){
      if(m_debug) cout << " Fail NJets " << endl;
      continue;
    }

    cutname = "2 jets";
    if(isSecondary) {
      m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(m_weight));
    } else {
      m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(m_weight));
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
      jet_clean_passLooseBad = m_secJet_clean_passLooseBad;
    }
    
    if(m_debug) cout << jet_pt->size() << " jets in collection" << endl;


    // select on relevant triggers (want to base trigger choice on relevant jet collection)
    // (this way it doesn't matter whether I run together or separately)

    if(m_doTrigger || m_doTrigger_j110){
      if(m_debug) Info("execute()", "Doing Trigger ");
      
      bool m_dumpTrig = false ;
      if(m_dumpTrig){
	cout << " --------" << endl;
	for(std::string& thisTrig: *m_passedTriggers)
	  cout << thisTrig << endl;
      }//end dump trigger
      
      std::string trig;
      
      // this is from 2015
      // if (m_jet_pt->at(0) < 85) return EL::StatusCode::SUCCESS; 
      // if (m_jet_pt->at(0) < 116) { trig = "HLT_j60"; }
      // else if (m_jet_pt->at(0) < 172) { trig = "HLT_j85"; }
      // else if (m_jet_pt->at(0) < 240) { trig = "HLT_j110"; }
      // else if (m_jet_pt->at(0) < 318) { trig = "HLT_j200"; }
      // else if (m_jet_pt->at(0) < 350) { trig = "HLT_j260"; }
      // else if (m_jet_pt->at(0) < 410) { trig = "HLT_j320"; }
      // else trig = "HLT_j360";
      
      
      // taking triggers active for all of 2016
      // j110 21
      // j150 26 <- ignore since barely any less prescaled than j110
      // j175 154
      // j260 950
      // j380 38074
      // full 38409
      
      // linear fit to 2015 turnons: t = (1.03*p + 36)
      // tried 10% increase but looked a tad worse re not-on-threshold
      
      if (m_jet_pt->at(0) < 216) { trig = "HLT_j110"; }
      else if (m_jet_pt->at(0) < 304) { trig = "HLT_j175"; }
      else if (m_jet_pt->at(0) < 427) { trig = "HLT_j260"; }
      else trig = "HLT_j380";
      
      if(m_doTrigger_j110)
	trig = "HLT_j110";
      
      std::vector<string>::iterator trigIt = std::find(m_passedTriggers->begin(), m_passedTriggers->end(), trig);
      if (trigIt == m_passedTriggers->end()) {
	continue;
      }
      else {
	prescaleWeight = m_triggerPrescales->at(std::distance(m_passedTriggers->begin(), trigIt));
	//std::cout << "trig: " << trig << ", distance from front of vector" << (std::distance(m_passedTriggers->begin(), trigIt)) << std::endl;
	//std::cout << "prescale: " << prescaleWeight << std::endl;
      }//end do trigger with prescales
      
      //bool passHLT_j360 = (find(m_passedTriggers->begin(), m_passedTriggers->end(), "HLT_j360" ) != m_passedTriggers->end());
      
      cutname = "trigger";
      if(isSecondary) {
	m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
	m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(m_weight));
      } else {
	m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
	m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(m_weight));
      }


      cutname = "trigger prescale";
      if(isSecondary) {
	m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(m_weight*prescaleWeight));
      } else {
	m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(m_weight*prescaleWeight));
      }
      
    }// end of do trigger
    

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

    if(m_debug) cout << " Weight: " << eventWeight << endl;

    cutname = "final weight";
    if(isSecondary) {
      m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    } else {
      m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    }

    
    // kinematic selection
    if(jet_pt->at(0) < m_leadJetPtCut) {
      if(m_debug) cout << " Fail LeadJetPt " << endl;
      continue;
    }

    cutname = "lead jet pT > "+to_string(m_leadJetPtCut);
    if(isSecondary) {
      m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    } else {
      m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    }


    if(jet_pt->at(1) < m_subleadJetPtCut) {
      if(m_debug) cout << " Fail subLeadJetPt " << endl;
      continue;
    }

    cutname = "sublead jet pT > "+to_string(m_subleadJetPtCut);
    if(isSecondary) {
      m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    } else {
      m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    }


    if(fabs(jet_eta->at(0)) > m_etaCut) {
      if(m_debug) cout << " Fail LeadJetEta " << endl;
      continue;
    }

    cutname = "lead jet eta < "+to_string(m_etaCut);
    if(isSecondary) {
      m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    } else {
      m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    }


    if(fabs(jet_eta->at(1)) > m_etaCut) {
      if(m_debug) cout << " Fail subLeadJetEta " << endl;
      continue;
    }

    cutname = "sublead jet eta < "+to_string(m_etaCut);
    if(isSecondary) {
      m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    } else {
      m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
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

    cutname = "yStar < "+to_string(m_YStarCut);
    if(isSecondary) {
      m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    } else {
      m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    }

    float yBoost = ( jet1.Rapidity() + jet2.Rapidity() ) / 2.0;
    if(fabs(yBoost) > m_YBoostCut){
      if(m_debug) cout << " Fail YBoost " << endl;
      continue;
    }

    cutname = "yBoost < "+to_string(m_YBoostCut);
    if(isSecondary) {
      m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    } else {
      m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    }


    /////////////////////////////////////////////////
    // doing cleaning on all jets above **60** GeV //
    /////////////////////////////////////////////////

    // first, set all other variables according to primary / secondary
    // jet_pt set earlier
    // jet_eta set earlier
    vector<float>*	jet_phi	                 = ( isSecondary ? m_secJet_phi          : m_jet_phi );
    vector<float>*	jet_E			 = ( isSecondary ? m_secJet_E            : m_jet_E );
    vector<float>*	jet_muonSegments	 = ( isSecondary ? m_secJet_muonSegments : m_jet_muonSegments );
    vector<float>*	jet_EMFrac		 = ( isSecondary ? m_secJet_EMFrac       : m_jet_EMFrac );
    vector<float>*	jet_HECFrac		 = ( isSecondary ? m_secJet_HECFrac      : m_jet_HECFrac );
    vector<float>*	jet_timing		 = ( isSecondary ? m_secJet_timing       : m_jet_timing );
    vector<float>*	jet_negativeE		 = ( isSecondary ? m_secJet_negativeE    : m_jet_negativeE );
    // jet_clean_passLooseBad set earlier (and is maybe about to be reset)
    vector<float>*	jet_LArQuality	               = ( isSecondary ? m_secJet_LArQuality                 : m_jet_LArQuality );
    vector<float>*	jet_AverageLArQF	       = ( isSecondary ? m_secJet_AverageLArQF               : m_jet_AverageLArQF );
    vector<float>*	jet_HECQuality		       = ( isSecondary ? m_secJet_HECQuality                 : m_jet_HECQuality );
    vector<float>*	jet_FracSamplingMax	       = ( isSecondary ? m_secJet_FracSamplingMax            : m_jet_FracSamplingMax );
    vector<int>*	jet_FracSamplingMaxIndex       = ( isSecondary ? m_secJet_FracSamplingMaxIndex       : m_jet_FracSamplingMaxIndex );
    vector<float>*	jet_LeadingClusterPt           = ( isSecondary ? m_secJet_LeadingClusterPt           : m_jet_LeadingClusterPt );
    vector<float>*	jet_LeadingClusterSecondLambda = ( isSecondary ? m_secJet_LeadingClusterSecondLambda : m_jet_LeadingClusterSecondLambda );
    vector<float>*	jet_LeadingClusterCenterLambda = ( isSecondary ? m_secJet_LeadingClusterCenterLambda : m_jet_LeadingClusterCenterLambda );
    vector<float>*	jet_LeadingClusterSecondR      = ( isSecondary ? m_secJet_LeadingClusterSecondR      : m_jet_LeadingClusterSecondR );


    if(m_debug) cout << " about to do cleaning" << endl;
    bool  passCleaning   = true;

    // redo cleaning??
    bool redo = m_recalculateJetCleaning;
    if(!redo) { // check it's actually filled, if not redo
      if(jet_clean_passLooseBad->at(0) != 0 && jet_clean_passLooseBad->at(0) != 1)
	redo = true;
    }
    if(m_debug) {
      cout << "jet_clean_passLooseBad->at(0) = " << jet_clean_passLooseBad->at(0) << endl;
      cout << "redo is " << redo << endl;
    }

    // now redo if requested / required
    if(redo) {
      jet_clean_passLooseBad->clear();
      
      for(unsigned int i = 0;  i < jet_pt->size(); ++i){

	// need to apply chf for jets below 60 -> we should ignore them
	if(jet_pt->at(i) < 60) {
	  jet_clean_passLooseBad->push_back(0);
	  continue;
	}

	bool jetPassCleaning = true;
	if(m_debug) cout << "  welcome to the cleaning of jet "<<i<<endl;
	// if (m_jet_clean_passLooseBad->at(i)==0)
	  // cout << "jet "<<i<<" failed cleaning according to AOD value";

	
	if ( jet_HECFrac->at(i) > 0.5 && 
	     fabs(jet_HECQuality->at(i)) > 0.5 && 
	     jet_AverageLArQF->at(i) > 0.8 ) {
	  if(m_debug) cout << "    HECfrac failed" << endl;
	  jetPassCleaning = false;
	}
	  
	if ( fabs(jet_negativeE->at(i)) > 60 ) {
	  if(m_debug) cout << "    negativeE failed" << endl;
	  jetPassCleaning = false;
	}

	if ( jet_EMFrac->at(i) > 0.95 && 
	     jet_LArQuality->at(i) > 0.8 && 
	     jet_AverageLArQF->at(i) > 0.8 && 
	     fabs(jet_eta->at(i)) < 2.8 ) {
	  if(m_debug) cout << "    EMfrac failed" << endl;
	  jetPassCleaning = false;
	}
	
	if ( jet_FracSamplingMax->at(i) > 0.99 &&
	     fabs(jet_eta->at(i)) < 2.0 ) {
	  if(m_debug) cout << "    high fracSamplingMax failed" << endl;
	  jetPassCleaning = false;
	}
	    
	// if (jet_EMFrac->at(i) < 0.05 &&
	    // jet_chFrac->at(i) < 0.05 && 
	    // fabs(jet_eta->at(i)) < 2.0 ) {
	  // if(m_debug) cout << "    low fracSamplingMax (low eta plus low chf) failed" << endl;
	  // jetPassCleaning = false;
	// }
	    
	if ( jet_EMFrac->at(i) < 0.05 && 
	     fabs(jet_eta->at(i)) >= 2.0 ) {
	  if(m_debug) cout << "    low fracSamplingMax (high eta) failed" << endl;
	  jetPassCleaning = false;
	}

	if (!jetPassCleaning) {
	  jet_clean_passLooseBad->push_back(0);
	  if(m_debug) cout <<"  did not pass recalcuated cleaning" << endl;
	}
	else {
	  if(m_debug) cout <<"  passed recalcuated cleaning" << endl;
	  jet_clean_passLooseBad->push_back(1);
	}

	// if(jet_clean_passLooseBad->at(i) != m_jet_clean_passLooseBad->at(i) ) {
	  // cout << "  WARNING in fill loop! Jet " << i << " has discrepant cleaning decision... original says passed = "<< m_jet_clean_passLooseBad->at(i) << " but recalculated version says " << jet_clean_passLooseBad->at(i) << endl;
	// }	
      }

    }

    // assign to the relevant part
    if (isSecondary) m_secJet_clean_passLooseBad_recalc = jet_clean_passLooseBad;
    else             m_jet_clean_passLooseBad_recalc = jet_clean_passLooseBad;


    // fill comparison histogram, print out decision if different
    // don't do for secondary since I know that's broken (at least, secJets in NTUP don't have it - you could have swapped them in the HIST config)
    if(!isSecondary) {

      for (unsigned int i = 0;  i < jet_pt->size(); ++i){ 
	if(jet_pt->at(i) < 60) continue;

	m_h2_jetCleaning->Fill(m_jet_clean_passLooseBad->at(i), m_jet_clean_passLooseBad_recalc->at(i));
	
	if( m_jet_clean_passLooseBad_recalc->at(i) != m_jet_clean_passLooseBad->at(i) ) {
	  cout << "WARNING! Jet " << i << " has discrepant cleaning decision... original says passed = "<< m_jet_clean_passLooseBad->at(i) << " but recalculated version says " << m_jet_clean_passLooseBad_recalc->at(i) << endl;
	  
	  cout << "here are the relevant cleaning variables..." << endl;
	  
	  cout << "  HECfrac: \t" << m_jet_HECFrac->at(i) << endl;
	  cout << "  abs(HECQual): " << fabs(m_jet_HECQuality->at(i)) << endl; 
	  cout << "  <LArQF>: \t" << m_jet_AverageLArQF->at(i) << endl;
	  cout << "  abs(negE): \t" << fabs(m_jet_negativeE->at(i)) << endl;
	  cout << "  EMfrac: \t" << m_jet_EMFrac->at(i) << endl;
	  cout << "  LARQual: \t" << m_jet_LArQuality->at(i) << endl;
	  cout << "  fracMax: \t" << m_jet_FracSamplingMax->at(i) << endl;
	  cout << "  abs(eta): \t" << fabs(m_jet_eta->at(i)) << endl;
	  cout << "  pT: \t\t" << fabs(m_jet_pt->at(i)) << endl;
	  
	  cout << "What about the subsets?" << endl;
	  
	  cout << "  HECfrac? \t" << 
	    ( m_jet_HECFrac->at(i) > 0.5 && 
	      fabs(m_jet_HECQuality->at(i)) > 0.5 && 
	      m_jet_AverageLArQF->at(i) > 0.8 ) << endl;
	  
	  cout << "  negE? \t" <<
	    ( fabs(m_jet_negativeE->at(i)) > 60 ) << endl;
	  cout << "  EMfrac? \t" << 
	    ( m_jet_EMFrac->at(i) > 0.95 && 
	      m_jet_LArQuality->at(i) > 0.8 && 
	      m_jet_AverageLArQF->at(i) > 0.8 && 
	      fabs(m_jet_eta->at(i)) < 2.8 ) << endl;
	  cout << "  fracMax? \t" <<
	    ( m_jet_FracSamplingMax->at(i) > 0.99 &&
	      fabs(m_jet_eta->at(i)) < 2.0 ) << endl;
	  cout << "  low EMfrac? \t" <<
	    ( m_jet_EMFrac->at(i) < 0.05 && 
	      fabs(m_jet_eta->at(i)) >= 2.0 ) << endl;
	  cout << "  with CHF? \t" << 
	    (m_jet_EMFrac->at(i) < 0.05 && fabs(m_jet_eta->at(i)) < 2.0 ) << endl; // if true then if chf<0.05 would have been flagged as bad
	}
      }
    }
    

    if(!m_doTruthOnly){
      for(unsigned int i = 0;  i < jet_pt->size(); ++i){ // change from njets since that is hard-coded from primary 
	if (i > 2) continue; // If we only care about cleaning the leading three jets?
	if(jet_pt->at(i) > 60){ // threshold for JVT has increased to 60 for 2016
	  if(fabs(jet_eta->at(i)) < m_etaCut){
	    if(!jet_clean_passLooseBad->at(i)) {
	      passCleaning = false;
	      if(!m_doCleaning) {
		cout << "not doing cleaning, but jet "<< i <<" failed passLooseBad" << endl;
	      }
	    }
	  }
	}
      }
    }

    
    
    //
    //  Jet Cleaning 
    //
    if(m_doCleaning) {
      if ( !m_invertJetCleaning && !passCleaning ) {
	if (m_debug) cout << "Fail jet cleaning " << endl;
	continue;
      }
      if ( m_invertJetCleaning && passCleaning ) {
	if (m_debug) cout << "Pass jet cleaning (but I'm inverting, so skip event) " << endl;
	continue;
      }
    }

    cutname = "jet cleaning";
    if(isSecondary) {
      m_h_cutflow_secondary   -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_secondary_w -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    } else {
      m_h_cutflow_primary     -> Fill(cutname.c_str(), 1.0);
      m_h_cutflow_primary_w   -> Fill(cutname.c_str(), static_cast<double>(eventWeight*prescaleWeight));
    }


    if(m_debug) cout << " Pass All Cut " << endl;


    // fill LAr Event cleaning histo after event selection cuts
    if(!isSecondary) {
      m_h2_LArError_postSelection->Fill(LArError_offline, LArError_tool);
      m_h2_LArError_postSelection_w->Fill(LArError_offline, LArError_tool, eventWeight*prescaleWeight);
    }


    // make event data
    if(m_debug) cout << " filling eventData" << endl;
    eventData thisEvent = eventData(m_runNumber,
				    m_eventNumber,
				    jet_pt,
				    jet_eta,
				    jet_phi,
				    jet_E,
				    jet_muonSegments,
				    jet_EMFrac,
				    jet_HECFrac,
				    jet_timing,
				    jet_negativeE,
				    jet_clean_passLooseBad,
				    jet_LArQuality,
				    jet_AverageLArQF,
				    jet_HECQuality,
				    jet_FracSamplingMax,
				    jet_FracSamplingMaxIndex,
				    jet_LeadingClusterPt,
				    jet_LeadingClusterSecondLambda,
				    jet_LeadingClusterCenterLambda,
				    jet_LeadingClusterSecondR,
				    m_MHT,
				    eventWeight,
				    prescaleWeight,
				    m_avgIntPerX);
    
    if(m_debug) cout << "done extracting event data to structs" << endl;
  
    // fill mjj-agnostic histograms
    if(isSecondary) {
      hSecIncl->Fill(thisEvent);
      if(m_debug) cout << "filled hSecIncl" << endl;
      if(m_plotEtaSlices) {
        //central: leading jet within 1.2
        if (fabs(m_secJet_eta->at(0))<1.2) hSecCentral->Fill(thisEvent);
        //crack: leading jet within 1.2-1.6
        else if (fabs(m_secJet_eta->at(0))>=1.2 && fabs(m_secJet_eta->at(0))<1.6) hSecCrack->Fill(thisEvent);
        //endcap: leading jet within 1.6-2.8
        else if (fabs(m_secJet_eta->at(0))>=1.6 && fabs(m_secJet_eta->at(0))<2.8) hSecEndcap->Fill(thisEvent);
        if(m_debug) cout << "filled hSecCentral, hSecCrack and hSecEndcap" << endl;
      }
      if(m_plotPtSlices) {
        if (m_secJet_pt->at(0) > 200) hSecPt200->Fill(thisEvent);
        if (m_secJet_pt->at(0) > 210) hSecPt210->Fill(thisEvent);
        if (m_secJet_pt->at(0) > 220) hSecPt220->Fill(thisEvent);
        if (m_secJet_pt->at(0) > 230) hSecPt230->Fill(thisEvent);
        if (m_secJet_pt->at(0) > 240) hSecPt240->Fill(thisEvent);
        if (m_secJet_pt->at(0) > 250) hSecPt250->Fill(thisEvent);
      }
    }
    else {
      hIncl->Fill(thisEvent);
      if(m_debug) cout << "filled hIncl" << endl;
      if(m_plotEtaSlices) {
        if (fabs(m_jet_eta->at(0))<1.2) hCentral->Fill(thisEvent);
        else if (fabs(m_jet_eta->at(0))>=1.2 && fabs(m_jet_eta->at(0))<1.6) hCrack->Fill(thisEvent);
        else if (fabs(m_jet_eta->at(0))>=1.6 && fabs(m_jet_eta->at(0))<2.8) hEndcap->Fill(thisEvent);
        if(m_debug) cout << "filled hCentral, hCrack and hEndcap" << endl;
      }
      if(m_plotPtSlices) {
        if (m_jet_pt->at(0) > 200) hPt200->Fill(thisEvent);
        if (m_jet_pt->at(0) > 210) hPt210->Fill(thisEvent);
        if (m_jet_pt->at(0) > 220) hPt220->Fill(thisEvent);
        if (m_jet_pt->at(0) > 230) hPt230->Fill(thisEvent);
        if (m_jet_pt->at(0) > 240) hPt240->Fill(thisEvent);
        if (m_jet_pt->at(0) > 250) hPt250->Fill(thisEvent);
      }
    }


    // mjj isn't in the new ntuples, calculate it here
    float mjj = 0;
    if(m_isDijetNtupleOffline || m_isDijetNtupleTrig || m_isTLANtupleTruth) {
      if(m_debug) cout << "calculating mjj" << endl;
      mjj = (thisEvent.jets.at(0).vec() + thisEvent.jets.at(1).vec()).M();
    }
    else {
      mjj = m_mjj;
    }

    // fill mjj-aware histograms
    if(isSecondary) {
      if(m_plotMjjWindow) {
        if ( mjj>394 && mjj<1236 ) {
          hSecIncl_mjjWindow->Fill(thisEvent);
          if(m_plotEtaSlices) {
            //central: leading jet within 1.2
            if (fabs(m_secJet_eta->at(0))<1.2) hSecCentral_mjjWindow->Fill(thisEvent);
            //crack: leading jet within 1.2-1.6
            else if (fabs(m_secJet_eta->at(0))>=1.2 && fabs(m_secJet_eta->at(0))<1.6) hSecCrack_mjjWindow->Fill(thisEvent);
            //endcap: leading jet within 1.6-2.8
            else if (fabs(m_secJet_eta->at(0))>=1.6 && fabs(m_secJet_eta->at(0))<2.8) hSecEndcap_mjjWindow->Fill(thisEvent);
          }
        }
      }
    }
    else {
      if(m_plotMjjWindow){
        if ( mjj>394 && mjj<1236 ) {
          hIncl_mjjWindow->Fill(thisEvent);
          if(m_plotEtaSlices) {
            if (fabs(m_jet_eta->at(0))<1.2) hCentral_mjjWindow->Fill(thisEvent);
            else if (fabs(m_jet_eta->at(0))>=1.2 && fabs(m_jet_eta->at(0))<1.6) hCrack_mjjWindow->Fill(thisEvent);
            else if (fabs(m_jet_eta->at(0))>=1.6 && fabs(m_jet_eta->at(0))<2.8) hEndcap_mjjWindow->Fill(thisEvent);
          }
        }
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


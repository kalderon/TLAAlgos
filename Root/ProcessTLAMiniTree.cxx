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
  m_applyGRL(true),
  m_GRLxml("$ROOTCOREBIN/data/xAODAnaHelpers/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml"),
  m_debug(false),
  m_doTrigger(false),
  m_doTruthOnly(false),
  m_detailLevel(0),
  m_YStarCut(99),
  m_YBoostCut(99),
  m_doBlind(false),
  m_doData(false),
  m_doPUReweight(false),
  m_useWeighted(false),
  m_doCleaning(false),
  m_etaCut(2.8),
  m_leadJetPtCut(200.),
  m_lumi(0),
  m_sampleEvents(0),
  m_grl(nullptr),
  m_jet_pt(0),
  m_jet_eta(0),
  m_jet_phi(0),
  m_jet_E(0),
  m_jet_Timing(0),
  m_jet_GhostMuonSegmentCount(0),
  m_jet_clean_passLooseBad(0),
  m_passedTriggers(nullptr),
  m_gamma_pt(0),
  m_gamma_eta(0),
  m_gamma_phi(0),
  m_gamma_E(0),
  hIncl(nullptr),
  hInclYStar(nullptr),
  hMjj200_400(nullptr),
  hMjj200_400_Ystar(nullptr)
{
  Info("ProcessTLAMiniTree()", "Calling constructor");
  m_applyGRL      = true;
  m_GRLxml   = "$ROOTCOREBIN/data/xAODAnaHelpers/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml";  //https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/GoodRunListsForAnalysis
}


EL::StatusCode  ProcessTLAMiniTree :: configure ()
{

//  //
//  // Read Input from .config file
//  //
//  m_debug                    = config->GetValue("Debug" ,          m_debug );
//  m_doTrigger                = config->GetValue("DoTrigger" ,      m_doTrigger );
//  m_detailLevel              = config->GetValue("DetailLevel" ,    m_detailLevel );
//  m_YStarCut                 = config->GetValue("YStarCut" ,       m_YStarCut);
//  m_YBoostCut                = config->GetValue("YBoostCut" ,      m_YBoostCut);
//  m_doBlind                  = config->GetValue("DoBlind" ,        m_doBlind );
//  m_doData                   = config->GetValue("DoData" ,         m_doData );
//  m_doPUReweight             = config->GetValue("DoPUReweight" ,   m_doPUReweight );
//  m_useWeighted              = config->GetValue("UseWeighted" ,    m_useWeighted );
//  m_doCleaning               = config->GetValue("DoCleaning" ,     m_doCleaning );
//  m_etaCut                   = config->GetValue("AbsEtaCut" ,      m_etaCut );
//  m_leadJetPtCut             = config->GetValue("LeadJetPtCut" ,   m_leadJetPtCut );
//  m_lumi                     = config->GetValue("Lumi",            m_lumi );
//  m_applyGRL                 = config->GetValue("ApplyGRL",        m_applyGRL);
//  m_GRLxml                   = config->GetValue("GRL",             m_GRLxml.c_str());

  hIncl          = new eventHists("Incl"  ,       wk(), m_etaCut, m_detailLevel);
  hMjj200_400    = new eventHists("Mjj200_400"  , wk(), m_etaCut, m_detailLevel);

  hInclYStar           = new eventHists("InclYStar"  ,        wk(), m_etaCut, m_detailLevel);
  hMjj200_400_Ystar    = new eventHists("Mjj200_400_YStar"  , wk(), m_etaCut, m_detailLevel);
  
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

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessTLAMiniTree :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ProcessTLAMiniTree :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  if(firstFile){
    if ( this->configure() == EL::StatusCode::FAILURE ) {
      Error("initialize()", "Failed to properly configure. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    TFile* inputFile = wk()->inputFile();
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
	  cout << "Setting Sample events to: " << m_sampleEvents << endl;
	  cout << "Setting Lumi  to: " << m_lumi << endl;
	  
	  std::string inFileName = inputFile->GetName();
	  std::size_t foundJZ5Pos = inFileName.find("JZ5");
	  bool foundJZ5 = (foundJZ5Pos!=std::string::npos);
	  if(foundJZ5){
	    m_sampleEvents = (m_sampleEvents - 352219559.60316455) ;
	    cout << "Correcting Sample events to: " << m_sampleEvents << endl;
	  }
	  
	} else if(!m_useWeighted && !foundWeighted) {
	  cout << "Getting NSample events from " << keyName << endl;
	  m_sampleEvents = ((TH1F*)key->ReadObj())->GetBinContent(1);
	  cout << "Setting Sample events to: " << m_sampleEvents << endl;
	  cout << "Setting Lumi  to: " << m_lumi << endl;
	}

      }//Found CF
    }//over Keys
    
  }// first file

  TTree *tree = wk()->tree();
  tree->SetBranchStatus ("*", 0);

  tree->SetBranchStatus  ("runNumber",    1);
  tree->SetBranchAddress ("runNumber",    &m_runNumber);

  tree->SetBranchStatus  ("eventNumber",    1);
  tree->SetBranchAddress ("eventNumber",    &m_eventNumber);

  if(!m_doTruthOnly){
    tree->SetBranchStatus  ("lumiBlock",    1);
    tree->SetBranchAddress ("lumiBlock",    &m_lumiBlock);

    tree->SetBranchStatus  ("NPV",    1);
    tree->SetBranchAddress ("NPV",    &m_NPV);

    tree->SetBranchStatus  ("actualInteractionsPerCrossing",    1);
    tree->SetBranchAddress ("actualInteractionsPerCrossing",    &m_actualInteractionsPerCrossing);

    tree->SetBranchStatus  ("averageInteractionsPerCrossing",    1);
    tree->SetBranchAddress ("averageInteractionsPerCrossing",    &m_averageInteractionsPerCrossing);

    tree->SetBranchStatus  ("weight_pileup", 1);
    tree->SetBranchAddress ("weight_pileup", &m_weight_pileup);

    if(m_detailLevel > 0){
      tree->SetBranchStatus  ("jet_Timing", 1);
      tree->SetBranchAddress ("jet_Timing", &m_jet_Timing);

      tree->SetBranchStatus  ("jet_GhostMuonSegmentCount", 1);
      tree->SetBranchAddress ("jet_GhostMuonSegmentCount", &m_jet_GhostMuonSegmentCount);
    }

    tree->SetBranchStatus  ("jet_clean_passLooseBad", 1);
    tree->SetBranchAddress ("jet_clean_passLooseBad", &m_jet_clean_passLooseBad);

    if(m_doTrigger){
      tree->SetBranchStatus  ("passedTriggers", 1);
      tree->SetBranchAddress ("passedTriggers", &m_passedTriggers);
    }
  }

  tree->SetBranchStatus  ("weight", 1);
  tree->SetBranchAddress ("weight", &m_weight);

  tree->SetBranchStatus  ("weight_xs", 1);
  tree->SetBranchAddress ("weight_xs", &m_weight_xs);


  tree->SetBranchStatus  ("jet_pt", 1);
  tree->SetBranchAddress ("jet_pt", &m_jet_pt);

  tree->SetBranchStatus  ("jet_eta", 1);
  tree->SetBranchAddress ("jet_eta", &m_jet_eta);

  tree->SetBranchStatus  ("jet_phi", 1);
  tree->SetBranchAddress ("jet_phi", &m_jet_phi);

  tree->SetBranchStatus  ("jet_E", 1);
  tree->SetBranchAddress ("jet_E", &m_jet_E);

  std::string photonName = "ph";
  if(m_doTruthOnly) photonName = "photons";

  tree->SetBranchStatus  ((photonName+"_pt" ).c_str(), 1);
  tree->SetBranchAddress ((photonName+"_pt" ).c_str(), &m_gamma_pt);

  tree->SetBranchStatus  ((photonName+"_eta").c_str(), 1);
  tree->SetBranchAddress ((photonName+"_eta").c_str(), &m_gamma_eta);

  tree->SetBranchStatus  ((photonName+"_phi").c_str(), 1);
  tree->SetBranchAddress ((photonName+"_phi").c_str(), &m_gamma_phi);
  
  tree->SetBranchStatus  ((photonName+"_E"  ).c_str(), 1);
  tree->SetBranchAddress ((photonName+"_E"  ).c_str(), &m_gamma_E);



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

  if(m_applyGRL){
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    std::vector<std::string> vecStringGRL;
    m_GRLxml = gSystem->ExpandPathName( m_GRLxml.c_str() );
    vecStringGRL.push_back(m_GRLxml);
    RETURN_CHECK("BasicEventSelection::initialize()", m_grl->setProperty( "GoodRunsListVec", vecStringGRL), "");
    RETURN_CHECK("BasicEventSelection::initialize()", m_grl->setProperty("PassThrough", false), "");
    RETURN_CHECK("BasicEventSelection::initialize()", m_grl->initialize(), "");
  }

  Info("initialize()", "Succesfully initialized! \n");
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

  if(m_doData){
   
    if ( m_applyGRL ) {
      if ( !m_grl->passRunLB( m_runNumber, m_lumiBlock ) ) {
	if(m_debug) cout << "GRL:: Fail Event " << endl;
        return EL::StatusCode::SUCCESS; // go to next event
      }else{
	//if(m_debug) cout << "GRL:: Pass Event " << endl;
      }
    }
  }

  //
  // Minimun selection is a dijet 
  //
  if(njets < 2){
    if(m_debug) cout << " Fail NJets " << endl;
    return EL::StatusCode::SUCCESS;
  }

  if(m_jet_pt->at(0) < m_leadJetPtCut) {
    if(m_debug) cout << " Fail LeadJetPt " << endl;
    return EL::StatusCode::SUCCESS;
  }


  TLorentzVector jet1 = TLorentzVector();
  jet1.SetPtEtaPhiE(m_jet_pt->at(0),m_jet_eta->at(0),m_jet_phi->at(0),m_jet_E->at(0));
  TLorentzVector jet2 = TLorentzVector();
  jet2.SetPtEtaPhiE(m_jet_pt->at(1),m_jet_eta->at(1),m_jet_phi->at(1),m_jet_E->at(1));
  float yStar = ( jet1.Rapidity() - jet2.Rapidity() ) / 2.0;
  
  if(fabs(yStar) > m_YStarCut){
    if(m_debug) cout << " Fail Ystar " << endl;
    return EL::StatusCode::SUCCESS;
  }

  float yBoost = ( jet1.Rapidity() + jet2.Rapidity() ) / 2.0;
  if(fabs(yBoost) > m_YBoostCut){
    if(m_debug) cout << " Fail YBoost " << endl;
    return EL::StatusCode::SUCCESS;
  }

  float mjj   = ( jet1 +  jet2 ).M();
//  if(mjj < m_mjjCut){
//    if(m_debug) cout << " Fail Mjj " << endl;
//    return EL::StatusCode::SUCCESS;
//  }


  //
  // doing cleaning
  //
  bool  passCleaning   = true;
  if(!m_doTruthOnly){
    for(unsigned int i = 0;  i< njets; ++i){
      if(m_jet_pt->at(i) > 50){
	if(fabs(m_jet_eta->at(i)) < m_etaCut){
	  if(!m_doCleaning && !m_jet_clean_passLooseBad->at(i)){
	    cout << "Skipping jet " << endl;
	    continue;
	  }

	  if(!m_jet_clean_passLooseBad->at(i)) passCleaning = false;
	}
      }else{
	break;
      }
    }
  }

  if(m_debug) cout << " Pass All Cut " << endl;
  float eventWeight = m_weight;
  
  if(!m_doData && m_doPUReweight) eventWeight *= m_weight_pileup;

  if(m_doData) eventWeight = 1.0;

  if(m_useWeighted && (eventWeight > 1000)){
    cout << "skipping event with Weight: " << eventWeight << " " << m_lumi << " " << m_weight << " " << m_weight_xs << endl;
    return EL::StatusCode::SUCCESS;    
  }

  //
  //  Jet Cleaning 
  //
  if(m_doCleaning && !passCleaning){
    //if(!passCleaning){
    cout << "Fail Cleaning " << endl;
    return EL::StatusCode::SUCCESS;
  }
  
  if(m_debug) cout << " Make EventData " << endl;
  eventData thisEvent = eventData(m_runNumber, m_eventNumber, 
				  m_jet_pt, m_jet_eta, m_jet_phi, m_jet_E,  m_jet_Timing, m_jet_GhostMuonSegmentCount,
				  m_gamma_pt, m_gamma_eta, m_gamma_phi, m_gamma_E,  
				  m_NPV, m_averageInteractionsPerCrossing, m_averageInteractionsPerCrossing, 
				  eventWeight);
  if(m_debug) cout << " Made EventDAta " << endl;				  
  
  hIncl->Fill(thisEvent);

  if(fabs(yStar) < 0.75) hInclYStar->Fill(thisEvent);
  if(m_debug) cout << " Filled inc.  event weight " << eventWeight << endl;				  

  if(mjj < 400 && mjj > 200){
    hMjj200_400->Fill(thisEvent);
    if(fabs(yStar) < 0.75) hMjj200_400_Ystar->Fill(thisEvent);  
  }

  if(m_doTrigger){
    if(m_debug) Info("execute()", "Doing Trigger ");

    bool m_dumpTrig = false ;
    if(m_dumpTrig){
      cout << " --------" << endl;
      for(std::string& thisTrig: *m_passedTriggers)
	cout << thisTrig << endl;
    }

    //bool passHLT_j360 = (find(m_passedTriggers->begin(), m_passedTriggers->end(), "HLT_j360" ) != m_passedTriggers->end());

  }// do trigger



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


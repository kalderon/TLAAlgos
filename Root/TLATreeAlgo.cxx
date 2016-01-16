#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>

#include <TLAAlgos/TLATreeAlgo.h>
#include <TLAAlgos/TLATreeHelper.h>

#include <xAODAnaHelpers/Algorithm.h>
#include <xAODAnaHelpers/HelperFunctions.h>
#include <xAODAnaHelpers/HelperClasses.h>
#include <xAODAnaHelpers/tools/ReturnCheck.h>

#include "TEnv.h"
#include "TSystem.h"

using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::ifstream;
using std::istringstream;


// this is needed to distribute the algorithm to the workers
ClassImp(TLATreeAlgo)

TLATreeAlgo :: TLATreeAlgo (std::string className) :
    TreeAlgo(className)
{
  this->SetName("TLATreeAlgo"); // needed if you want to retrieve this algo with wk()->getAlg(ALG_NAME) downstream
  //TODO: get from config file? Do the same in TreeAlgo and in TLATreeAlgo
  m_jetName = "";
  m_trigJetName = "";
  m_truthJetName = "";
    
  std::string m_jetContainerName;
  std::string m_truthJetContainerName;
  std::string m_trigJetContainerName;


}

//EL::StatusCode TLATreeAlgo :: setupJob (EL::Job& job)
//{
//    
//  job.useXAOD();
//  xAOD::Init("TLATreeAlgo").ignore();
//
//  EL::OutputStream outForTree("tree");
//  job.outputAdd (outForTree);
//
//  return EL::StatusCode::SUCCESS;
//}
//
EL::StatusCode TLATreeAlgo :: initialize ()
{
  Info("initialize()", m_name.c_str());
    
  m_eventCounter = -1;


  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  this->treeInitialize(); // better idea: call the base class
    
  //already fill the lumi weights for the first event
  const xAOD::EventInfo* eventInfo(nullptr);
  RETURN_CHECK("ResonanceAlgorithm::initialize()", HelperFunctions::retrieve(eventInfo, "EventInfo", m_event, m_store, m_debug), "");
  m_isMC = ( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) ? true : false;

  this->getLumiWeights(eventInfo);
    
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TLATreeAlgo :: treeInitialize ()
{
  Info("treeInitialize()", "%s", m_name.c_str() );
  // needed here and not in initalize since this is called first
  Info("treeInitialize()", "Attempting to configure using: %s", getConfig().c_str());

  TTree * outTree = new TTree(m_name.c_str(),m_name.c_str());
  if ( !outTree ) {
    Error("treeInitialize()","Failed to instantiate output tree!");
    return EL::StatusCode::FAILURE;
  }

  // get the input from user which determines which branches are created!
  if ( this->configure() != EL::StatusCode::SUCCESS ) {
    Error("treeInitialize()", "%s failed to properly configure. Exiting.", m_name.c_str() );
    return EL::StatusCode::FAILURE;
  } else {
    Info("treeInitialize()", "Succesfully configured! ");
  }

  // get the file we created already
  TFile* treeFile = wk()->getOutputFile ("tree");

  // CD: should delete this? Also m_units is hardwired in here so maybe not much point in having a getter in HelpTreeBase?
  if (!m_helpTree) {
    m_helpTree = new TLATreeHelper( m_event, outTree, treeFile, 1e3, m_debug, m_DC14 );
  }

  // tell the tree to go into the file
  outTree->SetDirectory( treeFile );
  // choose if want to add tree to same directory as ouput histograms
  if ( m_outHistDir ) {
    wk()->addOutput( outTree );
  }

  m_helpTree->AddEvent( m_evtDetailStr );

  if ( !m_trigDetailStr.empty() )       {   m_helpTree->AddTrigger    (m_trigDetailStr);    }
  if ( !m_muContainerName.empty() )     {   m_helpTree->AddMuons      (m_muDetailStr);      }
  if ( !m_elContainerName.empty() )     {   m_helpTree->AddElectrons  (m_elDetailStr);      }
  if ( !m_jetContainerName.empty() )    {
      m_helpTree->AddJets       (m_jetDetailStr, "jet");
      m_jetName = "jet";
  }
  if ( !m_trigJetContainerName.empty() ){
      m_helpTree->AddJets       (m_trigJetDetailStr, "trigJet");
      m_trigJetName = "trigJet";
  }
  if ( !m_truthJetContainerName.empty() ){
      m_helpTree->AddJets       (m_truthJetDetailStr, "truthJet");
      m_truthJetName = "truthJet";
  }
  if ( !m_fatJetContainerName.empty() ) {   m_helpTree->AddFatJets    (m_fatJetDetailStr);  }
  if ( !m_tauContainerName.empty() )    {   m_helpTree->AddTaus       (m_tauDetailStr);     }
  if ( !m_METContainerName.empty() )    {   m_helpTree->AddMET        (m_METDetailStr);     }
  if ( !m_photonContainerName.empty() ) {   m_helpTree->AddPhotons    (m_photonDetailStr);  }

  Info("treeInitialize()", "Successfully initialized output tree");

  return EL::StatusCode::SUCCESS;
}

//Function originally by J. Dandoy, from DijetResonanceAlgo
//This grabs cross section, acceptance, and eventNumber information from the respective text file
//text format:     147915 2.3793E-01 5.0449E-03 499000
EL::StatusCode TLATreeAlgo::getLumiWeights(const xAOD::EventInfo* eventInfo) {
    
    if(!m_isMC){
        m_mcChannelNumber = eventInfo->runNumber();
        m_xs = 1;
        m_filtEff = 1;
        m_numAMIEvents = 0;
        return EL::StatusCode::SUCCESS;
    }
    
    m_mcChannelNumber = eventInfo->mcChannelNumber();
    //if mcChannelNumber = 0 need to retrieve from runNumber
    if(eventInfo->mcChannelNumber()==0) m_mcChannelNumber = eventInfo->runNumber();
    //CD: fixme, hardwired
    ifstream fileIn(  gSystem->ExpandPathName( "$ROOTCOREBIN/data/TLAAlgos/XsAcc_13TeV.txt") );
    std::string runNumStr = std::to_string( m_mcChannelNumber );
    std::string line;
    std::string subStr;
    while (getline(fileIn, line)){
        istringstream iss(line);
        iss >> subStr;
        if (subStr.find(runNumStr) != string::npos){
            iss >> subStr;
            sscanf(subStr.c_str(), "%e", &m_xs);
            iss >> subStr;
            sscanf(subStr.c_str(), "%e", &m_filtEff);
            iss >> subStr;
            sscanf(subStr.c_str(), "%i", &m_numAMIEvents);
            cout << "Setting xs / acceptance / numAMIEvents to " << m_xs << ":" << m_filtEff << ":" << m_numAMIEvents << endl;
            continue;
        }
    }
    if( m_numAMIEvents == 0){
        cerr << "ERROR: Could not find proper file information for file number " << runNumStr << endl;
        return EL::StatusCode::FAILURE;
    }
    return EL::StatusCode::SUCCESS;
}

//%%%%%%%%functions coming straight from TreeAlgo

//EL::StatusCode TLATreeAlgo :: configure ()
//{
//  if (!getConfig().empty()) {
//
//    // the file exists, use TEnv to read it off
//    TEnv* config = new TEnv(getConfig(true).c_str());
//    m_evtDetailStr            = config->GetValue("EventDetailStr",       m_evtDetailStr.c_str());
//    m_trigDetailStr           = config->GetValue("TrigDetailStr",        m_trigDetailStr.c_str());
//    m_trigJetDetailStr        = config->GetValue("TrigJetDetailStr",     m_trigJetDetailStr.c_str());
//    m_muDetailStr             = config->GetValue("MuonDetailStr",        m_muDetailStr.c_str());
//    m_elDetailStr             = config->GetValue("ElectronDetailStr",    m_elDetailStr.c_str());
//    m_truthJetDetailStr       = config->GetValue("TruthJetDetailStr",    m_truthJetDetailStr.c_str());
//    m_jetDetailStr            = config->GetValue("JetDetailStr",         m_jetDetailStr.c_str());
//    m_fatJetDetailStr         = config->GetValue("FatJetDetailStr",      m_fatJetDetailStr.c_str());
//    m_tauDetailStr            = config->GetValue("TauDetailStr",         m_tauDetailStr.c_str());
//    m_METDetailStr            = config->GetValue("METDetailStr",         m_METDetailStr.c_str());
//    m_photonDetailStr         = config->GetValue("PhotonDetailStr",      m_photonDetailStr.c_str());
//
//    m_debug                   = config->GetValue("Debug" ,           m_debug);
//
//    m_outHistDir              = config->GetValue("SameHistsOutDir",  m_outHistDir);
//
//    m_muContainerName         = config->GetValue("MuonContainerName",       m_muContainerName.c_str());
//    m_elContainerName         = config->GetValue("ElectronContainerName",   m_elContainerName.c_str());
//    m_jetContainerName        = config->GetValue("JetContainerName",        m_jetContainerName.c_str());
//    m_truthJetContainerName   = config->GetValue("TruthJetContainerName",   m_truthJetContainerName.c_str());
//    m_trigJetContainerName    = config->GetValue("TrigJetContainerName",    m_trigJetContainerName.c_str());
//    m_fatJetContainerName     = config->GetValue("FatJetContainerName",     m_fatJetContainerName.c_str());
//    m_tauContainerName        = config->GetValue("TauContainerName",        m_tauContainerName.c_str());
//    m_METContainerName        = config->GetValue("METContainerName",        m_METContainerName.c_str());
//    m_photonContainerName     = config->GetValue("PhotonContainerName",     m_photonContainerName.c_str());
//
//    // DC14 switch for little things that need to happen to run
//    // for those samples with the corresponding packages
//    m_DC14                    = config->GetValue("DC14", m_DC14);
//
//    Info("configure()", "Loaded in configuration values");
//
//    // everything seems preliminarily ok, let's print config and say we were successful
//    config->Print();
//
//    delete config; config = nullptr;
//  }
//
//  return EL::StatusCode::SUCCESS;
//}

//EL::StatusCode TLATreeAlgo :: fileExecute () { return EL::StatusCode::SUCCESS; }
//EL::StatusCode TLATreeAlgo :: changeInput (bool /*firstFile*/) { return EL::StatusCode::SUCCESS; }

EL::StatusCode TLATreeAlgo::getJetVariables(std::string jetName, const xAOD::JetContainer* inJets, const xAOD::EventInfo* eventInfo) {

    //this only if 2 or more jets
    if (inJets->size() >= 2) {
        const xAOD::Jet* leadJet     = inJets->at(0);
        const xAOD::Jet* subLeadJet  = inJets->at(1);
        
        float mjj = ( leadJet->p4() + subLeadJet->p4() ).M();
        eventInfo->auxdecor< float >( ( jetName+"_mjj" ).c_str() ) = mjj/ m_units;
        /*if(m_applyResNLOKFactor){
         eventInfo->auxdecor< float >( "weight_resonanceKFactor" ) = 1.-m_ResNLOKFactorHist->GetBinContent( m_ResNLOKFactorHist->FindBin(mjj) );
         }*/
        
        eventInfo->auxdecor< float >( ( jetName+"_pTjj" ).c_str() ) = ( leadJet->p4() + subLeadJet->p4() ).Pt() / m_units;
        eventInfo->auxdecor< float >( ( jetName+"_yBoost" ).c_str() ) = ( leadJet->rapidity() + subLeadJet->rapidity() ) / 2.0;
        eventInfo->auxdecor< float >( ( jetName+"_yStar" ).c_str() ) = ( leadJet->rapidity() - subLeadJet->rapidity() ) / 2.0;
        eventInfo->auxdecor< float >( ( jetName+"_deltaPhi" ).c_str() ) = fabs( TVector2::Phi_mpi_pi( leadJet->phi() - subLeadJet->phi() ) );
        
        eventInfo->auxdecor< float >( ( jetName+"_pTBalance" ).c_str() ) = ( leadJet->pt() - subLeadJet->pt() ) / ( leadJet->pt() + subLeadJet->pt() );
        
        if(inJets->size() >= 3) {
            
            const xAOD::Jet* thirdLeadJet  = inJets->at(2);
            float m23 = ( thirdLeadJet->p4() + subLeadJet->p4() ).M();
            eventInfo->auxdecor< float >( ( jetName+"_m23" ).c_str() ) = m23 / m_units;
            eventInfo->auxdecor< float >( ( jetName+"_m3j" ).c_str() ) = ( inJets->at(0)->p4() + inJets->at(1)->p4() + inJets->at(2)->p4()).M() / m_units;
            
        } // 3 or more jets
    } // 2 or more jets
    
    TLorentzVector MHT = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    TLorentzVector MHTJVT = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    for( auto iJet : *inJets) {
        MHT -= iJet->p4();
        /*if( iJet->pt() < 50e3 &&
         fabs(iJet->getAttribute<xAOD::JetFourMom_t>("JetEMScaleMomentum").eta()) < 2.4 ) {
         if( iJet->getAttribute< float >( "Jvt" ) < 0.64 ) {
         continue;
         }
         }
         MHTJVT -= iJet->p4();*/
    } // end loop over all jet
    eventInfo->auxdecor< float >( ( jetName+"MHT" ).c_str() )    = ( MHT.Pt() / m_units );
    eventInfo->auxdecor< float >( ( jetName+"MHTPhi").c_str() )    = ( MHT.Phi() );
    //eventInfo->auxdecor< float >( "MHTJVT" )    = ( MHTJVT.Pt() / GeV );
    //eventInfo->auxdecor< float >( "MHTJVTPhi" ) = ( MHTJVT.Phi() );
 
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode TLATreeAlgo :: execute ()
{
    const xAOD::EventInfo* eventInfo(nullptr);
    RETURN_CHECK("TLATreeAlgo::execute()", HelperFunctions::retrieve(eventInfo, m_eventInfoContainerName, m_event, m_store, m_verbose) ,"");

    ++m_eventCounter;
    //fill the xsection variables
    if (m_eventCounter == 0 && m_isMC) {
        getLumiWeights(eventInfo);
    } else if (!m_isMC) {
        m_filtEff = m_xs = 1.0;
    }
    
    if(m_isMC)
      m_mcEventWeight = eventInfo->mcEventWeight();
    else
      m_mcEventWeight = 1;
    
    //float weight_pileup = 1.;
    //if( m_isMC && eventInfo->isAvailable< double >( "PileupWeight") )
    //weight_pileup = eventInfo->auxdecor< double >("PileupWeight");
    
    eventInfo->auxdecor< float >("weight_xs") = m_xs * m_filtEff;
    eventInfo->auxdecor< float >("weight") = m_mcEventWeight * m_xs * m_filtEff;
    
    //here, decorate the event with the event-level variables of the jet collections
    
    //get the various jets from the store
    
    //offline jets to start with
    const xAOD::JetContainer* inJets(nullptr);
    RETURN_CHECK("TLATreeAlgo::execute()", HelperFunctions::retrieve(inJets, m_jetContainerName, m_event, m_store, m_verbose) ,"");
    
    std::string jetName = "jet";
    getJetVariables(jetName, inJets, eventInfo);


/*
    //this only if 2 or more jets
    if (inJets->size() >= 2) {
      const xAOD::Jet* leadJet     = inJets->at(0);
      const xAOD::Jet* subLeadJet  = inJets->at(1);
    
      float mjj = ( leadJet->p4() + subLeadJet->p4() ).M() / m_units;
      eventInfo->auxdecor< float >( "mjj" ) = mjj;
    //if(m_applyResNLOKFactor){
    //    eventInfo->auxdecor< float >( "weight_resonanceKFactor" ) = 1.-m_ResNLOKFactorHist->GetBinContent( m_ResNLOKFactorHist->FindBin(mjj) );
    //  }
    
      eventInfo->auxdecor< float >( "pTjj" ) = ( leadJet->p4() + subLeadJet->p4() ).Pt() / m_units;
      eventInfo->auxdecor< float >( "yBoost" ) = ( leadJet->rapidity() + subLeadJet->rapidity() ) / 2.0;
      eventInfo->auxdecor< float >( "deltaPhi" ) = fabs( TVector2::Phi_mpi_pi( leadJet->phi() - subLeadJet->phi() ) );

      eventInfo->auxdecor< float >( "pTBalance" ) = ( leadJet->pt() - subLeadJet->pt() ) / ( leadJet->pt() + subLeadJet->pt() );

      if(inJets->size() >= 3) {
        
        const xAOD::Jet* thirdLeadJet  = inJets->at(2);
        float m23 = ( thirdLeadJet->p4() + subLeadJet->p4() ).M() / m_units;
        eventInfo->auxdecor< float >( "m23" ) = m23 / m_units;
        eventInfo->auxdecor< float >( "m3j" ) = ( inJets->at(0)->p4() + inJets->at(1)->p4() + inJets->at(2)->p4()).M() / m_units;
        
      } // 3 or more jets
    } // 2 or more jets

    TLorentzVector MHT = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    TLorentzVector MHTJVT = TLorentzVector(0.0, 0.0, 0.0, 0.0);
    for( auto iJet : *inJets) {
      MHT -= iJet->p4();
      // if( iJet->pt() < 50e3 &&
      //  fabs(iJet->getAttribute<xAOD::JetFourMom_t>("JetEMScaleMomentum").eta()) < 2.4 ) {
      //  if( iJet->getAttribute< float >( "Jvt" ) < 0.64 ) {
      //    continue;
      //  }
      //}
      //MHTJVT -= iJet->p4();
    } // end loop over all jets
    eventInfo->auxdecor< float >( "MHT"    )    = ( MHT.Pt() / 1000. );
    eventInfo->auxdecor< float >( "MHTPhi" )    = ( MHT.Phi() );
    //eventInfo->auxdecor< float >( "MHTJVT" )    = ( MHTJVT.Pt() / GeV );
    //eventInfo->auxdecor< float >( "MHTJVTPhi" ) = ( MHTJVT.Phi() );
*/
    //call the base class's execute method
    TreeAlgo::execute();

    //for the record, base class execute is below
    
//  // Get EventInfo and the PrimaryVertices
//  const xAOD::EventInfo* eventInfo(nullptr);
//  RETURN_CHECK("TLATreeAlgo::execute()", HelperFunctions::retrieve(eventInfo, m_eventInfoContainerName, m_event, m_store, m_verbose) ,"");
//  const xAOD::VertexContainer* vertices(nullptr);
//  RETURN_CHECK("TLATreeAlgo::execute()", HelperFunctions::retrieve(vertices, "PrimaryVertices", m_event, m_store, m_verbose) ,"");
//  // get the primaryVertex
//  const xAOD::Vertex* primaryVertex = HelperFunctions::getPrimaryVertex( vertices );
//
//  m_helpTree->FillEvent( eventInfo, m_event );
//
//  // Fill trigger information
//  if ( !m_trigDetailStr.empty() )    {
//    m_helpTree->FillTrigger( eventInfo );
//  }
//
//  // Fill jet trigger information - this can be used if with layer/cleaning info we need to turn off some variables?
//  /*if ( !m_trigJetDetailStr.empty() ) {
//    m_helpTree->FillJetTrigger();
//  }*/
//
//  // for the containers the were supplied, fill the appropriate vectors
//  if ( !m_muContainerName.empty() ) {
//    const xAOD::MuonContainer* inMuon(nullptr);
//    RETURN_CHECK("TLATreeAlgo::execute()", HelperFunctions::retrieve(inMuon, m_muContainerName, m_event, m_store, m_verbose) ,"");
//    m_helpTree->FillMuons( inMuon, primaryVertex );
//  }
//
//  if ( !m_elContainerName.empty() ) {
//    const xAOD::ElectronContainer* inElec(nullptr);
//    RETURN_CHECK("TLATreeAlgo::execute()", HelperFunctions::retrieve(inElec, m_elContainerName, m_event, m_store, m_verbose) ,"");
//    m_helpTree->FillElectrons( inElec, primaryVertex );
//  }
//  if ( !m_jetContainerName.empty() ) {
//    const xAOD::JetContainer* inJets(nullptr);
//    RETURN_CHECK("TLATreeAlgo::execute()", HelperFunctions::retrieve(inJets, m_jetContainerName, m_event, m_store, m_verbose) ,"");
//    m_helpTree->FillJets( inJets, HelperFunctions::getPrimaryVertexLocation(vertices), "jet" );
//  }
//  if ( !m_trigJetContainerName.empty() ) {
//    const xAOD::JetContainer* inTrigJets(nullptr);
//    RETURN_CHECK("TLATreeAlgo::execute()", HelperFunctions::retrieve(inTrigJets, m_trigJetContainerName, m_event, m_store, m_verbose) ,"");
//    m_helpTree->FillJets( inTrigJets, HelperFunctions::getPrimaryVertexLocation(vertices), "trigJet" );
//  }
//  if ( !m_truthJetContainerName.empty() ) {
//    const xAOD::JetContainer* inTruthJets(nullptr);
//    RETURN_CHECK("TLATreeAlgo::execute()", HelperFunctions::retrieve(inTruthJets, m_truthJetContainerName, m_event, m_store, m_verbose) ,"");
//        m_helpTree->FillJets( inTruthJets, HelperFunctions::getPrimaryVertexLocation(vertices), "truthJet" );
//  }
//  if ( !m_fatJetContainerName.empty() ) {
//    const xAOD::JetContainer* inFatJets(nullptr);
//    RETURN_CHECK("TLATreeAlgo::execute()", HelperFunctions::retrieve(inFatJets, m_fatJetContainerName, m_event, m_store, m_verbose) ,"");
//    m_helpTree->FillFatJets( inFatJets );
//  }
//  if ( !m_tauContainerName.empty() ) {
//    const xAOD::TauJetContainer* inTaus(nullptr);
//    RETURN_CHECK("HTopMultilepTLATreeAlgo::execute()", HelperFunctions::retrieve(inTaus, m_tauContainerName, m_event, m_store, m_verbose) , "");
//    m_helpTree->FillTaus( inTaus );
//  }
//  if ( !m_METContainerName.empty() ) {
//    const xAOD::MissingETContainer* inMETCont(nullptr);
//    RETURN_CHECK("HTopMultilepTLATreeAlgo::execute()", HelperFunctions::retrieve(inMETCont, m_METContainerName, m_event, m_store, m_debug) , "");
//    m_helpTree->FillMET( inMETCont );
//  }
//  if ( !m_photonContainerName.empty() ) {
//    const xAOD::PhotonContainer* inPhotons(nullptr);
//    RETURN_CHECK("TLATreeAlgo::execute()", HelperFunctions::retrieve(inPhotons, m_photonContainerName, m_event, m_store, m_verbose) ,"");
//    m_helpTree->FillPhotons( inPhotons );
//  }
//
//  // fill the tree
//  m_helpTree->Fill();

  return EL::StatusCode::SUCCESS;

}

//EL::StatusCode TLATreeAlgo :: postExecute () { return EL::StatusCode::SUCCESS; }
//
//EL::StatusCode TLATreeAlgo :: finalize () {
//
//  Info("finalize()", "Deleting tree instances...");
//
//  if ( m_helpTree ) { delete m_helpTree;   m_helpTree = nullptr; }
//
//  return EL::StatusCode::SUCCESS;
//}
//
//EL::StatusCode TLATreeAlgo :: treeFinalize () { return EL::StatusCode::SUCCESS; }

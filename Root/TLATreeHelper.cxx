#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODAnaHelpers/HelperFunctions.h"

#include "TLAAlgos/TLATreeHelper.h"

TLATreeHelper :: TLATreeHelper(xAOD::TEvent * event, TTree* tree, TFile* file, const float units, bool debug, bool DC14, xAOD::TStore* store, bool doJets, bool doTriggerJets, bool doTruthJets) :
HelpTreeBase(event, tree, file, units, debug, DC14, store)
{
    Info("TLATreeHelper", "Creating output TTree");
    m_firstEvent = true;
    m_store = store;
    m_doJets = doJets;
    m_doTriggerJets = doTriggerJets;
    m_doTruthJets = doTruthJets;
}

TLATreeHelper :: ~TLATreeHelper()
{
}
//////////////////// Connect Defined variables to branches here /////////////////////////////
void TLATreeHelper::AddEventUser(const std::string detailStr)
{
    if(m_debug) Info("AddEventUser()", "Adding event-level variables %s", detailStr.c_str());

    // weights
    m_tree->Branch("weight", &m_weight, "weight/F");
    m_tree->Branch("weight_xs", &m_weight_xs, "weight_xs/F");
    m_tree->Branch("distanceFromFront", &m_distanceFromFront, "distanceFromFront/I");
//    m_tree->Branch("weight_prescale", &m_weight_prescale, "weight_prescale/F");
//    m_tree->Branch("weight_resonanceKFactor", &m_weight_resonanceKFactor, "weight_resonanceKFactor/F");
    
    //here we need to find out which jet collections we have in this tree
    
    
    if (m_doJets) AddJetEvent("jet");
    if (m_doTruthJets) AddJetEvent("truthJet");
    if (m_doTriggerJets) AddJetEvent("trigJet");
    
}

void TLATreeHelper::AddJetEvent(const std::string jetName) {
    
    if(m_debug) Info("AddJetEvent()", "Adding event-level jet variables %s", jetName.c_str());
    m_jetEvent[jetName] = new jetEventInfo();
    jetEventInfo* thisJetEvent = m_jetEvent[jetName];
    m_tree->Branch((jetName+"_mjj").c_str(),  &thisJetEvent->m_mjj);
    m_tree->Branch((jetName+"_pTjj").c_str(),  &thisJetEvent->m_pTjj);
    m_tree->Branch((jetName+"_yStar").c_str(),  &thisJetEvent->m_yStar);
    m_tree->Branch((jetName+"_yBoost").c_str(),  &thisJetEvent->m_yBoost);
    m_tree->Branch((jetName+"_deltaPhi").c_str(),  &thisJetEvent->m_deltaPhi);
    m_tree->Branch((jetName+"_pTBalance").c_str(),  &thisJetEvent->m_pTBalance);
    m_tree->Branch((jetName+"_m23").c_str(),  &thisJetEvent->m_m23);
    m_tree->Branch((jetName+"_m3j").c_str(),  &thisJetEvent->m_m3j);
    m_tree->Branch((jetName+"_MHT").c_str(),  &thisJetEvent->m_MHT);
    m_tree->Branch((jetName+"_MHTPhi").c_str(),  &thisJetEvent->m_MHTPhi);
    
    
}

/*void TLATreeHelper::AddJetsUser(const std::string detailStr, const std::string jetName)
{
//    m_tree->Branch((jetName+"_constitScaleEta").c_str(), &m_jet_constitScaleEta);
//    m_tree->Branch((jetName+"_emScaleE").c_str(), &m_jet_emScaleE);
//    m_tree->Branch((jetName+"_emScaleEta").c_str(), &m_jet_emScaleEta);
//    m_tree->Branch((jetName+"_emScalePhi").c_str(), &m_jet_emScalePhi);
    m_tree->Branch((jetName+"_minDeltaR").c_str(), &m_jet_minDeltaR);
    m_tree->Branch((jetName+"_numberCloseJets").c_str(), &m_jet_numberCloseJets);
}*/

//////////////////// Clear any defined vectors here ////////////////////////////
void TLATreeHelper::ClearEventUser() {
    
//    m_weight_corr = -999;
    m_weight    = -999;
    m_weight_xs = -999;
    
    m_distanceFromFront = -999;
//    m_weight_prescale = -999;
//    m_weight_resonanceKFactor = -999;
    
    //clear jet event variables
    
    if (m_doJets) ClearJetEvent("jet");
    if (m_doTruthJets) ClearJetEvent("truthJet");
    if (m_doTriggerJets) ClearJetEvent("trigJet");

}

void TLATreeHelper::ClearJetEvent(const std::string jetName) {
    
    jetEventInfo* thisJetEvent = m_jetEvent[jetName];
    thisJetEvent->m_mjj = -999;
    thisJetEvent->m_pTjj = -999;
    thisJetEvent->m_yStar = -999;
    thisJetEvent->m_yBoost = -999;
    thisJetEvent->m_deltaPhi = -999;
    thisJetEvent->m_pTBalance = -999;
    thisJetEvent->m_m23 = -999;
    thisJetEvent->m_m3j = -999;
    thisJetEvent->m_MHT = -999;
    thisJetEvent->m_MHTPhi= -999;


}

/*void TLATreeHelper::ClearJetsUser(const std::string jetName ) {
//    m_jet_constitScaleEta.clear();
//    m_jet_emScaleE.clear();
//    m_jet_emScaleEta.clear();
//    m_jet_emScalePhi.clear();
    m_jet_minDeltaR.clear();
    m_jet_numberCloseJets.clear();
}*/

/////////////////// Assign values to defined event variables here ////////////////////////
void TLATreeHelper::FillEventUser( const xAOD::EventInfo* eventInfo ) {

//%%%%Left for later: variables that can be derived from what we have already
    
    if( eventInfo->isAvailable< float >( "weight" ) )
      m_weight = eventInfo->auxdecor< float >( "weight" );
    if( eventInfo->isAvailable< float >( "weight_xs" ) )
      m_weight_xs = eventInfo->auxdecor< float >( "weight_xs" );
    if( eventInfo->isAvailable< int >( "distanceFromFront" ) )
      m_distanceFromFront = eventInfo->auxdecor< int >( "distanceFromFront" );
//    if( eventInfo->isAvailable< float >( "weight_prescale" ) )
//    m_weight_prescale = eventInfo->auxdecor< float >( "weight_prescale" );
//    if( eventInfo->isAvailable< float >( "weight_resonanceKFactor" ) )
//    m_weight_resonanceKFactor = eventInfo->auxdecor< float >( "weight_resonanceKFactor" );
    

    if (m_doJets) FillJetEvent(eventInfo, "jet");
    if (m_doTruthJets) FillJetEvent(eventInfo, "truthJet");
    if (m_doTriggerJets) FillJetEvent(eventInfo, "trigJet");
    
}

void TLATreeHelper::FillJetEvent( const xAOD::EventInfo* eventInfo, const std::string jetName) {
    
    if(m_debug) Info("FillJetEvent()", "Filling event-level jet variables for branch '%s'", jetName.c_str());
    jetEventInfo* thisJetEvent = m_jetEvent[jetName];
    
    if (eventInfo->isAvailable< float > ((jetName+"_mjj").c_str()) )
        thisJetEvent->m_mjj = eventInfo->auxdecor< float > ((jetName+"_mjj").c_str());
    else thisJetEvent->m_mjj = -999;

    if (eventInfo->isAvailable< float > ((jetName+"_pTjj").c_str()) )
        thisJetEvent->m_pTjj = eventInfo->auxdecor< float > ((jetName+"_pTjj").c_str());
    else thisJetEvent->m_pTjj = -999;

    if (eventInfo->isAvailable< float > ((jetName+"_yStar").c_str()) )
        thisJetEvent->m_yStar = eventInfo->auxdecor< float > ((jetName+"_yStar").c_str());
    else thisJetEvent->m_yStar = -999;

    if (eventInfo->isAvailable< float > ((jetName+"_yBoost").c_str()) )
        thisJetEvent->m_yBoost = eventInfo->auxdecor< float > ((jetName+"_yBoost").c_str());
    else thisJetEvent->m_yBoost = -999;

    if (eventInfo->isAvailable< float > ((jetName+"_deltaPhi").c_str()) )
        thisJetEvent->m_deltaPhi = eventInfo->auxdecor< float > ((jetName+"_deltaPhi").c_str());
    else thisJetEvent->m_deltaPhi = -999;

    if (eventInfo->isAvailable< float > ((jetName+"_pTBalance").c_str()) )
        thisJetEvent->m_pTBalance = eventInfo->auxdecor< float > ((jetName+"_pTBalance").c_str());
    else thisJetEvent->m_pTBalance = -999;

    if (eventInfo->isAvailable< float > ((jetName+"_m23").c_str()) )
        thisJetEvent->m_m23 = eventInfo->auxdecor< float > ((jetName+"_m23").c_str());
    else thisJetEvent->m_m23 = -999;

    if (eventInfo->isAvailable< float > ((jetName+"_m3j").c_str()) )
        thisJetEvent->m_m3j = eventInfo->auxdecor< float > ((jetName+"_m3j").c_str());
    else thisJetEvent->m_m3j = -999;

    if (eventInfo->isAvailable< float > ((jetName+"_MHT").c_str()) )
        thisJetEvent->m_MHT = eventInfo->auxdecor< float > ((jetName+"_MHT").c_str());
    else thisJetEvent->m_MHT = -999;
    
    if (eventInfo->isAvailable< float > ((jetName+"_MHTPhi").c_str()) )
        thisJetEvent->m_MHTPhi = eventInfo->auxdecor< float > ((jetName+"_MHTPhi").c_str());
    else thisJetEvent->m_MHTPhi = -999;

}


//%%%%Left for later: jet variables

/////////////////// Assign values to defined jet variables here //////////////////
//void TLATreeHelper::FillJetsUser( const xAOD::Jet* jet, const std::string jetName ) {
//    if( jet->isAvailable< float >( (jetName+"_constitScaleEta").c_str() ) ) {
//        m_jet_constitScaleEta.push_back( jet->auxdata< float >((jetName+"_constitScaleEta").c_str()));
//    } else {
//        m_jet_constitScaleEta.push_back( -999 );
//    }
//    if( jet->isAvailable< float >( (jetName+"_emScaleE").c_str() ) ) {
//        m_jet_emScaleE.push_back( jet->auxdata< float >((jetName+"_emScaleE").c_str() ) );
//    } else {
//        m_jet_emScaleE.push_back( -999 );
//    }
//    if( jet->isAvailable< float >( (jetName+"_emScaleEta").c_str() ) ) {
//        m_jet_emScaleEta.push_back( jet->auxdata< float >((jetName+"_emScaleEta").c_str() ) );
//    } else {
//        m_jet_emScaleEta.push_back( -999 );
//    }
//    if( jet->isAvailable< float >( (jetName+"_emScalePhi").c_str() ) ) {
//        m_jet_emScalePhi.push_back( jet->auxdata< float >((jetName+"_emScalePhi").c_str() ) );
//    } else {
//        m_jet_emScalePhi.push_back( -999 );
//    }
//    
//    if( jet->isAvailable< float >( (jetName+"_minDeltaR").c_str() ) ) {
//        m_jet_minDeltaR.push_back( jet->auxdata< float >((jetName+"_minDeltaR").c_str() ) );
//    } else {
//        m_jet_minDeltaR.push_back( -999 );
//    }
//    
//    if( jet->isAvailable< int >( (jetName+"_numberCloseJets").c_str() ) ) {
//        m_jet_numberCloseJets.push_back( jet->auxdata< int >((jetName+"_numberCloseJets").c_str() ) );
//    } else {
//        m_jet_numberCloseJets.push_back( -999 );
//    }
//    
//}


//%%%%Left for later: b-tagging
/////////////////// Add variables for the b-tagged di-jet analysis //////////////////
//void TLATreeHelper::AddBtagHighEff() {
//    
//    m_tree->Branch("mbb_fix_7777",  &m_mbb_fix_7777,  "mbb_fix_7777/F");
//    m_tree->Branch("mbb_fix_8585",  &m_mbb_fix_8585,  "mbb_fix_8585/F");
//    m_tree->Branch("mbb_flt_7777",  &m_mbb_flt_7777,  "mbb_flt_7777/F");
//    m_tree->Branch("mbb_flt_8585",  &m_mbb_flt_8585,  "mbb_flt_8585/F");
//    
//    m_tree->Branch("mbj_fix_7777",  &m_mbj_fix_7777,  "mbj_fix_7777/F");
//    m_tree->Branch("mbj_fix_8585",  &m_mbj_fix_8585,  "mbj_fix_8585/F");
//    m_tree->Branch("mbj_flt_7777",  &m_mbj_flt_7777,  "mbj_flt_7777/F");
//    m_tree->Branch("mbj_flt_8585",  &m_mbj_flt_8585,  "mbj_flt_8585/F");
//    
//    m_tree->Branch("weight_btag_fix_77", &m_weight_btag_fix_77);
//    m_tree->Branch("weight_btag_fix_85", &m_weight_btag_fix_85);
//    m_tree->Branch("weight_btag_flt_77", &m_weight_btag_flt_77);
//    m_tree->Branch("weight_btag_flt_85", &m_weight_btag_flt_85);
//    
//    m_tree->Branch("systSF_btag_names",      &m_systSF_btag_names    );
//}

//void TLATreeHelper::AddBtagLowEff() {
//    
//    m_tree->Branch("mbb_fix_6060",  &m_mbb_fix_6060,  "mbb_fix_6060/F");
//    m_tree->Branch("mbb_fix_7070",  &m_mbb_fix_7070,  "mbb_fix_7070/F");
//    m_tree->Branch("mbb_flt_6060",  &m_mbb_flt_6060,  "mbb_flt_6060/F");
//    m_tree->Branch("mbb_flt_7070",  &m_mbb_flt_7070,  "mbb_flt_7070/F");
//    
//    m_tree->Branch("weight_btag_fix_60", &m_weight_btag_fix_60, "weight_btag_fix_60/F");
//    m_tree->Branch("weight_btag_fix_70", &m_weight_btag_fix_70, "weight_btag_fix_70/F");
//    m_tree->Branch("weight_btag_flt_60", &m_weight_btag_flt_60, "weight_btag_flt_60/F");
//    m_tree->Branch("weight_btag_flt_70", &m_weight_btag_flt_70, "weight_btag_flt_70/F");
//}

//void TLATreeHelper::FillBtagHighEff( const xAOD::EventInfo* eventInfo ) {
//    
//    this->ClearBtagHighEff();
//    
//    if( eventInfo->isAvailable< float >( "mbb_FixedCutBEff_77" ) )
//    m_mbb_fix_7777 = eventInfo->auxdecor< float >( "mbb_FixedCutBEff_77" );
//    if( eventInfo->isAvailable< float >( "mbb_FixedCutBEff_85" ) )
//    m_mbb_fix_8585 = eventInfo->auxdecor< float >( "mbb_FixedCutBEff_85" );
//    if( eventInfo->isAvailable< float >( "mbb_FlatBEff_77" ) )
//    m_mbb_flt_7777 = eventInfo->auxdecor< float >( "mbb_FlatBEff_77" );
//    if( eventInfo->isAvailable< float >( "mbb_FlatBEff_85" ) )
//    m_mbb_flt_8585 = eventInfo->auxdecor< float >( "mbb_FlatBEff_85" );
//    
//    if( eventInfo->isAvailable< float >( "mbj_FixedCutBEff_77" ) )
//    m_mbj_fix_7777 = eventInfo->auxdecor< float >( "mbj_FixedCutBEff_77" );
//    if( eventInfo->isAvailable< float >( "mbj_FixedCutBEff_85" ) )
//    m_mbj_fix_8585 = eventInfo->auxdecor< float >( "mbj_FixedCutBEff_85" );
//    if( eventInfo->isAvailable< float >( "mbj_FlatBEff_77" ) )
//    m_mbj_flt_7777 = eventInfo->auxdecor< float >( "mbj_FlatBEff_77" );
//    if( eventInfo->isAvailable< float >( "mbj_FlatBEff_85" ) )
//    m_mbj_flt_8585 = eventInfo->auxdecor< float >( "mbj_FlatBEff_85" );
//    
//    if( eventInfo->isAvailable< std::vector< float >  >( "weight_BTag_FixedCutBEff_77" ) )
//    m_weight_btag_fix_77 = eventInfo->auxdecor< std::vector< float >  >( "weight_BTag_FixedCutBEff_77" );
//    if( eventInfo->isAvailable< std::vector< float >  >( "weight_BTag_FixedCutBEff_85" ) )
//    m_weight_btag_fix_85 = eventInfo->auxdecor< std::vector< float >  >( "weight_BTag_FixedCutBEff_85" );
//    if( eventInfo->isAvailable< std::vector< float >  >( "weight_BTag_FlatBEff_77" ) )
//    m_weight_btag_flt_77 = eventInfo->auxdecor< std::vector< float >  >( "weight_BTag_FlatBEff_77" );
//    if( eventInfo->isAvailable< std::vector< float >  >( "weight_BTag_FlatBEff_85" ) )
//    m_weight_btag_flt_85 = eventInfo->auxdecor< std::vector< float >  >( "weight_BTag_FlatBEff_85" );
//    
//    
//    // get vector of string giving the syst names of the upstream algo from TStore
//    // (rememeber: 1st element is a blank string: nominal case!)
//    if( m_firstEvent){
//        std::vector< std::string >* systNames(nullptr);
//        if( m_store->contains<std::vector<std::string> >( "BJetEfficiency_Algo_FixedCutBEff_85" )){
//            HelperFunctions::retrieve(systNames, "BJetEfficiency_Algo_FixedCutBEff_85", 0, m_store, false); //.isSuccess()
//            std::cout << " Using " <<  systNames->size() << " systematics from BJetEfficiency_Algo_FixedCutBEff_85 " << std::endl;
//        }else if( m_store->contains<std::vector<std::string> >( "BJetEfficiency_Algo_FixedCutBEff_77" )){
//            HelperFunctions::retrieve(systNames, "BJetEfficiency_Algo_FixedCutBEff_77", 0, m_store, false); //.isSuccess()
//            std::cout << " Using " << systNames->size() << " systematics from BJetEfficiency_Algo_FixedCutBEff_77 " << std::endl;
//        } else {
//            std::cout << " No b-tagging systematics found !!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
//        }
//        
//        m_systSF_btag_names = *systNames;
//        if( m_systSF_btag_names.size() > 0 && m_systSF_btag_names.at(0).size() == 0)
//        m_systSF_btag_names.at(0) = "Nominal";
//        
//        m_firstEvent = false;
//    }
//    
//    
//    
//}
//
//void TLATreeHelper::FillBtagLowEff( const xAOD::EventInfo* eventInfo ) {
//    
//    this->ClearBtagLowEff();
//    
//    if( eventInfo->isAvailable< float >( "mbb_FixedCutBEff_60" ) )
//    m_mbb_fix_6060 = eventInfo->auxdecor< float >( "mbb_FixedCutBEff_60" );
//    if( eventInfo->isAvailable< float >( "mbb_FixedCutBEff_70" ) )
//    m_mbb_fix_7070 = eventInfo->auxdecor< float >( "mbb_FixedCutBEff_70" );
//    if( eventInfo->isAvailable< float >( "mbb_FlatBEff_60" ) )
//    m_mbb_flt_6060 = eventInfo->auxdecor< float >( "mbb_FlatBEff_60" );
//    if( eventInfo->isAvailable< float >( "mbb_FlatBEff_70" ) )
//    m_mbb_flt_7070 = eventInfo->auxdecor< float >( "mbb_FlatBEff_70" );
//    
//    if( eventInfo->isAvailable< float >( "weight_BTag_FixedCutBEff_60" ) )
//    m_weight_btag_fix_60 = eventInfo->auxdecor< float >( "weight_BTag_FixedCutBEff_60" );
//    if( eventInfo->isAvailable< float >( "weight_BTag_FixedCutBEff_70" ) )
//    m_weight_btag_fix_70 = eventInfo->auxdecor< float >( "weight_BTag_FixedCutBEff_70" );
//    if( eventInfo->isAvailable< float >( "weight_BTag_FlatBEff_60" ) )
//    m_weight_btag_flt_60 = eventInfo->auxdecor< float >( "weight_BTag_FlatBEff_60" );
//    if( eventInfo->isAvailable< float >( "weight_BTag_FlatBEff_70" ) )
//    m_weight_btag_flt_70 = eventInfo->auxdecor< float >( "weight_BTag_FlatBEff_70" );
//}
//
//void TLATreeHelper::ClearBtagHighEff() {
//    
//    m_mbb_fix_7777 = -999;
//    m_mbb_fix_8585 = -999;
//    m_mbb_flt_7777 = -999;
//    m_mbb_flt_8585 = -999;
//    
//    m_mbj_fix_7777 = -999;
//    m_mbj_fix_8585 = -999;
//    m_mbj_flt_7777 = -999;
//    m_mbj_flt_8585 = -999;
//    
//    m_weight_btag_fix_77.clear();
//    m_weight_btag_fix_85.clear();
//    m_weight_btag_flt_77.clear();
//    m_weight_btag_flt_85.clear();
//    
//    m_systSF_btag_names.clear();
//}
//
//void TLATreeHelper::ClearBtagLowEff() {
//    
//    m_mbb_fix_6060 = -999;
//    m_mbb_fix_7070 = -999;
//    m_mbb_flt_6060 = -999;
//    m_mbb_flt_7070 = -999;
//    
//    m_weight_btag_fix_60 = -999;
//    m_weight_btag_fix_70 = -999;
//    m_weight_btag_flt_60 = -999;
//    m_weight_btag_flt_70 = -999;
//}

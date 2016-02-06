#ifndef AnalysisExample_TLATreeHelper_H
#define AnalysisExample_TLATreeHelper_H

#include "xAODAnaHelpers/HelpTreeBase.h"
#include "TTree.h"

class TLATreeHelper : public HelpTreeBase
{
    
private:
    
    //struct to save event-level variables that need to belong to each jet collection
    //
    // jetEventInfo
    //
    struct jetEventInfo{
        
        float m_mjj;
        float m_pTjj;
        float m_yStar;
        float m_yBoost;
        float m_deltaPhi;
        float m_pTBalance;
        float m_m23;
        float m_m3j;
        float m_MHT;
        float m_MHTPhi;
        
        jetEventInfo(){ }
    };
    std::map<std::string, jetEventInfo*> m_jetEvent;

    bool m_firstEvent;
    
    bool m_doJets;
    bool m_doTriggerJets;
    bool m_doTruthJets;
    
    float m_weight;
//    float m_weight_corr;
    float m_weight_xs;
    int m_distanceFromFront;
//    float m_weight_prescale;
//    float m_weight_resonanceKFactor;
    
//
//    std::vector<float> m_jet_constitScaleEta;
//    std::vector<float> m_jet_emScaleE;
//    std::vector<float> m_jet_emScaleEta;
//    std::vector<float> m_jet_emScalePhi;
//    std::vector<float> m_jet_minDeltaR;
//    std::vector<int>   m_jet_numberCloseJets;
    
public:
    
    TLATreeHelper(xAOD::TEvent * event, TTree* tree, TFile* file, const float units, bool debug = false, bool DC14 = false, xAOD::TStore* store = nullptr, bool doJets = true, bool doTriggerJets = false, bool doTruthJets = false);
    ~TLATreeHelper();
    
    void AddEventUser( const std::string detailStr = "" );
    //void AddJetsUser( const std::string detailStr = "" , const std::string jetName = "jet");
    void FillEventUser( const xAOD::EventInfo* eventInfo );
    //void FillJetsUser( const xAOD::Jet* jet, const std::string jetName = "jet" );
    void ClearEventUser();
    
    void AddJetEvent( const std::string jetName);
    void FillJetEvent( const xAOD::EventInfo* eventInfo, const std::string jetName);
    void ClearJetEvent(const std::string jetName);
    

    //void ClearJetsUser(const std::string jetName = "jet");
    
//    void AddBtagHighEff();
//    void FillBtagHighEff( const xAOD::EventInfo* eventInfo );
//    void ClearBtagHighEff();
//    
//    void AddBtagLowEff();
//    void FillBtagLowEff( const xAOD::EventInfo* eventInfo );
//    void ClearBtagLowEff();
    
};
#endif

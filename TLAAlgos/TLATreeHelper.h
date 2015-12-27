#ifndef AnalysisExample_TLATreeHelper_H
#define AnalysisExample_TLATreeHelper_H

#include "xAODAnaHelpers/HelpTreeBase.h"
#include "TTree.h"

class TLATreeHelper : public HelpTreeBase
{
    
private:
    bool m_firstEvent;
    
    float m_yStar;
    float m_yBoost;
    float m_mjj;
    float m_pTjj;
    float m_m3j;
    float m_deltaPhi;
    
    float m_mbb_fix_6060;
    float m_mbb_fix_7070;
    float m_mbb_flt_6060;
    float m_mbb_flt_7070;
    
    float m_mbb_fix_7777;
    float m_mbb_fix_8585;
    float m_mbb_flt_7777;
    float m_mbb_flt_8585;
    
    float m_mbj_fix_7777;
    float m_mbj_fix_8585;
    float m_mbj_flt_7777;
    float m_mbj_flt_8585;
    
    float m_weight;
    float m_weight_corr;
    float m_weight_xs;
    float m_weight_prescale;
    float m_weight_resonanceKFactor;
    
    float m_weight_btag_fix_60;
    float m_weight_btag_fix_70;
    float m_weight_btag_flt_60;
    float m_weight_btag_flt_70;
    std::vector< float > m_weight_btag_fix_77;
    std::vector< float > m_weight_btag_fix_85;
    std::vector< float > m_weight_btag_flt_77;
    std::vector< float > m_weight_btag_flt_85;
    
    std::vector< std::string > m_systSF_btag_names;
    
    float m_pTbalance;
    float m_MHT;
    float m_MHTPhi;
    float m_MHTJVT;
    float m_MHTJVTPhi;
    
    float m_Insitu_Segs_response_E;
    float m_Insitu_Segs_response_pT;
    int m_punch_type_segs;
    
    std::vector<float> m_jet_constitScaleEta;
    std::vector<float> m_jet_emScaleE;
    std::vector<float> m_jet_emScaleEta;
    std::vector<float> m_jet_emScalePhi;
    std::vector<float> m_jet_minDeltaR;
    std::vector<int>   m_jet_numberCloseJets;
    
public:
    
    TLATreeHelper(xAOD::TEvent * event, TTree* tree, TFile* file, xAOD::TStore* store = nullptr);
    ~TLATreeHelper();
    
    void AddEventUser( const std::string detailStr = "" );
    void AddJetsUser( const std::string detailStr = "" , const std::string jetName = "jet");
    void FillEventUser( const xAOD::EventInfo* eventInfo );
    void FillJetsUser( const xAOD::Jet* jet, const std::string jetName = "jet" );
    void ClearEventUser();
    void ClearJetsUser(const std::string jetName = "jet");
    
    void AddBtagHighEff();
    void FillBtagHighEff( const xAOD::EventInfo* eventInfo );
    void ClearBtagHighEff();
    
    void AddBtagLowEff();
    void FillBtagLowEff( const xAOD::EventInfo* eventInfo );
    void ClearBtagLowEff();
    
};
#endif

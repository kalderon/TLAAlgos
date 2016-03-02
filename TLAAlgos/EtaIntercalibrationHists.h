#ifndef xAODAnaHelpers_EtaIntercalibrationHists_H
#define xAODAnaHelpers_EtaIntercalibrationHists_H

#include "xAODAnaHelpers/HistogramManager.h"
#include "xAODAnaHelpers/HelperClasses.h"
#include <xAODJet/JetContainer.h>

class EtaIntercalibrationHists : public HistogramManager
{
  public:


    EtaIntercalibrationHists(std::string name, std::string detailStr);
    virtual ~EtaIntercalibrationHists() ;

    bool m_debug;
    StatusCode initialize();
    StatusCode execute( const xAOD::JetContainer* jets, float eventWeight, int pvLoc = -1);
    StatusCode execute( const xAOD::Jet* jet, float eventWeight, int pvLoc = -1 );
    using HistogramManager::book; // make other overloaded version of book() to show up in subclass
    using HistogramManager::execute; // overload

  protected:

    // holds bools that control which histograms are filled
    HelperClasses::JetInfoSwitch* m_infoSwitch;

  private:
    
    std::map<std::string,TH1F*> m_hists1D;
//    std::map<std::string,TH2F*> m_hists2D;
//    std::map<std::string,TH3F*> m_hists3D;

    //basic
    TH1F* m_jetPt;                  //!
    TH1F* m_jetEta;                 //!
    TH1F* m_jetPhi;                 //!
    TH1F* m_jetM;                   //!
    TH1F* m_jetE;                   //!
    TH1F* m_jetRapidity;            //!

//    // kinematic
//    TH1F* m_jetPx;                  //!
//    TH1F* m_jetPy;                  //!
//    TH1F* m_jetPz;                  //!
//
//    //NLeadingJets
//    std::vector< TH1F* > m_NjetsPt;       //!
//    std::vector< TH1F* > m_NjetsEta;      //!
//    std::vector< TH1F* > m_NjetsPhi;      //!
//    std::vector< TH1F* > m_NjetsM;        //!
//    std::vector< TH1F* > m_NjetsE;        //!
//    std::vector< TH1F* > m_NjetsRapidity; //!
//
//    // clean
//    TH1F* m_jetTime;                //!
//    TH1F* m_LArQuality;             //!
//    TH1F* m_hecq;                   //!
//    TH1F* m_negE;                   //!
//    TH1F* m_avLArQF;                //!
//    TH1F* m_bchCorrCell;            //!
//    TH1F* m_N90Const;               //!
//
//    //layer
//    TH1F* m_PreSamplerB;
//    TH1F* m_EMB1;
//    TH1F* m_EMB2;
//    TH1F* m_EMB3;
//    TH1F* m_PreSamplerE;            //!
//    TH1F* m_EME1;                   //!
//    TH1F* m_EME2;                   //!
//    TH1F* m_EME3;                   //!
//    TH1F* m_HEC0;                   //!
//    TH1F* m_HEC1;                   //!
//    TH1F* m_HEC2;                   //!
//    TH1F* m_HEC3;                   //!
//    TH1F* m_TileBar0;               //!
//    TH1F* m_TileBar1;               //!
//    TH1F* m_TileBar2;               //!
//    TH1F* m_TileGap1;               //!
//    TH1F* m_TileGap2;               //!
//    TH1F* m_TileGap3;               //!
//    TH1F* m_TileExt0;               //!
//    TH1F* m_TileExt1;               //!
//    TH1F* m_TileExt2;               //!
//    TH1F* m_FCAL0;                  //!
//    TH1F* m_FCAL1;                  //!
//    TH1F* m_FCAL2;                  //!
//
//    // area
//    TH1F* m_actArea;                //!
//
//
//    TH1F* m_chf;                    //!
//
//    //energy
//    TH1F* m_HECf;                   //!
//    TH1F* m_EMf;                    //!
//    TH1F* m_centroidR;              //!
//    TH1F* m_fracSampMax;            //!
//    TH1F* m_fracSampMaxIdx;         //!
//    TH1F* m_lowEtFrac;              //!
//
//    // resolution
//    TH1F* m_jetGhostTruthPt;        //!
//    TH2F* m_jetPt_vs_resolution;    //!
//    TH2F* m_jetGhostTruthPt_vs_resolution; //!
//
//    // truth jets
//    TH1F* m_truthLabelID;          //!
//    TH1F* m_truthCount;            //!
//    TH1F* m_truthPt;               //!
//    TH1F* m_truthDr_B;             //!
//    TH1F* m_truthDr_C;             //!
//    TH1F* m_truthDr_T;             //!
//
//    // Detailed truth jet plots
//    TH1F* m_truthCount_BhadFinal;  //!
//    TH1F* m_truthCount_BhadInit ;  //!
//    TH1F* m_truthCount_BQFinal  ;  //!
//    TH1F* m_truthPt_BhadFinal;  //!
//    TH1F* m_truthPt_BhadInit ;  //!
//    TH1F* m_truthPt_BQFinal  ;  //!
//
//    TH1F* m_truthCount_ChadFinal;  //!
//    TH1F* m_truthCount_ChadInit ;  //!
//    TH1F* m_truthCount_CQFinal  ;  //!
//    TH1F* m_truthPt_ChadFinal;  //!
//    TH1F* m_truthPt_ChadInit ;  //!
//    TH1F* m_truthPt_CQFinal  ;  //!
//
//
//    TH1F* m_truthCount_TausFinal; //!
//    TH1F* m_truthPt_TausFinal   ; //!
//
//    // Flavor Tag
//    TH1F* m_MV1   ; //!
//    TH1F* m_MV2c00   ; //!
//    TH1F* m_MV2c10   ; //!
//    TH1F* m_MV2c20   ; //!
//    TH1F* m_COMB   ; //!
//    TH1F* m_SV0             ; //!
//    TH1F* m_SV1             ; //!
//    TH1F* m_IP2D            ; //!
//    TH1F* m_IP3D            ; //!
//    TH1F* m_JetFitter       ; //!
//    TH1F* m_JetFitterCombNN ; //!
//
//    TH1F* m_jf_nVTX           ; //!
//    TH1F* m_jf_nSingleTracks  ; //!
//    TH1F* m_jf_nTracksAtVtx   ; //!
//    TH1F* m_jf_mass           ; //!
//    TH1F* m_jf_energyFraction ; //!
//    TH1F* m_jf_significance3d ; //!
//    TH1F* m_jf_deltaeta       ; //!
//    TH1F* m_jf_deltaphi       ; //!
//    TH1F* m_jf_N2Tpar         ; //!
//    TH1F* m_jf_pb             ; //!
//    TH1F* m_jf_pc             ; //!
//    TH1F* m_jf_pu             ; //!
//
//    TH1F* m_sv0_NGTinSvx ; //!
//    TH1F* m_sv0_N2Tpair  ; //!
//    TH1F* m_sv0_massvx   ; //!
//    TH1F* m_sv0_efracsvx ; //!
//    TH1F* m_sv0_normdist ; //!
//
//    TH1F* m_sv1_NGTinSvx ; //!
//    TH1F* m_sv1_N2Tpair  ; //!
//    TH1F* m_sv1_massvx   ; //!
//    TH1F* m_sv1_efracsvx ; //!
//    TH1F* m_sv1_normdist ; //!
//
//    TH1F* m_nIP2DTracks              ; //!
//    TH1F* m_IP2D_gradeOfTracks       ; //!
//    TH1F* m_IP2D_flagFromV0ofTracks  ; //!
//    TH1F* m_IP2D_valD0wrtPVofTracks  ; //!
//    TH1F* m_IP2D_sigD0wrtPVofTracks  ; //!
//    TH1F* m_IP2D_sigD0wrtPVofTracks_l; //!
//    TH1F* m_IP2D_weightBofTracks     ; //!
//    TH1F* m_IP2D_weightCofTracks     ; //!
//    TH1F* m_IP2D_weightUofTracks     ; //!
//
//    TH1F* m_nIP3DTracks              ; //!
//    TH1F* m_IP3D_gradeOfTracks       ; //!
//    TH1F* m_IP3D_flagFromV0ofTracks  ; //!
//    TH1F* m_IP3D_valD0wrtPVofTracks  ; //!
//    TH1F* m_IP3D_sigD0wrtPVofTracks  ; //!
//    TH1F* m_IP3D_sigD0wrtPVofTracks_l; //!
//    TH1F* m_IP3D_valZ0wrtPVofTracks  ; //!
//    TH1F* m_IP3D_sigZ0wrtPVofTracks  ; //!
//    TH1F* m_IP3D_sigZ0wrtPVofTracks_l; //!
//    TH1F* m_IP3D_weightBofTracks     ; //!
//    TH1F* m_IP3D_weightCofTracks     ; //!
//    TH1F* m_IP3D_weightUofTracks     ; //!
//
//    // substructure
//    TH1F* m_tau1; //!
//    TH1F* m_tau2; //!
//    TH1F* m_tau3; //!
//    TH1F* m_tau21; //!
//    TH1F* m_tau32; //!
//    TH1F* m_tau1_wta; //!
//    TH1F* m_tau2_wta; //!
//    TH1F* m_tau3_wta; //!
//    TH1F* m_tau21_wta; //!
//    TH1F* m_tau32_wta; //!
//    TH1F* m_numConstituents; //!
};

#endif

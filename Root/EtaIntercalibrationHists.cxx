#include <TLAAlgos/EtaIntercalibrationHists.h>
#include <sstream>

#include "xAODAnaHelpers/tools/ReturnCheck.h"

using std::vector;

EtaIntercalibrationHists :: EtaIntercalibrationHists (std::string name, std::string detailStr) :
  HistogramManager(name, detailStr),
  m_infoSwitch(new HelperClasses::JetInfoSwitch(m_detailStr))
{
  m_debug = false;
}

EtaIntercalibrationHists :: ~EtaIntercalibrationHists () {
  if(m_infoSwitch) delete m_infoSwitch;
}

StatusCode EtaIntercalibrationHists::initialize() {

  if(m_debug) Info("EtaIntercalibrationHists::initialize()", m_name.c_str());

  // All the plots are always made
  m_jetPt          = book(m_name, "jetPt",  "jet p_{T} [GeV]", 120, 0, 3000.);
  m_jetEta         = book(m_name, "jetEta", "jet #eta",         80, -4, 4);
  m_jetPhi         = book(m_name, "jetPhi", "jet Phi",120, -TMath::Pi(), TMath::Pi() );
  m_jetM           = book(m_name, "jetMass", "jet Mass [GeV]",120, 0, 400);
  m_jetE           = book(m_name, "jetEnergy", "jet Energy [GeV]",120, 0, 4000.);
  m_jetRapidity    = book(m_name, "jetRapidity", "jet Rapidity",120, -10, 10);

    
    
  //individual eta intercalibration histograms
    std::cout <<"m_name" << m_name << std::endl;
    std::string prefix = m_name+"_";
    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"pTavg", book(m_name, prefix+"pTavg",";p_{T}^{avg};",2500,0,2500)) );

//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"jet2_pT", createHist1D(prefix+"jet2_pT",";Subleading jet p_{T} [GeV];",2500,0,2500)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"jet3_pT", createHist1D(prefix+"jet3_pT",";Third jet p_{T} [GeV];",2500,0,2500)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"jet1_eta", createHist1D(prefix+"jet1_eta",";Leading jet#eta;",100,-5,5)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"jet2_eta", createHist1D(prefix+"jet2_eta",";Subleading jet#eta;",100,-5,5)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"jet3_eta", createHist1D(prefix+"jet3_eta",";Third jet#eta;",100,-5,5)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"jet1_phi", createHist1D(prefix+"jet1_phi",";Leading jet #phi;",100,-TMath::Pi(),TMath::Pi())) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"jet2_phi", createHist1D(prefix+"jet2_phi",";Subleading jet #phi;",100,-TMath::Pi(),TMath::Pi())) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"jet3_phi", createHist1D(prefix+"jet3_phi",";Third jet #phi;",100,-TMath::Pi(),TMath::Pi())) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"pTavg", createHist1D(prefix+"pTavg",";p_{T}^{avg};",2500,0,2500)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"deltaPhi", createHist1D(prefix+"deltaPhi",";#Delta#phi(j1,j2);",100,0,TMath::Pi())) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"eta_right", createHist1D(prefix+"eta_right",";#eta_{right};",100,-5,5)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"eta_left", createHist1D(prefix+"eta_left",";#eta_{left};",100,-5,5)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"pT_right", createHist1D(prefix+"pT_right",";p_{T}^{right};",100,0,500)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"pT_left", createHist1D(prefix+"pT_left",";p_{T}^{left};",100,0,500)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"eta_ref", createHist1D(prefix+"eta_ref",";#eta_{ref};",100,-5,5)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"eta_probe", createHist1D(prefix+"eta_probe",";#eta_{probe};",100,-5,5)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"pT_ref", createHist1D(prefix+"pT_ref",";p_{T}^{ref};",100,0,500)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"pT_probe", createHist1D(prefix+"pT_probe",";p_{T}^{probe};",100,0,500)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"Asm_SM", createHist1D(prefix+"Asm_SM",";A_{SM};",100,-5,5)) );
//    m_hists1D.insert( std::pair<std::string,TH1F*>(prefix+"Asm_MM", createHist1D(prefix+"Asm_MM",";A_{MM};",100,-5,5)) );
  
    
  // details of the jet kinematics
//    if(m_debug) Info("EtaIntercalibrationHists::initialize()", "adding kinematic plots");
//    m_jetPx     = book(m_name, "jetPx",     "jet Px [GeV]",     120, 0, 1000);
//    m_jetPy     = book(m_name, "jetPy",     "jet Py [GeV]",     120, 0, 1000);
//    m_jetPz     = book(m_name, "jetPz",     "jet Pz [GeV]",     120, 0, 4000);


  // N leading jets
//  if( m_infoSwitch->m_numLeadingJets > 0 ){
//    std::stringstream jetNum;
//    std::stringstream jetTitle;
//    for(int iJet=0; iJet < m_infoSwitch->m_numLeadingJets; ++iJet){
//      jetNum << iJet;
//
//      jetTitle << iJet+1;
//      switch(iJet)
//	{
//	case 0:
//	  jetTitle << "^{st}";
//	  break;
//	case 1:
//	  jetTitle << "^{nd}";
//	  break;
//	case 2:
//	  jetTitle << "^{rd}";
//	  break;
//	default:
//	  jetTitle << "^{th}";
//	  break;
//	}
//
//      m_NjetsPt.push_back(       book(m_name, ("jetPt_jet"+jetNum.str()),       jetTitle.str()+" jet p_{T} [GeV]" ,120,            0,       3000. ) );
//      m_NjetsEta.push_back(      book(m_name, ("jetEta_jet"+jetNum.str()),      jetTitle.str()+" jet #eta"        , 80,           -4,           4 ) );
//      m_NjetsPhi.push_back(      book(m_name, ("jetPhi_jet"+jetNum.str()),      jetTitle.str()+" jet Phi"         ,120, -TMath::Pi(), TMath::Pi() ) );
//      m_NjetsM.push_back(        book(m_name, ("jetMass_jet"+jetNum.str()),     jetTitle.str()+" jet Mass [GeV]"  ,120,            0,         400 ) );
//      m_NjetsE.push_back(        book(m_name, ("jetEnergy_jet"+jetNum.str()),   jetTitle.str()+" jet Energy [GeV]",120,            0,       4000. ) );
//      m_NjetsRapidity.push_back( book(m_name, ("jetRapidity_jet"+jetNum.str()), jetTitle.str()+" jet Rapidity"    ,120,          -10,          10 ) );
//      jetNum.str("");
//      jetTitle.str("");
//    }//for iJet
//  }

  return StatusCode::SUCCESS;
}

StatusCode EtaIntercalibrationHists::execute( const xAOD::JetContainer* jets, float eventWeight, int pvLoc ) {
    
  //single jet histograms
  for( auto jet_itr : *jets ) {
    RETURN_CHECK("EtaIntercalibrationHists::execute()", this->execute( jet_itr, eventWeight, pvLoc ), "");
  }
   
    
  std::string prefix = m_name+"_";

  //get out if there are fewer than two jets
  if (jets->size()<2) return StatusCode::SUCCESS;
      
  const xAOD::Jet* jet1 = jets->at(0);
  const xAOD::Jet* jet2 = jets->at(1);
    
  //float Dphijj = deltaPhi(jet1,jet2);
  float pTavg = (jet1->pt()+jet2->pt())/2.0;
  m_hists1D[ prefix+"pTavg" ]->Fill(pTavg/1e3);

  //histograms jet by jet
//  if( m_infoSwitch->m_numLeadingJets > 0){
//
//    int numJets = std::min( m_infoSwitch->m_numLeadingJets, (int)jets->size() );
//    for(int iJet=0; iJet < numJets; ++iJet){
//      m_NjetsPt.at(iJet)->        Fill( jets->at(iJet)->pt()/1e3,   eventWeight);
//      m_NjetsEta.at(iJet)->       Fill( jets->at(iJet)->eta(),      eventWeight);
//      m_NjetsPhi.at(iJet)->       Fill( jets->at(iJet)->phi(),      eventWeight);
//      m_NjetsM.at(iJet)->         Fill( jets->at(iJet)->m()/1e3,    eventWeight);
//      m_NjetsE.at(iJet)->         Fill( jets->at(iJet)->e()/1e3,    eventWeight);
//      m_NjetsRapidity.at(iJet)->  Fill( jets->at(iJet)->rapidity(), eventWeight);
//    }
//  }

  return StatusCode::SUCCESS;
}

StatusCode EtaIntercalibrationHists::execute( const xAOD::Jet* jet, float eventWeight, int /*pvLoc*/ ) {

  if(m_debug) std::cout << "in execute jet by jet " <<std::endl;

  //basic
  m_jetPt ->        Fill( jet->pt()/1e3,    eventWeight );
  m_jetEta->        Fill( jet->eta(),       eventWeight );
  m_jetPhi->        Fill( jet->phi(),       eventWeight );
  m_jetM->          Fill( jet->m()/1e3,     eventWeight );
  m_jetE->          Fill( jet->e()/1e3,     eventWeight );
  m_jetRapidity->   Fill( jet->rapidity(),  eventWeight );

  return StatusCode::SUCCESS;
}

//  // kinematic
//  if( m_infoSwitch->m_kinematic ) {
//    m_jetPx->  Fill( jet->px()/1e3,  eventWeight );
//    m_jetPy->  Fill( jet->py()/1e3,  eventWeight );
//    m_jetPz->  Fill( jet->pz()/1e3,  eventWeight );
//  } // fillKinematic
//
//  // clean
//  if( m_infoSwitch->m_clean ) {
//
//    static SG::AuxElement::ConstAccessor<float> jetTime ("Timing");
//    if( jetTime.isAvailable( *jet ) ) {
//      m_jetTime ->  Fill( jetTime( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> LArQuality ("LArQuality");
//    if( LArQuality.isAvailable( *jet ) ) {
//      m_LArQuality ->  Fill( LArQuality( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> hecq ("HECQuality");
//    if( hecq.isAvailable( *jet ) ) {
//      m_hecq ->  Fill( hecq( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> negE ("NegativeE");
//    if( negE.isAvailable( *jet ) ) {
//      m_negE ->  Fill( negE( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> avLArQF ("AverageLArQF");
//    if( avLArQF.isAvailable( *jet ) ) {
//      m_avLArQF ->  Fill( avLArQF( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> bchCorrCell ("BchCorrCell");
//    if( bchCorrCell.isAvailable( *jet ) ) {
//      m_bchCorrCell ->  Fill( bchCorrCell( *jet ), eventWeight );
//    }
//
//    // 0062       N90Cells?
//    static SG::AuxElement::ConstAccessor<float> N90Const ("N90Constituents");
//    if( N90Const.isAvailable( *jet ) ) {
//      m_N90Const ->  Fill( N90Const( *jet ), eventWeight );
//    }
//
//
// // 0030       LArBadHVEnergy,
// // 0031       LArBadHVRatio,
// // 0035       NumTowers,
// // 0057       isBadLoose,
// // 0058       isBadMedium,
// // 0059       isBadTight,
// // 0060       isUgly,
// // 0063       OotFracClusters10,
// // 0064       OotFracClusters5,
// // 0065       OotFracCells5,
// // 0066       OotFracCells10,
//
//  } // fillClean
//
//  // energy
//  if( m_infoSwitch->m_energy ) {
//
//    static SG::AuxElement::ConstAccessor<float> HECf ("HECFrac");
//    if( HECf.isAvailable( *jet ) ) {
//      m_HECf ->  Fill( HECf( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> EMf ("EMFrac");
//    if( EMf.isAvailable( *jet ) ) {
//      m_EMf ->  Fill( EMf( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> centroidR ("CentroidR");
//    if( centroidR.isAvailable( *jet ) ) {
//      m_centroidR ->  Fill( centroidR( *jet ), eventWeight );
//    }
//
//    /*
//
//    static SG::AuxElement::ConstAccessor<float> samplingMax ("SamplingMax");
//    if( samplingMax.isAvailable( *jet ) ) {
//      m_samplingMax ->  Fill( samplingMax( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> ePerSamp ("EnergyPerSampling");
//    if( ePerSamp.isAvailable( *jet ) ) {
//      m_ePerSamp ->  Fill( ePerSamp( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> fracSampMax ("FracSamplingMax");
//    if( fracSampMax.isAvailable( *jet ) ) {
//      m_fracSampMax ->  Fill( fracSampMax( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> lowEtFrac ("LowEtConstituentsFrac");
//    if( lowEtFrac.isAvailable( *jet ) ) {
//      m_lowEtFrac ->  Fill( lowEtFrac( *jet ), eventWeight );
//    }
//
// // 0036       Offset,
// // 0037       OriginIndex    ,
//    */
//
//  }
//
//  if( m_infoSwitch->m_layer ){
//    static SG::AuxElement::ConstAccessor< vector<float> > ePerSamp ("EnergyPerSampling");
//    if( ePerSamp.isAvailable( *jet ) ) {
//      vector<float> ePerSampVals = ePerSamp( *jet );
//      float jetE = jet->e();
//      m_PreSamplerB -> Fill( ePerSampVals.at(0) / jetE );
//      m_EMB1        -> Fill( ePerSampVals.at(1) / jetE );
//      m_EMB2        -> Fill( ePerSampVals.at(2) / jetE );
//      m_EMB3        -> Fill( ePerSampVals.at(3) / jetE );
//      m_PreSamplerE -> Fill( ePerSampVals.at(4) / jetE );
//      m_EME1        -> Fill( ePerSampVals.at(5) / jetE );
//      m_EME2        -> Fill( ePerSampVals.at(6) / jetE );
//      m_EME3        -> Fill( ePerSampVals.at(7) / jetE );
//      m_HEC0        -> Fill( ePerSampVals.at(8) / jetE );
//      m_HEC1        -> Fill( ePerSampVals.at(9) / jetE );
//      m_HEC2        -> Fill( ePerSampVals.at(10) / jetE );
//      m_HEC3        -> Fill( ePerSampVals.at(11) / jetE );
//      m_TileBar0    -> Fill( ePerSampVals.at(12) / jetE );
//      m_TileBar1    -> Fill( ePerSampVals.at(13) / jetE );
//      m_TileBar2    -> Fill( ePerSampVals.at(14) / jetE );
//      m_TileGap1    -> Fill( ePerSampVals.at(15) / jetE );
//      m_TileGap2    -> Fill( ePerSampVals.at(16) / jetE );
//      m_TileGap3    -> Fill( ePerSampVals.at(17) / jetE );
//      m_TileExt0    -> Fill( ePerSampVals.at(18) / jetE );
//      m_TileExt1    -> Fill( ePerSampVals.at(19) / jetE );
//      m_TileExt2    -> Fill( ePerSampVals.at(20) / jetE );
//      m_FCAL0       -> Fill( ePerSampVals.at(21) / jetE );
//      m_FCAL1       -> Fill( ePerSampVals.at(22) / jetE );
//      m_FCAL2       -> Fill( ePerSampVals.at(23) / jetE );
//    }
//  }
//
//
//
//  // area
//  /*
//  if ( m_fillArea ) {
//
//    static SG::AuxElement::ConstAccessor<int> actArea ("ActiveArea");
//    if( actArea.isAvailable( *jet ) ) {
//      m_actArea ->  Fill( actArea( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<int> voroniA ("VoronoiArea");
//    if( voroniA.isAvailable( *jet ) ) {
//      m_voroniA ->  Fill( voroniA( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<int> voroniAE ("VoronoiAreaE");
//    if( voroniAE.isAvailable( *jet ) ) {
//      m_voroniAE ->  Fill( voroniAE( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<int> voroniAPx ("VoronoiAreaPx");
//    if( voroniAPx.isAvailable( *jet ) ) {
//      m_voroniAPx ->  Fill( voroniAPx( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<int> voroniAPy ("CentroidR");
//    if( voroniAPy.isAvailable( *jet ) ) {
//      m_voroniAPy ->  Fill( voroniAPy( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<int> voroniAPz ("CentroidR");
//    if( voroniAPz.isAvailable( *jet ) ) {
//      m_voroniAPz ->  Fill( voroniAPz( *jet ), eventWeight );
//    }
//
//  }
//  */
//
//
//  /*
//  // tracks
//  if ( m_fillTrks ) {
// // 0040       TrackMF,
// // 0041       TrackMFindex,
// // 0049       WIDTH,
// // 0067       NumTrkPt1000,
// // 0068       NumTrkPt500,
// // 0069       SumPtTrkPt1000,
// // 0070       SumPtTrkPt500,
// // 0071       TrackWidthPt1000,
// // 0072       TrackWidthPt500,
// // 0012       GhostTrackCount,
// // 0011       GhostMuonSegmentCount,
// // 0027       HighestJVFVtx,
// } */
//
//
//
//  /*
//  // isolation
//  if ( m_fillIso ) {
// // 0024       IsoKR20Par,
// // 0025       IsoKR20Perp,
//  }
//  */
//
//  /*
//  // substructure
//  if( m_fillSubstructure) {
//    // 0029       KtDR,
//    static SG::AuxElement::ConstAccessor<int> ktDR ("KtDR");
//    if( ktDR.isAvailable( *jet ) ) {
//      m_ktDR ->  Fill( ktDR( *jet ), eventWeight );
//    }
// // 0050       YFlip12,
// // 0051       YFlip13,
// // 0074       Tau1,
// // 0075       Tau2,
// // 0076       Tau3,
// // 0077       Dip12,
// // 0078       Dip13,
// // 0079       Dip23,
// // 0080       DipExcl12,
// //
// // 0081       ThrustMin,
// // 0082       ThrustMaj,
// //
// // 0083       FoxWolfram0,
// // 0084       FoxWolfram1,
// // 0085       FoxWolfram2,
// // 0086       FoxWolfram3,
// // 0087       FoxWolfram4,
// //
// // 0088       Sphericity,
// // 0089       Aplanarity,
//  }
//
//  */
//
//  // truth
// // 0073       PtTruth,
// // 0013       GhostTruthParticleCount,
// // 0028       JetLabel,
// // 0042       TruthMF,
// // 0043       TruthMFindex,
//
//  if( m_infoSwitch->m_truth ) {
//
//    static SG::AuxElement::ConstAccessor<int> TruthLabelID ("TruthLabelID");
//    if( TruthLabelID.isAvailable( *jet ) ) {
//      m_truthLabelID ->  Fill( TruthLabelID( *jet ), eventWeight );
//    }else{
//      static SG::AuxElement::ConstAccessor<int> PartonTruthLabelID ("PartonTruthLabelID");
//      if( PartonTruthLabelID.isAvailable( *jet ) ) {
//	m_truthLabelID ->  Fill( PartonTruthLabelID( *jet ), eventWeight );
//      }
//    }
//
//    static SG::AuxElement::ConstAccessor<int> TruthCount ("TruthCount");
//    if( TruthCount.isAvailable( *jet ) ) {
//      m_truthCount ->  Fill( TruthCount( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> TruthPt ("TruthPt");
//    if( TruthPt.isAvailable( *jet ) ) {
//      m_truthPt ->  Fill( TruthPt( *jet )/1000, eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> TruthLabelDeltaR_B ("TruthLabelDeltaR_B");
//    if( TruthLabelDeltaR_B.isAvailable( *jet ) ) {
//      m_truthDr_B ->  Fill( TruthLabelDeltaR_B( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> TruthLabelDeltaR_C ("TruthLabelDeltaR_C");
//    if( TruthLabelDeltaR_C.isAvailable( *jet ) ) {
//      m_truthDr_C ->  Fill( TruthLabelDeltaR_C( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> TruthLabelDeltaR_T ("TruthLabelDeltaR_T");
//    if( TruthLabelDeltaR_T.isAvailable( *jet ) ) {
//      m_truthDr_T ->  Fill( TruthLabelDeltaR_T( *jet ), eventWeight );
//    }
//
//  }
//
//
//  if( m_infoSwitch->m_truthDetails ) {
//
//    //
//    // B-Hadron Details
//    //
//    static SG::AuxElement::ConstAccessor<int> GhostBHadronsFinalCount ("GhostBHadronsFinalCount");
//    if( GhostBHadronsFinalCount.isAvailable( *jet ) ) {
//      m_truthCount_BhadFinal ->  Fill( GhostBHadronsFinalCount( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<int> GhostBHadronsInitialCount ("GhostBHadronsInitialCount");
//    if( GhostBHadronsInitialCount.isAvailable( *jet ) ) {
//      m_truthCount_BhadInit ->  Fill( GhostBHadronsInitialCount( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<int> GhostBQuarksFinalCount ("GhostBQuarksFinalCount");
//    if( GhostBQuarksFinalCount.isAvailable( *jet ) ) {
//      m_truthCount_BQFinal ->  Fill( GhostBQuarksFinalCount( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> GhostBHadronsFinalPt ("GhostBHadronsFinalPt");
//    if( GhostBHadronsFinalPt.isAvailable( *jet ) ) {
//      m_truthPt_BhadFinal ->  Fill( GhostBHadronsFinalPt( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> GhostBHadronsInitialPt ("GhostBHadronsInitialPt");
//    if( GhostBHadronsInitialPt.isAvailable( *jet ) ) {
//      m_truthPt_BhadInit ->  Fill( GhostBHadronsInitialPt( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> GhostBQuarksFinalPt ("GhostBQuarksFinalPt");
//    if( GhostBQuarksFinalPt.isAvailable( *jet ) ) {
//      m_truthPt_BQFinal ->  Fill( GhostBQuarksFinalPt( *jet ), eventWeight );
//    }
//
//
//    //
//    // C-Hadron Details
//    //
//    static SG::AuxElement::ConstAccessor<int> GhostCHadronsFinalCount ("GhostCHadronsFinalCount");
//    if( GhostCHadronsFinalCount.isAvailable( *jet ) ) {
//      m_truthCount_ChadFinal ->  Fill( GhostCHadronsFinalCount( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<int> GhostCHadronsInitialCount ("GhostCHadronsInitialCount");
//    if( GhostCHadronsInitialCount.isAvailable( *jet ) ) {
//      m_truthCount_ChadInit ->  Fill( GhostCHadronsInitialCount( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<int> GhostCQuarksFinalCount ("GhostCQuarksFinalCount");
//    if( GhostCQuarksFinalCount.isAvailable( *jet ) ) {
//      m_truthCount_CQFinal ->  Fill( GhostCQuarksFinalCount( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> GhostCHadronsFinalPt ("GhostCHadronsFinalPt");
//    if( GhostCHadronsFinalPt.isAvailable( *jet ) ) {
//      m_truthPt_ChadFinal ->  Fill( GhostCHadronsFinalPt( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> GhostCHadronsInitialPt ("GhostCHadronsInitialPt");
//    if( GhostCHadronsInitialPt.isAvailable( *jet ) ) {
//      m_truthPt_ChadInit ->  Fill( GhostCHadronsInitialPt( *jet ), eventWeight );
//    }
//
//    static SG::AuxElement::ConstAccessor<float> GhostCQuarksFinalPt ("GhostCQuarksFinalPt");
//    if( GhostCQuarksFinalPt.isAvailable( *jet ) ) {
//      m_truthPt_CQFinal ->  Fill( GhostCQuarksFinalPt( *jet ), eventWeight );
//    }
//
//
//    //
//    // Tau Details
//    //
//    static SG::AuxElement::ConstAccessor<int> GhostTausFinalCount ("GhostTausFinalCount");
//    if( GhostTausFinalCount.isAvailable( *jet ) ) {
//      m_truthCount_TausFinal ->  Fill( GhostTausFinalCount( *jet ), eventWeight );
//    }
//
//
//    static SG::AuxElement::ConstAccessor<float> GhostTausFinalPt ("GhostTausFinalPt");
//    if( GhostTausFinalPt.isAvailable( *jet ) ) {
//      m_truthPt_TausFinal ->  Fill( GhostTausFinalPt( *jet ), eventWeight );
//    }
//
//
//  }
//
//  //
//  // BTagging
//  //
//  if( m_infoSwitch->m_flavTag || m_infoSwitch->m_flavTagHLT ) {
//
//    const xAOD::BTagging *btag_info(0);
//    if(m_infoSwitch->m_flavTag){
//      btag_info = jet->btagging();
//    }else if(m_infoSwitch->m_flavTagHLT){
//      btag_info = jet->auxdata< const xAOD::BTagging* >("HLTBTag");
//    }
//
//    double MV2c00 = -99;
//    double MV2c10 = -99;
//    double MV2c20 = -99;
//    btag_info->MVx_discriminant("MV2c00", MV2c00);
//    btag_info->MVx_discriminant("MV2c10", MV2c10);
//    btag_info->MVx_discriminant("MV2c20", MV2c20);
//    m_MV2c00 ->  Fill( MV2c00, eventWeight );
//    m_MV2c10 ->  Fill( MV2c10, eventWeight );
//    m_MV2c20 ->  Fill( MV2c20, eventWeight );
//
//    static SG::AuxElement::ConstAccessor<double> SV0_significance3DAcc ("SV0_significance3D");
//    if ( SV0_significance3DAcc.isAvailable(*btag_info) ) {
//      m_SV0             ->  Fill( btag_info->SV0_significance3D() , eventWeight );
//      m_SV1             ->  Fill( btag_info->SV1_loglikelihoodratio() , eventWeight );
//      m_IP2D            ->  Fill( btag_info->IP2D_loglikelihoodratio() , eventWeight );
//      m_IP3D            ->  Fill( btag_info->IP3D_loglikelihoodratio() , eventWeight );
//      m_COMB            ->  Fill( btag_info->SV1_loglikelihoodratio() + btag_info->IP3D_loglikelihoodratio() , eventWeight );
//      m_JetFitter       ->  Fill( btag_info->JetFitter_loglikelihoodratio() , eventWeight );
//      m_JetFitterCombNN ->  Fill( btag_info->JetFitterCombNN_loglikelihoodratio() , eventWeight );
//    }
//
//    if(m_infoSwitch->m_jetFitterDetails){
//      static SG::AuxElement::ConstAccessor< int   > jf_nVTXAcc       ("JetFitter_nVTX");
//      static SG::AuxElement::ConstAccessor< int   > jf_nSingleTracks ("JetFitter_nSingleTracks");
//      static SG::AuxElement::ConstAccessor< int   > jf_nTracksAtVtx  ("JetFitter_nTracksAtVtx");
//      static SG::AuxElement::ConstAccessor< float > jf_mass          ("JetFitter_mass");
//      static SG::AuxElement::ConstAccessor< float > jf_energyFraction("JetFitter_energyFraction");
//      static SG::AuxElement::ConstAccessor< float > jf_significance3d("JetFitter_significance3d");
//      static SG::AuxElement::ConstAccessor< float > jf_deltaeta      ("JetFitter_deltaeta");
//      static SG::AuxElement::ConstAccessor< float > jf_deltaphi      ("JetFitter_deltaphi");
//      static SG::AuxElement::ConstAccessor< int   > jf_N2Tpar        ("JetFitter_N2Tpair");
//      static SG::AuxElement::ConstAccessor< double > jf_pb           ("JetFitter_pb");
//      static SG::AuxElement::ConstAccessor< double > jf_pc           ("JetFitter_pc");
//      static SG::AuxElement::ConstAccessor< double > jf_pu           ("JetFitter_pu");
//
//      if(jf_nVTXAcc.isAvailable       (*btag_info)) m_jf_nVTX           ->Fill(jf_nVTXAcc       (*btag_info), eventWeight);
//      if(jf_nSingleTracks.isAvailable (*btag_info)) m_jf_nSingleTracks  ->Fill(jf_nSingleTracks (*btag_info), eventWeight);
//      if(jf_nTracksAtVtx.isAvailable  (*btag_info)) m_jf_nTracksAtVtx   ->Fill(jf_nTracksAtVtx  (*btag_info), eventWeight);
//      if(jf_mass.isAvailable          (*btag_info)) m_jf_mass           ->Fill(jf_mass          (*btag_info)/1000, eventWeight);
//      if(jf_energyFraction.isAvailable(*btag_info)) m_jf_energyFraction ->Fill(jf_energyFraction(*btag_info), eventWeight);
//      if(jf_significance3d.isAvailable(*btag_info)) m_jf_significance3d ->Fill(jf_significance3d(*btag_info), eventWeight);
//      if(jf_deltaeta.isAvailable      (*btag_info)) m_jf_deltaeta       ->Fill(jf_deltaeta      (*btag_info), eventWeight);
//      if(jf_deltaphi.isAvailable      (*btag_info)) m_jf_deltaphi       ->Fill(jf_deltaphi      (*btag_info), eventWeight);
//      if(jf_N2Tpar.isAvailable        (*btag_info)) m_jf_N2Tpar         ->Fill(jf_N2Tpar        (*btag_info), eventWeight);
//      if(jf_pb.isAvailable            (*btag_info)) m_jf_pb             ->Fill(jf_pb            (*btag_info), eventWeight);
//      if(jf_pu.isAvailable            (*btag_info)) m_jf_pu             ->Fill(jf_pu            (*btag_info), eventWeight);
//    }
//
//
//    if(m_infoSwitch->m_svDetails){
//
//      //
//      // SV0
//      //
//
//      /// @brief SV0 : Number of good tracks in vertex
//      static SG::AuxElement::ConstAccessor< int   >   sv0_NGTinSvxAcc     ("SV0_NGTinSvx");
//      // @brief SV0 : Number of 2-track pairs
//      static SG::AuxElement::ConstAccessor< int   >   sv0_N2TpairAcc      ("SV0_N2Tpair");
//      /// @brief SV0 : vertex mass
//      static SG::AuxElement::ConstAccessor< float   > sv0_masssvxAcc      ("SV0_masssvx");
//      /// @brief SV0 : energy fraction
//      static SG::AuxElement::ConstAccessor< float   > sv0_efracsvxAcc     ("SV0_efracsvx");                                                                  	/// @brief SV0 : 3D vertex significance
//      static SG::AuxElement::ConstAccessor< float   > sv0_normdistAcc     ("SV0_normdist");
//
//      if(sv0_NGTinSvxAcc .isAvailable(*btag_info)) m_sv0_NGTinSvx -> Fill( sv0_NGTinSvxAcc (*btag_info), eventWeight);
//      if(sv0_N2TpairAcc  .isAvailable(*btag_info)) m_sv0_N2Tpair  -> Fill( sv0_N2TpairAcc  (*btag_info), eventWeight);
//      if(sv0_masssvxAcc  .isAvailable(*btag_info)) m_sv0_massvx   -> Fill( sv0_masssvxAcc  (*btag_info)/1000, eventWeight);
//      if(sv0_efracsvxAcc .isAvailable(*btag_info)) m_sv0_efracsvx -> Fill( sv0_efracsvxAcc (*btag_info), eventWeight);
//      if(sv0_normdistAcc .isAvailable(*btag_info)) m_sv0_normdist -> Fill( sv0_normdistAcc (*btag_info), eventWeight);
//
//      //
//      // SV1
//      //
//
//      /// @brief SV1 : Number of good tracks in vertex
//      static SG::AuxElement::ConstAccessor< int   >   sv1_NGTinSvxAcc     ("SV1_NGTinSvx");
//      // @brief SV1 : Number of 2-track pairs
//      static SG::AuxElement::ConstAccessor< int   >   sv1_N2TpairAcc      ("SV1_N2Tpair");
//      /// @brief SV1 : vertex mass
//      static SG::AuxElement::ConstAccessor< float   > sv1_masssvxAcc      ("SV1_masssvx");
//      /// @brief SV1 : energy fraction
//      static SG::AuxElement::ConstAccessor< float   > sv1_efracsvxAcc     ("SV1_efracsvx");                                                                  	/// @brief SV1 : 3D vertex significance
//      static SG::AuxElement::ConstAccessor< float   > sv1_normdistAcc     ("SV1_normdist");
//
//      if(sv1_NGTinSvxAcc .isAvailable(*btag_info)) m_sv1_NGTinSvx -> Fill( sv1_NGTinSvxAcc (*btag_info), eventWeight);
//      if(sv1_N2TpairAcc  .isAvailable(*btag_info)) m_sv1_N2Tpair  -> Fill( sv1_N2TpairAcc  (*btag_info), eventWeight);
//      if(sv1_masssvxAcc  .isAvailable(*btag_info)) m_sv1_massvx   -> Fill( sv1_masssvxAcc  (*btag_info)/1000, eventWeight);
//      if(sv1_efracsvxAcc .isAvailable(*btag_info)) m_sv1_efracsvx -> Fill( sv1_efracsvxAcc (*btag_info), eventWeight);
//      if(sv1_normdistAcc .isAvailable(*btag_info)) m_sv1_normdist -> Fill( sv1_normdistAcc (*btag_info), eventWeight);
//
//    }
//
//
//    if(m_infoSwitch->m_ipDetails){
//
//      //
//      // IP2D
//      //
//
//      /// @brief IP2D: track grade
//      static SG::AuxElement::ConstAccessor< vector<int>   >   IP2D_gradeOfTracksAcc     ("IP2D_gradeOfTracks");
//      /// @brief IP2D : tracks from V0
//      static SG::AuxElement::ConstAccessor< vector<bool>   >  IP2D_flagFromV0ofTracksAcc("IP2D_flagFromV0ofTracks");
//      /// @brief IP2D : d0 value with respect to primary vertex
//      static SG::AuxElement::ConstAccessor< vector<float>   > IP2D_valD0wrtPVofTracksAcc("IP2D_valD0wrtPVofTracks");
//      /// @brief IP2D : d0 significance with respect to primary vertex
//      static SG::AuxElement::ConstAccessor< vector<float>   > IP2D_sigD0wrtPVofTracksAcc("IP2D_sigD0wrtPVofTracks");
//      /// @brief IP2D : track contribution to B likelihood
//      static SG::AuxElement::ConstAccessor< vector<float>   > IP2D_weightBofTracksAcc   ("IP2D_weightBofTracks");
//      /// @brief IP2D : track contribution to C likelihood
//      static SG::AuxElement::ConstAccessor< vector<float>   > IP2D_weightCofTracksAcc   ("IP2D_weightCofTracks");
//      /// @brief IP2D : track contribution to U likelihood
//      static SG::AuxElement::ConstAccessor< vector<float>   > IP2D_weightUofTracksAcc   ("IP2D_weightUofTracks");
//
//      if(IP2D_gradeOfTracksAcc .isAvailable(*btag_info)){
//	unsigned int nIP2DTracks = IP2D_gradeOfTracksAcc(*btag_info).size();
//	m_nIP2DTracks -> Fill( nIP2DTracks, eventWeight);
//	for(int grade : IP2D_gradeOfTracksAcc(*btag_info))        m_IP2D_gradeOfTracks->Fill(grade, eventWeight);
//      }
//
//      if(IP2D_flagFromV0ofTracksAcc .isAvailable(*btag_info)){
//	for(bool flag : IP2D_flagFromV0ofTracksAcc(*btag_info))   m_IP2D_flagFromV0ofTracks->Fill(flag, eventWeight);
//      }
//
//      if(IP2D_valD0wrtPVofTracksAcc .isAvailable(*btag_info)){
//	for(float d0 : IP2D_valD0wrtPVofTracksAcc(*btag_info))    m_IP2D_valD0wrtPVofTracks->Fill(d0, eventWeight);
//      }
//
//      if(IP2D_sigD0wrtPVofTracksAcc .isAvailable(*btag_info)){
//	for(float d0Sig : IP2D_sigD0wrtPVofTracksAcc(*btag_info)) {
//	  m_IP2D_sigD0wrtPVofTracks  ->Fill(d0Sig, eventWeight);
//	  m_IP2D_sigD0wrtPVofTracks_l->Fill(d0Sig, eventWeight);
//	}
//      }
//
//      if(IP2D_weightBofTracksAcc .isAvailable(*btag_info)){
//	for(float weightB : IP2D_weightBofTracksAcc(*btag_info))  m_IP2D_weightBofTracks->Fill(weightB, eventWeight);
//      }
//
//      if(IP2D_weightCofTracksAcc .isAvailable(*btag_info)){
//	for(float weightC : IP2D_weightCofTracksAcc(*btag_info))  m_IP2D_weightCofTracks->Fill(weightC, eventWeight);
//      }
//
//      if(IP2D_weightUofTracksAcc .isAvailable(*btag_info)){
//	for(float weightU : IP2D_weightUofTracksAcc(*btag_info))  m_IP2D_weightUofTracks->Fill(weightU, eventWeight);
//      }
//
//      //
//      // IP3D
//      //
//
//      /// @brief IP3D: track grade
//      static SG::AuxElement::ConstAccessor< vector<int>   >   IP3D_gradeOfTracksAcc     ("IP3D_gradeOfTracks");
//      /// @brief IP3D : tracks from V0
//      static SG::AuxElement::ConstAccessor< vector<bool>   >  IP3D_flagFromV0ofTracksAcc("IP3D_flagFromV0ofTracks");
//      /// @brief IP3D : d0 value with respect to primary vertex
//      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_valD0wrtPVofTracksAcc("IP3D_valD0wrtPVofTracks");
//      /// @brief IP3D : d0 significance with respect to primary vertex
//      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_sigD0wrtPVofTracksAcc("IP3D_sigD0wrtPVofTracks");
//      /// @brief IP3D : z0 value with respect to primary vertex
//      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_valZ0wrtPVofTracksAcc("IP3D_valZ0wrtPVofTracks");
//      /// @brief IP3D : z0 significance with respect to primary vertex
//      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_sigZ0wrtPVofTracksAcc("IP3D_sigZ0wrtPVofTracks");
//      /// @brief IP3D : track contribution to B likelihood
//      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_weightBofTracksAcc   ("IP3D_weightBofTracks");
//      /// @brief IP3D : track contribution to C likelihood
//      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_weightCofTracksAcc   ("IP3D_weightCofTracks");
//      /// @brief IP3D : track contribution to U likelihood
//      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_weightUofTracksAcc   ("IP3D_weightUofTracks");
//
//      if(IP3D_gradeOfTracksAcc .isAvailable(*btag_info)){
//	unsigned int nIP3DTracks = IP3D_gradeOfTracksAcc(*btag_info).size();
//	m_nIP3DTracks -> Fill( nIP3DTracks, eventWeight);
//	for(int grade : IP3D_gradeOfTracksAcc(*btag_info))        m_IP3D_gradeOfTracks->Fill(grade, eventWeight);
//      }
//
//      if(IP3D_flagFromV0ofTracksAcc .isAvailable(*btag_info)){
//	for(bool flag : IP3D_flagFromV0ofTracksAcc(*btag_info))   m_IP3D_flagFromV0ofTracks->Fill(flag, eventWeight);
//      }
//
//      if(IP3D_valD0wrtPVofTracksAcc .isAvailable(*btag_info)){
//	for(float d0 : IP3D_valD0wrtPVofTracksAcc(*btag_info))    m_IP3D_valD0wrtPVofTracks->Fill(d0, eventWeight);
//      }
//
//      if(IP3D_sigD0wrtPVofTracksAcc .isAvailable(*btag_info)){
//	for(float d0Sig : IP3D_sigD0wrtPVofTracksAcc(*btag_info)){
//	  m_IP3D_sigD0wrtPVofTracks  ->Fill(d0Sig, eventWeight);
//	  m_IP3D_sigD0wrtPVofTracks_l->Fill(d0Sig, eventWeight);
//	}
//      }
//
//      if(IP3D_valZ0wrtPVofTracksAcc .isAvailable(*btag_info)){
//	for(float z0 : IP3D_valZ0wrtPVofTracksAcc(*btag_info))    m_IP3D_valZ0wrtPVofTracks->Fill(z0, eventWeight);
//      }
//
//      if(IP3D_sigZ0wrtPVofTracksAcc .isAvailable(*btag_info)){
//	for(float z0Sig : IP3D_sigZ0wrtPVofTracksAcc(*btag_info)){
//	  m_IP3D_sigZ0wrtPVofTracks  ->Fill(z0Sig, eventWeight);
//	  m_IP3D_sigZ0wrtPVofTracks_l->Fill(z0Sig, eventWeight);
//	}
//      }
//
//      if(IP3D_weightBofTracksAcc .isAvailable(*btag_info)){
//	for(float weightB : IP3D_weightBofTracksAcc(*btag_info))  m_IP3D_weightBofTracks->Fill(weightB, eventWeight);
//      }
//
//      if(IP3D_weightCofTracksAcc .isAvailable(*btag_info)){
//	for(float weightC : IP3D_weightCofTracksAcc(*btag_info))  m_IP3D_weightCofTracks->Fill(weightC, eventWeight);
//      }
//
//      if(IP3D_weightUofTracksAcc .isAvailable(*btag_info)){
//	for(float weightU : IP3D_weightUofTracksAcc(*btag_info))  m_IP3D_weightUofTracks->Fill(weightU, eventWeight);
//      }
//
//
//
//    }
//  }
//
//
//
//
//  /*
//  vector<float> chfs = jet->getAttribute< vector<float> >(xAOD::JetAttribute::SumPtTrkPt1000);
//  float chf(-1);
//  if( pvLoc >= 0 && pvLoc < (int)chfs.size() ) {
//    m_chf ->  Fill( chfs.at( pvLoc ) , eventWeight );
//  }
//  */
//
//
//  // testing
//  if( m_infoSwitch->m_resolution ) {
//    //float ghostTruthPt = jet->getAttribute( xAOD::JetAttribute::GhostTruthPt );
//    float ghostTruthPt = jet->auxdata< float >( "GhostTruthPt" );
//    m_jetGhostTruthPt -> Fill( ghostTruthPt/1e3, eventWeight );
//    float resolution = jet->pt()/ghostTruthPt - 1;
//    m_jetPt_vs_resolution -> Fill( jet->pt()/1e3, resolution, eventWeight );
//    m_jetGhostTruthPt_vs_resolution -> Fill( ghostTruthPt/1e3, resolution, eventWeight );
//  }
//
//  if( m_infoSwitch->m_substructure ){
//    static SG::AuxElement::ConstAccessor<float> Tau1("Tau1");
//    static SG::AuxElement::ConstAccessor<float> Tau2("Tau2");
//    static SG::AuxElement::ConstAccessor<float> Tau3("Tau3");
//    static SG::AuxElement::ConstAccessor<float> Tau1_wta("Tau1_wta");
//    static SG::AuxElement::ConstAccessor<float> Tau2_wta("Tau2_wta");
//    static SG::AuxElement::ConstAccessor<float> Tau3_wta("Tau3_wta");
//
//    if(Tau1.isAvailable(*jet)) m_tau1->Fill( Tau1(*jet), eventWeight );
//    if(Tau2.isAvailable(*jet)) m_tau2->Fill( Tau2(*jet), eventWeight );
//    if(Tau3.isAvailable(*jet)) m_tau3->Fill( Tau3(*jet), eventWeight );
//    if(Tau1.isAvailable(*jet) && Tau2.isAvailable(*jet)) m_tau21->Fill( Tau2(*jet)/Tau1(*jet), eventWeight );
//    if(Tau2.isAvailable(*jet) && Tau3.isAvailable(*jet)) m_tau32->Fill( Tau3(*jet)/Tau2(*jet), eventWeight );
//    if(Tau1_wta.isAvailable(*jet)) m_tau1_wta->Fill( Tau1_wta(*jet), eventWeight );
//    if(Tau2_wta.isAvailable(*jet)) m_tau2_wta->Fill( Tau2_wta(*jet), eventWeight );
//    if(Tau3_wta.isAvailable(*jet)) m_tau3_wta->Fill( Tau3_wta(*jet), eventWeight );
//    if(Tau1_wta.isAvailable(*jet) && Tau2_wta.isAvailable(*jet)) m_tau21_wta->Fill( Tau2_wta(*jet)/Tau1_wta(*jet), eventWeight );
//    if(Tau2_wta.isAvailable(*jet) && Tau3_wta.isAvailable(*jet)) m_tau32_wta->Fill( Tau3_wta(*jet)/Tau2_wta(*jet), eventWeight );
//
//    m_numConstituents->Fill( jet->numConstituents(), eventWeight );
//
//  }
//
//  return StatusCode::SUCCESS;
//}


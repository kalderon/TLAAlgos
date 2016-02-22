#ifndef TLAAlgos_ProcessTLAMiniTree_H
#define TLAAlgos_ProcessTLAMiniTree_H

#include <EventLoop/StatusCode.h>
#include <EventLoop/Algorithm.h>
#include <EventLoop/Worker.h>

// rootcore includes
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

//algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

// ROOT include(s):
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TLorentzVector.h"

#include <sstream>
#include <vector>

using namespace std;



class ProcessTLAMiniTree : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
  public:

    // GRL
    bool m_applyGRL;
    std::string m_GRLxml;

    //configuration variables
    bool m_debug;                     
    bool m_doTrigger;                 
    bool m_doTruthOnly;                 
    int  m_detailLevel;               
    float m_YStarCut;                 
    float m_YBoostCut;             
    bool m_doBlind;                
    bool m_doData;                 
    bool m_doPUReweight;           
    bool m_useWeighted;            
    bool m_doCleaning;             
    float m_etaCut;                
    float m_leadJetPtCut;          
    float m_lumi;                   // Lumi we are scaling to

    float m_sampleEvents;           //! MC Events in the sample we are processing


    // float cutValue;
    int m_eventCounter;     //!

    std::string m_name;
    float m_mcEventWeight;  //!

  private:

    GoodRunsListSelectionTool*   m_grl;       //!

    float m_ht; //!
    int m_NPV; //!
    float m_actualInteractionsPerCrossing; //!
    float m_averageInteractionsPerCrossing; //!
    int m_runNumber; //!
    int m_eventNumber; //!
    int m_lumiBlock; //!

    float m_weight; //!
    float m_weight_xs; //!
    float m_weight_pileup;    //!
    vector<float>* m_jet_pt; //!
    vector<float>* m_jet_eta; //!
    vector<float>* m_jet_phi; //!
    vector<float>* m_jet_E; //!
    vector<float>* m_jet_Timing; //!
    vector<float>* m_jet_GhostMuonSegmentCount; //!
    vector<int>*   m_jet_clean_passLooseBad; //!
    vector<string>* m_passedTriggers; //!

    vector<float>* m_gamma_pt; //!
    vector<float>* m_gamma_eta; //!
    vector<float>* m_gamma_phi; //!
    vector<float>* m_gamma_E; //!

    //
    // Jet Data
    //
    struct jetData{
      float pt;
      float eta;
      float phi;
      float E;
      float Timing;
      float MuonSegmentCount; 

      jetData(float m_pt, float m_eta, float m_phi, float m_E){ 
	      
	pt     = m_pt;
	eta    = m_eta;
	phi    = m_phi;
	E      = m_E;

	// Set directly out side of constructor
	Timing           = -99;
	MuonSegmentCount = -99;
      
      }

      TLorentzVector vec() const{
	TLorentzVector vec = TLorentzVector();
	vec.SetPtEtaPhiE(pt,eta,phi,E);
	return vec;
      }

    };

    //
    // Gamma Data
    //
    struct gammaData{
      float pt;
      float eta;
      float phi;
      float E;

      gammaData(float m_pt, float m_eta, float m_phi, float m_E){ 
	      
	pt     = m_pt;
	eta    = m_eta;
	phi    = m_phi;
	E      = m_E;
      
      }

    };


    //
    // Event Data
    //
    struct eventData{
      int             runNumber;
      int             eventNumber;
      vector<jetData> jets;
      vector<gammaData> gammas;
      int             NPV;
      float           mu_ave;
      float           mu_act;
      float           weight;
      vector<jetData>* l1_ROIs;      

      eventData(unsigned int m_runNumber, unsigned int m_eventNumber, 
		vector<float>* m_jet_pt, vector<float>* m_jet_eta, vector<float>* m_jet_phi, vector<float>* m_jet_E, 
		vector<float>* m_jet_Timing, vector<float>* m_jet_GhostMuonSegmentCount, 
		vector<float>* m_gamma_pt, vector<float>* m_gamma_eta, vector<float>* m_gamma_phi, vector<float>* m_gamma_E, 
		int m_NPV, float m_mu_ave, float m_mu_act, float& m_weight){
	runNumber = m_runNumber;
	eventNumber = m_eventNumber;

	for(unsigned int i =0; i < m_jet_pt->size(); ++i){
	  jetData thisJet = jetData(m_jet_pt->at(i), m_jet_eta->at(i), m_jet_phi->at(i), m_jet_E->at(i));

	  if(m_jet_Timing){
	    thisJet.Timing = m_jet_Timing->at(i);
	  }

	  if(m_jet_GhostMuonSegmentCount){
	    thisJet.MuonSegmentCount = m_jet_GhostMuonSegmentCount->at(i);
	  }

	    
	  jets.push_back(thisJet);
	}


	for(unsigned int i =0; i < m_gamma_pt->size(); ++i){
	  gammaData thisGamma = gammaData(m_gamma_pt->at(i), m_gamma_eta->at(i), m_gamma_phi->at(i), m_gamma_E->at(i));
	  gammas.push_back(thisGamma);
	}


	NPV	 = m_NPV;
	mu_ave   = m_mu_ave;
	mu_act   = m_mu_act;
	weight   = m_weight;
      }

      void dump() const{
	cout << "EventDump: run: " << runNumber << " event: " << eventNumber << endl;
	for(auto& jet : jets){
	  cout << "jet:  pt:" << jet.pt << " eta: " << jet.eta << " phi: " << jet.phi << " E:"<< jet.E <<endl;
	}
	for(auto& gamma : gammas){
	  cout << "jet:  pt:" << gamma.pt << " eta: " << gamma.eta << " phi: " << gamma.phi << " E:"<< gamma.E <<endl;
	}
      }

    };

    struct jetHists{
      
      int m_detailLevel;

      // Detail Level 0
      TH1F* h_pt;
      TH1F* h_pt_m;
      TH1F* h_pt_l;
      TH1F* h_eta;
      TH1F* h_phi;

      // Details
      TH2F* h_jettiming_vs_eta;
      TH1F* h_MuonSegments;
      TH1F* h_MuonSegments_l;
      TH1F* h_MuonSegments_vl;
      TProfile* h_MuonSegments_vs_eta;
      //TH2F* h_etaphi;
      //TH2F* h_etaphi_s;

      

      jetHists(std::string name, EL::Worker* wk, int detailLevel=0){
        h_pt  = new TH1F((name+"_Pt").c_str(), "jetPt;p_{T} [GeV];Entries",   100,  0, 500);
	wk->addOutput(h_pt);

        h_pt_m  = new TH1F((name+"_Pt_m").c_str(), "jetPt;p_{T} [GeV];Entries",   100,  0, 1000);
	wk->addOutput(h_pt_m);

        h_pt_l  = new TH1F((name+"_Pt_l").c_str(), "jetPt;p_{T} [GeV];Entries",   100,  0, 5000);
	wk->addOutput(h_pt_l);


	h_eta = new TH1F((name+"_Eta").c_str(),"jetEta;Eta;Entries",          100, -3,   3);
	wk->addOutput(h_eta);

	h_phi = new TH1F((name+"_Phi").c_str(),"jetPhi;Phi;Entries",          100, -3.2,   3.2);
	wk->addOutput(h_phi);

	m_detailLevel = detailLevel;
	
	if(m_detailLevel > 0){
	  h_jettiming_vs_eta = new TH2F((name+"_timing_vs_eta").c_str(),"jetTime_vs_eta;jet #eta;jetTime",     50, -3.2,   3.2, 100, -10,10);
	  wk->addOutput(h_jettiming_vs_eta);

	  h_MuonSegments = new TH1F((name+"_MuonSegs").c_str(),"jetMuonSegs;MuonSegments;Entries",          50, -0.5,   49.5);
	  wk->addOutput(h_MuonSegments);

	  h_MuonSegments_l = new TH1F((name+"_MuonSegs_l").c_str(),"jetMuonSegs;MuonSegments;Entries",          100, -0.5,   99.5);
	  wk->addOutput(h_MuonSegments_l);

	  h_MuonSegments_vl = new TH1F((name+"_MuonSegs_vl").c_str(),"jetMuonSegs;MuonSegments;Entries",          100, -0.5,   199.5);
	  wk->addOutput(h_MuonSegments_vl);

	  h_MuonSegments_vs_eta = new TProfile((name+"_MuonSegs_vs_eta").c_str(),"jetMuonSegs_vs_eta;jet #eta;MuonSegments",  50, -3.2, 3.2, -0.5, 49.5);
	  wk->addOutput(h_MuonSegments_vs_eta);
	}

	//h_etaphi = new TH2F((name+"_EtaPhi").c_str(),"jetEtaPhi;Eta;Phi",          100, -3,3,100,-3.2,   3.2);
	//wk->addOutput(h_etaphi);
	//
	//h_etaphi_s = new TH2F((name+"_EtaPhi_s").c_str(),"jetEtaPhi;Eta;Phi",          50, -3,3,50,-3.2,   3.2);
	//wk->addOutput(h_etaphi_s);
      }
      
      void Fill(const jetData& jet, const float& weight){
        h_pt       -> Fill(jet.pt,       weight);
        h_pt_l     -> Fill(jet.pt,       weight);
        h_pt_m     -> Fill(jet.pt,       weight);
	h_eta      -> Fill(jet.eta,      weight);
	h_phi      -> Fill(jet.phi,      weight);

	if(m_detailLevel > 0){
	  h_jettiming_vs_eta-> Fill(jet.eta, jet.Timing,   weight);
	  h_MuonSegments    -> Fill(jet.MuonSegmentCount,  weight);
	  h_MuonSegments_l  -> Fill(jet.MuonSegmentCount,  weight);
	  h_MuonSegments_vl -> Fill(jet.MuonSegmentCount,  weight);
	  h_MuonSegments_vs_eta  -> Fill(jet.eta, jet.MuonSegmentCount,  weight);
	}
	//h_etaphi   -> Fill(eta, phi, weight);
	//h_etaphi_s -> Fill(eta, phi, weight);
      }

    };

    struct gammaHists{
      
      int m_detailLevel;

      // Detail Level 0
      TH1F* h_pt;
      TH1F* h_pt_m;
      TH1F* h_pt_l;
      TH1F* h_eta;
      TH1F* h_phi;

      gammaHists(std::string name, EL::Worker* wk, int detailLevel=0){
        h_pt  = new TH1F((name+"_Pt").c_str(), "gammaPt;p_{T} [GeV];Entries",   100,  0, 500);
	wk->addOutput(h_pt);

        h_pt_m  = new TH1F((name+"_Pt_m").c_str(), "gammaPt;p_{T} [GeV];Entries",   100,  0, 1000);
	wk->addOutput(h_pt_m);

        h_pt_l  = new TH1F((name+"_Pt_l").c_str(), "gammaPt;p_{T} [GeV];Entries",   100,  0, 5000);
	wk->addOutput(h_pt_l);


	h_eta = new TH1F((name+"_Eta").c_str(),"gammaEta;Eta;Entries",          100, -3,   3);
	wk->addOutput(h_eta);

	h_phi = new TH1F((name+"_Phi").c_str(),"gammaPhi;Phi;Entries",          100, -3.2,   3.2);
	wk->addOutput(h_phi);

	m_detailLevel = detailLevel;
	
      }
      
      void Fill(const gammaData& gamma, const float& weight){
        h_pt       -> Fill(gamma.pt,       weight);
        h_pt_l     -> Fill(gamma.pt,       weight);
        h_pt_m     -> Fill(gamma.pt,       weight);
	h_eta      -> Fill(gamma.eta,      weight);
	h_phi      -> Fill(gamma.phi,      weight);
      }

    };



    struct eventHists{
      TH1F*     h_NPV;
      TH1F*     h_mu_ave;
      TH1F*     h_mu_act;
      TH1F*     h_nJet;
      TH1F*     h_nJet_l;
      TH1F*     h_nGamma;
      TH1F*     h_nGamma_l;

      //
      // jj
      //
      TH1F*     h_yStarjj;
      TH1F*     h_yBoostjj;
      TH1F*     h_dEtajj;
      TH1F*     h_dPhijj;
      TH1F*     h_dRjj;
      TH1F*     h_chijj;
      TH1F*     h_asymjj;
      TH1F*     h_mjj;
      TH1F*     h_ptjj;
      TH1F*     h_etajj;
      TH1F*     h_phijj;

      //
      // ISRjj
      //
      TH1F*     h_yStarISRjj;
      TH1F*     h_yBoostISRjj;
      TH1F*     h_dEtaISRjj;
      TH1F*     h_dPhiISRjj;
      TH1F*     h_dRISRjj;
      TH1F*     h_chiISRjj;
      TH1F*     h_asymISRjj;
      TH1F*     h_mISRjj;
      TH1F*     h_ptISRjj;
      TH1F*     h_etaISRjj;
      TH1F*     h_phiISRjj;


      //
      // LeadJet vs ISR
      //
      TH1F*     h_dRLeadJetISR;
      TH1F*     h_dRCloseJetISR;
      TH1F*     h_dPhiLeadJetISR;
      TH1F*     h_dEtaLeadJetISR;

      jetHists* h_jets;
      jetHists* h_leadJet;
      jetHists* h_sublJet;
      jetHists* h_thrdJet;
      float     m_etaCutVal;
      gammaHists* h_gammas;
      gammaHists* h_leadGamma;


      eventHists(std::string name, EL::Worker* wk, float etaCut = 2.8, int detailLevel = 0){
	

	//
	// Event Level
	//
	h_NPV      = book(wk, name, "NPV",       "NPV",           50,     -0.5,    49   );
	h_mu_ave   = book(wk, name, "mu_ave",    "mu_ave",        50,     -0.5,    49   );
	h_mu_act   = book(wk, name, "mu_act",    "mu_act",        50,     -0.5,    49   );
        h_nJet     = book(wk, name, "nJet",      "nJets",         10,     -0.5,     9.5 );
        h_nJet_l   = book(wk, name, "nJet_l",    "nJets",         20,     -0.5,    19.5 );
        h_nGamma   = book(wk, name, "nGamma",    "nGammas",       10,     -0.5,     9.5 );
        h_nGamma_l = book(wk, name, "nGamma_l",  "nGammas",       20,     -0.5,    19.5 );

	//
	// jj hists
	//
        h_yStarjj  = book(wk, name, "yStarjj",   "yStarjj",      100,     -3,       3   );
        h_yBoostjj = book(wk, name, "yBoostjj",  "yBoostjj",     100,     -3,       3   );
        h_dPhijj   = book(wk, name, "dPhijj",    "dPhijj",       100,     -3.2,     3.2 );
        h_dEtajj   = book(wk, name, "dEtajj",    "dEtajj",       100,     -5,       5   );
        h_dRjj     = book(wk, name, "dRjj",      "dRjj",         100,      0,       5   );
        h_asymjj   = book(wk, name, "asymjj",    "asymjj",       100,     -0.1,     1   );
        h_chijj    = book(wk, name, "chijj",     "chijj",        100,      0,      60   );
        h_mjj      = book(wk, name, "mjj",       "mjj",          100,      0,    5000   );
        h_ptjj     = book(wk, name, "ptjj",      "ptjj",         100,      0,     500   );
        h_etajj    = book(wk, name, "etajj",     "etajj",        100,     -4,       4   );
        h_phijj    = book(wk, name, "phijj",     "phijj",        100,     -3.2,     3.2 );

	//
	// ISR-jj hits
	//
        h_yStarISRjj  = book(wk, name, "yStarISRjj",   "yStarISRjj",      100,     -3,       3   );
        h_yBoostISRjj = book(wk, name, "yBoostISRjj",  "yBoostISRjj",     100,     -3,       3   );
        h_dPhiISRjj   = book(wk, name, "dPhiISRjj",    "dPhiISRjj",       100,     -3.2,     3.2 );
        h_dEtaISRjj   = book(wk, name, "dEtaISRjj",    "dEtaISRjj",       100,     -5,       5   );
        h_dRISRjj     = book(wk, name, "dRISRjj",      "dRISRjj",         100,      0,       5   );
        h_asymISRjj   = book(wk, name, "asymISRjj",    "asymISRjj",       100,     -0.1,     1   );
        h_chiISRjj    = book(wk, name, "chiISRjj",     "chiISRjj",        100,      0,      60   );
        h_mISRjj      = book(wk, name, "mISRjj",       "mISRjj",          100,      0,    5000   );
        h_ptISRjj     = book(wk, name, "ptISRjj",      "ptISRjj",         100,      0,     500   );
        h_etaISRjj    = book(wk, name, "etaISRjj",     "etaISRjj",        100,     -4,       4   );
        h_phiISRjj    = book(wk, name, "phiISRjj",     "phiISRjj",        100,     -3.2,     3.2 );

	//
	// ISR-Leadjet hists
	//
        h_dPhiLeadJetISR   = book(wk, name, "dPhiLeadJetISR",    "dPhiLeadJetISR",       100,     -3.2,     3.2 );
        h_dEtaLeadJetISR   = book(wk, name, "dEtaLeadJetISR",    "dEtaLeadJetISR",       100,     -5,       5   );
        h_dRLeadJetISR     = book(wk, name, "dRLeadJetISR",      "dRLeadJetISR",         100,      0,       5   );
        h_dRCloseJetISR    = book(wk, name, "dRCloseJetISR",     "dRCloseJetISR",        100,      0,       5   );


	h_jets    = new jetHists(name+"/jets",    wk, detailLevel);
	h_leadJet = new jetHists(name+"/leadJet", wk, detailLevel);
	h_sublJet = new jetHists(name+"/sublJet", wk, detailLevel);
	h_thrdJet = new jetHists(name+"/thrdJet", wk, detailLevel);

	h_gammas    = new gammaHists(name+"/gammas",    wk, detailLevel);
	h_leadGamma = new gammaHists(name+"/leadGamma", wk, detailLevel);

	m_etaCutVal   = etaCut;

      }

      TH1F* book(EL::Worker* wk, std::string name, std::string hname, std::string title, int nBins, float xmin, float xmax){
	TH1F* h_tmp = new TH1F((name+"/"+hname).c_str(),(hname+";"+title+";Entries").c_str(), nBins, xmin,   xmax);
	wk->addOutput(h_tmp);
	return h_tmp;
      }


      void Fill(const eventData& thisEvent, float jetPtCut = 50){

	float weight = thisEvent.weight;

	h_NPV     ->Fill(thisEvent.NPV, weight );
	h_mu_ave  ->Fill(thisEvent.mu_ave, weight );
	h_mu_act  ->Fill(thisEvent.mu_act, weight );
	
	jetData leadJet = thisEvent.jets.at(0);
	jetData sublJet = thisEvent.jets.at(1);

	TLorentzVector jet1 = leadJet.vec();
	TLorentzVector jet2 = sublJet.vec();

	float yStarjj  = ( jet1.Rapidity() - jet2.Rapidity() ) / 2.0;
	float yBoostjj = ( jet1.Rapidity() + jet2.Rapidity() ) / 2.0;

	h_yStarjj->Fill( yStarjj, weight);
	h_yBoostjj->Fill( yBoostjj, weight);
	h_dEtajj->Fill(jet1.Eta() - jet2.Eta(), weight);
	h_dPhijj->Fill(jet1.DeltaPhi(jet2), weight);
	h_dRjj  ->Fill(jet1.DeltaR(jet2), weight);
	h_chijj ->Fill(exp(2*fabs(yStarjj)), weight);
	double asymjj = (jet1.Pt()-jet2.Pt())/(jet1.Pt()+jet2.Pt());
	h_asymjj ->Fill(asymjj, weight);

	//
	//  jj kiniematics
	//
	TLorentzVector jjSystem = (jet1 + jet2);
	h_ptjj ->Fill(jjSystem.Pt(), weight);
	h_mjj  ->Fill(jjSystem.M(), weight);
	h_etajj  ->Fill(jjSystem.Eta(), weight);
	h_phijj  ->Fill(jjSystem.Phi(), weight);


	//
	//  ISR-jj kinematics
	//    (Assume photon for now...)
	TLorentzVector isrSystem = TLorentzVector();
	isrSystem.SetPtEtaPhiE(thisEvent.gammas.at(0).pt, 
			       thisEvent.gammas.at(0).eta, 
			       thisEvent.gammas.at(0).phi, 
			       thisEvent.gammas.at(0).E);

	float yStarISRjj  = ( jjSystem.Rapidity() - isrSystem.Rapidity() ) / 2.0;
	float yBoostISRjj = ( jjSystem.Rapidity() + isrSystem.Rapidity() ) / 2.0;
	h_yStarISRjj->Fill( yStarISRjj, weight);
	h_yBoostISRjj->Fill( yBoostISRjj, weight);
	h_dEtaISRjj->Fill(jjSystem.Eta() - isrSystem.Eta(), weight);
	h_dPhiISRjj->Fill(jjSystem.DeltaPhi(isrSystem), weight);
	h_dRISRjj  ->Fill(jjSystem.DeltaR(isrSystem), weight);
	h_chiISRjj ->Fill(exp(2*fabs(yStarISRjj)), weight);
	double asymISRjj = (jjSystem.Pt()-isrSystem.Pt())/(jjSystem.Pt()+isrSystem.Pt());
	h_asymISRjj ->Fill(asymISRjj, weight);

	TLorentzVector ISRjjSystem = (jjSystem + isrSystem);
	h_ptISRjj ->Fill(ISRjjSystem.Pt(), weight);
	h_mISRjj  ->Fill(ISRjjSystem.M(), weight);
	h_etaISRjj  ->Fill(ISRjjSystem.Eta(), weight);
	h_phiISRjj  ->Fill(ISRjjSystem.Phi(), weight);


	//
	//  LeadJet vs ISR 
	//
	h_dPhiLeadJetISR->Fill(isrSystem.DeltaPhi(jet1), weight);
	h_dEtaLeadJetISR->Fill(isrSystem.Eta() - jet1.Eta(), weight);
	h_dRLeadJetISR  ->Fill(isrSystem.DeltaR(jet1), weight);

	unsigned njets = thisEvent.jets.size();
	if(njets > 0)
	  h_leadJet -> Fill(leadJet, weight);
	if(njets > 1)
	  h_sublJet -> Fill(sublJet, weight);
	if(njets > 2)
	  h_thrdJet -> Fill(thisEvent.jets.at(2), weight);


	float dRClosestJetISR = 1e8;

	unsigned int njetsPTCut = 0;
	for(unsigned int i = 0;  i< njets; ++i){

	  float thisJetPt = thisEvent.jets.at(i).pt;
	  if(thisJetPt < jetPtCut) continue;

	  float thisJetAbsEta = fabs(thisEvent.jets.at(i).eta);
	  if(thisJetAbsEta > m_etaCutVal) continue;
	  ++njetsPTCut;

	  float thisDrISR = isrSystem.DeltaR(thisEvent.jets.at(i).vec());
	  if(thisDrISR  < dRClosestJetISR)
	    dRClosestJetISR = thisDrISR;

	  h_jets -> Fill(thisEvent.jets.at(i), weight);
	  
	}

	h_dRCloseJetISR  ->Fill(dRClosestJetISR, weight);
	
	h_nJet   ->Fill(njetsPTCut, weight  );
	h_nJet_l ->Fill(njetsPTCut, weight  );


	unsigned ngammas = thisEvent.gammas.size();
	if(ngammas > 0)
	  h_leadGamma -> Fill(thisEvent.gammas.at(0), weight);

	unsigned int ngammasPTCut = 0;
	for(unsigned int i = 0;  i< ngammas; ++i){

	  float thisJetAbsEta = fabs(thisEvent.gammas.at(i).eta);
	  if(thisJetAbsEta > 2.5) continue;
	  ++ngammasPTCut;

	  h_gammas -> Fill(thisEvent.gammas.at(i), weight);
	  
	}
	
	h_nGamma   ->Fill(ngammasPTCut, weight  );
	h_nGamma_l ->Fill(ngammasPTCut, weight  );


      }

    };


    eventHists*   hIncl; //!
    eventHists*   hInclYStar; //!
    eventHists*   hMjj200_400; //!
    eventHists*   hMjj200_400_Ystar; //!

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!

  // this is a standard constructor
  ProcessTLAMiniTree ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // these are the functions not inherited from Algorithm
  virtual EL::StatusCode configure ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(ProcessTLAMiniTree, 1);
};

#endif

#ifndef TLAAlgos_ProcessTLAMiniTree_H
#define TLAAlgos_ProcessTLAMiniTree_H

#include <EventLoop/StatusCode.h>
#include <EventLoop/Algorithm.h>
#include <EventLoop/Worker.h>

// rootcore includes
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

//algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

#include "TLAEventCleaning/TLALArEventVetoData.h"

// ROOT include(s):
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3D.h"
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

		std::string m_name;
		bool m_debug;

		// what am I running on
		bool m_is2015;
		bool m_doTruthOnly;
		bool m_isDijetNtupleTruth;
		bool m_isDijetNtupleDS;
		bool m_isDijetNtupleOffline;
		bool m_isTLANtupleTruth;
		bool m_isTLANtupleDS;
		bool m_isTLANtupleOffline;
		bool m_doData;

		// normalisation
		bool m_useCutflow;         // get normalisation from cutflow
		bool m_useWeighted;        // get normalisation from weighted cutflow (is m_useCutflow is True)
		float m_lumi;              // Lumi we are scaling to
		
		// GRL and pileup
		bool m_applyGRL;
		std::string m_GRLxml;
		bool m_doPileupFromMap;    // get pileup from map
		std::string m_pileupMap;
    
		// trigger selection
		bool m_doTrigger;          // this selects according to an OR of single jet triggers
		bool m_doTrigger_j110;     // this selects according to j110 only (overwrites above)
                bool m_useTriggerSF;       // apply scale factors? Default true
                std::string m_doTrigger_str;        // this requires the trigger given in the string (default "" does nothing)
                bool m_getTriggerFromMap;   // get the trigger decision base on the runNumber and lumiBlock from the pileup map
                bool m_getTriggerFromNTUP;  // for (most) DS it isn't there -> turn off or it will crash
                bool m_requireDStriggers;   // require one of the HLT_DS_perf chains, in selection and plotAllSRs

		// which jet collections to use
		std::string m_primaryJetInName;     // name in NTUP
		std::string m_primaryJetOutName;    // subdirectory name in histogram file
		bool m_doSecondaryJets;
		std::string m_secondaryJetInName;
		std::string m_secondaryJetOutName;
		
    		// Event cleaning configuration
		bool m_applyLArEventCleaning;        // veto events which fail LAr event cleaning
		bool m_invertLArEventCleaning;       // veto events which pass (only if above is true)
		bool m_applyTLALArEventVetoData;     // run tool
		std::string m_TLALArEventVetoFiles;
		TLALArEventVetoData * m_dataForLArEventVeto; //!

		// jet cleaning
		bool m_doCleaning;               // veto events which fail jet cleaning
		bool m_invertJetCleaning;        // veto events which pass (only if above is True)
		bool m_recalculateJetCleaning;   // don't take the value from the NTUP but recalculate based on saved variables. (will do this anyway if NTUP values != (0,1) )

		// event selection cuts
		float m_YStarCut;
		float m_YBoostCut;
		float m_etaCut;
		float m_leadJetPtCut;
		float m_subleadJetPtCut;

		// which hists to write
		bool m_plotCleaning;
                bool m_plotPtSlices;
		bool m_plotEtaSlices;
		bool m_plotMjjWindow;
		bool m_plotAllSRs;


		/*bool m_applySF;
		TString m_scaleFactorLocation;
		TString m_scaleFactorHistoName;*/
		// float cutValue;

		// counters etc
		float m_sampleEvents;   //! MC Events in the sample we are processing
		int m_eventCounter;     //!
		float m_mcEventWeight;  //!

	private:

		GoodRunsListSelectionTool*   m_grl;       //!

		Int_t     m_runNumber; //!
		Long64_t  m_eventNumber; //!
		int       m_lumiBlock; //!
		bool      m_LArError; //!
		uint32_t  m_LArFlags; //!
		uint32_t  m_timeStamp; //!
		uint32_t  m_timeStampNSOffset; //!
		int       m_NPV; //!

		float m_avgIntPerX; //!
		float m_avgIntPerX_fromMap; //!
		float m_avgIntPerX_fromAOD; //!

		float m_weight; //!
		float m_weight_xs; //!

		vector<float>* m_jet_pt; //!
		vector<float>* m_jet_eta; //!
		vector<float>* m_jet_phi; //!
		vector<float>* m_jet_E; //!
		vector<float>* m_jet_muonSegments; //!
		vector<float>* m_jet_EMFrac; //!
		vector<float>* m_jet_HECFrac; //!
		vector<float>* m_jet_timing; //!
		vector<float>* m_jet_negativeE; //!
		float m_MHT; //!
		float m_mjj; //!

		vector<int>*   m_jet_clean_passLooseBad; //!
		vector<int>*   m_jet_clean_passLooseBad_recalc;
		
		vector<float>* m_jet_LArQuality; //!
		vector<float>* m_jet_AverageLArQF; //!
		vector<float>* m_jet_HECQuality; //!
		vector<float>* m_jet_FracSamplingMax; //!
		vector<int>*   m_jet_FracSamplingMaxIndex; //!
		vector<float>* m_jet_LeadingClusterPt; //!
		vector<float>* m_jet_LeadingClusterSecondLambda; //!
		vector<float>* m_jet_LeadingClusterCenterLambda; //!
		vector<float>* m_jet_LeadingClusterSecondR; //!

                // these used to have //! but I want to modify them with my map
		vector<string>* m_passedTriggers; 
		vector<float>* m_triggerPrescales; 

		vector<float>* m_secJet_pt; //!
		vector<float>* m_secJet_eta; //!
		vector<float>* m_secJet_phi; //!
		vector<float>* m_secJet_E; //!
		vector<float>* m_secJet_muonSegments; //!
		vector<float>* m_secJet_EMFrac; //!
		vector<float>* m_secJet_HECFrac; //!
		vector<float>* m_secJet_timing; //!
		vector<float>* m_secJet_negativeE; //!
		
		vector<int>*   m_secJet_clean_passLooseBad; //!
		vector<int>*   m_secJet_clean_passLooseBad_recalc; //!

		vector<float>* m_secJet_LArQuality; //!
		vector<float>* m_secJet_AverageLArQF; //!
		vector<float>* m_secJet_HECQuality; //!
		vector<float>* m_secJet_FracSamplingMax; //!
		vector<int>*   m_secJet_FracSamplingMaxIndex; //!
		vector<float>* m_secJet_LeadingClusterPt; //!
		vector<float>* m_secJet_LeadingClusterSecondLambda; //!
		vector<float>* m_secJet_LeadingClusterCenterLambda; //!
		vector<float>* m_secJet_LeadingClusterSecondR; //!


		// for scale factors
		/*TH2D* m_hcalibration;//!
		double m_pt_freeze;//!
		double m_eta_freeze;//!*/

		TH2D* m_h2_LArError; //!
		TH2D* m_h2_LArFlags_Tool; //!
		TH2D* m_h2_Clean_Evt_Jet; //!
		TH2D* m_h2_Clean_Evt_Jet_w; //!
		TH2D* m_h2_jetCleaning; //!
		TH2D* m_h2_avgIntPerX_map_AOD; //!
		TH2F* m_h2_pileupMap; //!
                // could make this a vector and generic
                TH2F* m_h2_J75Map; //!
                TH2F* m_h2_J100Map; //!
		
		TH1D* m_h_cutflow_primary; //!
		TH1D* m_h_cutflow_secondary; //!
		TH1D* m_h_cutflow_primary_w; //!
		TH1D* m_h_cutflow_secondary_w; //!


		//
		// Jet Data
		//
		struct jetData{
		  float pt;
		  float eta;
		  float phi;
		  float E;
		  float muonSegments;
		  float timing;
		  float negativeE;
		  float EMFrac;
		  float HECFrac;
		  int   clean_passLooseBad;

		  float LArQuality;
		  float AverageLArQF;
		  float HECQuality;
		  float FracSamplingMax;
		  int   FracSamplingMaxIndex;
		  float LeadingClusterPt;
		  float LeadingClusterSecondLambda;
		  float LeadingClusterCenterLambda;
		  float LeadingClusterSecondR;

		  
		  jetData(float m_pt, float m_eta, float m_phi, float m_E, float m_muonSegments=0, float m_timing = 0, float m_negativeE = 0, 
			  float m_EMFrac = 0, float m_HECFrac = 0, int m_clean_passLooseBad = -9, float m_LArQuality = -9, float m_AverageLArQF = -9, float m_HECQuality = -9, 
			  float m_FracSamplingMax = -9, int m_FracSamplingMaxIndex = -9, float m_LeadingClusterPt = -9, float m_LeadingClusterSecondLambda = -9, 
			  float m_LeadingClusterCenterLambda = -9, float m_LeadingClusterSecondR = -9 ){
		    pt     = m_pt;
		    eta    = m_eta;
		    phi    = m_phi;
		    E      = m_E;
		    muonSegments = m_muonSegments;
		    timing	 = m_timing;
		    negativeE	 = m_negativeE;
		    EMFrac	 = m_EMFrac;
		    HECFrac      = m_HECFrac;
		    clean_passLooseBad = m_clean_passLooseBad;

		    LArQuality		       = m_LArQuality;
		    AverageLArQF	       = m_AverageLArQF;
		    HECQuality		       = m_HECQuality;
		    FracSamplingMax	       = m_FracSamplingMax;
		    FracSamplingMaxIndex       = m_FracSamplingMaxIndex;
		    LeadingClusterPt	       = m_LeadingClusterPt;
		    LeadingClusterSecondLambda = m_LeadingClusterSecondLambda;
		    LeadingClusterCenterLambda = m_LeadingClusterCenterLambda;
		    LeadingClusterSecondR      = m_LeadingClusterSecondR;

		    
		  }

		  TLorentzVector vec() const{
		    TLorentzVector vec = TLorentzVector();
		    vec.SetPtEtaPhiE(pt,eta,phi,E);
		    return vec;
		  }

		};

		//
		// Event Data
		//

		struct eventData{
		  int             runNumber;
		  int             eventNumber;
		  vector<jetData> jets;
		  float           yStar;
		  float           MHT;
		  float           MHTvec;
		  float           weight;
		  float           prescaleWeight = 1;
		  float           avgIntPerX;
                  int             LArError_offline;
                  int             LArError_tool;
                  int             LArFlags;
		  
		  eventData(unsigned int m_runNumber,
			    unsigned int m_eventNumber,
			    vector<float>* m_jet_pt,
			    vector<float>* m_jet_eta,
			    vector<float>* m_jet_phi,
			    vector<float>* m_jet_E,
			    vector<float>* m_jet_muonSegments,
			    vector<float>* m_jet_EMFrac,
			    vector<float>* m_jet_HECFrac,
			    vector<float>* m_jet_timing,
			    vector<float>* m_jet_negativeE,
			    vector<int>*   m_jet_clean_passLooseBad,

			    vector<float>* m_jet_LArQuality,
			    vector<float>* m_jet_AverageLArQF,
			    vector<float>* m_jet_HECQuality,
			    vector<float>* m_jet_FracSamplingMax,
			    vector<int>*   m_jet_FracSamplingMaxIndex,
			    vector<float>* m_jet_LeadingClusterPt,
			    vector<float>* m_jet_LeadingClusterSecondLambda,
			    vector<float>* m_jet_LeadingClusterCenterLambda,
			    vector<float>* m_jet_LeadingClusterSecondR,

			    float& m_MHT,
			    float& m_weight,
			    float& m_prescaleWeight,
			    float& m_avgIntPerX,
                            int m_LArError_offline,
                            int m_LArError_tool,
                            int m_LArFlags){
		    
		    runNumber       = m_runNumber;
		    eventNumber     = m_eventNumber;
		    prescaleWeight  = m_prescaleWeight;
		    avgIntPerX      = m_avgIntPerX;
                    LArError_offline = m_LArError_offline;
                    LArError_tool   = m_LArError_tool;
                    LArFlags        = m_LArFlags;

		    /* MHT = m_MHT; // it's not filled in Dijet NTUPs */		    
		    // calculate MHT
		    MHT=0;
		    TLorentzVector HTvec;
		    for(unsigned int i =0; i < m_jet_pt->size(); ++i){
		      if(m_jet_pt->at(i) < 60) continue;
		      MHT += m_jet_pt->at(i);
		      TLorentzVector tempJetTLV;
		      tempJetTLV.SetPtEtaPhiE(m_jet_pt->at(i), m_jet_eta->at(i), m_jet_phi->at(i), m_jet_E->at(i));
		      HTvec += tempJetTLV;
		    }		      
		    MHTvec = HTvec.Pt();

		    /* if(m_jet_clean_passLooseBad==NULL) { cout<< "whyyyyyy"<<endl;} */
		    /* cout<<"m_jet_clean_passLooseBad->size() = "<< m_jet_clean_passLooseBad->size() <<", m_jet_pt->size() = "<<m_jet_pt->size()<<endl; */

		    for(unsigned int i =0; i < m_jet_pt->size(); ++i){					       

		      /* cout << "filling jet " <<i<<endl; */

		      /* cout << "pt " << m_jet_pt->at(i) << endl; */
		      /* cout << "eta " << m_jet_eta->at(i) << endl; */
		      /* cout << "phi " << m_jet_phi->at(i) << endl; */
		      /* cout << "E " << m_jet_E->at(i) << endl; */
		      /* cout << "mS " << m_jet_muonSegments->at(i) << endl; */
		      /* cout << "tim " << m_jet_timing->at(i) << endl; */
		      /* cout << "negE " << m_jet_negativeE->at(i) << endl; */
		      /* cout << "EMf " << m_jet_EMFrac->at(i) << endl; */
		      /* cout << "HECf " << m_jet_HECFrac->at(i) << endl; */
		      /* cout << "cplb " << m_jet_clean_passLooseBad->at(i) << endl; */
		      
		      /* float LArQuality_i = -1; */
		      /* if(! m_jet_LArQuality==NULL) LArQuality_i = m_jet_LArQuality->at(i); */
		      /* cout << "LARq " << LArQuality_i << endl; */
		      /* cout << "LARqf " << m_jet_AverageLArQF->at(i) << endl; */
		      /* cout << "HECq " << m_jet_HECQuality->at(i) << endl; */
		      /* cout << "fsm " << m_jet_FracSamplingMax->at(i) << endl; */
		      /* cout << "fsmi " << m_jet_FracSamplingMaxIndex->at(i) << endl; */
		      /* cout << "lcp " << m_jet_LeadingClusterPt->at(i) << endl; */
		      /* cout << "lcsl " << m_jet_LeadingClusterSecondLambda->at(i) << endl; */
		      /* cout << "lccl " << m_jet_LeadingClusterCenterLambda->at(i) << endl; */
		      /* cout << "lcsr " << m_jet_LeadingClusterSecondR->at(i) << endl; */

		      /* cout << "all good?" << endl; */
	      
		      jetData thisJet = jetData(m_jet_pt->at(i),
						m_jet_eta->at(i),
						m_jet_phi->at(i),
						m_jet_E->at(i),
						m_jet_muonSegments->at(i),
						m_jet_timing->at(i),
						m_jet_negativeE->at(i),
						m_jet_EMFrac->at(i),
						m_jet_HECFrac->at(i),
						m_jet_clean_passLooseBad->at(i),
						m_jet_LArQuality->at(i),
						m_jet_AverageLArQF->at(i),
						m_jet_HECQuality->at(i),
						m_jet_FracSamplingMax->at(i),
						m_jet_FracSamplingMaxIndex->at(i),
						m_jet_LeadingClusterPt->at(i),
						m_jet_LeadingClusterSecondLambda->at(i),
						m_jet_LeadingClusterCenterLambda->at(i),
						m_jet_LeadingClusterSecondR->at(i)
						);
		      /* cout << "filled thisJet" << endl; */
		      jets.push_back(thisJet);
		    }
		    
		    weight   = m_weight * m_prescaleWeight;
		  }
		  
		  void dump() const{
		    cout << "EventDump: run: " << runNumber << " event: " << eventNumber << endl;
		    cout << "Prescale weight for this event:" << prescaleWeight << endl;
		    
		    for(auto& jet : jets){
		      cout << "jet:  pt:" << jet.pt << " eta: " << jet.eta << " phi: " << jet.phi << " E:"<< jet.E <<endl;
		    }
		  }
		  
		};


		struct eventHists{

		  //
		  //
		  //
		  
		  TH1D*     h_mjj;
		  TH1D*     h_mjj_uniform;
		  TH1D*     h_yStar;
		  TH1D*     h_MHT;
		  TH1D*     h_MHTvec;
		  TH1D*     h_pt_lead;
		  TH1D*     h_pt_sublead;
		  TH1D*     h_eta_lead;
		  TH1D*     h_eta_sublead;
		  TH1D*     h_phi_lead;
		  TH1D*     h_phi_sublead;

		  TH2F*     h_pt_eta_lead;
		  TH2F*     h_pt_eta_sublead;
		  TH3D*     h_pt_eta_muon;

		  TH1D*     h_EMFrac_lead;
		  TH1D*     h_EMFrac_sublead;
		  TH1D*     h_HECFrac_lead;
		  TH1D*     h_HECFrac_sublead;

		  TH1D*     h_LArQuality_lead;
		  TH1D*     h_AverageLArQF_lead;
		  TH1D*     h_HECQuality_lead;
		  TH2F*     h2_FracSamplingMax_FracSamplingMaxIndex_lead;
		  TH2F*     h2_EMfrac_phi_lead;
		  TH1D*     h_LeadingClusterPt_lead;
		  TH1D*     h_LeadingClusterSecondLambda_lead;
		  TH1D*     h_LeadingClusterCenterLambda_lead;
		  TH1D*     h_LeadingClusterSecondR_lead;

		  TH1D*     h_LArQuality_sublead;
		  TH1D*     h_AverageLArQF_sublead;
		  TH1D*     h_HECQuality_sublead;
		  TH2F*     h2_FracSamplingMax_FracSamplingMaxIndex_sublead;
		  TH2F*     h2_EMfrac_phi_sublead;
		  TH1D*     h_LeadingClusterPt_sublead;
		  TH1D*     h_LeadingClusterSecondLambda_sublead;
		  TH1D*     h_LeadingClusterCenterLambda_sublead;
		  TH1D*     h_LeadingClusterSecondR_sublead;

		  TH1D*     h_timing_lead;
		  TH1D*     h_timing_sublead;
		  TH1D*     h_negativeE_lead;
		  TH1D*     h_negativeE_sublead;

		  TH1D*     h_avgIntPerX;

		  TH2F*     h2_eta_pt_lead;
		  TH2F*     h2_eta_pt_sublead;
		  TH2F*     h2_eta_phi_lead;
		  TH2F*     h2_eta_phi_sublead;
                  TH2F*     h2_eta_lead_mjj_uniform;
                  TH2F*     h2_eta_sublead_mjj_uniform;
                  TH2F*     h2_yStar_mjj_uniform;
		  TH2F*     h2_avgIntPerX_timing_lead;
		  TH2F*     h2_avgIntPerX_negativeE_lead;
		  TH2F*     h2_avgIntPerX_pt_lead;
		  TH2F*     h2_avgIntPerX_eta_lead;

		  TH2F*    h2_timing_eta_lead;
		  TH2F*    h2_timing_eta_sublead;
		  TH2F*    h2_timing_mjj;

		  TH2F*    h2_LArError;
		  TH2F*    h2_LArError_w;
		  TH2F*    h2_LArFlags_Tool;

		  eventHists(std::string name, EL::Worker* wk){
		    
		    //
		    // hists
		    //

		    //earlier bins:
		    //Float_t bins[] = {203, 216, 229, 243, 257, 272, 287, 303, 319, 335, 352, 369, 387, 405, 424, 443, 462, 482, 502, 523, 544, 566, 588, 611, 634, 657, 681, 705, 730, 755, 781, 807, 834, 861, 889, 917, 946, 976, 1006, 1037, 1068, 1100, 1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165, 6285, 6407, 6531, 6658, 6787, 6918, 7052, 7188, 7326, 7467, 7610, 7756, 7904, 8055, 8208, 8364, 8523, 8685, 8850, 9019, 9191, 9366, 9544, 9726, 9911, 10100, 10292, 10488, 10688, 10892, 11100, 11312, 11528, 11748, 11972, 12200, 12432, 12669, 12910, 13156};
		    
		    //Latest bins by Edgar, after JES calibration
		    Float_t mjjbins[] = { 310, 323, 336, 350, 364, 379, 394, 410, 426, 443, 460, 478, 496, 515, 534, 554, 574, 595, 617, 639, 662, 685, 709, 733, 758, 783, 809, 836, 863, 891, 919, 948, 978, 1008, 1039, 1070, 1102, 1135, 1168, 1202, 1236, 1271, 1307, 1343, 1380, 1418, 1456, 1495, 1535, 1575, 1616, 1658, 1700, 1743, 1787, 1832, 1877, 1923, 1970, 2018, 2067, 2116, 2166, 2217, 2269, 2322, 2376, 2431, 2487, 2544, 2602, 2661, 2721, 2782, 2844, 2907, 2971, 3036, 3102, 3169, 3238, 3308, 3379, 3451, 3524, 3599, 3675, 3752, 3830, 3910, 3991, 4073, 4157, 4242, 4329, 4417, 4507, 4598, 4691, 4785, 4881, 4978, 5077, 5178, 5281, 5385, 5491, 5599, 5709, 5821, 5935, 6051, 6169, 6289, 6411, 6536, 6663, 6792, 6923, 7057, 7193, 7331, 7472, 7615, 7761, 7909, 8060, 8214, 8371, 8531, 8694, 8860, 9029, 9201, 9376, 9555, 9737, 9923, 10112, 10305, 10501, 10701, 10905, 11113, 11325, 11541, 11761, 11985, 12213, 12446, 12683, 12925, 13171 };

		    
		    Float_t ptbins[] =  {15. ,20. ,25. ,35. ,45. ,55. ,70. ,85. ,100. ,116. ,134. ,152. ,185. ,194. ,216. ,240. ,264. ,290. ,318. ,346.,376.,408.,442.,478.,516.,556.,598.,642.,688.,736.,786.,838.,894.,952.,1012.,1076.,1162.,1250.,1310.,1420.,1530.,1750.,1992.,2250.,2500.,2850.,3200.,3600.,4000.,4600.};
                
		    Float_t muonsegmentbins[] =  {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500};
		    
		    Float_t etabins[] =  {0,0.8,1.2,1.6,2.1,2.8,3.1,4.9};
		    
		    /*## Get an array which defines the eta bin boundaries for high pT jet checks
		      def getJetEtaBinsFine():
		      binLowE = [-4.8+0.08*x for x in range(0,121)]
		      return binLowE*/
		    
		    
		    Int_t  mjjbinnum = sizeof(mjjbins)/sizeof(Float_t) - 1;
		    Int_t  ptbinnum = sizeof(ptbins)/sizeof(Float_t) - 1;
		    Int_t  etabinnum = sizeof(etabins)/sizeof(Float_t) - 1;
		    Int_t  muonbinnum = sizeof(muonsegmentbins)/sizeof(Float_t) - 1;
		    
                    // kinematics
		    h_mjj             = book_variable(wk, name, "mjj"            , "mjj"                , mjjbinnum, mjjbins    );
		    h_yStar           = book         (wk, name, "yStar"          , "yStar"              , 100      , -2,    2    );
		    h_mjj_uniform     = book         (wk, name, "mjj_uniform"    , "mjj"                , 10000    , 0,    10000   );
		    h_MHT             = book_variable(wk, name, "MHT"            , "MHT"                , ptbinnum , ptbins   );
		    h_MHTvec          = book_variable(wk, name, "MHTvec"         , "MHTvec"             , ptbinnum , ptbins   );
		    h_pt_lead         = book_variable(wk, name, "pt_lead"        , "lead jet pt"        , ptbinnum , ptbins   );
		    h_pt_sublead      = book_variable(wk, name, "pt_sublead"     , "sublead jet pt"     , ptbinnum , ptbins   );
		    h_eta_lead        = book         (wk, name, "eta_lead"       , "lead jet eta"       , 100      , -4.5, 4.5   );
		    h_eta_sublead     = book         (wk, name, "eta_sublead"    , "sublead jet eta"    , 100      , -4.5, 4.5   );
		    h_phi_lead        = book         (wk, name, "phi_lead"       , "lead jet phi"       , 100      , -3.14, 3.14   );
		    h_phi_sublead     = book         (wk, name, "phi_sublead"    , "sublead jet phi"    , 100      , -3.14, 3.14  );
                    
                    // kinematic 2D
		    h_pt_eta_lead     = book         (wk, name, "pt_eta_lead"    , "lead jet pt;lead jet #eta"        , ptbinnum , ptbins, etabinnum, etabins );
		    h_pt_eta_sublead  = book         (wk, name, "pt_eta_sublead" , "sublead jet pt;sublead jet #eta"        , ptbinnum , ptbins, etabinnum, etabins );
                    // 3D
		    h_pt_eta_muon     = book         (wk, name, "pt_eta_muon"    , "pt;eta;muon segments"    , ptbinnum , ptbins, etabinnum, etabins, muonbinnum, muonsegmentbins );

                    // cleaning 1D
		    h_EMFrac_lead     = book         (wk, name, "EMFrac_lead"    , "lead jet EMFrac"    , 200      , -2, 2    );
		    h_EMFrac_sublead  = book         (wk, name, "EMFrac_sublead" , "sublead jet EMFrac" , 200      , -2, 2    );
		    h_HECFrac_lead    = book         (wk, name, "HECFrac_lead"   , "lead jet HECFrac"   , 150      , 0, 1.5   );
		    h_HECFrac_sublead = book         (wk, name, "HECFrac_sublead", "sublead jet HECFrac", 150      , 0, 1.5   );
		    h_timing_lead       = book (wk, name, "timing_lead",       "lead jet timing",        1000, -100, 100 );
		    h_timing_sublead    = book (wk, name, "timing_sublead",    "sublead jet timing",     1000, -100, 100 );
		    h_negativeE_lead    = book (wk, name, "negativeE_lead",    "lead jet negativeE",     210, -200, 10 );
		    h_negativeE_sublead = book (wk, name, "negativeE_sublead", "sublead jet negativeE",  210, -200, 10 );

		    h_LArQuality_lead                            = book (wk, name, "LArQuality_lead", "Lead jet LAr Quality",  200, -2, 8 );
		    h_AverageLArQF_lead                          = book (wk, name, "AverageLArQF_lead", "Lead jet average LAr QF",  140, 0, 70000 );
		    h_HECQuality_lead                            = book (wk, name, "HECQuality_lead", "Lead jet HEC Quality",  400, -10, 10 );
		    h_LeadingClusterPt_lead                      = book (wk, name, "LeadingClusterPt_lead", "Lead jet leading cluster p_{T}",  100, 0, 5000000 );
		    h_LeadingClusterSecondLambda_lead            = book (wk, name, "LeadingClusterSecondLambda_lead", "Lead jet leading cluster SecondLambda",  100, 0, 5000000 );
		    h_LeadingClusterCenterLambda_lead            = book (wk, name, "LeadingClusterCenterLambda_lead", "Lead jet leading cluster CenterLambda",  100, 0, 20000 );
		    h_LeadingClusterSecondR_lead                 = book (wk, name, "LeadingClusterSecondR_lead", "Lead jet leading cluster SecondR",  100, 0, 6000000 );

		    h_LArQuality_sublead                            = book (wk, name, "LArQuality_sublead", "Sublead jet LAr Quality",  200, -2, 8 );
		    h_AverageLArQF_sublead                          = book (wk, name, "AverageLArQF_sublead", "Sublead jet average LAr QF",  140, 0, 70000 );
		    h_HECQuality_sublead                            = book (wk, name, "HECQuality_sublead", "Sublead jet HEC Quality",  400, -10, 10 );
		    h_LeadingClusterPt_sublead                      = book (wk, name, "LeadingClusterPt_sublead", "LeadingClusterPt_sublead;Sublead jet leading cluster p_{T}",  100, 0, 5000000 );
		    h_LeadingClusterSecondLambda_sublead            = book (wk, name, "LeadingClusterSecondLambda_sublead", "Sublead jet leading cluster SecondLambda",  100, 0, 5000000 );
		    h_LeadingClusterCenterLambda_sublead            = book (wk, name, "LeadingClusterCenterLambda_sublead", "Sublead jet leading cluster CenterLambda",  100, 0, 20000 );
		    h_LeadingClusterSecondR_sublead                 = book (wk, name, "LeadingClusterSecondR_subleaed", "Sublead jet leading cluster SecondR",  100, 0, 6000000 );
		    
		    h_avgIntPerX        = book (wk, name, "avgIntPerX", "avgIntPerX",  100, 0, 100 );

		    // cleaning 2D
		    h2_FracSamplingMax_FracSamplingMaxIndex_lead = book (wk, name, "FracSamplingMax_FracSamplingMaxIndex_lead", "Lead jet FracSamplingMax;Lead jet FSM Index", 100, -2, -2, 25,0,25 );
		    h2_EMfrac_phi_lead                           = book (wk, name, "EMfrac_phi_lead_lead", "Lead jet EM frac;Lead jet #phi", 110,-0.1,1, 100,-3.14,3.14 );

		    h2_FracSamplingMax_FracSamplingMaxIndex_sublead = book (wk, name, "FracSamplingMax_FracSamplingMaxIndex_sublead", "Sublead jet FracSamplingMax;Sublead jet FSM Index", 100, -2, -2, 25,0,25 );
		    h2_EMfrac_phi_sublead                           = book (wk, name, "EMfrac_phi_sublead", "Sublead jet EM frac;Sublead jet #phi", 110, -0.1, 1, 100, -3.14, 3.14 );

		    h2_eta_pt_lead     = book (wk, name, "eta_pt_lead", "lead #eta;lead pt", 100, -4.5, 4.5, 100, 0, 2000 );
		    h2_eta_pt_sublead  = book (wk, name, "eta_pt_sublead", "sublead #eta;sublead pt", 100, -4.5, 4.5, 100, 0, 2000 );
		    h2_eta_phi_lead    = book (wk, name, "eta_phi_lead", "lead #eta;lead #phi", 100, -4.5, 4.5, 100, -3.14, 3.14 );
		    h2_eta_phi_sublead = book (wk, name, "eta_phi_sublead", "sublead #eta;sublead #phi", 100, -4.5, 4.5, 100, -3.14, 3.14 );

		    h2_eta_lead_mjj_uniform    = book (wk, name, "eta_lead_mjj_uniform", "lead #eta;m_{jj}",       60, -3, 3, 3600, 400, 4000 );
		    h2_eta_sublead_mjj_uniform = book (wk, name, "eta_sublead_mjj_uniform", "sublead #eta;m_{jj}", 60, -3, 3, 3600, 400, 4000 );
		    h2_yStar_mjj_uniform       = book (wk, name, "yStar_mjj_uniform", "yStar;m_{jj}",              50, -1, 1, 3600, 400, 4000 );

		    h2_avgIntPerX_timing_lead    = book (wk, name, "avgIntPerX_timing_lead", "avgIntPerX;timing_lead", 100, 0, 100, 100, -100, 100 );
		    h2_avgIntPerX_negativeE_lead = book (wk, name, "avgIntPerX_negativeE_lead", "avgIntPerX;negativeE_lead", 100, 0, 100, 105, -100, 10 );

		    h2_avgIntPerX_pt_lead  = book (wk, name, "avgIntPerX_pt_lead", "avgIntPerX;lead pt", 100, 0, 100, 100, 0, 2000 );
		    h2_avgIntPerX_eta_lead = book (wk, name, "avgIntPerX_eta_lead", "avgIntPerX;lead #eta", 100, 0, 100, 100, -4.5, 4.5 );

		    h2_timing_eta_lead    = book (wk, name, "timing_eta_lead", "lead jet timing;lead jet #eta", 120, -60, 60, 160, -4, 4 );
		    h2_timing_eta_sublead = book (wk, name, "timing_eta_sublead", "sublead jet timing;sublead jet #eta", 120, -60, 60, 160, -4, 4 );

		    h2_timing_mjj = book (wk, name, "timing_mjj", "lead two jets timing;mjj", 120, -60, 60, 200, 0, 10000 );
                    
                    // cleaning validation
                    h2_LArError = book (wk, name, "h2_LArError", "Offline: isArError;Tool: isLArError", 2, 0, 2, 4, 0, 4 );
                    h2_LArError_w = book (wk, name, "h2_LArError_w", "Offline: isArError;Tool: isLArError", 2, 0, 2, 4, 0, 4 );
                    h2_LArFlags_Tool = book (wk, name, "h2_LArFlags_Tool", "LAr Flags;Tool: isLArError", 500, 0, 500, 4, 0, 4 );
                    /* m_h2_LArError_postSelection->GetYaxis()->SetBinLabel(1,"None"); */
                    /* m_h2_LArError_postSelection->GetYaxis()->SetBinLabel(2,"NoiseBurst"); */
                    /* m_h2_LArError_postSelection->GetYaxis()->SetBinLabel(3,"MiniNoiseBurst"); */
                    /* m_h2_LArError_postSelection->GetYaxis()->SetBinLabel(4,"Other/Both"); */

		  }
		  
		  TH1D* book(EL::Worker* wk, std::string name, std::string hname, std::string title, int nBins, float xmin, float xmax){
		    TH1D* h_tmp = new TH1D((name+"/"+hname).c_str(),(hname+";"+title+";Entries").c_str(), nBins, xmin,   xmax);
		    wk->addOutput(h_tmp);
		    return h_tmp;
		  }
		  
		  
		  TH1D* book_variable(EL::Worker* wk, std::string name, std::string hname, std::string title, int nBins, const Float_t *xBins){
		    TH1D* h_tmp = new TH1D((name+"/"+hname).c_str(),(hname+";"+title+";Entries").c_str(), nBins, xBins);
		    wk->addOutput(h_tmp);
		    return h_tmp;
		  }
		  
		  TH2F* book(EL::Worker* wk, std::string name, std::string hname, std::string title, int nBinsX, float xmin, float xmax, int nBinsY, float ymin, float ymax){
		    TH2F* h_tmp = new TH2F((name+"/"+hname).c_str(),(hname+";"+title+";Entries").c_str(), nBinsX, xmin, xmax, nBinsY, ymin, ymax);
		    wk->addOutput(h_tmp);
		    return h_tmp;
		  }

		  TH2F* book(EL::Worker* wk, std::string name, std::string hname, std::string title, int nBinsX, const Float_t *xBins, int nBinsY, const Float_t *yBins){
		    TH2F* h_tmp = new TH2F((name+"/"+hname).c_str(),(hname+";"+title+";Entries").c_str(), nBinsX, xBins, nBinsY, yBins);
		    wk->addOutput(h_tmp);
		    return h_tmp;
		  }

		  TH3D* book(EL::Worker* wk, std::string name, std::string hname, std::string title, int nBinsX, const Float_t *xBins, int nBinsY, const Float_t *yBins, int nBinsZ, const Float_t *zBins){
		    TH3D* h_tmp = new TH3D((name+"/"+hname).c_str(),(hname+";"+title+";Entries").c_str(), nBinsX, xBins, nBinsY, yBins, nBinsZ, zBins);
		    wk->addOutput(h_tmp);
		    return h_tmp;
		  }
		  
		  void Fill(const eventData& thisEvent){
		    
		    float weight = thisEvent.weight;
		    
		    jetData leadJet = thisEvent.jets.at(0);
		    jetData sublJet = thisEvent.jets.at(1);
		    float leadJetMuonSegments = leadJet.muonSegments;
		    float sublJetMuonSegments = sublJet.muonSegments;
		    //some cleaning quantities
		    float leadJetEMFrac = leadJet.EMFrac;
		    float sublJetEMFrac = sublJet.EMFrac;
		    float leadJetHECFrac = leadJet.HECFrac;
		    float sublJetHECFrac = sublJet.HECFrac;
		    float MHT = thisEvent.MHT;
		    float MHTvec = thisEvent.MHTvec;

		    TLorentzVector jet1 = leadJet.vec();
		    TLorentzVector jet2 = sublJet.vec();
		    
		    float yStar = ( jet1.Rapidity() - jet2.Rapidity() ) / 2.0;
		    
		    //
		    //  simple kinematics
		    //
		    
		    
		    TLorentzVector jjSystem = (jet1 + jet2);
		    h_mjj             ->Fill(jjSystem.M(), weight);
		    h_yStar           ->Fill(yStar, weight);
		    h_MHT             ->Fill(MHT, weight);
		    h_MHTvec          ->Fill(MHTvec, weight);
		    h_mjj_uniform     ->Fill(jjSystem.M(), weight);
		    h_pt_lead         ->Fill(jet1.Pt(), weight);
		    h_pt_sublead      ->Fill(jet2.Pt(), weight);
		    h_eta_lead         ->Fill(jet1.Eta(), weight);
		    h_eta_sublead      ->Fill(jet2.Eta(), weight);
		    h_phi_lead         ->Fill(jet1.Phi(), weight);
		    h_phi_sublead      ->Fill(jet2.Phi(), weight);
		    h_pt_eta_lead      ->Fill(jet1.Pt(), jet1.Eta(), weight);
		    h_pt_eta_sublead   ->Fill(jet2.Pt(), jet2.Eta(), weight);
		    h_pt_eta_muon      ->Fill(jet1.Pt(), jet1.Eta(), leadJetMuonSegments, weight);
		    h_pt_eta_muon      ->Fill(jet2.Pt(), jet2.Eta(), sublJetMuonSegments, weight);

		    h_EMFrac_lead         ->Fill(leadJetEMFrac, weight);
		    h_EMFrac_sublead      ->Fill(sublJetEMFrac, weight);
		    h_HECFrac_lead         ->Fill(leadJetHECFrac, weight);
		    h_HECFrac_sublead      ->Fill(sublJetHECFrac, weight);
		    
		    h_timing_lead       ->Fill(leadJet.timing, weight);
		    h_timing_sublead    ->Fill(sublJet.timing, weight);
		    h_negativeE_lead    ->Fill(leadJet.negativeE, weight);
		    h_negativeE_sublead ->Fill(sublJet.negativeE, weight);
		    
		    h_avgIntPerX        -> Fill(thisEvent.avgIntPerX, weight); 
		    h2_eta_pt_lead     ->Fill(jet1.Eta(), jet1.Pt(), weight);
		    h2_eta_pt_sublead     ->Fill(jet1.Eta(), jet2.Pt(), weight);
		    h2_eta_phi_lead     ->Fill(jet1.Eta(), jet1.Phi(), weight);
		    h2_eta_phi_sublead     ->Fill(jet2.Eta(), jet2.Phi(), weight);

		    h2_eta_lead_mjj_uniform    -> Fill(jet1.Eta(), jjSystem.M(), weight);
		    h2_eta_sublead_mjj_uniform -> Fill(jet2.Eta(), jjSystem.M(), weight);
		    h2_yStar_mjj_uniform       -> Fill(yStar, jjSystem.M(), weight);

		    h2_avgIntPerX_timing_lead -> Fill(thisEvent.avgIntPerX, leadJet.timing, weight); 
		    h2_avgIntPerX_negativeE_lead -> Fill(thisEvent.avgIntPerX, leadJet.negativeE, weight); 

		    h2_avgIntPerX_pt_lead -> Fill(thisEvent.avgIntPerX, jet1.Pt(), weight); 
		    h2_avgIntPerX_eta_lead -> Fill(thisEvent.avgIntPerX, jet1.Eta(), weight); 

		    h_LArQuality_lead                            -> Fill (leadJet.LArQuality, weight);
		    h_AverageLArQF_lead                          -> Fill (leadJet.AverageLArQF, weight);
		    h_HECQuality_lead                            -> Fill (leadJet.HECQuality, weight);
		    h2_FracSamplingMax_FracSamplingMaxIndex_lead  -> Fill (leadJet.FracSamplingMax, leadJet.FracSamplingMaxIndex, weight);
		    h2_EMfrac_phi_lead                            -> Fill (leadJet.EMFrac, leadJet.phi, weight);
		    h_LeadingClusterPt_lead                      -> Fill (leadJet.LeadingClusterPt, weight);
		    h_LeadingClusterSecondLambda_lead            -> Fill (leadJet.LeadingClusterSecondLambda, weight);
		    h_LeadingClusterCenterLambda_lead            -> Fill (leadJet.LeadingClusterCenterLambda, weight);
		    h_LeadingClusterSecondR_lead                 -> Fill (leadJet.LeadingClusterSecondR, weight);

		    h_LArQuality_sublead                            -> Fill (sublJet.LArQuality, weight);
		    h_AverageLArQF_sublead                          -> Fill (sublJet.AverageLArQF, weight);
		    h_HECQuality_sublead                            -> Fill (sublJet.HECQuality, weight);
		    h2_FracSamplingMax_FracSamplingMaxIndex_sublead  -> Fill (sublJet.FracSamplingMax, sublJet.FracSamplingMaxIndex, weight);
		    h2_EMfrac_phi_sublead                            -> Fill (sublJet.EMFrac, sublJet.phi, weight);
		    h_LeadingClusterPt_sublead                      -> Fill (sublJet.LeadingClusterPt, weight);
		    h_LeadingClusterSecondLambda_sublead            -> Fill (sublJet.LeadingClusterSecondLambda, weight);
		    h_LeadingClusterCenterLambda_sublead            -> Fill (sublJet.LeadingClusterCenterLambda, weight);
		    h_LeadingClusterSecondR_sublead                 -> Fill (sublJet.LeadingClusterSecondR, weight);

		    h2_timing_eta_lead                              -> Fill (leadJet.timing, leadJet.eta, weight);
		    h2_timing_eta_sublead                           -> Fill (sublJet.timing, sublJet.eta, weight);
		    h2_timing_mjj                                   -> Fill (leadJet.timing, jjSystem.M(), weight);
		    h2_timing_mjj                                   -> Fill (sublJet.timing, jjSystem.M(), weight);

                    h2_LArError                                     -> Fill (thisEvent.LArError_offline, thisEvent.LArError_tool);
                    h2_LArError_w                                   -> Fill (thisEvent.LArError_offline, thisEvent.LArError_tool, weight);
                    h2_LArFlags_Tool                                -> Fill (thisEvent.LArFlags, thisEvent.LArError_tool);
		  }

		};
		
		
		//		eventHists*   hOffline; //!
		eventHists*   hIncl; //!
		eventHists*   hClean_PassEvt; //!
		eventHists*   hClean_FailEvt; //!
		eventHists*   hClean_PassJet; //!
		eventHists*   hClean_FailJet; //!
		eventHists*   hClean_PassEvt_PassJet; //!
		eventHists*   hClean_PassEvt_FailJet; //!
		eventHists*   hClean_FailEvt_PassJet; //!
		eventHists*   hClean_FailEvt_FailJet; //!
                eventHists*   hClean_PassToolFailEvt; //!

		eventHists*   hCentral; //!
		eventHists*   hCrack; //!
		eventHists*   hEndcap; //!
		eventHists*   hIncl_mjjWindow; //!
		eventHists*   hCentral_mjjWindow; //!
		eventHists*   hCrack_mjjWindow; //!
		eventHists*   hEndcap_mjjWindow; //!
		eventHists*   hIncl_mjjWindow2015; //!
		eventHists*   hCentral_mjjWindow2015; //!
		eventHists*   hCrack_mjjWindow2015; //!
		eventHists*   hEndcap_mjjWindow2015; //!
		//        eventHists*   hTrigger; //!
		//        eventHists*   hTruth; //!
		eventHists*   hPt200; //!
		eventHists*   hPt210; //!
		eventHists*   hPt220; //!
		eventHists*   hPt230; //!
		eventHists*   hPt240; //!
		eventHists*   hPt250; //!
                eventHists*   hJ100_yStar03; //!
                eventHists*   hJ100_yStar06; //!
                eventHists*   hJ75_yStar03; //!
                eventHists*   hJ75_yStar06; //!


		eventHists*   hSecIncl; //!
		eventHists*   hSecClean_PassEvt; //!
		eventHists*   hSecClean_FailEvt; //!
		eventHists*   hSecClean_PassJet; //!
		eventHists*   hSecClean_FailJet; //!
		eventHists*   hSecClean_PassEvt_PassJet; //!
		eventHists*   hSecClean_PassEvt_FailJet; //!
		eventHists*   hSecClean_FailEvt_PassJet; //!
		eventHists*   hSecClean_FailEvt_FailJet; //!
                eventHists*   hSecClean_PassToolFailEvt; //!

		eventHists*   hSecCentral; //!
		eventHists*   hSecCrack; //!
		eventHists*   hSecEndcap; //!
		eventHists*   hSecIncl_mjjWindow; //!
		eventHists*   hSecCentral_mjjWindow; //!
		eventHists*   hSecCrack_mjjWindow; //!
		eventHists*   hSecEndcap_mjjWindow; //!
		eventHists*   hSecIncl_mjjWindow2015; //!
		eventHists*   hSecCentral_mjjWindow2015; //!
		eventHists*   hSecCrack_mjjWindow2015; //!
		eventHists*   hSecEndcap_mjjWindow2015; //!
		eventHists*   hSecPt200; //!
		eventHists*   hSecPt210; //!
		eventHists*   hSecPt220; //!
		eventHists*   hSecPt230; //!
		eventHists*   hSecPt240; //!
		eventHists*   hSecPt250; //!
                eventHists*   hSecJ100_yStar03; //!
                eventHists*   hSecJ100_yStar06; //!
                eventHists*   hSecJ75_yStar03; //!
                eventHists*   hSecJ75_yStar06; //!

		
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

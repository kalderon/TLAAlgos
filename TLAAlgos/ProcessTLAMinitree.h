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
        bool m_useCutflow;
		bool m_doTruthOnly;
		bool m_isDijetNtupleTruth;
        bool m_isTLANtupleTruth;
        bool m_isTLANtupleTrig;
        bool m_isTLANtupleOffline;
        float m_YStarCut;
		float m_YBoostCut;             
		bool m_doBlind;                
		bool m_doData;                 
		bool m_useWeighted;            
		bool m_doCleaning;        
		float m_etaCut;       
		float m_leadJetPtCut;          
		float m_subleadJetPtCut;          
		float m_lumi;                   // Lumi we are scaling to

		float m_sampleEvents;           //! MC Events in the sample we are processing


		// float cutValue;
		int m_eventCounter;     //!

		std::string m_name;
		float m_mcEventWeight;  //!

	private:

		GoodRunsListSelectionTool*   m_grl;       //!

		long int m_runNumber; //!
		int m_eventNumber; //!
		int m_lumiBlock; //!

		float m_weight; //!
		float m_weight_xs; //!

		vector<float>* m_jet_pt; //!
		vector<float>* m_jet_eta; //!
		vector<float>* m_jet_phi; //!
		vector<float>* m_jet_E; //!

		vector<int>*   m_jet_clean_passLooseBad; //!

		vector<string>* m_passedTriggers; //!

		//
		// Jet Data
		//
		struct jetData{
			float pt;
			float eta;
			float phi;
			float E;

			jetData(float m_pt, float m_eta, float m_phi, float m_E){ 

				pt     = m_pt;
				eta    = m_eta;
				phi    = m_phi;
				E      = m_E;

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
			float           weight;

			eventData(unsigned int m_runNumber, unsigned int m_eventNumber, 
					vector<float>* m_jet_pt, vector<float>* m_jet_eta, vector<float>* m_jet_phi, vector<float>* m_jet_E, float& m_weight){

				runNumber = m_runNumber;
				eventNumber = m_eventNumber;

				for(unsigned int i =0; i < m_jet_pt->size(); ++i){
					jetData thisJet = jetData(m_jet_pt->at(i), m_jet_eta->at(i), m_jet_phi->at(i), m_jet_E->at(i));
					jets.push_back(thisJet);
				}

				weight   = m_weight;
			}

			void dump() const{
				cout << "EventDump: run: " << runNumber << " event: " << eventNumber << endl;
				for(auto& jet : jets){
					cout << "jet:  pt:" << jet.pt << " eta: " << jet.eta << " phi: " << jet.phi << " E:"<< jet.E <<endl;
				}
			}

		};


		struct eventHists{

			//
			// jj
			//

			TH1F*     h_mjj;
			TH1F*     h_mjj_uniform;
			TH1F*     h_ptjj;
			TH1F*     h_etajj;
			TH1F*     h_phijj;

			eventHists(std::string name, EL::Worker* wk){

				//
				// jj hists
				//

				Float_t bins[] = {203, 216, 229, 243, 257, 272, 287, 303, 319, 335, 352, 369, 387, 405, 424, 443, 462, 482, 502, 523, 544, 566, 588, 611, 634, 657, 681, 705, 730, 755, 781, 807, 834, 861, 889, 917, 946, 976, 1006, 1037, 1068, 1100, 1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165, 6285, 6407, 6531, 6658, 6787, 6918, 7052, 7188, 7326, 7467, 7610, 7756, 7904, 8055, 8208, 8364, 8523, 8685, 8850, 9019, 9191, 9366, 9544, 9726, 9911, 10100, 10292, 10488, 10688, 10892, 11100, 11312, 11528, 11748, 11972, 12200, 12432, 12669, 12910, 13156};

				Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1; 

				h_mjj          = book_variable(wk, name, "mjj"        ,       "mjj"        ,        binnum,   bins         );
				h_mjj_uniform  = book(wk, name, "mjj_uniform",       "mjj_uniform",        100,      0,    5000   );
				h_ptjj         = book(wk, name, "ptjj"       ,       "ptjj"       ,        100,      0,     500   );
				h_etajj        = book(wk, name, "etajj"      ,       "etajj"      ,        100,     -4,       4   );
				h_phijj        = book(wk, name, "phijj"      ,       "phijj"      ,        100,     -3.2,     3.2 );


			}

			TH1F* book(EL::Worker* wk, std::string name, std::string hname, std::string title, int nBins, float xmin, float xmax){
				TH1F* h_tmp = new TH1F((name+"/"+hname).c_str(),(hname+";"+title+";Entries").c_str(), nBins, xmin,   xmax);
				wk->addOutput(h_tmp);
				return h_tmp;
			}


			TH1F* book_variable(EL::Worker* wk, std::string name, std::string hname, std::string title, int nBins, const Float_t *xBins){
				TH1F* h_tmp = new TH1F((name+"/"+hname).c_str(),(hname+";"+title+";Entries").c_str(), nBins, xBins);
				wk->addOutput(h_tmp);
				return h_tmp;
			}

			void Fill(const eventData& thisEvent){

				float weight = thisEvent.weight;

				jetData leadJet = thisEvent.jets.at(0);
				jetData sublJet = thisEvent.jets.at(1);

				TLorentzVector jet1 = leadJet.vec();
				TLorentzVector jet2 = sublJet.vec();


				//
				//  jj kiniematics
				//
				TLorentzVector jjSystem = (jet1 + jet2);
				h_mjj          ->Fill(jjSystem.M(), weight);
				h_mjj_uniform  ->Fill(jjSystem.M(), weight);
				h_ptjj         ->Fill(jjSystem.Pt(), weight);
				h_etajj        ->Fill(jjSystem.Eta(), weight);
				h_phijj        ->Fill(jjSystem.Phi(), weight);

			}

		};


		eventHists*   hIncl; //!

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

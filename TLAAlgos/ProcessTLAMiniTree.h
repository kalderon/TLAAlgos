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
        /*bool m_applySF;
        TString m_scaleFactorLocation;
        TString m_scaleFactorHistoName;*/
    
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
        vector<float>* m_jet_muonSegments; //!

		vector<int>*   m_jet_clean_passLooseBad; //!

		vector<string>* m_passedTriggers; //!
        vector<float>* m_triggerPrescales; //!

		// for scale factors
		/*TH2D* m_hcalibration;//!
		double m_pt_freeze;//!
		double m_eta_freeze;//!*/

		//
		// Jet Data
		//
		struct jetData{
			float pt;
			float eta;
			float phi;
			float E;
            float muonSegments;

			jetData(float m_pt, float m_eta, float m_phi, float m_E, float m_muonSegments=0){

				pt     = m_pt;
				eta    = m_eta;
				phi    = m_phi;
				E      = m_E;
                muonSegments      = m_muonSegments;

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
            float prescaleWeight = 1;

			eventData(unsigned int m_runNumber, unsigned int m_eventNumber, 
					vector<float>* m_jet_pt, vector<float>* m_jet_eta, vector<float>* m_jet_phi, vector<float>* m_jet_E, vector<float>* m_jet_muonSegments, float& m_weight, float& m_prescaleWeight){

				runNumber = m_runNumber;
				eventNumber = m_eventNumber;
                prescaleWeight = m_prescaleWeight;

				for(unsigned int i =0; i < m_jet_pt->size(); ++i){
					jetData thisJet = jetData(m_jet_pt->at(i), m_jet_eta->at(i), m_jet_phi->at(i), m_jet_E->at(i), m_jet_muonSegments->at(i));
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
            TH1D*     h_pt_lead;
            TH1D*     h_pt_sublead;
            TH1D*     h_eta_lead;
            TH1D*     h_eta_sublead;
            TH1D*     h_phi_lead;
            TH1D*     h_phi_sublead;
            TH3D*     h_pt_eta_muon;

			eventHists(std::string name, EL::Worker* wk){

				//
				// hists
				//

               //earlier bins:
               //Float_t bins[] = {203, 216, 229, 243, 257, 272, 287, 303, 319, 335, 352, 369, 387, 405, 424, 443, 462, 482, 502, 523, 544, 566, 588, 611, 634, 657, 681, 705, 730, 755, 781, 807, 834, 861, 889, 917, 946, 976, 1006, 1037, 1068, 1100, 1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165, 6285, 6407, 6531, 6658, 6787, 6918, 7052, 7188, 7326, 7467, 7610, 7756, 7904, 8055, 8208, 8364, 8523, 8685, 8850, 9019, 9191, 9366, 9544, 9726, 9911, 10100, 10292, 10488, 10688, 10892, 11100, 11312, 11528, 11748, 11972, 12200, 12432, 12669, 12910, 13156};

                //Latest bins by Edgar, after JES calibration
                Float_t mjjbins[] = { 310, 323, 336, 350, 364, 379, 394, 410, 426, 443, 460, 478, 496, 515, 534, 554, 574, 595, 617, 639, 662, 685, 709, 733, 758, 783, 809, 836, 863, 891, 919, 948, 978, 1008, 1039, 1070, 1102, 1135, 1168, 1202, 1236, 1271, 1307, 1343, 1380, 1418, 1456, 1495, 1535, 1575, 1616, 1658, 1700, 1743, 1787, 1832, 1877, 1923, 1970, 2018, 2067, 2116, 2166, 2217, 2269, 2322, 2376, 2431, 2487, 2544, 2602, 2661, 2721, 2782, 2844, 2907, 2971, 3036, 3102, 3169, 3238, 3308, 3379, 3451, 3524, 3599, 3675, 3752, 3830, 3910, 3991, 4073, 4157, 4242, 4329, 4417, 4507, 4598, 4691, 4785, 4881, 4978, 5077, 5178, 5281, 5385, 5491, 5599, 5709, 5821, 5935, 6051, 6169, 6289, 6411, 6536, 6663, 6792, 6923, 7057, 7193, 7331, 7472, 7615, 7761, 7909, 8060, 8214, 8371, 8531, 8694, 8860, 9029, 9201, 9376, 9555, 9737, 9923, 10112, 10305, 10501, 10701, 10905, 11113, 11325, 11541, 11761, 11985, 12213, 12446, 12683, 12925, 13171 };

                
                Float_t ptbins[] =  {15. ,20. ,25. ,35. ,45. ,55. ,70. ,85. ,100. ,116. ,134. ,152. ,172. ,194. ,216. ,240. ,264. ,290. ,318. ,346.,376.,408.,442.,478.,516.,556.,598.,642.,688.,736.,786.,838.,894.,952.,1012.,1076.,1162.,1250.,1310.,1420.,1530.,1750.,1992.,2250.,2500.,2850.,3200.,3600.,4000.,4600.};
                
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

				h_mjj           = book_variable(wk, name, "mjj"           ,       "mjj"            ,        mjjbinnum    ,      mjjbins    );
                h_yStar         = book         (wk, name, "yStar"         ,       "yStar"          ,        100          ,      -2,    2    );
				h_mjj_uniform   = book         (wk, name, "mjj_uniform"   ,       "mjj_uniform"    ,        100          ,      0,    5000   );
                h_pt_lead       = book_variable(wk, name, "pt_lead"       ,       "pt_lead"        ,        ptbinnum     ,      ptbins   );
                h_pt_sublead    = book_variable(wk, name, "pt_sublead"    ,       "pt_sublead"     ,        ptbinnum     ,      ptbins   );
                h_eta_lead      = book         (wk, name, "eta_lead"      ,       "eta_lead"       ,        100          ,      -4.5, 4.5   );
                h_eta_sublead   = book         (wk, name, "eta_sublead"   ,       "eta_sublead"    ,        100          ,      -4.5, 4.5   );
                h_phi_lead      = book         (wk, name, "phi_lead"      ,       "phi_lead"       ,        100          ,      -3.14, 3.14   );
                h_phi_sublead   = book         (wk, name, "phi_sublead"   ,       "phi_sublead"    ,        100          ,      -3.14, 3.14  );
                h_pt_eta_muon   = book         (wk, name, "pt_eta_muon"   ,       "pt_eta_muon"    ,        ptbinnum          ,      ptbins, etabinnum, etabins, muonbinnum, muonsegmentbins );


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

				TLorentzVector jet1 = leadJet.vec();
				TLorentzVector jet2 = sublJet.vec();
                
                float yStar = ( jet1.Rapidity() - jet2.Rapidity() ) / 2.0;

				//
				//  simple kinematics
				//
				TLorentzVector jjSystem = (jet1 + jet2);
				h_mjj             ->Fill(jjSystem.M(), weight);
                h_yStar           ->Fill(yStar, weight);
				h_mjj_uniform     ->Fill(jjSystem.M(), weight);
                h_pt_lead         ->Fill(jet1.Pt(), weight);
                h_pt_sublead      ->Fill(jet2.Pt(), weight);
                h_eta_lead         ->Fill(jet1.Eta(), weight);
                h_eta_sublead      ->Fill(jet2.Eta(), weight);
                h_phi_lead         ->Fill(jet1.Phi(), weight);
                h_phi_sublead      ->Fill(jet2.Phi(), weight);
                h_pt_eta_muon      ->Fill(jet1.Pt(), jet1.Eta(), leadJetMuonSegments, weight);
                h_pt_eta_muon      ->Fill(jet2.Pt(), jet2.Eta(), sublJetMuonSegments, weight);

			}

		};


//		eventHists*   hOffline; //!
        eventHists*   hIncl; //!
//        eventHists*   hTrigger; //!
//        eventHists*   hTruth; //!

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

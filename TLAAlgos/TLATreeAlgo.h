#ifndef xAODAnaHelpers_TLATreeAlgo_H
#define xAODAnaHelpers_TLATreeAlgo_H

#include <xAODAnaHelpers/TreeAlgo.h>
#include <sstream>
#include <fstream>



class TLATreeAlgo : public TreeAlgo
{

public:

  // this is a standard constructor
  TLATreeAlgo (std::string className = "TLATreeAlgo");                                           //!

  // these are the functions inherited from Algorithm
//  virtual EL::StatusCode setupJob (EL::Job& job);           //!
//  virtual EL::StatusCode fileExecute ();                    //!
  virtual EL::StatusCode treeInitialize ();                 //!
//  virtual EL::StatusCode changeInput (bool firstFile);      //!
  virtual EL::StatusCode initialize ();                     //!
  virtual EL::StatusCode execute ();                        //!
//  virtual EL::StatusCode postExecute ();                    //!
//  virtual EL::StatusCode finalize ();                       //!
//  virtual EL::StatusCode treeFinalize ();                   //!

  // these are the functions not inherited from Algorithm
  //virtual EL::StatusCode configure ();                      //!

  /// @cond
  // this is needed to distribute the algorithm to the workers
  ClassDef(TLATreeAlgo, 1);                                 //!
  /// @endcond
    
private :
    
  bool m_isMC;//!
  float m_xs; //!
  float m_filtEff; //!
  int m_eventCounter; //!
  int m_numAMIEvents; //!
  int m_mcChannelNumber; //!
  double m_mcEventWeight; //!
  std::stringstream m_ss; //!

  
  //these variables are there to keep track of jet collection index in event-level variables, in case of multiple jet collections
  //kept for later when event-level jet-based quantities are needed
  std::string m_jetName;//!
  std::string m_trigJetName;//!
  std::string m_truthJetName;//!
    
  //CD: not sure why this is here? This function is needed to grab lumi weights from txt files
#ifndef __CINT__
    EL::StatusCode getLumiWeights(const xAOD::EventInfo* eventInfo);
#endif // not __CINT__

};

#endif

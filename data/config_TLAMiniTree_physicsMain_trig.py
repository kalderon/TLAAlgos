import ROOT
from xAH_config import xAH_config

c = xAH_config()

GRL = "$ROOTCOREBIN/data/TLAAlgos/data15_13TeV.periodAllYear_DetStatus-v73-pro19-08_DQDefects-00-01-02_PHYS_StandardGRL_All_Good_25ns_tolerable_IBLSTANDBY-DISABLE.xml"

#
#  Data Config
#
if not args.is_MC:
    applyGRL           = True

#
#  MC Config
#
else:
    applyGRL           = False


doTruthOnly = False

#
# Process Ntuple
#

# Cuts: The first thing (probably most urgent) we're going to do is run on truth, so for the moment no cleaning + trigger requirements. 
# The cuts we do on truth are usually y*<0.6 and pTLead>150 GeV and pTSublead>85 GeV and then we start with mjj at a reasonably high value to avoid the turn-on.

# TODO mjj binning
# TODO compare with Cat
# TODO run multi files!

c.setalg("ProcessTLAMiniTree", { "m_name"                   : "TLAAlgo",
                                 "m_debug"                  : False,
                                 "m_doTrigger"              : False,
                                 "m_doTruthOnly"            : doTruthOnly,
                                 "m_applyTLALArEventVetoData" : True,
                                 "m_useCutflow"             : False,
                                 "m_isTLANtupleTrig"        : True,
                                 "m_useWeighted"            : False,
                                 "m_doCleaning"             : True, 
                                 "m_etaCut"                 : 2.8,
                                 "m_leadJetPtCut"           : 185,
                                 "m_subleadJetPtCut"        : 85,
                                 "m_doTrigger"              : True,
                                 "m_YStarCut"               : 0.6,
                                 "m_lumi"                   : 1.0, 
                                 "m_applyGRL"               : False,
                                 "m_doData"                 : (not args.is_MC),
                                 "m_GRLxml"                 : GRL,
                               } )


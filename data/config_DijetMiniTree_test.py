import ROOT
from xAH_config import xAH_config

c = xAH_config()

GRL = "$ROOTCOREBIN/data/TLAAlgos/data16_13TeV.periodAllYear_DetStatus-v83-pro20-14_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml" # circulated 19.10
PUmap = "$ROOTCOREBIN/data/TLAAlgos/PRW/pileup_map_None_297730-309759_OflLumi-13TeV-005.root" # using relevant ilumicalc file from this, and make_pileup_map.py

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
                                 "m_useCutflow"             : False,
                                 "m_isDijetNtuple"          : True,
                                 "m_useWeighted"            : True,
                                 "m_doCleaning"             : False,
                                 "m_etaCut"                 : 2.8, # only applied if not Truth only? FIXME
                                 # "m_leadJetPtCut"           : 150,
                                 "m_leadJetPtCut"           : 200,
                                 "m_subleadJetPtCut"        : 85,
                                 "m_YStarCut"               : 0.6,
                                 "m_lumi"                   : 1.0, 
                                 "m_applyGRL"               : applyGRL,
                                 "m_doData"                 : (not args.is_MC),
                                 "m_GRLxml"                 : GRL,
                                 "m_doSecondaryJets"        : True,
                                 "m_secondaryJetName"       : "uncalibJet",
                                 "m_doPileupFromMap"        : True,
                                 "m_pileupMap"              : PUmap,
                               } )


import ROOT
from xAH_config import xAH_config

c = xAH_config()

# circulated 19.10
# GRL = "$ROOTCOREBIN/data/TLAAlgos/data16_13TeV.periodAllYear_DetStatus-v83-pro20-14_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml"
# PUmap = "$ROOTCOREBIN/data/TLAAlgos/PRW/pileup_map_None_297730-309759_OflLumi-13TeV-005.root" # made with dijet-TLA/lumi_and_pileup/make_pileup_map.py

# final 2016 w/ toroid off runs included circulated 31.10 / 01.11
GRL = "$ROOTCOREBIN/data/TLAAlgos/data16_13TeV.periodAllYear_DetStatus-v83-pro20-15_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns_ignore_TOROID_STATUS.xml"
PUmap = "$ROOTCOREBIN/data/TLAAlgos/PRW/pileup_map_None_297730-311481_OflLumi-13TeV-005.root"

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

c.setalg("ProcessTLAMiniTree",
         { "m_name"                   : "TLAAlgo",
           "m_debug"                  : False,
           
           # what am I running on
           "m_doTruthOnly"            : doTruthOnly,
           "m_isDijetNtupleOffline"   : True,
           "m_isDijetNtupleTrig"      : False,
           "m_doData"                 : (not args.is_MC),
           
           # normalisation
           "m_useCutflow"             : False, # get normalisation from cutflow
           "m_useWeighted"            : True,  # get normalisation from weighted cutflow (is m_useCutflow is True)
           "m_lumi"                   : 1.0, 
           
           # GRL and pileup
           "m_applyGRL"               : applyGRL,
           "m_GRLxml"                 : GRL,
           "m_doPileupFromMap"        : True, # get pileup from map
           "m_pileupMap"              : PUmap,
           
           # trigger selection
           "m_doTrigger"              : True, # this selects according to an OR of single jet triggers
           "m_doTrigger_j110"         : False, # this selects according to j110 only (overwrites above)
           
           # which jet collections to use
           "m_primaryJetInName"       : "jet",
           "m_primaryJetOutName"      : "OfflineJets", # the name of the histogram file
           # "m_primaryJetInName"       : "trigJet",
           # "m_primaryJetOutName"      : "TriggerJets", # the name of the histogram file
           
           "m_doSecondaryJets"        : True,
           "m_secondaryJetInName"     : "trigJet",
           "m_secondaryJetOutName"    : "TriggerJets",
           # "m_secondaryJetInName"     : "jet",
           # "m_secondaryJetOutName"    : "OfflineJets", # the name of the histogram file
           
           
           # event cleaning
           "m_applyLArEventCleaning"    : True, # veto events which fail LAr event cleaning, from AOD if isDijetNtupleOffline, from tool if ....Trig
           "m_invertLArEventCleaning"   : False, # veto events which pass event cleaning (only if above is true)
           "m_applyTLALArEventVetoData" : False, # run tool
           "m_TLALArEventVetoFiles"   : "$ROOTCOREBIN/data/TLAEventCleaning/event-veto-data/",
           
           # jet cleaning
           "m_recalculateJetCleaning" : False, # don't take the value from the NTUP but recalculate based on saved variables. (will do this anyway if NTUP values != (0,1) )
           "m_doCleaning"             : True, # veto events which fail jet cleaning
           "m_invertJetCleaning"      : False, # veto events which pass jet cleaning (only if above is True)
           
           # event selection cuts
           "m_etaCut"                 : 2.8, # only applied if not Truth only? FIXME
           "m_leadJetPtCut"           : 200, # was 150, can cull bonus selections
           "m_subleadJetPtCut"        : 85,
           "m_YStarCut"               : 0.6,
           
           # which hists to print - TO DO
           # ptslices, mjj slices, eta slices, etc
           
           } )

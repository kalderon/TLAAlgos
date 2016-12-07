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

c.setalg("ProcessTLAMiniTree", { "m_name"                   : "TLAAlgo",
                                 "m_debug"                  : False,
                                 "m_doTrigger"              : False, # this selects according to an OR of single jet triggers
                                 "m_doTrigger_j110"         : False, # this selects according to j110 only
                                 "m_doTruthOnly"            : doTruthOnly,
                                 "m_isDijetNtupleOffline"   : True,
                                 "m_isDijetNtupleTrig"      : False,
                                 "m_useCutflow"             : False, # get normalisation from cutflow
                                 "m_useWeighted"            : True,  # get normalisation from weighted cutflow

                                 # jet cleaning
                                 "m_doCleaning"             : False, # veto events which fail jet cleaning
                                 "m_invertJetCleaning"      : True,  # veto events which pass jet cleaning (only if above is True)

                                 # event cleaning
                                 "m_applyLArEventCleaning"    : True, # veto events which fail LAr event cleaning, from AOD if isDijetNtupleOffline, from tool if ....Trig
                                 "m_invertLArEventCleaning"   : True, # veto events which pass event cleaning (only if above is true)
                                 "m_applyTLALArEventVetoData" : True, # run tool
                                 "m_TLALArEventVetoFiles"   : "$ROOTCOREBIN/data/TLAEventCleaning/event-veto-data/",

                                 "m_etaCut"                 : 2.8, # only applied if not Truth only? FIXME
                                 # "m_leadJetPtCut"           : 150,
                                 "m_leadJetPtCut"           : 200,
                                 "m_subleadJetPtCut"        : 85,
                                 "m_YStarCut"               : 0.6,
                                 "m_lumi"                   : 1.0, 
                                 "m_applyGRL"               : applyGRL,
                                 "m_doData"                 : (not args.is_MC),
                                 "m_GRLxml"                 : GRL,
                                 "m_doSecondaryJets"        : False,
                                 "m_secondaryJetName"       : "trigJet",
                                 "m_doPileupFromMap"        : True, # get pileup from map
                                 "m_pileupMap"              : PUmap,
                               } )


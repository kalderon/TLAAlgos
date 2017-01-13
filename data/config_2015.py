import ROOT
from xAH_config import xAH_config

c = xAH_config()


GRL = "$ROOTCOREBIN/data/TLAAlgos/data15_13TeV.periodAllYear_DetStatus-v73-pro19-08_DQDefects-00-01-02_PHYS_StandardGRL_All_Good_25ns_tolerable_IBLSTANDBY-DISABLE.xml"
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
           "m_is2015"                   : True,
           "m_doTruthOnly"            : doTruthOnly,
           "m_isDijetNtupleOffline"   : False,
           "m_isDijetNtupleTrig"      : True,
           "m_doData"                 : (not args.is_MC),
           
           # normalisation
           "m_useCutflow"             : False, # get normalisation from cutflow
           "m_useWeighted"            : True,  # get normalisation from weighted cutflow (is m_useCutflow is True)
           "m_lumi"                   : 1.0, 
           
           # GRL and pileup
           "m_applyGRL"               : applyGRL,
           "m_GRLxml"                 : GRL,
           "m_doPileupFromMap"        : False, # get pileup from map
           "m_pileupMap"              : "",
           
           # trigger selection
           "m_doTrigger"              : False, # this selects according to an OR of single jet triggers
           "m_doTrigger_j110"         : False, # this selects according to j110 only (overwrites above)
           
           # which jet collections to use
           "m_primaryJetInName"       : "trigJet",
           "m_primaryJetOutName"      : "TriggerJets", # subdirectory name in histogram file
           "m_doSecondaryJets"        : False,
           
           # event cleaning
           "m_applyLArEventCleaning"    : True, # veto events which fail LAr event cleaning, from AOD if isDijetNtupleOffline, from tool if ....Trig
           "m_invertLArEventCleaning"   : False, # veto events which pass event cleaning (only if above is true)
           "m_applyTLALArEventVetoData" : True, # run tool
           "m_TLALArEventVetoFiles"   : "$ROOTCOREBIN/data/TLAEventCleaning/event-veto-data-2015/",


           # jet cleaning
           "m_doCleaning"             : True, # veto events which fail jet cleaning
           "m_invertJetCleaning"      : False, # veto events which pass jet cleaning (only if above is True)
           "m_recalculateJetCleaning" : False, # don't take the value from the NTUP but recalculate based on saved variables. (will do this anyway if NTUP values != (0,1) )

           # event selection cuts
           "m_etaCut"                 : 2.8, # only applied if not Truth only? FIXME
           "m_leadJetPtCut"           : 220, # was 150, can cull bonus selections
           "m_subleadJetPtCut"        : 85,
           "m_YStarCut"               : 0.6,
           
           # which hists to write
           "m_plotCleaning"           : True,
           "m_plotPtSlices"           : False,
           "m_plotEtaSlices"          : False,
           "m_plotMjjWindow"          : False,

           } )


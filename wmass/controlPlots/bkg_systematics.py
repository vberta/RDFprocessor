# bkg_systematics = {
#      # "puWeight"  : ["puWeightUp", "puWeightDown"],
#      "Muon_Trigger_BCDEF_SF": ["Muon_Trigger_BCDEF_SFstatUp", "Muon_Trigger_BCDEF_SFstatDown", "Muon_Trigger_BCDEF_SFsystUp", "Muon_Trigger_BCDEF_SFsystDown"],
#      "Muon_ID_BCDEF_SF"     : ["Muon_ID_BCDEF_SFstatUp", "Muon_ID_BCDEF_SFstatDown", "Muon_ID_BCDEF_SFsystUp", "Muon_ID_BCDEF_SFsystDown"],
#      "Muon_ISO_BCDEF_SF"    : ["Muon_ISO_BCDEF_SFstatUp", "Muon_ISO_BCDEF_SFstatDown", "Muon_ISO_BCDEF_SFsystUp", "Muon_ISO_BCDEF_SFsystDown"],
#      "corrected" : ["correctedUp", "correctedDown"],
#      "nom"       : ["jerUp", "jerDown", "jesTotalUp", "jesTotalDown", "unclustEnUp","unclustEnDown"],
# }
# bkg_systematics = {
#     "nom"       : ["jerUp", "jerDown"],
# }
bkg_systematics = {
     "ISO" : ["ISOreco_mc_eigen0Up","ISOreco_mc_eigen0Down","ISOreco_mc_eigen1Up","ISOreco_mc_eigen1Down","ISOreco_mc_eigen2Up","ISOreco_mc_eigen2Down"],
     "LHEScaleWeight" : ["LHEScaleWeight_muR0p5_muF0p5", "LHEScaleWeight_muR0p5_muF1p0", "LHEScaleWeight_muR1p0_muF0p5","LHEScaleWeight_muR1p0_muF2p0","LHEScaleWeight_muR2p0_muF1p0","LHEScaleWeight_muR2p0_muF2p0"],
     "Trigger" : ["Triggertrigger_mc_eigen0Up","Triggertrigger_mc_eigen0Down","Triggertrigger_mc_eigen1Up","Triggertrigger_mc_eigen1Down","Triggertrigger_mc_eigen2Up","Triggertrigger_mc_eigen2Down"],
     "corrected" : ["correctedUp", "correctedDown"],
     # "nom"       : ["jerUp", "jerDown", "jesTotalUp", "jesTotalDown", "unclustEnUp","unclustEnDown"],
     "nom"       : ["jesTotalUp", "jesTotalDown", "unclustEnUp","unclustEnDown"],
     # "nominal" : [''] #NEEDED FOR PREPARING HISTOS FILES
}
# bkg_systematics = {}
 
#  
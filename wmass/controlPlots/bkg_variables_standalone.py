import math

bkg_variables_standalone = {
    'appliesTo' : ['bkg_Signal*', 'bkg_Sideband*', "Dimuon"],
    'prefix' : "bkgSel",
    'variables' :{
        'Muon_corrected_MET_nom_mt':   ('M_{T} (Rochester corr./smear MET)',  120, 0, 120, 'Muon_corrected_MET_nom_mt[Idx_mu1]',True),
        'MET_pt':  ('MET P_{T}',  120, 0, 120, 'MET_pt',True),
        'PV_npvsGood' :  ('Number of good primary vertices',  100, 0, 100, 'PV_npvsGood',True),
        'Muon_pfRelIso04_all_corrected_pt': ('muon pfAbsIso04',100, 0., 40,'Muon_pfRelIso04_all[Idx_mu1]*Muon_corrected_pt[Idx_mu1]',False),
        'Muon_pfRelIso04_all': ('muon pfRelIso04', 100, 0., 0.5, 'Muon_pfRelIso04_all[Idx_mu1]',False),
        'Muon_eta':            ('muon eta', 100, -2.5, 2.5, 'Muon_eta[Idx_mu1]',False),
        'Muon_corrected_pt':   ('muon p_{T} (Rochester corr.)',  100, 25, 65, 'Muon_corrected_pt[Idx_mu1]',False),

        
    },
    'D2variables':{
        'Muon_pfRelIso04_all_corrected_MET_nom_mt':   ('M_{T} (Rochester corr./smear MET) VS muon pfRelIso04', 60, 0, 120, 50, 0., 0.5,"Muon_corrected_MET_nom_mt", "Muon_pfRelIso04_all","Muon_eta"),      
        # 'Muon_pfRelIso04_all_corrected_pt_corrected_MET_nom_mt': ('M_{T} (Rochester corr./smear MET) VS muon pfAbsIso04',60, 0, 120, 40, 0., 40, "Muon_corrected_MET_nom_mt","Muon_pfRelIso04_all_corrected_pt","Muon_eta"),         
        # 'Muon_pfRelIso04_all_MET_pt':   ('MET p_{T} VS muon pfRelIso04', 60, 0, 120,50, 0., 0.5,"MET_pt","Muon_pfRelIso04_all", "Muon_eta"),        
        # 'Muon_pfRelIso04_all_corrected_pt_MET_pt':   ('MET p_{T} VS muon pfIso04',  60, 0, 120,40, 0., 40,"MET_pt", "Muon_pfRelIso04_all_corrected_pt","Muon_eta"),
    },
    'ClousureVariables' : {
        'Muon_corrected_MET_nom_mt_VS_pt_VS_eta' :  ("M_{T} (Rochester corr./smear MET)", 120, 0, 120, "Muon_corrected_MET_nom_mt","Muon_corrected_pt","Muon_eta"),
        # 'MET_pt_VS_pt_VS_eta' :     ('MET P_{T}', 120, 0, 120, "MET_pt","Muon_corrected_pt","Muon_eta"),
        # 'Muon_eta_VS_pt_VS_eta':    ('muon eta', 100, -2.5, 2.5, "Muon_eta","Muon_corrected_pt","Muon_eta"),
        'Muon_corrected_pt_VS_pt_VS_eta':   ('muon p_{T} (Rochester corr.)',100, 25, 65, "Muon_corrected_pt","Muon_corrected_pt","Muon_eta"),
        },
    'WptVariables' : {
            'Muon_eta':            ('muon eta', 100, -2.5, 2.5, 'Muon_eta[Idx_mu1]',False),
            'Muon_corrected_pt':   ('muon p_{T} (Rochester corr.)',  100, 25, 65, 'Muon_corrected_pt[Idx_mu1]',False),
            'RecoZ_Muon_mass' :  ("Reco Z mass [GeV]", 16, 50, 130,  "RecoZ_Muon_mass",100, 0, 200,'RecoZ_Muon_corrected_pt'),
            'RecoZ_Muon_corrected_pt' :  ("Reco Z p_{T} [GeV]", 100, 0, 200, "RecoZ_Muon_corrected_pt","Muon_corrected_pt","Muon_eta"),
            }
}

# bkg_variables = {
#     'Muon1': {
#         'appliesTo' : ['bkg_Signal*', 'bkg_Sideband*'],
#         'inputCollection' : 'Muon',
#         'newCollection': 'bkgSelMuon1',
#         'index': 'Idx_mu1',
#         'variables': {
#             'corrected_pt':   ('muon p_{T} (Rochester corr.)',  100, 25, 65),
#             'eta':            ('muon eta', 100, -2.5, 2.5),
#             'corrected_MET_nom_mt':   ('M_{T} (Rochester corr./smear MET)',  120, 0, 120),
#             'pfRelIso04_all': ('muon pfRelIso04', 400, 0., 0.5),
#             'dxy':            ('muon dxy', 100, -0.01, 0.01),
#             'dz':             ('muon dz', 100, -0.05, 0.05),
#             },
#         'newvariables':{
#             'pfRelIso04_all_corrected_pt': ('muon pfAbsIso04',800, 0., 40,'Muon_pfRelIso04_all*Muon_corrected_pt'),
#         },
#         'D2variables':{
#             'pfRelIso04_all_corrected_MET_nom_mt':   ('M_{T} (Rochester corr./smear MET) VS muon pfRelIso04', 400, 0., 0.5,120, 0, 120, "Muon_pfRelIso04_all[Idx_mu1]","Muon_corrected_MET_nom_mt[Idx_mu1]", "Muon_eta[Idx_mu1]"),
#             'pfRelIso04_all_corrected_pt_corrected_MET_nom_mt':   ('M_{T} (Rochester corr./smear MET) VS muon pfIso04', 800, 0., 40,120, 0, 120, "bkgSelMuon1_pfRelIso04_all_corrected_pt[Idx_mu1]","Muon_corrected_MET_nom_mt[Idx_mu1]", "Muon_eta[Idx_mu1]"),            
#             # 'pfRelIso04_all_corrected_pt_corrected_MET_nom_mt':   ('M_{T} (Rochester corr./smear MET) VS muon pfIso04', 800, 0., 40,120, 0, 120, "Muon_pfRelIso04_all_corrected_pt[Idx_mu1]","Muon_corrected_MET_nom_mt[Idx_mu1]", "Muon_eta[Idx_mu1]"), #INTEGRATED WITH CONTROL PLOTS IMPLEMENTATION
# 
#             # 'pfRelIso04_all_MET_pt':   ('MET p_{T} VS muon pfRelIso04', 400, 0., 0.5,120, 0, 120, "Muon_pfRelIso04_all[Idx_mu1]","MET_pt"),
#             # 'pfRelIso04_all_corrected_pt_MET_pt':   ('MET p_{T} VS muon pfIso04', 800, 0., 40, 120, 0, 120,"Muon_pfRelIso04_all_corrected_pt[Idx_mu1]","MET_pt"),
#         },
#     },
#     'MET' :{
#         'appliesTo' : ['bkg_Signal*','bkg_Sideband*'],
#         'inputCollection' : 'MET',
#         'newCollection': 'bkgSelMET',#out of controlplots only implementation
#         'variables': {
#             'pt':  ('MET P_{T}',  120, 0, 120),
#         },
#     },
#         
#     'PV' : {
#         'appliesTo' : ['Signal*','Sideband*'],
#         'inputCollection' : 'PV',
#         'newCollection': 'bkgSelPV',#out of controlplots only implementation
#         'variables': {
#             'npvsGood' :  ('Number of good primary vertices',  100, 0, 100),
#             },
#     }
# }

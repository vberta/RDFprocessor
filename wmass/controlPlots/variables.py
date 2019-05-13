import math

# respect nanoAOD structure: Collection_modifiers_variable
variables =  {
    'PV' : {
        'appliesTo' : ['Signal*','Sideband*','Dimuon'],
        'inputCollection' : 'PV',
        'variables': {
            'npvsGood' :  ('Number of good primary vertexes',  100, 0, 100),
            },
        },
    'RecoZ' : {
        'appliesTo' : ['Dimuon'],
        'inputCollection' : 'RecoZ',
         'variables': {
            'Muon_corrected_mass' :  ('Dimuon mass (Rochester corr.)',  100, 50, 150),
            'Muon_corrected_pt' :  ('Dimuon p_{T} (Rochester corr.)',  100, 00, 100),
            'Muon_corrected_y' :  ('Dimuon y (Rochester corr.)',  100, -3, 3),
            'Muon_corrected_CStheta' :  ('Dimuon CS cosTheta (Rochester corr.)',  100, -1, 1),
            'Muon_corrected_CSphi' :  ('Dimuon CS phi (Rochester corr.)',  100, -math.pi, math.pi),
            'Muon_corrected_MET_nom_uPar' :  ('uPar on Z (Rochester corr./smear MET)',  100, -100, 100),
            'Muon_corrected_MET_nom_uPer' :  ('uPer on Z (Rochester corr./smear MET)',  100, -100, 100),
            },
        },
    'Muon1': {
        'appliesTo' : ['Signal*','Sideband*','Dimuon'],
        'inputCollection' : 'Muon',
        'newCollection': 'SelectedMuon1',
        'index': 'Idx_mu1',
        'variables': {
            'corrected_pt':   ('muon p_{T} (Rochester corr.)',  100, 25, 65),
            'corrected_MET_uPar':   ('uPar (Rochester corr./smear MET)',  100, -100, 100),
            'corrected_MET_uPer':   ('uPer (Rochester corr./smear MET)',  100, -100, 100),
            'corrected_MET_nom_Wlikemt':   ('W-like M_{T} (Rochester corr./smear MET)',  100, 0, 200),
            'corrected_MET_nom_mt':   ('M_{T} (Rochester corr./smear MET)',  100, 0, 100),
            'pfRelIso04_all': ('muon pfRelIso04', 100, 0., 0.5),
            'eta':            ('muon eta', 100, -2.5, 2.5),
            'dxy':            ('muon dxy', 100, -0.01, 0.01),
            'dz':             ('muon dz', 100, -0.05, 0.05),
            'sip3d':          ('sip3D', 100, 0, 5),
            'corrected_MET_nom_hpt': ('recoil p_{T} (Rochester corr./smear MET)',  100, 0, 200),
            'corrected_MET_nom_Wlikehpt': ('W-like recoil p_{T} (Rochester corr./smear MET)',  100, 0, 200),

            },
    },
    'Muon2': {
        'appliesTo' : ['Dimuon'],
        'inputCollection' : 'Muon',
        'newCollection': 'SelectedMuon2',
        'index': 'Idx_mu2',
        'variables': {
            'corrected_pt':   ('muon p_{T} (Rochester corr.)',  100, 25, 65),
            'corrected_MET_nom_Wlikemt':   ('W-like M_{T} (Rochester corr./smear MET)',  100, 0, 200),
            'corrected_MET_uPar':   ('uPar (Rochester corr./smear MET)',  100, -100, 100),
            'corrected_MET_uPer':   ('uPer (Rochester corr./smear MET)',  100, -100, 100),
            'eta':            ('muon eta', 100, -2.5, 2.5),
            'dxy':            ('muon dxy', 100, -0.01, 0.01),
            'dz':             ('muon dz', 100, -0.05, 0.05),
            'pfRelIso04_all': ('muon pfRelIso04', 100, 0., 0.5),
            'sip3d':          ('sip3D', 100, 0, 5),
            'corrected_MET_nom_Wlikehpt': ('W-like recoil p_{T} (Rochester corr./smear MET)',  100, 0, 200),
            },
    },
}


bkg_variables = {
    'Muon1': {
        'appliesTo' : ['bkg_Signal*', 'bkg_Sideband*'],
        'inputCollection' : 'Muon',
        'newCollection': 'bkgSelMuon1',
        'index': 'Idx_mu1',
        'variables': {
            'corrected_pt':   ('muon p_{T} (Rochester corr.)',  100, 25, 65),
            'corrected_MET_nom_mt':   ('M_{T} (Rochester corr./smear MET)',  120, 0, 120),
            'pfRelIso04_all': ('muon pfRelIso04', 100, 0., 0.5),
            'pfRelIso04_all_TIMES_corrected_pt': ('muon pfIso04', 200, 0., 40),
            'dxy':            ('muon dxy', 100, -0.01, 0.01),
            'dz':             ('muon dz', 100, -0.05, 0.05),
            'pfRelIso04_all_VS_corrected_MET_nom_mt':   ('M_{T} (Rochester corr./smear MET) VS muon pfRelIso04',  120, 0, 120,100, 0., 0.5),
            'pfRelIso04_all_TIMES_corrected_pt_VS_corrected_MET_nom_mt':   ('M_{T} (Rochester corr./smear MET) VS muon pfIso04',  120, 0, 120, 200, 0., 40),
            'pfRelIso04_all_VS_MET_pt':   ('M_{T} (Rochester corr./smear MET)  VS muon pfRelIso04',  120, 0, 120, 100, 0., 0.5),
            'pfRelIso04_all_TIMES_corrected_pt_VS_MET_pt':   ('M_{T} (Rochester corr./smear MET) VS muon pfIso04',  120, 0, 120, 200, 0., 40),

            },
    },
    'MET' :{
        'appliesTo' : ['bkg_Signal*','bkg_Sideband*'],
        'inputCollection' : 'MET',
        'variables': {
            'pt':  ('MET P_{T} (smear MET)',  120, 0, 120),
        },
    }
}

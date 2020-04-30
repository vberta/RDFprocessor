# corrected = Rochester
muon = '_corrected'

v0 = False 
if(v0) : version = ''
else :version = 'BCDEF_'

# nom = MET w/ jet smearing
met = '_nom'

mz_over_mw='1.135'
mw_over_mz = '0.881'
WptWeight = '(1.131-0.0105*GenV_preFSR_qt*%s+0.000145*(GenV_preFSR_qt*GenV_preFSR_qt*%s*%s))*' % (mw_over_mz,mw_over_mz,mw_over_mz) 

bkg_selections = {
    'bkg_Signal' : {
        'MC' : {
            'cut': \
                'Vtype==0 && ' + \
                'HLT_SingleMu24 && '+ \
                ('Muon%s_pt[Idx_mu1]>25. && ' % muon) + \
                ('Muon%s_MET%s_mt[Idx_mu1]>0. && ' % (muon, met) ) + \
                'MET_filters==1 && ' + \
                'nVetoElectrons==0 && ' + \
                '1',
            'weight' : 'Generator_weight*' +\
                'puWeight*' + \
                ('Muon_Trigger_%sSF[Idx_mu1]*'% version) + \
                ('Muon_ID_%sSF[Idx_mu1]*'% version) + \
                ('Muon_ISO_%sSF[Idx_mu1]'% version),
            },
        'DATA' : {
            'cut': \
                'Vtype==0 && ' + \
                'HLT_SingleMu24 && ' + \
                ('Muon%s_pt[Idx_mu1]>25. && ' % muon) + \
                ('Muon%s_MET%s_mt[Idx_mu1]>0. && ' % (muon, met)) +\
                'MET_filters==1 && ' + \
                'nVetoElectrons==0 && ' + \
                '1',
            'weight' : '',
            },
        },
    'bkg_Sideband' : {
        'MC' : {
            'cut': \
                'Vtype==1 && ' +\
                'HLT_SingleMu24 && ' +\
                ('Muon%s_pt[Idx_mu1]>25. && ' % muon) + \
                ('Muon%s_MET%s_mt[Idx_mu1]>0. && ' % (muon, met)) + \
                'MET_filters==1 && ' +\
                'nVetoElectrons==0 && ' +\
                '1',
            'weight' : 'Generator_weight*' +\
                'puWeight*' + \
                ('Muon_Trigger_%sSF[Idx_mu1]*'% version) + \
                ('Muon_ID_%sSF[Idx_mu1]*'% version) + \
                ('Muon_ISO_%sSF[Idx_mu1]'% version),
            },
        'DATA' : {
            'cut': \
                'Vtype==1 && ' +\
                'HLT_SingleMu24 && '+ \
                ('Muon%s_pt[Idx_mu1]>25. && ' % muon) + \
                ('Muon%s_MET%s_mt[Idx_mu1]>0. && ' % (muon, met)) + \
                'MET_filters==1 && ' +\
                'nVetoElectrons==0 && ' + \
                '1',
            'weight' : '',
            },
        },
    'Dimuon' : {
        'MC' : {
            'cut': \
                'Vtype==2 && '+ \
                'HLT_SingleMu24 && '+ \
                ('Muon%s_pt[Idx_mu1]>25. && Muon%s_pt[Idx_mu2]>20. && ' % (muon, muon) ) + \
                'MET_filters==1 && '+ \
                'nVetoElectrons==0 && '+ \
                '1',
            'weight' : 'Generator_weight*' +\
                'puWeight*' +\
                'Muon_Trigger_BCDEF_SF[Idx_mu1]*Muon_ISO_BCDEF_SF[Idx_mu1]*Muon_ID_BCDEF_SF[Idx_mu1]*' +\
                'Muon_ISO_BCDEF_SF[Idx_mu2]*Muon_ID_BCDEF_SF[Idx_mu2]',
            },
        'DATA' : {
            'cut': \
                'Vtype==2 && '+ \
                'HLT_SingleMu24 && '+ \
                ('Muon%s_pt[Idx_mu1]>25. && Muon%s_pt[Idx_mu2]>20. && ' % (muon, muon))+ \
                'MET_filters==1 && '+ \
                'nVetoElectrons==0 && '+ \
                '1',
            'weight' : '1',
            },
        },
        
    # 'bkg_all' : {
    #     'MC' : {
    #         'cut': \
    #             'Vtype==1 || Vtype==0 ' +\
    #             'HLT_SingleMu24 && ' +\
    #             ('Muon%s_pt[Idx_mu1]>25. && ' % muon) + \
    #             ('Muon%s_MET%s_mt[Idx_mu1]>0. && ' % (muon, met)) + \
    #             'MET_filters==1 && ' +\
    #             'nVetoElectrons==0 && ' +\
    #             '1',
    #         'weight' : 'Generator_weight*' +\
    #             'puWeight*' + \
    #             ('Muon_Trigger_%sSF[Idx_mu1]*'% version) + \
    #             ('Muon_ID_%sSF[Idx_mu1]*'% version) + \
    #             ('Muon_ISO_%sSF[Idx_mu1]'% version),
    #         },
    #     'DATA' : {
    #         'cut': \
    #             'Vtype==1 || Vtype==0' +\
    #             'HLT_SingleMu24 && '+ \
    #             ('Muon%s_pt[Idx_mu1]>25. && ' % muon) + \
    #             ('Muon%s_MET%s_mt[Idx_mu1]>0. && ' % (muon, met)) + \
    #             'MET_filters==1 && ' +\
    #             'nVetoElectrons==0 && ' + \
    #             '1',
    #         'weight' : '',
    #         },
    #     },    

}

#marc
# ptBinning = [30,31,32,33,34,35,36,37,38,40,42,44,46,48,50,52,54,57,60,65]
# etaBinning = [-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4]

#fine final binning (for fast analisys)
# ptBinning = [25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65]
# ptBinning = [26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65] #lore ana
ptBinning = [26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55] #lore ana bis
etaBinning = [-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4]
mtBinning = [0,30,200]

looseCutDict = {
    'A' : [0,15],
    'B' : [15,30],
    'C' : [30,40],
}
    

# etaBinning = [0.0,0.1]


#SYSTEMATIC BNNINGS
# etaBinning = [-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4] #largeEtaBins analysis for systematics
# etaBinning = [-2.4,-1.8,-1.2,-0.6,0,0.6,1.2,1.8,2.4] #VeryLargeEtaBins analysis for systematics
# [26,28,30,32,34,36,38,40,42,44,47,50,53,56,59,62,65] #used in correlated Fit
# etaBinning = [-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.0,2.4]#LargeEtaPtBins
# ptBinning = [25,29,33,37,41,45,49,53,57,61,65] #LargeEtaPtBins

#TESTING ONES
# ptBinning = [25,26,27,28,29,30]
# etaBinning = [-2.4,-2.3,-2.2]
# etaBinning = [0.0,0.1,0.2]#,0.3,0.4,0.5]
# etaBinning = [-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0]
#ptstudy fine
# ptBinning = [30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65]
# etaBinning = [0.0,0.1,0.2]

#RAW
# ptBinning = [30,32,34,36,38,41,45,50,65]
# etaBinning = [0.0,0.2,0.4,0.8,1.2,1.6,2.0,2.4]

#ptstudy
# ptBinning = [30,31,32,33,34,35,36,37,38,40,42,44,46,48,50,52,54,57,60,65]
# etaBinning = [0.0,0.2,0.4]

# 'PV_npvsGood<35 && ' + \
# 'PV_npvsGood<35 && ' + \
# 'PV_npvsGood<35 && ' + \
# 'PV_npvsGood<35 && ' + \
# (WptWeight) +\
# (WptWeight) +\
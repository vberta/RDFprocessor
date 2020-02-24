import os
import copy
import math
import argparse
import sys
sys.path.append('../../framework')
import ROOT
from RDFtreeV2 import RDFtree

# from bkg_histos import *
from bkg_histos_standalone import *
from bkg_histos_fast import *
from plotter import *

from sampleParser import *
from bkg_selections import *
# from bkg_variables import *
from bkg_variables_standalone import *
# from systematics import *
from bkg_fakerateAnalyzer import *
from bkg_systematics import *
from bkg_WptReweighter import *


ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(True);

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

parser = argparse.ArgumentParser("")
parser.add_argument('-tag', '--tag', type=str, default="TEST",      help="")
parser.add_argument('-dataYear', '--dataYear',type=int, default=2016, help="")
parser.add_argument('-hadd', '--hadd',type=int, default=False, help="")
parser.add_argument('-plot', '--plot',type=int, default=False, help="")
parser.add_argument('-fakerate', '--fakerate',type=int, default=False, help="")
parser.add_argument('-fakerateVar', '--fakeVar',type=int, default=False, help="variation of tight and loose cut for fakerate")
parser.add_argument('-fakeAna', '--fakeAna',type=int, default=False, help="complete bkg anaysis, with systematics")
parser.add_argument('-fakeFinal', '--fakeFinal',type=int, default=False, help="comparison plots of systematic analysis")
parser.add_argument('-fakeAnaCorr', '--fakeAnaCorr',type=int, default=False, help="complete bkg anaysis, with systematics with correlated Fit")
parser.add_argument('-fakeSS', '--fakeSS',type=int, default=False, help="strategy systemtatics plots")
parser.add_argument('-fakeSyst', '--fakeSyst',type=int, default=False, help="input systemtatics plots")
parser.add_argument('-rdf', '--rdf',type=int, default=True, help="")
parser.add_argument('-pretend', '--pretend',type=int, default=False, help="")
parser.add_argument('-restrict', '--restrict',type=str, default="", help="")
parser.add_argument('-wpt', '--wpt',type=int, default=False, help="W pt rewight option studies (to evaluate weight)")
parser.add_argument('-wptRew', '--wptRew',type=int, default=False, help="W pt rewight option")
parser.add_argument('-clousure', '--clousure',type=int, default=False, help="produce clousure histos (RDF)")




args = parser.parse_args()
tag = args.tag
dataYear = args.dataYear
hadd = args.hadd
plot = args.plot
fakerate = args.fakerate
fakerateVar = args.fakeVar
fakerateAnalysis = args.fakeAna
fakerateFinalPlots = args.fakeFinal
fakerateCorrelatedFit = args.fakeAnaCorr
fakerateStraSyst = args.fakeSS
fakerateSyst = args.fakeSyst
rdf = args.rdf
pretend = args.pretend
wpt =args.wpt
wptRew =args.wptRew
clousureFlag = args.clousure
restrictDataset = [ x for x in args.restrict.split(',') ]

print "tag =", bcolors.OKGREEN, tag, bcolors.ENDC, \
    ", dataYear =", bcolors.OKGREEN, str(dataYear), bcolors.ENDC


# def filterVariables(variables={}, selection='Signal', verbose=False):
#     if verbose: print '>>', variables
#     new_variables = copy.deepcopy(variables)
#     delete_vars = []
#     for ivar,var in new_variables.iteritems():
#         match = False
#         appliesTo = var['appliesTo']
#         for applyTo in appliesTo:
#             if applyTo[-1]=="*":
#                 if applyTo[0:-1] in selection:
#                     match = True
#             else:
#                 if applyTo==selection:  match = True
#         if not match:
#             delete_vars.append(ivar)
#     for ivar in delete_vars:
#         del new_variables[ivar]
#     if verbose: print '<<', new_variables
#     return new_variables

def filterVariables(variables={}, selection='Signal', verbose=False):
    if verbose: print '>>', variables
    new_variables = copy.deepcopy(variables)
    delete_vars = []
    appliesTo = new_variables['appliesTo']
    for applyTo in appliesTo:
        if applyTo[-1]=="*":
            if applyTo[0:-1] in selection:
                match = True
        else:
            if applyTo==selection:  match = True
    if not match:
        new_variables.clear()
    return new_variables

def RDFprocess(outDir, inputFile, selections, sample, clousureFlag):

    outDir = outDir
    inputFile = inputFile
    myselections = selections
    sample = sample
    clousureFlag = clousureFlag

    outputFile = "%s.root" % (sample_key)
    QCDFlag = False
    WFlag = False
    if 'QCD' in sample_key :
        QCDFlag = clousureFlag
    if 'WJets' in sample_key or 'WToMuNu' in sample_key :
        WFlag = wptRew


    p = RDFtree(outputDir=outDir, outputFile = outputFile,inputFile=inputFile,pretend = pretend)

      # create branches
    # for subsel_key, subsel in sample['subsel'].iteritems():
    outputFiles.append("%s" % (sample_key))
    for sel_key, sel in myselections.iteritems():
            # if len(sample['subsel'])>1 and subsel_key=='none': continue
            myvariables = filterVariables(bkg_variables_standalone, sel_key)
            # print '\tBranching: subselection', bcolors.OKBLUE, subsel_key, bcolors.ENDC, 'with selection' , bcolors.OKBLUE, sel_key, bcolors.ENDC
            print '\tBranching selection' , bcolors.OKBLUE, sel_key, bcolors.ENDC
            # print '\tAdding variables for collections', bcolors.OKBLUE, myvariables.keys(), bcolors.ENDC
            print '\tAdding variables for collections', bcolors.OKBLUE, myvariables['variables'].keys(), myvariables['D2variables'].keys(), bcolors.ENDC

            myselection = copy.deepcopy(sel)
            # myselection[dataType]['cut'] += subsel if subsel_key!='none' else ''
            # subsel_str= subsel if subsel_key!='none' else ''

            # for sKind, sList in bkg_systematics.iteritems():
            #     for sName in sList :
                    # print "systematic kind", bcolors.OKBLUE, sKind, bcolors.ENDC, "name=",bcolors.OKBLUE, sName, bcolors.ENDC
            p.branch(nodeToStart='input',
                        nodeToEnd=sel_key,#+'_'+sKind+'_'+sName,
                        outputFile=outputFile,
                        # modules = [bkg_histos_standalone(selections=myselection, variables=myvariables, dataType=dataType, xsec=sample['xsec'], inputFile=inputFile,ptBins=ptBinning, etaBins=etaBinning, systDict=bkg_systematics,clousureFlag=QCDFlag, wpt=wpt, wptFunc=ROOT.wptFunc, wptRew=WFlag)])
                        modules = [bkg_histos_fast(selections=myselection, variables=myvariables, dataType=dataType, xsec=sample['xsec'], inputFile=inputFile,ptBins=ptBinning, etaBins=etaBinning, systDict=bkg_systematics,clousureFlag=QCDFlag, wpt=wpt, wptFunc=ROOT.wptFunc, wptRew=WFlag)])

    p.getOutput()


myselections = {}

if not wpt :
    for cut in ['bkg_Signal', 'bkg_Sideband']:
        myselections['%sPlus' % cut]  = copy.deepcopy(bkg_selections['%s' % cut])
        myselections['%sMinus' % cut] = copy.deepcopy(bkg_selections['%s' % cut])
        for d in ['MC','DATA']:
            myselections['%sPlus' % cut][d]['cut']    += ' && Muon_charge[Idx_mu1]>0'
            myselections['%sMinus' % cut][d]['cut']   += ' && Muon_charge[Idx_mu1]<0'
else :
    myselections['Dimuon'] = copy.deepcopy(bkg_selections['Dimuon'])


# inputDir = ('/scratch/bertacch/NanoAOD%s-%s/' % (str(dataYear), tag))
inputDir = ('/scratch/sroychow/NanoAOD%s-%s/' % (str(dataYear), tag))


# outDir =  'NanoAOD%s-%s_RAW_multivar/' % (str(dataYear), tag)
# outDir =  'NanoAOD%s-%s_RAW_clousure/' % (str(dataYear), tag)
# outDir =  'NanoAOD%s-%sRAW_ptStudy_fine_FIXED_highPU/' % (str(dataYear), tag)
# outDir =  'NanoAOD%s-%sRAW_ptStudy_fine_FIXED_lowPU/' % (str(dataYear), tag)
# outDir =  'NanoAOD%s-%sRAW_ptStudy_fine_FIXED_rewighted/' % (str(dataYear), tag)  #buona per fare i plot di Zpt
# outDir =  'NanoAOD%s-%sRAW_ptStudy_fine_FIXED_rewighted_2/' % (str(dataYear), tag)  #davvero ripesato
# outDir =  'NanoAOD%s-%sRAW_ptStudy_fine_FIXED_rewighted_2_ratiomass_noPUcut/' % (str(dataYear), tag)  #davvero ripesato con rapporto masse

#outDir =  'NanoAOD%s-%sRAW_ptStudy_fine_FIXED/' % (str(dataYear), tag)
# outDir =  'NanoAOD%s-%s/' % (str(dataYear), tag)
# outDir =  'NanoAOD%s-%s_ptStudy_clousureQCD/' % (str(dataYear), tag)

# outDir =  'NanoAOD%s-%s_ptStudy_TEST/' % (str(dataYear), tag)
# outDir =  'NanoAOD%s-%s_ptStudy_reWApplied/' % (str(dataYear), tag)
# outDir =  'NanoAOD%s-%s_ptStudy_reWApplied_check/' % (str(dataYear), tag)

# outDir =  'NanoAOD%s-%s_ptStudy_syst/' % (str(dataYear), tag)

# outDir =  'NanoAOD%s-%s_fast_syst/' % (str(dataYear), tag) #correlated fit in nom only
# outDir =  'NanoAOD%s-%s_fast_DoubleCorrelated/' % (str(dataYear), tag) #correlated fit in all syst
# outDir =  'NanoAOD%s-%s_fast_LargeEtaBins/' % (str(dataYear), tag) #deltaEra=0.3
# outDir =  'NanoAOD%s-%s_fast_VeryLargeEtaBins/' % (str(dataYear), tag) #deltaeta=0.6
# outDir =  'NanoAOD%s-%s_fast_LinearUncorr/' % (str(dataYear), tag) #deltaeta=0.6
# outDir =  'NanoAOD%s-%s_fast_LargeEtaPtBins/' % (str(dataYear), tag) #deltaeta=0.4,deltaPt=4
# outDir =  'NanoAOD%s-%s_fast_TEST_largeBins/' % (str(dataYear), tag) #deltaeta=0.4,deltaPt=4 for testing purpose (EW asymmetry)
# outDir =  'NanoAOD%s-%s_LoreHistos/' % (str(dataYear), tag) #new histos from lorenzo
# outDir =  'NanoAOD%s-%s_LoreHistos_MetFix/' % (str(dataYear), tag) #new histos from lorenzo with nom met.
# outDir =  'NanoAOD%s-%s_LoreHistos_AntiIso_MetFix/' % (str(dataYear), tag) #new histos from lorenzo with nom met and antiIso SF evaluated
# outDir =  'NanoAOD%s-%s_LoreHistos_trigger_AntiIso_MetFix/' % (str(dataYear), tag) #new histos from lorenzo with nom met and antiIso SF evaluated, trigger SF from Marc
# outDir =  'NanoAOD%s-%s_LoreHistos_fix3_trigger_AntiIso_MetFix/' % (str(dataYear), tag) # FINAL GOOD LORE HISTOS, W/O SYST
outDir =  'NanoAOD%s-%s_LoreHistos_syst/' % (str(dataYear), tag) #lore histos with systs


if not os.path.isdir(outDir): os.system('mkdir '+outDir)

outputFiles = []
# restricted = ['SingleMuon', 'TTJets', 'DY', 'WJets']
# restricted = ['DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root','DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root','SingleMuon_Run2016C','SingleMuon_Run2016D','SingleMuon_Run2016E','SingleMuon_Run2016F','SingleMuon_Run2016G','SingleMuon_Run2016H','SingleMuon_Run2016B_ver2']

parser = sampleParser(restrict = restrictDataset, tag=tag, inputDir=inputDir, production_file ='/scratch/sroychow/mcsamples_'+str(dataYear)+'-'+str(tag)+'.txt', exclude = ['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'])
# parser = sampleParser(restrict = restricted, tag=tag, inputDir=inputDir, production_file ='/scratch/sroychow/mcsamples_'+str(dataYear)+'-'+str(tag)+'.txt', exclude = ['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'])


samples_dict = parser.getSampleDict()

for sample_key, sample in samples_dict.iteritems():
    # for subsel_key, subsel in sample['subsel'].iteritems():
        outputFiles.append("%s" % (sample_key))

if rdf:
    print "CAMBIA PARAMETRI SELF DEL MONDULO!!!!!!!!!!!!!!!!!!!!!!!!!!"
    from multiprocessing import Process

    procs = []


    ROOT.gInterpreter.Declare('''TSpline3 *wptFunc; ''')
    wptFile = ROOT.TFile.Open(outDir+"/WptReweight.root")
    wptFuncpy = wptFile.Get("splineRatio")
    ROOT.wptFunc = wptFuncpy


    for sample_key, sample in samples_dict.iteritems():

        print 'doing multiprocessing'

        if not sample['multiprocessing']: continue

        dataType = 'MC' if 'Run' not in sample_key else 'DATA'

        print 'Analysing sample', bcolors.OKBLUE, sample_key, bcolors.ENDC
        print '\tdirectories =', bcolors.OKBLUE, sample['dir'] , bcolors.ENDC
        print '\txsec = '+'{:0.3f}'.format(sample['xsec'])+' pb', \
        ', (data type is', bcolors.OKBLUE, dataType, bcolors.ENDC, ')'
        print '\tsubselections =', bcolors.OKBLUE, sample['subsel'] , bcolors.ENDC

        inputFile = ROOT.std.vector("std::string")()
        for x in sample['dir']: inputFile.push_back(inputDir+x+"/tree*.root")




        p = Process(target=RDFprocess, args=(outDir, inputFile, myselections,sample, clousureFlag))

        p.start()

        procs.append(p)

    for p in procs:
        p.join()


    ROOT.ROOT.EnableImplicitMT(32)#48

    for sample_key, sample in samples_dict.iteritems():

        print 'doing multithreading'

        if sample['multiprocessing']: continue

        dataType = 'MC' if 'Run' not in sample_key else 'DATA'

        print 'Analysing sample', bcolors.OKBLUE, sample_key, bcolors.ENDC
        print '\tdirectories =', bcolors.OKBLUE, sample['dir'] , bcolors.ENDC
        print '\txsec = '+'{:0.3f}'.format(sample['xsec'])+' pb', \
        ', (data type is', bcolors.OKBLUE, dataType, bcolors.ENDC, ')'
        print '\tsubselections =', bcolors.OKBLUE, sample['subsel'] , bcolors.ENDC

        inputFile = ROOT.std.vector("std::string")()
        for x in sample['dir']: inputFile.push_back(inputDir+x+"/tree*.root")

        RDFprocess(outDir, inputFile, myselections, sample,clousureFlag)



samples_merging = {
    'WToMuNu'  : [x for x in outputFiles if ('WJets' and 'WToMuNu') in x],###
    'WToETauNu'  : [x for x in outputFiles if ('WJets' and 'WToETauNu') in x],
    'DYJets' : [x for x in outputFiles if 'DYJetsToLL' in x],
    'QCD' : [x for x in outputFiles if 'QCD_' in x],
    'Top' : [x for x in outputFiles if ('TTJets' in x or 'ST_' in x)],
    'Diboson' : [x for x in outputFiles if ('WW_' in x or 'WZ_' in x or 'ZZ_' in x)],
    'Data' : [x for x in outputFiles if 'Run' in x],
}

print 'Samples to be merged:'
print bcolors.OKBLUE, samples_merging, bcolors.ENDC

outputMergedFiles = []
for sample_merging_key, sample_merging_value in samples_merging.iteritems():
        if len(sample_merging_value)>0:
            outputMergedFiles.append( '%s.root' % (sample_merging_key))
            cmd = 'hadd -f -k %s/%s.root' % (outDir,sample_merging_key)
            for isample in sample_merging_value:
                cmd += ' %s/%s.root' % (outDir,isample)
            if hadd:
                print bcolors.OKGREEN, cmd, bcolors.ENDC
                os.system(cmd)

print 'Final samples:'
print bcolors.OKBLUE, outputMergedFiles, bcolors.ENDC

def fakerate_analysis(systName='nom',systKind='nom', outdir=outDir+'/bkg',folder=outDir,norm = 35.922, tightCut = 0.15, looseCut=40, fitOnTemplate=True, ptBinning=ptBinning, etaBinning=etaBinning, onData=True, slicing=True,template=True,variations=fakerateVar, tightCutList=[0.02, 0.05, 0.10, 0.15, 0.2, 0.3, 0.5],looseCutList=[10,20,30,35,40,45,50,60,70,80], parabolaFit = False, EWSF='Iso_pt',correlatedFit=False,fastHistos=False, extraValidationPlots = False,statAna=False, nameSuff='') :
    # possible EWSF option: Fit, Fit_pt, Iso_pt, Mt, Mt_pt, 1, 0
    #if fastHistos=True also: Iso_global
    # if LoreHistos = True also: Mt_global
    outdirBkg = outdir+'/bkg_'+systName
    if not os.path.isdir(outdirBkg): os.system('mkdir '+outdirBkg)
    # if not os.path.isdir(outdirBkg+'_plot'): os.system('mkdir '+outdirBkg+'_plot')
    
    fake = bkg_fakerateAnalyzer(outdir=outdirBkg, folder=folder,norm = norm, fitOnTemplate=fitOnTemplate, ptBinning=ptBinning, etaBinning=etaBinning, onData=onData, slicing=slicing,systName=systName,systKind=systKind, parabolaFit=parabolaFit, EWSF=EWSF,tightCut = tightCut, looseCut=looseCut,fastHistos=fastHistos, extraValidationPlots=extraValidationPlots,correlatedFit=correlatedFit,statAna=statAna,nameSuff=nameSuff)#, tightCut = 5, varFake='pfRelIso04_all_corrected_pt_corrected_MET_nom_mt')
    # fake.integrated_preliminary()
    fake.differential_preliminary(template=template, correlatedFit=correlatedFit, produce_ext_output=True)

    if fakerateVar  :
        print "->starting variations..."
        for lcut in looseCutList :
                for tcut in tightCutList :
                    print "variation (l, t)",lcut, tcut
                    fake_i = bkg_fakerateAnalyzer(outdir=outdirBkg, folder=folder,norm = norm, fitOnTemplate=fitOnTemplate, ptBinning=ptBinning, etaBinning=etaBinning, onData=True, tightCut = tcut, looseCut=lcut, nameSuff="_"+str(lcut)+"_"+str(tcut),systName=systName, systKind=systKind, parabolaFit=parabolaFit,EWSF=EWSF)
                    fake_i.differential_preliminary(fakerate=True, produce_ext_output=True, template=template)

    # if not correlatedFit : #usa parametro "Final plots per tutto"
    if template :
        fake.fakerate_plots(variations=variations,tightCutList=tightCutList,looseCutList=looseCutList, parabolaFit = parabolaFit)


if fakerate :
    print "----> FakeRate analyzer:"
    selected = [s for s in outputMergedFiles]
    if not rdf : #if rdf=false --->expected already created the rootfiles with proper histos (done previously by rdf)
        skipHisto = True
    if not os.path.isdir(outDir+'/bkg'): os.system('mkdir '+outDir+'/bkg')
    # tightCutList = [0.02, 0.05, 0.10, 0.15, 0.2, 0.3, 0.5]
    tightCutList = [0.15]
    looseCutList = [10,20,30,35,40,45,50,60,70,80]

    #generalPars:
    TCUT = 0.15
    LCUT = 30 
    PARFIT = True
    NORM = 35.922
    FASTHISTOS=True
    EXTRAVALPLOT = False
    CORRFIT = False # this option is for nom only correlated fit
    STATANA = True
    CORRFITFINAL= True #correlated fit in the final plots
    EWKSF='1'#'Mt_global'
    SF_EVAL = False
    TEMPLATE = True
    NOM = ['nominal','']
    # NOM = ['nom','nom']

    if fakerateAnalysis :
        print "--> Not correlated analysis path..."
        fakerate_analysis(systName=NOM[1],systKind=NOM[0], norm=NORM, tightCutList=tightCutList,looseCutList=looseCutList, parabolaFit=PARFIT, looseCut=LCUT, tightCut = TCUT, slicing = True, fastHistos=FASTHISTOS, extraValidationPlots=EXTRAVALPLOT,  correlatedFit=False, EWSF=EWKSF, template = TEMPLATE)
        if SF_EVAL :
            fakerate_analysis(systName=NOM[1],systKind=NOM[0], norm=NORM, tightCutList=tightCutList,looseCutList=looseCutList, parabolaFit=PARFIT, looseCut=LCUT, tightCut = TCUT, slicing = True, fastHistos=FASTHISTOS, extraValidationPlots=EXTRAVALPLOT,  correlatedFit=False, EWSF=EWKSF, nameSuff='_noIsoSF_',template = TEMPLATE)
        for sKind, sList in bkg_systematics.iteritems():
            for sName in sList :
                print "-> systematatic:", sKind,sName
                if fakerateSyst :
                    fakerate_analysis(systKind=sKind,systName=sName,norm=NORM, tightCutList=tightCutList,looseCutList=looseCutList, parabolaFit=PARFIT, looseCut=LCUT, tightCut = TCUT,slicing = True,fastHistos=FASTHISTOS, extraValidationPlots=EXTRAVALPLOT,correlatedFit=False, EWSF=EWKSF, template = TEMPLATE)
    
    if fakerateCorrelatedFit :
        print "--> Correlated analysis path..."
        fakerate_analysis(systName=NOM[1],systKind=NOM[0], norm=NORM, tightCutList=tightCutList,looseCutList=looseCutList, parabolaFit=PARFIT, looseCut=LCUT, tightCut = TCUT, slicing = False, fastHistos=FASTHISTOS, extraValidationPlots=EXTRAVALPLOT, correlatedFit=True,EWSF=EWKSF)
        if STATANA :
            print "-> statistical uncertainity analysis..."
            fakerate_analysis(systName=NOM[1],systKind=NOM[0], norm=NORM, tightCutList=tightCutList,looseCutList=looseCutList, parabolaFit=PARFIT, looseCut=LCUT, tightCut = TCUT, slicing = False, fastHistos=FASTHISTOS, extraValidationPlots=EXTRAVALPLOT, correlatedFit=True, statAna=STATANA,EWSF=EWKSF)
        for sKind, sList in bkg_systematics.iteritems():
            for sName in sList :
                print "-> systematatic:", sKind,sName
                if fakerateSyst :
                    fakerate_analysis(systKind=sKind,systName=sName,norm=NORM, tightCutList=tightCutList,looseCutList=looseCutList, parabolaFit=PARFIT, looseCut=LCUT, tightCut = TCUT,slicing = False, fastHistos=FASTHISTOS, extraValidationPlots=EXTRAVALPLOT,correlatedFit=True, EWSF=EWKSF)
    
    if fakerateFinalPlots :
        print "--> Final plots..."
        fakeFinal = bkg_fakerateAnalyzer(outdir=outDir+'/bkg/bkg_'+NOM[1], folder=outDir,norm = NORM, fitOnTemplate=True, ptBinning=ptBinning, etaBinning=etaBinning, onData=True, slicing=False, systName=NOM[1],systKind=NOM[0], parabolaFit=PARFIT,EWSF=EWKSF, tightCut=TCUT, looseCut=LCUT, fastHistos=FASTHISTOS, extraValidationPlots=EXTRAVALPLOT, correlatedFit=CORRFITFINAL)#,  nameSuff='_noIsoSF_')
        fakeFinal.finalPlots(systDict=bkg_systematics, SymBands=True, outDir=outDir+'/bkg/',correlatedFit=CORRFIT, noratio=False, correlatedFit_alreadyDone=True, statAna=STATANA)

        
    #OLD FINAL PLOTS, fix folder sistem and remove me! 
    # if fakerateFinalPlots :
    #     print "-> Final plots..."
    #     fakeFinal = bkg_fakerateAnalyzer(outdir=outDir+'/bkg/bkg_nom', folder=outDir,norm = NORM, fitOnTemplate=True, ptBinning=ptBinning, etaBinning=etaBinning, onData=True, slicing=False, systName=NOM[1],systKind=NOM[0], parabolaFit=PARFIT,EWSF='Iso_pt', tightCut=TCUT, looseCut=LCUT, fastHistos=FASTHISTOS, extraValidationPlots=EXTRAVALPLOT, correlatedFit=CORRFIT)
    #     fakeFinal.differential_preliminary(fakerate=fakerate, correlatedFit=CORRFIT, produce_ext_output=True)
    #     fakeFinal.fakerate_plots(variations=fakerateVar,tightCutList=tightCutList,looseCutList=looseCutList, parabolaFit = PARFIT,correlatedFit=CORRFIT)
    #     fakeFinal.finalPlots(systDict=bkg_systematics, SymBands=True, outDir=outDir+'/bkg/',correlatedFit=CORRFIT, noratio=False, correlatedFit_alreadyDone=True)
        
        # fakeFinal.clousure_plots(outdir=outDir+'/bkg/', MC=False)

    if fakerateStraSyst :
        print "--> Strategy syst plots..."
        fakeFinal.strategy_syst(preOutDir=outDir)


if wpt :
    if not os.path.isdir(outDir+'/bkg'): os.system('mkdir '+outDir+'/bkg')
    wpt_rewighter= bkg_WptReweighter(outdir=outDir+'/bkg/', folder=outDir, norm = NORM, ptBinning=ptBinning, etaBinning=etaBinning, systName=NOM[1], systKind=NOM[0],interpol=True)
    wpt_rewighter.ratio_fitter()
            # if not os.path.isdir(outDir+'/bkg'): os.system('mkdir '+outDir+'/bkg_'+sName)
            # if not os.path.isdir(outDir+'/bkg/bkg_plot'): os.system('mkdir '+outDir+'/bkg/bkg_plot')
            # fake = bkg_fakerateAnalyzer(outdir=outDir+'/bkg', folder=outDir,norm = 35.922, fitOnTemplate=True, ptBinning=ptBinning, etaBinning=etaBinning, onData=True, slicing=False)#, tightCut = 5, varFake='pfRelIso04_all_corrected_pt_corrected_MET_nom_mt')
            # fake.integrated_preliminary()
            # fake.differential_preliminary(fakerate=True)
            #
            # print "starting variation"
            # if(fakerateVar) :
            #     for lcut in looseCutList :
            #             for tcut in tightCutList :
            #                 print "---variation (l, t)",lcut, tcut
            #                 fake_i = bkg_fakerateAnalyzer(outdir=outDir+'/bkg', folder=outDir,norm = 35.922, fitOnTemplate=False, ptBinning=ptBinning, etaBinning=etaBinning, onData=True, tightCut = tcut, looseCut=lcut, nameSuff="_"+str(lcut)+"_"+str(tcut))
            #                 fake_i.differential_preliminary(fakerate=True)
            # print "plotting"
            # fake.fakerate_plots(variations=fakerateVar,tightCutList=tightCutList,looseCutList=looseCutList )

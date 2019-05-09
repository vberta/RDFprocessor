import os
import copy
import math
import argparse
import sys
sys.path.append('../../framework')
import ROOT
from RDFtreeV2 import RDFtree

from controlPlots import *
from plotter import *

from sampleParser import *
from selections import *
from variables import *
from systematics import *

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
parser.add_argument('-rdf', '--rdf',type=int, default=True, help="")
parser.add_argument('-pretend', '--pretend',type=bool, default=False, help="")
args = parser.parse_args()
tag = args.tag
dataYear = args.dataYear
hadd = args.hadd
plot = args.plot
rdf = args.rdf
pretend = args.pretend
print "tag =", bcolors.OKGREEN, tag, bcolors.ENDC, \
    ", dataYear =", bcolors.OKGREEN, str(dataYear), bcolors.ENDC


def filterVariables(variables={}, selection='Signal', verbose=False):
    if verbose: print '>>', variables
    new_variables = copy.deepcopy(variables)
    delete_vars = []
    for ivar,var in new_variables.iteritems():
        match = False
        appliesTo = var['appliesTo']
        for applyTo in appliesTo:
            if applyTo[-1]=="*":
                if applyTo[0:-1] in selection:
                    match = True
            else:
                if applyTo==selection:  match = True
        if not match:
            delete_vars.append(ivar)
    for ivar in delete_vars:
        del new_variables[ivar]
    if verbose: print '<<', new_variables
    return new_variables

def RDFprocess(outDir, inputFile, selections, sample):

    outDir = outDir
    inputFile = inputFile
    myselections = selections
    sample = sample

    outputFile = "%s.root" % (sample_key)


    p = RDFtree(outputDir=outDir, outputFile = outputFile,inputFile=inputFile,pretend = pretend, syst = systematics)

      # create branches
    for subsel_key, subsel in sample['subsel'].iteritems():
        print "!!!!DEBUG!!!!!! IN subsel_key=", subsel_key, "subsel=", subsel
        outputFiles.append("%s%s" % (sample_key, ('_'+subsel_key if subsel_key!='none' else '')) )
        print "!!!!DEBUG!!!!!! outputfile in %s%s" % (sample_key, ('_'+subsel_key if subsel_key!='none' else ''))
        for sel_key, sel in myselections.iteritems():
            if len(sample['subsel'])>1 and subsel_key=='none': continue
            myvariables = filterVariables(variables, sel_key)
            print '\tBranching: subselection', bcolors.OKBLUE, subsel_key, bcolors.ENDC, 'with selection' , bcolors.OKBLUE, sel_key, bcolors.ENDC
            print '\tAdding variables for collections', bcolors.OKBLUE, myvariables.keys(), bcolors.ENDC

            myselection = copy.deepcopy(sel)
            myselection[dataType]['cut'] += subsel if subsel_key!='none' else ''
            subsel_str= subsel if subsel_key!='none' else ''
            p.branch(nodeToStart='input',
                        nodeToEnd='controlPlots'+sel_key+subsel_str,
                        outputFile=outputFile,
                        modules = [controlPlots(selections=myselection, variables=myvariables, dataType=dataType, xsec=sample['xsec'], inputFile=inputFile)])

    p.getOutput()


myselections = {}

for cut in ['Signal', 'Sideband', 'Dimuon']:
    if cut=='Dimuon':
        myselections['%s' % cut] = copy.deepcopy(selections['%s' % cut])
        continue
    myselections['%sPlus' % cut]  = copy.deepcopy(selections['%s' % cut])
    myselections['%sMinus' % cut] = copy.deepcopy(selections['%s' % cut])
    for d in ['MC','DATA']:
        myselections['%sPlus' % cut][d]['cut']    += ' && Muon_charge[Idx_mu1]>0'
        myselections['%sMinus' % cut][d]['cut']   += ' && Muon_charge[Idx_mu1]<0'


inputDir = ('/scratch/bertacch/NanoAOD%s-%s/' % (str(dataYear), tag))

outDir =  'NanoAOD%s-%s/' % (str(dataYear), tag)
if not os.path.isdir(outDir): os.system('mkdir '+outDir)

outputFiles = []

parser = sampleParser(restrict= ['QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8', 'QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'])
# parser = sampleParser(restrict=['WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'])
# parser = sampleParser()
samples_dict = parser.getSampleDict()

for sample_key, sample in samples_dict.iteritems():
    print "!!!!DEBUG!!!!!! OUT sample_key=", sample_key, "sample=", sample
    for subsel_key, subsel in sample['subsel'].iteritems():
        print "!!!!DEBUG!!!!!! OUT subsel_key=", subsel_key, "subsel=", subsel
        outputFiles.append("%s%s" % (sample_key, ('_'+subsel_key if subsel_key!='none' else '')) )
        print "!!!!DEBUG!!!!!! outputfile out %s%s" % (sample_key, ('_'+subsel_key if subsel_key!='none' else ''))


if rdf:

    from multiprocessing import Process

    procs = []

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

        p = Process(target=RDFprocess, args=(outDir, inputFile, myselections,sample))
        p.start()

        procs.append(p)

    for p in procs:
        p.join()


    # ROOT.ROOT.EnableImplicitMT(24)

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

        RDFprocess(outDir, inputFile, myselections, sample)



samples_merging = {
    'WToMuNu'  : [x for x in outputFiles if ('WJets' and 'WToMuNu') in x],
    'WToETauNu'  : [x for x in outputFiles if ('WJets' and 'WToETauNu') in x],
    'DYJets' : [x for x in outputFiles if 'DYJetsToLL' in x],
    'QCD' : [x for x in outputFiles if 'QCD_' in x],
    'Top' : [x for x in outputFiles if ('TTJets' in x or 'ST_' in x)],
    'Diboson' : [x for x in outputFiles if ('WW_' in x or 'WZ_' in x or 'ZZ_' in x)],
    'Data' : [x for x in outputFiles if 'Run' in x],
}

print 'Samples to be merged:'
print bcolors.OKBLUE, samples_merging, bcolors.ENDC

# outputMergedFiles = []
# for sel_key, sel in myselections.iteritems():
#     print '!!!!DEBUG!!!! ----' , sel_key, sel
#     for sample_merging_key, sample_merging in samples_merging.iteritems():
#         if len(sample_merging)>0:
#             outputMergedFiles.append( '%s_%s.root' % (sample_merging_key,sel_key))
#             # print '!!!!DEBUG!!!! ----File name:  %s_%s.root' % (sample_merging_key,sel_key)
#             cmd = 'hadd -f -k %s/%s_%s.root' % (outDir,sample_merging_key,sel_key)
#             # print '!!!!DEBUG!!!! ----initial File name: %s/%s_%s.root' % (outDir,sample_merging_key,sel_key)
#             for isample in sample_merging:
#
#                 cmd += ' %s/%s_%s.root' % (outDir,isample,sel_key)
#                 # cmd += ' %s/%s_%s.root' % (outDir,sample_merging_key,sel_key)
#                 # print '!!!!DEBUG!!!! ----isample', isample
#                 # print '!!!!DEBUG!!!! ----File dir:  %s/%s_%s.root'% (outDir,isample,sel_key)
#             if hadd:
#                 print bcolors.OKGREEN, cmd, bcolors.ENDC
#                 os.system(cmd)
#
# print 'Final samples:'
# print bcolors.OKBLUE, outputMergedFiles, bcolors.ENDC


outputMergedFiles = []
for sample_merging_key, sample_merging in samples_merging.iteritems():
        if len(sample_merging)>0:
            outputMergedFiles.append( '%s.root' % (sample_merging_key))
            cmd = 'hadd -f -k %s/%s.root' % (outDir,sample_merging_key)
            for isample in sample_merging:
                cmd += ' %s/%s.root' % (outDir,isample)
                print '!!!!DEBUG!!!! ----isample', isample

            if hadd:
                print bcolors.OKGREEN, cmd, bcolors.ENDC
                os.system(cmd)

print 'Final samples:'
print bcolors.OKBLUE, outputMergedFiles, bcolors.ENDC


# if plot:
#
#     for sel_key, sel in myselections.iteritems():
#
#         print sel_key
#         selected = [s for s in outputMergedFiles if sel_key in s]
#
#         plt = plotter(outdir=outDir+'/'+sel_key, folder=outDir, tag = sel_key, fileList=selected, norm = 35.922)
#         plt.plotStack()

if plot:

    for sel_key, sel in myselections.iteritems():

        print sel_key
        selected = [s for s in outputMergedFiles]

        plt = plotter(outdir=outDir+'/'+sel_key, folder=outDir, tag = sel_key, syst=systematics , fileList=selected, norm = 35.922)
        plt.plotStack()

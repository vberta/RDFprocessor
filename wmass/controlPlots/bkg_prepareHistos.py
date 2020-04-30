# import math
import sys
# sys.path.append('../../framework')
# from header import *
# from module import *
from bkg_systematics import *
# import numpy as np
# import ctypes
import ROOT
# ROOT.gROOT.SetBatch(True)
import os
import copy
from bkg_selections import *




class bkg_prepareHistos:
    def __init__(self, outDir, inputDir, ultraFast=True,extrap=True) :
        self.outDir = outDir
        self.inputDir = inputDir
        self.ultraFast = ultraFast
        self.extrap = extrap
        
        self.systDict = copy.deepcopy(bkg_systematics)
        self.systDict.update({'nominal':['']})
    
    def prepare(self) :
        fileList = ['TTbar', 'DYJets', 'DiBoson', 'ST', 'QCD','SingleMuon','WJets']
        # outDir = 'NanoAOD2016-V1MCFinal_LoreHistos_syst'
        varName = '__SelMuon1_eta_SelMuon1_corrected_pt_SelMuon1_charge__'
        regList = ["SIGNAL", "AISO","QCD","SIDEBAND"]
        
        if self.extrap :
            for lcut, lbin in looseCutDict.iteritems() :
                regList.append('QCD'+lcut)
                regList.append('SIDEBAND'+lcut)
                
        dirDict = {}

        for f in fileList : 
            inFile =  ROOT.TFile.Open("./"+self.inputDir+f+".root")
            outFile =  ROOT.TFile("./"+self.outDir+f+".root", "recreate")
            for r in regList :
                regExtrapFlag = self.isExtrapReg(r)
                for sKind, sList in self.systDict.iteritems():  
                    if sKind!='nominal' and regExtrapFlag : #extrap region only for nominal
                        continue                      
                    dirDict[f+r+sKind] =    outFile.mkdir(r+'_'+sKind)
                    for sName in sList : 
                        inFile.cd()
                        if ROOT.gDirectory.Get(r+'_'+sKind+'/'+r+varName+sName)==None : #this syst is not present
                            # if sKind=='LHEScaleWeight' and f=='WJets' :
                            print "no syst in:", f, r, sKind, sName
                            
                            if sKind == 'ISO_stat' : # if ISO_stat the trigger sist is used (these systs have the same magnitude)
                                isoKind = 'Trigger'
                                isoName = sName.replace('ISOreco','Triggertrigger')
                                if ROOT.gDirectory.Get(r+'_'+isoKind+'/'+r+varName+isoName)!=None : 
                                    print "iso stat built in:", f, r
                                    sNameNom = isoName
                                    sKindNom = isoKind
                                else : #if trigger is missing use the nominal
                                    sNameNom = ''
                                    sKindNom = 'nominal'
                            else : 
                                sNameNom = ''
                                sKindNom = 'nominal'
                            
                            h = inFile.Get(r+'_'+sKindNom+'/'+r+varName+sNameNom)
                            h.SetName(r+varName+sName)
                            
                        else :
                            h = inFile.Get(r+'_'+sKind+'/'+r+varName+sName)
                        dirDict[f+r+sKind].cd()
                        h.Write()
                        
    def h_add(self) :
        cmdList = []
        d = './'+self.outDir
        if self.ultraFast :
             cmdList.append('hadd '+d+'WToMuNu_ultra.root '+d+'TTbar.root '+d+'DYJets.root '+d+'DiBoson.root '+d+'ST.root '+d+'WJets.root')
             cmdList.append('cp '+d+'SingleMuon.root '+d+'Data_ultra.root')
        
        cmdList.append('hadd '+d+'EWKbkg.root '+d+'TTbar.root '+d+'DYJets.root '+d+'DiBoson.root '+d+'ST.root')
        cmdList.append('rm '+d+'TTbar.root '+d+'DYJets.root '+d+'DiBoson.root '+d+'ST.root')
        cmdList.append('mv '+d+'WJets.root '+d+'WToMuNu.root')
        cmdList.append('mv '+d+'SingleMuon.root '+d+'Data.root')        
        for i in cmdList :
            os.system(i)

    def isExtrapReg(self, reg) :
        FlagOut = False
        for lcut, lbin in looseCutDict.iteritems() :
                if reg == 'QCD'+lcut:
                    FlagOut = True
                if reg == 'SIDEBAND'+lcut:
                    FlagOut = True
        return  FlagOut 



#debugging
outDir =  'NanoAOD2016-V1MCFinal_LoreHistos_newSyst/'
preparator = bkg_prepareHistos(outDir=outDir+'/', inputDir=outDir+'/rawInput/',ultraFast = True,extrap=True)
preparator.prepare()
preparator.h_add()


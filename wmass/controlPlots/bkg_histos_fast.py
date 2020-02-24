import math
import ROOT
import re
import copy

import sys
sys.path.append('../../framework')
from module import *
from header import *
from array import array
import numpy as np
from bkg_selections import *

varDict_fast = {
    'Muon_eta':            'Muon_eta[Idx_mu1]',
    'Muon_corrected_pt':   'Muon_corrected_pt[Idx_mu1]',
    'Muon_pfRelIso04_all': 'Muon_pfRelIso04_all[Idx_mu1]',
    'Muon_corrected_MET_nom_mt':   'Muon_corrected_MET_nom_mt[Idx_mu1]'
}


class bkg_histos_fast(module):

    def __init__(self, selections, variables, dataType, xsec, inputFile, ptBins,etaBins,systDict, wptFunc, targetLumi = 1., clousureFlag=False, wpt=False, wptRew=False):

        # TH lists
        self.myTH1 = []
        self.myTH2 = []
        self.myTH3 = []

        # MC or DATA
        self.dataType = dataType

        self.selections = selections
        self.variablesUnvaried = variables
        self.ptBins = ptBins
        self.etaBins = etaBins
        self.systDict = systDict

        # pb to fb conversion
        self.xsec = xsec / 0.001
        self.targetLumi = targetLumi
        self.inputFile = inputFile
        self.clousureFlag = clousureFlag
        self.wpt= wpt #wpt calibration of Z ONLY
        self.wptFunc = wptFunc #function to rewight wpt
        self.wptRew = wptRew #activation of rewight wpt
        
        
    def dictReplacerVar(self, old,new) :
        self.old = old
        self.new = new
        
        mod_variables = {}
        mod_D2variables = {}
        
        for keys, value in self.variables["variables"].iteritems():
            mod_keys = keys.replace(self.old,self.new)
            temp = []
            for i in range(0,len(value)) :
                if (i==0 or i==4) :
                    temp.append(value[i].replace(self.old, self.new))
                else : 
                    temp.append(value[i])
            temp = tuple(temp)
            mod_variables[mod_keys] = temp

        for keys, value in self.variables["D2variables"].iteritems():
            mod_keys = keys.replace(self.old,self.new)
            temp = []
            for i in range(0,len(value)) :
                if (i>0 and i<7)  : 
                    temp.append(value[i])
                else :
                    temp.append(value[i].replace(self.old, self.new))
            temp = tuple(temp)
            mod_D2variables[mod_keys] = temp
        
        del self.variables["variables"]
        self.variables["variables"] =  mod_variables  
        
        del self.variables["D2variables"]
        self.variables["D2variables"] =  mod_D2variables  
    
    def bkg_histos(self,sKind,sName, selection,weight,colWeightName,wpt,wptFunc,wptRew) :

        self.systKind = sKind
        self.systName = sName
        self.colWeightName= colWeightName
        self.selection= selection
        self.weight= weight
        self.wpt = wpt
        self.wptFunc = wptFunc #function to rewight wpt
        self.wptRew = wptRew #activation of rewight wpt

        self.variables = copy.deepcopy(self.variablesUnvaried)

        
        varMT='Muon_corrected_MET_nom_mt[Idx_mu1]'
        varPt='Muon_corrected_pt[Idx_mu1]'
        
        #Apply sistematics            
        if self.dataType == 'MC':
            
            if "SF" in self.systKind or "Weight" in self.systKind:
                newWeight = weight.replace(self.systKind,self.systName)
                self.colWeightName = 'totweight_{}'.format(self.systName)
                self.d = self.d.Define(self.colWeightName, 'lumiweight*{}'.format(newWeight))
                if self.wptRew :
                    self.d = self.d.Define(self.colWeightName+'_wptRew', self.colWeightName+'*wptRew')
                    

            elif "nom" in self.systKind or "corrected" in self.systKind :
                print 'ERROR: there is a bug here (look to the code). nom and corrected not applied'
                # for some reason this part doesn't work, nom and corrected are not applied in hMain (but in hEWKSF yes), solve this issue!!
                self.dictReplacerVar(self.systKind,self.systName)
                self.selection = self.selection.replace(self.systKind,self.systName)            
                varMT = varMT.replace(self.systKind,self.systName)
                varPt = varPt.replace(self.systKind,self.systName)                          
            
            if self.wptRew :
                 self.colWeightName = self.colWeightName+'_wptRew'
        
        collName =  self.variables['prefix']

        if not self.wpt: 

            for var,tools in varDict_fast.iteritems():
                if ((not self.systName in var) and self.systName!='nom'): continue
                self.d = self.d.Define(collName+'_'+var,tools)
            
            
            hMain =  self.d.Filter(self.selection).Histo3D(('ptVSetaVSmt_'+self.systName, " ; {}; ".format('p_{T} VS Eta VS M_{T}'),  len(self.ptBins)-1, array('f',self.ptBins), len(self.etaBins)-1, array('f',self.etaBins),len(mtBinning)-1, array('f',mtBinning)),collName+'_Muon_corrected_pt',collName+'_Muon_eta', collName+'_Muon_corrected_MET_nom_mt', self.colWeightName)
            # hMain_highMt =  self.d.Filter(self.selection+'&& '+varMT+'>30').Histo2D(('ptVSeta_highMt'+self.systName, " ; {}; ".format('p_{T} VS Eta'),  len(self.etaBins)-1, array('f',self.etaBins), len(self.ptBins)-1, array('f',self.ptBins)),collName+'_Muon_corrected_pt',collName+'_Muon_eta', self.colWeightName)
            # hMain_lowMt = self.d.Filter(self.selection+'&& '+varMT+'<30').Histo2D(('ptVSeta_highMt'+self.systName, " ; {}; ".format('p_{T} VS Eta'),  len(self.etaBins)-1, array('f',self.etaBins), len(self.ptBins)-1, array('f',self.ptBins)),collName+'_Muon_corrected_pt',collName+'_Muon_eta', self.colWeightName)
            hEWKSF = self.d.Filter(self.selection+'&& '+varMT+'>30 && Muon_pfRelIso04_all[Idx_mu1]<0.01').Histo2D(('ptVSeta_IsoEWKSF_'+self.systName, " ; {}; ".format('p_{T} VS Eta'),  len(self.ptBins)-1, array('f',self.ptBins), len(self.etaBins)-1, array('f',self.etaBins)),collName+'_Muon_corrected_pt',collName+'_Muon_eta', self.colWeightName)
              
            # self.myTH2.append(hMain_highMt)
            # self.myTH2.append(hMain_lowMt)
            
            self.myTH3.append(hMain)
            self.myTH2.append(hEWKSF)

            if self.clousureFlag :
                for Cvar, tools in self.variables['ClousureVariables'].iteritems():
                        h3 = self.d.Filter(self.selection).Histo3D((collName+'_'+Cvar+'_'+self.systName, " ; {}; ".format(tools[0]),  tools[1],tools[2],tools[3],h_fake.GetNbinsY(), self.ptBins[0],self.ptBins[len(self.ptBins)-1], h_fake.GetNbinsX(), self.etaBins[0],self.etaBins[len(self.etaBins)-1]),collName+'_'+tools[4],collName+'_'+tools[5],collName+'_'+tools[6],self.colWeightName)
                        self.myTH3.append(h3)    
        
        else : #wpt=True
            for var,tools in self.variables['WptVariables'].iteritems():
                if ((not self.systName in var) and self.systName!='nom'): continue
                self.d = self.d.Define(collName+'_'+var,tools[4])   
                    
            h_fake = ROOT.TH2F("h_fake", "h_fake", len(self.etaBins)-1, array('f',self.etaBins), len(self.ptBins)-1, array('f',self.ptBins)) #fake hiso for binning only
            
            wptvar='RecoZ_Muon_corrected_pt' #z-pt and z-pt VS eta-mu VS pt-mu
            tools = self.variables['WptVariables'][wptvar]
            h = self.d.Filter(self.selection).Histo1D((collName+'_'+wptvar+'_'+self.systName, " ; {}; ".format(tools[0]), tools[1],tools[2], tools[3]), collName+'_'+wptvar, self.colWeightName)
            h3 = self.d.Filter(self.selection).Histo3D((collName+'_'+wptvar+'_'+self.systName+'_VS_eta_VS_pt', " ; {}; ".format(tools[0]),  tools[1],tools[2],tools[3],h_fake.GetNbinsY(), self.ptBins[0],self.ptBins[len(self.ptBins)-1], h_fake.GetNbinsX(), self.etaBins[0],self.etaBins[len(self.etaBins)-1]),collName+'_'+tools[4],collName+'_'+tools[5],collName+'_'+tools[6],self.colWeightName)
            self.myTH1.append(h) 
            self.myTH3.append(h3)
            
            wptvar='RecoZ_Muon_mass' #z-pt VS z-mass
            tools = self.variables['WptVariables'][wptvar]
            h2 = self.d.Filter(self.selection).Histo2D((collName+'_'+wptvar+'_'+self.systName+'_VS_Wpt', " ; {}; ".format(tools[0]),  tools[1],tools[2],tools[3],tools[5],tools[6],tools[7]),collName+'_'+tools[4],collName+'_'+tools[8],self.colWeightName)
            self.myTH2.append(h2) 
            
                
           
    def run(self,d):
        self.d = d
        # print "-"

        RDF = ROOT.ROOT.RDataFrame
        runs = RDF('Runs', self.inputFile)

        if self.dataType == 'MC':
            genEventSumw = runs.Sum("genEventSumw").GetValue()
        selection = self.selections[self.dataType]['cut']
        weight = self.selections[self.dataType]['weight']
       
        if self.dataType == 'MC':
            self.d = self.d.Define('lumiweight', '({L}*{xsec})/({genEventSumw})'.format(L=self.targetLumi, genEventSumw = genEventSumw, xsec = self.xsec)) 
            colWeightName = 'totweight'
            self.d = self.d.Define(colWeightName, 'lumiweight*{}'.format(weight))
            if self.wptRew :
                self.d = self.d.Define('wptRew', 'wptFunc->Eval(GenV_preFSR_qt)')
                self.d = self.d.Define(colWeightName+'_wptRew','wptRew*'+colWeightName)

        else:
            colWeightName = 'totweight'
            self.d = self.d.Define('totweight', '1')   
        
        self.bkg_histos("nom","nom",selection,weight,colWeightName, self.wpt, self.wptFunc,self.wptRew)    
        
        for sKind, sList in self.systDict.iteritems():
            for sName in sList :
                self.bkg_histos(sKind,sName,selection,weight,colWeightName, self.wpt, self.wptFunc,self.wptRew)
                
        return self.d

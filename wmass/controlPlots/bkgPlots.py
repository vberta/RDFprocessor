import math
import ROOT

import sys
sys.path.append('../../framework')
from module import *
from header import *

class bkgPlots(module):

    def __init__(self, selections, variables, dataType, xsec, inputFile, targetLumi = 1.):

        # TH lists
        self.myTH1 = []
        self.myTH2 = []
        self.myTH3 = []

        # MC or DATA
        self.dataType = dataType
        self.selections = selections
        self.variables = variables

        # pb to fb conversion
        self.xsec = xsec / 0.001
        self.targetLumi = targetLumi
        self.inputFile = inputFile

    def getSyst(self, syst):

        self.syst = syst # this is a dictionary. if empty, it corresponds to nominal

    def run(self,d):

        self.d = d

        RDF = ROOT.ROOT.RDataFrame
        runs = RDF('Runs', self.inputFile)

        if self.dataType == 'MC':
            genEventSumw = runs.Sum("genEventSumw").GetValue()
            print 'genEventSumw : '+'{:1.1f}'.format(genEventSumw)+' weighted events'
            print 'xsec         : '+'{:1.1f}'.format(self.xsec)+' pb'
            print 'lumiweight   : '+'{:1.8f}'.format((1.*self.xsec)/genEventSumw)+' (|Generator_weight| not accounted for)'

        selection = self.selections[self.dataType]['cut']
        weight = self.selections[self.dataType]['weight']

        # define mc specific weights (nominal)

        if self.dataType == 'MC':
            self.d = self.d.Define('lumiweight', '({L}*{xsec})/({genEventSumw})'.format(L=self.targetLumi, genEventSumw = genEventSumw, xsec = self.xsec)) \
                    .Define('totweight', 'lumiweight*{}'.format(weight))
        else:
            self.d = self.d.Define('totweight', '1')

        for nom,variations in self.syst.iteritems():
            if "SF" in nom or "Weight" in nom: #if this is a systematic of type "weight variations"

                print nom, "this is a systematic of type weight variations"
                if not self.dataType == 'MC': break

                for v in variations:
                    newWeight = weight.replace(nom,v)
                    print weight, newWeight

                    # define mc specific weights
                    if self.dataType == 'MC':
                        self.d = self.d.Define('totweight_{}'.format(v), 'lumiweight*{}[0]'.format(v))
                    else:
                        self.d = self.d.Define('totweight', '1') # to be checked what to do with data



                # loop over variables
                for Collection,dic in self.variables.iteritems():
                    collectionName = ''
                    if dic.has_key('newCollection') and dic['newCollection'] != '':
                        if 'index' in dic:
                            # define a new subcollection with all the columns of the original collection
                            self.d = self.defineSubcollectionFromIndex(dic['inputCollection'], dic['newCollection'], dic['index'], self.d)
                            collectionName = dic['newCollection']
                    else:
                        collectionName = dic['inputCollection']

                    for var,tools in dic['variables'].iteritems():


                        defined_columns = list(self.d.GetColumnNames())
                        defined_columns.extend(self.d.GetDefinedColumnNames())
                        if not(dic['inputCollection']+'_'+var in defined_columns) :
                            if ('_VS_' in var) and ('_TIMES_' in var)  : #a*b:c*d variable or a*b:c or a:b*c)
                                    splittedVarName = var.split('_VS_')
                                    if '_TIMES_' in splittedVarName[0] : #a*b:something
                                        secondSplitVarNameY = splittedVarName[0].split('_TIMES_')
                                        mergedVarNameY = dic['inputCollection']+'_'+secondSplitVarNameY[0]+'['+dic['index']+']'+'*'+dic['inputCollection']+'_'+secondSplitVarNameY[1]+'['+dic['index']+']'
                                    else : mergedVarNameY= dic['inputCollection']+'_'+splittedVarName[0]+'['+dic['index']+']'
                                    self.d = self.d.Define(collectionName+'_'+var+'_Y', '{vec}'.format(vec=mergedVarNameY))
                                    if '_TIMES_' in splittedVarName[1] : #something:c*d
                                        secondSplitVarNameX = splittedVarName[1].split('_TIMES_')
                                        mergedVarNameX = dic['inputCollection']+'_'+secondSplitVarNameX[0]+'['+dic['index']+']'+'*'+dic['inputCollection']+'_'+secondSplitVarNameX[1]+'['+dic['index']+']'
                                    else : mergedVarNameX= dic['inputCollection']+'_'+splittedVarName[1]+'['+dic['index']+']'
                                    self.d = self.d.Define(collectionName+'_'+var+'_X', '{vec}'.format(vec=mergedVarNameX))
                            elif '_VS_' in var : #a:b variables
                                    splittedVarName = var.split('_VS_')
                                    varNameY = dic['inputCollection']+'_'+splittedVarName[0]+'['+dic['index']+']'
                                    varNameX = dic['inputCollection']+'_'+splittedVarName[1]+'['+dic['index']+']'
                                    self.d = self.d.Define(collectionName+'_'+var+'_Y', '{vec}'.format(vec=varNameY))
                                    self.d = self.d.Define(collectionName+'_'+var+'_X', '{vec}'.format(vec=varNameX))
                            elif '_TIMES_' in var : #a*b variable
                                    splittedVarName = var.split('_TIMES_')
                                    mergedVarName = dic['inputCollection']+'_'+splittedVarName[0]+'['+dic['index']+']'+'*'+dic['inputCollection']+'_'+splittedVarName[1]+'['+dic['index']+']'
                                    self.d = self.d.Define(collectionName+'_'+var, '{vec}'.format(vec=mergedVarName))

                        for nom, variations in self.syst.iteritems():
                            for v in variations:
                                print Collection+'_'+var+'_'+v,collectionName+'_'+var,'totweight_{}'.format(v)


                                self.d = self.d.Filter(selection)
                                #cols = ROOT.vector('string')(); cols.push_back(collectionName+'_'+var);
                                #d2 = self.d.Display(cols)
                                #d2.Print()
                                if not('_VS_' in var) :
                                    h =self.d.Histo1D((Collection+'_'+var+'_'+v, " ; {}; ".format(tools[0]), tools[1],tools[2], tools[3]), collectionName+'_'+var, 'totweight_{}'.format(v))
                                    self.myTH1.append(h)
                                else :
                                    h =self.d.Histo2D((Collection+'_'+var+'_'+v, " ; {}; ".format(tools[0]), tools[1],tools[2], tools[3],tools[4],tools[5], tools[6]), collectionName+'_'+var+'_X',collectionName+'_'+var+'_Y','totweight_{}'.format(v))
                                    self.myTH2.append(h)


            else:
                print nom, "this is a systematic of type Up/Down variations"

                # loop over variables
                for Collection,dic in self.variables.iteritems():
                    collectionName = ''
                    if dic.has_key('newCollection') and dic['newCollection'] != '':
                        if 'index' in dic:
                            # define a new subcollection with all the columns of the original collection
                            if self.dataType == 'MC':
                                self.d = self.defineSubcollectionFromIndex(dic['inputCollection'], dic['newCollection'], dic['index'], self.d, self.syst)
                            else:
                                self.d = self.defineSubcollectionFromIndex(dic['inputCollection'], dic['newCollection'], dic['index'], self.d)

                            collectionName = dic['newCollection']
                    else:
                        collectionName = dic['inputCollection']

                    for var,tools in dic['variables'].iteritems():
                        # print "!!!!DEBUG!!!!!! inputCollection", dic['inputCollection']
                        # print "!!!!DEBUG!!!!!! var, tools,", var, tools
                        # print "!!!!DEBUG!!!!!! collectionName,Collection,", collectionName,Collection
                        defined_columns = list(self.d.GetColumnNames())
                        defined_columns.extend(self.d.GetDefinedColumnNames())
                        # print "!!!!DEBUG!!!!!! old  branch=", Collection+'_'+var
                        # print "!!!!DEBUG!!!!!! cadidato branch (s)=", collectionName+'_'+var
                        # print "!!!!DEBUG!!!!!! altro candidato nome=", dic['inputCollection']+'_'+var


                        # print "!!!!DEBUG!!!!!! defined_columns=", defined_columns
                        if not(dic['inputCollection']+'_'+var in defined_columns) :
                            if ('_VS_' in var) and ('_TIMES_' in var)  : #a*b:c*d variable or a*b:c or a:b*c)
                                    splittedVarName = var.split('_VS_')
                                    if '_TIMES_' in splittedVarName[0] : #a*b:something
                                        secondSplitVarNameY = splittedVarName[0].split('_TIMES_')
                                        mergedVarNameY = dic['inputCollection']+'_'+secondSplitVarNameY[0]+'['+dic['index']+']'+'*'+dic['inputCollection']+'_'+secondSplitVarNameY[1]+'['+dic['index']+']'
                                    else : mergedVarNameY= dic['inputCollection']+'_'+splittedVarName[0]+'['+dic['index']+']'
                                    self.d = self.d.Define(collectionName+'_'+var+'_Y', '{vec}'.format(vec=mergedVarNameY))
                                    if '_TIMES_' in splittedVarName[1] : #something:c*d
                                        secondSplitVarNameX = splittedVarName[1].split('_TIMES_')
                                        mergedVarNameX = dic['inputCollection']+'_'+secondSplitVarNameX[0]+'['+dic['index']+']'+'*'+dic['inputCollection']+'_'+secondSplitVarNameX[1]+'['+dic['index']+']'
                                    else : mergedVarNameX= dic['inputCollection']+'_'+splittedVarName[1]+'['+dic['index']+']'
                                    if splittedVarName[1]=='MET_pt' : mergedVarNameX=splittedVarName[1] #TEMPORARY: ONLY FOR MET (DIFFERENT COLLECTION)
                                    self.d = self.d.Define(collectionName+'_'+var+'_X', '{vec}'.format(vec=mergedVarNameX))
                            elif '_VS_' in var : #a:b variables
                                    splittedVarName = var.split('_VS_')
                                    varNameY = dic['inputCollection']+'_'+splittedVarName[0]+'['+dic['index']+']'
                                    varNameX = dic['inputCollection']+'_'+splittedVarName[1]+'['+dic['index']+']'
                                    if splittedVarName[1]=='MET_pt' : varNameX=splittedVarName[1] #TEMPORARY: ONLY FOR MET (DIFFERENT COLLECTION)
                                    self.d = self.d.Define(collectionName+'_'+var+'_Y', '{vec}'.format(vec=varNameY))
                                    self.d = self.d.Define(collectionName+'_'+var+'_X', '{vec}'.format(vec=varNameX))
                            elif '_TIMES_' in var : #a*b variable
                                    splittedVarName = var.split('_TIMES_')
                                    mergedVarName = dic['inputCollection']+'_'+splittedVarName[0]+'['+dic['index']+']'+'*'+dic['inputCollection']+'_'+splittedVarName[1]+'['+dic['index']+']'
                                    self.d = self.d.Define(collectionName+'_'+var, '{vec}'.format(vec=mergedVarName))


                        if not self.dataType == 'MC':
                            if not('_VS_' in var) :
                                h = self.d.Filter(selection).Histo1D((Collection+'_'+var, " ; {}; ".format(tools[0]), tools[1],tools[2], tools[3]), collectionName+'_'+var, 'totweight')
                                # print "il mio histo"
                                print h.GetName()
                                self.myTH1.append(h)
                            else :
                                h = self.d.Filter(selection).Histo2D((Collection+'_'+var, " ; {}; ".format(tools[0]), tools[1],tools[2], tools[3], tools[4],tools[5], tools[6]), collectionName+'_'+var+'_X',collectionName+'_'+var+'_Y', 'totweight')
                                # print "il mio histo"
                                print h.GetName()
                                self.myTH2.append(h)

                        else:
                            for nom, variations in self.syst.iteritems():

                                if len(variations)==0:
                                    if not('_VS_' in var) :
                                        h = self.d.Filter(selection).Histo1D((Collection+'_'+var, " ; {}; ".format(tools[0]), tools[1],tools[2], tools[3]), collectionName+'_'+var, 'totweight')
                                        self.myTH1.append(h)
                                    else :
                                        h = self.d.Filter(selection).Histo2D((Collection+'_'+var, " ; {}; ".format(tools[0]), tools[1],tools[2], tools[3],tools[4],tools[5], tools[6]), collectionName+'_'+var+'_X',collectionName+'_'+var+'_Y', 'totweight')
                                        self.myTH2.append(h)
                                else:

                                    for v in variations:
                                        if not nom in var: continue
                                        if not('_VS_' in var) :
                                            h = self.d.Filter(selection.replace(nom,v)).Histo1D((Collection+'_'+var.replace(nom,v), " ; {}; ".format(tools[0]), tools[1],tools[2], tools[3]), collectionName+'_'+var.replace(nom,v), 'totweight')
                                            self.myTH1.append(h)
                                        else :
                                            h = self.d.Filter(selection.replace(nom,v)).Histo2D((Collection+'_'+var.replace(nom,v), " ; {}; ".format(tools[0]), tools[1],tools[2], tools[3],tools[4],tools[5], tools[6]), collectionName+'_'+var+'_X',collectionName+'_'+var+'_Y', 'totweight')
                                            self.myTH2.append(h)

        return self.d

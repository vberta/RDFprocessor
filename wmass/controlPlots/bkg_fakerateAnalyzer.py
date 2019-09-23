import ROOT
from array import array
import math
import sys
import copy as copy
sys.path.append('../../framework')
from module import *
from header import *
from bkg_variables_standalone import *
from bkg_selections import *
from bkg_systematics import *




# class bkg_fakerateAnalyzer(module):
class bkg_fakerateAnalyzer:
    def __init__(self, ptBinning, etaBinning, outdir='./bkg', folder='./', norm = 1, varFake = 'Muon_pfRelIso04_all_corrected_MET_nom_mt', tightCut = 0.15, looseCut=40, fitOnTemplate=False, onData=True, nameSuff = '',slicing=True,systKind='nom',systName='nom',parabolaFit=False, EWSFfit=True)  :

        self.outdir = outdir
        self.folder = folder
        self.norm = norm
        self.ptBinning = ptBinning
        self.etaBinning = etaBinning
        self.varFake = varFake
        self.tightCut = tightCut
        self.looseCut = looseCut
        self.nameSuff = nameSuff
        self.systKind = systKind
        self.systName = systName
        self.parabolaFit = parabolaFit
        self.EWSFfit =EWSFfit

        self.fitOnTemplate = fitOnTemplate
        self.onData = onData
        self.slicing = slicing

        self.rootFilesRaw = []
        # self.rootFiles = []
        # self.relisoCUT = 0.15
        self.isoCUT = 5 # used for ratios VS Mt in preliminary studies only
        self.QCDmult = 1. #multiplication factor to QCD bkg, not implemented

        self.sampleList = ['WToMuNu','QCD','EWKbkg','Data','DataLike']
        # self.sampleList = ['WToMuNu','QCD','EWKbkg','DataLike']
        # self.sampleList = ['WToMuNu','QCD','Data','DataLike']

        self.signList = ['Plus','Minus']
        self.regionList = ['Signal','Sideband', 'Tot']
        self.varList = []
        for var in bkg_variables_standalone['D2variables'] : self.varList.append(var)
        # self.varList = ["pfRelIso04_all_VS_corrected_MET_nom_mt","pfRelIso04_all_TIMES_corrected_pt_VS_corrected_MET_nom_mt","pfRelIso04_all_VS_MET_pt","pfRelIso04_all_TIMES_corrected_pt_VS_MET_pt"]
        # self.varName = ["relIso_vs_Mt", "absIso_vs_Mt","relIso_vs_MET", "absIso_vs_MET"]
        # self.varName = ["relIso_vs_Mt", "absIso_vs_Mt"]
        self.varName = ["relIso_vs_Mt"]

        # self.ptBinningS = ['{:.2g}'.format(x) for x in self.ptBinning[1:]]
        # self.etaBinningS = ['{:.2g}'.format(x) for x in self.etaBinning[1:]]
        self.ptBinningS = ['{:.2g}'.format(x) for x in self.ptBinning[:-1]]
        self.etaBinningS = ['{:.2g}'.format(x) for x in self.etaBinning[:-1]]

        self.dataOpt = 'fake'
        if not self.onData : self.dataOpt = 'fakeMC'

        #open all the useful rootfile
        # for f in fileList
        #     rootFiles.append(ROOT.TFile.Open(self.folder+'/'+f))
        for f in range(len(self.sampleList)-1) :
            if (self.sampleList[f]!='DataLike') : self.rootFilesRaw.append(ROOT.TFile.Open(self.folder+'/'+self.sampleList[f]+'.root'))

        if self.slicing  :
            self.slicer()

        self.rootFiles = []
        for f in range(len(self.sampleList)-1) :
                if (self.sampleList[f]!='DataLike') : self.rootFiles.append(ROOT.TFile.Open(self.outdir+"/"+self.sampleList[f]+"_sliced"+".root"))


    def slicer(self) :
        print "-->Slincing:"
        # outlist = []
        print "files raw", self.rootFilesRaw
        for f in range(len(self.rootFilesRaw)) :

            var1DList = []
            for var,tools in bkg_variables_standalone['variables'].iteritems() :
                if not tools[5] : continue
                var1DList.append(var)

            output = ROOT.TFile(self.outdir+"/"+self.sampleList[f]+"_sliced"+".root","recreate")
            dirOutDict = {}
            histo1DDict ={}
            histo2DDict ={}
            histo3DDict = {}
            for s in self.signList :
                for r in self.regionList :
                    print "filename, region,sign", self.rootFilesRaw[f], r, s
                    if r== 'Tot': continue
                    dirOutDict[s+r] = output.mkdir('bkg_'+r+s+'/nom')
                    dirOutDict[s+r+'nom'] = dirOutDict[s+r].GetDirectory("nom")
                    dirOutDict[s+r+'nom'].cd()

                    #COPY 1D
                    for v in var1DList :
                        v_var = v
                        if "nom" in self.systKind or "corrected" in self.systKind :
                            if  not "Data.root" in self.rootFilesRaw[f].GetName() : v_var = v.replace(self.systKind,self.systName)
                        histo1DDict[s+r+v] = self.rootFilesRaw[f].Get('bkg_'+r+s+'/nom/bkgSel_'+v_var+'_'+self.systName)
                        print "histo name", 'bkgSel_'+v_var+'_'+self.systName
                        histo1DDict[s+r+v].SetName('bkgSel_'+v)
                        histo1DDict[s+r+v].Write()

                    #SLICE 2D
                    name2D = 'Muon_corrected_MET_nom_mt_VS_eta'
                    # name2D = 'MET_pt_VS_eta'
                    name2D_var = name2D
                    if "nom" in self.systKind or "corrected" in self.systKind :
                        if not "Data.root" in self.rootFilesRaw[f].GetName() : name2D_var = name2D.replace(self.systKind,self.systName)
                    histo2DDict[s+r] = self.rootFilesRaw[f].Get('bkg_'+r+s+'/nom/bkgSel_'+name2D_var+'_'+self.systName)
                    for eta in range(1, histo2DDict[s+r].GetNbinsY()+1) :
                        histo2DDict[s+r+str(eta)] = histo2DDict[s+r].ProjectionX("bkgSel_Muon_corrected_MET_nom_mt_"+self.etaBinningS[eta-1],eta,eta+1,"e")
                        histo2DDict[s+r+str(eta)].Write()
                    histo2DDict[s+r+"int"] = histo2DDict[s+r].ProjectionX("bkgSel_Muon_corrected_MET_nom_mt",1,histo2DDict[s+r].GetNbinsY(),"e")
                    histo2DDict[s+r+"int"].Write()

                    #SLICE 3D
                    for v in self.varList :
                        v_var = v
                        if "nom" in self.systKind or "corrected" in self.systKind :
                            if not "Data.root" in self.rootFilesRaw[f].GetName() : v_var = v.replace(self.systKind,self.systName)
                        histo3DDict[s+r+v+"int"] =  self.rootFilesRaw[f].Get('bkg_'+r+s+'/nom/bkgSel_'+v_var+'_VS_eta_'+self.ptBinningS[0]+'_'+self.systName).Project3D("yxe")
                        for p in self.ptBinningS :
                            histo3DDict[s+r+v] = self.rootFilesRaw[f].Get('bkg_'+r+s+'/nom/bkgSel_'+v_var+'_VS_eta_'+p+'_'+self.systName)
                            for eta in range(0, histo3DDict[s+r+v].GetNbinsZ()+1) :
                                histo3DDict[s+r+v].GetZaxis().SetRange(eta,eta+1)
                                histo3DDict[s+r+v+str(eta)] = histo3DDict[s+r+v].Project3D("yxe")
                                histo3DDict[s+r+v+str(eta)].SetName('bkgSel_'+v+'_'+self.etaBinningS[eta-1]+'_'+p)
                                histo3DDict[s+r+v+str(eta)].Write()
                                if(p!=self.ptBinningS[0]) :
                                    histo3DDict[s+r+v+"int"].Add(histo3DDict[s+r+v+str(eta)])
                        histo3DDict[s+r+v+"int"].SetName('bkgSel_'+v)
                        histo3DDict[s+r+v+"int"].Write()

            output.Close()
            # outlist.append(output)
        # return outlist


    def ratio_2Dto1D(self,histo,isocut =0.15,name = "histoRate") : #histo = 2D isto iso:Mt, isocut=tight def., name=output histo name
        #this func. produce an histogram of fake or prompt rate in fuction of Mt (to verify ABCD assumption)
        self.histo = histo
        self.isocut = isocut
        self.name = name

        isoMin= self.histo.GetYaxis().GetBinCenter(1)-self.histo.GetYaxis().GetBinWidth(1)/2
        binsize = self.histo.GetYaxis().GetBinWidth(1)
        Ncut=(self.isocut-isoMin)/binsize
        # print "DEBUG RATIO", isoMin, binsize, Ncut
        Ncut = int(Ncut)
        # print isocut, Ncut, isoMin
        # histoDen = self.histo.ProjectionX("histoDen",Ncut,-1)

        histoDen = self.histo.ProjectionX("histoDen")
        histoNum = self.histo.ProjectionX("histoNum",0,Ncut-1)
        histoRate = histoNum.Clone(self.name)
        histoRate.Divide(histoNum,histoDen,1,1)
        # print "preliminary, NUM, DEN", histoNum.GetEntries(),histoDen.GetEntries()
        return histoRate

    def Fit4ScaleFactorEW(self,mtDict, sign, eta,datakind,pt='') :
        self.mtDict = mtDict
        self.sign = sign
        self.datakind = datakind
        self.eta = eta
        self.pt = pt

        outlist =[] #sig,bkg, chi2,chi2err

        class linearHistoFit:
            def __call__(self, x, parameters):
                s = parameters[0] # weight signal
                b = parameters[1] # weight bkg
                # c = parameters[2]
                x = x[0]
                ysig = hsig.GetBinContent(hsig.GetXaxis().FindFixBin(x));
                ybkg = hbkg.GetBinContent(hbkg.GetXaxis().FindFixBin(x));
                y = s*ysig+b*ybkg
                # print "value y",y
                return y

        hsig = mtDict[self.pt+self.eta+self.sign+'WToMuNuTot'].Clone()
        hsig.Add(mtDict[self.pt+self.eta+self.sign+'EWKbkgTot'])
        hbkg = mtDict[self.pt+self.eta+self.sign+'QCDTot'].Clone()
        hsig.Rebin(4)
        hbkg.Rebin(4)
        fitFunc = ROOT.TF1("fitFunc", linearHistoFit(),0,120,2)
        # fitFunc = ROOT.TF1("fitFunc", linearHistoFit,0,120,2)
        fitFunc.SetParameters(1,1,0)
        fitFunc.SetParNames("sig","bkg","const")
        hdata = mtDict[self.pt+self.eta+self.sign+self.datakind+'Tot'].Clone()
        hdata.Rebin(4)
        hdata.Fit(fitFunc,"Q","",0,120)
        
        sumDiffPre = 0 
        sumDiffPost = 0 
        for x in range(0,120) :
            btot = mtDict[self.pt+self.eta+self.sign+self.datakind+'Tot'].GetBinContent(hsig.GetXaxis().FindFixBin(x))
            bsig = hsig.GetBinContent(hsig.GetXaxis().FindFixBin(x));
            bbkg = hbkg.GetBinContent(hbkg.GetXaxis().FindFixBin(x));
            ss = fitFunc.GetParameter(0)
            bb = fitFunc.GetParameter(1)
            sumDiffPre = sumDiffPre+ (bsig+bbkg-btot)**2
            sumDiffPost = sumDiffPost+ (ss*bsig+bb*bbkg-btot)**2

        # print "FIT RESULTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Post fit, values:", fitFunc.GetParameter(0), fitFunc.GetParameter(1), fitFunc.GetChisquare()/fitFunc.GetNDF(), "residual prefit/postfit =",sumDiffPre/sumDiffPost #PRINT THIS ONE
        # print "Pre fit values: (tot,w,QCD)   ", mtDict[self.sign+self.datakind+'Tot'].GetBinContent(10), mtDict[self.sign+'WToMuNuTot'].GetBinContent(10)+mtDict[self.sign+'QCDTot'].GetBinContent(10)+mtDict[self.sign+'EWKbkgTot'].GetBinContent(10)

        outlist = [fitFunc.GetParameter(0), fitFunc.GetParError(0),fitFunc.GetParameter(1), fitFunc.GetParError(1),fitFunc.GetChisquare()/fitFunc.GetNDF(),math.sqrt(2*fitFunc.GetNDF())/fitFunc.GetNDF() ]
        # return fitFunc.GetParameter(1)
        return outlist


    def isolationAna(self, hdict,  loosecut=40, varname = 'Muon_pfRelIso04_all_corrected_MET_nom_mt', kind = 'fake') :
        self.loosecut = loosecut
        self.varname = varname
        self.kind = kind # fake = calculate the fakerate (measurement ABCD), prompt = the promptrate (from MC), validation = the fakerate on MC QCD, fakeMC = fakerate from dataLike (MC) SEE DICT BELOW
        self.hdct = hdict

        kindDict = {
            'fake' : 'Data',
            'validation' : 'QCD',
            'fakeMC' : 'DataLike' ,
            'prompt' : 'WToMuNu',
            'EWKbkg' : 'EWKbkg',
        }
        datakind = kindDict[self.kind]

        hIsos = {}

        for s in self.signList :
            for e in self.etaBinningS :
                for p in self.ptBinningS :
                    # hIso = ROOT.TH1F("hIso_{kind}_{sign}_{eta}_{pt}".format(kind=self.kind,sign=s,eta=e, pt=p),"hIso_{kind}_{sign}_{eta}_{pt}".format(kind=self.kind,sign=s,eta=e, pt=p),400,0,0.5)
                    mtMin= hdict[p+e+s+varname+datakind+'Tot'].GetXaxis().GetBinCenter(1)-hdict[p+e+s+varname+datakind+'Tot'].GetXaxis().GetBinWidth(1)/2
                    binsizeLoose = hdict[p+e+s+varname+datakind+'Tot'].GetXaxis().GetBinWidth(1)
                    NcutLoose=(self.loosecut-mtMin)/binsizeLoose
                    NcutLoose = int(NcutLoose)

                    hIso = hdict[p+e+s+varname+datakind+'Tot'].ProjectionY("Iso_{kind}_{sign}_{eta}_{pt}".format(kind=self.kind,sign=s,eta=e, pt=p),0,NcutLoose-1)

                    hIsos[p+e+s] = (hIso)
        return hIsos

    def differential_fakerate(self, hdict, mtDict, tightcut = 0.15, loosecut=40, varname = 'pfRelIso04_all_corrected_MET_nom_mt', kind = 'fake', EWSFfit=True, highMtCut=90,parabolaFit=False ) :
        self.loosecut = loosecut
        self.tightcut = tightcut
        self.varname = varname
        self.kind = kind # fake = calculate the fakerate (measurement ABCD), prompt = the promptrate (from MC), validation = the fakerate on MC QCD, fakeMC = fakerate from dataLike (MC) SEE DICT BELOW
        self.hdct = hdict
        self.mtDict =mtDict
        self.EWSFfit = EWSFfit
        self.highMtCut = highMtCut
        self.parabolaFit = parabolaFit

        kindDict = {
            'fake' : 'Data',
            'validation' : 'QCD',
            'fakeMC' : 'DataLike' ,
            'prompt' : 'WToMuNu',
            'validationSigReg' : 'QCD',
            'promptSideband' : 'WToMuNu',
        }
        datakind = kindDict[self.kind]

        hFakes = {}
        h2Fakes = {}
        hEWSF_Fit = {}
        hTempl_Fit = {}
        # TH2F h2Fakes = TH2F("h2Fakes","h2Fakes",len(self.etaBinning)-1, array('f',self.etaBinning), len(self.ptBinning)-1, array('f',self.ptBinning) )
        # h2Fakes[0] = TH2F("h2Fakes_plus","h2Fakes_plus",len(self.etaBinning)-1, array('f',self.etaBinning), len(self.ptBinning)-1, array('f',self.ptBinning) )
        # h2Fakes[1] = TH2F("h2Fakes_minus","h2Fakes_minus",len(self.etaBinning)-1, array('f',self.etaBinning), len(self.ptBinning)-1, array('f',self.ptBinning) )
        for s in self.signList :
            h2Fakes_sign = ROOT.TH2F("h2Fakes_{kind}_{sign}".format(kind=self.kind,sign=s),"h2Fakes_{kind}_{sign}".format(kind=self.kind,sign=s),len(self.etaBinning)-1, array('f',self.etaBinning), len(self.ptBinning)-1, array('f',self.ptBinning) )
            if self.kind == 'fake' or self.kind == 'fakeMC' :
                hEWSF_chi2 = ROOT.TH1F("hEWSF_chi2_{kind}_{sign}".format(kind=self.kind,sign=s),"hEWSF_chi2_{kind}_{sign}".format(kind=self.kind,sign=s), len(self.etaBinning)-1, array('f',self.etaBinning), )
                hEWSF_bkg = ROOT.TH1F("hEWSF_bkg_{kind}_{sign}".format(kind=self.kind,sign=s),"hEWSF_bkg_{kind}_{sign}".format(kind=self.kind,sign=s), len(self.etaBinning)-1, array('f',self.etaBinning), )
                hEWSF_sig = ROOT.TH1F("hEWSF_sig_{kind}_{sign}".format(kind=self.kind,sign=s),"hEWSF_sig_{kind}_{sign}".format(kind=self.kind,sign=s), len(self.etaBinning)-1, array('f',self.etaBinning), )
            hTempl_chi2 = ROOT.TH1F("hTempl_chi2_{kind}_{sign}".format(kind=self.kind,sign=s),"hTempl_chi2_{kind}_{sign}".format(kind=self.kind,sign=s), len(self.etaBinning)-1, array('f',self.etaBinning), )
            hTempl_slope = ROOT.TH1F("hTempl_slope_{kind}_{sign}".format(kind=self.kind,sign=s),"hTempl_slope_{kind}_{sign}".format(kind=self.kind,sign=s), len(self.etaBinning)-1, array('f',self.etaBinning), )
            hTempl_offset = ROOT.TH1F("hTempl_offset_{kind}_{sign}".format(kind=self.kind,sign=s),"hTempl_offset_{kind}_{sign}".format(kind=self.kind,sign=s), len(self.etaBinning)-1, array('f',self.etaBinning), )
            hTempl_2deg = ROOT.TH1F("hTempl_2deg_{kind}_{sign}".format(kind=self.kind,sign=s),"hTempl_2deg_{kind}_{sign}".format(kind=self.kind,sign=s), len(self.etaBinning)-1, array('f',self.etaBinning), )


            for e in self.etaBinningS :
                hFakes_pt = ROOT.TH1F("hFakes_pt_{kind}_{sign}_{eta}".format(kind=self.kind,sign=s,eta=e),"hFakes_pt_{kind}_{sign}_{eta}".format(kind=self.kind,sign=s,eta=e), len(self.ptBinning)-1, array('f',self.ptBinning) )
                scaleFactorEW =1
                if self.kind == 'fake' or self.kind == 'fakeMC' :
                    # print "PRE FIT (sign, eta, kind))", s, e, datakind
                    # scaleFactorEW=self.Fit4ScaleFactorEW(mtDict=mtDict,sign=s,eta=e,datakind=datakind)
                    scaleFactorEWPars = self.Fit4ScaleFactorEW(mtDict=self.mtDict,sign=s,eta=e,datakind=datakind)
                    if(self.EWSFfit) :
                        scaleFactorEW=scaleFactorEWPars[0]
                        # scaleFactorEW = 1
                        # scaleFactorEW = scaleFactorEW-0.1*scaleFactorEW #variation of EWSF
                    else :
                        minBin = mtDict[e+s+'WToMuNuTot'].GetXaxis().FindBin(self.highMtCut)
                        maxBin = mtDict[e+s+'WToMuNuTot'].GetSize()-1 #number of bins (overflow included, -1 is to the underflow)
                        EWKInt = mtDict[e+s+'WToMuNuTot'].Integral(minBin,maxBin)
                        EWKInt = EWKInt + mtDict[e+s+'EWKbkgTot'].Integral(minBin,maxBin)
                        dataInt = mtDict[e+s+datakind+'Tot'].Integral(minBin,maxBin)
                        scaleFactorEW = EWKInt/dataInt
                        # print "SCALE FACTOR (eta,sign)",e, s, ",   VALUE=", scaleFactorEW, "  fit one=",scaleFactorEWPars[0], "   ratio (int/fit)=",scaleFactorEW/scaleFactorEWPars[0] #PRINT THIS ONE
                    hEWSF_bkg.SetBinContent(self.etaBinningS.index(e)+1,scaleFactorEWPars[0])
                    hEWSF_bkg.SetBinError(self.etaBinningS.index(e)+1,scaleFactorEWPars[1])
                    hEWSF_sig.SetBinContent(self.etaBinningS.index(e)+1,scaleFactorEWPars[2])
                    hEWSF_sig.SetBinError(self.etaBinningS.index(e)+1,scaleFactorEWPars[3])
                    hEWSF_chi2.SetBinContent(self.etaBinningS.index(e)+1,scaleFactorEWPars[4])
                    hEWSF_chi2.SetBinError(self.etaBinningS.index(e)+1,scaleFactorEWPars[5])

                
                    hEWSF_chi2_pt = ROOT.TH1F("hEWSF_chi2_{kind}_{sign}_{eta}".format(kind=self.kind,sign=s,eta=e),"hEWSF_chi2_{kind}_{sign}_{eta}".format(kind=self.kind,sign=s,eta=e), len(self.ptBinning)-1, array('f',self.ptBinning), )
                    hEWSF_bkg_pt = ROOT.TH1F("hEWSF_bkg_{kind}_{sign}_{eta}".format(kind=self.kind,sign=s,eta=e),"hEWSF_bkg_{kind}_{sign}_{eta}".format(kind=self.kind,sign=s,eta=e), len(self.ptBinning)-1, array('f',self.ptBinning), )
                    hEWSF_sig_pt = ROOT.TH1F("hEWSF_sig_{kind}_{sign}_{eta}".format(kind=self.kind,sign=s,eta=e),"hEWSF_sig_{kind}_{sign}_{eta}".format(kind=self.kind,sign=s,eta=e), len(self.ptBinning)-1, array('f',self.ptBinning), )
                    hEWSF_count_pt = ROOT.TH1F("hEWSF_count_{kind}_{sign}_{eta}".format(kind=self.kind,sign=s,eta=e),"hEWSF_count_{kind}_{sign}_{eta}".format(kind=self.kind,sign=s,eta=e), len(self.ptBinning)-1, array('f',self.ptBinning), )
                for p in self.ptBinningS :
                    scaleFactorEWPars_p = self.Fit4ScaleFactorEW(mtDict=self.mtDict,sign=s,eta=e,datakind=datakind,pt=p)
                    if(self.EWSFfit) :
                        scaleFactorEW=scaleFactorEWPars_p[0]
                    if self.kind == 'fake' or self.kind == 'fakeMC' :
                        hEWSF_bkg_pt.SetBinContent(self.ptBinningS.index(p)+1,scaleFactorEWPars_p[0])
                        hEWSF_bkg_pt.SetBinError(self.ptBinningS.index(p)+1,scaleFactorEWPars_p[1])
                        hEWSF_sig_pt.SetBinContent(self.ptBinningS.index(p)+1,scaleFactorEWPars_p[2])
                        hEWSF_sig_pt.SetBinError(self.ptBinningS.index(p)+1,scaleFactorEWPars_p[3])
                        hEWSF_chi2_pt.SetBinContent(self.ptBinningS.index(p)+1,scaleFactorEWPars_p[4])
                        hEWSF_chi2_pt.SetBinError(self.ptBinningS.index(p)+1,scaleFactorEWPars_p[5])       
                        # print "e,p, s=", e,p,s, ", value=", scaleFactorEWPars_p[0]            
                        
                    hsubtract= hdict[p+e+s+varname+datakind+'Tot'].Clone(p+'_'+e+'_'+'datakind'+s+'_'+varname)
                    # if self.kind == 'fake' or self.kind == 'fakeMC' :
                    #     hsubtract.Add(hdict[p+e+s+varname+'EWKbkg'+'Tot'],-1)
                    #     hsubtract.Add(hdict[p+e+s+varname+'WToMuNu'+'Tot'],-1)

                    isoMin= hsubtract.GetYaxis().GetBinCenter(1)-hsubtract.GetYaxis().GetBinWidth(1)/2
                    binsizeTight = hsubtract.GetYaxis().GetBinWidth(1)
                    NcutTight=(self.tightcut-isoMin)/binsizeTight
                    NcutTight = int(NcutTight)


                    mtMin= hsubtract.GetXaxis().GetBinCenter(1)-hsubtract.GetXaxis().GetBinWidth(1)/2
                    binsizeLoose = hsubtract.GetXaxis().GetBinWidth(1)
                    NcutLoose=(self.loosecut-mtMin)/binsizeLoose
                    NcutLoose = int(NcutLoose)
                    # print "cuts=, ", NcutTight, NcutLoose

                    numErr = ROOT.Double(0)
                    denErr = ROOT.Double(0)
                    antiNumErr = ROOT.Double(0) #antinum = b region (not tight)
                    fake_err =0
                    if (self.kind=='fake' or self.kind=='validation' or self.kind=='fakeMC' or self.kind=='promptSideband') :
                        # print "here: cuts (t,l)", NcutTight, NcutLoose
                        den = hsubtract.ProjectionX("histoDen",0,-1,"e").IntegralAndError(0,NcutLoose-1,denErr)
                        num = hsubtract.ProjectionX("histoNum",0,NcutTight-1, "e").IntegralAndError(0,NcutLoose-1,numErr)
                        antiNum = hsubtract.ProjectionX("histoNum",NcutTight,-1, "e").IntegralAndError(0,NcutLoose-1,antiNumErr)

                        # scaleFactorLumi = 1
                        # if(num!= 0 and den!=0) :
                        #     scaleFactorLumi= numErr*numErr/num
                        #     num4err = num /scaleFactorLumi
                        #     den4err = den /scaleFactorLumi
                        #     fake_err = 1/den4err*math.sqrt(num4err*(1-num4err/den4err))*scaleFactorLumi #standard eff error rewighted on scale factor
                        # print "SCALE FACTOR LUMI:", scaleFactorLumi


                        if self.kind == 'fake' or self.kind == 'fakeMC' : #not QCD
                            numErr_EWKbkg = ROOT.Double(0)
                            denErr_EWKbkg = ROOT.Double(0)
                            numErr_WToMuNu = ROOT.Double(0)
                            denErr_WToMuNu = ROOT.Double(0)
                            den_EWKbkg = hdict[p+e+s+varname+'EWKbkg'+'Tot'].ProjectionX("histoDen",0,-1,"e").IntegralAndError(0,NcutLoose-1,denErr_EWKbkg)
                            num_EWKbkg = hdict[p+e+s+varname+'EWKbkg'+'Tot'].ProjectionX("histoNum",0,NcutTight-1, "e").IntegralAndError(0,NcutLoose-1,numErr_EWKbkg)
                            den_WToMuNu = hdict[p+e+s+varname+'WToMuNu'+'Tot'].ProjectionX("histoDen",0,-1,"e").IntegralAndError(0,NcutLoose-1,denErr_WToMuNu)
                            num_WToMuNu = hdict[p+e+s+varname+'WToMuNu'+'Tot'].ProjectionX("histoNum",0,NcutTight-1, "e").IntegralAndError(0,NcutLoose-1,numErr_WToMuNu)

                            # print "PRE FIT (sign, kind))", s, datakind
                            # scaleFactorEW = self.Fit4ScaleFactorEW(mtDict,s,datakind)
                            # scaleFactorEW=1
                            
                            #evaluaiton of the correction pt-dependent of the scaleFactorEW
                            w_only_region_data_Err_w = ROOT.Double(0)
                            w_only_region_data_Err_qcd = ROOT.Double(0)
                            w_only_region_MC_Err_w = ROOT.Double(0)
                            w_only_region_MC_Err_ewk = ROOT.Double(0)
                            w_only_region_data= hsubtract.ProjectionX("histoNum",0,1, "e").IntegralAndError(NcutLoose-1,-1,w_only_region_data_Err_w)-hdict[p+e+s+varname+'QCD'+'Tot'].ProjectionX("histoNum",0,1, "e").IntegralAndError(NcutLoose-1,-1,w_only_region_data_Err_qcd)
                            w_only_region_MC = hdict[p+e+s+varname+'WToMuNu'+'Tot'].ProjectionX("histoNum",0,1, "e").IntegralAndError(NcutLoose-1,-1,w_only_region_MC_Err_w)+hdict[p+e+s+varname+'EWKbkg'+'Tot'].ProjectionX("histoNum",0,1, "e").IntegralAndError(NcutLoose-1,-1,w_only_region_MC_Err_ewk)
                            EWKSF_pt = w_only_region_data/w_only_region_MC
                            EWKSF_pt_err_num=math.sqrt(w_only_region_data_Err_w**2+w_only_region_data_Err_qcd**2)
                            EWKSF_pt_err_den=math.sqrt(w_only_region_MC_Err_w**2+w_only_region_MC_Err_ewk**2)
                            EWKSF_pt_err = 1/(w_only_region_MC**2)*math.sqrt((w_only_region_MC**2)*(EWKSF_pt_err_num**2)+(w_only_region_data**2)*(EWKSF_pt_err_den**2))
                            hEWSF_count_pt.SetBinContent(self.ptBinningS.index(p)+1,EWKSF_pt)
                            hEWSF_count_pt.SetBinError(self.ptBinningS.index(p)+1,EWKSF_pt_err)

                            # print "pt dependent EWSF=", EWKSF_pt, "pt=",p, "MC QCD % = ", hdict[p+e+s+varname+'QCD'+'Tot'].ProjectionX("histoNum",0,1, "e").Integral(NcutLoose-1,-1)/hsubtract.ProjectionX("histoNum",0,1, "e").Integral(NcutLoose-1,-1)
                            # EWKSF_pt =1 
                            
                            # den = den - (den_EWKbkg + den_WToMuNu)*scaleFactorEW
                            # num = num - (num_EWKbkg + num_WToMuNu)*scaleFactorEW
                            den = den - (den_EWKbkg + den_WToMuNu)*EWKSF_pt
                            num = num - (num_EWKbkg + num_WToMuNu)*EWKSF_pt
                            # denErr = math.sqrt(denErr**2 + denErr_EWKbkg**2 + denErr_WToMuNu**2)
                            # numErr = math.sqrt(numErr**2 + numErr_EWKbkg**2 + numErr_WToMuNu**2)

                            #ERROR ON SUM EVALUATION
                            antiNumErr_EWKbkg = ROOT.Double(0)
                            antiNumErr_WToMuNu = ROOT.Double(0)
                            antiNum_EWKbkg = hdict[p+e+s+varname+'EWKbkg'+'Tot'].ProjectionX("histoNum",NcutTight,-1, "e").IntegralAndError(0,NcutLoose-1,antiNumErr_EWKbkg)
                            antiNum_WToMuNu = hdict[p+e+s+varname+'WToMuNu'+'Tot'].ProjectionX("histoNum",NcutTight,-1, "e").IntegralAndError(0,NcutLoose-1,antiNumErr_WToMuNu)
                            
                        
                            antiNum = antiNum - (antiNum_EWKbkg + antiNum_WToMuNu)*scaleFactorEW

                            # if(num!=0 and den!=0 and antiNum!=0) :
                            #     num4err = num /scaleFactorLumi
                            #     den4err = den /scaleFactorLumi
                            #     antiNum4err = antiNum /scaleFactorLumi

                                # d_num2 = (scaleFactorLumi**2) * ((num4err**2)+(scaleFactorEW**2)*((num4err_EWKbkg**2)+(num4err_WToMuNu**2))) #error SQUARE on num4err
                                # d_antiNum2 = (scaleFactorLumi**2) * ((antiNum4err**2)+(scaleFactorEW**2)*((antiNum4err_EWKbkg**2)+(antiNum4err_WToMuNu**2))) #error SQUARE on antinum4err
                                # fake_err=(1/(den4err**2))*math.sqrt(d_num2*(antiNum4err**2)+d_antiNum2*(num4err**2))


                        # for x in range(hsubtract.GetNbinsX()) :
                        #     for y in range(hsubtract.GetNbinsY()) :
                        #             if(hsubtract.GetBinContent(x,y)>0) : radice = math.sqrt(hsubtract.GetBinContent(x,y))
                        #             else : radice = 0
                        #             print "esempio bin ------------------------------", hsubtract.GetBinContent(x,y),  hsubtract.GetBinError(x,y),radice


                    else : #prompt calulated in signal region (or valdiationSigReg)
                        den = hsubtract.ProjectionX("histoDen",0,-1,"e").IntegralAndError(NcutLoose-1,-1,denErr)
                        num = hsubtract.ProjectionX("histoNum",0,NcutTight-1, "e").IntegralAndError(NcutLoose-1,-1,numErr)
                        # print"bin 3,val + content,", hsubtract.ProjectionX("histoDen",0,-1,"e").GetBinContent(3), hsubtract.ProjectionX("histoDen",0,-1,"e").GetBinError(3)
                        # print"bin 3,val + content, no e opt", hsubtract.ProjectionX("histoDen").GetBinContent(3), hsubtract.ProjectionX("histoDen").GetBinError(3)
                        antiNum = hsubtract.ProjectionX("histoNum",NcutTight,-1, "e").IntegralAndError(NcutLoose-1,-1,antiNumErr)
                        # scaleFactorLumi = 1
                        # if(num!= 0 and den!=0) :
                        #     scaleFactorLumi= numErr*numErr/num
                        #     num = num /scaleFactorLumi
                        #     den = den /scaleFactorLumi
                        #     fake_err = 1/den*math.sqrt(num*(1-num/den))*scaleFactorLumi #standard eff error rewighted on scal factor
                        # print "SCALE FACTOR LUMI :", scaleFactorLumi

                    # print("kind", self.kind, "eta,pt", e, p, "num,den", num, den, "s",s)
                    if(den == 0) :
                        fake = 0
                        print "WARNING: fake rate den = 0 --> num=", num, "data kind=",self.kind, "(pt,eta,sign)=",p,e,s
                    if(num==0) :
                        fake = 0
                        print "WARNING: fake rate num = 0 --> den=", den, "data kind=",self.kind, "(pt,eta,sign)=",p,e,s
                    else:
                        # print "Ok: fake rate --> num/den=", num, den, num/den, "data kind=",self.kind, "(pt,eta,sign)=",p,e,s
                        fake = num/den
                        # print num, numErr, math.sqrt(num),den, denErr, math.sqrt(den)
                        # fake_err = 1/(den**2)*math.sqrt((numErr**2)*(den**2)+(denErr**2)*(num**2)) #UNCORRELATED!!!!
                        fake_err=(1/(den**2))*math.sqrt((numErr**2)*(antiNum**2)+(antiNumErr**2)*(num**2))

                    hFakes_pt.SetBinContent(self.ptBinningS.index(p)+1,fake)
                    h2Fakes_sign.SetBinContent(self.etaBinningS.index(e)+1, self.ptBinningS.index(p)+1,fake)
                    hFakes_pt.SetBinError(self.ptBinningS.index(p)+1,fake_err)
                    h2Fakes_sign.SetBinError(self.etaBinningS.index(e)+1, self.ptBinningS.index(p)+1,fake_err)
                
                if self.kind == 'fake' or self.kind == 'fakeMC' :
                    hEWSF_Fit["EWSF_chi2"+s+e] = hEWSF_chi2_pt
                    hEWSF_Fit["EWSF_sig"+s+e] = hEWSF_sig_pt
                    hEWSF_Fit["EWSF_bkg"+s+e] = hEWSF_bkg_pt
                    hEWSF_Fit["EWSF_count"+s+e] = hEWSF_count_pt
                    # print "EWSF_count"+s+e
                
                hFakes[s+e] = hFakes_pt

                if not self.parabolaFit :
                    fitFake = ROOT.TF1("fitFake", 'pol1',30,65,2)
                    fitFake.SetParameters(0.5,0.1)
                    fitFake.SetParNames("offset","slope")
                else :
                    fitFake = ROOT.TF1("fitFake", 'pol2',30,65,2)
                    fitFake.SetParameters(0.5,0.1,-0.1)
                    fitFake.SetParNames("offset","slope",'2deg')                    
                # fitFake = ROOT.TF1("fitFake", 'pol2',30,65,2)
                # fitFake.SetParameters(0.5,0.1,0.1)
                # fitFake.SetParNames("offset","slope","sndDeg")
                hFakes_pt.Fit(fitFake,"Q","",0,120)
                hFakes[s+e+'offset']=fitFake.GetParameter(0)
                hFakes[s+e+'slope']=fitFake.GetParameter(1)#+0.1*fitFake.GetParameter(1) #systemtatic
                hFakes[s+e+'offsetErr']=fitFake.GetParError(0)
                hFakes[s+e+'slopeErr']=fitFake.GetParError(1)
                if self.parabolaFit :
                    hFakes[s+e+'2deg']=fitFake.GetParameter(2)
                    hFakes[s+e+'2degErr']=fitFake.GetParError(2)
                if(fitFake.GetNDF()>0) :
                    hTempl_chi2.SetBinContent(self.etaBinningS.index(e)+1,fitFake.GetChisquare()/fitFake.GetNDF())
                    hTempl_chi2.SetBinError(self.etaBinningS.index(e)+1,math.sqrt(2*fitFake.GetNDF())/fitFake.GetNDF())
                    hTempl_slope.SetBinContent(self.etaBinningS.index(e)+1,fitFake.GetParameter(1))
                    hTempl_slope.SetBinError(self.etaBinningS.index(e)+1,fitFake.GetParError(1))
                    hTempl_offset.SetBinContent(self.etaBinningS.index(e)+1,fitFake.GetParameter(0))
                    hTempl_offset.SetBinError(self.etaBinningS.index(e)+1,fitFake.GetParError(0))
                    if self.parabolaFit :
                        hTempl_2deg.SetBinContent(self.etaBinningS.index(e)+1,fitFake.GetParameter(2))
                        hTempl_2deg.SetBinError(self.etaBinningS.index(e)+1,fitFake.GetParError(2))
                    
            h2Fakes[s] = h2Fakes_sign
            if self.kind == 'fake' or self.kind == 'fakeMC' :
                hEWSF_Fit["EWSF_chi2"+s] = hEWSF_chi2
                hEWSF_Fit["EWSF_sig"+s] = hEWSF_sig
                hEWSF_Fit["EWSF_bkg"+s] = hEWSF_bkg
            hTempl_Fit["Templ_chi2"+s] = hTempl_chi2
            hTempl_Fit["Templ_slop"+s] = hTempl_slope # the index of the dict is not a typo, it is to avoid error in continue during writing.
            hTempl_Fit["Templ_offse"+s] = hTempl_offset # the index of the dict is not a typo, it is to avoid error in continue during writing.
            if self.parabolaFit :
                hTempl_Fit["Templ_2de"+s] = hTempl_2deg # the index of the dict is not a typo, it is to avoid error in continue during writing.
            

        hFakes.update(h2Fakes)
        hFakes.update(hTempl_Fit)
        hFakes.update(hEWSF_Fit)
        
        return hFakes


    def bkg_template(self, kind, fakedict, promptdict, hdict, fit = False, tightcut = 0.15, loosecut=40, varname = 'pfRelIso04_all_corrected_MET_nom_mt', parabolaFit=False) :
        self.kind = kind
        self.fakedict = fakedict
        self.promptdict = promptdict
        self.hdict = hdict
        self.fit = fit
        self.varname = varname
        self.loosecut = loosecut
        self.tightcut = tightcut
        self.parabolaFit = parabolaFit

        kindDict = {
            'fake' : 'Data',
            'validation' : 'QCD',
            'fakeMC' : 'DataLike' ,
            'prompt' : 'WToMuNu',
            'EWKbkg' : 'EWKbkg',
        }

        datakind = kindDict[self.kind]

        htempl = {}
        h2templ = {}


        for s in self.signList :
            h2templ_sign = ROOT.TH2F("h2templ_{kind}_{sign}".format(kind=self.kind, sign=s),"h2templ_{kind}_{sign}".format(kind=self.kind, sign=s),len(self.etaBinning)-1, array('f',self.etaBinning), len(self.ptBinning)-1, array('f',self.ptBinning) )
            for e in self.etaBinningS :
                htempl_pt = ROOT.TH1F("htempl_pt_{kind}_{sign}_{eta}".format(kind=self.kind,sign=s,eta=e),"htempl_pt_{kind}_{sign}_{eta}".format(kind=self.kind,sign=s,eta=e), len(self.ptBinning)-1, array('f',self.ptBinning) )
                for p in self.ptBinningS :
                    # print "----------------------------"
                    # print "-----------------------------"
                    # print "sign, eta, pt, kind::", s, e, p, self.kind
                    htemp = hdict[p+e+s+varname+datakind+'Tot'].Clone('templ'+p+'_'+e+'_'+datakind+s+'_'+varname)

                    isoMin= htemp.GetYaxis().GetBinCenter(1)-htemp.GetYaxis().GetBinWidth(1)/2
                    binsizeTight = htemp.GetYaxis().GetBinWidth(1)
                    NcutTight=(self.tightcut-isoMin)/binsizeTight
                    NcutTight = int(NcutTight)

                    mtMin= htemp.GetXaxis().GetBinCenter(1)-htemp.GetXaxis().GetBinWidth(1)/2
                    binsizeLoose = htemp.GetXaxis().GetBinWidth(1)
                    NcutLoose=(self.loosecut-mtMin)/binsizeLoose
                    NcutLoose = int(NcutLoose)

                    tightErr = ROOT.Double(0)
                    notTightErr = ROOT.Double(0)
                    # print "cut, l, t, ", NcutLoose,NcutTight
                    Ntight = htemp.ProjectionX("htight",0,NcutTight, "e").IntegralAndError(NcutLoose-1,-1,tightErr)
                    NnotTight = htemp.ProjectionX("hNotTigth",NcutTight-1,-1, "e").IntegralAndError(NcutLoose-1,-1,notTightErr)

                    pr = 1
                    fr = 0
                    dpr =0
                    dfr =0

                    # print "Ntight,NnotTight",  Ntight, NnotTight
                    if(self.fit) :
                        fr = fakedict[s+e+'offset']+fakedict[s+e].GetBinCenter(self.ptBinningS.index(p)+1)*fakedict[s+e+'slope']
                        pr = promptdict[s+e+'offset']+promptdict[s+e].GetBinCenter(self.ptBinningS.index(p)+1)*promptdict[s+e+'slope']
                        dx_p =(promptdict[s+e].GetBinWidth(self.ptBinningS.index(p)+1))
                        dx_f = (fakedict[s+e].GetBinWidth(self.ptBinningS.index(p)+1))
                        dfr = math.sqrt(fakedict[s+e+'offsetErr']**2+(dx_f**2)*(fakedict[s+e+'slope']**2)+((fakedict[s+e].GetBinCenter(self.ptBinningS.index(p)+1))**2)*(fakedict[s+e+'slopeErr']**2))
                        dpr = math.sqrt(promptdict[s+e+'offsetErr']**2+(dx_p**2)*(promptdict[s+e+'slope']**2)+((promptdict[s+e].GetBinCenter(self.ptBinningS.index(p)+1))**2)*(promptdict[s+e+'slopeErr']**2))
                        
                        if self.parabolaFit :
                            fr = fr+(fakedict[s+e].GetBinCenter(self.ptBinningS.index(p)+1))**2*fakedict[s+e+'2deg']
                            pr = pr+(promptdict[s+e].GetBinCenter(self.ptBinningS.index(p)+1))**2*promptdict[s+e+'2deg']
                            dfr =math.sqrt(dfr**2+(4*fakedict[s+e+'2deg']**2*(fakedict[s+e].GetBinCenter(self.ptBinningS.index(p)+1))**2+4*fakedict[s+e+'2deg']*fakedict[s+e].GetBinCenter(self.ptBinningS.index(p)+1)*fakedict[s+e+'slope'])*dx_f**2)
                            dpr =math.sqrt(dpr**2+(4*promptdict[s+e+'2deg']**2*(promptdict[s+e].GetBinCenter(self.ptBinningS.index(p)+1))**2+4*promptdict[s+e+'2deg']*promptdict[s+e].GetBinCenter(self.ptBinningS.index(p)+1)*promptdict[s+e+'slope'])*dx_f**2)
                                                        

                        if pr>1 or fr>1 :
                            print "WARNING!!!!!!, pr>1 or fr>1,", pr, fr, ", with eta, pt, sign =", e, p, s
                            fr = 0
                            pr = 1
                    else :
                        fr = fakedict[s].GetBinContent(self.etaBinningS.index(e)+1, self.ptBinningS.index(p)+1)
                        pr = promptdict[s].GetBinContent(self.etaBinningS.index(e)+1, self.ptBinningS.index(p)+1)
                        dfr = fakedict[s].GetBinError(self.etaBinningS.index(e)+1, self.ptBinningS.index(p)+1)
                        dpr = promptdict[s].GetBinError(self.etaBinningS.index(e)+1, self.ptBinningS.index(p)+1)
                        if pr==0 and fr==0 :
                            print "WARNING!!!!!!, pr=fr=0, with eta, pt, sign =", e, p, s
                            fr = 0
                            pr = 1

                    # print "fake rate=", fr, "prompt rate=",pr, "(eta,pt)=",e,p, "sign=", s
                    # pr = 1
                    # print "WARNING!!!!!! pr=1"
                    scaleTight = -fr*(1-pr)/(pr-fr)
                    scaleNotTight = fr*pr/(pr-fr)
                    # print "scale tight", scaleTight, "scale not tight", scaleNotTight, fr, pr, p
                    # print "scale (tight, not tight)", scaleTight,scaleNotTight
                    scaleTightErr = math.sqrt((pr**4)*(dfr**2)-2*(pr**3)*(dfr**2)+(pr**2)*(dfr**2)+(fr**2)*((fr-1)**2)*(dpr**2))/((pr-fr)**2)
                    scaleNotTightErr = math.sqrt((pr**4)*(dfr**2)+(fr**4)*(dpr**2))/((pr-fr)**2)
                    NQCD=Ntight*scaleTight+NnotTight*scaleNotTight
                    NQCDErr = math.sqrt((scaleTight**2)*(tightErr**2)+(scaleTightErr**2)*(Ntight**2)+(scaleNotTight**2)*(notTightErr**2)+(scaleNotTightErr**2)*(NnotTight**2)   )
                    if self.kind == "prompt" or self.kind=="validation" or self.kind=='EWKbkg' :
                        # print "prompt and fake rate not applied"
                        NQCD=Ntight
                        NQCDErr = tightErr

                    htempl_pt.SetBinContent(self.ptBinningS.index(p)+1,NQCD)
                    # print "NQCD", NQCD, "kind=", self.kind,"---not tight=", NnotTight, "tight=", Ntight
                    h2templ_sign.SetBinContent(self.etaBinningS.index(e)+1, self.ptBinningS.index(p)+1,NQCD)
                    htempl_pt.SetBinError(self.ptBinningS.index(p)+1,NQCDErr)
                    # print "PROBLEMA: ERRORE NON ASSEGNATO", NQCDErr
                    ERRBIS = NQCDErr #error solving without any sense
                    h2templ_sign.SetBinError(self.etaBinningS.index(e)+1, self.ptBinningS.index(p)+1,ERRBIS)

                h2templ[s+e] = htempl_pt
            htempl[s] = h2templ_sign

        htempl.update(h2templ)
        return htempl


    def integrated_preliminary(self) :

        # print "getting histos"
        histoDict = {}
        for s in self.signList :
            for v,name in map(None,self.varList,self.varName) :
                for f,rootFile in map(None,self.sampleList,self.rootFiles) :
                    if(f!='DataLike') :
                        for r in self.regionList :
                            if(r!='Tot') :
                                # print "file", f, rootFile
                                # print "Get histo:", 'bkg_'+r+s+'/nom/bkgSel_'+v
                                # print "key=",s+v+f+r
                                histoDict[s+v+f+r] =  rootFile.Get('bkg_'+r+s+'/nom/bkgSel_'+v)
                                # print histoDict[s+v+f+r].GetName()
                                if(f!='Data') : histoDict[s+v+f+r].Scale(self.norm)
                            else:
                                # print "clone histo:", f+'_'+s+'_'+name+'_Tot'
                                histoDict[s+v+f+r] = histoDict[s+v+f+'Sideband'].Clone(f+'_'+s+'_'+name+'_Tot')
                                for rr in self.regionList :
                                    if (rr =="Sideband" or rr=="Tot") : continue
                                    # print "Added histo:", histoDict[s+v+f+rr].GetName()
                                    # print "key=",s+v+f+r
                                    # print "key added=",s+v+f+rr
                                    histoDict[s+v+f+r].Add(histoDict[s+v+f+rr])

                    else :
                        # print "cloned (datalike) histo:", 'DataLike_'+s+'_'+name+'_Tot'
                        # print "key=",s+v+f+'Tot'
                        histoDict[s+v+f+'Tot']= histoDict[s+v+'WToMuNuTot'].Clone('DataLike_'+s+'_'+name+'_Tot')
                        for ff in self.sampleList :
                            if (ff =="WToMuNu" or ff=="Data" or ff=="DataLike"):  continue
                            # print "Added (data like) histo:", histoDict[s+v+ff+'Tot'].GetName()
                            # print "key summed=",s+v+ff+'Tot'
                            histoDict[s+v+f+'Tot'].Add(histoDict[s+v+ff+'Tot'])
                # print "======================================================================="

        # print histos
        # print "ratios (integrated, preliminary)"
        ratios = []
        for s in self.signList :
            for v,name in map(None,self.varList,self.varName) :
                for f in self.sampleList :
                    # print 'ratio_'+histoDict[s+v+f+'Tot'].GetName()
                    if 'relIso' in name : cutvalue =  self.tightCut #self.relisoCUT
                    else : cutvalue = self.isoCUT
                    ratios.append(self.ratio_2Dto1D(histoDict[s+v+f+'Tot'],cutvalue,'ratio_'+histoDict[s+v+f+'Tot'].GetName()))

        output = ROOT.TFile(self.outdir+"/bkg_integrated_preliminary"+self.nameSuff+".root","recreate")
        for h in range(len(ratios)):
            ratios[h].Write()



    def differential_preliminary(self, fakerate = False) :

        self.fakerate = fakerate #if true calculate the fakerate in bin of eta, pt

        print "getting histos"
        histoDict = {}
        mtDict = {}
        PVDict = {}

        for p in self.ptBinningS :
            for e in self.etaBinningS :
                for s in self.signList :
                    for v,name in map(None,self.varList,self.varName) :
                        PVCond = bool(self.ptBinningS.index(p)==0 and self.etaBinningS.index(e)==0 and self.varList.index(v)==0)
                        MtCond = bool(self.ptBinningS.index(p)==0 and self.varList.index(v)==0)
                        for f,rootFile in map(None,self.sampleList,self.rootFiles) :
                            if(f!='DataLike') :
                                for r in self.regionList :
                                    if(r!='Tot') :
                                        # print "file", f, rootFile
                                        # print "Get histo:" 'bkg_'+r+s+'/nom/bkgSel_'+v+'_'+e+'_'+p
                                        # print "key=",p+e+s+v+f+r
                                        histoDict[p+e+s+v+f+r] =  rootFile.Get('bkg_'+r+s+'/nom/bkgSel_'+v+'_'+e+'_'+p)
                                        if(f!='Data') : histoDict[p+e+s+v+f+r].Scale(self.norm)
                                        if(MtCond) :
                                            # print "Get histo Mt:" 'bkg_'+r+s+'/nom/bkgSel_Muon_corrected_MET_nom_mt_'+e
                                            mtDict[e+s+f+r] = rootFile.Get('bkg_'+r+s+'/nom/bkgSel_Muon_corrected_MET_nom_mt_'+e)
                                            if(f!='Data') : mtDict[e+s+f+r].Scale(self.norm)
                                            # for pp in self.ptBinningS :
                                            #     mtDict[pp+e+s+f+r] =rootFile.Get('bkg_'+r+s+'/nom/bkgSel_'+v+'_'+e+'_'+pp).ProjectionX('bkgSel_Muon_corrected_MET_nom_mt_'+e+'_'+pp,0,-1,"e")
                                            #     if(f!='Data') : mtDict[pp+e+s+f+r].Scale(self.norm)
                                        if(PVCond) :
                                            PVDict[s+f+r] = rootFile.Get('bkg_'+r+s+'/nom/bkgSel_PV_npvsGood')
                                            if(f!='Data') : PVDict[s+f+r].Scale(self.norm)

                                    else:
                                        # print "clone histo:", p+e+f+'_'+s+'_'+name+'_Tot'
                                        histoDict[p+e+s+v+f+r] = histoDict[p+e+s+v+f+'Sideband'].Clone(p+'_'+e+'_'+f+'_'+s+'_'+name+'_Tot')
                                        if(MtCond) : 
                                            mtDict[e+s+f+r] = mtDict[e+s+f+'Sideband'].Clone('Mt_'+e+'_'+f+'_'+s+'_Tot') #PUT SIGNAL PER FAR FUNZIONARE
                                            # for pp in self.ptBinningS :
                                                # mtDict[pp+e+s+f+r] = mtDict[pp+e+s+f+'Sideband'].Clone('Mt_'+pp+'_'+e+'_'+f+'_'+s+'_Tot')
                                        if(PVCond) : PVDict[s+f+r] = PVDict[s+f+'Sideband'].Clone('PV_'+f+'_'+s+'_Tot')
                                        for rr in self.regionList :
                                            if (rr =="Sideband" or rr=="Tot") : continue
                                            # print "Added histo:", histoDict[p+e+s+v+f+rr].GetName()
                                            # print "key=",p+e+s+v+f+r
                                            # print "key added=",p+e+s+v+f+rr
                                            histoDict[p+e+s+v+f+r].Add(histoDict[p+e+s+v+f+rr])
                                            if(MtCond) :
                                                # print "skipped sideband adding in Mt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                                                # continue
                                                mtDict[e+s+f+r].Add(mtDict[e+s+f+rr])
                                                # for pp in self.ptBinningS :
                                                #     mtDict[pp+e+s+f+r].Add(mtDict[pp+e+s+f+rr])
                                            if(PVCond) :
                                                PVDict[s+f+r].Add(PVDict[s+f+rr])



                            else :
                                # print "cloned (datalike) histo:", 'DataLike_'+s+'_'+name+'_Tot'
                                # print "key=",p+e+s+v+f+'Tot'
                                histoDict[p+e+s+v+f+'Tot']= histoDict[p+e+s+v+'WToMuNuTot'].Clone(p+'_'+e+'_'+'DataLike_'+s+'_'+name+'_Tot')
                                if(MtCond) : 
                                    mtDict[e+s+f+'Tot'] = mtDict[e+s+'WToMuNuTot'].Clone('Mt_'+e+'_'+'DataLike_'+s+'_'+'_Tot')
                                    # for pp in self.ptBinningS :
                                        # mtDict[pp+e+s+f+'Tot'] = mtDict[pp+e+s+'WToMuNuTot'].Clone('Mt_'+pp+'_'+e+'_'+'DataLike_'+s+'_'+'_Tot')
                                if(PVCond) : PVDict[s+f+'Tot'] = PVDict[s+'WToMuNuTot'].Clone('Mt_'+'DataLike_'+s+'_'+'_Tot')
                                # print "><>>><<>><<>><<>>>>>>><<>><>>>>>>>>"
                                for ff in self.sampleList :
                                    if (ff =="WToMuNu" or ff=="Data" or ff=="DataLike"):  continue
                                    # print "Added (data like) histo:", histoDict[p+e+s+v+ff+'Tot'].GetName()
                                    # print "key summed=",p+e+s+v+ff+'Tot'
                                    histoDict[p+e+s+v+f+'Tot'].Add(histoDict[p+e+s+v+ff+'Tot'])
                                    if(MtCond) : 
                                        mtDict[e+s+f+'Tot'].Add(mtDict[e+s+ff+'Tot'])
                                        # for pp in self.ptBinningS :
                                            # mtDict[pp+e+s+f+'Tot'].Add(mtDict[pp+e+s+ff+'Tot'])
                                    if(PVDict) : PVDict[s+f+'Tot'].Add(PVDict[s+ff+'Tot'])

                                #DEBUUUG
                                errore = ROOT.Double(0)
                                binsizeTight = histoDict[p+e+s+v+f+'Tot'].GetYaxis().GetBinWidth(1)
                                NcutTight=(self.tightCut)/binsizeTight
                                NcutTight = int(NcutTight)

                                binsizeLoose = histoDict[p+e+s+v+f+'Tot'].GetXaxis().GetBinWidth(1)
                                NcutLoose=(self.looseCut)/binsizeLoose
                                NcutLoose = int(NcutLoose)
                                # print "cut, l, t, ", NcutLoose,NcutTight
                                # print "ENTRIES=", histoDict[p+e+s+v+f+'Tot'].ProjectionX("htight",NcutTight,-1, "e").IntegralAndError(NcutLoose-1,-1,errore), histoDict[p+e+s+v+f+'Tot'].ProjectionX("htight",-1,NcutTight, "e").IntegralAndError(NcutLoose-1,-1,errore)
                                #END OF DEBUG
                        # print "======================================================================="
        for p in self.ptBinningS :
            for e in self.etaBinningS :
                for s in self.signList :    
                    for v,name in map(None,self.varList,self.varName) :
                        if  self.varList.index(v)!=0 : continue
                        for f,rootFile in map(None,self.sampleList,self.rootFiles) :
                            mtDict[p+e+s+f+'Tot']= histoDict[p+e+s+v+f+'Tot'].ProjectionX('Mt_'+p+'_'+e+'_'+f+'_'+s+'_Tot',0,-1,"e")
                            # mtDict[p+e+s+f+'Tot'].SetName()



        # print histos
        # print "ratios (integrated, preliminary)"
        ratios = []
        for p in self.ptBinningS :
            for e in self.etaBinningS :
                for s in self.signList :
                    for v,name in map(None,self.varList,self.varName) :
                        for f in self.sampleList :
                            # print 'ratio_'+histoDict[p+e+s+v+f+'Tot'].GetName()
                            if 'relIso' in name : cutvalue =  self.tightCut #self.relisoCUT
                            else : cutvalue = self.isoCUT
                            # if(p+e == '30-2.4') : print "DEBUG KEY", p+e+s+v+f+'Tot'
                            ratios.append(self.ratio_2Dto1D(histoDict[p+e+s+v+f+'Tot'],cutvalue,'ratio_'+histoDict[p+e+s+v+f+'Tot'].GetName()))

        output = ROOT.TFile(self.outdir+"/bkg_differential_fakerate"+self.nameSuff+".root","recreate")
        preliminary_dir = output.mkdir("RatiosVSMt")
        preliminary_dir.cd()
        for h in range(len(ratios)):
            ratios[h].Write()
        # output.Close()
        
        histo_dir = output.mkdir("IsoVSMt")
        histo_dir.cd()
        for p in self.ptBinningS :
            for e in self.etaBinningS :
                for s in self.signList :    
                    for v,name in map(None,self.varList,self.varName) :
                        for f,rootFile in map(None,self.sampleList,self.rootFiles) :
                            histoDict[p+e+s+v+f+'Tot'].Write()
                            
        # ptVSeta_map(hdict = histoDict)
        # ABCD_hypo(hdict = histoDict)
        Mt_dir = output.mkdir("Mt")
        Mt_dir.cd()
        for e in self.etaBinningS :
            for s in self.signList :
                if(self.onData) : dataNameMt = 'Data'
                else : dataNameMt = 'DataLike'
                mtDict[e+s+'WToMuNuTot'].Write()
                mtDict[e+s+'EWKbkgTot'].Write()
                mtDict[e+s+'QCDTot'].Write()
                mtDict[e+s+dataNameMt+'Tot'].Write()
                for pp in self.ptBinningS :
                    mtDict[pp+e+s+'WToMuNuTot'].Write()
                    mtDict[pp+e+s+'EWKbkgTot'].Write()
                    mtDict[pp+e+s+'QCDTot'].Write()
                    mtDict[pp+e+s+dataNameMt+'Tot'].Write()

        # ABCD_dir = output.mkdir("ABCD_checks")
        # ABCD_dir.cd()

        hIsoMC = self.isolationAna(kind=self.dataOpt,hdict=histoDict,varname = self.varFake, loosecut = self.looseCut)
        hIsoValidation = self.isolationAna(kind = 'validation',hdict=histoDict,varname = self.varFake, loosecut = self.looseCut)
        hIsoPrompt = self.isolationAna(kind = 'prompt',hdict=histoDict,varname = self.varFake, loosecut = self.looseCut)
        hIsoEWKbkg = self.isolationAna(kind = 'EWKbkg',hdict=histoDict,varname = self.varFake, loosecut = self.looseCut)
        fakerate_dir = output.mkdir("Isolation")
        fakerate_dir.cd()
        for a,b,c,d in map(None,hIsoMC,hIsoValidation,hIsoPrompt,hIsoEWKbkg):
                hIsoMC[a].Write()
                hIsoValidation[b].Write()
                hIsoPrompt[c].Write()
                hIsoEWKbkg[d].Write()

        PV_dir = output.mkdir("PV")
        PV_dir.cd()
        for s in self.signList :
                if(self.onData) : dataNameMt = 'Data'
                else : dataNameMt = 'DataLike'
                PVDict[s+'WToMuNuTot'].Write()
                PVDict[s+'EWKbkgTot'].Write()
                PVDict[s+'QCDTot'].Write()
                PVDict[s+dataNameMt+'Tot'].Write()

        if(self.fakerate) :
            print "evaluating fakerates"
            hfakesMC = self.differential_fakerate(kind = self.dataOpt, hdict = histoDict, mtDict = mtDict, varname = self.varFake, tightcut = self.tightCut, loosecut = self.looseCut, EWSFfit=self.EWSFfit, highMtCut=60,parabolaFit = self.parabolaFit)
            hprompt =self.differential_fakerate(kind = 'prompt', hdict = histoDict, mtDict = mtDict, varname = self.varFake, tightcut = self.tightCut, loosecut = self.looseCut,parabolaFit = self.parabolaFit)
            hvalidation =self.differential_fakerate(kind = 'validation', hdict = histoDict, mtDict = mtDict, varname = self.varFake, tightcut = self.tightCut, loosecut = self.looseCut,parabolaFit = self.parabolaFit)
            hvalidationSigReg =self.differential_fakerate(kind = 'validationSigReg', hdict = histoDict, mtDict = mtDict, varname = self.varFake, tightcut = self.tightCut, loosecut = self.looseCut,parabolaFit = self.parabolaFit)
            hpromptSideband =self.differential_fakerate(kind = 'promptSideband', hdict = histoDict, mtDict = mtDict, varname = self.varFake, tightcut = self.tightCut, loosecut = self.looseCut,parabolaFit = self.parabolaFit)

            # outputFake = ROOT.TFile(self.outdir+"/bkg_differential_fakerate.root","recreate")
            fakerate_dir = output.mkdir("Fakerate")
            fakerate_dir.cd()

            print ("writing fakerates")
            for b,c,d,e in map(None,hprompt,hvalidation,hvalidationSigReg,hpromptSideband):
                if 'offset' in (b or c or d or e) or 'slope' in (b or c or d or e) or '2deg' in (b or c or d or e): continue
                # if 'offset' in (b or c or d or e) or 'slope' in (b or c or d or e) : continue
                hprompt[b].Write()
                hvalidation[c].Write()
                hvalidationSigReg[d].Write()
                hpromptSideband[e].Write()
            for a in hfakesMC:
                if 'offset' in a or 'slope' in a  or '2deg' in a : continue
                # if 'offset' in a or 'slope' in a: continue
                hfakesMC[a].Write()

            template_dir = output.mkdir("Template")
            template_dir.cd()
            print "evaluating templates"
            bkg_templateMC = self.bkg_template(kind = self.dataOpt, fakedict=hfakesMC, promptdict=hprompt, hdict =histoDict, fit =self.fitOnTemplate, tightcut = self.tightCut, loosecut=self.looseCut, varname = self.varFake, parabolaFit = self.parabolaFit)
            bkg_templatePrompt = self.bkg_template(kind = 'prompt', fakedict=hfakesMC, promptdict=hprompt, hdict =histoDict, fit =self.fitOnTemplate, tightcut = self.tightCut, loosecut=self.looseCut, varname = self.varFake, parabolaFit = self.parabolaFit)
            bkg_templateValidation = self.bkg_template(kind = 'validation', fakedict=hfakesMC, promptdict=hprompt, hdict =histoDict, fit =self.fitOnTemplate, tightcut = self.tightCut, loosecut=self.looseCut, varname = self.varFake, parabolaFit = self.parabolaFit)
            bkg_templateEWKbkg = self.bkg_template(kind = 'EWKbkg', fakedict=hfakesMC, promptdict=hprompt, hdict =histoDict, fit =self.fitOnTemplate, tightcut = self.tightCut, loosecut=self.looseCut, varname = self.varFake, parabolaFit = self.parabolaFit)

            print "Writing templates"
            for a,b,c in map(None, bkg_templateMC,bkg_templatePrompt,bkg_templateValidation) :
                bkg_templateMC[a].Write()
                bkg_templatePrompt[b].Write()
                bkg_templateValidation[c].Write()

        output.Close()

    def fakerate_plots(self, variations=False,tightCutList=[0.15],looseCutList=[40], parabolaFit = False) :
        print "plotting fakerate"

        self.variations = variations
        self.tightCutList= tightCutList
        self.looseCutList= looseCutList
        self.parabolaFit = parabolaFit

        inputFile = ROOT.TFile.Open(self.outdir+"/bkg_differential_fakerate"+self.nameSuff+".root")

        canvasList = []
        legDict = {}
        stackDict ={}

        for s in self.signList :
            #---------------------------------------------PV PLOTS ---------------------------------------------#
            # print "PV plots"

            c_PV = ROOT.TCanvas("c_PV_{sign}".format(sign=s),"c_PV_{sign}".format(sign=s),800,600)
            c_PV.cd()
            c_PV.SetGridx()
            c_PV.SetGridy()

            if(self.onData) : dataNameMt = 'Data'
            else : dataNameMt = 'DataLike'

            h_PV_data = inputFile.Get("PV/PV_"+dataNameMt+"_"+s+"_Tot")
            h_PV_bkg = inputFile.Get("PV/PV_QCD_"+s+"_Tot")
            h_PV_sig = inputFile.Get("PV/PV_WToMuNu_"+s+"_Tot")
            h_PV_EWKbkg = inputFile.Get("PV/PV_EWKbkg_"+s+"_Tot")
            h_PV_sig.Add(h_PV_EWKbkg)

            # h_PV_sig.Rebin(3)
            # h_PV_bkg.Rebin(3)
            # h_PV_data.Rebin(3)

            h_PV_data.SetLineWidth(3)
            h_PV_bkg.SetLineWidth(3)
            h_PV_sig.SetLineWidth(3)
            # h_fake.SetLineWidth(3)

            # h_PV_data.SetLineColor(632+2) #red
            h_PV_data.SetLineColor(1) #black
            h_PV_bkg.SetLineColor(600-4) #blue
            h_PV_sig.SetLineColor(416+2) #green
            h_PV_bkg.SetFillColor(600-4) #blue
            h_PV_sig.SetFillColor(416+2) #green
            # h_fake.SetLineColor(1) #black

            h_PV_data.Sumw2()
            h_PV_data.SetMarkerStyle(20)
            h_PV_data.Draw()


            stackDict[s+"PV"] = ROOT.THStack("PV"+s,"")
            stackDict[s+"PV"].Add(h_PV_sig)
            stackDict[s+"PV"].Add(h_PV_bkg)
            stackDict[s+"PV"].Draw("SAME HIST")
            h_PV_data.DrawCopy("SAME")
            # h_PV_bkg.Draw("SAME")
            # h_PV_sig.Draw("SAME")
            # h_fake.Draw("SAME")

            h_PV_data.GetYaxis().SetTitle("dN/dPV/{size} [1/GeV]".format(size=h_PV_data.GetBinWidth(1)))
            h_PV_data.GetYaxis().SetTitleOffset(1)
            h_PV_data.GetXaxis().SetTitle("N PV good")
            # h_PV_data.SetTitle("Fake Rates, {min}<#eta<{max}, W {sign}".format(min=self.etaBinning[self.etaBinning.index(float(e))-1], max=e, sign=s))
            h_PV_data.SetTitle("Number of Primary vertices, W {sign}".format(sign=s))

            legDict[s+"PV"] = ROOT.TLegend(0.1,0.7,0.48,0.9)
            legDict[s+"PV"].AddEntry(h_PV_data,"Data")
            legDict[s+"PV"].AddEntry(h_PV_bkg, "QCD MC")
            legDict[s+"PV"].AddEntry(h_PV_sig, "EWK MC")
            # legDict[e+s].AddEntry(h_fake, "Data")
            legDict[s+"PV"].Draw("SAME")

            canvasList.append(c_PV)
            # c_comparison.SaveAs(self.outdir+"/plot/"+c_Mt_EWSF.GetName()+'.png')

            for e in self.etaBinningS :

                    #---------------------------------------------COMPARISON PLOTS ---------------------------------------------#

                    # print "comparison plots"

                    c_comparison = ROOT.TCanvas("c_comparison_{sign}_{eta}".format(sign=s,eta=e),"c_comparison_{sign}_{eta}".format(sign=s,eta=e),800,600)
                    c_comparison.cd()
                    c_comparison.SetGridx()
                    c_comparison.SetGridy()


                    h_fakeMC = inputFile.Get("Fakerate/hFakes_pt_"+self.dataOpt+"_"+s+"_"+e)
                    h_prompt = inputFile.Get("Fakerate/hFakes_pt_prompt_"+s+"_"+e)
                    h_validation = inputFile.Get("Fakerate/hFakes_pt_validation_"+s+"_"+e)

                    # h_fake = inputFile.Get("Fakerate/hFakes_pt_fake_"+s+"_"+e)

                    h_fakeMC.SetLineWidth(3)
                    h_prompt.SetLineWidth(3)
                    h_validation.SetLineWidth(3)
                    # h_fake.SetLineWidth(3)

                    h_fakeMC.SetLineColor(632+2) #red
                    h_prompt.SetLineColor(600-4) #blue
                    h_validation.SetLineColor(416+2) #green
                    # h_fake.SetLineColor(1) #black

                    h_fakeMC.Draw()
                    h_prompt.Draw("SAME")
                    h_validation.Draw("SAME")
                    # h_fake.Draw("SAME")

                    h_fakeMC.GetYaxis().SetRangeUser(0,1.1)
                    h_fakeMC.GetYaxis().SetTitle("Isolation Cut Efficiency")
                    h_fakeMC.GetYaxis().SetTitleOffset(1)
                    h_fakeMC.GetXaxis().SetTitle("p_{T}^{#mu} [GeV]")
                    # h_fakeMC.SetTitle("Fake Rates, {min}<#eta<{max}, W {sign}".format(min=self.etaBinning[self.etaBinning.index(float(e))-1], max=e, sign=s))
                    h_fakeMC.SetTitle("Fake Rates, {min}<#eta<{max}, W {sign}".format(min=e, max=self.etaBinning[self.etaBinning.index(float(e))+1], sign=s))

                    legDict[e+s] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                    legDict[e+s].AddEntry(h_fakeMC,"Data")
                    legDict[e+s].AddEntry(h_prompt, "W MC")
                    legDict[e+s].AddEntry(h_validation, "QCD MC")
                    # legDict[e+s].AddEntry(h_fake, "Data")
                    legDict[e+s].Draw("SAME")

                    canvasList.append(c_comparison)
                    # c_comparison.SaveAs(self.outdir+"/plot/"+c_comparison.GetName()+'.png')

                    #---------------------------------------------ABCD cehcks PLOTS ---------------------------------------------#
                    # print "ABCD check plots"

                    c_ABCDcheck = ROOT.TCanvas("c_ABCDcheck_{sign}_{eta}".format(sign=s,eta=e),"c_ABCDcheck_{sign}_{eta}".format(sign=s,eta=e),800,600)
                    c_ABCDcheck.cd()
                    c_ABCDcheck.SetGridx()
                    c_ABCDcheck.SetGridy()


                    h_prompt = inputFile.Get("Fakerate/hFakes_pt_prompt_"+s+"_"+e)
                    h_validation = inputFile.Get("Fakerate/hFakes_pt_validation_"+s+"_"+e)
                    h_promptSideband = inputFile.Get("Fakerate/hFakes_pt_promptSideband_"+s+"_"+e)
                    h_validationSigReg = inputFile.Get("Fakerate/hFakes_pt_validationSigReg_"+s+"_"+e)

                    # h_fake = inputFile.Get("Fakerate/hFakes_pt_fake_"+s+"_"+e)

                    h_promptSideband.SetLineWidth(3)
                    h_prompt.SetLineWidth(3)
                    h_validation.SetLineWidth(3)
                    h_validationSigReg.SetLineWidth(3)

                    h_promptSideband.SetLineColor(632+2) #red
                    h_prompt.SetLineColor(600-4) #blue
                    h_validation.SetLineColor(416+2) #green
                    h_validationSigReg.SetLineColor(1) #black

                    h_promptSideband.Draw()
                    h_prompt.Draw("SAME")
                    h_validation.Draw("SAME")
                    h_validationSigReg.Draw("SAME")

                    h_promptSideband.GetYaxis().SetRangeUser(0,1.1)
                    h_promptSideband.GetYaxis().SetTitle("Isolation Cut Efficiency")
                    h_promptSideband.GetYaxis().SetTitleOffset(1)
                    h_promptSideband.GetXaxis().SetTitle("p_{T}^{#mu} [GeV]")
                    h_promptSideband.SetTitle("Fake Rates, {min}<#eta<{max}, ABCD check, W {sign}".format(min=e, max=self.etaBinning[self.etaBinning.index(float(e))+1], sign=s))

                    legDict[e+s+'ABCD'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                    legDict[e+s+'ABCD'].AddEntry(h_promptSideband,"W, Sideband")
                    legDict[e+s+'ABCD'].AddEntry(h_prompt, "W Signal Region")
                    legDict[e+s+'ABCD'].AddEntry(h_validation, "QCD Sideband")
                    legDict[e+s+'ABCD'].AddEntry(h_validationSigReg, "QCD Signal Region")
                    legDict[e+s+'ABCD'].Draw("SAME")

                    canvasList.append(c_ABCDcheck)
                    # c_comparison.SaveAs(self.outdir+"/plot/"+c_ABCDcheck.GetName()+'.png')


                    #---------------------------------------------TEMPLATE PLOTS ---------------------------------------------#
                    # print "template plots"

                    c_template = ROOT.TCanvas("c_template_{sign}_{eta}".format(sign=s,eta=e),"c_template_{sign}_{eta}".format(sign=s,eta=e),800,600)
                    c_template.cd()
                    c_template.SetGridx()
                    c_template.SetGridy()


                    h_template = inputFile.Get("Template/htempl_pt_"+self.dataOpt+"_"+s+"_"+e)
                    h_W = inputFile.Get("Template/htempl_pt_prompt_"+s+"_"+e)
                    h_qcd = inputFile.Get("Template/htempl_pt_validation_"+s+"_"+e)

                    # h_fake = inputFile.Get("Fakerate/hFakes_pt_fake_"+s+"_"+e)

                    h_template.SetLineWidth(3)
                    h_W.SetLineWidth(3)
                    h_qcd.SetLineWidth(3)
                    # h_fake.SetLineWidth(3)

                    h_template.SetLineColor(632+2) #red
                    h_W.SetLineColor(600-4) #blue
                    h_qcd.SetLineColor(416+2) #green
                    # h_fake.SetLineColor(1) #black

                    h_W.Draw()
                    h_template.Draw("SAME")
                    h_qcd.Draw("SAME")
                    # h_fake.Draw("SAME")

                    h_W.GetYaxis().SetRangeUser(0,1500000)
                    h_W.GetYaxis().SetTitle("Events")
                    h_W.GetYaxis().SetTitleOffset(1)
                    h_W.GetXaxis().SetTitle("p_{T}^{#mu} [GeV]")
                    h_W.SetTitle("Compared yelds, {min}<#eta<{max}, W {sign}".format(min=e, max=self.etaBinning[self.etaBinning.index(float(e))+1], sign=s))

                    legDict[e+s+"templ"] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                    legDict[e+s+"templ"].AddEntry(h_template,"Template bkg MC")
                    legDict[e+s+"templ"].AddEntry(h_W, "W MC")
                    legDict[e+s+"templ"].AddEntry(h_qcd, "QCD MC")
                    # legDict[e+s].AddEntry(h_fake, "Data")
                    legDict[e+s].Draw("SAME")

                    canvasList.append(c_template)
                    # c_comparison.SaveAs(self.outdir+"/plot/"+c_template.GetName()+'.png')


                    #---------------------------------------------Mt EWSF PLOTS ---------------------------------------------#
                    # print "Mt EWSF plots"

                    c_Mt_EWSF = ROOT.TCanvas("c_Mt_EWSF_{sign}_{eta}".format(sign=s,eta=e),"c_Mt_EWSF_{sign}_{eta}".format(sign=s,eta=e),800,600)
                    c_Mt_EWSF.cd()
                    c_Mt_EWSF.SetGridx()
                    c_Mt_EWSF.SetGridy()

                    if(self.onData) : dataNameMt = 'Data'
                    else : dataNameMt = 'DataLike'

                    h_Mt_data = inputFile.Get("Mt/Mt_"+e+"_"+dataNameMt+"_"+s+"_Tot")
                    h_Mt_bkg = inputFile.Get("Mt/Mt_"+e+"_QCD_"+s+"_Tot")
                    h_Mt_sig = inputFile.Get("Mt/Mt_"+e+"_WToMuNu_"+s+"_Tot")
                    h_Mt_EWKbkg = inputFile.Get("Mt/Mt_"+e+"_EWKbkg_"+s+"_Tot")
                    h_Mt_sig.Add(h_Mt_EWKbkg)

                    # h_Mt_sig.Rebin(3)
                    # h_Mt_bkg.Rebin(3)
                    # h_Mt_data.Rebin(3)

                    h_Mt_data.SetLineWidth(3)
                    h_Mt_bkg.SetLineWidth(3)
                    h_Mt_sig.SetLineWidth(3)
                    # h_fake.SetLineWidth(3)

                    # h_Mt_data.SetLineColor(632+2) #red
                    h_Mt_data.SetLineColor(1) #black
                    h_Mt_bkg.SetLineColor(600-4) #blue
                    h_Mt_sig.SetLineColor(416+2) #green
                    h_Mt_bkg.SetFillColor(600-4) #blue
                    h_Mt_sig.SetFillColor(416+2) #green
                    # h_fake.SetLineColor(1) #black

                    # h_Mt_data.Sumw2()
                    h_Mt_data.SetMarkerStyle(20)
                    h_Mt_data.Draw()


                    stackDict[e+s+"Mt"] = ROOT.THStack("Mt"+e+s,"")
                    stackDict[e+s+"Mt"].Add(h_Mt_sig)
                    stackDict[e+s+"Mt"].Add(h_Mt_bkg)
                    stackDict[e+s+"Mt"].Draw("SAME HIST")
                    h_Mt_data.DrawCopy("SAME")
                    # h_Mt_bkg.Draw("SAME")
                    # h_Mt_sig.Draw("SAME")
                    # h_fake.Draw("SAME")

                    h_Mt_data.GetYaxis().SetTitle("dN/dMt/{size} [1/GeV]".format(size=h_Mt_data.GetBinWidth(1)))
                    h_Mt_data.GetYaxis().SetTitleOffset(1)
                    h_Mt_data.GetXaxis().SetTitle("M_{T} [GeV]")
                    # h_Mt_data.SetTitle("Fake Rates, {min}<#eta<{max}, W {sign}".format(min=self.etaBinning[self.etaBinning.index(float(e))-1], max=e, sign=s))
                    h_Mt_data.SetTitle("Transverse Mass, {min}<#eta<{max}, W {sign}".format(min=e, max=self.etaBinning[self.etaBinning.index(float(e))+1], sign=s))

                    legDict[e+s+"Mt"] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                    legDict[e+s+"Mt"].AddEntry(h_Mt_data,"Data")
                    legDict[e+s+"Mt"].AddEntry(h_Mt_bkg, "QCD MC")
                    legDict[e+s+"Mt"].AddEntry(h_Mt_sig, "EWK MC")
                    # legDict[e+s].AddEntry(h_fake, "Data")
                    legDict[e+s+"Mt"].Draw("SAME")

                    canvasList.append(c_Mt_EWSF)
                    # c_comparison.SaveAs(self.outdir+"/plot/"+c_Mt_EWSF.GetName()+'.png')

                    notIsoPlots = False
                    for p in self.ptBinningS :
                        if(notIsoPlots) : continue
                        #---------------------------------------------ISO PLOTS ---------------------------------------------#
                        # print "ISO plots"

                        c_Iso = ROOT.TCanvas("c_Iso_{sign}_{eta}_{pt}".format(sign=s,eta=e, pt=p),"c_Iso_{sign}_{eta}_{pt}".format(sign=s,eta=e, pt=p),800,600)
                        c_Iso.cd()
                        c_Iso.SetGridx()
                        c_Iso.SetGridy()

                        if(self.onData) : dataNameIso = 'fake'
                        else : dataNameIso = 'FakeMC'

                        h_iso_data = inputFile.Get("Isolation/Iso_"+dataNameIso+"_"+s+"_"+e+"_"+p)
                        h_iso_bkg = inputFile.Get("Isolation/Iso_validation_"+s+"_"+e+"_"+p)
                        h_iso_sig = inputFile.Get("Isolation/Iso_prompt_"+s+"_"+e+"_"+p)
                        h_iso_EWKbkg = inputFile.Get("Isolation/Iso_EWKbkg_"+s+"_"+e+"_"+p)
                        # print "Isolation/iso_EWKbkg_"+s+"_"+e+"_"+p
                        # print "Isolation/iso_prompt_"+s+"_"+e+"_"+p
                        h_iso_sig.Add(h_iso_EWKbkg)

                        # h_iso_sig.Rebin(20)
                        # h_iso_bkg.Rebin(20)
                        # h_iso_data.Rebin(20)

                        h_iso_data.SetLineWidth(3)
                        h_iso_bkg.SetLineWidth(3)
                        h_iso_sig.SetLineWidth(3)
                        # h_fake.SetLineWidth(3)

                        # h_iso_data.SetLineColor(632+2) #red
                        h_iso_data.SetLineColor(1) #black
                        h_iso_bkg.SetLineColor(600-4) #blue
                        h_iso_sig.SetLineColor(416+2) #green
                        h_iso_bkg.SetFillColor(600-4) #blue
                        h_iso_sig.SetFillColor(416+2) #green
                        # h_fake.SetLineColor(1) #black

                        h_iso_data.Sumw2()
                        h_iso_data.SetMarkerStyle(20)
                        h_iso_data.Draw()

                        stackDict[e+s+p+"Iso"] = ROOT.THStack("Iso"+e+s+p,"")
                        stackDict[e+s+p+"Iso"].Add(h_iso_sig)
                        stackDict[e+s+p+"Iso"].Add(h_iso_bkg)
                        stackDict[e+s+p+"Iso"].Draw("SAME HIST")
                        h_iso_data.DrawCopy("SAME")
                        # h_iso_bkg.Draw("SAME")
                        # h_iso_sig.Draw("SAME")
                        # h_fake.Draw("SAME")

                        h_iso_data.GetYaxis().SetTitle("dN/dIso/{size} [1/GeV]".format(size=h_iso_data.GetBinWidth(1)))
                        h_iso_data.GetYaxis().SetTitleOffset(1)
                        h_iso_data.GetXaxis().SetTitle("RelIso_{T} ")
                        # h_iso_data.SetTitle("Fake Rates, {min}<#eta<{max}, W {sign}".format(min=self.etaBinning[self.etaBinning.index(float(e))-1], max=e, sign=s))
                        h_iso_data.SetTitle("Relaitve Isolation 04, {min}<#eta<{max},{pmin}<Pt<{pmax}, W {sign}".format(min=e, max=self.etaBinning[self.etaBinning.index(float(e))+1],pmin=p, pmax=self.ptBinning[self.ptBinning.index(float(p))+1], sign=s))

                        legDict[e+s+p+"Iso"] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                        legDict[e+s+p+"Iso"].AddEntry(h_iso_data,"Data")
                        legDict[e+s+p+"Iso"].AddEntry(h_iso_bkg, "QCD MC")
                        legDict[e+s+p+"Iso"].AddEntry(h_iso_sig, "EWK MC")
                        # legDict[e+s].AddEntry(h_fake, "Data")
                        legDict[e+s+p+"Iso"].Draw("SAME")

                        canvasList.append(c_Iso)
                        # c_comparison.SaveAs(self.outdir+"/plot/"+c_Iso.GetName()+'.png')
                        
                        
                        #---------------------------------------------Mt EWSF PLOTS (PT dependednt)---------------------------------------------#
                        # print "Mt EWSF plots"

                        c_Mt_EWSF = ROOT.TCanvas("c_Mt_EWSF_{sign}_{eta}_{pt}".format(sign=s,eta=e,pt=p),"c_Mt_EWSF_{sign}_{eta}_{pt}".format(sign=s,eta=e,pt=p),800,600)
                        c_Mt_EWSF.cd()
                        c_Mt_EWSF.SetGridx()
                        c_Mt_EWSF.SetGridy()

                        if(self.onData) : dataNameMt = 'Data'
                        else : dataNameMt = 'DataLike'

                        h_Mt_data = inputFile.Get("Mt/Mt_"+p+"_"+e+"_"+dataNameMt+"_"+s+"_Tot")
                        h_Mt_bkg = inputFile.Get("Mt/Mt_"+p+"_"+e+"_QCD_"+s+"_Tot")
                        h_Mt_sig = inputFile.Get("Mt/Mt_"+p+"_"+e+"_WToMuNu_"+s+"_Tot")
                        h_Mt_EWKbkg = inputFile.Get("Mt/Mt_"+p+"_"+e+"_EWKbkg_"+s+"_Tot")
                        h_Mt_sig.Add(h_Mt_EWKbkg)

                        # h_Mt_sig.Rebin(3)
                        # h_Mt_bkg.Rebin(3)
                        # h_Mt_data.Rebin(3)

                        h_Mt_data.SetLineWidth(3)
                        h_Mt_bkg.SetLineWidth(3)
                        h_Mt_sig.SetLineWidth(3)
                        # h_fake.SetLineWidth(3)

                        # h_Mt_data.SetLineColor(632+2) #red
                        h_Mt_data.SetLineColor(1) #black
                        h_Mt_bkg.SetLineColor(600-4) #blue
                        h_Mt_sig.SetLineColor(416+2) #green
                        h_Mt_bkg.SetFillColor(600-4) #blue
                        h_Mt_sig.SetFillColor(416+2) #green
                        # h_fake.SetLineColor(1) #black

                        # h_Mt_data.Sumw2()
                        h_Mt_data.SetMarkerStyle(20)
                        h_Mt_data.Draw()


                        stackDict[p+e+s+"Mt"] = ROOT.THStack("Mt"+e+s,"")
                        stackDict[p+e+s+"Mt"].Add(h_Mt_sig)
                        stackDict[p+e+s+"Mt"].Add(h_Mt_bkg)
                        stackDict[p+e+s+"Mt"].Draw("SAME HIST")
                        h_Mt_data.DrawCopy("SAME")
                        # h_Mt_bkg.Draw("SAME")
                        # h_Mt_sig.Draw("SAME")
                        # h_fake.Draw("SAME")

                        h_Mt_data.GetYaxis().SetTitle("dN/dMt/{size} [1/GeV]".format(size=h_Mt_data.GetBinWidth(1)))
                        h_Mt_data.GetYaxis().SetTitleOffset(1)
                        h_Mt_data.GetXaxis().SetTitle("M_{T} [GeV]")
                        # h_Mt_data.SetTitle("Fake Rates, {min}<#eta<{max}, W {sign}".format(min=self.etaBinning[self.etaBinning.index(float(e))-1], max=e, sign=s))
                        h_Mt_data.SetTitle("Transverse Mass, {min}<#eta<{max}, {pmin}<Pt<{pmax}, W {sign}".format(min=e, max=self.etaBinning[self.etaBinning.index(float(e))+1], sign=s,pmin=p, pmax=self.ptBinning[self.ptBinning.index(float(p))+1]))

                        legDict[p+e+s+"Mt"] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                        legDict[p+e+s+"Mt"].AddEntry(h_Mt_data,"Data")
                        legDict[p+e+s+"Mt"].AddEntry(h_Mt_bkg, "QCD MC")
                        legDict[p+e+s+"Mt"].AddEntry(h_Mt_sig, "EWK MC")
                        # legDict[e+s].AddEntry(h_fake, "Data")
                        legDict[p+e+s+"Mt"].Draw("SAME")

                        canvasList.append(c_Mt_EWSF)
                        # c_comparison.SaveAs(self.outdir+"/plot/"+c_Mt_EWSF.GetName()+'.png')                        
                        
                        
                        
                        




        #---------------------------------------------FIT CHECK PLOTS ---------------------------------------------#
        # print "Fit check plots"
        
        if not self.parabolaFit :
            typeFitDict = {
                'EWSF' : ['chi2', 'bkg', 'sig'],
                'Templ' : ['chi2','slope', 'offset'],
            }
        else :
            typeFitDict = {
                'EWSF' : ['chi2', 'bkg', 'sig' ],
                'Templ' : ['chi2','slope', 'offset', '2deg'],
            }            
        
        

        c_fitDict = {}
        hFitDict_Plus = {}
        hFitDict_Minus = {}


        for ty in typeFitDict :
            if ty == "Templ" : kind = [self.dataOpt, 'prompt']
            else : kind = [self.dataOpt]
            for var in typeFitDict[ty] :

                for ki in kind :
                    c_fitDict[ty+var+ki] = ROOT.TCanvas("c_"+ty+"_"+var+'_'+ki,"c_"+ty+"_"+var+'_'+ki,800,600)
                    c_fitDict[ty+var+ki].cd()
                    c_fitDict[ty+var+ki].SetGridx()
                    c_fitDict[ty+var+ki].SetGridy()
                    # print "Fakerate/h"+ty+"_"+var+"_"+ki+"_Plus"
                    hFitDict_Plus[ty+var+ki] = inputFile.Get("Fakerate/h"+ty+"_"+var+"_"+ki+"_Plus")
                    hFitDict_Minus[ty+var+ki] = inputFile.Get("Fakerate/h"+ty+"_"+var+"_"+ki+"_Minus")
                    # hEWSF_chi2_fake_Plus = inputFile.Get("Fakerate/hEWSF_chi2_fake_Plus")
                    # hEWSF_chi2_fake_Minus = inputFile.Get("Fakerate/hEWSF_chi2_fake_Minus")
                    hFitDict_Plus[ty+var+ki].SetLineWidth(3)
                    hFitDict_Minus[ty+var+ki].SetLineWidth(3)
                    hFitDict_Plus[ty+var+ki].SetLineColor(632+2) #red
                    hFitDict_Minus[ty+var+ki].SetLineColor(600-4) #blue
                    hFitDict_Plus[ty+var+ki].Draw()
                    hFitDict_Minus[ty+var+ki].Draw("SAME")
                    hFitDict_Plus[ty+var+ki].GetYaxis().SetTitleOffset(1)
                    hFitDict_Plus[ty+var+ki].GetXaxis().SetTitle("#eta^{#mu}")

                    if(var=="chi2") :
                        hFitDict_Plus[ty+var+ki].GetYaxis().SetTitle("Reduced #chi^{2}")
                    else :
                        hFitDict_Plus[ty+var+ki].GetYaxis().SetTitle("Par. Value")
                    hFitDict_Plus[ty+var+ki].SetTitle("{FitKind} Fit, {var}, {kind} ".format(FitKind=ty,var=var,kind=ki))
                    legDict[ty+var+ki] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                    legDict[ty+var+ki].AddEntry(hFitDict_Plus[ty+var+ki],"W plus")
                    legDict[ty+var+ki].AddEntry(hFitDict_Minus[ty+var+ki], "W minus")
                    legDict[ty+var+ki].Draw("SAME")
                    # canvasList.append(c_EWSF_chi2)
                    canvasList.append(c_fitDict[ty+var+ki])
                    # c_comparison.SaveAs(self.outdir+"/plot/"+c_fitDict[ty+var+ki].GetName()+'.png')
        ty='EWSF'
        ki = self.dataOpt
        for var in ['chi2', 'bkg', 'sig', 'count'] :
            for e in self.etaBinningS :
                            c_fitDict[ty+var+ki+e] = ROOT.TCanvas("c_"+ty+"_"+var+'_'+ki+'_'+e,"c_"+ty+"_"+var+'_'+ki+'_'+e,800,600)
                            c_fitDict[ty+var+ki+e].cd()
                            c_fitDict[ty+var+ki+e].SetGridx()
                            c_fitDict[ty+var+ki+e].SetGridy()
                            hFitDict_Plus[ty+var+ki+e] = inputFile.Get("Fakerate/h"+ty+"_"+var+"_"+ki+"_Plus_"+e)
                            hFitDict_Minus[ty+var+ki+e] = inputFile.Get("Fakerate/h"+ty+"_"+var+"_"+ki+"_Minus_"+e)
                            hFitDict_Plus[ty+var+ki+e].SetLineWidth(3)
                            hFitDict_Minus[ty+var+ki+e].SetLineWidth(3)
                            hFitDict_Plus[ty+var+ki+e].SetLineColor(632+2) #red
                            hFitDict_Minus[ty+var+ki+e].SetLineColor(600-4) #blue
                            hFitDict_Plus[ty+var+ki+e].Draw()
                            hFitDict_Minus[ty+var+ki+e].Draw("SAME")
                            hFitDict_Plus[ty+var+ki+e].GetYaxis().SetTitleOffset(1)
                            hFitDict_Plus[ty+var+ki+e].GetXaxis().SetTitle("P_{T}^{#mu}")

                            if(var=="chi2") :
                                hFitDict_Plus[ty+var+ki+e].GetYaxis().SetTitle("Reduced #chi^{2}")
                            else :
                                hFitDict_Plus[ty+var+ki+e].GetYaxis().SetTitle("Par. Value")
                            hFitDict_Plus[ty+var+ki+e].SetTitle("{FitKind} Fit, {var}, {kind} ".format(FitKind=ty,var=var,kind=ki))
                            legDict[ty+var+ki+e] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                            legDict[ty+var+ki+e].AddEntry(hFitDict_Plus[ty+var+ki+e],"W plus")
                            legDict[ty+var+ki+e].AddEntry(hFitDict_Minus[ty+var+ki+e], "W minus")
                            legDict[ty+var+ki+e].Draw("SAME")
                            # canvasList.append(c_EWSF_chi2)
                            canvasList.append(c_fitDict[ty+var+ki+e])                            

                    
                    


        #---------------------------------------------VARIATION OF CUTS PLOTS ---------------------------------------------#

        if(self.variations) :
            # print "Cut Variation plots"

            fileVarDict = {}
            graphVarDict = {}

            for lcut in self.looseCutList :
                for tcut in self.tightCutList :
                    fileVarDict[str(lcut)+str(tcut)] = ROOT.TFile.Open(self.outdir+"/bkg_differential_fakerate"+"_"+str(lcut)+"_"+str(tcut)+".root")

            for s in self.signList :
                for e in self.etaBinningS :
                    ptPoint = 1
                    for p in self.ptBinningS :
                        c_variation = ROOT.TCanvas("c_variations_{sign}_{eta}_{pt}".format(sign=s,eta=e,pt=p),"c_variations_{sign}_{eta}_{pt}".format(sign=s,eta=e, pt=p),800,600)
                        c_variation.cd()
                        legDict[e+s+p+"Variation"] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                        tightPoint =0

                        #line of W QCD value
                        h_templW = fileVarDict[str(lcut)+str(tcut)].Get("Template/htempl_pt_prompt_"+s+"_"+e)
                        Wvalue = ROOT.TLine(4,h_templW.GetBinContent(ptPoint),86,h_templW.GetBinContent(ptPoint))
                        Wvalue.SetLineColor(880)
                        Wvalue.SetLineStyle(2)
                        Wvalue.SetLineWidth(3)
                        # print "W VALUE", h_templW.GetBinContent(ptPoint), "sep", s,e, p

                        for tcut in self.tightCutList :
                            graphVarDict[s+e+p+str(tcut)] = ROOT.TGraphErrors()
                            loosePoint = 0
                            for lcut in self.looseCutList :
                                h_template = fileVarDict[str(lcut)+str(tcut)].Get("Template/htempl_pt_"+self.dataOpt+"_"+s+"_"+e)
                                graphVarDict[s+e+p+str(tcut)].SetPoint(loosePoint,lcut,h_template.GetBinContent(ptPoint))
                                graphVarDict[s+e+p+str(tcut)].SetPointError(loosePoint,0,h_template.GetBinError(ptPoint))
                                # print "s,e,p", s,e,p,"_____ l-t=",lcut, tcut, loosePoint, tightPoint, "point=",h_template.GetBinContent(ptPoint)
                                loosePoint = loosePoint+1
                            graphVarDict[s+e+p+str(tcut)].SetMarkerColor(tightPoint+1)
                            graphVarDict[s+e+p+str(tcut)].SetLineColor(tightPoint+1)
                            if(tightPoint==0) :
                                graphVarDict[s+e+p+str(tcut)].Draw()
                            else :
                                graphVarDict[s+e+p+str(tcut)].Draw("SAME")  #"P"
                            graphVarDict[s+e+p+str(tcut)].GetXaxis().SetTitle("M_{T} cut [GeV] ")
                            graphVarDict[s+e+p+str(tcut)].GetYaxis().SetTitle("N QCD tight")
                            graphVarDict[s+e+p+str(tcut)].GetYaxis().SetTitleOffset(1)
                            graphVarDict[s+e+p+str(tcut)].SetTitle("QCD tight varying the cuts, {min}<#eta<{max},{pmin}<Pt<{pmax}, W {sign}".format(min=e, max=self.etaBinning[self.etaBinning.index(float(e))+1],pmin=p, pmax=self.ptBinning[self.ptBinning.index(float(p))+1], sign=s))
                            if(self.tightCutList.index(tcut)==len(self.tightCutList)-1) :
                                valuemaxX = ROOT.Double(0)
                                valuemaxY = ROOT.Double(0)
                                graphVarDict[s+e+p+str(tcut)].GetPoint(0,valuemaxX,valuemaxY)
                                graphVarDict[s+e+p+str(self.tightCutList[0])].GetYaxis().SetRangeUser(0,valuemaxY+valuemaxY/20)
                            legDict[e+s+p+"Variation"].AddEntry(graphVarDict[s+e+p+str(tcut)], "Iso<"+str(tcut))
                            tightPoint = tightPoint+1
                        Wvalue.Draw("SAME")
                        legDict[e+s+p+"Variation"].AddEntry(Wvalue, "W MC={val:.0g}".format(val=h_templW.GetBinContent(ptPoint)))
                        legDict[e+s+p+"Variation"].Draw("SAME")
                        ptPoint = ptPoint+1
                        canvasList.append(c_variation)
                        # c_comparison.SaveAs(self.outdir+"/plot/"+c_variation.GetName()+'.png')
            #---------------------------------------------QCD TRENDS with VARIAION OF CUTS PLOTS ---------------------------------------------#
            # print "QCD trends plots"

            for s in self.signList :
                for e in self.etaBinningS :
                    for p in self.ptBinningS :

                        h2qcd = ROOT.TH2F("h2qcd_{sign}_{eta}_{pt}".format(sign=s,eta=e, pt=p),"h2qcd_{sign}_{eta}_{pt}".format(sign=s,eta=e, pt=p),len(self.looseCutList)-1, array('f',self.looseCutList), len(self.tightCutList)-1, array('f',self.tightCutList) )
                        # tightPoint =0
                        for tcut in self.tightCutList :
                            # loosePoint = 0
                            if self.tightCutList.index(tcut)== len(self.tightCutList)-1 : continue #skip the lat
                            y_t_down = graphVarDict[s+e+p+str(tcut)].GetY()
                            tcutUP = self.tightCutList[self.tightCutList.index(tcut)+1]
                            y_t_up = graphVarDict[s+e+p+str(tcutUP)].GetY()
                            diff_list= []
                            for ll in range(len(self.looseCutList)) :
                                diff_list.append(y_t_up[ll]-y_t_down[ll])
                            for ll in range(len(self.looseCutList)) :
                                if ll == len(self.looseCutList)-1 : continue #skip the last
                                valBin = diff_list[ll]-diff_list[ll+1]
                                h2qcd.SetBinContent(ll+1,self.tightCutList.index(tcut)+1,valBin)
                                # print ll+1, self.tightCutList.index(tcut)+1, valBin


                        h2qcd.GetXaxis().SetTitle("M_{T} [GeV]")
                        h2qcd.GetYaxis().SetTitle("Relative Isolation")
                        h2qcd.GetYaxis().SetTitleOffset(1)
                        h2qcd.SetTitle("QCD trends, {min}<#eta<{max},{pmin}<Pt<{pmax}, W {sign}".format(min=e, max=self.etaBinning[self.etaBinning.index(float(e))+1],pmin=p, pmax=self.ptBinning[self.ptBinning.index(float(p))+1], sign=s))
                        canvasList.append(h2qcd)
                        c_QCDtrend = ROOT.TCanvas("c_QCDtrend_{sign}_{eta}_{pt}".format(sign=s,eta=e,pt=p),"c_QCDtrend_{sign}_{eta}_{pt}".format(sign=s,eta=e, pt=p),800,600)
                        c_QCDtrend.cd()
                        h2qcd.Draw("lego2z")
                        canvasList.append(c_QCDtrend)
                        # c_comparison.SaveAs(self.outdir+"/plot/"+c_QCDtrend.GetName()+'.png')


        outputFake = ROOT.TFile(self.outdir+"/bkg_plots"+self.nameSuff+".root","recreate")
        for h in range(len(canvasList)) :
            if "h2qcd" in canvasList[h].GetName() :
                 continue
            else :
                canvasList[h].Write()
                # canvasList[h].SaveAs(self.outdir+"/bkg_plot/"+canvasList[h].GetName()+'.png')


    def finalPlots(self, systDict =bkg_systematics, sum2Bands=False) :
        self.systDict = systDict
        self.sum2Bands =sum2Bands

        histoNameDict = {
        'comparison' : {
            'Fakes' : ['fake', 'prompt', 'validation']
            },
        'ABCDcheck' : {
            'Fakes' : ['promptSideband','prompt', 'validation','validationSigReg']
            },
        'template' : {
            'templ' : ['fake', 'prompt', 'validation']
            }
        }

        #getting canvas and histoss
        finalCanvasDict = {}
        finalPlotFileDict = {}
        finalHistoDict = {}
        finalLegDict = {}
        # finalPlotFileDict['nom']=ROOT.TFile.Open(self.outdir+"/bkg_"+"nom"+"/bkg_plots"+self.nameSuff+".root")
        # for s in self.signList :
        #     for e in self.etaBinningS :
        #         for canvas in histoNameDict :
        #             finalCanvasDict['nom'+canvas+s+e] = finalPlotFileDict['nom'].Get('c_'+canvas+'_'+s+'_'+e).Clone()
        #             # finalCanvasDict['nom'+canvas+s+e].SetDirectory(0)
        #             for name in histoNameDict[canvas] :
        #                 for histo in histoNameDict[canvas][name] :
        #                     finalHistoDict['nom'+canvas+histo+s+e] = finalCanvasDict['nom'+canvas+s+e].GetPrimitive('h'+name+'_pt_'+histo+'_'+s+'_'+e)
        #                     # finalHistoDict['nom'+canvas+histo+s+e+'error'] =finalHistoDict['nom'+canvas+histo+s+e].Clone(finalHistoDict['nom'+canvas+histo+s+e].GetName()+'_error')
        #                     print "DEUBUG", type(finalHistoDict['nom'+canvas+histo+s+e])
        #                     finalHistoDict['nom'+canvas+histo+s+e+'error'] = ROOT.TGraphAsymmErrors()#finalCanvasDict['nom'+canvas+s+e]
        #
        #                     # finalHistoDict['nom'+canvas+histo+s+e+'error'] =finalCanvasDict['nom'+canvas+s+e].GetPrimitive('h'+name+'_pt_'+histo+'_'+s+'_'+e)
        #                     finalHistoDict['nom'+canvas+histo+s+e+'error'].SetName(finalHistoDict['nom'+canvas+histo+s+e].GetName()+'_error')


                # finalCanvasDict['nom'+'comparison'+s+e] = finalPlotFileDict['nom'].Get('c_comparison_'+s+'_'+e+'.root')
                # finalCanvasDict['nom'+'ABCDcheck'+s+e] =finalPlotFileDict['nom'].Get('c_ABCDcheck_'+s+'_'+e+'.root')
                # finalCanvasDict['nom'+'template'+s+e] =finalPlotFileDict['nom'].Get('c_template_'+s+'_'+e+'.root')
                # finalHistoDict['nom'+'comparison'+'fake'+s+e] = finalCanvasDict['nom'+'comparison'+s+e].GetPrimitive('hFakes_pt_'+'fake'+s+e)
                # finalHistoDict['nom'+'comparison'+'prompt'+s+e] = finalCanvasDict['nom'+'comparison'+s+e].GetPrimitive('hFakes_pt_'+'prompt'+s+e)
                # finalHistoDict['nom'+'comparison'+'validation'+s+e] = finalCanvasDict['nom'+'comparison'+s+e].GetPrimitive('hFakes_pt_'+'validation'+s+e)
                # finalHistoDict['nom'+'template'+'fake'+s+e] = finalCanvasDict['nom'+'template'+s+e].GetPrimitive('htempl_pt_'+'fake'+s+e)
                # finalHistoDict['nom'+'template'+'prompt'+s+e] = finalCanvasDict['nom'+'template'+s+e].GetPrimitive('htempl_pt_'+'prompt'+s+e)
                # finalHistoDict['nom'+'template'+'validation'+s+e] = finalCanvasDict['nom'+'template'+s+e].GetPrimitive('htempl_pt_'+'validation'+s+e)
        for sKind, sList in self.systDict.iteritems():
            for sName in sList :
                finalPlotFileDict[sName]=ROOT.TFile.Open(self.outdir+"/bkg_"+sName+"/bkg_plots"+self.nameSuff+".root")
                for s in self.signList :
                    for e in self.etaBinningS :
                        for canvas in histoNameDict :
                            # finalCanvasDict[sName+canvas+s+e] = finalPlotFileDict[sName].Get('c_'+canvas+'_'+s+'_'+e)
                            for name in histoNameDict[canvas] :
                                for histo in histoNameDict[canvas][name] :
                                    # finalHistoDict[sName+canvas+histo+s+e] = finalCanvasDict[sName+canvas+s+e].GetPrimitive('h'+name+'_pt_'+histo+'_'+s+'_'+e)
                                    finalHistoDict[sName+canvas+histo+s+e] =   finalPlotFileDict[sName].Get('c_'+canvas+'_'+s+'_'+e).GetPrimitive('h'+name+'_pt_'+histo+'_'+s+'_'+e)
                            # finalCanvasDict[sName+'comparison'+s+e] = finalPlotFileDict[sKind+sName].Get('c_comparison_'+s+'_'+e+'.root')
                            # finalCanvasDict[sName+'ABCDcheck'+s+e] =finalPlotFileDict[sKind+sName].Get('c_ABCDcheck_'+s+'_'+e+'.root')
                            # finalCanvasDict[sName+'template'+s+e] =finalPlotFileDict[sKind+sName].Get('c_template_'+s+'_'+e+'.root')

        #building of error bands on the "nom" hisotgrams
        # outputFinal = ROOT.TFile(self.outdir+"/final_plots"+self.nameSuff+".root","recreate")
        # outputFinal.cd()
        finalPlotFileDict['nom']=ROOT.TFile.Open(self.outdir+"/bkg_"+"nom"+"/bkg_plots"+self.nameSuff+".root")
        for s in self.signList :
            for e in self.etaBinningS :
                for canvas in histoNameDict :
                    finalCanvasDict['nom'+canvas+s+e] = finalPlotFileDict['nom'].Get('c_'+canvas+'_'+s+'_'+e).Clone()
                    finalCanvasDict['nom'+canvas+s+e].cd()
                    fillSyle =0
                    for name in histoNameDict[canvas] :
                        for histo in histoNameDict[canvas][name] :
                            finalHistoDict['nom'+canvas+histo+s+e] = finalCanvasDict['nom'+canvas+s+e].GetPrimitive('h'+name+'_pt_'+histo+'_'+s+'_'+e)
                            finalHistoDict['nom'+canvas+histo+s+e+'error'] = ROOT.TGraphAsymmErrors()#finalCanvasDict['nom'+canvas+s+e]
                            finalHistoDict['nom'+canvas+histo+s+e+'error'].SetName(finalHistoDict['nom'+canvas+histo+s+e].GetName()+'_error')

                            for p in self.ptBinning :
                                if(p==self.ptBinning[-1]) : continue
                                varErr = []
                                varErr.append(finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))
                                varErrSum2Up = 0
                                varErrSum2Down = 0.
                                for sKind, sList in self.systDict.iteritems():
                                    for sName in sList :
                                        varErr.append(finalHistoDict[sName+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))
                                        if "Up" in sName :
                                            varErrSum2Up = varErrSum2Up+ (finalHistoDict[sName+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)-finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))**2
                                        if "Down" in sName :
                                            varErrSum2Down = varErrSum2Down+ (finalHistoDict[sName+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)-finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))**2
                                varErrSum2Up = math.sqrt(varErrSum2Up)
                                varErrSum2Down = math.sqrt(varErrSum2Down)
                                varErr.sort()
                                minBand = varErr[0]
                                maxBand = varErr[len(varErr)-1]
                                errLow = finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)-minBand
                                errHigh = maxBand-finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)
                                if self.sum2Bands :
                                    errLow = varErrSum2Down
                                    errHigh = varErrSum2Up
                                # if histo == 'fake' :
                                    # print "DEBUG", name, "s,e,p", s,e,p,"___   err=",errLow, errHigh, minBand, maxBand, finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1)
                                # print "DEBUG:",histo, self.ptBinning.index(float(p)),finalHistoDict['nom'+canvas+histo+s+e].GetBinCenter(self.ptBinning.index(float(p))+1), errLow, errHigh
                                # symBand = (maxBand-minBand)/2
                                # print "WARNING: SYMMETRIC BANDS FOR SYST USED!!!"
                                # finalHistoDict['nom'+canvas+histo+s+e+'error'].SetBinError(self.ptBinning.index(float(p))+1,symBand)
                                finalHistoDict['nom'+canvas+histo+s+e+'error'].SetPoint(self.ptBinning.index(float(p)),finalHistoDict['nom'+canvas+histo+s+e].GetBinCenter(self.ptBinning.index(float(p))+1),finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))
                                finalHistoDict['nom'+canvas+histo+s+e+'error'].SetPointEYhigh(self.ptBinning.index(float(p)),errHigh)
                                finalHistoDict['nom'+canvas+histo+s+e+'error'].SetPointEYlow(self.ptBinning.index(float(p)),errLow)
                                finalHistoDict['nom'+canvas+histo+s+e+'error'].SetPointEXhigh(self.ptBinning.index(float(p)),finalHistoDict['nom'+canvas+histo+s+e].GetBinWidth(self.ptBinning.index(float(p))+1)/2)
                                finalHistoDict['nom'+canvas+histo+s+e+'error'].SetPointEXlow(self.ptBinning.index(float(p)),finalHistoDict['nom'+canvas+histo+s+e].GetBinWidth(self.ptBinning.index(float(p))+1)/2)
                                # print "DEBUG, p=", finalHistoDict['nom'+canvas+histo+s+e].GetBinCenter(self.ptBinning.index(float(p))+1), "|||||", p, "|||||", self.ptBinning.index(float(p))+1

                            finalHistoDict['nom'+canvas+histo+s+e+'error'].SetLineColor(finalHistoDict['nom'+canvas+histo+s+e].GetLineColor()-3)
                            finalHistoDict['nom'+canvas+histo+s+e+'error'].SetMarkerColor(finalHistoDict['nom'+canvas+histo+s+e].GetLineColor()-3)
                            finalHistoDict['nom'+canvas+histo+s+e+'error'].SetFillColor(finalHistoDict['nom'+canvas+histo+s+e].GetLineColor()-3)
                            finalHistoDict['nom'+canvas+histo+s+e+'error'].SetFillStyle(3003+fillSyle)
                            finalHistoDict['nom'+canvas+histo+s+e+'error'].SetLineWidth(1)

                            finalHistoDict['nom'+canvas+histo+s+e+'error'].Draw("SAME 0P5")
                            finalCanvasDict['nom'+canvas+s+e].Update()
                            finalCanvasDict['nom'+canvas+s+e].Modified()
                            fillSyle =fillSyle+1

                            # c1 = finalCanvasDict['nom'+canvas+s+e].Clone()
                            # c1.cd()
                            # c1.Update()
                            # c1.Modified()
                            # c1.Write()
                            # finalCanvasDict['nom'+canvas+s+e].Write()
                            # outputFinal.cd()
                            # print "DEBUG >>>>", canvas, s, e, histo
                            # print "name", finalCanvasDict['nom'+canvas+s+e].GetName()
                            # print "primitive:", list(finalCanvasDict['nom'+canvas+s+e].GetListOfPrimitives())
                            # print "error name", finalHistoDict['nom'+canvas+histo+s+e+'error'].GetName()

        #ratio plot creation
        for s in self.signList :
            for e in self.etaBinningS :
                for canvas in histoNameDict :
                    for name in histoNameDict[canvas] :
                        for histo in histoNameDict[canvas][name] :
                            c_ratioSyst = ROOT.TCanvas("c_ratioSyst_{sign}_{eta}_{canvas}_{histo}".format(sign=s,eta=e,canvas=canvas,histo=histo),"c_ratioSyst_{sign}_{eta}_{canvas}_{histo}".format(sign=s,eta=e,canvas=canvas,histo=histo),800,600)
                            c_ratioSyst.cd()
                            c_ratioSyst.SetGridx()
                            c_ratioSyst.SetGridy()
                            finalLegDict[e+s+canvas+histo+"ratioSyst"] = ROOT.TLegend(0.1,0.7,0.48,0.9)

                            sameFlag = True
                            colorNumber = 1
                            colorList = [600,616,416,632,432,800,900]
                            colorCounter = 0

                            for sKind, sList in self.systDict.iteritems():
                                colorNumber = colorList[colorCounter]
                                colorCounter = colorCounter+1
                                for sName in sList :
                                    colorNumber = colorNumber-2
                                    if colorNumber < colorList[colorCounter-1]-10 :
                                        colorNumber = colorList[colorCounter]+2
                                    # hRatio = finalHistoDict['nom'+canvas+histo+s+e].Clone(finalHistoDict[sName+canvas+histo+s+e].GetName()+'_'+sName+'_ratio')
                                    # print finalHistoDict[sName+canvas+histo+s+e].GetName()+'_ratio'
                                    # print "DEBUG memory leak", sName, canvas, histo, s,e
                                    # print "name ", canvas+'_'+finalHistoDict[sName+canvas+histo+s+e].GetName()+'_'+sName+'_ratio'
                                    finalHistoDict[sName+canvas+histo+s+e+'ratio'] = ROOT.TH1F(canvas+'_'+finalHistoDict[sName+canvas+histo+s+e].GetName()+'_'+sName+'_ratio',canvas+'_'+finalHistoDict[sName+canvas+histo+s+e].GetName()+'_'+sName+'_ratio',len(self.ptBinning)-1, array('f',self.ptBinning))

                                    #finalHistoDict[sName+canvas+histo+s+e+'ratio'].SetName(finalHistoDict[sName+canvas+histo+s+e].GetName()+'_'+sName+'_ratio')
                                    # finalHistoDict[sName+canvas+histo+s+e+'ratio'].SetTitle(finalHistoDict[sName+canvas+histo+s+e].GetName()+'_'+sName+'_ratio')
                                    # print sName, finalHistoDict[sName+canvas+histo+s+e].GetNbinsX(), "nom", finalHistoDict['nom'+canvas+histo+s+e].GetNbinsX(), finalHistoDict[sName+canvas+histo+s+e+'ratio'].GetNbinsX()
                                    finalHistoDict[sName+canvas+histo+s+e+'ratio'].Divide(finalHistoDict[sName+canvas+histo+s+e],finalHistoDict['nom'+canvas+histo+s+e],1,1)
                                    finalHistoDict[sName+canvas+histo+s+e+'ratio'].SetLineColor(colorNumber)
                                    c_ratioSyst.cd()
                                    if sameFlag :
                                        finalHistoDict[sName+canvas+histo+s+e+'ratio'].Draw()
                                        finalHistoDict[sName+canvas+histo+s+e+'ratio'].GetYaxis().SetRangeUser(0.8,1.2)
                                        # print "quante volte entro qui?"
                                    else :
                                        finalHistoDict[sName+canvas+histo+s+e+'ratio'].Draw("SAME")
                                        finalHistoDict[sName+canvas+histo+s+e+'ratio'].GetYaxis().SetRangeUser(0,3)
                                        # print "same!", sName
                                    finalHistoDict[sName+canvas+histo+s+e+'ratio'].SetLineWidth(3)
                                    sameFlag=False
                                    finalLegDict[e+s+canvas+histo+"ratioSyst"].AddEntry(finalHistoDict[sName+canvas+histo+s+e+'ratio'], sName)

                            finalLegDict[e+s+canvas+histo+"ratioSyst"].Draw("SAME")

                            finalCanvasDict['ratio'+canvas+histo+s+e] = c_ratioSyst
                            # print "DEBUG", list(finalCanvasDict['ratio'+canvas+histo+s+e].GetListOfPrimitives())



        #unrolled plot creation

        #binning
        unrolledPtEta= list(self.ptBinning)
        # lastPtValue = self.ptBinning[-1] 
        # lastPtValue = 100 #here set to 100 the unroll limit to read clearly the pt 
        intervalPtBin = []
        for p in self.ptBinning[:-1] :
            intervalPtBin.append(self.ptBinning[self.ptBinning.index(p)+1]-self.ptBinning[self.ptBinning.index(p)])
        print "intervals", intervalPtBin
        
        for e in range(len(self.etaBinning)-2) :
             # tempPt = []
            # tempPt = [p+lastPtValue*(self.etaBinning.index(e)+1) for p in self.ptBinning]
            # tempPt = [p-self.ptBinning[0]+self.ptBinning[-1]*(self.etaBinning.index(e)+1) for p in self.ptBinning[1:]]
            for p in intervalPtBin :
                unrolledPtEta.append(unrolledPtEta[-1]+p)
                
            # unrolledPtEta.extend(tempPt)
        print "unrolledPtEta", unrolledPtEta
        # print "eta binning", self.etaBinning
        # print "pt binning", self.ptBinning

        #final plot and syst
        for s in self.signList :
            for canvas in histoNameDict :
                c_unrolled = ROOT.TCanvas("c_unrolled_{canvas}_{sign}".format(canvas=canvas,sign=s),"c_unrolled_{canvas}_{sign}".format(canvas=canvas,sign=s),800,600)
                c_unrolled.cd()
                c_unrolled.SetGridx()
                c_unrolled.SetGridy()
                finalLegDict[s+canvas+'unrolled'] = ROOT.TLegend(0.1,0.7,0.48,0.9)
                sameFlagUNR=True
                for name in histoNameDict[canvas] :
                    for histo in histoNameDict[canvas][name] :
                        finalHistoDict[s+canvas+histo+'unrolled'] = ROOT.TH1F(canvas+'_'+finalHistoDict['nom'+canvas+histo+s+'0'].GetName()+'_unrolled',canvas+'_'+finalHistoDict['nom'+canvas+histo+s+'0'].GetName()+'_unrolled',len(unrolledPtEta)-1, array('f',unrolledPtEta))
                        finalHistoDict[s+canvas+histo+'unrolled'].SetLineColor(finalHistoDict['nom'+canvas+histo+s+'0'].GetLineColor())
                        finalHistoDict[s+canvas+histo+'unrolled'].SetLineWidth(finalHistoDict['nom'+canvas+histo+s+'0'].GetLineWidth())
                        finalHistoDict[s+canvas+histo+'unrolled'].GetXaxis().SetTitle(finalHistoDict['nom'+canvas+histo+s+'0'].GetXaxis().GetTitle())
                        finalHistoDict[s+canvas+histo+'unrolled'].GetYaxis().SetTitle(finalHistoDict['nom'+canvas+histo+s+'0'].GetYaxis().GetTitle())
                        

                        finalHistoDict[s+canvas+histo+'unrolled'+'error'] = ROOT.TGraphAsymmErrors()
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetName(finalHistoDict[s+canvas+histo+'unrolled'].GetName()+'_error')
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetLineColor(finalHistoDict['nom'+canvas+histo+s+'0'+'error'].GetLineColor())
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetMarkerColor(finalHistoDict['nom'+canvas+histo+s+'0'+'error'].GetLineColor())
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetFillColor(finalHistoDict['nom'+canvas+histo+s+'0'+'error'].GetLineColor())
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetFillStyle(finalHistoDict['nom'+canvas+histo+s+'0'+'error'].GetFillStyle())
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetLineWidth(finalHistoDict['nom'+canvas+histo+s+'0'+'error'].GetLineWidth())



                        for e in self.etaBinningS :
                            for p in self.ptBinningS :
                                indexUNR = self.etaBinning.index(float(e))*len(self.etaBinning)+self.ptBinning.index(float(p))
                                # print "index", indexUNR
                                finalHistoDict[s+canvas+histo+'unrolled'].SetBinContent(indexUNR+1,finalHistoDict['nom'+canvas+histo+s+e].GetBinContent(self.ptBinning.index(float(p))+1))
                                finalHistoDict[s+canvas+histo+'unrolled'].SetBinError(indexUNR+1,finalHistoDict['nom'+canvas+histo+s+e].GetBinError(self.ptBinning.index(float(p))+1))

                                finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetPoint(indexUNR,finalHistoDict[s+canvas+histo+'unrolled'].GetBinCenter(indexUNR+1),finalHistoDict[s+canvas+histo+'unrolled'].GetBinContent(indexUNR+1))
                                finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetPointEYhigh(indexUNR,finalHistoDict['nom'+canvas+histo+s+e+'error'].GetErrorYhigh(self.ptBinning.index(float(p))))
                                finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetPointEYlow(indexUNR,finalHistoDict['nom'+canvas+histo+s+e+'error'].GetErrorYlow(self.ptBinning.index(float(p))))
                                finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetPointEXhigh(indexUNR,finalHistoDict[s+canvas+histo+'unrolled'].GetBinWidth(indexUNR+1)/2)
                                finalHistoDict[s+canvas+histo+'unrolled'+'error'].SetPointEXlow(indexUNR,finalHistoDict[s+canvas+histo+'unrolled'].GetBinWidth(indexUNR+1)/2)

                        if sameFlagUNR :
                            finalHistoDict[s+canvas+histo+'unrolled'].Draw()
                            if(canvas!="template") :
                                finalHistoDict[s+canvas+histo+'unrolled'].GetYaxis().SetRangeUser(0,1.1)
                        else :
                            finalHistoDict[s+canvas+histo+'unrolled'].Draw("SAME")
                        sameFlagUNR=False
                        finalLegDict[s+canvas+'unrolled'].AddEntry(finalHistoDict[s+canvas+histo+'unrolled'], histo)
                        finalHistoDict[s+canvas+histo+'unrolled'+'error'].Draw("SAME 0P5")
                finalLegDict[s+canvas+'unrolled'].Draw("SAME")
                finalCanvasDict['unrolled'+canvas+histo+s] = c_unrolled


        #ratios unrolled
        for s in self.signList :
            for canvas in histoNameDict :
                for name in histoNameDict[canvas] :
                    for histo in histoNameDict[canvas][name] :
                        c_ratioSyst_unrolled = ROOT.TCanvas("c_ratioSyst_unrolled_{sign}_{canvas}_{histo}".format(sign=s,canvas=canvas,histo=histo),"c_ratioSyst_unrolled_{sign}_{canvas}_{histo}".format(sign=s,canvas=canvas,histo=histo),800,600)
                        c_ratioSyst_unrolled.cd()
                        c_ratioSyst_unrolled.SetGridx()
                        c_ratioSyst_unrolled.SetGridy()
                        finalLegDict[s+canvas+histo+"ratioSyst_unrolled"] = ROOT.TLegend(0.1,0.7,0.48,0.9)

                        sameFlagUNR = True
                        # colorNumber = 1
                        # colorList = [600,616,416,632,432,800,900]
                        # colorCounter = 0

                        for sKind, sList in self.systDict.iteritems():
                            # colorNumber = colorList[colorCounter]
                            # colorCounter = colorCounter+1
                            for sName in sList :
                                # colorNumber = colorNumber-2
                                # if colorNumber < colorList[colorCounter-1]-10 :
                                    # colorNumber = colorList[colorCounter]+2
                                # print "DEBUG momory leak", finalHistoDict[sName+canvas+histo+s+'0ratio'].GetName()+'_ratio_unrolled', s, canvas, sName, histo
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'] = ROOT.TH1F(finalHistoDict[sName+canvas+histo+s+'0ratio'].GetName()+'_ratio_unrolled',finalHistoDict[sName+canvas+histo+s+'0ratio'].GetName()+'_ratio_unrolled',len(unrolledPtEta)-1, array('f',unrolledPtEta))
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].SetLineColor(finalHistoDict[sName+canvas+histo+s+'0ratio'].GetLineColor())
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].SetLineWidth(finalHistoDict[sName+canvas+histo+s+'0ratio'].GetLineWidth())
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].GetXaxis().SetTitle(finalHistoDict[sName+canvas+histo+s+'0ratio'].GetXaxis().GetTitle())
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].GetYaxis().SetTitle(finalHistoDict[sName+canvas+histo+s+'0ratio'].GetYaxis().GetTitle())

                                for e in self.etaBinningS :
                                    for p in self.ptBinningS :
                                        indexUNR = self.etaBinning.index(float(e))*len(self.etaBinning)+self.ptBinning.index(float(p))
                                        finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].SetBinContent(indexUNR,finalHistoDict[sName+canvas+histo+s+e+'ratio'].GetBinContent(self.ptBinning.index(float(p))+1))
                                        finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].SetBinError(indexUNR,finalHistoDict[sName+canvas+histo+s+e+'ratio'].GetBinError(self.ptBinning.index(float(p))+1))
                                # finalHistoDict[sName+canvas+histo+s+e+'ratio'] = ROOT.TH1F(finalHistoDict[sName+canvas+histo+s+e].GetName()+'_'+sName+'_ratio',finalHistoDict[sName+canvas+histo+s+e].GetName()+'_'+sName+'_ratio',len(self.ptBinning)-1, array('f',self.ptBinning))
                                # finalHistoDict[sName+canvas+histo+s+e+'ratio'].Divide(finalHistoDict[sName+canvas+histo+s+e],finalHistoDict['nom'+canvas+histo+s+e],1,1)
                                # finalHistoDict[sName+canvas+histo+s+e+'ratio'].SetLineColor(colorNumber)
                                c_ratioSyst_unrolled.cd()
                                if sameFlagUNR :
                                    finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].Draw()
                                    finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].GetYaxis().SetRangeUser(0.8,1.2)
                                    # print "quante volte entro qui?"
                                else :
                                    finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].Draw("SAME")
                                    finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].GetYaxis().SetRangeUser(0,3)
                                    # print "same!", sName
                                finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'].SetLineWidth(3)
                                sameFlagUNR=False
                                finalLegDict[s+canvas+histo+"ratioSyst_unrolled"].AddEntry(finalHistoDict[sName+s+canvas+histo+'ratio_unrolled'], sName)

                        finalLegDict[s+canvas+histo+"ratioSyst_unrolled"].Draw("SAME")

                        finalCanvasDict['ratio_unrolled'+canvas+histo+s] = c_ratioSyst_unrolled
                        # print "DEBUG", list(finalCanvasDict['ratio'+canvas+histo+s+e].GetListOfPrimitives())



        # outputFinal.Close()
        outputFinal = ROOT.TFile(self.outdir+"/final_plots"+self.nameSuff+".root","recreate")
        # outputFinal.cd()
        # preliminary_dir = output.mkdir("RatiosVSMt")
        # preliminary_dir.cd()
        dirFinalDict = {}
        for s in self.signList :
            for e in self.etaBinningS :
                dirFinalDict[s+e] =    outputFinal.mkdir(s+"_eta"+e)
                # dirFinalDict[s+e+'bis'] =    dirFinalDict[s+e].GetDirectory()
                # dirFinalDict[s+e].cd()   #DECOMMENTA QUESTO SE NON FUNZIONANO LE DIRECTORY
            dirFinalDict[s+'unrolled'] =    outputFinal.mkdir(s+'_unrolled')
        for ind, obj in finalCanvasDict.iteritems():
                # print ind, obj.GetName()
                for s in self.signList :
                    
                    if 'unrolled' in ind and s in ind:
                        # print "inside unrolled"
                        dirFinalDict[s+'unrolled'].cd()
                        obj.Write()
                    else :
                        for e in self.etaBinningS :
                            if ind.endswith(s+e) :
                            # print "ind, e", ind, e
                            # if "ratio" in ind : continue
                                dirFinalDict[s+e].cd()
                                obj.Write()
        # for ind, obj in finalHistoDict.iteritems():
        #         if 'error' in ind : obj.Write()
    
    def strategy_syst(self, preOutDir) :
        self.preOutDir = preOutDir
        
        systlist = ['mT_fit', 'EWSFdown10','EWSFup10', 'slopeDown10','slopeUp10', 'mT_fit_2deg', 'mT_countHigh', 'MET_fit', 'MET_countHigh']
        colorList = [600,616,416,632,432,800,900,881,402]
        straSystDict_file = {}
        straSystDict_histo = {}
        straSystDict_leg = {}
        straSystDict_canvas = {}
        straSystDict_pad = {}
        
        unrolledPtEta= list(self.ptBinning)
        intervalPtBin = []
        for p in self.ptBinning[:-1] :
            intervalPtBin.append(self.ptBinning[self.ptBinning.index(p)+1]-self.ptBinning[self.ptBinning.index(p)])
        for e in range(len(self.etaBinning)-2) :
            for p in intervalPtBin :
                unrolledPtEta.append(unrolledPtEta[-1]+p)        
        
        for var in systlist :
            straSystDict_file[var]=ROOT.TFile.Open(self.preOutDir+"/bkg_"+var+"/final_plots"+self.nameSuff+".root")
            for s in self.signList : 
                 
                straSystDict_histo[var+s] =   straSystDict_file[var].Get(s+'_unrolled/c_unrolled_template'+'_'+s).GetPrimitive('template_htempl_pt_fake_'+s+'_0_unrolled')
                straSystDict_histo[var+s].SetLineColor(colorList[systlist.index(var)])
                if var == 'mT_fit' :
                    straSystDict_histo[var+s+'error'] =   straSystDict_file[var].Get(s+'_unrolled/c_unrolled_template'+'_'+s).GetPrimitive('template_htempl_pt_fake_'+s+'_0_unrolled_error') 
                    straSystDict_histo[var+s+'error'].SetLineColor(colorList[systlist.index(var)])

                if var != 'mT_fit' :
                    straSystDict_histo[var+s+'ratio'] = ROOT.TH1F(straSystDict_histo[var+s].GetName()+'_ratio',straSystDict_histo[var+s].GetName()+'_ratio',len(unrolledPtEta)-1, array('f',unrolledPtEta))
                    straSystDict_histo[var+s+'ratio'].Divide(straSystDict_histo[var+s],straSystDict_histo['mT_fit'+s],1,1)
                    straSystDict_histo[var+s+'ratio'].SetLineColor(colorList[systlist.index(var)])
                    straSystDict_histo[var+s+'ratio'].SetLineWidth(2)
                    straSystDict_histo[var+s+'ratio'].GetYaxis().SetTitle("Syst/Nom")
    
        
        for s in self.signList :
            c_straSyst = ROOT.TCanvas("c_straSyst_{sign}".format(sign=s),"c_straSyst_{sign}".format(sign=s),800,700)
            c_straSyst.cd()
            c_straSyst.SetGridx()
            c_straSyst.SetGridy()
            
            straSystDict_leg[s] = ROOT.TLegend(0.30,0.10,0.80,0.4)
            
            p_straSyst_histo = ROOT.TPad("p_straSyst_histo_"+s, "c_straSyst_"+s,0,0.5,1,1)
            p_straSyst_ratio = ROOT.TPad("p_straSyst_ratio_"+s, "c_straSyst_"+s,0,0,1,0.5)
            p_straSyst_histo.SetBottomMargin(0.02)
            p_straSyst_histo.Draw()
            p_straSyst_ratio.SetTopMargin(0)
            p_straSyst_ratio.SetBottomMargin(0.25)
            p_straSyst_ratio.Draw()
            
            sameFlagSS=True
            sameFlagSS_ratio = True
            for var in systlist :
                p_straSyst_histo.cd()
                if sameFlagSS :
                    straSystDict_histo[var+s].Draw()
                else : 
                    straSystDict_histo[var+s].Draw("SAME")
                straSystDict_leg[s].AddEntry(straSystDict_histo[var+s],var)
                sameFlagSS=False
                
                if var!= "mT_fit" :
                    p_straSyst_ratio.cd()
                    if sameFlagSS_ratio :
                        straSystDict_histo[var+s+'ratio'].Draw()
                    else :
                        straSystDict_histo[var+s+'ratio'].Draw("SAME")
                    sameFlagSS_ratio = False                    
            p_straSyst_histo.cd()
            straSystDict_leg[s].Draw("SAME")
                
            straSystDict_canvas[s] = c_straSyst
        
        outputraSyst= ROOT.TFile(self.preOutDir+"/bkg_mT_fit/strategy_syst"+self.nameSuff+".root","recreate")
        for s in self.signList :
            straSystDict_canvas[s].Write()
        
        

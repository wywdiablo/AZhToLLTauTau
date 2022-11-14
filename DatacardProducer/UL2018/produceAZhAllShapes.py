#! /usr/bin/env python
# Author: Yiwen Wen (Oct 2021)
# Description: Create datacards for combine
import sys, os,array
from CombineHarvester.CombineTools import ch
import ROOT as R
import math

Zchannels = ['ee','mm']
hchannels = ['em','et','mt','tt']
AZhCats = ['0btag','btag']
Syst = ['Up','Down']
Uncerts =['e','m','t','eft','mft','unclmet']
#TauIDSyst = ['VSjetSF']
TauIDSyst = ['VSjetSF','VSeSF','VSmuSF','Prefire']
FakeSyst = ['JetToEle','JetToMu','JetToTau']
Era = 'UL2018'
for iCat in AZhCats :
    for iZch in Zchannels :
        for ihch in hchannels :
            cut = ""
            if iCat is '0btag' : cut = "*(nbtag==0 && m_sv>90 && m_sv<180)"
            else : cut = "*(nbtag>=1 && m_sv>90 && m_sv<180)"
            if iCat is '0btag' and ihch is 'tt' : cut = "*(nbtag==0 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
            if iCat is 'btag' and ihch is 'tt' : cut = "*(nbtag>=1 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
            produceCentralShape_cmd = 'root -l -q \'ProduceAZhDataCards.C("%s","Mass_svc","m^{c}_{ll#tau#tau} [GeV]","Events",true,"%s","%s","%s","%s")\'' %(Era,iZch,ihch,cut,iCat)
            #os.system(produceCentralShape_cmd)
            #Produce Systematic Shapes
            for iuncer in Uncerts :
                for updown in Syst :
                    produceSystShape_cmd = 'root -l -q \'ProduceAZhDataCards_Syst.C("%s","Mass_svc_%sscale%s","CMS_%sshape_2018%s","m^{c}_{ll#tau#tau} [GeV]","Events",true,"%s","%s","%s","%s")\'' %(Era,iuncer,updown,iuncer,updown,iZch,ihch,cut,iCat)
                    #os.system(produceSystShape_cmd)
            #Produce Weight Systematic Shapes
            for iTauIDSyst in TauIDSyst :
                for updown in Syst :
                    if iTauIDSyst is 'VSjetSF': systname = 'CMS_eff_t_2018'
                    elif iTauIDSyst is 'VSeSF': systname = 'CMS_fake_e_2018'
                    elif iTauIDSyst is 'VSmuSF': systname = 'CMS_fake_m_2018'
                    elif iTauIDSyst is 'Prefire': systname = 'CMS_prefire'
                    produceTauIDSystShape_cmd = 'root -l -q \'ProduceAZhDataCards_Syst_TauIDWeight.C("%s","Mass_svc","%s","%s","%s","m^{c}_{ll#tau#tau} [GeV]","Events",true,"%s","%s","%s","%s")\''%(Era,systname,updown,iTauIDSyst,iZch,ihch,cut,iCat)
                    #os.system(produceTauIDSystShape_cmd)
            #Produce Fake Systematic Shapes
            for iFakeSyst in FakeSyst :
                for updown in Syst :
                    if iFakeSyst is 'JetToEle' : systname = 'CMS_AZh_fr_JetFakeE_2018'
                    elif iFakeSyst is 'JetToMu' : systname = 'CMS_AZh_fr_JetFakeM_2018'
                    elif iFakeSyst is 'JetToTau' : systname = 'CMS_AZh_fr_JetFakeT_2018'
                    produceFakeSystShape_cmd = 'root -l -q \'ProduceAZhDataCards_Syst_FakeRates.C("%s","Mass_svc","%s","%s","%s","m^{c}_{ll#tau#tau} [GeV]","Events",true,"%s","%s","%s","%s")\''%(Era,systname,updown,iFakeSyst,iZch,ihch,cut,iCat)
                    os.system(produceFakeSystShape_cmd)
                    
            #Produce B-tag Systematics Shapes(cut induced)
            if iCat is '0btag' :
                BtagUpcut = "*(nbtag_btagUp==0 && m_sv>90 && m_sv<180)"
                BtagDowncut = "*(nbtag_btagDown==0 && m_sv>90 && m_sv<180)"
                MisTagUpcut = "*(nbtag_mistagUp==0 && m_sv>90 && m_sv<180)"
                MisTagDowncut = "*(nbtag_mistagDown==0 && m_sv>90 && m_sv<180)"
                JECUpcut = "*(nbtagUp==0 && m_sv>90 && m_sv<180)"
                JECDowncut = "*(nbtagDown==0 && m_sv>90 && m_sv<180)"
            else :
                BtagUpcut = "*(nbtag_btagUp>=1 && m_sv>90 && m_sv<180)"
                BtagDowncut = "*(nbtag_btagDown>=1 && m_sv>90 && m_sv<180)"
                MisTagUpcut = "*(nbtag_mistagUp>=1 && m_sv>90 && m_sv<180)"
                MisTagDowncut = "*(nbtag_mistagDown>=1 && m_sv>90 && m_sv<180)"
                JECUpcut = "*(nbtagUp>=1 && m_sv>90 && m_sv<180)"
                JECDowncut = "*(nbtagDown>=1 && m_sv>90 && m_sv<180)"
                
            if iCat is '0btag' and ihch is 'tt' :
                BtagUpcut = "*(nbtag_btagUp==0 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
                BtagDowncut = "*(nbtag_btagDown==0 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
                MisTagUpcut = "*(nbtag_mistagUp==0 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
                MisTagDowncut = "*(nbtag_mistagDown==0 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
                JECUpcut = "*(nbtagUp==0 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
                JECDowncut = "*(nbtagDown==0 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
                
            if iCat is 'btag' and ihch is 'tt' :
                BtagUpcut = "*(nbtag_btagUp>=1 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
                BtagDowncut = "*(nbtag_btagDown>=1 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
                MisTagUpcut = "*(nbtag_mistagUp>=1 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
                MisTagDowncut = "*(nbtag_mistagDown>=1 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
                JECUpcut = "*(nbtagUp>=1 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
                JECDowncut = "*(nbtagDown>=1 && (pt_3+pt_4)>60 && m_sv>90 && m_sv<180)"
 
            produceBtagSystShapeBtagUp_cmd = 'root -l -q \'ProduceAZhDataCards_Syst_Btag.C("%s","Mass_svc","CMS_DeepFlavBtag_2018","Up","Btag","m^{c}_{ll#tau#tau} [GeV]","Events",true,"%s","%s","%s","%s")\''%(Era,iZch,ihch,BtagUpcut,iCat)
            produceBtagSystShapeBtagDown_cmd = 'root -l -q \'ProduceAZhDataCards_Syst_Btag.C("%s","Mass_svc","CMS_DeepFlavBtag_2018","Down","Btag","m^{c}_{ll#tau#tau} [GeV]","Events",true,"%s","%s","%s","%s")\''%(Era,iZch,ihch,BtagDowncut,iCat)
            produceBtagSystShapeMistagUp_cmd = 'root -l -q \'ProduceAZhDataCards_Syst_Btag.C("%s","Mass_svc","CMS_DeepFlavMistag_2018","Up","MisTag","m^{c}_{ll#tau#tau} [GeV]","Events",true,"%s","%s","%s","%s")\''%(Era,iZch,ihch,MisTagUpcut,iCat)
            produceBtagSystShapeMistagDown_cmd = 'root -l -q \'ProduceAZhDataCards_Syst_Btag.C("%s","Mass_svc","CMS_DeepFlavMistag_2018","Down","MisTag","m^{c}_{ll#tau#tau} [GeV]","Events",true,"%s","%s","%s","%s")\''%(Era,iZch,ihch,MisTagDowncut,iCat)
            produceJECSystShapeUp_cmd = 'root -l -q \'ProduceAZhDataCards_Syst_Btag.C("%s","Mass_svc","CMS_JEC_2018","Up","JEC","m^{c}_{ll#tau#tau} [GeV]","Events",true,"%s","%s","%s","%s")\''%(Era,iZch,ihch,JECUpcut,iCat)
            produceJECSystShapeDown_cmd = 'root -l -q \'ProduceAZhDataCards_Syst_Btag.C("%s","Mass_svc","CMS_JEC_2018","Down","JEC","m^{c}_{ll#tau#tau} [GeV]","Events",true,"%s","%s","%s","%s")\''%(Era,iZch,ihch,JECDowncut,iCat)
            
            #os.system(produceBtagSystShapeBtagUp_cmd)
            #os.system(produceBtagSystShapeBtagDown_cmd)
            #os.system(produceBtagSystShapeMistagUp_cmd)
            #os.system(produceBtagSystShapeMistagDown_cmd)
            #os.system(produceJECSystShapeUp_cmd)
            #os.system(produceJECSystShapeDown_cmd)

            rm_cmd = 'rm ./Combineworkspace/shapes/AZh_%s_%s%s_Mass_svc_%s.root'%(Era,iZch,ihch,iCat)
            hadd_cmd = 'hadd ./Combineworkspace/shapes/AZh_%s_%s%s_Mass_svc_%s.root ./Combineworkspace/shapes/produceArea/AZh_%s_%s%s_Mass_svc_*_%s.root' %(Era,iZch,ihch,iCat,Era,iZch,ihch,iCat)
            os.system(rm_cmd)
            os.system(hadd_cmd)


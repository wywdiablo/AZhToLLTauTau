#! /usr/bin/env python
# Author: Yiwen Wen (Feb 2021)
# Description: Create datacards for combine
import sys,os,array
from CombineHarvester.CombineTools import ch
import ROOT as R
import math

R.PyConfig.IgnoreCommandLineOptions = True
R.gROOT.SetBatch(R.kTRUE)
R.gStyle.SetOptStat(0)
canvas = R.TCanvas()
canvas2 = R.TCanvas()

wps = 'Medium_VLoose_Tight'
pts = ['20to30','30to40','40to50','50to60','60to70','70to80','80to90','Gt90']
dms = ['0','1','10','11']
bins = ['Num','Den']
nptbins = 8
h2_FR_prefit = R.TH2D("JetToTauPrefitFR", "Jet To #tau_{h} Prefit Fake Rate",nptbins,20,100,4,0,4)
h2_FR_postfit = R.TH2D("JetToTauPostfitFR", "Jet To #tau_{h} Postfit Fake Rate",nptbins,20,100,4,0,4)

outfile = R.TFile("JetTauFakeRate_%s.root"%(wps), "RECREATE")
ndm = 1

for idm in dms :
    print '<<<<<<< decay mode: ', idm
    h_FR_prefit = R.TH1D("PrefitFR_DM%s"%idm,"Prefit Fake Rate DM%s"%idm, nptbins, 20, 100)
    h_FR_postfit = R.TH1D("PostfitFR_DM%s"%idm,"Postfit Fake Rate DM%s"%idm, nptbins, 20, 100)
    nptbin = 1
    for ipt in pts :
        print '<<<<<<<<<<<<<< pt range : ', ipt
        
        Num_prefit_value = R.Double()
        Den_prefit_value = R.Double()
        Num_postfit_value = R.Double()
        Den_postfit_value = R.Double()
        
        Num_data_prefit_err = R.Double()
        Den_data_prefit_err = R.Double()
        Num_prompt_prefit_err = R.Double()
        Den_prompt_prefit_err = R.Double()
        Num_data_postfit_err = R.Double()
        Den_data_postfit_err = R.Double()
        Num_prompt_postfit_err = R.Double()
        Den_prompt_postfit_err = R.Double()
        for bin in bins :
            print '<<<<<<<<<<<<<< bin  : ', bin

            cb = ch.CombineHarvester()
            #cb.SetFlag('workspaces-use-clone', True)
            signals = ['fake']
            mc_backgrounds = ['prompt']
            backgrounds = mc_backgrounds

            categories = {
                'llt' : [( 1,'JetToTau')],
            }

            cb.AddObservations(['*'], ['JetTauFR'], ['2018'], ['llt'], categories['llt']) # adding observed data
            cb.AddProcesses(   ['*'], ['JetTauFR'], ['2018'], ['llt'], signals, categories['llt'], True) # adding signals
            cb.AddProcesses(   ['*'], ['JetTauFR'], ['2018'], ['llt'], backgrounds, categories['llt'], False) # adding backgrounds
            
            cb.cp().process(mc_backgrounds+signals).AddSyst(cb,'lumi_2018', 'lnN', ch.SystMap()(1.026))
            #cb.cp().process(['fake']).AddSyst(cb, 'Norm_fake', 'lnN', ch.SystMap()(2.0))
            cb.cp().process(['prompt']).AddSyst(cb, 'Norm_prompt', 'lnN', ch.SystMap()(1.1))
            filepath = os.path.join("/nfs/dust/cms/user/ywen/AZh/workSpacev9/FRworkspace/shapes/pfmt_3_JetToTauFR_%s_%s_%s_%s.root")%(wps,ipt,idm,bin)
            processName = '$BIN/$PROCESS'
            systematicName = '$BIN/$PROCESS_$SYSTEMATIC'
            cb.cp().backgrounds().ExtractShapes(filepath, processName, systematicName)
            cb.cp().signals().ExtractShapes(filepath, processName, systematicName)
            ch.SetStandardBinNames(cb, '$BIN') # Define the name of the category names
            #cb.SetAutoMCStats(cb, 0.0) # Introducing statistical uncertainties on the total background for each histogram bin (Barlow-Beeston lite approach)
            
            bbb = ch.BinByBinFactory()
            bbb.SetAddThreshold(0.1).SetMergeThreshold(0.5).SetFixNorm(True)
            bbb.MergeBinErrors(cb.cp().backgrounds())
            bbb.AddBinByBin(cb.cp().backgrounds(), cb)
            
            datacardPath = 'JetTauFRDatacard/%s_PT%s_DM%s_%s.txt'%(wps,ipt,idm,bin)
            shapePath = 'JetTauFRDatacard/common/%s_PT%s_DM%s_%s.root'%(wps,ipt,idm,bin)
            writer = ch.CardWriter(datacardPath,shapePath)
            writer.SetWildcardMasses([])
            writer.WriteCards('cmb', cb) # writing all datacards into one folder for combination
            #cb.PrintAll()
            #writer.WriteCards(channel, cb.cp().channel([channel])) # writing datacards for each final state in a corresponding folder to be able to perform the measurement individually in each final state
            
            #Fit commands
            text2workspace_cmd = 'text2workspace.py JetTauFRDatacard/%s_PT%s_DM%s_%s.txt -o JetTauFRWorkspace/%s_PT%s_DM%s_%s.root' %(wps,ipt,idm,bin,wps,ipt,idm,bin)
            comine_cmd = 'combine -M FitDiagnostics JetTauFRWorkspace/%s_PT%s_DM%s_%s.root --expectSignal=1 --saveWithUncertainties --saveNormalizations --saveShapes'%(wps,ipt,idm,bin)
            producePostFit_cmd = 'PostFitShapesFromWorkspace -o JetTauFRPostFitShapes/PostFitShapes_%s_PT%s_DM%s_%s.root -f fitDiagnostics.root:fit_s --postfit --sampling --print -d JetTauFRDatacard/%s_PT%s_DM%s_%s.txt -w JetTauFRWorkspace/%s_PT%s_DM%s_%s.root'%(wps,ipt,idm,bin,wps,ipt,idm,bin,wps,ipt,idm,bin)
            
            os.system(text2workspace_cmd)
            os.system(comine_cmd)
            os.system(producePostFit_cmd)
            
            #Get FR
            shapefile = './JetTauFRPostFitShapes/PostFitShapes_%s_PT%s_DM%s_%s.root'%(wps,ipt,idm,bin)
            currentfile = R.TFile.Open(shapefile, "read")
            h_data_prefit = currentfile.Get("JetToTau_prefit/data_obs")
            h_prompt_prefit = currentfile.Get("JetToTau_prefit/prompt")
            h_data_postfit = currentfile.Get("JetToTau_postfit/data_obs")
            h_prompt_postfit = currentfile.Get("JetToTau_postfit/prompt")
            if bin is 'Num':
                Num_prefit_value = h_data_prefit.IntegralAndError(0,1000,Num_data_prefit_err) - h_prompt_prefit.IntegralAndError(0,1000,Num_prompt_prefit_err)
                Num_postfit_value = h_data_postfit.IntegralAndError(0,1000,Num_data_postfit_err) - h_prompt_postfit.IntegralAndError(0,1000,Num_prompt_postfit_err)
            else:
                Den_prefit_value = h_data_prefit.IntegralAndError(0,1000,Den_data_prefit_err) - h_prompt_prefit.IntegralAndError(0,1000,Den_prompt_prefit_err)
                Den_postfit_value = h_data_postfit.IntegralAndError(0,1000,Den_data_postfit_err) - h_prompt_postfit.IntegralAndError(0,1000,Den_prompt_postfit_err)
            currentfile.Close()
            
        FR_prefit_value = R.Double()
        FR_postfit_value = R.Double()
        FR_prefit_error = R.Double()
        FR_postfit_error = R.Double()
        
        FR_prefit_value = Num_prefit_value/Den_prefit_value
        FR_postfit_value = Num_postfit_value/Den_postfit_value
        
        #Binominal distribution error
        #Num_prefit_error = math.sqrt(Num_data_prefit_err*Num_data_prefit_err + Num_prompt_prefit_err*Num_prompt_prefit_err)
        Num_prefit_error = math.sqrt(Den_prefit_value*(1-FR_prefit_value)*FR_prefit_value)
        #Num_postfit_error = math.sqrt(Num_data_postfit_err*Num_data_postfit_err + Num_prompt_postfit_err*Num_prompt_postfit_err)
        Num_postfit_error = math.sqrt(Den_postfit_value*(1-FR_postfit_value)*FR_postfit_value)
        
        Den_prefit_error = math.sqrt(Den_data_prefit_err*Den_data_prefit_err + Den_prompt_prefit_err*Den_prompt_prefit_err)
        Den_postfit_error = math.sqrt(Den_data_postfit_err*Den_data_postfit_err + Den_prompt_postfit_err*Den_prompt_postfit_err)
        FR_prefit_error = math.sqrt(((Num_prefit_error*Num_prefit_error)/(Num_prefit_value*Num_prefit_value)) + ((Den_prefit_error*Den_prefit_error)/(Den_prefit_value*Den_prefit_value)))*FR_prefit_value
        FR_postfit_error = math.sqrt(((Num_postfit_error*Num_postfit_error)/(Num_postfit_value*Num_postfit_value)) + ((Den_postfit_error*Den_postfit_error)/(Den_postfit_value*Den_postfit_value)))*FR_postfit_value
        
        h_FR_prefit.SetBinContent(nptbin,FR_prefit_value)
        h_FR_prefit.SetBinError(nptbin,FR_prefit_error)
        
        h_FR_postfit.SetBinContent(nptbin,FR_postfit_value)
        h_FR_postfit.SetBinError(nptbin,FR_postfit_error)
        
        h2_FR_prefit.SetBinContent(nptbin,ndm,FR_prefit_value)
        h2_FR_prefit.SetBinError(nptbin,ndm,FR_prefit_error)
        h2_FR_prefit.GetYaxis().SetBinLabel(ndm,"DM%s"%idm)
        
        h2_FR_postfit.SetBinContent(nptbin,ndm,FR_postfit_value)
        h2_FR_postfit.SetBinError(nptbin,ndm,FR_postfit_error)
        h2_FR_postfit.GetYaxis().SetBinLabel(ndm,"DM%s"%idm)
        print 'FR prefit:', FR_prefit_value, 'error:', FR_prefit_error
        print 'FR postfit:', FR_postfit_value, 'error:', FR_postfit_error
        nptbin = nptbin + 1
        
    #h_FR_prefit.Draw("e1")
    #h_FR_postfit.Draw("e1same")
    #canvas.SaveAs("prefit.png")
    outfile.cd()
    h_FR_prefit.Write()
    h_FR_postfit.Write()
    ndm = ndm + 1

outfile.cd()
h2_FR_prefit.Write()
h2_FR_postfit.Write()
canvas.cd()
h2_FR_prefit.Draw("COLZ L")
canvas.SaveAs("JetTauFR_%s_Prefit.png"%(wps))
canvas2.cd()
h2_FR_postfit.Draw("COLZ L")
canvas2.SaveAs("JetTauFR_%s_Postfit.png"%(wps))
outfile.Close()

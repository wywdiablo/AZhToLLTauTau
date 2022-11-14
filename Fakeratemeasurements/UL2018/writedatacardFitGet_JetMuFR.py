#! /usr/bin/env python
# Author: Yiwen Wen (Oct 2021)
# Description: Create datacards for combine
import sys, os,array
from CombineHarvester.CombineTools import ch
import ROOT as R
import math
from array import array
R.PyConfig.IgnoreCommandLineOptions = True
R.gROOT.SetBatch(R.kTRUE)
R.gStyle.SetOptStat(0)
canvas = R.TCanvas()
canvas2 = R.TCanvas()

wps = 'Medium_Iso0p15'
pts = ['10to20','20to30','30to40','40to60','Gt60']
etas = ['barrel','endcap']
bins = ['Num','Den']
netabins = 2
etabin = [0,1,2]
nptbins = 5
ptbins = [10,20,30,40,60,80]
h2_FR_prefit = R.TH2D("JetToMuPrefitFR", "Jet To #mu Prefit Fake Rate",nptbins,array('d',ptbins),netabins,array('d',etabin))
h2_FR_postfit = R.TH2D("JetToMuPostfitFR", "Jet To #mu Postfit Fake Rate",nptbins,array('d',ptbins),netabins,array('d',etabin))
outfile = R.TFile("JetMuFakeRate_%s.root"%(wps), "RECREATE")
neta = 1
print '<<<<<<< Numerator working point  : ', wps
for ieta in etas :
    print '<<<<<<<<<<<<<< eta  : ', ieta
    #gr_err = R.TGraphAsymmErrors(nptbins)
    h_FR_prefit = R.TH1D("PrefitFR_%s"%ieta,"Prefit Fake Rate %s"%ieta, nptbins, array('d',ptbins))
    h_FR_postfit = R.TH1D("PostfitFR_%s"%ieta,"Postfit Fake Rate %s"%ieta, nptbins, array('d',ptbins))
    nptbin = 1
    for ipt in pts :
        print '<<<<<<< pt range  : ', ipt
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
            #looping numerator and demoninator
            print '<<<<<<<<<<<<<< bin  : ', bin
            cb = ch.CombineHarvester()
            #cb.SetFlag('workspaces-use-clone', True)
            signals = ['fake']
            mc_backgrounds = ['prompt']
            backgrounds = mc_backgrounds
            all = signals + backgrounds

            categories = {
                'llm' : [( 1,'JetToMu')],
            }

            cb.AddObservations(['*'], ['JetMuFR'], ['2018'], ['llm'], categories['llm']) # adding observed data
            cb.AddProcesses(   ['*'], ['JetMuFR'], ['2018'], ['llm'], backgrounds, categories['llm'], False) # adding backgrounds
            cb.AddProcesses(   ['*'], ['JetMuFR'], ['2018'], ['llm'], signals, categories['llm'], True) # adding signal

            cb.cp().process(all).AddSyst(cb,'lumi_2018', 'lnN', ch.SystMap()(1.026))
            #cb.cp().process(['fake']).AddSyst(cb, 'Norm_fake', 'lnN', ch.SystMap()(2.0))
            cb.cp().process(['prompt']).AddSyst(cb, 'Norm_prompt', 'lnN', ch.SystMap()(1.05))
            filepath = os.path.join("/nfs/dust/cms/user/ywen/AZh/workSpacev9/FRworkspace/shapes/pfmt_3_JetToMuFR_%s_%s_%s_%s.root")%(wps,ipt,ieta,bin)
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
            
            datacardPath = 'JetMuFRDatacard/%s_PT%s_%s_%s.txt'%(wps,ipt,ieta,bin)
            shapePath = 'JetMuFRDatacard/common/%s_PT%s_%s_%s.root'%(wps,ipt,ieta,bin)
            writer = ch.CardWriter(datacardPath,shapePath)
            writer.SetWildcardMasses([])
            writer.WriteCards('cmb', cb) # writing all datacards into one folder for combination
            #cb.PrintAll()
            #writer.WriteCards(channel, cb.cp().channel([channel])) # writing datacards for each final state in a corresponding folder to be able to perform the measurement individually in each final state
                        
            ###Fitting started
            text2workspace_cmd = 'text2workspace.py JetMuFRDatacard/%s_PT%s_%s_%s.txt -o JetMuFRWorkspace/%s_PT%s_%s_%s.root' %(wps,ipt,ieta,bin,wps,ipt,ieta,bin)
            comine_cmd = 'combine -M FitDiagnostics JetMuFRWorkspace/%s_PT%s_%s_%s.root --expectSignal=1 --saveWithUncertainties --saveNormalizations --saveShapes'%(wps,ipt,ieta,bin)
            producePostFit_cmd = 'PostFitShapesFromWorkspace -o JetMuFRPostFitShapes/PostFitShapes_%s_PT%s_%s_%s.root -f fitDiagnostics.root:fit_s --postfit --sampling --print -d JetMuFRDatacard/%s_PT%s_%s_%s.txt -w JetMuFRWorkspace/%s_PT%s_%s_%s.root'%(wps,ipt,ieta,bin,wps,ipt,ieta,bin,wps,ipt,ieta,bin)
            
            os.system(text2workspace_cmd)
            os.system(comine_cmd)
            os.system(producePostFit_cmd)
            ###Getting FR started
            shapefile = './JetMuFRPostFitShapes/PostFitShapes_%s_PT%s_%s_%s.root'%(wps,ipt,ieta,bin)
            currentfile = R.TFile.Open(shapefile, "read")
            h_data_prefit = currentfile.Get("JetToMu_prefit/data_obs")
            h_prompt_prefit = currentfile.Get("JetToMu_prefit/prompt")
            h_data_postfit = currentfile.Get("JetToMu_postfit/data_obs")
            h_prompt_postfit = currentfile.Get("JetToMu_postfit/prompt")
            if bin is 'Num':
                Num_prefit_value = h_data_prefit.IntegralAndError(0,1000,Num_data_prefit_err) - h_prompt_prefit.IntegralAndError(0,1000,Num_prompt_prefit_err)
                Num_postfit_value = h_data_postfit.IntegralAndError(0,1000,Num_data_postfit_err) - h_prompt_postfit.IntegralAndError(0,1000,Num_prompt_postfit_err)
            else:
                Den_prefit_value = h_data_prefit.IntegralAndError(0,1000,Den_data_prefit_err) - h_prompt_prefit.IntegralAndError(0,1000,Den_prompt_prefit_err)
                Den_postfit_value = h_data_postfit.IntegralAndError(0,1000,Den_data_postfit_err) - h_prompt_postfit.IntegralAndError(0,1000,Den_prompt_postfit_err)
            currentfile.Close()
            #End looping numerator and demoninator
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
        
        h2_FR_prefit.SetBinContent(nptbin,neta,FR_prefit_value)
        h2_FR_prefit.SetBinError(nptbin,neta,FR_prefit_error)
        h2_FR_prefit.GetYaxis().SetBinLabel(neta,"%s"%ieta)
        
        h2_FR_postfit.SetBinContent(nptbin,neta,FR_postfit_value)
        h2_FR_postfit.SetBinError(nptbin,neta,FR_postfit_error)
        h2_FR_postfit.GetYaxis().SetBinLabel(neta,"%s"%ieta)
        print 'FR prefit:', FR_prefit_value, 'error:', FR_prefit_error
        print 'FR postfit:', FR_postfit_value, 'error:', FR_postfit_error
        nptbin = nptbin + 1
        
    #h_FR_prefit.Draw("e1")
    #h_FR_postfit.Draw("e1same")
    #canvas.SaveAs("prefit.png")
    outfile.cd()
    h_FR_prefit.Write()
    h_FR_postfit.Write()
    neta = neta + 1

outfile.cd()
h2_FR_prefit.Write()
h2_FR_postfit.Write()
canvas.cd()
h2_FR_prefit.Draw("COLZ L")
canvas.SaveAs("JetMuFR_%s_Prefit.png"%(wps))
canvas2.cd()
h2_FR_postfit.Draw("COLZ L")
canvas2.SaveAs("JetMuFR_%s_Postfit.png"%(wps))
outfile.Close()

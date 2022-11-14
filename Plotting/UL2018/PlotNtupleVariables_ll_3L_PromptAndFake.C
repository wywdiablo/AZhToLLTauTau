//----------------------Version 1.0-------------------------//
//Ploting Ntuple variables for AZh to lltautau analysis
//Author: Yiwen Wen
//DESY
//----------------------------------------------------------//
//v1.0
#include "HttStylesNew.cc"
#include "HtoH.h"
#include "TColor.h"
#include "CMS_lumi.C"
#include "Samples.C"
void PlotNtupleVariables_ll_3L_PromptAndFake(
                                 TString varName= "pfmt_3",
                                 TString xtitle = "Transverse Mass [GeV]",
                                 TString ytitle = "Events",
                                 float xLower = 0,
                                 float xUpper = 300,
                                 int numberofbins = 15,
                                 bool logY = false,
                                 bool legLeft = false,
                                 TString finalState = "e",
                                 TString VSjetWP = "Medium",
                                 TString VSeWP = "VLoose",
                                 TString VSmuWP = "VLoose",
                                 bool getDenominator = true,
                                 TString extraCuts = "",
                                 TString tags = "")
{
	SetStyle();
    //gStyle->SetOptStat(1111);


	double xsec[34] = { //in pico barn (pb)
        1, // (0)data
        6077.22,// (1)Drell-Yan Z
        0.00319*1.7,// (2)GGContinToZZTo2e2mu
        0.00319*1.7,// (3)
        0.00319*1.7,// (4)
        0.00159*1.7,// (5)
        0.00159*1.7,//(6)
        0.00159*1.7,//(7)
        1.256*1.1,//(8)
        3.22,//(9)ZZTo2L2Q
        0.78,//(10)TTZ
        0.2043,//(11)TTZW
        88.29, //(12)TTTo2L2Nu
        377.96,//(13)TT to hadronic
        365.35, //(14)TT to semileptonic
        0.2086,
        0.1651,
        0.01398,
        0.05565,
        4.708,
        5.595,//(20) WZTo2L2Q
        48.58*0.0627,
        3.782*0.0627,
        0.5328*0.0627,
        0.8400*0.0627,
        0.884*0.0627,
        0.507*0.0627,
        48.58*0.2137*0.1061,
        3.925*0.2137*0.1061,
        0.114,//HW- H-> WW
        0.180,//HW+ H-> WW
        0.163,//HZJ H -> WW
        0.0262*0.1,//GGZH H -> WW
        0.0128//GG H -> ZZ -> 4L
    };
    
        float lumi = 59740;
        //Deal with bins and binning
        float xMin = xLower;
        float xMax = xUpper;
        int nBins = numberofbins;
        float bins[100];
        float binWidth = (xMax-xMin)/float(nBins);
            for (int iB=0; iB<=nBins; ++iB)
                bins[iB] = xMin + float(iB)*binWidth;

        int nSamples = 34;
        TH1D * hist_prompt[34];   //This is "hist_mumu"
        TH1D * hist_ee_prompt[34];

        TH1D * hist_fake[34];   //This is "hist_mumu"
        TH1D * hist_ee_fake[34];
    
        //inistiating cuts
        TString cuts[34];
        TString ZeeCuts = "(electronIDWP90v2_1>0.5 && iso_1<0.15 && electronIDWP90v2_2>0.5 && iso_2<0.15)*";
        TString ZmumuCuts = "(muonLoose_1>0.5 && iso_1<0.15 && muonLoose_2>0.5 && iso_2<0.15)*";
        
        TString MCCommonWeight = "IdSF_1*IdSF_2*TrigSF_1*TrigSF_2*puweight*mcweight*";
        TString PromptCut = "(gen_match_3!=6)*";
        TString FakeCut = "(gen_match_3==6)*";
    
        TString NumeratorCuts = "";
        TString DenominatorCuts = "";
        TString FRCuts = "";

        if(finalState == "t")
        {
            if(VSjetWP == "Medium" && VSeWP == "VLoose" && VSmuWP == "Tight")
            {
                NumeratorCuts = "(VSjetMedium_3>0.5&&VSeVLoose_3>0.5&&VSmuTight_3>0.5)";
            }
            if(VSjetWP == "Medium" && VSeWP == "Tight" && VSmuWP == "VLoose")
            {
                NumeratorCuts = "(VSjetMedium_3>0.5&&VSeTight_3>0.5&&VSmuVLoose_3>0.5)";
            }
            if(VSjetWP == "Medium" && VSeWP == "VLoose" && VSmuWP == "VLoose")
            {
                NumeratorCuts = "(VSjetMedium_3>0.5&&VSeVLoose_3>0.5&&VSmuVLoose_3>0.5)";
            }
            
            DenominatorCuts = "(VSjetVVVLoose_3>0.5)";
        }
    
        if(finalState == "e")
        {
            NumeratorCuts = "(electronIDWP90v2_3>0.5&&iso_3<0.15)";
            //NumeratorCuts = "(electronIDWP80v2_3>0.5&&iso_3<0.15)";
            DenominatorCuts = "(iso_3<9999)";
        }
        if(finalState == "m")
        {
            //NumeratorCuts = "(muonLoose_3>0.5&&iso_3<0.15)";
            NumeratorCuts = "(muonMedium_3>0.5&&iso_3<0.15)";
            DenominatorCuts = "(iso_3<9999)";
        }
        if(finalState == "em")
        {
            MCCommonWeight = "IdSF_1*IdSF_2*IdSF_3*IdSF_4*TrigSF_1*TrigSF_2*puweight*mcweight*";
            NumeratorCuts = "(electronIDWP90v2_3>0.5 && iso_3<0.15 && muonLoose_4>0.5 && iso_4<0.15&&H_SS<0.5)";
            DenominatorCuts = "(H_SS<0.5)";
        }
        if(finalState == "et")
        {
            MCCommonWeight = "IdSF_1*IdSF_2*IdSF_3*tauVSjetSF_4*tauVSeSF_4*tauVSmuSF_4*TrigSF_1*TrigSF_2*puweight*mcweight*";
            NumeratorCuts = "(electronIDWP90v2_3>0.5 && iso_3<0.15 && H_SS<0.5)";
            DenominatorCuts = "(H_SS<0.5)";
        }
        if(finalState == "mt")
        {
            MCCommonWeight = "IdSF_1*IdSF_2*IdSF_3*tauVSjetSF_4*tauVSeSF_4*tauVSmuSF_4*TrigSF_1*TrigSF_2*puweight*mcweight*";
            NumeratorCuts = "(muonLoose_3>0.5 && iso_3<0.15 && H_SS<0.5)";
            DenominatorCuts = "(H_SS<0.5)";
        }
        if(finalState == "tt")
        {
            MCCommonWeight = "IdSF_1*IdSF_2*tauVSjetSF_3*tauVSeSF_3*tauVSmuSF_3*tauVSjetSF_4*tauVSeSF_4*tauVSmuSF_4*TrigSF_1*TrigSF_2*puweight*mcweight*";
            NumeratorCuts = "(H_SS<0.5)";
            DenominatorCuts = "(H_SS<0.5)";
        }

        if(getDenominator)
        {
            FRCuts = DenominatorCuts;
        }
        else
        {
            FRCuts = NumeratorCuts;
        }

            for (int i=0; i<nSamples; ++i)
            {
                cuts[i] = MCCommonWeight + FRCuts;
            }
            cuts[0] = FRCuts;
    
            for (int i=0; i<nSamples; ++i)
            {
                //ee
                std::cout << "EE channel " <<i<<":"+samples_ee_UL2018_FR[i]<<std::endl;
                TFile * file_ee = new TFile("./NTuples_ee"+finalState+"/"+samples_ee_UL2018_FR[i]+".root");
                TH1D * histWeightsH_ee = (TH1D*)file_ee->Get("histWeightsH");
                TTree * tree_ee = (TTree*)file_ee->Get("SynTree");
                double normaliza_ee = xsec[i]*lumi/histWeightsH_ee->GetSumOfWeights();
                TString histName_ee_prompt = samples_ee_UL2018_FR[i] + "_"+varName + "_ee_prompt";
                hist_ee_prompt[i] = new TH1D(histName_ee_prompt,"",nBins,xMin,xMax);
                hist_ee_prompt[i]->Sumw2();
                if(i > 0)
                {
                    tree_ee->Draw(varName+">>"+histName_ee_prompt,PromptCut+ZeeCuts+cuts[i]+extraCuts);

                }
                else
                {
                    tree_ee->Draw(varName+">>"+histName_ee_prompt,ZeeCuts+cuts[i]+extraCuts);
                }
                if(i > 0)
                {
                    for (int iB=1; iB<=nBins; ++iB)
                    {
                        double x = hist_ee_prompt[i]->GetBinContent(iB);
                        double e = hist_ee_prompt[i]->GetBinError(iB);
                        hist_ee_prompt[i]->SetBinContent(iB,normaliza_ee*x);
                        hist_ee_prompt[i]->SetBinError(iB,normaliza_ee*e);
                    }
                }
                std::cout << hist_ee_prompt[i]->GetSumOfWeights() << std::endl;
                
                TString histName_ee_fake = samples_ee_UL2018_FR[i] + "_"+varName + "_ee_fake";
                hist_ee_fake[i] = new TH1D(histName_ee_fake,"",nBins,xMin,xMax);
                hist_ee_fake[i]->Sumw2();
                if(i > 0)
                {
                        tree_ee->Draw(varName+">>"+histName_ee_fake,FakeCut+ZeeCuts+cuts[i]+extraCuts);
                        for (int iB=1; iB<=nBins; ++iB)
                        {
                            double x = hist_ee_fake[i]->GetBinContent(iB);
                            double e = hist_ee_fake[i]->GetBinError(iB);
                            hist_ee_fake[i]->SetBinContent(iB,normaliza_ee*x);
                            hist_ee_fake[i]->SetBinError(iB,normaliza_ee*e);
                        }
                        std::cout << hist_ee_fake[i]->GetSumOfWeights() << std::endl;
                }
                //mumu
                std::cout << "MuMu channel " <<i<<":"+samples_mumu_UL2018_FR[i]<<std::endl;
                TFile * file_mumu = new TFile("./NTuples_mm"+finalState+"/"+samples_mumu_UL2018_FR[i]+".root");
                TH1D * histWeightsH_mumu = (TH1D*)file_mumu->Get("histWeightsH");
                TTree * tree_mumu = (TTree*)file_mumu->Get("SynTree");
                double normaliza_mumu = xsec[i]*lumi/histWeightsH_mumu->GetSumOfWeights();
                TString histName_mumu_prompt = samples_mumu_UL2018_FR[i] + "_"+varName+ "_mumu_prompt";
                hist_prompt[i] = new TH1D(histName_mumu_prompt,"",nBins,xMin,xMax);
                hist_prompt[i]->Sumw2();
                if(i > 0)
                {
                    tree_mumu->Draw(varName+">>"+histName_mumu_prompt,PromptCut+ZmumuCuts+cuts[i]+extraCuts);

                }
                else
                {
                    tree_mumu->Draw(varName+">>"+histName_mumu_prompt,ZmumuCuts+cuts[i]+extraCuts);
                }
                if(i > 0)
                {
                    for (int iB=1; iB<=nBins; ++iB)
                    {
                        double x = hist_prompt[i]->GetBinContent(iB);
                        double e = hist_prompt[i]->GetBinError(iB);
                        hist_prompt[i]->SetBinContent(iB,normaliza_mumu*x);
                        hist_prompt[i]->SetBinError(iB,normaliza_mumu*e);
                    }
                }
                std::cout << hist_prompt[i]->GetSumOfWeights() << std::endl;
                hist_prompt[i]->Add(hist_prompt[i],hist_ee_prompt[i]);
                
                TString histName_mumu_fake = samples_mumu_UL2018_FR[i] + "_"+varName+ "_mumu_fake";
                hist_fake[i] = new TH1D(histName_mumu_fake,"",nBins,xMin,xMax);
                hist_fake[i]->Sumw2();
                if(i > 0)
                {
                    tree_mumu->Draw(varName+">>"+histName_mumu_fake,FakeCut+ZmumuCuts+cuts[i]+extraCuts);
                    for (int iB=1; iB<=nBins; ++iB)
                    {
                        double x = hist_fake[i]->GetBinContent(iB);
                        double e = hist_fake[i]->GetBinError(iB);
                        hist_fake[i]->SetBinContent(iB,normaliza_mumu*x);
                        hist_fake[i]->SetBinError(iB,normaliza_mumu*e);
                    }
                    std::cout << hist_fake[i]->GetSumOfWeights() << std::endl;
                    hist_fake[i]->Add(hist_fake[i],hist_ee_fake[i]);
                }
            }
    
    // *******************************
    // ***** DY+Jets samples *******
    // *******************************
    
    TH1D * histZ_prompt[9];
    TH1D * histZ_ee_prompt[9];
    
    TH1D * histZ_fake[9];
    TH1D * histZ_ee_fake[9];

    double dyRefXSec[5] = {6077.22,
        1.137*877.8,
        1.137*304.4,
        1.137*111.5,
        1.137*44.05};//2018 values
    double dyRefEvents_ee[5];
    double dyRefEvents_mumu[5];

    for (int iDY=0; iDY<5; ++iDY) {
        TFile * file_ee = new TFile("./NTuples_ee"+finalState+"/"+dyRefSamples_UL2018[iDY]+".root");
        TH1D * histWeightsH_ee = (TH1D*)file_ee->Get("histWeightsH");
        dyRefEvents_ee[iDY] = histWeightsH_ee->GetSumOfWeights();
        
        TFile * file_mumu = new TFile("./NTuples_mm"+finalState+"/"+dyRefSamples_UL2018[iDY]+".root");
        TH1D * histWeightsH_mumu = (TH1D*)file_mumu->Get("histWeightsH");
        dyRefEvents_mumu[iDY] = histWeightsH_mumu->GetSumOfWeights();
    }

    double dyNorm_ee[9];
    dyNorm_ee[0] = lumi*dyRefXSec[0]/dyRefEvents_ee[0];
    dyNorm_ee[1] = lumi/(dyRefEvents_ee[0]/dyRefXSec[0]+dyRefEvents_ee[1]/dyRefXSec[1]);
    dyNorm_ee[2] = lumi/(dyRefEvents_ee[0]/dyRefXSec[0]+dyRefEvents_ee[2]/dyRefXSec[2]);
    dyNorm_ee[3] = lumi/(dyRefEvents_ee[0]/dyRefXSec[0]+dyRefEvents_ee[3]/dyRefXSec[3]);
    dyNorm_ee[4] = lumi/(dyRefEvents_ee[0]/dyRefXSec[0]+dyRefEvents_ee[4]/dyRefXSec[4]);
    dyNorm_ee[5] = lumi/(dyRefEvents_ee[0]/dyRefXSec[0]+dyRefEvents_ee[1]/dyRefXSec[1]);
    dyNorm_ee[6] = lumi/(dyRefEvents_ee[0]/dyRefXSec[0]+dyRefEvents_ee[2]/dyRefXSec[2]);
    dyNorm_ee[7] = lumi/(dyRefEvents_ee[0]/dyRefXSec[0]+dyRefEvents_ee[3]/dyRefXSec[3]);
    dyNorm_ee[8] = lumi/(dyRefEvents_ee[0]/dyRefXSec[0]+dyRefEvents_ee[4]/dyRefXSec[4]);
    
    double dyNorm_mumu[9];
    dyNorm_mumu[0] = lumi*dyRefXSec[0]/dyRefEvents_mumu[0];
    dyNorm_mumu[1] = lumi/(dyRefEvents_mumu[0]/dyRefXSec[0]+dyRefEvents_mumu[1]/dyRefXSec[1]);
    dyNorm_mumu[2] = lumi/(dyRefEvents_mumu[0]/dyRefXSec[0]+dyRefEvents_mumu[2]/dyRefXSec[2]);
    dyNorm_mumu[3] = lumi/(dyRefEvents_mumu[0]/dyRefXSec[0]+dyRefEvents_mumu[3]/dyRefXSec[3]);
    dyNorm_mumu[4] = lumi/(dyRefEvents_mumu[0]/dyRefXSec[0]+dyRefEvents_mumu[4]/dyRefXSec[4]);
    dyNorm_mumu[5] = lumi/(dyRefEvents_mumu[0]/dyRefXSec[0]+dyRefEvents_mumu[1]/dyRefXSec[1]);
    dyNorm_mumu[6] = lumi/(dyRefEvents_mumu[0]/dyRefXSec[0]+dyRefEvents_mumu[2]/dyRefXSec[2]);
    dyNorm_mumu[7] = lumi/(dyRefEvents_mumu[0]/dyRefXSec[0]+dyRefEvents_mumu[3]/dyRefXSec[3]);
    dyNorm_mumu[8] = lumi/(dyRefEvents_mumu[0]/dyRefXSec[0]+dyRefEvents_mumu[4]/dyRefXSec[4]);
    
    
    TString npartonCuts[9] = {"*(npartons==0||npartons>4)",
        "*(npartons==1)",
        "*(npartons==2)",
        "*(npartons==3)",
        "*(npartons==4)",
        "",
        "",
        "",
        ""
    };
    TString cutsZ[9];

        for (int iDY=0; iDY<9; ++iDY)
        {
            cutsZ[iDY]   = MCCommonWeight + "zptweight*"+FRCuts+npartonCuts[iDY];
        }
    
    // filling histograms for DYJets samples
    for (int i=0; i<9; ++i) { // run over samples
        TFile * file_ee = new TFile("./NTuples_ee"+finalState+"/"+dySampleNames_UL2018[i]+".root");
        TTree * tree_ee = (TTree*)file_ee->Get("SynTree");
        double normaliza_ee = dyNorm_ee[i];
        
        TString histNameZ_ee_prompt = dySampleNames_UL2018[i] + "_"+varName+"_Z_ee_prompt";
        histZ_ee_prompt[i] = new TH1D(histNameZ_ee_prompt,"",nBins,xMin,xMax);
        histZ_ee_prompt[i]->Sumw2();
        tree_ee->Draw(varName+">>"+histNameZ_ee_prompt, PromptCut + ZeeCuts+cutsZ[i] + extraCuts);

        for (int iB=1; iB<=nBins; ++iB) {
            
            double x = histZ_ee_prompt[i]->GetBinContent(iB);
            double e = histZ_ee_prompt[i]->GetBinError(iB);
            histZ_ee_prompt[i]->SetBinContent(iB,normaliza_ee*x);
            histZ_ee_prompt[i]->SetBinError(iB,normaliza_ee*e);
        }
        std::cout << dySampleNames_UL2018[i] << " -> DY Z = " << histZ_ee_prompt[i]->GetEntries() << " : " << histZ_ee_prompt[i]->GetSumOfWeights() << std::endl;
        
        TString histNameZ_ee_fake = dySampleNames_UL2018[i] + "_"+varName+"_Z_ee_fake";
        histZ_ee_fake[i] = new TH1D(histNameZ_ee_fake,"",nBins,xMin,xMax);
        histZ_ee_fake[i]->Sumw2();
        tree_ee->Draw(varName+">>"+histNameZ_ee_fake, FakeCut + ZeeCuts+cutsZ[i] + extraCuts);

        for (int iB=1; iB<=nBins; ++iB) {
            
            double x = histZ_ee_fake[i]->GetBinContent(iB);
            double e = histZ_ee_fake[i]->GetBinError(iB);
            histZ_ee_fake[i]->SetBinContent(iB,normaliza_ee*x);
            histZ_ee_fake[i]->SetBinError(iB,normaliza_ee*e);
        }
        std::cout << dySampleNames_UL2018[i] << " -> DY Z = " << histZ_ee_fake[i]->GetEntries() << " : " << histZ_ee_fake[i]->GetSumOfWeights() << std::endl;
        
        
        TFile * file_mumu = new TFile("./NTuples_mm"+finalState+"/"+dySampleNames_UL2018[i]+".root");
        TTree * tree_mumu = (TTree*)file_mumu->Get("SynTree");
        double normaliza_mumu = dyNorm_mumu[i];
        
        TString histNameZ_mumu_prompt = dySampleNames_UL2018[i] + "_"+varName+"_Z_mumu_prompt";
        histZ_prompt[i] = new TH1D(histNameZ_mumu_prompt,"",nBins,xMin,xMax);
        histZ_prompt[i]->Sumw2();
        tree_mumu->Draw(varName+">>"+histNameZ_mumu_prompt, PromptCut + ZmumuCuts+cutsZ[i] + extraCuts);
        for (int iB=1; iB<=nBins; ++iB) {
            double x = histZ_prompt[i]->GetBinContent(iB);
            double e = histZ_prompt[i]->GetBinError(iB);
            histZ_prompt[i]->SetBinContent(iB,normaliza_mumu*x);
            histZ_prompt[i]->SetBinError(iB,normaliza_mumu*e);
        }
        std::cout << dySampleNames_UL2018[i] << " -> DY Z = " << histZ_prompt[i]->GetEntries() << " : " << histZ_prompt[i]->GetSumOfWeights() << std::endl;
        histZ_prompt[i]->Add(histZ_prompt[i],histZ_ee_prompt[i]);
        
        TString histNameZ_mumu_fake = dySampleNames_UL2018[i] + "_"+varName+"_Z_mumu_fake";
        histZ_fake[i] = new TH1D(histNameZ_mumu_fake,"",nBins,xMin,xMax);
        histZ_fake[i]->Sumw2();
        tree_mumu->Draw(varName+">>"+histNameZ_mumu_fake, FakeCut + ZmumuCuts+cutsZ[i] + extraCuts);
        for (int iB=1; iB<=nBins; ++iB) {
            double x = histZ_fake[i]->GetBinContent(iB);
            double e = histZ_fake[i]->GetBinError(iB);
            histZ_fake[i]->SetBinContent(iB,normaliza_mumu*x);
            histZ_fake[i]->SetBinError(iB,normaliza_mumu*e);
        }
        std::cout << dySampleNames_UL2018[i] << " -> DY Z = " << histZ_fake[i]->GetEntries() << " : " << histZ_fake[i]->GetSumOfWeights() << std::endl;
        histZ_fake[i]->Add(histZ_fake[i],histZ_ee_fake[i]);
    }
    
    hist_prompt[1] = histZ_prompt[0];
    hist_fake[1] = histZ_fake[0];
    std::cout << " -> DY Adding = " << hist_prompt[1]->GetEntries() << " : " << hist_prompt[1]->GetSumOfWeights() << std::endl;
    std::cout << " -> DY Adding = " << hist_fake[1]->GetEntries() << " : " << hist_fake[1]->GetSumOfWeights() << std::endl;

    for (int iDY=1; iDY<9; ++iDY)
    {
        hist_prompt[1]->Add(hist_prompt[1],histZ_prompt[iDY]);
        std::cout << " -> DY Adding = " << hist_prompt[1]->GetEntries() << " : " << hist_prompt[1]->GetSumOfWeights() << std::endl;
        
        hist_fake[1]->Add(hist_fake[1],histZ_fake[iDY]);
        std::cout << " -> DY Adding = " << hist_fake[1]->GetEntries() << " : " << hist_fake[1]->GetSumOfWeights() << std::endl;
    }
    
    std::cout << "end: DY stitching" << std::endl;
    
        //  adding up ZZ continue backgrounds
        for (int iH=3; iH<10; ++iH)
        {
            hist_prompt[2]->Add(hist_prompt[2],hist_prompt[iH]);
            hist_fake[2]->Add(hist_fake[2],hist_fake[iH]);
        }
        
        // adding up TT backgrounds
        for (int iH=11; iH<15; ++iH)
        {
            hist_prompt[10]->Add(hist_prompt[10],hist_prompt[iH]);
            hist_fake[10]->Add(hist_fake[10],hist_fake[iH]);
        }
    
        // adding up VVV backgrounds
        for (int iH=16; iH<19; ++iH)
        {
            hist_prompt[15]->Add(hist_prompt[15],hist_prompt[iH]);
            hist_fake[15]->Add(hist_fake[15],hist_fake[iH]);
        }
    
        // (18+19) is WZ
        for (int iH=20; iH<21; ++iH)
        {
            hist_prompt[19]->Add(hist_prompt[19],hist_prompt[iH]);
            hist_fake[19]->Add(hist_fake[19],hist_fake[iH]);
        }
    
        // adding up SM H 125 backgrounds
        for(int iH=22; iH<34; ++iH)
        {
            hist_prompt[21]->Add(hist_prompt[21],hist_prompt[iH]);
            hist_fake[21]->Add(hist_fake[21],hist_fake[iH]);
        }

        TH1D * data_obs = (TH1D*)hist_prompt[0]->Clone("data_obs");
        TH1D * DY_prompt = (TH1D*)hist_prompt[1]->Clone("DY_prompt");
        TH1D * ZZ_prompt = (TH1D*)hist_prompt[2]->Clone("ZZ_prompt");
        TH1D * TT_prompt = (TH1D*)hist_prompt[10]->Clone("TT_prompt");
        TH1D * VVV_prompt = (TH1D*)hist_prompt[15]->Clone("VVV_prompt");
        TH1D * WZ_prompt = (TH1D*)hist_prompt[19]->Clone("WZ_prompt");
        TH1D * SMH_prompt = (TH1D*)hist_prompt[21]->Clone("SMH_prompt");
        float totData = data_obs->GetSumOfWeights();
        float totDY_prompt = DY_prompt->GetSumOfWeights();
        float totZZ_prompt = ZZ_prompt->GetSumOfWeights();
        float totTT_prompt = TT_prompt->GetSumOfWeights();
        float totVVV_prompt = VVV_prompt->GetSumOfWeights();
        float totWZ_prompt = WZ_prompt->GetSumOfWeights();
        float totSMH_prompt = SMH_prompt->GetSumOfWeights();
        float totBkg_prompt = totDY_prompt+totZZ_prompt+totTT_prompt+totVVV_prompt+totWZ_prompt+totSMH_prompt;

        std::cout << "data : " << totData << std::endl;
        std::cout << "Total bkg prompt: " << totBkg_prompt << std::endl;
        std::cout << "DY prompt : " << DY_prompt->GetSumOfWeights() << std::endl;
        std::cout << "ZZ prompt  : " << ZZ_prompt->GetSumOfWeights() << std::endl;
        std::cout << "TT prompt  : " << TT_prompt->GetSumOfWeights() << std::endl;
        std::cout << "VVV prompt : " << VVV_prompt->GetSumOfWeights() << std::endl;
        std::cout << "WZ prompt : " << WZ_prompt->GetSumOfWeights() << std::endl;
        std::cout << "SM H prompt : " << SMH_prompt->GetSumOfWeights() << std::endl;
    
        VVV_prompt->Add(VVV_prompt,WZ_prompt);
        TT_prompt->Add(TT_prompt,VVV_prompt);
        ZZ_prompt->Add(ZZ_prompt,TT_prompt);
        SMH_prompt->Add(SMH_prompt,ZZ_prompt);
        DY_prompt->Add(DY_prompt,SMH_prompt);
    
        TH1D * DY_fake = (TH1D*)hist_fake[1]->Clone("DY_fake");
        TH1D * ZZ_fake = (TH1D*)hist_fake[2]->Clone("ZZ_fake");
        TH1D * TT_fake = (TH1D*)hist_fake[10]->Clone("TT_fake");
        TH1D * VVV_fake = (TH1D*)hist_fake[15]->Clone("VVV_fake");
        TH1D * WZ_fake = (TH1D*)hist_fake[19]->Clone("WZ_fake");
        TH1D * SMH_fake = (TH1D*)hist_fake[21]->Clone("SMH_fake");
        float totDY_fake = DY_fake->GetSumOfWeights();
        float totZZ_fake = ZZ_fake->GetSumOfWeights();
        float totTT_fake = TT_fake->GetSumOfWeights();
        float totVVV_fake = VVV_fake->GetSumOfWeights();
        float totWZ_fake = WZ_fake->GetSumOfWeights();
        float totSMH_fake = SMH_fake->GetSumOfWeights();
        float totBkg_fake = totDY_fake+totZZ_fake+totTT_fake+totVVV_fake+totWZ_fake+totSMH_fake;

        std::cout << "Total bkg fake: " << totBkg_fake << std::endl;
        std::cout << "DY fake : " << DY_fake->GetSumOfWeights() << std::endl;
        std::cout << "ZZ fake  : " << ZZ_fake->GetSumOfWeights() << std::endl;
        std::cout << "TT fake  : " << TT_fake->GetSumOfWeights() << std::endl;
        std::cout << "VVV fake : " << VVV_fake->GetSumOfWeights() << std::endl;
        std::cout << "WZ fake : " << WZ_fake->GetSumOfWeights() << std::endl;
        std::cout << "SM H fake : " << SMH_fake->GetSumOfWeights() << std::endl;
    
        VVV_fake->Add(VVV_fake,WZ_fake);
        TT_fake->Add(TT_fake,VVV_fake);
        ZZ_fake->Add(ZZ_fake,TT_fake);
        SMH_fake->Add(SMH_fake,ZZ_fake);
        DY_fake->Add(DY_fake,SMH_fake);
    
        TString finalStateTag = "";
        if(finalState=="t")
        {
            finalStateTag = "Tau";
        }
        if(finalState=="m")
        {
            finalStateTag = "Mu";
        }
        if(finalState=="e")
        {
            finalStateTag = "Ele";
        }
        TFile * file = new TFile("./FRworkspace/shapes/"+varName+"_JetTo"+finalStateTag+"FR_"+tags+".root","recreate");
        file->mkdir("JetTo"+finalStateTag);
        file->cd("JetTo"+finalStateTag);
        TH1D * fake = (TH1D*)DY_fake->Clone("fake");
        TH1D * prompt = (TH1D*)DY_prompt->Clone("prompt");
        TH1D * data = (TH1D*)data_obs->Clone("data_obs");
        file->Write();
        file->Close();

        DY_fake->Add(DY_fake,DY_prompt);

        TH1D * dummy = (TH1D*)DY_fake->Clone("dummy");
        float errLumi = 0.03;
        float errDY=0.1;
        float errVVV = 0.15;
        float errTT = 0.15;
        float errWZ = 0.15;

        for (int iB=1; iB<=nBins; ++iB)
        {
            float eVVV = errVVV*(VVV_prompt->GetBinContent(iB)+VVV_fake->GetBinContent(iB));
            float eDY = errDY*(DY_prompt->GetBinContent(iB)+DY_fake->GetBinContent(iB));
            float eWZ = errWZ*(WZ_prompt->GetBinContent(iB)+WZ_fake->GetBinContent(iB));
            float eTT = errTT*(TT_prompt->GetBinContent(iB)+TT_fake->GetBinContent(iB));
            float err2 = eVVV*eVVV + eDY*eDY + eTT*eTT+eWZ*eWZ;
            float errTot = TMath::Sqrt(err2);
            dummy->SetBinError(iB,errTot);
        }
    
        TH1D * bkgdErr = (TH1D*)DY_fake->Clone("bkgdErr");
        bkgdErr->SetFillStyle(3013);
        bkgdErr->SetFillColor(1);
        bkgdErr->SetMarkerStyle(21);
        bkgdErr->SetMarkerSize(0);

        for (int iB=1; iB<=nBins; ++iB)
        {
            DY_prompt->SetBinError(iB,0);
            DY_fake->SetBinError(iB,0);
            float eStat =  bkgdErr->GetBinError(iB);
            float X = bkgdErr->GetBinContent(iB);
            float eLumi = errLumi * X;
            float eBkg = dummy->GetBinError(iB);
            float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eBkg*eBkg);
            bkgdErr->SetBinError(iB,Err);
        }
        //Colors
        Int_t colorDY_fake = TColor::GetColor("#ffcc66");
        Int_t colorDY_prompt = TColor::GetColor("#DE5A6A");

        InitData(data_obs);
        InitHist(DY_prompt,"","",colorDY_prompt,1001);
        InitHist(DY_fake,"","",colorDY_fake,1001);

            data_obs->GetXaxis()->SetTitle(xtitle);
            data_obs->GetYaxis()->SetTitle(ytitle);
            data_obs->GetYaxis()->SetTitleOffset(1.5);
            data_obs->GetYaxis()->SetTitleSize(0.06);
            data_obs->GetXaxis()->SetRangeUser(xLower,xUpper);
            float yUpper = data_obs->GetMaximum();
            if (logY)
            data_obs->GetYaxis()->SetRangeUser(0.5,2*yUpper);
            else
            data_obs->GetYaxis()->SetRangeUser(0,2*yUpper);
            data_obs->SetMarkerSize(1.5);
            data_obs->GetXaxis()->SetLabelSize(0);
            data_obs->GetYaxis()->SetLabelSize(0.06);

        TCanvas * canv1 = MakeCanvas("canv1", "", 700, 800);
        
        TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
        upper->Draw();
        upper->cd();
        upper->SetFillColor(0);
        upper->SetBorderMode(0);
        upper->SetBorderSize(10);
        upper->SetTickx(1);
        upper->SetTicky(1);
        upper->SetLeftMargin(0.17);
        upper->SetRightMargin(0.05);
        upper->SetBottomMargin(0.02);
        upper->SetFrameFillStyle(0);
        upper->SetFrameLineStyle(0);
        upper->SetFrameLineWidth(2);
        upper->SetFrameBorderMode(0);
        upper->SetFrameBorderSize(10);
        upper->SetFrameFillStyle(0);
        upper->SetFrameLineStyle(0);
        upper->SetFrameLineWidth(2);
        upper->SetFrameBorderMode(0);
        upper->SetFrameBorderSize(10);

        //Drawing histogram
            data_obs->Draw("e1");
            DY_fake->Draw("sameHIST");
            DY_prompt->Draw("sameHIST");
            data_obs->Draw("e1same");
        
        bkgdErr->Draw("e2same");
    
    
        //Calculating chi2
        float chi2 = 0;
        for (int iB=1; iB<=nBins; ++iB)
        {
            float xData = data_obs->GetBinContent(iB);
            float xMC = DY_fake->GetBinContent(iB);
            if (xMC>1e-1)
            {
                    float diff2 = (xData-xMC)*(xData-xMC);
                    chi2 += diff2/xMC;
            }
        }
        std::cout << std::endl;
        std::cout << "Chi2 = " << chi2 << std::endl;
        std::cout << std::endl;

        float x1Leg = 0.50;
        float x2Leg = 0.70;
        if (legLeft)
        {
                x1Leg = 0.20;
                x2Leg = 0.45;
        }
        TLegend * leg = new TLegend(x1Leg,0.6,x2Leg,0.88);
        SetLegendStyle(leg);
        leg->SetTextSize(0.04);
        leg->AddEntry(data_obs,"Data","lp");
        leg->AddEntry(DY_fake,"MC Fakes","f");
        leg->AddEntry(DY_prompt,"MC prompts","f");
        leg->Draw();
        TLatex * cms = new TLatex(0.20,0.94,"CMS Preliminary L = 59.7 fb^{-1} at #sqrt{s} = 13 TeV, ll"+finalState);

        cms->SetNDC();
        cms->SetTextSize(0.05);
        cms->Draw();
        TLatex * workinprogress = new TLatex(x1Leg,0.55,tags);
        workinprogress->SetNDC();
        workinprogress->SetTextSize(0.05);
        workinprogress->Draw();
        
        if (logY) upper->SetLogy(true);

        upper->Draw("SAME");
        upper->RedrawAxis();
        upper->Modified();
        upper->Update();
        canv1->cd();
        
        TH1D * ratioH = (TH1D*)data_obs->Clone("ratioH");
        TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
        ratioH->SetMarkerColor(1);
        ratioH->SetMarkerStyle(20);
        ratioH->SetMarkerSize(1.5);
        ratioH->SetLineColor(1);
        ratioH->GetYaxis()->SetRangeUser(0.3,1.8);
        ratioH->GetYaxis()->SetNdivisions(505);
        ratioH->GetXaxis()->SetLabelFont(42);
        ratioH->GetXaxis()->SetLabelOffset(0.04);
        ratioH->GetXaxis()->SetLabelSize(0.14);
        ratioH->GetXaxis()->SetTitleSize(0.13);
        ratioH->GetXaxis()->SetTitleOffset(1.2);
        ratioH->GetXaxis()->SetTitle(xtitle);
        ratioH->GetYaxis()->SetTitle("obs/exp");
        ratioH->GetYaxis()->SetLabelFont(42);
        ratioH->GetYaxis()->SetLabelOffset(0.015);
        ratioH->GetYaxis()->SetLabelSize(0.13);
        ratioH->GetYaxis()->SetTitleSize(0.14);
        ratioH->GetYaxis()->SetTitleOffset(0.5);
        ratioH->GetXaxis()->SetTickLength(0.07);
        ratioH->GetYaxis()->SetTickLength(0.04);
        ratioH->GetYaxis()->SetLabelOffset(0.01);
        
        for (int iB=1; iB<=nBins; ++iB)
        {
                float x1 = data_obs->GetBinContent(iB);
                float x2 = DY_fake->GetBinContent(iB);
                ratioErrH->SetBinContent(iB,1.0);
                ratioErrH->SetBinError(iB,0.0);
            float xBkg = bkgdErr->GetBinContent(iB);
            float errBkg = bkgdErr->GetBinError(iB);
            if (xBkg>0)
            {
                    float relErr = errBkg/xBkg;
                    ratioErrH->SetBinError(iB,relErr);
            }
            if (x1>0&&x2>0)
            {
                    float e1 = data_obs->GetBinError(iB);
                    float ratio = x1/x2;
                    float eratio = e1/x2;
                    ratioH->SetBinContent(iB,ratio);
                    ratioH->SetBinError(iB,eratio);
            }
            else
            {
                    ratioH->SetBinContent(iB,1000);
            }
        }

        TPad *lower = new TPad("lower", "pad",0,0,1,0.30);
        lower->Draw();
        lower->cd();
        lower->SetFillColor(0);
        lower->SetBorderMode(0);
        lower->SetBorderSize(10);
        lower->SetGridy();
        lower->SetTickx(1);
        lower->SetTicky(1);
        lower->SetLeftMargin(0.17);
        lower->SetRightMargin(0.05);
        lower->SetTopMargin(0.026);
        lower->SetBottomMargin(0.35);
        lower->SetFrameFillStyle(0);
        lower->SetFrameLineStyle(0);
        lower->SetFrameLineWidth(2);
        lower->SetFrameBorderMode(0);
        lower->SetFrameBorderSize(10);
        lower->SetFrameFillStyle(0);
        lower->SetFrameLineStyle(0);
        lower->SetFrameLineWidth(2);
        lower->SetFrameBorderMode(0);
        lower->SetFrameBorderSize(10);

        ratioH->Draw("e1");
        ratioErrH->Draw("e2same");
        lower->Modified();
        lower->RedrawAxis();
        canv1->cd();
        canv1->Modified();
        canv1->cd();
        canv1->SetSelected(canv1);
        if(getDenominator)
            tags = tags + "_Den";
        else
            tags = tags + "_Num";
        canv1->Print("./FRworkspace/figures/"+varName+"_ll"+finalState+"_"+tags+"_PromptAndFake.pdf","Portrait pdf");
    }

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

void ProduceAZhDataCards_Syst(
                                 TString era = "UL2018",
                                 TString varName = "Mass_svc_escaleUp",
                                 TString UncerSuffix = "CMS_13TeV_2018_eshapeUp",
                                 TString xtitle = "m^{c}_{ll#tau#tau} [GeV]",
                                 TString ytitle = "Events",
                                 bool blindData = true,
                                 TString ZfinalStates = "ee", //options: "ll", "ee", "mm"
                                 TString finalState = "mt",
                                 TString extraCuts = "*(nbtag==0)",
                                 TString CatName = "0btag")
{
    bool checkSS = false;//swtich to check SS region
    bool UseFRShape = false;
	SetStyle();

	double xsec[76] = { //in pico barn (fb)
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
        377.96, //(13)TT to hadronic
        365.35, //(14)TT to semileptonic
        0.2086,//(15) WWW
        0.1651,//(16) WWZ
        0.01398,//(17) ZZZ
        0.05565,//(18) WZZ
        4.708,//(18) WZTo3LNu
        5.595,//(20) WZTo2L2Q
        48.58*0.0627, //(21) GGH tau tau
        3.782*0.0627, //(22) VBF tau tau
        0.5328*0.0627,
        0.8400*0.0627,
        0.884*0.0627,
        0.507*0.0627,
        48.58*0.2137*0.1061,
        3.925*0.2137*0.1061,
        0.114,//HW- H-> WW
        0.180,//HW+ H-> WW
        0.163,//HZJ H -> WW
        0.0262*0.1,//()GGZH H -> WW,Z->LL
        0.0128,//GG H -> ZZ -> 4L
        0.001, //AZh (225)
        0.001, //
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001, //AZh
        0.001, //AZh
        0.001, //AZh
        0.001, //AZh
        0.001,
        0.001,
        0.001,
        0.001, //BBA (225)
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001
    };
    
        float lumi = 59740;
        //Deal with bins and binning
        //int nBins = numberofbins;
        int nBins = 15;
        float nBinsX_AZh[16] = {200,220,240,260,280,300,320,340,360,380,400,450,550,700,1000,2400};

        float bins[100];
        for (int iB=0; iB<=nBins; ++iB)
        {
            bins[iB] = nBinsX_AZh[iB];
        }
        
        int nSamples = 76;
        TH1D * hist[nSamples];   //This is "hist_mumu"
        TH1D * hist_ee[nSamples];
    
        //inistiating cuts
        TString cuts[nSamples];
    
        TString ZeeCuts = "(electronIDWP90v2_1>0.5 && iso_1<0.15 && electronIDWP90v2_2>0.5 && iso_2<0.15)*";
        TString ZmumuCuts = "(muonLoose_1>0.5 && iso_1<0.15 && muonLoose_2>0.5 && iso_2<0.15)*";
        
        TString MCCommonWeight = "IdSF_1*IdSF_2*TrigSF_1*TrigSF_2*puweight*mcweight*prefiringweight*";
        TString SRMCWeight = "";
        TString HiggsSRCuts = "*(H_SS<0.5)";
        extraCuts += "*(passMETFilters==1)";

        if(checkSS == true)
        {
            HiggsSRCuts = "*(H_SS>0.5)";
        }
        TString HiggsSSCuts = "*(H_SS>0.5)";
        TString eeLeg1ID = "electronIDWP90v2_1>0.5 && iso_1<0.15";
        TString eeLeg2ID = "electronIDWP90v2_2>0.5 && iso_2<0.15";
        TString mumuLeg1ID = "muonLoose_1>0.5 && iso_1<0.15";
        TString mumuLeg2ID = "muonLoose_2>0.5 && iso_2<0.15";
        TString Leg3ID = "";
        TString Leg4ID = "";
    
        if(finalState == "em")
        {
            SRMCWeight = "IdSF_3*IdSF_4*";
            Leg3ID = "electronIDWP90v2_3>0.5 && iso_3<0.15";
            Leg4ID = "muonLoose_4>0.5 && iso_4<0.15";
        }
    
        if(finalState == "et")
        {
            SRMCWeight = "IdSF_3*tauVSjetSF_4*tauVSeSF_4*tauVSmuSF_4*";
            Leg3ID = "electronIDWP90v2_3>0.5 && iso_3<0.15";
            Leg4ID = "VSjetMedium_4>0.5 && VSmuVLoose_4>0.5 && VSeTight_4>0.5";
        }
    
        if(finalState == "mt")
        {
            SRMCWeight = "IdSF_3*tauVSjetSF_4*tauVSeSF_4*tauVSmuSF_4*";
            Leg3ID = "muonLoose_3>0.5 && iso_3<0.15";
            Leg4ID = "VSjetMedium_4>0.5 && VSmuTight_4>0.5 && VSeVLoose_4>0.5";
        }
    
        if(finalState == "tt")
        {
            SRMCWeight = "tauVSjetSF_3*tauVSeSF_3*tauVSmuSF_3*tauVSjetSF_4*tauVSeSF_4*tauVSmuSF_4*";
            Leg3ID = "VSjetMedium_3>0.5 && VSmuVLoose_3>0.5 && VSeVLoose_3>0.5";
            Leg4ID = "VSjetMedium_4>0.5 && VSmuVLoose_4>0.5 && VSeVLoose_4>0.5";
        }
    
        for (int i=0; i<nSamples; ++i)
        {
            if(i<34)//for background
                cuts[i] = MCCommonWeight + SRMCWeight +"(" + Leg3ID + " && " + Leg4ID + "&& gen_match_1!=6 && gen_match_2!=6 && gen_match_3!=6 && gen_match_4!=6)";
            else//for signals
                cuts[i] = MCCommonWeight + SRMCWeight +"(" + Leg3ID + " && " + Leg4ID + ")";
        }
        cuts[0] = "("+ Leg3ID + "&&" + Leg4ID + ")";
    

            for (int i=0; i<nSamples; ++i)
            {
                //ee
                std::cout << "EE channel " <<i<<":"+samples_ee_UL2018_DC[i]<<std::endl;
                //TFile * file = new TFile("../NTuples/v7/"+samples[i]+".root");
                TFile * file_ee = new TFile("./NTuples_ee"+finalState+"/"+samples_ee_UL2018_DC[i]+".root");
                TH1D * histWeightsH_ee = (TH1D*)file_ee->Get("histWeightsH");
                TTree * tree_ee = (TTree*)file_ee->Get("SynTree");
                double normaliza_ee = xsec[i]*lumi/histWeightsH_ee->GetSumOfWeights();
                TString histName_ee = samples_ee_UL2018_DC[i] + "_"+varName + "_ee";
                hist_ee[i] = new TH1D(histName_ee,"",nBins,bins);
                hist_ee[i]->Sumw2();
                tree_ee->Draw(varName+">>"+histName_ee,ZeeCuts+cuts[i]+HiggsSRCuts+extraCuts);
                //std::cout << "ee Cuts " <<i<<":"+ ZeeCuts+cuts[i]+HiggsSRCuts+extraCuts <<std::endl;

                if(i > 0)
                {
                    for (int iB=1; iB<=nBins; ++iB)
                    {
                        double x = hist_ee[i]->GetBinContent(iB);
                        double e = hist_ee[i]->GetBinError(iB);
                        hist_ee[i]->SetBinContent(iB,normaliza_ee*x);
                        hist_ee[i]->SetBinError(iB,normaliza_ee*e);
                    }
                }
                std::cout << hist_ee[i]->GetSumOfWeights() << std::endl;
                
                //mumu
                std::cout << "MuMu channel " <<i<<":"+samples_mumu_UL2018_DC[i]<<std::endl;
                TFile * file_mumu = new TFile("./NTuples_mm"+finalState+"/"+samples_mumu_UL2018_DC[i]+".root");
                TH1D * histWeightsH_mumu = (TH1D*)file_mumu->Get("histWeightsH");
                TTree * tree_mumu = (TTree*)file_mumu->Get("SynTree");
                double normaliza_mumu = xsec[i]*lumi/histWeightsH_mumu->GetSumOfWeights();
                TString histName_mumu = samples_mumu_UL2018_DC[i] + "_"+varName+ "_mumu";
                hist[i] = new TH1D(histName_mumu,"",nBins,bins);
                hist[i]->Sumw2();
                tree_mumu->Draw(varName+">>"+histName_mumu,ZmumuCuts+cuts[i]+HiggsSRCuts+extraCuts);
                //std::cout << "MuMu Cuts " <<i<<":" +ZmumuCuts+cuts[i]+HiggsSRCuts+extraCuts <<std::endl;

                if(i > 0)
                {
                    for (int iB=1; iB<=nBins; ++iB)
                    {
                        double x = hist[i]->GetBinContent(iB);
                        double e = hist[i]->GetBinError(iB);
                        hist[i]->SetBinContent(iB,normaliza_mumu*x);
                        hist[i]->SetBinError(iB,normaliza_mumu*e);
                    }
                }
                std::cout << hist[i]->GetSumOfWeights() << std::endl;
                
                //SUM of ee and mumu
                if(ZfinalStates == "ee")
                    hist[i] = hist_ee[i];
                else
                    if(ZfinalStates == "mm")
                        hist[i] = hist[i];
                    else
                        hist[i]->Add(hist[i],hist_ee[i]);
            }
    
    // *******************************
    // ***** DY+Jets samples *******
    // *******************************
    
    TH1D * histZ[9];
    TH1D * histZ_ee[9];

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
    
    
    TString npartonCuts[9] = {"&&(npartons==0||npartons>4)",
        "&&npartons==1",
        "&&npartons==2",
        "&&npartons==3",
        "&&npartons==4",
        "",
        "",
        "",
        ""
    };
    TString cutsZ[9];

    if (finalState == "em")
    {
        for (int iDY=0; iDY<9; ++iDY)
        {
            cutsZ[iDY]   = MCCommonWeight + SRMCWeight + "zptweight*(" + Leg3ID + " && " + Leg4ID + " && gen_match_1!=6 && gen_match_2!=6 && gen_match_3!=6 && gen_match_4!=6"+npartonCuts[iDY]+")";
        }
    }
    if (finalState == "et")
    {
        for (int iDY=0; iDY<9; ++iDY)
        {
            cutsZ[iDY]   = MCCommonWeight + SRMCWeight + "zptweight*(" + Leg3ID + " && " + Leg4ID + " && gen_match_1!=6 && gen_match_2!=6 && gen_match_3!=6 && gen_match_4!=6"+npartonCuts[iDY]+")";
        }
    }
    if (finalState == "mt")
    {
        for (int iDY=0; iDY<9; ++iDY)
        {
            cutsZ[iDY]   = MCCommonWeight + SRMCWeight + "zptweight*(" + Leg3ID + " && " + Leg4ID + " && gen_match_1!=6 && gen_match_2!=6 && gen_match_3!=6 && gen_match_4!=6"+npartonCuts[iDY]+")";
        }
    }
    if (finalState == "tt")
    {
        for (int iDY=0; iDY<9; ++iDY)
        {
            cutsZ[iDY]   = MCCommonWeight + SRMCWeight + "zptweight*(" + Leg3ID + " && " + Leg4ID + " && gen_match_1!=6 && gen_match_2!=6 && gen_match_3!=6 && gen_match_4!=6"+npartonCuts[iDY]+")";
        }
    }
    
    // filling histograms for DYJets samples
    for (int i=0; i<9; ++i) { // run over samples
        TFile * file_ee = new TFile("./NTuples_ee"+finalState+"/"+dySampleNames_UL2018[i]+".root");
        TTree * tree_ee = (TTree*)file_ee->Get("SynTree");
        double normaliza_ee = dyNorm_ee[i];
        TString histNameZ_ee = dySampleNames_UL2018[i] + "_"+varName+"_Z_ee";

        histZ_ee[i] = new TH1D(histNameZ_ee,"",nBins,bins);
        histZ_ee[i]->Sumw2();

        tree_ee->Draw(varName+">>"+histNameZ_ee,ZeeCuts + cutsZ[i] + HiggsSRCuts + extraCuts);

        for (int iB=1; iB<=nBins; ++iB) {
            
            double x = histZ_ee[i]->GetBinContent(iB);
            double e = histZ_ee[i]->GetBinError(iB);
            histZ_ee[i]->SetBinContent(iB,normaliza_ee*x);
            histZ_ee[i]->SetBinError(iB,normaliza_ee*e);
        }
        std::cout << dySampleNames_UL2018[i] << " -> DY Z = " << histZ_ee[i]->GetEntries() << " : " << histZ_ee[i]->GetSumOfWeights() << std::endl;
        
        TFile * file_mumu = new TFile("./NTuples_mm"+finalState+"/"+dySampleNames_UL2018[i]+".root");
        TTree * tree_mumu = (TTree*)file_mumu->Get("SynTree");
        double normaliza_mumu = dyNorm_mumu[i];
        TString histNameZ_mumu = dySampleNames_UL2018[i] + "_"+varName+"_Z_mumu";

        histZ[i] = new TH1D(histNameZ_mumu,"",nBins,bins);
        histZ[i]->Sumw2();

        tree_mumu->Draw(varName+">>"+histNameZ_mumu,ZmumuCuts + cutsZ[i] + HiggsSRCuts + extraCuts);

        for (int iB=1; iB<=nBins; ++iB) {
            double x = histZ[i]->GetBinContent(iB);
            double e = histZ[i]->GetBinError(iB);
            histZ[i]->SetBinContent(iB,normaliza_mumu*x);
            histZ[i]->SetBinError(iB,normaliza_mumu*e);
        }
        std::cout << dySampleNames_UL2018[i] << " -> DY Z = " << histZ[i]->GetEntries() << " : " << histZ[i]->GetSumOfWeights() << std::endl;
        
        //SUM of ee and mumu
        if(ZfinalStates == "ee")
            histZ[i] = histZ_ee[i];
        else
        {
            if(ZfinalStates == "mm")
                histZ[i] = histZ[i];
            else
                histZ[i]->Add(histZ[i],histZ_ee[i]);
        }
    }
    
    hist[1] = histZ[0];
    std::cout << " -> DY Adding = " << hist[1]->GetEntries() << " : " << hist[1]->GetSumOfWeights() << std::endl;
    for (int iDY=1; iDY<9; ++iDY)
    {
        hist[1]->Add(hist[1],histZ[iDY]);
        std::cout << " -> DY Adding = " << hist[1]->GetEntries() << " : " << hist[1]->GetSumOfWeights() << std::endl;

    }
    
    std::cout << "end: DY stitching" << std::endl;
    
    //  adding up ZZ continue backgrounds
        for (int iH=3; iH<8; ++iH)// ggZZ
        {
            hist[2]->Add(hist[2],hist[iH]);
        }
    
        hist[8]->Add(hist[8],hist[9]);
        //ttZ (10)
        //ttW (11)
        // adding up TT backgrounds
        for (int iH=13; iH<15; ++iH)
        {
            hist[12]->Add(hist[12],hist[iH]);
        }
    
        // adding up VVV backgrounds
        for (int iH=16; iH<19; ++iH)
        {
            hist[15]->Add(hist[15],hist[iH]);
        }
    
        // (19+20) is WZ
        hist[19]->Add(hist[19],hist[20]);
        // (23+24) is WH, H->TauTau
        hist[23]->Add(hist[23],hist[24]);
        // (29+30) is WH, H->WW
        hist[29]->Add(hist[29],hist[30]);
    
        //set MC negative bin to 0
        for(int i=0; i<nSamples;i++)
        {
            for (int iB=1; iB<=nBins; ++iB)
            {
                float ySS = hist[i]->GetBinContent(iB);
                if (ySS<0)
                {
                    hist[i]->SetBinContent(iB,0.);
                    hist[i]->SetBinError(iB,0.);
                }
            }
        }


        TFile * file = new TFile("./Combineworkspace/shapes/produceArea/AZh_"+era+"_"+ZfinalStates+finalState+"_"+varName+"_"+CatName+".root","recreate");
    
        file->mkdir(ZfinalStates+finalState);
        file->cd(ZfinalStates+finalState);
        //TH1D * data_obs = (TH1D*)hist[0]->Clone("data_obs_"+UncerSuffix);
        TH1D * DY = (TH1D*)hist[1]->Clone("DY_"+UncerSuffix);
        TH1D * ggZZ = (TH1D*)hist[2]->Clone("ggZZ_"+UncerSuffix);
        TH1D * ZZ = (TH1D*)hist[8]->Clone("ZZ_"+UncerSuffix);
        TH1D * TTZ = (TH1D*)hist[10]->Clone("TTZ_"+UncerSuffix);
        TH1D * TTW = (TH1D*)hist[11]->Clone("TTW_"+UncerSuffix);
        TH1D * TT = (TH1D*)hist[12]->Clone("TT_"+UncerSuffix);
        TH1D * VVV = (TH1D*)hist[15]->Clone("VVV_"+UncerSuffix);
        TH1D * WZ = (TH1D*)hist[19]->Clone("WZ_"+UncerSuffix);
        TH1D * ggHtt = (TH1D*)hist[21]->Clone("ggHtt_"+UncerSuffix);
        TH1D * VBFHtt = (TH1D*)hist[22]->Clone("VBFHtt_"+UncerSuffix);
        TH1D * WHtt = (TH1D*)hist[23]->Clone("WHtt_"+UncerSuffix);
        TH1D * ZHtt = (TH1D*)hist[25]->Clone("ZHtt_"+UncerSuffix);
        TH1D * TTHtt = (TH1D*)hist[26]->Clone("TTHtt_"+UncerSuffix);
        TH1D * ggHWW = (TH1D*)hist[27]->Clone("ggHWW_"+UncerSuffix);
        TH1D * VBFHWW = (TH1D*)hist[28]->Clone("VBFHWW_"+UncerSuffix);
        TH1D * WHWW = (TH1D*)hist[29]->Clone("WHWW_"+UncerSuffix);
        TH1D * ZHWW = (TH1D*)hist[31]->Clone("ZHWW_"+UncerSuffix);
        TH1D * ggZHWW = (TH1D*)hist[32]->Clone("ggZHWW_"+UncerSuffix);
        TH1D * ggHZZ = (TH1D*)hist[33]->Clone("ggHZZ_"+UncerSuffix);

        TH1D * AZh225 = (TH1D*)hist[34]->Clone("AZh225_"+UncerSuffix);
        TH1D * AZh250 = (TH1D*)hist[35]->Clone("AZh250_"+UncerSuffix);
        TH1D * AZh275 = (TH1D*)hist[36]->Clone("AZh275_"+UncerSuffix);
        TH1D * AZh300 = (TH1D*)hist[37]->Clone("AZh300_"+UncerSuffix);
        TH1D * AZh325 = (TH1D*)hist[38]->Clone("AZh325_"+UncerSuffix);
        TH1D * AZh350 = (TH1D*)hist[39]->Clone("AZh350_"+UncerSuffix);
        TH1D * AZh375 = (TH1D*)hist[40]->Clone("AZh375_"+UncerSuffix);
        TH1D * AZh400 = (TH1D*)hist[41]->Clone("AZh400_"+UncerSuffix);
        TH1D * AZh450 = (TH1D*)hist[42]->Clone("AZh450_"+UncerSuffix);
        TH1D * AZh500 = (TH1D*)hist[43]->Clone("AZh500_"+UncerSuffix);
        TH1D * AZh600 = (TH1D*)hist[44]->Clone("AZh600_"+UncerSuffix);
        TH1D * AZh700 = (TH1D*)hist[45]->Clone("AZh700_"+UncerSuffix);
        TH1D * AZh750 = (TH1D*)hist[46]->Clone("AZh750_"+UncerSuffix);
        TH1D * AZh800 = (TH1D*)hist[47]->Clone("AZh800_"+UncerSuffix);
        TH1D * AZh900 = (TH1D*)hist[48]->Clone("AZh900_"+UncerSuffix);
        TH1D * AZh1000 = (TH1D*)hist[49]->Clone("AZh1000_"+UncerSuffix);
        TH1D * AZh1200 = (TH1D*)hist[50]->Clone("AZh1200_"+UncerSuffix);
        TH1D * AZh1400 = (TH1D*)hist[51]->Clone("AZh1400_"+UncerSuffix);
        TH1D * AZh1600 = (TH1D*)hist[52]->Clone("AZh1600_"+UncerSuffix);
        TH1D * AZh1800 = (TH1D*)hist[53]->Clone("AZh1800_"+UncerSuffix);
        TH1D * AZh2000 = (TH1D*)hist[54]->Clone("AZh2000_"+UncerSuffix);

        TH1D * BBAZh225 = (TH1D*)hist[55]->Clone("BBAZh225_"+UncerSuffix);
        TH1D * BBAZh250 = (TH1D*)hist[56]->Clone("BBAZh250_"+UncerSuffix);
        TH1D * BBAZh275 = (TH1D*)hist[57]->Clone("BBAZh275_"+UncerSuffix);
        TH1D * BBAZh300 = (TH1D*)hist[58]->Clone("BBAZh300_"+UncerSuffix);
        TH1D * BBAZh325 = (TH1D*)hist[59]->Clone("BBAZh325_"+UncerSuffix);
        TH1D * BBAZh350 = (TH1D*)hist[60]->Clone("BBAZh350_"+UncerSuffix);
        TH1D * BBAZh375 = (TH1D*)hist[61]->Clone("BBAZh375_"+UncerSuffix);
        TH1D * BBAZh400 = (TH1D*)hist[62]->Clone("BBAZh400_"+UncerSuffix);
        TH1D * BBAZh450 = (TH1D*)hist[63]->Clone("BBAZh450_"+UncerSuffix);
        TH1D * BBAZh500 = (TH1D*)hist[64]->Clone("BBAZh500_"+UncerSuffix);
        TH1D * BBAZh600 = (TH1D*)hist[65]->Clone("BBAZh600_"+UncerSuffix);
        TH1D * BBAZh700 = (TH1D*)hist[66]->Clone("BBAZh700_"+UncerSuffix);
        TH1D * BBAZh750 = (TH1D*)hist[67]->Clone("BBAZh750_"+UncerSuffix);
        TH1D * BBAZh800 = (TH1D*)hist[68]->Clone("BBAZh800_"+UncerSuffix);
        TH1D * BBAZh900 = (TH1D*)hist[69]->Clone("BBAZh900_"+UncerSuffix);
        TH1D * BBAZh1000 = (TH1D*)hist[70]->Clone("BBAZh1000_"+UncerSuffix);
        TH1D * BBAZh1200 = (TH1D*)hist[71]->Clone("BBAZh1200_"+UncerSuffix);
        TH1D * BBAZh1400 = (TH1D*)hist[72]->Clone("BBAZh1400_"+UncerSuffix);
        TH1D * BBAZh1600 = (TH1D*)hist[73]->Clone("BBAZh1600_"+UncerSuffix);
        TH1D * BBAZh1800 = (TH1D*)hist[74]->Clone("BBAZh1800_"+UncerSuffix);
        TH1D * BBAZh2000 = (TH1D*)hist[75]->Clone("BBAZh2000_"+UncerSuffix);

        float totData = 0;
        float totDY = DY->GetSumOfWeights();
        float totggZZ = ggZZ->GetSumOfWeights();
        float totZZ = ZZ->GetSumOfWeights();
        float totTTZ = TTZ->GetSumOfWeights();
        float totTTW = TTW->GetSumOfWeights();
        float totTT = TT->GetSumOfWeights();
        float totVVV = VVV->GetSumOfWeights();
        float totWZ = WZ->GetSumOfWeights();
        float totggHtt = ggHtt->GetSumOfWeights();
        float totVBFHtt = VBFHtt->GetSumOfWeights();
        float totWHtt = WHtt->GetSumOfWeights();
        float totZHtt = ZHtt->GetSumOfWeights();
        float totTTHtt = TTHtt->GetSumOfWeights();
        float totggHWW = ggHWW->GetSumOfWeights();
        float totVBFHWW = VBFHWW->GetSumOfWeights();
        float totWHWW = WHWW->GetSumOfWeights();
        float totZHWW = ZHWW->GetSumOfWeights();
        float totggZHWW = ggZHWW->GetSumOfWeights();
        float totggHZZ = ggHZZ->GetSumOfWeights();
        float totSMH = totggHtt + totVBFHtt + totWHtt + totZHtt + totTTHtt + totggHWW + totVBFHWW + totWHWW + totZHWW + totggZHWW + totggHZZ;
    
        float totBkg = totDY+totggZZ+totZZ+totTTW+totTTZ+totTT+totVVV+totWZ+totSMH;

        if(!blindData)
            std::cout << "data : " << totData << std::endl;
        std::cout << "Total bkg : " << totBkg << std::endl;
        std::cout << "GGAZh(225) : " << AZh225->GetSumOfWeights() << std::endl;
        std::cout << "GGAZh(300) : " << AZh300->GetSumOfWeights() << ", " << AZh300->GetEntries() << std::endl;
        std::cout << "GGAZh(400) : " << AZh400->GetSumOfWeights() << std::endl;
        std::cout << "GGAZh(700) : " << AZh700->GetSumOfWeights() << std::endl;
        std::cout << "GGAZh(1000) : " << AZh1000->GetSumOfWeights() << ", " << AZh1000->GetEntries() << std::endl;
        std::cout << "BBAZh(225) : " << BBAZh225->GetSumOfWeights() << std::endl;
        std::cout << "BBAZh(300) : " << BBAZh300->GetSumOfWeights() << ", " << BBAZh300->GetEntries() << std::endl;
        std::cout << "BBAZh(400) : " << BBAZh400->GetSumOfWeights() << std::endl;
        std::cout << "BBAZh(700) : " << BBAZh700->GetSumOfWeights() << std::endl;
        std::cout << "BBAZh(1000) : " << BBAZh1000->GetSumOfWeights() << ", "  << AZh1000->GetEntries()<< std::endl;
        std::cout << "DY : " << DY->GetSumOfWeights() << std::endl;
        std::cout << "ggZZ  : " << ggZZ->GetSumOfWeights() << std::endl;
        std::cout << "ZZ  : " << ZZ->GetSumOfWeights() << std::endl;
        std::cout << "TTZ  : " << TTZ->GetSumOfWeights() << std::endl;
        std::cout << "TTW  : " << TTW->GetSumOfWeights() << std::endl;
        std::cout << "TT  : " << TT->GetSumOfWeights() << std::endl;
        std::cout << "VVV : " << VVV->GetSumOfWeights() << std::endl;
        std::cout << "WZ : " << WZ->GetSumOfWeights() << std::endl;
        std::cout << "SM H : " << totSMH << std::endl;

        file->Write();
        file->Close();
    }

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

void ProduceAZhDataCards(
                                 TString era = "UL2016preVFP",
                                 TString varName = "Mass_svc",
                                 TString xtitle = "m^{c}_{ll#tau#tau} [GeV]",
                                 TString ytitle = "Events",
                                 bool blindData = true,
                                 TString ZfinalStates = "ee", //options: "ll", "ee", "mm"
                                 TString finalState = "mt",
                                 TString extraCuts = "",
                                 TString CatName = "inclusive")
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
        0.163,//HZJ H -> WW ->2L2Nu, Z->LL
        0.0262*0.1,//()GGZH H -> WW, Z->LL
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
    
        float lumi = 19500;
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
    
        //Data-driven Normalization Hist
        TH1D * hist_DR[nSamples];
        TH1D * hist_DR_ee[nSamples];

        //Data-driven Shape Hist -- Same Sign region
        TH1D * hist_SS;
        TH1D * hist_SS_ee;

        //inistiating cuts
        TString cuts[nSamples];
    
        TString ZeeCuts = "(electronIDWP90v2_1>0.5 && iso_1<0.15 && electronIDWP90v2_2>0.5 && iso_2<0.15)*";
        TString ZmumuCuts = "(muonLoose_1>0.5 && iso_1<0.15 && muonLoose_2>0.5 && iso_2<0.15)*";
        
        TString MCCommonWeight = "IdSF_1*IdSF_2*TrigSF_1*TrigSF_2*puweight*mcweight*prefiringweight*";
        TString SRMCWeight = "";
        extraCuts += "*(passMETFilters)";

        TString Cat0Cuts = "";
        TString Cat1Cuts = "";
        TString Cat2Cuts = "";

        TString eeCat1Cuts = "";
        TString eeCat2Cuts = "";
        TString eeCat3Cuts = "";
        TString eeCat4Cuts = "";
        TString eeCat12Cuts = "";
        TString eeCat13Cuts = "";
        TString eeCat14Cuts = "";
        TString eeCat23Cuts = "";
        TString eeCat24Cuts = "";
        TString eeCat34Cuts = "";
        TString eeCat123Cuts = "";
        TString eeCat124Cuts = "";
        TString eeCat234Cuts = "";
        TString eeCat1234Cuts = "";
        TString mumuCat1Cuts = "";
        TString mumuCat2Cuts = "";
        TString mumuCat3Cuts = "";
        TString mumuCat4Cuts = "";
        TString mumuCat12Cuts = "";
        TString mumuCat13Cuts = "";
        TString mumuCat14Cuts = "";
        TString mumuCat23Cuts = "";
        TString mumuCat24Cuts = "";
        TString mumuCat34Cuts = "";
        TString mumuCat123Cuts = "";
        TString mumuCat124Cuts = "";
        TString mumuCat234Cuts = "";
        TString mumuCat1234Cuts = "";
        TString cutsRelaxed = "";
        TString HiggsSRCuts = "*(H_SS<0.5)";
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
            cutsRelaxed = "(iso_3<9999 && iso_4<9999)";
        }
    
        if(finalState == "et")
        {
            SRMCWeight = "IdSF_3*tauVSjetSF_4*tauVSeSF_4*tauVSmuSF_4*";
            Leg3ID = "electronIDWP90v2_3>0.5 && iso_3<0.15";
            Leg4ID = "VSjetMedium_4>0.5 && VSmuVLoose_4>0.5 && VSeTight_4>0.5";
            cutsRelaxed = "(VSjetVVVLoose_4>0.5 && VSmuVLoose_4>0.5 && VSeVLoose_4>0.5)";
        }
    
        if(finalState == "mt")
        {
            SRMCWeight = "IdSF_3*tauVSjetSF_4*tauVSeSF_4*tauVSmuSF_4*";
            Leg3ID = "muonLoose_3>0.5 && iso_3<0.15";
            Leg4ID = "VSjetMedium_4>0.5 && VSmuTight_4>0.5 && VSeVLoose_4>0.5";
            cutsRelaxed = "(VSjetVVVLoose_4>0.5 && VSmuVLoose_4>0.5 && VSeVLoose_4>0.5)";
        }
    
        if(finalState == "tt")
        {
            SRMCWeight = "tauVSjetSF_3*tauVSeSF_3*tauVSmuSF_3*tauVSjetSF_4*tauVSeSF_4*tauVSmuSF_4*";
            Leg3ID = "VSjetMedium_3>0.5 && VSmuVLoose_3>0.5 && VSeVLoose_3>0.5";
            Leg4ID = "VSjetMedium_4>0.5 && VSmuVLoose_4>0.5 && VSeVLoose_4>0.5";
            cutsRelaxed = "(VSjetVVVLoose_3>0.5 && VSmuVLoose_3>0.5 && VSeVLoose_3 >0.5 && VSjetVVVLoose_4>0.5 && VSmuVLoose_4>0.5 && VSeVLoose_4>0.5)";
        }
    
        for (int i=0; i<nSamples; ++i)
        {
            if(i<34)
                cuts[i] = MCCommonWeight + SRMCWeight +"(" + Leg3ID + " && " + Leg4ID + "&& gen_match_1!=6 && gen_match_2!=6 && gen_match_3!=6 && gen_match_4!=6)";
            else
                cuts[i] = MCCommonWeight + SRMCWeight +"(" + Leg3ID + " && " + Leg4ID + ")";
        }
        cuts[0] = "("+ Leg3ID + "&&" + Leg4ID + ")";
    
        eeCat1Cuts = "fakeRateWeight_Z1*(!(" + eeLeg1ID + ") && ("+ eeLeg2ID + ") && (" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat2Cuts = "fakeRateWeight_Z2*((" + eeLeg1ID + ") && !("+ eeLeg2ID + ") && (" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat3Cuts = "fakeRateWeight_1*((" + eeLeg1ID + ") && (" + eeLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat4Cuts = "fakeRateWeight_2*((" + eeLeg1ID + ") && (" + eeLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        eeCat12Cuts = "fakeRateWeight_Z1*fakeRateWeight_Z2*(!(" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && (" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat13Cuts = "fakeRateWeight_Z1*fakeRateWeight_1*(!(" + eeLeg1ID + ") && (" + eeLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat14Cuts = "fakeRateWeight_Z1*fakeRateWeight_2*(!(" + eeLeg1ID + ") && (" + eeLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        eeCat23Cuts = "fakeRateWeight_Z2*fakeRateWeight_1*((" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat24Cuts = "fakeRateWeight_Z2*fakeRateWeight_2*((" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        eeCat34Cuts = "fakeRateWeight_1*fakeRateWeight_2*((" + eeLeg1ID + ") && (" + eeLeg2ID + ") && !(" + Leg3ID + ") && !(" + Leg4ID + "))";
        eeCat123Cuts = "fakeRateWeight_Z1*fakeRateWeight_Z2*fakeRateWeight_1*(!(" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat124Cuts = "fakeRateWeight_Z1*fakeRateWeight_Z2*fakeRateWeight_2*(!(" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        eeCat234Cuts = "fakeRateWeight_Z2*fakeRateWeight_1*fakeRateWeight_2*((" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && !(" + Leg3ID + ") && !(" + Leg4ID + "))";
        eeCat1234Cuts = "fakeRateWeight_Z1*fakeRateWeight_Z2*fakeRateWeight_1*fakeRateWeight_2*(!(" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && !(" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat1Cuts = "fakeRateWeight_Z1*(!(" + mumuLeg1ID + ") && ("+ mumuLeg2ID + ") && (" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat2Cuts = "fakeRateWeight_Z2*((" + mumuLeg1ID + ") && !("+ mumuLeg2ID + ") && (" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat3Cuts = "fakeRateWeight_1*((" + mumuLeg1ID + ") && (" + mumuLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat4Cuts = "fakeRateWeight_2*((" + mumuLeg1ID + ") && (" + mumuLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat12Cuts = "fakeRateWeight_Z1*fakeRateWeight_Z2*(!(" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && (" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat13Cuts = "fakeRateWeight_Z1*fakeRateWeight_1*(!(" + mumuLeg1ID + ") && (" + mumuLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat14Cuts = "fakeRateWeight_Z1*fakeRateWeight_2*(!(" + mumuLeg1ID + ") && (" + mumuLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat23Cuts = "fakeRateWeight_Z2*fakeRateWeight_1*((" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat24Cuts = "fakeRateWeight_Z2*fakeRateWeight_2*((" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat34Cuts = "fakeRateWeight_1*fakeRateWeight_2*((" + mumuLeg1ID + ") && (" + mumuLeg2ID + ") && !(" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat123Cuts = "fakeRateWeight_Z1*fakeRateWeight_Z2*fakeRateWeight_1*(!(" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat124Cuts = "fakeRateWeight_Z1*fakeRateWeight_Z2*fakeRateWeight_2*(!(" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat234Cuts = "fakeRateWeight_Z2*fakeRateWeight_1*fakeRateWeight_2*((" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && !(" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat1234Cuts = "fakeRateWeight_Z1*fakeRateWeight_Z2*fakeRateWeight_1*fakeRateWeight_2*(!(" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && !(" + Leg3ID + ") && !(" + Leg4ID + "))";

            for (int i=0; i<nSamples; ++i)
            {
                //ee
                std::cout << "EE channel " <<i<<":"+samples_ee_UL2016preVFP_DC[i]<<std::endl;
                //TFile * file = new TFile("../NTuples/v7/"+samples[i]+".root");
                TFile * file_ee = new TFile("./NTuples_ee"+finalState+"/"+samples_ee_UL2016preVFP_DC[i]+".root");
                TH1D * histWeightsH_ee = (TH1D*)file_ee->Get("histWeightsH");
                TTree * tree_ee = (TTree*)file_ee->Get("SynTree");
                double normaliza_ee = xsec[i]*lumi/histWeightsH_ee->GetSumOfWeights();
                TString histName_ee = samples_ee_UL2016preVFP_DC[i] + "_"+varName + "_ee";
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
                std::cout << "MuMu channel " <<i<<":"+samples_mumu_UL2016preVFP_DC[i]<<std::endl;
                TFile * file_mumu = new TFile("./NTuples_mm"+finalState+"/"+samples_mumu_UL2016preVFP_DC[i]+".root");
                TH1D * histWeightsH_mumu = (TH1D*)file_mumu->Get("histWeightsH");
                TTree * tree_mumu = (TTree*)file_mumu->Get("SynTree");
                double normaliza_mumu = xsec[i]*lumi/histWeightsH_mumu->GetSumOfWeights();
                TString histName_mumu = samples_mumu_UL2016preVFP_DC[i] + "_"+varName+ "_mumu";
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
        TFile * file_ee = new TFile("./NTuples_ee"+finalState+"/"+dyRefSamples_UL2016preVFP[iDY]+".root");
        TH1D * histWeightsH_ee = (TH1D*)file_ee->Get("histWeightsH");
        dyRefEvents_ee[iDY] = histWeightsH_ee->GetSumOfWeights();
        
        TFile * file_mumu = new TFile("./NTuples_mm"+finalState+"/"+dyRefSamples_UL2016preVFP[iDY]+".root");
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
        TFile * file_ee = new TFile("./NTuples_ee"+finalState+"/"+dySampleNames_UL2016preVFP[i]+".root");
        TTree * tree_ee = (TTree*)file_ee->Get("SynTree");
        double normaliza_ee = dyNorm_ee[i];
        TString histNameZ_ee = dySampleNames_UL2016preVFP[i] + "_"+varName+"_Z_ee";

        histZ_ee[i] = new TH1D(histNameZ_ee,"",nBins,bins);
        histZ_ee[i]->Sumw2();

        tree_ee->Draw(varName+">>"+histNameZ_ee,ZeeCuts + cutsZ[i] + HiggsSRCuts + extraCuts);

        for (int iB=1; iB<=nBins; ++iB) {
            
            double x = histZ_ee[i]->GetBinContent(iB);
            double e = histZ_ee[i]->GetBinError(iB);
            histZ_ee[i]->SetBinContent(iB,normaliza_ee*x);
            histZ_ee[i]->SetBinError(iB,normaliza_ee*e);
        }
        std::cout << dySampleNames_UL2016preVFP[i] << " -> DY Z = " << histZ_ee[i]->GetEntries() << " : " << histZ_ee[i]->GetSumOfWeights() << std::endl;
        
        TFile * file_mumu = new TFile("./NTuples_mm"+finalState+"/"+dySampleNames_UL2016preVFP[i]+".root");
        TTree * tree_mumu = (TTree*)file_mumu->Get("SynTree");
        double normaliza_mumu = dyNorm_mumu[i];
        TString histNameZ_mumu = dySampleNames_UL2016preVFP[i] + "_"+varName+"_Z_mumu";

        histZ[i] = new TH1D(histNameZ_mumu,"",nBins,bins);
        histZ[i]->Sumw2();

        tree_mumu->Draw(varName+">>"+histNameZ_mumu,ZmumuCuts + cutsZ[i] + HiggsSRCuts + extraCuts);

        for (int iB=1; iB<=nBins; ++iB) {
            double x = histZ[i]->GetBinContent(iB);
            double e = histZ[i]->GetBinError(iB);
            histZ[i]->SetBinContent(iB,normaliza_mumu*x);
            histZ[i]->SetBinError(iB,normaliza_mumu*e);
        }
        std::cout << dySampleNames_UL2016preVFP[i] << " -> DY Z = " << histZ[i]->GetEntries() << " : " << histZ[i]->GetSumOfWeights() << std::endl;
        
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
    //Start Data driven method
    
    TFile * file_data = new TFile("./NTuples_mm"+finalState+"/"+samples_mumu_UL2016preVFP_DC[0]+".root");
    TTree * tree_data = (TTree*)file_data->Get("SynTree");
    
    TString histName_dataCat1 = samples_mumu_UL2016preVFP_DC[0] + "_Cat1_"+varName;
    TString histName_dataCat2 = samples_mumu_UL2016preVFP_DC[0] + "_Cat2_"+varName;
    TString histName_dataCat3 = samples_mumu_UL2016preVFP_DC[0] + "_Cat3_"+varName;
    TString histName_dataCat4 = samples_mumu_UL2016preVFP_DC[0] + "_Cat4_"+varName;
    TString histName_dataCat12 = samples_mumu_UL2016preVFP_DC[0] + "_Cat12_"+varName;
    TString histName_dataCat13 = samples_mumu_UL2016preVFP_DC[0] + "_Cat13_"+varName;
    TString histName_dataCat14 = samples_mumu_UL2016preVFP_DC[0] + "_Cat14_"+varName;
    TString histName_dataCat23 = samples_mumu_UL2016preVFP_DC[0] + "_Cat23_"+varName;
    TString histName_dataCat24 = samples_mumu_UL2016preVFP_DC[0] + "_Cat24_"+varName;
    TString histName_dataCat34 = samples_mumu_UL2016preVFP_DC[0] + "_Cat34_"+varName;
    TString histName_dataCat123 = samples_mumu_UL2016preVFP_DC[0] + "_Cat123_"+varName;
    TString histName_dataCat124 = samples_mumu_UL2016preVFP_DC[0] + "_Cat124_"+varName;
    TString histName_dataCat234 = samples_mumu_UL2016preVFP_DC[0] + "_Cat234_"+varName;
    TString histName_dataCat1234 = samples_mumu_UL2016preVFP_DC[0] + "_Cat1234_"+varName;
    
    hist_DR[0] = new TH1D(histName_dataCat1,"",nBins,bins);
    hist_DR[1] = new TH1D(histName_dataCat2,"",nBins,bins);
    hist_DR[2] = new TH1D(histName_dataCat3,"",nBins,bins);
    hist_DR[3] = new TH1D(histName_dataCat4,"",nBins,bins);
    hist_DR[4] = new TH1D(histName_dataCat12,"",nBins,bins);
    hist_DR[5] = new TH1D(histName_dataCat13,"",nBins,bins);
    hist_DR[6] = new TH1D(histName_dataCat14,"",nBins,bins);
    hist_DR[7] = new TH1D(histName_dataCat23,"",nBins,bins);
    hist_DR[8] = new TH1D(histName_dataCat24,"",nBins,bins);
    hist_DR[9] = new TH1D(histName_dataCat34,"",nBins,bins);
    hist_DR[10] = new TH1D(histName_dataCat123,"",nBins,bins);
    hist_DR[11] = new TH1D(histName_dataCat124,"",nBins,bins);
    hist_DR[12] = new TH1D(histName_dataCat234,"",nBins,bins);
    hist_DR[13] = new TH1D(histName_dataCat1234,"",nBins,bins);
    hist_DR[0]->Sumw2();
    hist_DR[1]->Sumw2();
    hist_DR[2]->Sumw2();
    hist_DR[3]->Sumw2();
    hist_DR[4]->Sumw2();
    hist_DR[5]->Sumw2();
    hist_DR[6]->Sumw2();
    hist_DR[7]->Sumw2();
    hist_DR[8]->Sumw2();
    hist_DR[9]->Sumw2();
    hist_DR[10]->Sumw2();
    hist_DR[11]->Sumw2();
    hist_DR[12]->Sumw2();
    hist_DR[13]->Sumw2();
    tree_data->Draw(varName+">>"+histName_dataCat1,mumuCat1Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat2,mumuCat2Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat3,mumuCat3Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat4,mumuCat4Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat12,mumuCat12Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat13,mumuCat13Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat14,mumuCat14Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat23,mumuCat23Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat24,mumuCat24Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat34,mumuCat34Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat123,mumuCat123Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat124,mumuCat124Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat234,mumuCat234Cuts + HiggsSRCuts + extraCuts);
    tree_data->Draw(varName+">>"+histName_dataCat1234,mumuCat1234Cuts + HiggsSRCuts + extraCuts);
    
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat1 = " << hist_DR[0]->GetEntries() << " : " << hist_DR[0]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat2 = " << hist_DR[1]->GetEntries() << " : " << hist_DR[1]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat3 = " << hist_DR[2]->GetEntries() << " : " << hist_DR[2]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat4 = " << hist_DR[3]->GetEntries() << " : " << hist_DR[3]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat12 = " << hist_DR[4]->GetEntries() << " : " << hist_DR[4]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat13 = " << hist_DR[5]->GetEntries() << " : " << hist_DR[5]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat14 = " << hist_DR[6]->GetEntries() << " : " << hist_DR[6]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat23 = " << hist_DR[7]->GetEntries() << " : " << hist_DR[7]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat24 = " << hist_DR[8]->GetEntries() << " : " << hist_DR[8]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat34 = " << hist_DR[9]->GetEntries() << " : " << hist_DR[9]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat123 = " << hist_DR[10]->GetEntries() << " : " << hist_DR[10]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat124 = " << hist_DR[11]->GetEntries() << " : " << hist_DR[11]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat234 = " << hist_DR[12]->GetEntries() << " : " << hist_DR[12]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven Cat1234 = " << hist_DR[13]->GetEntries() << " : " << hist_DR[13]->GetSumOfWeights() << std::endl;

    float reducibleYield_mumu = hist_DR[0]->GetSumOfWeights() + hist_DR[1]->GetSumOfWeights() + hist_DR[2]->GetSumOfWeights() + hist_DR[3]->GetSumOfWeights() - hist_DR[4]->GetSumOfWeights() - hist_DR[5]->GetSumOfWeights() - hist_DR[6]->GetSumOfWeights() - hist_DR[7]->GetSumOfWeights() - hist_DR[8]->GetSumOfWeights() - hist_DR[9]->GetSumOfWeights() + hist_DR[10]->GetSumOfWeights() + hist_DR[11]->GetSumOfWeights() + hist_DR[12]->GetSumOfWeights() - hist_DR[13]->GetSumOfWeights();
    
    float reducibleYield_mumu_old = hist_DR[2]->GetSumOfWeights() + hist_DR[3]->GetSumOfWeights() - hist_DR[9]->GetSumOfWeights();
    hist_DR[13]->Add(hist_DR[13],hist_DR[12],-1,1);

    std::cout << " Reducible Yields(mumu) = " << reducibleYield_mumu << std::endl;

    TFile * file_data_ee = new TFile("./NTuples_ee"+finalState+"/"+samples_ee_UL2016preVFP_DC[0]+".root");
    TTree * tree_data_ee = (TTree*)file_data_ee->Get("SynTree");
    
    TString histName_dataCat1_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat1_"+varName;
    TString histName_dataCat2_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat2_"+varName;
    TString histName_dataCat3_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat3_"+varName;
    TString histName_dataCat4_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat4_"+varName;
    TString histName_dataCat12_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat12_"+varName;
    TString histName_dataCat13_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat13_"+varName;
    TString histName_dataCat14_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat14_"+varName;
    TString histName_dataCat23_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat23_"+varName;
    TString histName_dataCat24_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat24_"+varName;
    TString histName_dataCat34_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat34_"+varName;
    TString histName_dataCat123_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat123_"+varName;
    TString histName_dataCat124_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat124_"+varName;
    TString histName_dataCat234_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat234_"+varName;
    TString histName_dataCat1234_ee = samples_ee_UL2016preVFP_DC[0] + "_Cat1234_"+varName;
    
    hist_DR_ee[0] = new TH1D(histName_dataCat1_ee,"",nBins,bins);
    hist_DR_ee[1] = new TH1D(histName_dataCat2_ee,"",nBins,bins);
    hist_DR_ee[2] = new TH1D(histName_dataCat3_ee,"",nBins,bins);
    hist_DR_ee[3] = new TH1D(histName_dataCat4_ee,"",nBins,bins);
    hist_DR_ee[4] = new TH1D(histName_dataCat12_ee,"",nBins,bins);
    hist_DR_ee[5] = new TH1D(histName_dataCat13_ee,"",nBins,bins);
    hist_DR_ee[6] = new TH1D(histName_dataCat14_ee,"",nBins,bins);
    hist_DR_ee[7] = new TH1D(histName_dataCat23_ee,"",nBins,bins);
    hist_DR_ee[8] = new TH1D(histName_dataCat24_ee,"",nBins,bins);
    hist_DR_ee[9] = new TH1D(histName_dataCat34_ee,"",nBins,bins);
    hist_DR_ee[10] = new TH1D(histName_dataCat123_ee,"",nBins,bins);
    hist_DR_ee[11] = new TH1D(histName_dataCat124_ee,"",nBins,bins);
    hist_DR_ee[12] = new TH1D(histName_dataCat234_ee,"",nBins,bins);
    hist_DR_ee[13] = new TH1D(histName_dataCat1234_ee,"",nBins,bins);
    hist_DR_ee[0]->Sumw2();
    hist_DR_ee[1]->Sumw2();
    hist_DR_ee[2]->Sumw2();
    hist_DR_ee[3]->Sumw2();
    hist_DR_ee[4]->Sumw2();
    hist_DR_ee[5]->Sumw2();
    hist_DR_ee[6]->Sumw2();
    hist_DR_ee[7]->Sumw2();
    hist_DR_ee[8]->Sumw2();
    hist_DR_ee[9]->Sumw2();
    hist_DR_ee[10]->Sumw2();
    hist_DR_ee[11]->Sumw2();
    hist_DR_ee[12]->Sumw2();
    hist_DR_ee[13]->Sumw2();

    tree_data_ee->Draw(varName+">>"+histName_dataCat1_ee,eeCat1Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat2_ee,eeCat2Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat3_ee,eeCat3Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat4_ee,eeCat4Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat12_ee,eeCat12Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat13_ee,eeCat13Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat14_ee,eeCat14Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat23_ee,eeCat23Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat24_ee,eeCat24Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat34_ee,eeCat34Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat123_ee,eeCat123Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat124_ee,eeCat124Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat234_ee,eeCat234Cuts + HiggsSRCuts + extraCuts);
    tree_data_ee->Draw(varName+">>"+histName_dataCat1234_ee,eeCat1234Cuts + HiggsSRCuts + extraCuts);
    
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat1 = " << hist_DR_ee[0]->GetEntries() << " : " << hist_DR_ee[0]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat2 = " << hist_DR_ee[1]->GetEntries() << " : " << hist_DR_ee[1]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat3 = " << hist_DR_ee[2]->GetEntries() << " : " << hist_DR_ee[2]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat4 = " << hist_DR_ee[3]->GetEntries() << " : " << hist_DR_ee[3]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat12 = " << hist_DR_ee[4]->GetEntries() << " : " << hist_DR_ee[4]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat13 = " << hist_DR_ee[5]->GetEntries() << " : " << hist_DR_ee[5]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat14 = " << hist_DR_ee[6]->GetEntries() << " : " << hist_DR_ee[6]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat23 = " << hist_DR_ee[7]->GetEntries() << " : " << hist_DR_ee[7]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat24 = " << hist_DR_ee[8]->GetEntries() << " : " << hist_DR_ee[8]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat34 = " << hist_DR_ee[9]->GetEntries() << " : " << hist_DR_ee[9]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat123 = " << hist_DR_ee[10]->GetEntries() << " : " << hist_DR_ee[10]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat124 = " << hist_DR_ee[11]->GetEntries() << " : " << hist_DR_ee[11]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat234 = " << hist_DR_ee[12]->GetEntries() << " : " << hist_DR_ee[12]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven Cat1234 = " << hist_DR_ee[13]->GetEntries() << " : " << hist_DR_ee[13]->GetSumOfWeights() << std::endl;

    float reducibleYield_ee = hist_DR_ee[0]->GetSumOfWeights() + hist_DR_ee[1]->GetSumOfWeights() + hist_DR_ee[2]->GetSumOfWeights() + hist_DR_ee[3]->GetSumOfWeights() - hist_DR_ee[4]->GetSumOfWeights() - hist_DR_ee[5]->GetSumOfWeights() - hist_DR_ee[6]->GetSumOfWeights() - hist_DR_ee[7]->GetSumOfWeights() - hist_DR_ee[8]->GetSumOfWeights() - hist_DR_ee[9]->GetSumOfWeights() + hist_DR_ee[10]->GetSumOfWeights() + hist_DR_ee[11]->GetSumOfWeights() + hist_DR_ee[12]->GetSumOfWeights() - hist_DR_ee[13]->GetSumOfWeights();
    std::cout << " Reducible Yields(ee) = " << reducibleYield_ee << std::endl;
    
    float reducibleYield_ee_old = hist_DR_ee[2]->GetSumOfWeights() + hist_DR_ee[3]->GetSumOfWeights() - hist_DR_ee[9]->GetSumOfWeights();
    
    float reducibleYield;

    //SUM of ee and mumu
    if(ZfinalStates == "ee")
        reducibleYield = reducibleYield_ee;
    else
    {
        if(ZfinalStates == "mm")
            reducibleYield = reducibleYield_mumu;
        else
            reducibleYield = reducibleYield_mumu + reducibleYield_ee;
    }
        
    std::cout << " Reducible Yields(total) = " << reducibleYield << std::endl;
    
    std::cout << " Reducible Yields(total, old) = " << reducibleYield_mumu_old + reducibleYield_ee_old << std::endl;

    
    //Data driven shapes: option 1 from SS relaxed region
    TString histName_data_SS = samples_mumu_UL2016preVFP_DC[0] + "_SS_"+varName;
    hist_SS = new TH1D(histName_data_SS,"",nBins,bins);
    hist_SS->Sumw2();
    tree_data->Draw(varName+">>"+histName_data_SS, cutsRelaxed + HiggsSSCuts + extraCuts);
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven SS mumu (before Nor.) = " << hist_SS->GetEntries() << " : " << hist_SS->GetSumOfWeights() << std::endl;
    hist_SS->Scale(reducibleYield_mumu/hist_SS->GetSumOfWeights());

    TString histName_data_SS_ee = samples_ee_UL2016preVFP_DC[0] + "_SS_"+varName;
    hist_SS_ee = new TH1D(histName_data_SS_ee,"",nBins,bins);
    hist_SS_ee->Sumw2();
    tree_data_ee->Draw(varName+">>"+histName_data_SS_ee, cutsRelaxed + HiggsSSCuts + extraCuts);
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven SS ee (before Nor.) = " << hist_SS_ee->GetEntries() << " : " << hist_SS_ee->GetSumOfWeights() << std::endl;
    hist_SS_ee->Scale(reducibleYield_ee/hist_SS_ee->GetSumOfWeights());
       
    if(ZfinalStates == "ee")
        hist_SS = hist_SS_ee;
    else
    {
        if(ZfinalStates == "mm")
            hist_SS = hist_SS;
        else
            hist_SS->Add(hist_SS,hist_SS_ee);
    }
    
    //Data Driven shapes: option 2 from FR region
    TH1D * hist_FR;
    TH1D * hist_FR_ee;

    TString cuts_FR_ee = "!("+ eeLeg1ID + "&&"+ eeLeg2ID+ "&&" + Leg3ID + "&&" + Leg4ID + ")";
    TString cuts_FR_mumu = "!("+ mumuLeg1ID + "&&"+ mumuLeg2ID+ "&&" + Leg3ID + "&&" + Leg4ID + ")";
    
    TString histName_data_FR = samples_mumu_UL2016preVFP_DC[0] + "_FR_"+varName;
    hist_FR = new TH1D(histName_data_FR,"",nBins,bins);
    hist_FR->Sumw2();
    tree_data->Draw(varName+">>"+histName_data_FR, cuts_FR_mumu + HiggsSRCuts);
    std::cout << samples_mumu_UL2016preVFP_DC[0] << " -> Data Driven FR mumu (before Nor.) = " << hist_FR->GetEntries() << " : " << hist_FR->GetSumOfWeights() << std::endl;
    hist_FR->Scale(reducibleYield_mumu/hist_FR->GetSumOfWeights());

    TString histName_data_FR_ee = samples_ee_UL2016preVFP_DC[0] + "_FR_"+varName;
    hist_FR_ee = new TH1D(histName_data_FR_ee,"",nBins,bins);
    hist_FR_ee->Sumw2();
    tree_data_ee->Draw(varName+">>"+histName_data_FR_ee, cuts_FR_ee + HiggsSRCuts);
    std::cout << samples_ee_UL2016preVFP_DC[0] << " -> Data Driven FR ee (before Nor.) = " << hist_FR_ee->GetEntries() << " : " << hist_FR_ee->GetSumOfWeights() << std::endl;
    hist_FR_ee->Scale(reducibleYield_ee/hist_FR_ee->GetSumOfWeights());
    if(ZfinalStates == "ee")
        hist_FR = hist_FR_ee;
    else
    {
        if(ZfinalStates == "mm")
            hist_FR = hist_FR;
        else
            hist_FR->Add(hist_FR,hist_FR_ee);
    }
    
        //Data Driven shapes: option weighted FR from FR region
        bool weightFRshape = false;
        //updated take weighted shapes
        hist_DR[1]->Add(hist_DR[1],hist_DR[0]);
        hist_DR[2]->Add(hist_DR[2],hist_DR[1]);
        hist_DR[3]->Add(hist_DR[3],hist_DR[2]);
        hist_DR[4]->Add(hist_DR[4],hist_DR[3],-1,1);
        hist_DR[5]->Add(hist_DR[5],hist_DR[4],-1,1);
        hist_DR[6]->Add(hist_DR[6],hist_DR[5],-1,1);
        hist_DR[7]->Add(hist_DR[7],hist_DR[6],-1,1);
        hist_DR[8]->Add(hist_DR[8],hist_DR[7],-1,1);
        hist_DR[9]->Add(hist_DR[9],hist_DR[8],-1,1);
        hist_DR[10]->Add(hist_DR[10],hist_DR[9]);
        hist_DR[11]->Add(hist_DR[11],hist_DR[10]);
        hist_DR[12]->Add(hist_DR[12],hist_DR[11]);
        hist_DR_ee[1]->Add(hist_DR_ee[1],hist_DR_ee[0]);
        hist_DR_ee[2]->Add(hist_DR_ee[2],hist_DR_ee[1]);
        hist_DR_ee[3]->Add(hist_DR_ee[3],hist_DR_ee[2]);
        hist_DR_ee[4]->Add(hist_DR_ee[4],hist_DR_ee[3],-1,1);
        hist_DR_ee[5]->Add(hist_DR_ee[5],hist_DR_ee[4],-1,1);
        hist_DR_ee[6]->Add(hist_DR_ee[6],hist_DR_ee[5],-1,1);
        hist_DR_ee[7]->Add(hist_DR_ee[7],hist_DR_ee[6],-1,1);
        hist_DR_ee[8]->Add(hist_DR_ee[8],hist_DR_ee[7],-1,1);
        hist_DR_ee[9]->Add(hist_DR_ee[9],hist_DR_ee[8],-1,1);
        hist_DR_ee[10]->Add(hist_DR_ee[10],hist_DR_ee[9]);
        hist_DR_ee[11]->Add(hist_DR_ee[11],hist_DR_ee[10]);
        hist_DR_ee[12]->Add(hist_DR_ee[12],hist_DR_ee[11]);
        hist_DR_ee[13]->Add(hist_DR_ee[13],hist_DR_ee[12],-1,1);
        
        if(weightFRshape)
        {
            if(ZfinalStates == "ee")
                hist_FR = hist_DR_ee[13];
            else
            {
                if(ZfinalStates == "mm")
                    hist_FR = hist_DR[13];
                else
                    hist_FR->Add(hist_DR[13],hist_DR_ee[13]);
            }
        }
        std::cout << "reducible entries and yield: " << hist_FR->GetEntries() << " : " << hist_FR->GetSumOfWeights() << std::endl;
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

        TFile * file = new TFile("./Combineworkspace/shapes/produceArea/AZh_"+era+"_"+ZfinalStates+finalState+"_"+varName+"_Central_"+CatName+".root","recreate");
    
        file->mkdir(ZfinalStates+finalState);
        file->cd(ZfinalStates+finalState);
        TH1D * data_obs = (TH1D*)hist[0]->Clone("data_obs");
        TH1D * DY = (TH1D*)hist[1]->Clone("DY");
        TH1D * ggZZ = (TH1D*)hist[2]->Clone("ggZZ");
        TH1D * ZZ = (TH1D*)hist[8]->Clone("ZZ");
        TH1D * TTZ = (TH1D*)hist[10]->Clone("TTZ");
        TH1D * TTW = (TH1D*)hist[11]->Clone("TTW");
        TH1D * TT = (TH1D*)hist[12]->Clone("TT");
        TH1D * VVV = (TH1D*)hist[15]->Clone("VVV");
        TH1D * WZ = (TH1D*)hist[19]->Clone("WZ");
        TH1D * ggHtt = (TH1D*)hist[21]->Clone("ggHtt");
        TH1D * VBFHtt = (TH1D*)hist[22]->Clone("VBFHtt");
        TH1D * WHtt = (TH1D*)hist[23]->Clone("WHtt");
        TH1D * ZHtt = (TH1D*)hist[25]->Clone("ZHtt");
        TH1D * TTHtt = (TH1D*)hist[26]->Clone("TTHtt");
        TH1D * ggHWW = (TH1D*)hist[27]->Clone("ggHWW");
        TH1D * VBFHWW = (TH1D*)hist[28]->Clone("VBFHWW");
        TH1D * WHWW = (TH1D*)hist[29]->Clone("WHWW");
        TH1D * ZHWW = (TH1D*)hist[31]->Clone("ZHWW");
        TH1D * ggZHWW = (TH1D*)hist[32]->Clone("ggZHWW");
        TH1D * ggHZZ = (TH1D*)hist[33]->Clone("ggHZZ");

        TH1D * Reducible = (TH1D*)hist_SS->Clone("Reducible");
        if(UseFRShape)
            Reducible = (TH1D*)hist_FR->Clone("Reducible");
        TH1D * AZh225 = (TH1D*)hist[34]->Clone("AZh225");
        TH1D * AZh250 = (TH1D*)hist[35]->Clone("AZh250");
        TH1D * AZh275 = (TH1D*)hist[36]->Clone("AZh275");
        TH1D * AZh300 = (TH1D*)hist[37]->Clone("AZh300");
        TH1D * AZh325 = (TH1D*)hist[38]->Clone("AZh325");
        TH1D * AZh350 = (TH1D*)hist[39]->Clone("AZh350");
        TH1D * AZh375 = (TH1D*)hist[40]->Clone("AZh375");
        TH1D * AZh400 = (TH1D*)hist[41]->Clone("AZh400");
        TH1D * AZh450 = (TH1D*)hist[42]->Clone("AZh450");
        TH1D * AZh500 = (TH1D*)hist[43]->Clone("AZh500");
        TH1D * AZh600 = (TH1D*)hist[44]->Clone("AZh600");
        TH1D * AZh700 = (TH1D*)hist[45]->Clone("AZh700");
        TH1D * AZh750 = (TH1D*)hist[46]->Clone("AZh750");
        TH1D * AZh800 = (TH1D*)hist[47]->Clone("AZh800");
        TH1D * AZh900 = (TH1D*)hist[48]->Clone("AZh900");
        TH1D * AZh1000 = (TH1D*)hist[49]->Clone("AZh1000");
        TH1D * AZh1200 = (TH1D*)hist[50]->Clone("AZh1200");
        TH1D * AZh1400 = (TH1D*)hist[51]->Clone("AZh1400");
        TH1D * AZh1600 = (TH1D*)hist[52]->Clone("AZh1600");
        TH1D * AZh1800 = (TH1D*)hist[53]->Clone("AZh1800");
        TH1D * AZh2000 = (TH1D*)hist[54]->Clone("AZh2000");

        TH1D * BBAZh225 = (TH1D*)hist[55]->Clone("BBAZh225");
        TH1D * BBAZh250 = (TH1D*)hist[56]->Clone("BBAZh250");
        TH1D * BBAZh275 = (TH1D*)hist[57]->Clone("BBAZh275");
        TH1D * BBAZh300 = (TH1D*)hist[58]->Clone("BBAZh300");
        TH1D * BBAZh325 = (TH1D*)hist[59]->Clone("BBAZh325");
        TH1D * BBAZh350 = (TH1D*)hist[60]->Clone("BBAZh350");
        TH1D * BBAZh375 = (TH1D*)hist[61]->Clone("BBAZh375");
        TH1D * BBAZh400 = (TH1D*)hist[62]->Clone("BBAZh400");
        TH1D * BBAZh450 = (TH1D*)hist[63]->Clone("BBAZh450");
        TH1D * BBAZh500 = (TH1D*)hist[64]->Clone("BBAZh500");
        TH1D * BBAZh600 = (TH1D*)hist[65]->Clone("BBAZh600");
        TH1D * BBAZh700 = (TH1D*)hist[66]->Clone("BBAZh700");
        TH1D * BBAZh750 = (TH1D*)hist[67]->Clone("BBAZh750");
        TH1D * BBAZh800 = (TH1D*)hist[68]->Clone("BBAZh800");
        TH1D * BBAZh900 = (TH1D*)hist[69]->Clone("BBAZh900");
        TH1D * BBAZh1000 = (TH1D*)hist[70]->Clone("BBAZh1000");
        TH1D * BBAZh1200 = (TH1D*)hist[71]->Clone("BBAZh1200");
        TH1D * BBAZh1400 = (TH1D*)hist[72]->Clone("BBAZh1400");
        TH1D * BBAZh1600 = (TH1D*)hist[73]->Clone("BBAZh1600");
        TH1D * BBAZh1800 = (TH1D*)hist[74]->Clone("BBAZh1800");
        TH1D * BBAZh2000 = (TH1D*)hist[75]->Clone("BBAZh2000");

        float totData = data_obs->GetSumOfWeights();
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
    
        float totReducible = Reducible->GetSumOfWeights();
        float totBkg = totDY+totggZZ+totZZ+totTTW+totTTZ+totTT+totVVV+totWZ+totSMH+totReducible;

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
        std::cout << "Reducible: " << Reducible->GetSumOfWeights() << std::endl;

        file->Write();
        file->Close();
    }

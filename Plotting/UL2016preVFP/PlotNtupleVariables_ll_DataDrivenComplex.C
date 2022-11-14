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
void PlotNtupleVariables_ll_DataDrivenComplex(
                                 TString varName= "Mass_svc",
                                 TString xtitle = "m^{c}_{ll#tau#tau} [GeV]",
                                 TString ytitle = "Events",
                                 float xLower = 0,
                                 float xUpper = 200,
                                 int numberofbins = 20,
                                 bool showSignal = true,
                                 bool logY = false,
                                 bool legLeft = false,
                                 bool blindData = true,
                                 TString ZfinalStates = "", //options: "all", "ee", "mm"
                                 TString finalState = "tt",
                                 TString extraCuts = "")
{
    bool checkSS = false;//swtich to check SS region
    bool UseFRShape = false;
    bool logX = false;
    bool UseDataCardBinning = false;
	SetStyle();

	double xsec[44] = { //in pico barn (pb)
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
        0.0262*0.1,//()GGZH H -> WW, Z->LL
        0.0128,//GG H -> ZZ -> 4L
        0.026, //AZh (225)
        0.01, //AZh (300)
        0.1, //AZh (400)
        0.1, //AZh (700)
        0.1, //AZh (1000)
        0.004, //BBAZh (225)
        0.0025, //BBAZh (300)
        0.1, //BBAZh (400)
        0.1, //BBAZh (700)
        0.1 //BBAZh (1000)
    };
    
        float lumi = 19500;
        //Deal with bins and binning
        float xMin = xLower;
        float xMax = xUpper;
        int nBins = numberofbins;
        float bins[100];
        float nBinsX_AZh[16] = {200,220,240,260,280,300,320,340,360,380,400,450,550,700,1000,2400};

        if(UseDataCardBinning)
        {
            nBins = 15;
            for (int iB=0; iB<=nBins; ++iB)
            {
                bins[iB] = nBinsX_AZh[iB];
            }
        }
        else
        {
            float binWidth = (xMax-xMin)/float(nBins);
            for (int iB=0; iB<=nBins; ++iB)
                bins[iB] = xMin + float(iB)*binWidth;
        }
        int nSamples = 44;
        TH1D * hist[44];   //This is "hist_mumu"
        TH1D * hist_ee[44];
    
        //Data-driven Normalization Hist
        TH1D * hist_DR[14];
        TH1D * hist_DR_ee[14];

        //Data-driven Shape Hist -- Same Sign region
        TH1D * hist_SS;
        TH1D * hist_SS_ee;

        //inistiating cuts
        TString cuts[44];
    
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
            if(i<34)//for background
                cuts[i] = MCCommonWeight + SRMCWeight +"(" + Leg3ID + " && " + Leg4ID + "&& gen_match_1!=6 && gen_match_2!=6 && gen_match_3!=6 && gen_match_4!=6)";
            else//for signals
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
                std::cout << "EE channel " <<i<<":"+samples_ee_UL2016preVFP[i]<<std::endl;
                //TFile * file = new TFile("../NTuples/v7/"+samples[i]+".root");
                TFile * file_ee = new TFile("./NTuples_ee"+finalState+"/"+samples_ee_UL2016preVFP[i]+".root");
                TH1D * histWeightsH_ee = (TH1D*)file_ee->Get("histWeightsH");
                TTree * tree_ee = (TTree*)file_ee->Get("SynTree");
                double normaliza_ee = xsec[i]*lumi/histWeightsH_ee->GetSumOfWeights();
                TString histName_ee = samples_ee_UL2016preVFP[i] + "_"+varName + "_ee";
                if(UseDataCardBinning)
                    hist_ee[i] = new TH1D(histName_ee,"",nBins,bins);
                else
                    hist_ee[i] = new TH1D(histName_ee,"",nBins,xMin,xMax);
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
                std::cout << "MuMu channel " <<i<<":"+samples_mumu_UL2016preVFP[i]<<std::endl;
                TFile * file_mumu = new TFile("./NTuples_mm"+finalState+"/"+samples_mumu_UL2016preVFP[i]+".root");
                TH1D * histWeightsH_mumu = (TH1D*)file_mumu->Get("histWeightsH");
                TTree * tree_mumu = (TTree*)file_mumu->Get("SynTree");
                double normaliza_mumu = xsec[i]*lumi/histWeightsH_mumu->GetSumOfWeights();
                TString histName_mumu = samples_mumu_UL2016preVFP[i] + "_"+varName+ "_mumu";
                if(UseDataCardBinning)
                    hist[i] = new TH1D(histName_mumu,"",nBins,bins);
                else
                    hist[i] = new TH1D(histName_mumu,"",nBins,xMin,xMax);
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
        1.227*1012.5,
        1.227*332.8,
        1.227*101.8,
        1.227*54.8};//2016 values
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

        if(UseDataCardBinning)
            histZ_ee[i] = new TH1D(histNameZ_ee,"",nBins,bins);
        else
            histZ_ee[i] = new TH1D(histNameZ_ee,"",nBins,xMin,xMax);
        
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

        if(UseDataCardBinning)
            histZ[i] = new TH1D(histNameZ_mumu,"",nBins,bins);
        else
            histZ[i] = new TH1D(histNameZ_mumu,"",nBins,xMin,xMax);
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
    
    TFile * file_data = new TFile("./NTuples_mm"+finalState+"/"+samples_mumu_UL2016preVFP[0]+".root");
    TTree * tree_data = (TTree*)file_data->Get("SynTree");
    
    TString histName_dataCat1 = samples_mumu_UL2016preVFP[0] + "_Cat1_"+varName;
    TString histName_dataCat2 = samples_mumu_UL2016preVFP[0] + "_Cat2_"+varName;
    TString histName_dataCat3 = samples_mumu_UL2016preVFP[0] + "_Cat3_"+varName;
    TString histName_dataCat4 = samples_mumu_UL2016preVFP[0] + "_Cat4_"+varName;
    TString histName_dataCat12 = samples_mumu_UL2016preVFP[0] + "_Cat12_"+varName;
    TString histName_dataCat13 = samples_mumu_UL2016preVFP[0] + "_Cat13_"+varName;
    TString histName_dataCat14 = samples_mumu_UL2016preVFP[0] + "_Cat14_"+varName;
    TString histName_dataCat23 = samples_mumu_UL2016preVFP[0] + "_Cat23_"+varName;
    TString histName_dataCat24 = samples_mumu_UL2016preVFP[0] + "_Cat24_"+varName;
    TString histName_dataCat34 = samples_mumu_UL2016preVFP[0] + "_Cat34_"+varName;
    TString histName_dataCat123 = samples_mumu_UL2016preVFP[0] + "_Cat123_"+varName;
    TString histName_dataCat124 = samples_mumu_UL2016preVFP[0] + "_Cat124_"+varName;
    TString histName_dataCat234 = samples_mumu_UL2016preVFP[0] + "_Cat234_"+varName;
    TString histName_dataCat1234 = samples_mumu_UL2016preVFP[0] + "_Cat1234_"+varName;
    
    if(UseDataCardBinning)
    {
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
    }
    else
    {
        hist_DR[0] = new TH1D(histName_dataCat1,"",nBins,xMin,xMax);
        hist_DR[1] = new TH1D(histName_dataCat2,"",nBins,xMin,xMax);
        hist_DR[2] = new TH1D(histName_dataCat3,"",nBins,xMin,xMax);
        hist_DR[3] = new TH1D(histName_dataCat4,"",nBins,xMin,xMax);
        hist_DR[4] = new TH1D(histName_dataCat12,"",nBins,xMin,xMax);
        hist_DR[5] = new TH1D(histName_dataCat13,"",nBins,xMin,xMax);
        hist_DR[6] = new TH1D(histName_dataCat14,"",nBins,xMin,xMax);
        hist_DR[7] = new TH1D(histName_dataCat23,"",nBins,xMin,xMax);
        hist_DR[8] = new TH1D(histName_dataCat24,"",nBins,xMin,xMax);
        hist_DR[9] = new TH1D(histName_dataCat34,"",nBins,xMin,xMax);
        hist_DR[10] = new TH1D(histName_dataCat123,"",nBins,xMin,xMax);
        hist_DR[11] = new TH1D(histName_dataCat124,"",nBins,xMin,xMax);
        hist_DR[12] = new TH1D(histName_dataCat234,"",nBins,xMin,xMax);
        hist_DR[13] = new TH1D(histName_dataCat1234,"",nBins,xMin,xMax);
    }
    
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
    
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat1 = " << hist_DR[0]->GetEntries() << " : " << hist_DR[0]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat2 = " << hist_DR[1]->GetEntries() << " : " << hist_DR[1]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat3 = " << hist_DR[2]->GetEntries() << " : " << hist_DR[2]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat4 = " << hist_DR[3]->GetEntries() << " : " << hist_DR[3]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat12 = " << hist_DR[4]->GetEntries() << " : " << hist_DR[4]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat13 = " << hist_DR[5]->GetEntries() << " : " << hist_DR[5]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat14 = " << hist_DR[6]->GetEntries() << " : " << hist_DR[6]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat23 = " << hist_DR[7]->GetEntries() << " : " << hist_DR[7]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat24 = " << hist_DR[8]->GetEntries() << " : " << hist_DR[8]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat34 = " << hist_DR[9]->GetEntries() << " : " << hist_DR[9]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat123 = " << hist_DR[10]->GetEntries() << " : " << hist_DR[10]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat124 = " << hist_DR[11]->GetEntries() << " : " << hist_DR[11]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat234 = " << hist_DR[12]->GetEntries() << " : " << hist_DR[12]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven Cat1234 = " << hist_DR[13]->GetEntries() << " : " << hist_DR[13]->GetSumOfWeights() << std::endl;

    float reducibleYield_mumu = hist_DR[0]->GetSumOfWeights() + hist_DR[1]->GetSumOfWeights() + hist_DR[2]->GetSumOfWeights() + hist_DR[3]->GetSumOfWeights() - hist_DR[4]->GetSumOfWeights() - hist_DR[5]->GetSumOfWeights() - hist_DR[6]->GetSumOfWeights() - hist_DR[7]->GetSumOfWeights() - hist_DR[8]->GetSumOfWeights() - hist_DR[9]->GetSumOfWeights() + hist_DR[10]->GetSumOfWeights() + hist_DR[11]->GetSumOfWeights() + hist_DR[12]->GetSumOfWeights() - hist_DR[13]->GetSumOfWeights();
    
    float reducibleYield_mumu_old = hist_DR[2]->GetSumOfWeights() + hist_DR[3]->GetSumOfWeights() - hist_DR[9]->GetSumOfWeights();
    hist_DR[13]->Add(hist_DR[13],hist_DR[12],-1,1);

    std::cout << " Reducible Yields(mumu) = " << reducibleYield_mumu << std::endl;

    TFile * file_data_ee = new TFile("./NTuples_ee"+finalState+"/"+samples_ee_UL2016preVFP[0]+".root");
    TTree * tree_data_ee = (TTree*)file_data_ee->Get("SynTree");
    
    TString histName_dataCat1_ee = samples_ee_UL2016preVFP[0] + "_Cat1_"+varName;
    TString histName_dataCat2_ee = samples_ee_UL2016preVFP[0] + "_Cat2_"+varName;
    TString histName_dataCat3_ee = samples_ee_UL2016preVFP[0] + "_Cat3_"+varName;
    TString histName_dataCat4_ee = samples_ee_UL2016preVFP[0] + "_Cat4_"+varName;
    TString histName_dataCat12_ee = samples_ee_UL2016preVFP[0] + "_Cat12_"+varName;
    TString histName_dataCat13_ee = samples_ee_UL2016preVFP[0] + "_Cat13_"+varName;
    TString histName_dataCat14_ee = samples_ee_UL2016preVFP[0] + "_Cat14_"+varName;
    TString histName_dataCat23_ee = samples_ee_UL2016preVFP[0] + "_Cat23_"+varName;
    TString histName_dataCat24_ee = samples_ee_UL2016preVFP[0] + "_Cat24_"+varName;
    TString histName_dataCat34_ee = samples_ee_UL2016preVFP[0] + "_Cat34_"+varName;
    TString histName_dataCat123_ee = samples_ee_UL2016preVFP[0] + "_Cat123_"+varName;
    TString histName_dataCat124_ee = samples_ee_UL2016preVFP[0] + "_Cat124_"+varName;
    TString histName_dataCat234_ee = samples_ee_UL2016preVFP[0] + "_Cat234_"+varName;
    TString histName_dataCat1234_ee = samples_ee_UL2016preVFP[0] + "_Cat1234_"+varName;
    
    if(UseDataCardBinning)
    {
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
    }
    else
    {
        hist_DR_ee[0] = new TH1D(histName_dataCat1_ee,"",nBins,xMin,xMax);
        hist_DR_ee[1] = new TH1D(histName_dataCat2_ee,"",nBins,xMin,xMax);
        hist_DR_ee[2] = new TH1D(histName_dataCat3_ee,"",nBins,xMin,xMax);
        hist_DR_ee[3] = new TH1D(histName_dataCat4_ee,"",nBins,xMin,xMax);
        hist_DR_ee[4] = new TH1D(histName_dataCat12_ee,"",nBins,xMin,xMax);
        hist_DR_ee[5] = new TH1D(histName_dataCat13_ee,"",nBins,xMin,xMax);
        hist_DR_ee[6] = new TH1D(histName_dataCat14_ee,"",nBins,xMin,xMax);
        hist_DR_ee[7] = new TH1D(histName_dataCat23_ee,"",nBins,xMin,xMax);
        hist_DR_ee[8] = new TH1D(histName_dataCat24_ee,"",nBins,xMin,xMax);
        hist_DR_ee[9] = new TH1D(histName_dataCat34_ee,"",nBins,xMin,xMax);
        hist_DR_ee[10] = new TH1D(histName_dataCat123_ee,"",nBins,xMin,xMax);
        hist_DR_ee[11] = new TH1D(histName_dataCat124_ee,"",nBins,xMin,xMax);
        hist_DR_ee[12] = new TH1D(histName_dataCat234_ee,"",nBins,xMin,xMax);
        hist_DR_ee[13] = new TH1D(histName_dataCat1234_ee,"",nBins,xMin,xMax);
    }
    
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
    
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat1 = " << hist_DR_ee[0]->GetEntries() << " : " << hist_DR_ee[0]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat2 = " << hist_DR_ee[1]->GetEntries() << " : " << hist_DR_ee[1]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat3 = " << hist_DR_ee[2]->GetEntries() << " : " << hist_DR_ee[2]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat4 = " << hist_DR_ee[3]->GetEntries() << " : " << hist_DR_ee[3]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat12 = " << hist_DR_ee[4]->GetEntries() << " : " << hist_DR_ee[4]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat13 = " << hist_DR_ee[5]->GetEntries() << " : " << hist_DR_ee[5]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat14 = " << hist_DR_ee[6]->GetEntries() << " : " << hist_DR_ee[6]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat23 = " << hist_DR_ee[7]->GetEntries() << " : " << hist_DR_ee[7]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat24 = " << hist_DR_ee[8]->GetEntries() << " : " << hist_DR_ee[8]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat34 = " << hist_DR_ee[9]->GetEntries() << " : " << hist_DR_ee[9]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat123 = " << hist_DR_ee[10]->GetEntries() << " : " << hist_DR_ee[10]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat124 = " << hist_DR_ee[11]->GetEntries() << " : " << hist_DR_ee[11]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat234 = " << hist_DR_ee[12]->GetEntries() << " : " << hist_DR_ee[12]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven Cat1234 = " << hist_DR_ee[13]->GetEntries() << " : " << hist_DR_ee[13]->GetSumOfWeights() << std::endl;

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
    TString histName_data_SS = samples_mumu_UL2016preVFP[0] + "_SS_"+varName;
    if(UseDataCardBinning)
        hist_SS = new TH1D(histName_data_SS,"",nBins,bins);
    else
        hist_SS = new TH1D(histName_data_SS,"",nBins,xMin,xMax);
    hist_SS->Sumw2();
    tree_data->Draw(varName+">>"+histName_data_SS, cutsRelaxed + HiggsSSCuts + extraCuts);
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven SS mumu (before Nor.) = " << hist_SS->GetEntries() << " : " << hist_SS->GetSumOfWeights() << std::endl;
    hist_SS->Scale(reducibleYield_mumu/hist_SS->GetSumOfWeights());

    TString histName_data_SS_ee = samples_ee_UL2016preVFP[0] + "_SS_"+varName;
    if(UseDataCardBinning)
        hist_SS_ee = new TH1D(histName_data_SS_ee,"",nBins,bins);
    else
        hist_SS_ee = new TH1D(histName_data_SS_ee,"",nBins,xMin,xMax);
    hist_SS_ee->Sumw2();
    tree_data_ee->Draw(varName+">>"+histName_data_SS_ee, cutsRelaxed + HiggsSSCuts + extraCuts);
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven SS ee (before Nor.) = " << hist_SS_ee->GetEntries() << " : " << hist_SS_ee->GetSumOfWeights() << std::endl;
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
    
    TString histName_data_FR = samples_mumu_UL2016preVFP[0] + "_FR_"+varName;
    if(UseDataCardBinning)
        hist_FR = new TH1D(histName_data_FR,"",nBins,bins);
    else
        hist_FR = new TH1D(histName_data_FR,"",nBins,xMin,xMax);
    hist_FR->Sumw2();
    tree_data->Draw(varName+">>"+histName_data_FR, cuts_FR_mumu + HiggsSRCuts + extraCuts);
    std::cout << samples_mumu_UL2016preVFP[0] << " -> Data Driven FR mumu (before Nor.) = " << hist_FR->GetEntries() << " : " << hist_FR->GetSumOfWeights() << std::endl;
    hist_FR->Scale(reducibleYield_mumu/hist_FR->GetSumOfWeights());

    TString histName_data_FR_ee = samples_ee_UL2016preVFP[0] + "_FR_"+varName;
    if(UseDataCardBinning)
        hist_FR_ee = new TH1D(histName_data_FR_ee,"",nBins,bins);
    else
        hist_FR_ee = new TH1D(histName_data_FR_ee,"",nBins,xMin,xMax);
    hist_FR_ee->Sumw2();
    tree_data_ee->Draw(varName+">>"+histName_data_FR_ee, cuts_FR_ee + HiggsSRCuts + extraCuts);
    std::cout << samples_ee_UL2016preVFP[0] << " -> Data Driven FR ee (before Nor.) = " << hist_FR_ee->GetEntries() << " : " << hist_FR_ee->GetSumOfWeights() << std::endl;
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
        for (int iH=3; iH<10; ++iH)
        {
            hist[2]->Add(hist[2],hist[iH]);
        }
        
        // adding up TT backgrounds
        for (int iH=11; iH<15; ++iH)
        {
            hist[10]->Add(hist[10],hist[iH]);
        }
    
        // adding up VVV backgrounds
        for (int iH=16; iH<19; ++iH)
        {
            hist[15]->Add(hist[15],hist[iH]);
        }
    
        // (18+19) is WZ
        for (int iH=20; iH<21; ++iH)
        {
            hist[19]->Add(hist[19],hist[iH]);
        }
    
        // adding up SM H 125 backgrounds
        for(int iH=22; iH<34; ++iH)
        {
            hist[21]->Add(hist[21],hist[iH]);
        }
    
        TString PlotShapeName = "";
        if(checkSS)
            PlotShapeName = "_SS";
        if(extraCuts.Contains("nbtag>=1"))
            PlotShapeName += "_btag";
        if(extraCuts.Contains("nbtag==0"))
            PlotShapeName += "_0btag";

        TFile * file = new TFile("./PlottingShapes/AZh_"+ZfinalStates+finalState+"_"+varName+PlotShapeName+".root","recreate");
        file->mkdir(ZfinalStates+finalState);
        file->cd(ZfinalStates+finalState);

        TH1D * data_obs = (TH1D*)hist[0]->Clone("data_obs");
        TH1D * DY = (TH1D*)hist[1]->Clone("DY");
        TH1D * ZZ = (TH1D*)hist[2]->Clone("ZZ");
        TH1D * TT = (TH1D*)hist[10]->Clone("TT");
        TH1D * VVV = (TH1D*)hist[15]->Clone("VVV");
        TH1D * WZ = (TH1D*)hist[19]->Clone("WZ");
        TH1D * SMH = (TH1D*)hist[21]->Clone("SMH");
        TH1D * Reducible = (TH1D*)hist_SS->Clone("Reducible");
        if(UseFRShape)
            Reducible = (TH1D*)hist_FR->Clone("Reducible");
        TH1D * AZh = (TH1D*)hist[34]->Clone("AZh");
        TH1D * AZh300 = (TH1D*)hist[35]->Clone("AZh300");
        TH1D * AZh400 = (TH1D*)hist[36]->Clone("AZh400");
        TH1D * AZh700 = (TH1D*)hist[37]->Clone("AZh700");
        TH1D * AZh1000 = (TH1D*)hist[38]->Clone("AZh1000");
        TH1D * BBAZh225 = (TH1D*)hist[39]->Clone("BBAZh225");
        TH1D * BBAZh300 = (TH1D*)hist[40]->Clone("BBAZh300");
        TH1D * BBAZh400 = (TH1D*)hist[41]->Clone("BBAZh400");
        TH1D * BBAZh700 = (TH1D*)hist[42]->Clone("BBAZh700");
        TH1D * BBAZh1000 = (TH1D*)hist[43]->Clone("BBAZh1000");
    
        file->Write();

        float totData = data_obs->GetSumOfWeights();
        float totDY = DY->GetSumOfWeights();
        float totZZ = ZZ->GetSumOfWeights();
        float totTT = TT->GetSumOfWeights();
        float totVVV = VVV->GetSumOfWeights();
        float totWZ = WZ->GetSumOfWeights();
        float totSMH = SMH->GetSumOfWeights();
        float totReducible = Reducible->GetSumOfWeights();
        float totBkg = totDY+totZZ+totTT+totVVV+totWZ+totSMH+totReducible;

        if(!blindData)
            std::cout << "data : " << totData << std::endl;
        std::cout << "Total bkg : " << totBkg << std::endl;
        std::cout << "GGAZh(225) : " << AZh->GetSumOfWeights() << std::endl;
        std::cout << "GGAZh(300) : " << AZh300->GetSumOfWeights() << std::endl;
        std::cout << "GGAZh(400) : " << AZh400->GetSumOfWeights() << std::endl;
        std::cout << "GGAZh(700) : " << AZh700->GetSumOfWeights() << std::endl;
        std::cout << "GGAZh(1000) : " << AZh1000->GetSumOfWeights() << std::endl;
        std::cout << "BBAZh(225) : " << BBAZh225->GetSumOfWeights() << std::endl;
        std::cout << "BBAZh(300) : " << BBAZh300->GetSumOfWeights() << std::endl;
        std::cout << "BBAZh(400) : " << BBAZh400->GetSumOfWeights() << std::endl;
        std::cout << "BBAZh(700) : " << BBAZh700->GetSumOfWeights() << std::endl;
        std::cout << "BBAZh(1000) : " << BBAZh1000->GetSumOfWeights() << std::endl;
        std::cout << "DY : " << DY->GetSumOfWeights() << std::endl;
        std::cout << "ZZ  : " << ZZ->GetSumOfWeights() << std::endl;
        std::cout << "TT  : " << TT->GetSumOfWeights() << std::endl;
        std::cout << "VVV : " << VVV->GetSumOfWeights() << std::endl;
        std::cout << "WZ : " << WZ->GetSumOfWeights() << std::endl;
        std::cout << "SM H : " << SMH->GetSumOfWeights() << std::endl;
        std::cout << "Reducible: " << Reducible->GetSumOfWeights() << std::endl;
    
        if(checkSS)//for getting systematics from non-closure
        {
            float data_error = 0;
            float data_error_r = 0;
            for(int iB=1; iB<=nBins; ++iB)
            {
                data_error += data_obs->GetBinError(iB)*data_obs->GetBinError(iB);
            }
            data_error = sqrt(data_error);
            float datasubtractedMC = totData - (totDY+totZZ+totTT+totVVV+totWZ+totSMH);
            float difference = abs(datasubtractedMC - totReducible)/((totReducible + datasubtractedMC)/2);
            data_error_r = data_error/totData;
            std::cout << "data statistics error: " << data_error_r*100 << "%" << std::endl;
            std::cout << "Non closure difference: " << difference*100 << "%" << std::endl;
            std::cout << "Non closure error: " << sqrt(data_error_r*data_error_r+difference*difference) *100  << "%" << std::endl;
        }

        TH1D * dummy = (TH1D*)DY->Clone("dummy");
        float errLumi = 0.03;
        float errDY=0.1;
        float errVVV = 0.15;
        float errTT = 0.15;
        float errWZ = 0.15;
        float errReducible = 0.05;

        for (int iB=1; iB<=nBins; ++iB)
        {
            float eVVV = errVVV*VVV->GetBinContent(iB);
            float eDY = errDY*DY->GetBinContent(iB);
            float eWZ = errWZ*WZ->GetBinContent(iB);
            float eTT = errTT*TT->GetBinContent(iB);
            float eReducible = errReducible*Reducible->GetBinContent(iB);
            float err2 = eVVV*eVVV + eDY*eDY + eTT*eTT+eWZ*eWZ+eReducible*eReducible;
            float errTot = TMath::Sqrt(err2);
            dummy->SetBinError(iB,errTot);
            AZh300->SetBinError(iB,0);
            BBAZh300->SetBinError(iB,0);
        }
        VVV->Add(VVV,WZ);
        TT->Add(TT,VVV);
        ZZ->Add(ZZ,TT);
        SMH->Add(SMH,ZZ);
        DY->Add(DY,SMH);
        Reducible->Add(Reducible,DY);

        TH1D * bkgdErr = (TH1D*)Reducible->Clone("bkgdErr");
        bkgdErr->SetFillStyle(3013);
        bkgdErr->SetFillColor(1);
        bkgdErr->SetMarkerStyle(21);
        bkgdErr->SetMarkerSize(0);

        for (int iB=1; iB<=nBins; ++iB)
        {
            SMH->SetBinError(iB,0);
            DY->SetBinError(iB,0);
            ZZ->SetBinError(iB,0);
            TT->SetBinError(iB,0);
            VVV->SetBinError(iB,0);
            WZ->SetBinError(iB,0);
            Reducible->SetBinError(iB,0);
            float eStat =  bkgdErr->GetBinError(iB);
            float X = bkgdErr->GetBinContent(iB);
            float eLumi = errLumi * X;
            float eBkg = dummy->GetBinError(iB);
            float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eBkg*eBkg);
            bkgdErr->SetBinError(iB,Err);
        }
        //Colors
        Int_t colorDY = TColor::GetColor("#ffcc66");
        Int_t colorSMH = TColor::GetColor("#FFCCFF");
        Int_t colorZZ = TColor::GetColor("#4496C8");
        Int_t colorTT = TColor::GetColor("#9999CC");
        Int_t colorVVV = TColor::GetColor("#6F2D35");
        Int_t colorWZ = TColor::GetColor("#DE5A6A");
        Int_t colorReducible = TColor::GetColor("#07752e");

        InitData(data_obs);
        InitHist(WZ,"","",colorWZ,1001);
        InitHist(VVV,"","",colorVVV,1001);
        InitHist(TT,"","",colorTT,1001);
        InitHist(ZZ,"","",colorZZ,1001);
        InitHist(SMH,"","",colorSMH,1001);
        InitHist(DY,"","",colorDY,1001);
        InitHist(Reducible,"","",colorReducible,1001);

            data_obs->GetXaxis()->SetTitle(xtitle);
            data_obs->GetYaxis()->SetTitle(ytitle);
            data_obs->GetYaxis()->SetTitleOffset(1.5);
            data_obs->GetYaxis()->SetTitleSize(0.06);
            data_obs->GetXaxis()->SetRangeUser(xLower,xUpper);
            float yUpper = data_obs->GetMaximum();
            if (logY)
            data_obs->GetYaxis()->SetRangeUser(0.1,2*yUpper);
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
    
        if(blindData)
        {
            for (int iB=1; iB<=nBins; ++iB)
            {
                float bkgbin = Reducible->GetBinContent(iB);
                float sigbin = AZh300->GetBinContent(iB)+BBAZh300->GetBinContent(iB);
                if((sigbin/(TMath::Sqrt(bkgbin+(0.09*bkgbin)*(0.09*bkgbin))))>=0.5)
                {
                    data_obs->SetBinError(iB,0);
                    data_obs->SetBinContent(iB,0);
                }
            }
        }
        if(blindData && UseDataCardBinning)
        {
            for(int iB=1; iB<=nBins; ++iB)
            {
                data_obs->SetBinError(iB,0);
                data_obs->SetBinContent(iB,0);
            }
        }

        //Drawing histogram
            data_obs->Draw("e1");
            Reducible->Draw("sameHIST");
            DY->Draw("sameHIST");
            SMH->Draw("sameHIST");
            ZZ->Draw("sameHIST");
            TT->Draw("sameHIST");
            VVV->Draw("sameHIST");
            WZ->Draw("sameHIST");
            data_obs->Draw("e1same");
        bkgdErr->Draw("e2same");
    
        InitSignal(AZh300);
        InitSignal(BBAZh300);
    
        AZh300->SetLineColor(2);
        BBAZh300->SetLineColor(4);
    
        if (showSignal)
        {
            AZh300->Draw("hsame");
            BBAZh300->Draw("hsame");
        }
    
        //Calculating chi2
        float chi2 = 0;
        for (int iB=1; iB<=nBins; ++iB)
        {
            float xData = 0;
            if(!blindData)
                xData = data_obs->GetBinContent(iB);
            float xMC = Reducible->GetBinContent(iB);
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
        leg->AddEntry(Reducible,"Reducible","f");
        leg->AddEntry(DY,"DY","f");
        leg->AddEntry(SMH,"SM H(125)","f");
        leg->AddEntry(ZZ,"ZZ","f");
        leg->AddEntry(TT,"t#bar{t}","f");
        leg->AddEntry(VVV,"VVV","f");
        leg->AddEntry(WZ,"WZ","f");
        if(showSignal)
        {
            leg->AddEntry(AZh300,"GGA(300) #rightarrow Zh, #sigma = 10 fb","l");
            leg->AddEntry(BBAZh300,"BBA(300) #rightarrow Zh, #sigma = 2.5 fb","l");
        }
        leg->Draw();
        TString ZfinalStatesTag = "";
        if(ZfinalStates == "ee")
            ZfinalStatesTag = "ee";
        else
        {
            if(ZfinalStates == "mm")
                ZfinalStatesTag = "mm";
            else
                ZfinalStatesTag = "ll";
        }
        TLatex * cms = new TLatex(0.20,0.94,"CMS L = 19.5 fb^{-1} at #sqrt{s} = 13 TeV, "+ZfinalStatesTag+finalState);

        cms->SetNDC();
        cms->SetTextSize(0.05);
        cms->Draw();
        TLatex * workinprogress = new TLatex(x1Leg,0.55,"Work in progress");
        workinprogress->SetNDC();
        workinprogress->SetTextSize(0.05);
        workinprogress->Draw();
        
        if (logY) upper->SetLogy(true);
        if (logX) upper->SetLogx(true);

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
                float x2 = Reducible->GetBinContent(iB);
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
        if (logX) lower->SetLogx(true);
        canv1->cd();
        canv1->Modified();
        canv1->cd();
        canv1->SetSelected(canv1);

        if(checkSS)
            canv1->Print("./figures/"+varName+"_"+ZfinalStatesTag+finalState+"_SS_DataDriven_2016.pdf","Portrait pdf");
        else
            canv1->Print("./figures/"+varName+"_"+ZfinalStatesTag+finalState+"_DataDriven_2016.pdf","Portrait pdf");

        file->Close();
    }

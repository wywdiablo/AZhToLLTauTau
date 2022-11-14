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

void ProduceAZhDataCards_Syst_FakeRates(
                                 TString era = "UL2018",
                                 TString varName = "Mass_svc",
                                 TString UncerSuffix = "CMS_azh_fr_JetFakeM_2018",
                                 TString Shift = "Up",
                                 TString Measurement = "JetToMu",
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
    
        TString FRZ1 = "fakeRateWeight_Z1";
        TString FRZ2 = "fakeRateWeight_Z2";
        TString FRH1 = "fakeRateWeight_1";
        TString FRH2 = "fakeRateWeight_2";
    
        if(Measurement == "JetToEle")
        {
            if(ZfinalStates == "ee")
            {
                if(Shift == "Up")
                {
                    FRZ1 = "fakeRateWeightUp_Z1";
                    FRZ2 = "fakeRateWeightUp_Z2";
                }
                if(Shift == "Down")
                {
                    FRZ1 = "fakeRateWeightDown_Z1";
                    FRZ2 = "fakeRateWeightDown_Z2";
                }
            }
            if(finalState == "em")
            {
                if(Shift == "Up")
                {
                    FRH1 = "fakeRateWeightUp_1";
                }
                if(Shift == "Down")
                {
                    FRH1 = "fakeRateWeightDown_1";
                }
            }
            if(finalState == "et")
            {
                if(Shift == "Up")
                {
                    FRH1 = "fakeRateWeightUp_1";
                }
                if(Shift == "Down")
                {
                    FRH1 = "fakeRateWeightDown_1";
                }
            }
        }
        if(Measurement == "JetToMu")
        {
            if(ZfinalStates == "mm")
            {
                if(Shift == "Up")
                {
                    FRZ1 = "fakeRateWeightUp_Z1";
                    FRZ2 = "fakeRateWeightUp_Z2";
                }
                if(Shift == "Down")
                {
                    FRZ1 = "fakeRateWeightDown_Z1";
                    FRZ2 = "fakeRateWeightDown_Z2";
                }
            }
            if(finalState == "em")
            {
                if(Shift == "Up")
                {
                    FRH2 = "fakeRateWeightUp_2";
                }
                if(Shift == "Down")
                {
                    FRH2 = "fakeRateWeightDown_2";
                }
            }
            if(finalState == "mt")
            {
                if(Shift == "Up")
                {
                    FRH1 = "fakeRateWeightUp_1";
                }
                if(Shift == "Down")
                {
                    FRH1 = "fakeRateWeightDown_1";
                }
            }
        }
        if(Measurement == "JetToTau")
        {
            if(finalState == "et")
            {
                if(Shift == "Up")
                {
                    FRH2 = "fakeRateWeightUp_2";
                }
                if(Shift == "Down")
                {
                    FRH2 = "fakeRateWeightDown_2";
                }
            }
            if(finalState == "mt")
            {
                if(Shift == "Up")
                {
                    FRH2 = "fakeRateWeightUp_2";
                }
                if(Shift == "Down")
                {
                    FRH2 = "fakeRateWeightDown_2";
                }
            }
            if(finalState == "tt")
            {
                if(Shift == "Up")
                {
                    FRH1 = "fakeRateWeightUp_1";
                    FRH2 = "fakeRateWeightUp_2";
                }
                if(Shift == "Down")
                {
                    FRH1 = "fakeRateWeightDown_1";
                    FRH2 = "fakeRateWeightDown_2";
                }
            }
        }
        eeCat1Cuts = FRZ1+"*(!(" + eeLeg1ID + ") && ("+ eeLeg2ID + ") && (" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat2Cuts = FRZ2+"*((" + eeLeg1ID + ") && !("+ eeLeg2ID + ") && (" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat3Cuts = FRH1+"*((" + eeLeg1ID + ") && (" + eeLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat4Cuts = FRH2+"*((" + eeLeg1ID + ") && (" + eeLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        eeCat12Cuts = FRZ1+"*"+FRZ2+"*(!(" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && (" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat13Cuts = FRZ1+"*"+FRH1+"*(!(" + eeLeg1ID + ") && (" + eeLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat14Cuts = FRZ1+"*"+FRH2+"*(!(" + eeLeg1ID + ") && (" + eeLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        eeCat23Cuts = FRZ2+"*"+FRH1+"*((" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat24Cuts = FRZ2+"*"+FRH2+"*((" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        eeCat34Cuts = FRH1+"*"+FRH2+"*((" + eeLeg1ID + ") && (" + eeLeg2ID + ") && !(" + Leg3ID + ") && !(" + Leg4ID + "))";
        eeCat123Cuts = FRZ1+"*"+FRZ2+"*"+FRH1+"*(!(" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        eeCat124Cuts = FRZ1+"*"+FRZ2+"*"+FRH2+"*(!(" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        eeCat234Cuts = FRZ2+"*"+FRH1+"*"+FRH2+"*((" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && !(" + Leg3ID + ") && !(" + Leg4ID + "))";
        eeCat1234Cuts = FRZ1+"*"+FRZ2+"*"+FRH1+"*"+FRH2+"*(!(" + eeLeg1ID + ") && !(" + eeLeg2ID + ") && !(" + Leg3ID + ") && !(" + Leg4ID + "))";
    
        mumuCat1Cuts = FRZ1+"*(!(" + mumuLeg1ID + ") && ("+ mumuLeg2ID + ") && (" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat2Cuts = FRZ2+"*((" + mumuLeg1ID + ") && !("+ mumuLeg2ID + ") && (" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat3Cuts = FRH1+"*((" + mumuLeg1ID + ") && (" + mumuLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat4Cuts = FRH2+"*((" + mumuLeg1ID + ") && (" + mumuLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat12Cuts = FRZ1+"*"+FRZ2+"*(!(" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && (" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat13Cuts = FRZ1+"*"+FRH1+"*(!(" + mumuLeg1ID + ") && (" + mumuLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat14Cuts = FRZ1+"*"+FRH2+"*(!(" + mumuLeg1ID + ") && (" + mumuLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat23Cuts = FRZ2+"*"+FRH1+"*((" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat24Cuts = FRZ2+"*"+FRH2+"*((" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat34Cuts = FRH1+"*"+FRH2+"*((" + mumuLeg1ID + ") && (" + mumuLeg2ID + ") && !(" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat123Cuts = FRZ1+"*"+FRZ2+"*"+FRH1+"*(!(" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && !(" + Leg3ID + ") && (" + Leg4ID + "))";
        mumuCat124Cuts = FRZ1+"*"+FRZ2+"*"+FRH2+"*(!(" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && (" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat234Cuts = FRZ2+"*"+FRH1+"*"+FRH2+"*((" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && !(" + Leg3ID + ") && !(" + Leg4ID + "))";
        mumuCat1234Cuts = FRZ1+"*"+FRZ2+"*"+FRH1+"*"+FRH2+"*(!(" + mumuLeg1ID + ") && !(" + mumuLeg2ID + ") && !(" + Leg3ID + ") && !(" + Leg4ID + "))";
            

    //Start Data driven method
    TFile * file_data = new TFile("./NTuples_mm"+finalState+"/"+samples_mumu_UL2018_DC[0]+".root");
    TTree * tree_data = (TTree*)file_data->Get("SynTree");
    
    TString histName_dataCat1 = samples_mumu_UL2018_DC[0] + "_Cat1_"+varName;
    TString histName_dataCat2 = samples_mumu_UL2018_DC[0] + "_Cat2_"+varName;
    TString histName_dataCat3 = samples_mumu_UL2018_DC[0] + "_Cat3_"+varName;
    TString histName_dataCat4 = samples_mumu_UL2018_DC[0] + "_Cat4_"+varName;
    TString histName_dataCat12 = samples_mumu_UL2018_DC[0] + "_Cat12_"+varName;
    TString histName_dataCat13 = samples_mumu_UL2018_DC[0] + "_Cat13_"+varName;
    TString histName_dataCat14 = samples_mumu_UL2018_DC[0] + "_Cat14_"+varName;
    TString histName_dataCat23 = samples_mumu_UL2018_DC[0] + "_Cat23_"+varName;
    TString histName_dataCat24 = samples_mumu_UL2018_DC[0] + "_Cat24_"+varName;
    TString histName_dataCat34 = samples_mumu_UL2018_DC[0] + "_Cat34_"+varName;
    TString histName_dataCat123 = samples_mumu_UL2018_DC[0] + "_Cat123_"+varName;
    TString histName_dataCat124 = samples_mumu_UL2018_DC[0] + "_Cat124_"+varName;
    TString histName_dataCat234 = samples_mumu_UL2018_DC[0] + "_Cat234_"+varName;
    TString histName_dataCat1234 = samples_mumu_UL2018_DC[0] + "_Cat1234_"+varName;
    
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
    
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat1 = " << hist_DR[0]->GetEntries() << " : " << hist_DR[0]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat2 = " << hist_DR[1]->GetEntries() << " : " << hist_DR[1]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat3 = " << hist_DR[2]->GetEntries() << " : " << hist_DR[2]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat4 = " << hist_DR[3]->GetEntries() << " : " << hist_DR[3]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat12 = " << hist_DR[4]->GetEntries() << " : " << hist_DR[4]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat13 = " << hist_DR[5]->GetEntries() << " : " << hist_DR[5]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat14 = " << hist_DR[6]->GetEntries() << " : " << hist_DR[6]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat23 = " << hist_DR[7]->GetEntries() << " : " << hist_DR[7]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat24 = " << hist_DR[8]->GetEntries() << " : " << hist_DR[8]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat34 = " << hist_DR[9]->GetEntries() << " : " << hist_DR[9]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat123 = " << hist_DR[10]->GetEntries() << " : " << hist_DR[10]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat124 = " << hist_DR[11]->GetEntries() << " : " << hist_DR[11]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat234 = " << hist_DR[12]->GetEntries() << " : " << hist_DR[12]->GetSumOfWeights() << std::endl;
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven Cat1234 = " << hist_DR[13]->GetEntries() << " : " << hist_DR[13]->GetSumOfWeights() << std::endl;

    float reducibleYield_mumu = hist_DR[0]->GetSumOfWeights() + hist_DR[1]->GetSumOfWeights() + hist_DR[2]->GetSumOfWeights() + hist_DR[3]->GetSumOfWeights() - hist_DR[4]->GetSumOfWeights() - hist_DR[5]->GetSumOfWeights() - hist_DR[6]->GetSumOfWeights() - hist_DR[7]->GetSumOfWeights() - hist_DR[8]->GetSumOfWeights() - hist_DR[9]->GetSumOfWeights() + hist_DR[10]->GetSumOfWeights() + hist_DR[11]->GetSumOfWeights() + hist_DR[12]->GetSumOfWeights() - hist_DR[13]->GetSumOfWeights();
    
    float reducibleYield_mumu_old = hist_DR[2]->GetSumOfWeights() + hist_DR[3]->GetSumOfWeights() - hist_DR[9]->GetSumOfWeights();
    hist_DR[13]->Add(hist_DR[13],hist_DR[12],-1,1);

    std::cout << " Reducible Yields(mumu) = " << reducibleYield_mumu << std::endl;

    TFile * file_data_ee = new TFile("./NTuples_ee"+finalState+"/"+samples_ee_UL2018_DC[0]+".root");
    TTree * tree_data_ee = (TTree*)file_data_ee->Get("SynTree");
    
    TString histName_dataCat1_ee = samples_ee_UL2018_DC[0] + "_Cat1_"+varName;
    TString histName_dataCat2_ee = samples_ee_UL2018_DC[0] + "_Cat2_"+varName;
    TString histName_dataCat3_ee = samples_ee_UL2018_DC[0] + "_Cat3_"+varName;
    TString histName_dataCat4_ee = samples_ee_UL2018_DC[0] + "_Cat4_"+varName;
    TString histName_dataCat12_ee = samples_ee_UL2018_DC[0] + "_Cat12_"+varName;
    TString histName_dataCat13_ee = samples_ee_UL2018_DC[0] + "_Cat13_"+varName;
    TString histName_dataCat14_ee = samples_ee_UL2018_DC[0] + "_Cat14_"+varName;
    TString histName_dataCat23_ee = samples_ee_UL2018_DC[0] + "_Cat23_"+varName;
    TString histName_dataCat24_ee = samples_ee_UL2018_DC[0] + "_Cat24_"+varName;
    TString histName_dataCat34_ee = samples_ee_UL2018_DC[0] + "_Cat34_"+varName;
    TString histName_dataCat123_ee = samples_ee_UL2018_DC[0] + "_Cat123_"+varName;
    TString histName_dataCat124_ee = samples_ee_UL2018_DC[0] + "_Cat124_"+varName;
    TString histName_dataCat234_ee = samples_ee_UL2018_DC[0] + "_Cat234_"+varName;
    TString histName_dataCat1234_ee = samples_ee_UL2018_DC[0] + "_Cat1234_"+varName;
    
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
    
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat1 = " << hist_DR_ee[0]->GetEntries() << " : " << hist_DR_ee[0]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat2 = " << hist_DR_ee[1]->GetEntries() << " : " << hist_DR_ee[1]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat3 = " << hist_DR_ee[2]->GetEntries() << " : " << hist_DR_ee[2]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat4 = " << hist_DR_ee[3]->GetEntries() << " : " << hist_DR_ee[3]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat12 = " << hist_DR_ee[4]->GetEntries() << " : " << hist_DR_ee[4]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat13 = " << hist_DR_ee[5]->GetEntries() << " : " << hist_DR_ee[5]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat14 = " << hist_DR_ee[6]->GetEntries() << " : " << hist_DR_ee[6]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat23 = " << hist_DR_ee[7]->GetEntries() << " : " << hist_DR_ee[7]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat24 = " << hist_DR_ee[8]->GetEntries() << " : " << hist_DR_ee[8]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat34 = " << hist_DR_ee[9]->GetEntries() << " : " << hist_DR_ee[9]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat123 = " << hist_DR_ee[10]->GetEntries() << " : " << hist_DR_ee[10]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat124 = " << hist_DR_ee[11]->GetEntries() << " : " << hist_DR_ee[11]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat234 = " << hist_DR_ee[12]->GetEntries() << " : " << hist_DR_ee[12]->GetSumOfWeights() << std::endl;
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven Cat1234 = " << hist_DR_ee[13]->GetEntries() << " : " << hist_DR_ee[13]->GetSumOfWeights() << std::endl;

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
    TString histName_data_SS = samples_mumu_UL2018_DC[0] + "_SS_"+varName;
    hist_SS = new TH1D(histName_data_SS,"",nBins,bins);
    hist_SS->Sumw2();
    tree_data->Draw(varName+">>"+histName_data_SS, cutsRelaxed + HiggsSSCuts + extraCuts);
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven SS mumu (before Nor.) = " << hist_SS->GetEntries() << " : " << hist_SS->GetSumOfWeights() << std::endl;
    hist_SS->Scale(reducibleYield_mumu/hist_SS->GetSumOfWeights());

    TString histName_data_SS_ee = samples_ee_UL2018_DC[0] + "_SS_"+varName;
    hist_SS_ee = new TH1D(histName_data_SS_ee,"",nBins,bins);
    hist_SS_ee->Sumw2();
    tree_data_ee->Draw(varName+">>"+histName_data_SS_ee, cutsRelaxed + HiggsSSCuts + extraCuts);
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven SS ee (before Nor.) = " << hist_SS_ee->GetEntries() << " : " << hist_SS_ee->GetSumOfWeights() << std::endl;
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
    
    TString histName_data_FR = samples_mumu_UL2018_DC[0] + "_FR_"+varName;
    hist_FR = new TH1D(histName_data_FR,"",nBins,bins);
    hist_FR->Sumw2();
    tree_data->Draw(varName+">>"+histName_data_FR, cuts_FR_mumu + HiggsSRCuts);
    std::cout << samples_mumu_UL2018_DC[0] << " -> Data Driven FR mumu (before Nor.) = " << hist_FR->GetEntries() << " : " << hist_FR->GetSumOfWeights() << std::endl;
    hist_FR->Scale(reducibleYield_mumu/hist_FR->GetSumOfWeights());

    TString histName_data_FR_ee = samples_ee_UL2018_DC[0] + "_FR_"+varName;
    hist_FR_ee = new TH1D(histName_data_FR_ee,"",nBins,bins);
    hist_FR_ee->Sumw2();
    tree_data_ee->Draw(varName+">>"+histName_data_FR_ee, cuts_FR_ee + HiggsSRCuts);
    std::cout << samples_ee_UL2018_DC[0] << " -> Data Driven FR ee (before Nor.) = " << hist_FR_ee->GetEntries() << " : " << hist_FR_ee->GetSumOfWeights() << std::endl;
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


        TFile * file = new TFile("./Combineworkspace/shapes/produceArea/AZh_"+era+"_"+ZfinalStates+finalState+"_"+varName+"_"+Measurement+Shift+"_"+CatName+".root","recreate");
    
        file->mkdir(ZfinalStates+finalState);
        file->cd(ZfinalStates+finalState);

        TH1D * Reducible = (TH1D*)hist_SS->Clone("Reducible_"+UncerSuffix+Shift);
        if(UseFRShape)
            Reducible = (TH1D*)hist_FR->Clone("Reducible_"+UncerSuffix+Shift);

        float totReducible = Reducible->GetSumOfWeights();
        float totBkg = totReducible;
        std::cout << "Reducible: " << Reducible->GetSumOfWeights() << std::endl;

        file->Write();
        file->Close();
    }

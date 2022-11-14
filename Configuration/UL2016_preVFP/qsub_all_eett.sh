#!/bin/sh

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_Data.conf NTuples_SingleElectron_Run2016B-ver1 50
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_Data.conf NTuples_SingleElectron_Run2016B-ver2 50
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_Data.conf NTuples_SingleElectron_Run2016C 50
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_Data.conf NTuples_SingleElectron_Run2016D 50
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_Data.conf NTuples_SingleElectron_Run2016E 50
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_Data.conf NTuples_SingleElectron_Run2016F_preVFP 50

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_DY.conf NTuples_DY1JetsToLL_M50 200
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_DY.conf NTuples_DY2JetsToLL_M50 200
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_DY.conf NTuples_DY3JetsToLL_M50 200
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_DY.conf NTuples_DY4JetsToLL_M50 200
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_DY.conf NTuples_DYJetsToLL_M50 200

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToContinToZZTo2e2mu 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToContinToZZTo2e2tau 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToContinToZZTo2mu2tau 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToContinToZZTo4tau 10

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_TTWJetsToLNu 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_TTTo2L2Nu 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_TTToHadronic 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_TTToSemiLeptonic 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8 40

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_WWW 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_WWW-ext 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_WWZ 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_WZZ 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_WZZ_ext 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_ZZZ 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_ZZZ-ext 40

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8 40

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_ZZTo4L_TuneCP5_13TeV_powheg_pythia8 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8 40

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_VBFHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_WminusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_WplusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_ZHToTauTau_M125_CP5_13TeV-powheg-pythia8_ext1 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluHToWWTo2L2Nu 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_VBFHToWWTo2L2Nu 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_HWminusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_HWplusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_HZJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8 40
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluZH_HToWW_ZTo2L_M125_13TeV_powheg_pythia8_TuneCP5_PSweights 40


./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M225_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M250_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M275_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M300_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M325_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M350_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M375_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M400_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M450_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M500_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M600_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M700_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M750_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M800_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M900_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M1000_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M1200_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M1400_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M1600_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M1800_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_GluGluToAToZhToLLTauTau_M2000_TuneCP5_13TeV-madgraph-pythia8 10

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M225_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M250_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M275_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M300_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M325_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M350_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M375_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M400_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M450_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M500_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M600_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M700_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M750_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M800_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M900_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M1000_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M1200_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M1400_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M1600_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M1800_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eett_MC.conf NTuples_BBAToZhToLLTauTau_M2000_TuneCP5_13TeV-madgraph-pythia8 10

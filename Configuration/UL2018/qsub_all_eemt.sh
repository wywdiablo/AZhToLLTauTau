#!/bin/sh

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_Data.conf NTuples_EGamma-Run2018A-UL2018 50
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_Data.conf NTuples_EGamma-Run2018B-UL2018 50
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_Data.conf NTuples_EGamma-Run2018C-UL2018 50
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_Data.conf NTuples_EGamma-Run2018D-UL2018 50

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_DY.conf NTuples_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 300
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_DY.conf NTuples_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext 300

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_DY.conf NTuples_DY1JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8 300
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_DY.conf NTuples_DY2JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8 300
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_DY.conf NTuples_DY3JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8 300
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_DY.conf NTuples_DY4JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8 300

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8 15
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8 15
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8 15
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8 15
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8 15
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8 15

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8 150
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_TTToHadronic_TuneCP5_13TeV-powheg-pythia8 150
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8 150
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8 60

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_ext1 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_WZZ_TuneCP5_13TeV-amcatnlo-pythia8 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_ZZZ_TuneCP5_13TeV-amcatnlo-pythia8 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_ZZTo4L_TuneCP5_13TeV_powheg_pythia8 60
#./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_WZTo3LNu_TuneCP5_13TeV-powheg-pythia8 100

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluHToTauTau 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_VBFHToTauTau 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_WminusHToTauTau 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_WplusHToTauTau 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_ZHToTauTau_M125_CP5_13TeV-powheg-pythia8_ext1 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_ttHToTauTau 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluHToWWTo2L2Nu 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_VBFHToWWTo2L2Nu 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_HWminusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_HWplusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_HZJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluZH_HToWW_ZTo2L_M125_13TeV_powheg_pythia8_TuneCP5_PSweights 60
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8 60

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M225_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M250_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M275_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M300_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M325_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M350_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M375_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M400_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M450_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M500_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M600_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M700_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M750_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M800_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M900_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M1000_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M1200_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M1400_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M1600_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M1800_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_GluGluToAToZhToLLTauTau_M2000_TuneCP5_13TeV-madgraph-pythia8 10

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M225_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M250_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M275_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M300_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M325_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M350_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M375_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M400_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M450_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M500_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M600_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M700_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M750_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M800_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M900_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M1000_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M1200_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M1400_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M1600_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M1800_TuneCP5_13TeV-madgraph-pythia8 10
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eemt_MC.conf NTuples_BBAToZhToLLTauTau_M2000_TuneCP5_13TeV-madgraph-pythia8 10

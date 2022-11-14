#!/bin/sh

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_Data.conf NTuples_SingleElectron_Run2016F 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_Data.conf NTuples_SingleElectron_Run2016G 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_Data.conf NTuples_SingleElectron_Run2016H 100

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_DY.conf NTuples_DY1JetsToLL_M-50 200
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_DY.conf NTuples_DY2JetsToLL_M-50 200
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_DY.conf NTuples_DY3JetsToLL_M-50 200
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_DY.conf NTuples_DY4JetsToLL_M-50 200
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_DY.conf NTuples_DYJetsToLL_M-50 200

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_GluGluToContinToZZTo2e2mu 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_GluGluToContinToZZTo2e2tau 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_GluGluToContinToZZTo2mu2tau 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_GluGluToContinToZZTo4e 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_GluGluToContinToZZTo4mu 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_GluGluToContinToZZTo4tau 100

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_TTWJetsToLNu 200
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_TTTo2L2Nu 200
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_TTToHadronic 200
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_TTToSemiLeptonic 200
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8 200

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_WWW_4F 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_WWW_4F_ext1 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_WWZ_4F 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_WWZ_4F_ext1 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_WZZ 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_WZZ_ext1 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_ZZZ 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_ZZZ_ext1 100

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_WZTo3LNu 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8 100

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_ZZTo4L_TuneCP5_13TeV_powheg_pythia8 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8 100

./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_GluGluHToTauTau_M125 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_VBFHToTauTau_M125 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_WminusHToTauTau_M125 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_WplusHToTauTau_M125 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_ZHToTauTau_M125 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_ttHToTauTau_M125 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_GluGluHToWWTo2L2Nu 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_VBFHToWWTo2L2Nu 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_HWminusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_HWplusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_HZJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8 100
./HTC_submit_seq.sh SynchNtupleProducer_AZh_ee analysisMacro_AZh_eee_MC.conf NTuples_GluGluZH_HToWW_ZTo2L_M125_13TeV_powheg_pythia8_TuneCP5_PSweights 100

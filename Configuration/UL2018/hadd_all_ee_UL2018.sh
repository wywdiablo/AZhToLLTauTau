#!/bin/sh
./hadd.sh NTuples_EGamma-Run2018A-UL2018
./hadd.sh NTuples_EGamma-Run2018B-UL2018
./hadd.sh NTuples_EGamma-Run2018C-UL2018
./hadd.sh NTuples_EGamma-Run2018D-UL2018

rm NTuples_EGamma_UL2018.root
hadd NTuples_EGamma_UL2018.root NTuples_EGamma-Run2018A-UL2018.root NTuples_EGamma-Run2018B-UL2018.root NTuples_EGamma-Run2018C-UL2018.root NTuples_EGamma-Run2018D-UL2018.root

./hadd.sh NTuples_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
./hadd.sh NTuples_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext

rm NTuples_DYJetsToLL_M-50.root
hadd NTuples_DYJetsToLL_M-50.root NTuples_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root NTuples_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext.root

./hadd.sh NTuples_DY1JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8
./hadd.sh NTuples_DY2JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8
./hadd.sh NTuples_DY3JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8
./hadd.sh NTuples_DY4JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8

./hadd.sh NTuples_GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8
./hadd.sh NTuples_GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8
./hadd.sh NTuples_GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8
./hadd.sh NTuples_GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8
./hadd.sh NTuples_GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8
./hadd.sh NTuples_GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8

./hadd.sh NTuples_TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8
./hadd.sh NTuples_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8
./hadd.sh NTuples_TTToHadronic_TuneCP5_13TeV-powheg-pythia8
./hadd.sh NTuples_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8
./hadd.sh NTuples_ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8

./hadd.sh NTuples_WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8
./hadd.sh NTuples_WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_ext1

rm NTuples_WWW.root
hadd NTuples_WWW.root NTuples_WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8.root NTuples_WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_ext1.root

./hadd.sh NTuples_WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8
./hadd.sh NTuples_WZZ_TuneCP5_13TeV-amcatnlo-pythia8
./hadd.sh NTuples_ZZZ_TuneCP5_13TeV-amcatnlo-pythia8

./hadd.sh NTuples_WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8
./hadd.sh NTuples_ZZTo4L_TuneCP5_13TeV_powheg_pythia8

./hadd.sh NTuples_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8
./hadd.sh NTuples_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8

./hadd.sh NTuples_GluGluHToTauTau
./hadd.sh NTuples_VBFHToTauTau
./hadd.sh NTuples_WminusHToTauTau
./hadd.sh NTuples_WplusHToTauTau
./hadd.sh NTuples_ZHToTauTau_M125_CP5_13TeV-powheg-pythia8_ext1
./hadd.sh NTuples_ttHToTauTau
./hadd.sh NTuples_GluGluHToWWTo2L2Nu
./hadd.sh NTuples_VBFHToWWTo2L2Nu
./hadd.sh NTuples_HWplusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8
./hadd.sh NTuples_HWminusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8
./hadd.sh NTuples_HZJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8
./hadd.sh NTuples_GluGluZH_HToWW_ZTo2L_M125_13TeV_powheg_pythia8_TuneCP5_PSweights
./hadd.sh NTuples_GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8

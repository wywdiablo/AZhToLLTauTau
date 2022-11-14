#!/bin/sh
./hadd.sh NTuples_SingleElectron_Run2016B-ver1
./hadd.sh NTuples_SingleElectron_Run2016B-ver2
./hadd.sh NTuples_SingleElectron_Run2016C
./hadd.sh NTuples_SingleElectron_Run2016D
./hadd.sh NTuples_SingleElectron_Run2016E
./hadd.sh NTuples_SingleElectron_Run2016F_preVFP

rm NTuples_SingleElectron_UL2016preVFP.root
hadd NTuples_SingleElectron_UL2016preVFP.root NTuples_SingleElectron_Run2016B-ver1.root NTuples_SingleElectron_Run2016B-ver2.root NTuples_SingleElectron_Run2016C.root NTuples_SingleElectron_Run2016D.root NTuples_SingleElectron_Run2016E.root NTuples_SingleElectron_Run2016F_preVFP.root

./hadd.sh NTuples_DYJetsToLL_M50
./hadd.sh NTuples_DY1JetsToLL_M50
./hadd.sh NTuples_DY2JetsToLL_M50
./hadd.sh NTuples_DY3JetsToLL_M50
./hadd.sh NTuples_DY4JetsToLL_M50

./hadd.sh NTuples_GluGluToContinToZZTo2e2mu
./hadd.sh NTuples_GluGluToContinToZZTo2e2tau
./hadd.sh NTuples_GluGluToContinToZZTo2mu2tau
./hadd.sh NTuples_GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8
./hadd.sh NTuples_GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8
./hadd.sh NTuples_GluGluToContinToZZTo4tau

./hadd.sh NTuples_TTWJetsToLNu
./hadd.sh NTuples_TTTo2L2Nu
./hadd.sh NTuples_TTToHadronic
./hadd.sh NTuples_TTToSemiLeptonic
./hadd.sh NTuples_ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8

./hadd.sh NTuples_WWW
./hadd.sh NTuples_WWW-ext
rm NTuples_WWW_all.root
hadd NTuples_WWW_all.root NTuples_WWW.root NTuples_WWW-ext.root

./hadd.sh NTuples_WWZ

./hadd.sh NTuples_WZZ
./hadd.sh NTuples_WZZ_ext
rm NTuples_WZZ_all.root
hadd NTuples_WZZ_all.root NTuples_WZZ.root NTuples_WZZ_ext.root

./hadd.sh NTuples_ZZZ
./hadd.sh NTuples_ZZZ-ext
rm NTuples_ZZZ_all.root
hadd NTuples_ZZZ_all.root NTuples_ZZZ.root NTuples_ZZZ-ext.root

./hadd.sh NTuples_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8
./hadd.sh NTuples_WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8

./hadd.sh NTuples_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8
./hadd.sh NTuples_ZZTo4L_TuneCP5_13TeV_powheg_pythia8

./hadd.sh NTuples_GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8
./hadd.sh NTuples_VBFHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8
./hadd.sh NTuples_WminusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8
./hadd.sh NTuples_WplusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8
./hadd.sh NTuples_ZHToTauTau_M125_CP5_13TeV-powheg-pythia8_ext1
./hadd.sh NTuples_ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8
./hadd.sh NTuples_GluGluHToWWTo2L2Nu
./hadd.sh NTuples_VBFHToWWTo2L2Nu
./hadd.sh NTuples_HWminusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8
./hadd.sh NTuples_HWplusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8
./hadd.sh NTuples_HZJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8
./hadd.sh NTuples_GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8
./hadd.sh NTuples_GluGluZH_HToWW_ZTo2L_M125_13TeV_powheg_pythia8_TuneCP5_PSweights

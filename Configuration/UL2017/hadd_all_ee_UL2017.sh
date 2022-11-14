#!/bin/sh
./hadd.sh NTuples_SingleEleRun2017B
./hadd.sh NTuples_SingleElectron_2017C
./hadd.sh NTuples_SingleElectron_2017D
./hadd.sh NTuples_SingleElectron_2017E
./hadd.sh NTuples_SingleElectron_2017F

rm NTuples_SingleElectron_UL2017.root
hadd NTuples_SingleElectron_UL2017.root NTuples_SingleEleRun2017B.root NTuples_SingleElectron_2017C.root NTuples_SingleElectron_2017D.root NTuples_SingleElectron_2017E.root NTuples_SingleElectron_2017F.root

./hadd.sh NTuples_DYJetsToLL_M50
./hadd.sh NTuples_DYJetsToLL_M50_ext1

rm NTuples_DYJetsToLL_M50_all.root
hadd NTuples_DYJetsToLL_M50_all.root NTuples_DYJetsToLL_M50.root NTuples_DYJetsToLL_M50_ext1.root

./hadd.sh NTuples_DY1JetsToLL_M50
./hadd.sh NTuples_DY2JetsToLL_M50
./hadd.sh NTuples_DY3JetsToLL_M50
./hadd.sh NTuples_DY4JetsToLL_M50

./hadd.sh NTuples_GluGluToContinToZZTo2e2mu
./hadd.sh NTuples_GluGluToContinToZZTo2e2tau
./hadd.sh NTuples_GluGluToContinToZZTo2mu2tau
./hadd.sh NTuples_GluGluToContinToZZTo4e
./hadd.sh NTuples_GluGluToContinToZZTo4mu
./hadd.sh NTuples_GluGluToContinToZZTo4tau

./hadd.sh NTuples_TTWJetsToLNu
./hadd.sh NTuples_TTTo2L2Nu
./hadd.sh NTuples_TTToHadronic
./hadd.sh NTuples_TTToSemiLeptonic
./hadd.sh NTuples_TTZJets

./hadd.sh NTuples_WWW_4F
./hadd.sh NTuples_WWW_4F_ext1

rm NTuples_WWW.root
hadd NTuples_WWW.root NTuples_WWW_4F.root NTuples_WWW_4F_ext1.root

./hadd.sh NTuples_WWZ
./hadd.sh NTuples_WZZ
./hadd.sh NTuples_WZZ_ext

rm NTuples_WZZ_all.root
hadd NTuples_WZZ_all.root NTuples_WZZ.root NTuples_WZZ_ext.root

./hadd.sh NTuples_ZZZ
./hadd.sh NTuples_WZTo3LNu
./hadd.sh NTuples_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8
./hadd.sh NTuples_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8

./hadd.sh NTuples_ZZTo4L

./hadd.sh NTuples_GluGluHToTauTau
./hadd.sh NTuples_VBFHToTauTau
./hadd.sh NTuples_WminusHToTauTau
./hadd.sh NTuples_WplusHToTauTau
./hadd.sh NTuples_ZHToTauTau_M125_CP5_13TeV-powheg-pythia8
./hadd.sh NTuples_ZHToTauTau_M125_CP5_13TeV-powheg-pythia8_ext1

rm NTuples_ZHToTauTau.root
hadd NTuples_ZHToTauTau.root NTuples_ZHToTauTau_M125_CP5_13TeV-powheg-pythia8.root NTuples_ZHToTauTau_M125_CP5_13TeV-powheg-pythia8_ext1.root

./hadd.sh NTuples_ttHToTauTau
./hadd.sh NTuples_GluGluHToWWTo2L2Nu
./hadd.sh NTuples_VBFHToWWTo2L2Nu
./hadd.sh NTuples_HWminusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8
./hadd.sh NTuples_HWplusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8
./hadd.sh NTuples_HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8
./hadd.sh NTuples_GluGluZH_HToWW_ZTo2L_M125_13TeV_powheg_pythia8_TuneCP5_PSweights
./hadd.sh NTuples_GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8

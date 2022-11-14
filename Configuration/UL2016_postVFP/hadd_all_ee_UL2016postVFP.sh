#!/bin/sh
./hadd.sh NTuples_SingleElectron_Run2016F
./hadd.sh NTuples_SingleElectron_Run2016G
./hadd.sh NTuples_SingleElectron_Run2016H

rm NTuples_SingleElectron_UL2016postVFP.root
hadd NTuples_SingleElectron_UL2016postVFP.root NTuples_SingleElectron_Run2016F.root NTuples_SingleElectron_Run2016G.root NTuples_SingleElectron_Run2016H.root

./hadd.sh NTuples_DYJetsToLL_M-50
./hadd.sh NTuples_DY1JetsToLL_M-50
./hadd.sh NTuples_DY2JetsToLL_M-50
./hadd.sh NTuples_DY3JetsToLL_M-50
./hadd.sh NTuples_DY4JetsToLL_M-50

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
./hadd.sh NTuples_ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8

./hadd.sh NTuples_WWW_4F
./hadd.sh NTuples_WWW_4F_ext1
rm NTuples_WWW.root
hadd NTuples_WWW.root NTuples_WWW_4F.root NTuples_WWW_4F_ext1.root

./hadd.sh NTuples_WWZ_4F
./hadd.sh NTuples_WWZ_4F_ext1
rm NTuples_WWZ.root
hadd NTuples_WWZ.root NTuples_WWZ_4F.root NTuples_WWZ_4F_ext1.root

./hadd.sh NTuples_WZZ
./hadd.sh NTuples_WZZ_ext1
rm NTuples_WZZ_all.root
hadd NTuples_WZZ_all.root NTuples_WZZ.root NTuples_WZZ_ext1.root

./hadd.sh NTuples_ZZZ
./hadd.sh NTuples_ZZZ_ext1
rm NTuples_ZZZ_all.root
hadd NTuples_ZZZ_all.root NTuples_ZZZ.root NTuples_ZZZ_ext1.root


./hadd.sh NTuples_WZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8
./hadd.sh NTuples_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8
./hadd.sh NTuples_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8_ext1

rm NTuples_WZTo3LNu.root
hadd NTuples_WZTo3LNu.root NTuples_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root NTuples_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8_ext1.root

./hadd.sh NTuples_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8
./hadd.sh NTuples_WZTo3LNu
./hadd.sh NTuples_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8
./hadd.sh NTuples_ZZTo4L_TuneCP5_13TeV_powheg_pythia8

./hadd.sh NTuples_GluGluHToTauTau_M125
./hadd.sh NTuples_VBFHToTauTau_M125
./hadd.sh NTuples_WminusHToTauTau_M125
./hadd.sh NTuples_WplusHToTauTau_M125
./hadd.sh NTuples_ZHToTauTau_M125
./hadd.sh NTuples_ttHToTauTau_M125
./hadd.sh NTuples_GluGluHToWWTo2L2Nu
./hadd.sh NTuples_VBFHToWWTo2L2Nu
./hadd.sh NTuples_HWminusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8
./hadd.sh NTuples_HWplusJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8
./hadd.sh NTuples_HZJ_HToWW_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8
./hadd.sh NTuples_GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8
./hadd.sh NTuples_GluGluZH_HToWW_ZTo2L_M125_13TeV_powheg_pythia8_TuneCP5_PSweights

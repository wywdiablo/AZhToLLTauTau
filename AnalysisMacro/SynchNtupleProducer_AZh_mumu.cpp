#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LepTauFakeRate.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauPOG/TauIDSFs/interface/TauIDSFTool.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TMath.h"
#include "AZhTools/FakeRateTools/interface/FakeRateWeightTools.h"
#include "DesyTauAnalyses/NTupleMaker/interface/btagSF.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AZh_functions.h"

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

    const string era = cfg.get<string>("Era");
    const bool isData = cfg.get<bool>("IsData");
    const bool isDY = cfg.get<bool>("IsDY");
    const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection");

    const string jsonFile = cfg.get<string>("jsonFile");
    const bool applyRecoilCorrections = cfg.get<bool>("ApplyRecoilCorrections");
    const string recoilFileName   = cfg.get<string>("RecoilFileName");
    const string dataPUFile = cfg.get<string>("DataPUFile");
    const string mcPUFile = cfg.get<string>("MCPUFile");
    const string sampleName = cfg.get<string>("sampleName");
    // pile up reweighting
    const bool applyPUreweighting = cfg.get<bool>("ApplyPUreweighting");
    const bool ApplySVFit       = cfg.get<bool>("ApplySVFit");
    const bool ApplyFastMTT     = cfg.get<bool>("ApplyFastMTT");
    const bool ApplyPuppiMet     = cfg.get<bool>("ApplyPuppiMet");
    const bool applySystShift    = cfg.get<bool>("ApplySystShift");

    // which Higgs final states: TT, MT, ET, EM, MM, EE
    const string HiggsFinalState = cfg.get<string>("HiggsFinalState");
    
    // kinematic cuts on electrons
    const float ptEleCut = cfg.get<float>("ptEleCut");
    const float etaEleCut = cfg.get<float>("etaEleCut");
    const float dxyEleCut = cfg.get<float>("dxyEleCut");
    const float dzEleCut = cfg.get<float>("dzEleCut");

    // kinematic cuts on muons
    const float ptMuonCut  = cfg.get<float>("ptMuonCut");
    const float etaMuonCut = cfg.get<float>("etaMuonCut");
    const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
    const float dzMuonCut      = cfg.get<float>("dzMuonCut");

    //kinematic cuts on taus
    const float ptTauCut = cfg.get<float>("ptTauCut");
    const float etaTauCut = cfg.get<float>("etaTauCut");
    const float dzTauCut = cfg.get<float>("dzTauCut");
    
    // jets
    const float jetPtLowCut = cfg.get<float>("JetPtLowCut");
    const float jetPtHighCut = cfg.get<float>("JetPtHighCut");
    const float jetEtaCut = cfg.get<float>("JetEtaCut");
    const float dRJetLeptonCut = cfg.get<float>("dRJetLeptonCut");
    const float btagCut = cfg.get<float>("btagCut");
    const string bTagAlgorithm = cfg.get<string>("BTagAlgorithm");
    const string bTagFile = cfg.get<string>("BTagFile");
    const string bTagEffFile = cfg.get<string>("BTagEffFile");
    const string bTagDiscriminator1 = cfg.get<string>("BTagDiscriminator1");
    const string bTagDiscriminator2 = cfg.get<string>("BTagDiscriminator2");
    const string bTagDiscriminator3 = cfg.get<string>("BTagDiscriminator3");
    const float bJetEtaCut = cfg.get<float>("bJetEtaCut");
    
    // topological cuts
    const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
    const float dPhileptonsCut = cfg.get<float>("dPhileptonsCut");
    const float DRTrigMatch    = cfg.get<float>("DRTrigMatch");

    // trigger
    const bool applyTrigger = cfg.get<bool>("ApplyTrigger");
    const bool applySingleMuTriggerOnly = cfg.get<bool>("ApplySingleMuTriggerOnly");
    const string muonTriggerName_1  = cfg.get<string>("MuonTriggerName_1");
    const string muonTriggerName_2  = cfg.get<string>("MuonTriggerName_2");
    const string muonTriggerName_3  = cfg.get<string>("MuonTriggerName_3");
    const string muonTriggerName_4  = cfg.get<string>("MuonTriggerName_4");

    const string singleMuonFilterName_1_a = cfg.get<string>("SingleMuonFilterName_1_a");
    const string singleMuonFilterName_1_b = cfg.get<string>("SingleMuonFilterName_1_b");
    const float singleMuonTriggerPtCut_1 = cfg.get<float>("SingleMuonTriggerPtCut_1");
    
    const string singleMuonFilterName_2 = cfg.get<string>("SingleMuonFilterName_2");
    const float singleMuonTriggerPtCut_2 = cfg.get<float>("SingleMuonTriggerPtCut_2");
    
    const string doubleMuonLeg1FilterName_3 = cfg.get<string>("DoubleMuonLeg1FilterName_3");
    const float doubleMuonLeg1TriggerPtCut_3 = cfg.get<float>("DoubleMuonLeg1TriggerPtCut_3");
    const string doubleMuonLeg2FilterName_3_a = cfg.get<string>("DoubleMuonLeg2FilterName_3_a");
    const string doubleMuonLeg2FilterName_3_b = cfg.get<string>("DoubleMuonLeg2FilterName_3_b");
    const float doubleMuonLeg2TriggerPtCut_3 = cfg.get<float>("DoubleMuonLeg2TriggerPtCut_3");
    const string doubleMuonDzFilterName_3 = cfg.get<string>("DoubleMuonDzFilterName_3");
    const string doubleMuonMassFilterName_3 = cfg.get<string>("DoubleMuonMassFilterName_3");

    const string doubleMuonLeg1FilterName_4 = cfg.get<string>("DoubleMuonLeg1FilterName_4");
    const float doubleMuonLeg1TriggerPtCut_4 = cfg.get<float>("DoubleMuonLeg1TriggerPtCut_4");
    const string doubleMuonLeg2FilterName_4 = cfg.get<string>("DoubleMuonLeg2FilterName_4");
    const float doubleMuonLeg2TriggerPtCut_4 = cfg.get<float>("DoubleMuonLeg2TriggerPtCut_4");
    const string doubleMuonDzFilterName_4 = cfg.get<string>("DoubleMuonDzFilterName_4");

    TString MuonTriggerName_1(muonTriggerName_1);
    TString MuonTriggerName_2(muonTriggerName_2);
    TString MuonTriggerName_3(muonTriggerName_3);
    TString MuonTriggerName_4(muonTriggerName_4);

    TString SingleMuonFilterName_1_a(singleMuonFilterName_1_a);
    TString SingleMuonFilterName_1_b(singleMuonFilterName_1_b);
    TString SingleMuonFilterName_2(singleMuonFilterName_2);
    TString DoubleMuonLeg1FilterName_3(doubleMuonLeg1FilterName_3);
    TString DoubleMuonLeg2FilterName_3_a(doubleMuonLeg2FilterName_3_a);
    TString DoubleMuonLeg2FilterName_3_b(doubleMuonLeg2FilterName_3_b);
    TString DoubleMuonDzFilterName_3(doubleMuonDzFilterName_3);
    TString DoubleMuonMassFilterName_3(doubleMuonMassFilterName_3);

    TString DoubleMuonLeg1FilterName_4(doubleMuonLeg1FilterName_4);
    TString DoubleMuonLeg2FilterName_4(doubleMuonLeg2FilterName_4);
    TString DoubleMuonDzFilterName_4(doubleMuonDzFilterName_4);
    
    TString BTagDiscriminator1(bTagDiscriminator1);
    TString BTagDiscriminator2(bTagDiscriminator2);
    TString BTagDiscriminator3(bTagDiscriminator3);

    // vertex cuts
    const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");
    const float zVertexCut     = cfg.get<float>("ZVertexCut");
    const float dVertexCut     = cfg.get<float>("DVertexCut");
    
    //scale factor
    const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoFile");
    const string MuonTriggerFile = cfg.get<string>("MuonTriggerFile");
    const string EleIdIsoFile = cfg.get<string>("EleIdIsoFile");
    
    const float MuonScale = cfg.get<float>("MuonScale");
    const float EleScaleBarrel = cfg.get<float>("EleScaleBarrel");
    const float EleScaleEndcap = cfg.get<float>("EleScaleEndcap");
    
    
    //svfit
    const string svFitPtResFile = cfg.get<string>("svFitPtResFile").data();
    const string FakeRateFilesZ = cfg.get<string>("FakeRateFilesZ");
    const string FakeRateFiles1 = cfg.get<string>("FakeRateFiles1");
    const string FakeRateFiles2 = cfg.get<string>("FakeRateFiles2");

    // Z pt weight
    const string ZptweightFile = cfg.get<string>("ZptweightFile");
    
    // **** end of configuration

    string cmsswBase = (getenv ("CMSSW_BASE"));
    string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;

    // Run-lumi selector
    std::vector<Period> periods;
    if (isData)
    { // read the good runs
	  std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
  	  if (inputFileStream.fail() )
      {
            std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
            std::cout << "please check" << std::endl;
            std::cout << "quitting program" << std::endl;
	    exit(-1);
      }
  
        for(std::string s; std::getline(inputFileStream, s); )
        {
           periods.push_back(Period());
           std::stringstream ss(s);
           ss >> periods.back();
        }
    }

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  TH1::SetDefaultSumw2(true);
    
    RecoilCorrector* recoilPFMetCorrector = (RecoilCorrector*) malloc(sizeof(*recoilPFMetCorrector));
    recoilPFMetCorrector = new RecoilCorrector(recoilFileName);

  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
    
    //declare the tree
    TTree * SynTree = new TTree("SynTree","SynTree");
    
    TH1D * histWeightsH = new TH1D("histWeightsH","",1,-0.5,0.5);
    TH1D * histWeightsTriggerH = new TH1D("histWeightsTriggerH","",1,-0.5,0.5);
    TH1D * histWeightsVertexH = new TH1D("histWeightsVertexH","",1,-0.5,0.5);
    TH1D * histWeightsZFoundH = new TH1D("histWeightsZFoundH","",1,-0.5,0.5);
    TH1D * histWeightsHiggsFoundH = new TH1D("histWeightsHiggsFoundH","",1,-0.5,0.5);
    TH1D * histWeightsDRLepsH = new TH1D("histWeightsDRLepsH","",1,-0.5,0.5);
    TH1D * histWeightsExtraVetoH = new TH1D("histWeightsExtraVetoH","",1,-0.5,0.5);
    TH1D * histWeightsOSZH = new TH1D("histWeightsOSZH","",1,-0.5,0.5);
    TH1D * histWeightsTriggerMatchedH = new TH1D("histWeightsTriggerMatchedH","",1,-0.5,0.5);
    TH1D * histWeightsmllH = new TH1D("histWeightsmllH","",1,-0.5,0.5);
    
    Float_t         puweight;
    Float_t         mcweight;
    Float_t         zptweight;

    //Event ID variable
    ULong64_t run;
    UInt_t lumi;
    ULong64_t evt;
    
    //Pile up variable
    Int_t           npu;
    Int_t           npv;
    Float_t         rho;
    
    //Jet related variables
    Int_t   nbtag;
    Int_t   nbtagUp;
    Int_t   nbtagDown;
    Int_t   nbtagHEM2018Up;
    Int_t   nbtagHEM2018Down;
    Int_t   nbtag_mistagUp;
    Int_t   nbtag_mistagDown;
    Int_t   nbtag_btagUp;
    Int_t   nbtag_btagDown;
    Int_t   nbtag_noSF;
    
    Int_t   njets;
    Int_t   njetsUp;
    Int_t   njetsDown;
    Int_t   njetsHEM2018Up;
    Int_t   njetsHEM2018Down;
    Int_t   njetspt20;
    Int_t   njetspt20Up;
    Int_t   njetspt20Down;
    Int_t   njetspt20HEM2018Up;
    Int_t   njetspt20HEM2018Down;
    
    Float_t jpt_1;
    Float_t jeta_1;
    Float_t jphi_1;
    Float_t jpt_1_Up;
    Float_t jeta_1_Up;
    Float_t jphi_1_Up;
    Float_t jpt_1_Down;
    Float_t jeta_1_Down;
    Float_t jphi_1_Down;
    Float_t jpt_2;
    Float_t jeta_2;
    Float_t jphi_2;
    Float_t jpt_2_Up;
    Float_t jeta_2_Up;
    Float_t jphi_2_Up;
    Float_t jpt_2_Down;
    Float_t jeta_2_Down;
    Float_t jphi_2_Down;
    Float_t bpt_1;
    Float_t beta_1;
    Float_t bphi_1;
    Float_t bpt_1_Up;
    Float_t beta_1_Up;
    Float_t bphi_1_Up;
    Float_t bpt_1_Down;
    Float_t beta_1_Down;
    Float_t bphi_1_Down;
    Float_t bhadronFlavor_1;
    Float_t bpt_2;
    Float_t beta_2;
    Float_t bphi_2;
    Float_t bpt_2_Up;
    Float_t beta_2_Up;
    Float_t bphi_2_Up;
    Float_t bpt_2_Down;
    Float_t beta_2_Down;
    Float_t bphi_2_Down;
    Float_t bhadronFlavor_2;
    
    //Leg 1 (highest pt elec/mu from Z candidate)
    Float_t pt_1;
    Float_t phi_1;
    Float_t eta_1;
    Float_t iso_1;
    //Float_t iso_ea_1;
    Int_t gen_match_1;
    Float_t IdRawMva_1;
    Float_t pfmt_1;
    Float_t muonSF_1;
    Float_t IdSF_1;
    Float_t TrigSF_1;
    Float_t electronSF_1;
    Bool_t muonLoose_1;
    Bool_t muonMedium_1;
    Bool_t muonTight_1;

    //Leg 2 (subleading pt elec/mu from Z candidate)
    Float_t pt_2;
    Float_t phi_2;
    Float_t eta_2;
    Float_t iso_2;
    //Float_t iso_ea_2;
    Int_t gen_match_2;
    Float_t IdRawMva_2;
    Float_t pfmt_2;
    Float_t muonSF_2;
    Float_t IdSF_2;
    Float_t TrigSF_2;
    Float_t electronSF_2;
    Bool_t muonLoose_2;
    Bool_t muonMedium_2;
    Bool_t muonTight_2;

    //Leg 3 (Based off of Higgs legs: ET -> E, MT -> M, EM -> E, TT -> highest pt tau)
    Float_t pt_3;
    Float_t phi_3;
    Float_t eta_3;
    Float_t iso_3; //if tau, then it is VSjet raw value
    //Float_t iso_ea_3;
    Int_t gen_match_3;
    Float_t IdRawMva_3;
    Float_t m_3;
    Int_t dm_3;
    Bool_t VSjetVVTight_3;
    Bool_t VSjetVTight_3;
    Bool_t VSjetTight_3;
    Bool_t VSjetMedium_3;
    Bool_t VSjetLoose_3;
    Bool_t VSjetVLoose_3;
    Bool_t VSjetVVLoose_3;
    Bool_t VSjetVVVLoose_3;
    Bool_t VSeVVTight_3;
    Bool_t VSeVTight_3;
    Bool_t VSeTight_3;
    Bool_t VSeMedium_3;
    Bool_t VSeLoose_3;
    Bool_t VSeVLoose_3;
    Bool_t VSeVVLoose_3;
    Bool_t VSeVVVLoose_3;
    Bool_t VSmuTight_3;
    Bool_t VSmuMedium_3;
    Bool_t VSmuLoose_3;
    Bool_t VSmuVLoose_3;
    Float_t pfmt_3;
    Float_t muonSF_3;
    Float_t IdSF_3;
    Float_t electronSF_3;
    Float_t tauSF_3;
    Bool_t electronIDWP90v2_3;
    Bool_t electronIDWP80v2_3;
    Bool_t muonLoose_3;
    Bool_t muonMedium_3;
    Bool_t muonTight_3;
    Float_t tauVSjetSF_3;
    Float_t tauVSjetSFUp_3;
    Float_t tauVSjetSFDown_3;
    Float_t tauVSeSF_3;
    Float_t tauVSeSFUp_3;
    Float_t tauVSeSFDown_3;
    Float_t tauVSmuSF_3;
    Float_t tauVSmuSFUp_3;
    Float_t tauVSmuSFDown_3;
    Float_t tauVSjetVVVLooseSF_3;
    
    //Leg 4 (Based off of Higgs legs: ET -> T, MT -> T, EM -> M, TT -> subleading pt tau)
    Float_t pt_4;
    Float_t phi_4;
    Float_t eta_4;
    Float_t iso_4; //if tau, then it is VSjet raw value
    //Float_t iso_ea_4;
    Int_t gen_match_4;
    Float_t IdRawMva_4;
    Float_t m_4;
    Int_t dm_4;
    Bool_t VSjetVVTight_4;
    Bool_t VSjetVTight_4;
    Bool_t VSjetTight_4;
    Bool_t VSjetMedium_4;
    Bool_t VSjetLoose_4;
    Bool_t VSjetVLoose_4;
    Bool_t VSjetVVLoose_4;
    Bool_t VSjetVVVLoose_4;
    Bool_t VSeVVTight_4;
    Bool_t VSeVTight_4;
    Bool_t VSeTight_4;
    Bool_t VSeMedium_4;
    Bool_t VSeLoose_4;
    Bool_t VSeVLoose_4;
    Bool_t VSeVVLoose_4;
    Bool_t VSeVVVLoose_4;
    Bool_t VSmuTight_4;
    Bool_t VSmuMedium_4;
    Bool_t VSmuLoose_4;
    Bool_t VSmuVLoose_4;
    Float_t pfmt_4;
    Float_t muonSF_4;
    Float_t IdSF_4;
    Float_t electronSF_4;
    Float_t tauSF_4;
    Bool_t electronIDWP90v2_4;
    Bool_t electronIDWP80v2_4;
    Bool_t muonLoose_4;
    Bool_t muonMedium_4;
    Bool_t muonTight_4;
    Float_t tauVSjetSF_4;
    Float_t tauVSjetSFUp_4;
    Float_t tauVSjetSFDown_4;
    Float_t tauVSeSF_4;
    Float_t tauVSeSFUp_4;
    Float_t tauVSeSFDown_4;
    Float_t tauVSmuSF_4;
    Float_t tauVSmuSFUp_4;
    Float_t tauVSmuSFDown_4;
    
    //met variables
    Float_t met;
    Float_t metphi;
    Float_t metcov00;
    Float_t metcov01;
    Float_t metcov10;
    Float_t metcov11;
    
    Bool_t passMETFilters;

    Float_t DR_12;
    Float_t DR_13;
    Float_t DR_14;
    
    Float_t DR_23;
    Float_t DR_24;
    Float_t DR_34;

    Float_t mll;
    Float_t Z_Pt;
    Float_t Z_DR;
    Bool_t  Z_SS;
    
    Float_t m_vis;
    Float_t pt_tt;
    Float_t H_DR;
    Bool_t  H_SS;
    Float_t m_sv;
    Float_t pt_sv;
    Float_t mt_sv;
    Float_t eta_sv;
    Float_t phi_sv;
    Float_t m_svc;
    Float_t pt_svc;
    Float_t mt_svc;
    Float_t eta_svc;
    Float_t phi_svc;
    Float_t m_fast;
    Float_t pt_fast;
    Float_t mt_fast;
    Float_t eta_fast;
    Float_t phi_fast;
    Float_t m_fastc;
    Float_t pt_fastc;
    Float_t mt_fastc;
    Float_t eta_fastc;
    Float_t phi_fastc;
    
    Float_t Mass_gen;
    Float_t Mass;
    Float_t Mass_sv;
    Float_t res_sv;

    Float_t Mass_svc;
    Float_t res_svc;

    Float_t Mass_fast;
    Float_t res_fast;
    Float_t Mass_fastc;
    Float_t res_fastc;
    
    //Systematics
    Float_t m_svc_mscaleUp;
    Float_t pt_svc_mscaleUp;
    Float_t mt_svc_mscaleUp;
    Float_t eta_svc_mscaleUp;
    Float_t phi_svc_mscaleUp;
    Float_t Mass_svc_mscaleUp;
    Float_t res_svc_mscaleUp;
    
    Float_t m_svc_mscaleDown;
    Float_t pt_svc_mscaleDown;
    Float_t mt_svc_mscaleDown;
    Float_t eta_svc_mscaleDown;
    Float_t phi_svc_mscaleDown;
    Float_t Mass_svc_mscaleDown;
    Float_t res_svc_mscaleDown;
    
    Float_t m_svc_escaleUp;
    Float_t pt_svc_escaleUp;
    Float_t mt_svc_escaleUp;
    Float_t eta_svc_escaleUp;
    Float_t phi_svc_escaleUp;
    Float_t Mass_svc_escaleUp;
    Float_t res_svc_escaleUp;
    
    Float_t m_svc_escaleDown;
    Float_t pt_svc_escaleDown;
    Float_t mt_svc_escaleDown;
    Float_t eta_svc_escaleDown;
    Float_t phi_svc_escaleDown;
    Float_t Mass_svc_escaleDown;
    Float_t res_svc_escaleDown;
    
    Float_t m_svc_tscaleUp;
    Float_t pt_svc_tscaleUp;
    Float_t mt_svc_tscaleUp;
    Float_t eta_svc_tscaleUp;
    Float_t phi_svc_tscaleUp;
    Float_t Mass_svc_tscaleUp;
    Float_t res_svc_tscaleUp;
    
    Float_t m_svc_tscaleDown;
    Float_t pt_svc_tscaleDown;
    Float_t mt_svc_tscaleDown;
    Float_t eta_svc_tscaleDown;
    Float_t phi_svc_tscaleDown;
    Float_t Mass_svc_tscaleDown;
    Float_t res_svc_tscaleDown;
    
    Float_t m_svc_eftscaleUp;
    Float_t pt_svc_eftscaleUp;
    Float_t mt_svc_eftscaleUp;
    Float_t eta_svc_eftscaleUp;
    Float_t phi_svc_eftscaleUp;
    Float_t Mass_svc_eftscaleUp;
    Float_t res_svc_eftscaleUp;
    
    Float_t m_svc_eftscaleDown;
    Float_t pt_svc_eftscaleDown;
    Float_t mt_svc_eftscaleDown;
    Float_t eta_svc_eftscaleDown;
    Float_t phi_svc_eftscaleDown;
    Float_t Mass_svc_eftscaleDown;
    Float_t res_svc_eftscaleDown;
    
    Float_t m_svc_mftscaleUp;
    Float_t pt_svc_mftscaleUp;
    Float_t mt_svc_mftscaleUp;
    Float_t eta_svc_mftscaleUp;
    Float_t phi_svc_mftscaleUp;
    Float_t Mass_svc_mftscaleUp;
    Float_t res_svc_mftscaleUp;
    
    Float_t m_svc_mftscaleDown;
    Float_t pt_svc_mftscaleDown;
    Float_t mt_svc_mftscaleDown;
    Float_t eta_svc_mftscaleDown;
    Float_t phi_svc_mftscaleDown;
    Float_t Mass_svc_mftscaleDown;
    Float_t res_svc_mftscaleDown;
    
    Float_t m_svc_unclmetscaleUp;
    Float_t pt_svc_unclmetscaleUp;
    Float_t mt_svc_unclmetscaleUp;
    Float_t eta_svc_unclmetscaleUp;
    Float_t phi_svc_unclmetscaleUp;
    Float_t Mass_svc_unclmetscaleUp;
    Float_t res_svc_unclmetscaleUp;
    
    Float_t m_svc_unclmetscaleDown;
    Float_t pt_svc_unclmetscaleDown;
    Float_t mt_svc_unclmetscaleDown;
    Float_t eta_svc_unclmetscaleDown;
    Float_t phi_svc_unclmetscaleDown;
    Float_t Mass_svc_unclmetscaleDown;
    Float_t res_svc_unclmetscaleDown;
    
    Float_t m_svc_jecscaleUp;
    Float_t pt_svc_jecscaleUp;
    Float_t mt_svc_jecscaleUp;
    Float_t eta_svc_jecscaleUp;
    Float_t phi_svc_jecscaleUp;
    Float_t Mass_svc_jecscaleUp;
    Float_t res_svc_jecscaleUp;
    
    Float_t m_svc_jecscaleDown;
    Float_t pt_svc_jecscaleDown;
    Float_t mt_svc_jecscaleDown;
    Float_t eta_svc_jecscaleDown;
    Float_t phi_svc_jecscaleDown;
    Float_t Mass_svc_jecscaleDown;
    Float_t res_svc_jecscaleDown;
    
    Float_t m_svc_HEM2018scaleUp;
    Float_t pt_svc_HEM2018scaleUp;
    Float_t mt_svc_HEM2018scaleUp;
    Float_t eta_svc_HEM2018scaleUp;
    Float_t phi_svc_HEM2018scaleUp;
    Float_t Mass_svc_HEM2018scaleUp;
    Float_t res_svc_HEM2018scaleUp;
    
    Float_t m_svc_HEM2018scaleDown;
    Float_t pt_svc_HEM2018scaleDown;
    Float_t mt_svc_HEM2018scaleDown;
    Float_t eta_svc_HEM2018scaleDown;
    Float_t phi_svc_HEM2018scaleDown;
    Float_t Mass_svc_HEM2018scaleDown;
    Float_t res_svc_HEM2018scaleDown;


    UInt_t  npartons;
    
    Float_t fakeRateWeight_1;
    Float_t fakeRateWeightUp_1;
    Float_t fakeRateWeightDown_1;

    Float_t fakeRateWeight_2;
    Float_t fakeRateWeightUp_2;
    Float_t fakeRateWeightDown_2;
    
    Float_t fakeRateWeight_Z1;
    Float_t fakeRateWeightUp_Z1;
    Float_t fakeRateWeightDown_Z1;
    
    Float_t fakeRateWeight_Z2;
    Float_t fakeRateWeightUp_Z2;
    Float_t fakeRateWeightDown_Z2;
    
    Float_t prefiringweight;
    Float_t prefiringweightUp;
    Float_t prefiringweightDown;

    SynTree->Branch("run",&run,"run/l");
    SynTree->Branch("lumi",&lumi,"lumi/i");
    SynTree->Branch("evt",&evt,"evt/l");
    SynTree->Branch("npv",&npv,"npv/I");
    SynTree->Branch("npu",&npu,"npu/I");
    SynTree->Branch("rho",&rho,"rho/F");
    SynTree->Branch("mcweight",&mcweight,"mcweight/F");
    SynTree->Branch("puweight", &puweight, "puweight/F");
    SynTree->Branch("zptweight",&zptweight,"zptweight/F");
    SynTree->Branch("nbtag",&nbtag,"nbtag/I");
    SynTree->Branch("nbtagUp",&nbtagUp,"nbtagUp/I");
    SynTree->Branch("nbtagDown",&nbtagDown,"nbtagDown/I");
    SynTree->Branch("nbtagHEM2018Up",&nbtagHEM2018Up,"nbtagHEM2018Up/I");
    SynTree->Branch("nbtagHEM2018Down",&nbtagHEM2018Down,"nbtagHEM2018Down/I");
    SynTree->Branch("nbtag_mistagUp", &nbtag_mistagUp, "nbtag_mistagUp/I");
    SynTree->Branch("nbtag_mistagDown", &nbtag_mistagDown, "nbtag_mistagDown/I");
    SynTree->Branch("nbtag_btagUp", &nbtag_btagUp, "nbtag_btagUp/I");
    SynTree->Branch("nbtag_btagDown", &nbtag_btagDown, "nbtag_btagDown/I");
    SynTree->Branch("nbtag_noSF", &nbtag_noSF, "nbtag_noSF/I");
    SynTree->Branch("njets",&njets,"njets/I");
    SynTree->Branch("njetsUp",&njetsUp,"njetsUp/I");
    SynTree->Branch("njetsDown",&njetsDown,"njetsDown/I");
    SynTree->Branch("njetsHEM2018Up",&njetsHEM2018Up,"njetsHEM2018Up/I");
    SynTree->Branch("njetsHEM2018Down",&njetsHEM2018Down,"njetsHEM2018Down/I");
    SynTree->Branch("njetspt20",&njetspt20,"njetspt20/I");
    SynTree->Branch("njetspt20Up",&njetspt20Up,"njetspt20Up/I");
    SynTree->Branch("njetspt20Down",&njetspt20Down,"njetspt20Down/I");
    SynTree->Branch("njetspt20HEM2018Up",&njetspt20HEM2018Up,"njetspt20HEM2018Up/I");
    SynTree->Branch("njetspt20HEM2018Down",&njetspt20HEM2018Down,"njetspt20HEM2018Down/I");
    
    SynTree->Branch("jpt_1",&jpt_1,"jpt_1/F");
    SynTree->Branch("jeta_1",&jeta_1,"jeta_1/F");
    SynTree->Branch("jphi_1",&jphi_1,"jphi_1/F");
    SynTree->Branch("jpt_1_Up",&jpt_1_Up,"jpt_1_Up/F");
    SynTree->Branch("jeta_1_Up",&jeta_1_Up,"jeta_1_Up/F");
    SynTree->Branch("jphi_1_Up",&jphi_1_Up,"jphi_1_Up/F");
    SynTree->Branch("jpt_1_Down",&jpt_1_Down,"jpt_1_Down/F");
    SynTree->Branch("jeta_1_Down",&jeta_1_Down,"jeta_1_Down/F");
    SynTree->Branch("jphi_1_Down",&jphi_1_Down,"jphi_1_Down/F");
    SynTree->Branch("jpt_2",&jpt_2,"jpt_2/F");
    SynTree->Branch("jeta_2",&jeta_2,"jeta_2/F");
    SynTree->Branch("jphi_2",&jphi_2,"jphi_2/F");
    SynTree->Branch("jpt_2_Up",&jpt_2_Up,"jpt_2_Up/F");
    SynTree->Branch("jeta_2_Up",&jeta_2_Up,"jeta_2_Up/F");
    SynTree->Branch("jphi_2_Up",&jphi_2_Up,"jphi_2_Up/F");
    SynTree->Branch("jpt_2_Down",&jpt_2_Down,"jpt_2_Down/F");
    SynTree->Branch("jeta_2_Down",&jeta_2_Down,"jeta_2_Down/F");
    SynTree->Branch("jphi_2_Down",&jphi_2_Down,"jphi_2_Down/F");
    SynTree->Branch("bpt_1",&bpt_1,"bpt_1/F");
    SynTree->Branch("beta_1",&beta_1,"beta_1/F");
    SynTree->Branch("bphi_1",&bphi_1,"bphi_1/F");
    SynTree->Branch("bpt_1_Up",&bpt_1_Up,"bpt_1_Up/F");
    SynTree->Branch("beta_1_Up",&beta_1_Up,"beta_1_Up/F");
    SynTree->Branch("bphi_1_Up",&bphi_1_Up,"bphi_1_Up/F");
    SynTree->Branch("bpt_1_Down",&bpt_1_Down,"bpt_1_Down/F");
    SynTree->Branch("beta_1_Down",&beta_1_Down,"beta_1_Down/F");
    SynTree->Branch("bphi_1_Down",&bphi_1_Down,"bphi_1_Down/F");
    SynTree->Branch("bhadronFlavor_1",&bhadronFlavor_1,"bhadronFlavor_1/F");
    SynTree->Branch("bpt_2",&bpt_2,"bpt_2/F");
    SynTree->Branch("beta_2",&beta_2,"beta_2/F");
    SynTree->Branch("bphi_2",&bphi_2,"bphi_2/F");
    SynTree->Branch("bpt_2_Up",&bpt_2_Up,"bpt_2_Up/F");
    SynTree->Branch("beta_2_Up",&beta_2_Up,"beta_2_Up/F");
    SynTree->Branch("bphi_2_Up",&bphi_2_Up,"bphi_2_Up/F");
    SynTree->Branch("bpt_2_Down",&bpt_2_Down,"bpt_2_Down/F");
    SynTree->Branch("beta_2_Down",&beta_2_Down,"beta_2_Down/F");
    SynTree->Branch("bphi_2_Down",&bphi_2_Down,"bphi_2_Down/F");
    SynTree->Branch("bhadronFlavor_2",&bhadronFlavor_2,"bhadronFlavor_2/F");
    
    SynTree->Branch("pt_1",&pt_1,"pt_1/F");
    SynTree->Branch("phi_1",&phi_1,"phi_1/F");
    SynTree->Branch("eta_1",&eta_1,"eta_1/F");
    SynTree->Branch("iso_1",&iso_1,"iso_1/F");
    //SynTree->Branch("iso_ea_1",&iso_ea_1,"iso_ea_1/F");
    SynTree->Branch("gen_match_1",&gen_match_1,"gen_match_1/I");
    SynTree->Branch("IdRawMva_1",&IdRawMva_1,"IdRawMva_1/F");
    SynTree->Branch("pfmt_1",&pfmt_1,"pfmt_1/F");
    SynTree->Branch("muonSF_1",&muonSF_1,"muonSF_1/F");
    SynTree->Branch("IdSF_1",&IdSF_1,"IdSF_1/F");
    SynTree->Branch("TrigSF_1",&TrigSF_1,"TrigSF_1/F");
    SynTree->Branch("electronSF_1",&electronSF_1,"electronSF_1/F");
    SynTree->Branch("muonLoose_1",&muonLoose_1,"muonLoose_1/O");
    SynTree->Branch("muonMedium_1",&muonMedium_1,"muonMedium_1/O");
    SynTree->Branch("muonTight_1",&muonTight_1,"muonTight_1/O");

    SynTree->Branch("pt_2",&pt_2,"pt_2/F");
    SynTree->Branch("phi_2",&phi_2,"phi_2/F");
    SynTree->Branch("eta_2",&eta_2,"eta_2/F");
    SynTree->Branch("iso_2",&iso_2,"iso_2/F");
    //SynTree->Branch("iso_ea_2",&iso_ea_2,"iso_ea_2/F");
    SynTree->Branch("gen_match_2",&gen_match_2,"gen_match_2/I");
    SynTree->Branch("IdRawMva_2",&IdRawMva_2,"IdRawMva_2/F");
    SynTree->Branch("pfmt_2",&pfmt_2,"pfmt_2/F");
    SynTree->Branch("muonSF_2",&muonSF_2,"muonSF_2/F");
    SynTree->Branch("IdSF_2",&IdSF_2,"IdSF_2/F");
    SynTree->Branch("TrigSF_2",&TrigSF_2,"TrigSF_2/F");
    SynTree->Branch("electronSF_2",&electronSF_2,"electronSF_2/F");
    SynTree->Branch("muonLoose_2",&muonLoose_2,"muonLoose_2/O");
    SynTree->Branch("muonMedium_2",&muonMedium_2,"muonMedium_2/O");
    SynTree->Branch("muonTight_2",&muonTight_2,"muonTight_2/O");

    SynTree->Branch("pt_3",&pt_3,"pt_3/F");
    SynTree->Branch("phi_3",&phi_3,"phi_3/F");
    SynTree->Branch("eta_3",&eta_3,"eta_3/F");
    SynTree->Branch("iso_3",&iso_3,"iso_3/F");
    //SynTree->Branch("iso_ea_3",&iso_ea_3,"iso_ea_3/F");
    SynTree->Branch("gen_match_3",&gen_match_3,"gen_match_3/I");
    SynTree->Branch("IdRawMva_3",&IdRawMva_3,"IdRawMva_3/F");
    SynTree->Branch("m_3",&m_3,"m_3/F");
    SynTree->Branch("dm_3",&dm_3,"dm_3/I");
    SynTree->Branch("VSjetVVTight_3",&VSjetVVTight_3,"VSjetVVTight_3/O");
    SynTree->Branch("VSjetVTight_3",&VSjetVTight_3,"VSjetVTight_3/O");
    SynTree->Branch("VSjetTight_3",&VSjetTight_3,"VSjetTight_3/O");
    SynTree->Branch("VSjetMedium_3",&VSjetMedium_3,"VSjetMedium_3/O");
    SynTree->Branch("VSjetLoose_3",&VSjetLoose_3,"VSjetLoose_3/O");
    SynTree->Branch("VSjetVLoose_3",&VSjetVLoose_3,"VSjetVLoose_3/O");
    SynTree->Branch("VSjetVVLoose_3",&VSjetVVLoose_3,"VSjetVVLoose_3/O");
    SynTree->Branch("VSjetVVVLoose_3",&VSjetVVVLoose_3,"VSjetVVVLoose_3/O");
    SynTree->Branch("VSeVVTight_3",&VSeVVTight_3,"VSeVVTight_3/O");
    SynTree->Branch("VSeVTight_3",&VSeVTight_3,"VSeVTight_3/O");
    SynTree->Branch("VSeTight_3",&VSeTight_3,"VSeTight_3/O");
    SynTree->Branch("VSeMedium_3",&VSeMedium_3,"VSeMedium_3/O");
    SynTree->Branch("VSeLoose_3",&VSeLoose_3,"VSeLoose_3/O");
    SynTree->Branch("VSeVLoose_3",&VSeVLoose_3,"VSeVLoose_3/O");
    SynTree->Branch("VSeVVLoose_3",&VSeVVLoose_3,"VSeVVLoose_3/O");
    SynTree->Branch("VSeVVVLoose_3",&VSeVVVLoose_3,"VSeVVVLoose_3/O");
    SynTree->Branch("VSmuTight_3",&VSmuTight_3,"VSmuTight_3/O");
    SynTree->Branch("VSmuMedium_3",&VSmuMedium_3,"VSmuMedium_3/O");
    SynTree->Branch("VSmuLoose_3",&VSmuLoose_3,"VSmuLoose_3/O");
    SynTree->Branch("VSmuVLoose_3",&VSmuVLoose_3,"VSmuVLoose_3/O");
    SynTree->Branch("pfmt_3",&pfmt_3,"pfmt_3/F");
    SynTree->Branch("muonSF_3",&muonSF_3,"muonSF_3/F");
    SynTree->Branch("IdSF_3",&IdSF_3,"IdSF_3/F");
    SynTree->Branch("electronSF_3",&electronSF_3,"electronSF_3/F");
    SynTree->Branch("tauSF_3",&tauSF_3,"tauSF_3/F");
    SynTree->Branch("muonLoose_3",&muonLoose_3,"muonLoose_3/O");
    SynTree->Branch("muonMedium_3",&muonMedium_3,"muonMedium_3/O");
    SynTree->Branch("muonTight_3",&muonTight_3,"muonTight_3/O");
    SynTree->Branch("electronSF_3",&electronSF_3,"electronSF_3/F");
    SynTree->Branch("electronIDWP90v2_3",&electronIDWP90v2_3,"electronIDWP90v2_3/O");
    SynTree->Branch("electronIDWP80v2_3",&electronIDWP80v2_3,"electronIDWP80v2_3/O");
    SynTree->Branch("tauVSjetSF_3",&tauVSjetSF_3,"tauVSjetSF_3/F");
    SynTree->Branch("tauVSjetSFUp_3",&tauVSjetSFUp_3,"tauVSjetSFUp_3/F");
    SynTree->Branch("tauVSjetSFDown_3",&tauVSjetSFDown_3,"tauVSjetSFDown_3/F");
    SynTree->Branch("tauVSeSF_3",&tauVSeSF_3,"tauVSeSF_3/F");
    SynTree->Branch("tauVSeSFUp_3",&tauVSeSFUp_3,"tauVSeSFUp_3/F");
    SynTree->Branch("tauVSeSFDown_3",&tauVSeSFDown_3,"tauVSeSFDown_3/F");
    SynTree->Branch("tauVSmuSF_3",&tauVSmuSF_3,"tauVSmuSF_3/F");
    SynTree->Branch("tauVSmuSFUp_3",&tauVSmuSFUp_3,"tauVSmuSFUp_3/F");
    SynTree->Branch("tauVSmuSFDown_3",&tauVSmuSFDown_3,"tauVSmuSFDown_3/F");
    SynTree->Branch("tauVSjetVVVLooseSF_3",&tauVSjetVVVLooseSF_3,"tauVSjetVVVLooseSF_3/F");
    
    SynTree->Branch("pt_4",&pt_4,"pt_4/F");
    SynTree->Branch("phi_4",&phi_4,"phi_4/F");
    SynTree->Branch("eta_4",&eta_4,"eta_4/F");
    SynTree->Branch("iso_4",&iso_4,"iso_4/F");
    //SynTree->Branch("iso_ea_4",&iso_ea_4,"iso_ea_4/F");
    SynTree->Branch("gen_match_4",&gen_match_4,"gen_match_4/I");
    SynTree->Branch("IdRawMva_4",&IdRawMva_4,"IdRawMva_4/F");
    SynTree->Branch("m_4",&m_4,"m_4/F");
    SynTree->Branch("dm_4",&dm_4,"dm_4/I");
    SynTree->Branch("VSjetVVTight_4",&VSjetVVTight_4,"VSjetVVTight_4/O");
    SynTree->Branch("VSjetVTight_4",&VSjetVTight_4,"VSjetVTight_4/O");
    SynTree->Branch("VSjetTight_4",&VSjetTight_4,"VSjetTight_4/O");
    SynTree->Branch("VSjetMedium_4",&VSjetMedium_4,"VSjetMedium_4/O");
    SynTree->Branch("VSjetLoose_4",&VSjetLoose_4,"VSjetLoose_4/O");
    SynTree->Branch("VSjetVLoose_4",&VSjetVLoose_4,"VSjetVLoose_4/O");
    SynTree->Branch("VSjetVVLoose_4",&VSjetVVLoose_4,"VSjetVVLoose_4/O");
    SynTree->Branch("VSjetVVVLoose_4",&VSjetVVVLoose_4,"VSjetVVVLoose_4/O");
    SynTree->Branch("VSeVVTight_4",&VSeVVTight_4,"VSeVVTight_4/O");
    SynTree->Branch("VSeVTight_4",&VSeVTight_4,"VSeVTight_4/O");
    SynTree->Branch("VSeTight_4",&VSeTight_4,"VSeTight_4/O");
    SynTree->Branch("VSeMedium_4",&VSeMedium_4,"VSeMedium_4/O");
    SynTree->Branch("VSeLoose_4",&VSeLoose_4,"VSeLoose_4/O");
    SynTree->Branch("VSeVLoose_4",&VSeVLoose_4,"VSeVLoose_4/O");
    SynTree->Branch("VSeVVLoose_4",&VSeVVLoose_4,"VSeVVLoose_4/O");
    SynTree->Branch("VSeVVVLoose_4",&VSeVVVLoose_4,"VSeVVVLoose_4/O");
    SynTree->Branch("VSmuTight_4",&VSmuTight_4,"VSmuTight_4/O");
    SynTree->Branch("VSmuMedium_4",&VSmuMedium_4,"VSmuMedium_4/O");
    SynTree->Branch("VSmuLoose_4",&VSmuLoose_4,"VSmuLoose_4/O");
    SynTree->Branch("VSmuVLoose_4",&VSmuVLoose_4,"VSmuVLoose_4/O");
    SynTree->Branch("pfmt_4",&pfmt_4,"pfmt_4/F");
    SynTree->Branch("muonSF_4",&muonSF_4,"muonSF_4/F");
    SynTree->Branch("IdSF_4",&IdSF_4,"IdSF_4/F");
    SynTree->Branch("electronSF_4",&electronSF_4,"electronSF_4/F");
    SynTree->Branch("tauSF_4",&tauSF_4,"tauSF_4/F");
    SynTree->Branch("muonLoose_4",&muonLoose_4,"muonLoose_4/O");
    SynTree->Branch("muonMedium_4",&muonMedium_4,"muonMedium_4/O");
    SynTree->Branch("muonTight_4",&muonTight_4,"muonTight_4/O");
    SynTree->Branch("electronIDWP90v2_4",&electronIDWP90v2_4,"electronIDWP90v2_4/O");
    SynTree->Branch("electronIDWP80v2_4",&electronIDWP80v2_4,"electronIDWP80v2_4/O");
    SynTree->Branch("tauVSjetSF_4",&tauVSjetSF_4,"tauVSjetSF_4/F");
    SynTree->Branch("tauVSjetSFUp_4",&tauVSjetSFUp_4,"tauVSjetSFUp_4/F");
    SynTree->Branch("tauVSjetSFDown_4",&tauVSjetSFDown_4,"tauVSjetSFDown_4/F");
    SynTree->Branch("tauVSeSF_4",&tauVSeSF_4,"tauVSeSF_4/F");
    SynTree->Branch("tauVSeSFUp_4",&tauVSeSFUp_4,"tauVSeSFUp_4/F");
    SynTree->Branch("tauVSeSFDown_4",&tauVSeSFDown_4,"tauVSeSFDown_4/F");
    SynTree->Branch("tauVSmuSF_4",&tauVSmuSF_4,"tauVSmuSF_4/F");
    SynTree->Branch("tauVSmuSFUp_4",&tauVSmuSFUp_4,"tauVSmuSFUp_4/F");
    SynTree->Branch("tauVSmuSFDown_4",&tauVSmuSFDown_4,"tauVSmuSFDown_4/F");
    
    SynTree->Branch("met",&met,"met/F");
    SynTree->Branch("metphi",&metphi,"metphi/F");
    SynTree->Branch("metcov00",&metcov00,"metcov00/F");
    SynTree->Branch("metcov01",&metcov01,"metcov01/F");
    SynTree->Branch("metcov10",&metcov10,"metcov10/F");
    SynTree->Branch("metcov11",&metcov11,"metcov11/F");
    
    SynTree->Branch("passMETFilters",&passMETFilters,"passMETFilters/O");

    SynTree->Branch("DR_12",&DR_12,"DR_12/F");
    SynTree->Branch("DR_13",&DR_13,"DR_13/F");
    SynTree->Branch("DR_14",&DR_14,"DR_14/F");
    SynTree->Branch("DR_23",&DR_23,"DR_23/F");
    SynTree->Branch("DR_24",&DR_24,"DR_24/F");
    SynTree->Branch("DR_34",&DR_34,"DR_34/F");

    SynTree->Branch("mll",&mll,"mll/F");
    SynTree->Branch("Z_Pt",&Z_Pt,"Z_Pt/F");
    SynTree->Branch("Z_DR",&Z_DR,"Z_DR/F");
    SynTree->Branch("Z_SS",&Z_SS,"Z_SS/O");

    SynTree->Branch("m_vis",&m_vis,"m_vis/F");
    SynTree->Branch("pt_tt",&pt_tt,"pt_tt/F");
    SynTree->Branch("H_DR",&H_DR,"H_DR/F");
    SynTree->Branch("H_SS",&H_SS,"H_SS/O");
    SynTree->Branch("m_sv",&m_sv,"m_sv/F");
    SynTree->Branch("pt_sv",&pt_sv,"pt_sv/F");
    SynTree->Branch("mt_sv",&mt_sv,"mt_sv/F");
    SynTree->Branch("eta_sv",&eta_sv,"eta_sv/F");
    SynTree->Branch("phi_sv",&phi_sv,"phi_sv/F");
    SynTree->Branch("m_svc",&m_svc,"m_svc/F");
    SynTree->Branch("pt_svc",&pt_svc,"pt_svc/F");
    SynTree->Branch("mt_svc",&mt_svc,"mt_svc/F");
    SynTree->Branch("eta_svc",&eta_svc,"eta_svc/F");
    SynTree->Branch("phi_svc",&phi_svc,"phi_svc/F");
    SynTree->Branch("m_fast",&m_fast,"m_fast/F");
    SynTree->Branch("pt_fast",&pt_fast,"pt_fast/F");
    SynTree->Branch("mt_fast",&mt_fast,"mt_fast/F");
    SynTree->Branch("eta_fast",&eta_fast,"eta_fast/F");
    SynTree->Branch("phi_fast",&phi_fast,"phi_fast/F");
    SynTree->Branch("m_fastc",&m_fastc,"m_fastc/F");
    SynTree->Branch("pt_fastc",&pt_fastc,"pt_fastc/F");
    SynTree->Branch("mt_fastc",&mt_fastc,"mt_fastc/F");
    SynTree->Branch("eta_fastc",&eta_fastc,"eta_fastc/F");
    SynTree->Branch("phi_fastc",&phi_fastc,"phi_fastc/F");

    SynTree->Branch("Mass_gen",&Mass_gen,"Mass_gen/F");
    SynTree->Branch("Mass",&Mass,"Mass/F");
    SynTree->Branch("Mass_sv",&Mass_sv,"Mass_sv/F");
    SynTree->Branch("res_sv",&res_sv,"res_sv/F");

    SynTree->Branch("Mass_svc",&Mass_svc,"Mass_svc/F");
    SynTree->Branch("res_svc",&res_svc,"res_svc/F");

    SynTree->Branch("Mass_fast",&Mass_fast,"Mass_fast/F");
    SynTree->Branch("res_fast",&res_fast,"res_fast/F");
    
    SynTree->Branch("Mass_fastc",&Mass_fastc,"Mass_fastc/F");
    SynTree->Branch("res_fastc",&res_fastc,"res_fastc/F");

    //Systematics
    SynTree->Branch("m_svc_mscaleUp",&m_svc_mscaleUp,"m_svc_mscaleUp/F");
    SynTree->Branch("pt_svc_mscaleUp",&pt_svc_mscaleUp,"pt_svc_mscaleUp/F");
    SynTree->Branch("mt_svc_mscaleUp",&mt_svc_mscaleUp,"mt_svc_mscaleUp/F");
    SynTree->Branch("eta_svc_mscaleUp",&eta_svc_mscaleUp,"eta_svc_mscaleUp/F");
    SynTree->Branch("phi_svc_mscaleUp",&phi_svc_mscaleUp,"phi_svc_mscaleUp/F");
    SynTree->Branch("Mass_svc_mscaleUp",&Mass_svc_mscaleUp,"Mass_svc_mscaleUp/F");
    SynTree->Branch("res_svc_mscaleUp",&res_svc_mscaleUp,"res_svc_mscaleUp/F");
    
    SynTree->Branch("m_svc_mscaleDown",&m_svc_mscaleDown,"m_svc_mscaleDown/F");
    SynTree->Branch("pt_svc_mscaleDown",&pt_svc_mscaleDown,"pt_svc_mscaleDown/F");
    SynTree->Branch("mt_svc_mscaleDown",&mt_svc_mscaleDown,"mt_svc_mscaleDown/F");
    SynTree->Branch("eta_svc_mscaleDown",&eta_svc_mscaleDown,"eta_svc_mscaleDown/F");
    SynTree->Branch("phi_svc_mscaleDown",&phi_svc_mscaleDown,"phi_svc_mscaleDown/F");
    SynTree->Branch("Mass_svc_mscaleDown",&Mass_svc_mscaleDown,"Mass_svc_mscaleDown/F");
    SynTree->Branch("res_svc_mscaleDown",&res_svc_mscaleDown,"res_svc_mscaleDown/F");
    
    SynTree->Branch("m_svc_escaleUp",&m_svc_escaleUp,"m_svc_escaleUp/F");
    SynTree->Branch("pt_svc_escaleUp",&pt_svc_escaleUp,"pt_svc_escaleUp/F");
    SynTree->Branch("mt_svc_escaleUp",&mt_svc_escaleUp,"mt_svc_escaleUp/F");
    SynTree->Branch("eta_svc_escaleUp",&eta_svc_escaleUp,"eta_svc_escaleUp/F");
    SynTree->Branch("phi_svc_escaleUp",&phi_svc_escaleUp,"phi_svc_escaleUp/F");
    SynTree->Branch("Mass_svc_escaleUp",&Mass_svc_escaleUp,"Mass_svc_escaleUp/F");
    SynTree->Branch("res_svc_escaleUp",&res_svc_escaleUp,"res_svc_escaleUp/F");
    
    SynTree->Branch("m_svc_escaleDown",&m_svc_escaleDown,"m_svc_escaleDown/F");
    SynTree->Branch("pt_svc_escaleDown",&pt_svc_escaleDown,"pt_svc_escaleDown/F");
    SynTree->Branch("mt_svc_escaleDown",&mt_svc_escaleDown,"mt_svc_escaleDown/F");
    SynTree->Branch("eta_svc_escaleDown",&eta_svc_escaleDown,"eta_svc_escaleDown/F");
    SynTree->Branch("phi_svc_escaleDown",&phi_svc_escaleDown,"phi_svc_escaleDown/F");
    SynTree->Branch("Mass_svc_escaleDown",&Mass_svc_escaleDown,"Mass_svc_escaleDown/F");
    SynTree->Branch("res_svc_escaleDown",&res_svc_escaleDown,"res_svc_escaleDown/F");
    
    SynTree->Branch("m_svc_tscaleUp",&m_svc_tscaleUp,"m_svc_tscaleUp/F");
    SynTree->Branch("pt_svc_tscaleUp",&pt_svc_tscaleUp,"pt_svc_tscaleUp/F");
    SynTree->Branch("mt_svc_tscaleUp",&mt_svc_tscaleUp,"mt_svc_tscaleUp/F");
    SynTree->Branch("eta_svc_tscaleUp",&eta_svc_tscaleUp,"eta_svc_tscaleUp/F");
    SynTree->Branch("phi_svc_tscaleUp",&phi_svc_tscaleUp,"phi_svc_tscaleUp/F");
    SynTree->Branch("Mass_svc_tscaleUp",&Mass_svc_tscaleUp,"Mass_svc_tscaleUp/F");
    SynTree->Branch("res_svc_tscaleUp",&res_svc_tscaleUp,"res_svc_tscaleUp/F");
    
    SynTree->Branch("m_svc_tscaleDown",&m_svc_tscaleDown,"m_svc_tscaleDown/F");
    SynTree->Branch("pt_svc_tscaleDown",&pt_svc_tscaleDown,"pt_svc_tscaleDown/F");
    SynTree->Branch("mt_svc_tscaleDown",&mt_svc_tscaleDown,"mt_svc_tscaleDown/F");
    SynTree->Branch("eta_svc_tscaleDown",&eta_svc_tscaleDown,"eta_svc_tscaleDown/F");
    SynTree->Branch("phi_svc_tscaleDown",&phi_svc_tscaleDown,"phi_svc_tscaleDown/F");
    SynTree->Branch("Mass_svc_tscaleDown",&Mass_svc_tscaleDown,"Mass_svc_tscaleDown/F");
    SynTree->Branch("res_svc_tscaleDown",&res_svc_tscaleDown,"res_svc_tscaleDown/F");
    
    SynTree->Branch("m_svc_eftscaleUp",&m_svc_eftscaleUp,"m_svc_eftscaleUp/F");
    SynTree->Branch("pt_svc_eftscaleUp",&pt_svc_eftscaleUp,"pt_svc_eftscaleUp/F");
    SynTree->Branch("mt_svc_eftscaleUp",&mt_svc_eftscaleUp,"mt_svc_eftscaleUp/F");
    SynTree->Branch("eta_svc_eftscaleUp",&eta_svc_eftscaleUp,"eta_svc_eftscaleUp/F");
    SynTree->Branch("phi_svc_eftscaleUp",&phi_svc_eftscaleUp,"phi_svc_eftscaleUp/F");
    SynTree->Branch("Mass_svc_eftscaleUp",&Mass_svc_eftscaleUp,"Mass_svc_eftscaleUp/F");
    SynTree->Branch("res_svc_eftscaleUp",&res_svc_eftscaleUp,"res_svc_eftscaleUp/F");
    
    SynTree->Branch("m_svc_eftscaleDown",&m_svc_eftscaleDown,"m_svc_eftscaleDown/F");
    SynTree->Branch("pt_svc_eftscaleDown",&pt_svc_eftscaleDown,"pt_svc_eftscaleDown/F");
    SynTree->Branch("mt_svc_eftscaleDown",&mt_svc_eftscaleDown,"mt_svc_eftscaleDown/F");
    SynTree->Branch("eta_svc_eftscaleDown",&eta_svc_eftscaleDown,"eta_svc_eftscaleDown/F");
    SynTree->Branch("phi_svc_eftscaleDown",&phi_svc_eftscaleDown,"phi_svc_eftscaleDown/F");
    SynTree->Branch("Mass_svc_eftscaleDown",&Mass_svc_eftscaleDown,"Mass_svc_eftscaleDown/F");
    SynTree->Branch("res_svc_eftscaleDown",&res_svc_eftscaleDown,"res_svc_eftscaleDown/F");
    
    SynTree->Branch("m_svc_mftscaleUp",&m_svc_mftscaleUp,"m_svc_mftscaleUp/F");
    SynTree->Branch("pt_svc_mftscaleUp",&pt_svc_mftscaleUp,"pt_svc_mftscaleUp/F");
    SynTree->Branch("mt_svc_mftscaleUp",&mt_svc_mftscaleUp,"mt_svc_mftscaleUp/F");
    SynTree->Branch("eta_svc_mftscaleUp",&eta_svc_mftscaleUp,"eta_svc_mftscaleUp/F");
    SynTree->Branch("phi_svc_mftscaleUp",&phi_svc_mftscaleUp,"phi_svc_mftscaleUp/F");
    SynTree->Branch("Mass_svc_mftscaleUp",&Mass_svc_mftscaleUp,"Mass_svc_mftscaleUp/F");
    SynTree->Branch("res_svc_mftscaleUp",&res_svc_mftscaleUp,"res_svc_mftscaleUp/F");
    
    SynTree->Branch("m_svc_mftscaleDown",&m_svc_mftscaleDown,"m_svc_mftscaleDown/F");
    SynTree->Branch("pt_svc_mftscaleDown",&pt_svc_mftscaleDown,"pt_svc_mftscaleDown/F");
    SynTree->Branch("mt_svc_mftscaleDown",&mt_svc_mftscaleDown,"mt_svc_mftscaleDown/F");
    SynTree->Branch("eta_svc_mftscaleDown",&eta_svc_mftscaleDown,"eta_svc_mftscaleDown/F");
    SynTree->Branch("phi_svc_mftscaleDown",&phi_svc_mftscaleDown,"phi_svc_mftscaleDown/F");
    SynTree->Branch("Mass_svc_mftscaleDown",&Mass_svc_mftscaleDown,"Mass_svc_mftscaleDown/F");
    SynTree->Branch("res_svc_mftscaleDown",&res_svc_mftscaleDown,"res_svc_mftscaleDown/F");
    
    SynTree->Branch("m_svc_unclmetscaleUp",&m_svc_unclmetscaleUp,"m_svc_unclmetscaleUp/F");
    SynTree->Branch("pt_svc_unclmetscaleUp",&pt_svc_unclmetscaleUp,"pt_svc_unclmetscaleUp/F");
    SynTree->Branch("mt_svc_unclmetscaleUp",&mt_svc_unclmetscaleUp,"mt_svc_unclmetscaleUp/F");
    SynTree->Branch("eta_svc_unclmetscaleUp",&eta_svc_unclmetscaleUp,"eta_svc_unclmetscaleUp/F");
    SynTree->Branch("phi_svc_unclmetscaleUp",&phi_svc_unclmetscaleUp,"phi_svc_unclmetscaleUp/F");
    SynTree->Branch("Mass_svc_unclmetscaleUp",&Mass_svc_unclmetscaleUp,"Mass_svc_unclmetscaleUp/F");
    SynTree->Branch("res_svc_unclmetscaleUp",&res_svc_unclmetscaleUp,"res_svc_unclmetscaleUp/F");
    
    SynTree->Branch("m_svc_unclmetscaleDown",&m_svc_unclmetscaleDown,"m_svc_unclmetscaleDown/F");
    SynTree->Branch("pt_svc_unclmetscaleDown",&pt_svc_unclmetscaleDown,"pt_svc_unclmetscaleDown/F");
    SynTree->Branch("mt_svc_unclmetscaleDown",&mt_svc_unclmetscaleDown,"mt_svc_unclmetscaleDown/F");
    SynTree->Branch("eta_svc_unclmetscaleDown",&eta_svc_unclmetscaleDown,"eta_svc_unclmetscaleDown/F");
    SynTree->Branch("phi_svc_unclmetscaleDown",&phi_svc_unclmetscaleDown,"phi_svc_unclmetscaleDown/F");
    SynTree->Branch("Mass_svc_unclmetscaleDown",&Mass_svc_unclmetscaleDown,"Mass_svc_unclmetscaleDown/F");
    SynTree->Branch("res_svc_unclmetscaleDown",&res_svc_unclmetscaleDown,"res_svc_unclmetscaleDown/F");
    
    SynTree->Branch("m_svc_jecscaleUp",&m_svc_jecscaleUp,"m_svc_jecscaleUp/F");
    SynTree->Branch("pt_svc_jecscaleUp",&pt_svc_jecscaleUp,"pt_svc_jecscaleUp/F");
    SynTree->Branch("mt_svc_jecscaleUp",&mt_svc_jecscaleUp,"mt_svc_jecscaleUp/F");
    SynTree->Branch("eta_svc_jecscaleUp",&eta_svc_jecscaleUp,"eta_svc_jecscaleUp/F");
    SynTree->Branch("phi_svc_jecscaleUp",&phi_svc_jecscaleUp,"phi_svc_jecscaleUp/F");
    SynTree->Branch("Mass_svc_jecscaleUp",&Mass_svc_jecscaleUp,"Mass_svc_jecscaleUp/F");
    SynTree->Branch("res_svc_jecscaleUp",&res_svc_jecscaleUp,"res_svc_jecscaleUp/F");
    
    SynTree->Branch("m_svc_jecscaleDown",&m_svc_jecscaleDown,"m_svc_jecscaleDown/F");
    SynTree->Branch("pt_svc_jecscaleDown",&pt_svc_jecscaleDown,"pt_svc_jecscaleDown/F");
    SynTree->Branch("mt_svc_jecscaleDown",&mt_svc_jecscaleDown,"mt_svc_jecscaleDown/F");
    SynTree->Branch("eta_svc_jecscaleDown",&eta_svc_jecscaleDown,"eta_svc_jecscaleDown/F");
    SynTree->Branch("phi_svc_jecscaleDown",&phi_svc_jecscaleDown,"phi_svc_jecscaleDown/F");
    SynTree->Branch("Mass_svc_jecscaleDown",&Mass_svc_jecscaleDown,"Mass_svc_jecscaleDown/F");
    SynTree->Branch("res_svc_jecscaleDown",&res_svc_jecscaleDown,"res_svc_jecscaleDown/F");
    
    SynTree->Branch("m_svc_HEM2018scaleUp",&m_svc_HEM2018scaleUp,"m_svc_HEM2018scaleUp/F");
    SynTree->Branch("pt_svc_HEM2018scaleUp",&pt_svc_HEM2018scaleUp,"pt_svc_HEM2018scaleUp/F");
    SynTree->Branch("mt_svc_HEM2018scaleUp",&mt_svc_HEM2018scaleUp,"mt_svc_HEM2018scaleUp/F");
    SynTree->Branch("eta_svc_HEM2018scaleUp",&eta_svc_HEM2018scaleUp,"eta_svc_HEM2018scaleUp/F");
    SynTree->Branch("phi_svc_HEM2018scaleUp",&phi_svc_HEM2018scaleUp,"phi_svc_HEM2018scaleUp/F");
    SynTree->Branch("Mass_svc_HEM2018scaleUp",&Mass_svc_HEM2018scaleUp,"Mass_svc_HEM2018scaleUp/F");
    SynTree->Branch("res_svc_HEM2018scaleUp",&res_svc_HEM2018scaleUp,"res_svc_HEM2018scaleUp/F");
    
    SynTree->Branch("m_svc_HEM2018scaleDown",&m_svc_HEM2018scaleDown,"m_svc_HEM2018scaleDown/F");
    SynTree->Branch("pt_svc_HEM2018scaleDown",&pt_svc_HEM2018scaleDown,"pt_svc_HEM2018scaleDown/F");
    SynTree->Branch("mt_svc_HEM2018scaleDown",&mt_svc_HEM2018scaleDown,"mt_svc_HEM2018scaleDown/F");
    SynTree->Branch("eta_svc_HEM2018scaleDown",&eta_svc_HEM2018scaleDown,"eta_svc_HEM2018scaleDown/F");
    SynTree->Branch("phi_svc_HEM2018scaleDown",&phi_svc_HEM2018scaleDown,"phi_svc_HEM2018scaleDown/F");
    SynTree->Branch("Mass_svc_HEM2018scaleDown",&Mass_svc_HEM2018scaleDown,"Mass_svc_HEM2018scaleDown/F");
    SynTree->Branch("res_svc_HEM2018scaleDown",&res_svc_HEM2018scaleDown,"res_svc_HEM2018scaleDown/F");

    SynTree->Branch("npartons",&npartons,"npartons/i");
    
    SynTree->Branch("fakeRateWeight_1",&fakeRateWeight_1,"fakeRateWeight_1/F");
    SynTree->Branch("fakeRateWeightUp_1",&fakeRateWeightUp_1,"fakeRateWeightUp_1/F");
    SynTree->Branch("fakeRateWeightDown_1",&fakeRateWeightDown_1,"fakeRateWeightDown_1/F");
    
    SynTree->Branch("fakeRateWeight_2",&fakeRateWeight_2,"fakeRateWeight_2/F");
    SynTree->Branch("fakeRateWeightUp_2",&fakeRateWeightUp_2,"fakeRateWeightUp_2/F");
    SynTree->Branch("fakeRateWeightDown_2",&fakeRateWeightDown_2,"fakeRateWeightDown_2/F");
    
    SynTree->Branch("fakeRateWeight_Z1",&fakeRateWeight_Z1,"fakeRateWeight_Z1/F");
    SynTree->Branch("fakeRateWeightUp_Z1",&fakeRateWeightUp_Z1,"fakeRateWeightUp_Z1/F");
    SynTree->Branch("fakeRateWeightDown_Z1",&fakeRateWeightDown_Z1,"fakeRateWeightDown_Z1/F");
    
    SynTree->Branch("fakeRateWeight_Z2",&fakeRateWeight_Z2,"fakeRateWeight_Z2/F");
    SynTree->Branch("fakeRateWeightUp_Z2",&fakeRateWeightUp_Z2,"fakeRateWeightUp_Z2/F");
    SynTree->Branch("fakeRateWeightDown_Z2",&fakeRateWeightDown_Z2,"fakeRateWeightDown_Z2/F");
    
    SynTree->Branch("prefiringweight",&prefiringweight,"prefiringweight/F");
    SynTree->Branch("prefiringweightUp",&prefiringweightUp,"prefiringweightUp/F");
    SynTree->Branch("prefiringweightDown",&prefiringweightDown,"prefiringweightDown/F");
    
    
    TH1D * PUweightsOfficialH = new TH1D("PUweightsOfficialH","PU weights w/ official reweighting",1000, 0, 10);
    TH1D * nTruePUInteractionsH = new TH1D("nTruePUInteractionsH","",50,-0.5,49.5);

    // initialize pile up object
    PileUp * PUofficial = new PileUp();
    TString puHistName("pileup");
    if(era == "2017")
    {
        puHistName = sampleName + "_pileup";//for 2017 only
    }
    
    if (applyPUreweighting)
    {
        TFile * filePUdistribution_data = new TFile(TString(cmsswBase)+"/src/"+TString(dataPUFile),"read");
        TFile * filePUdistribution_MC = new TFile (TString(cmsswBase)+"/src/"+TString(mcPUFile), "read");
        TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
        TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get(puHistName);
        PUofficial->set_h_data(PU_data);
        PUofficial->set_h_MC(PU_mc);
    }
    // Muon/electron scale factors
    ScaleFactor * SF_muonIdIso = new ScaleFactor();
    SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
    ScaleFactor * SF_eleIdIso = new ScaleFactor();
    SF_eleIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(EleIdIsoFile));
    ScaleFactor * SF_muonTrig = new ScaleFactor();
    SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTriggerFile));
    
    // set MET filters =================================================================================================================================
    // for recommendations see : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
    std::vector<TString> metFlags; metFlags.clear();
    if (era=="2016"){
       metFlags.push_back("Flag_HBHENoiseFilter");
       metFlags.push_back("Flag_HBHENoiseIsoFilter");
       metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
       metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
       metFlags.push_back("Flag_goodVertices");
       //if (isData) metFlags.push_back("Flag_eeBadScFilter");
       metFlags.push_back("Flag_BadPFMuonFilter");
       // metFlags.push_back("Flag_BadChargedCandidateFilter");  // currently not recommended, under review
    }
    else if (era=="2017"){
       metFlags.push_back("Flag_HBHENoiseFilter");
       metFlags.push_back("Flag_HBHENoiseIsoFilter");
       metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
       metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
       metFlags.push_back("Flag_goodVertices");
       //if (isData) metFlags.push_back("Flag_eeBadScFilter");// not suggested
       metFlags.push_back("Flag_BadPFMuonFilter");
       //metFlags.push_back("Flag_BadChargedCandidateFilter");  // currently not recommended, under review
       metFlags.push_back("ecalBadCalibReducedMINIAODFilter");//version we currently use is outdated, will be updated with next ntuple production campaign
    }
    else if (era == "2018"){
       metFlags.push_back("Flag_goodVertices");
       metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
       metFlags.push_back("Flag_HBHENoiseFilter");
       metFlags.push_back("Flag_HBHENoiseIsoFilter");
       metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
       metFlags.push_back("Flag_BadPFMuonFilter");
       //metFlags.push_back("Flag_BadChargedCandidateFilter");  // currently not recommended, under review
       //if (isData) metFlags.push_back("Flag_eeBadScFilter"); // not suggested
       metFlags.push_back("ecalBadCalibReducedMINIAODFilter");
    }
    else if ((era == "UL2018") || (era == "UL2017")){
       metFlags.push_back("Flag_goodVertices");
       metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
       metFlags.push_back("Flag_HBHENoiseFilter");
       metFlags.push_back("Flag_HBHENoiseIsoFilter");
       metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
       metFlags.push_back("Flag_BadPFMuonFilter");
       metFlags.push_back("Flag_BadPFMuonDzFilter");
       metFlags.push_back("Flag_eeBadScFilter");
       metFlags.push_back("Flag_ecalBadCalibFilter");
       //metFlags.push_back("Flag_BadChargedCandidateFilter");  // currently not recommended, under review
       //if (isData) metFlags.push_back("Flag_eeBadScFilter"); // not suggested
       //metFlags.push_back("ecalBadCalibReducedMINIAODFilter");
    }
    else if ((era == "UL2016postVFP") || (era == "UL2016preVFP")){
       metFlags.push_back("Flag_goodVertices");
       metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
       metFlags.push_back("Flag_HBHENoiseFilter");
       metFlags.push_back("Flag_HBHENoiseIsoFilter");
       metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
       metFlags.push_back("Flag_BadPFMuonFilter");
       metFlags.push_back("Flag_BadPFMuonDzFilter");
       metFlags.push_back("Flag_eeBadScFilter");
       //metFlags.push_back("Flag_BadChargedCandidateFilter");  // currently not recommended, under review
       //if (isData) metFlags.push_back("Flag_eeBadScFilter"); // not suggested
       //metFlags.push_back("ecalBadCalibReducedMINIAODFilter");
    }
    else {
       std::cout << "MET filters not defined for era "<<era<<std::endl;
       exit(-1);
    }
    
    
    // Zpt reweighting for LO DY samples
    TFile * f_zptweight = new TFile(TString(cmsswBase)+"/src/"+ZptweightFile,"read");
    TH2D * h_zptweight = (TH2D*)f_zptweight->Get("zptmass_histo");
    
    //TauIDSF
    std::string year_label;
    if (era == "2016") year_label = "2016Legacy";
    else if (era == "2017") year_label = "2017ReReco";
    else if (era == "2018") year_label = "2018ReReco";
    else if (era == "UL2018") year_label = "UL2018";
    else if (era == "UL2017") year_label = "UL2017";
    else if (era == "UL2016postVFP") year_label = "UL2016_postVFP";
    else if (era == "UL2016preVFP") year_label = "UL2016_preVFP";
    else {std::cout << "year is not 2016, 2017, 2018 or UL - exiting" << '\n'; exit(-1);}
    TauIDSFTool *tauIDSF_VSjetMedium = new TauIDSFTool(year_label, "DeepTau2017v2p1VSjet", "Medium", false);
    TauIDSFTool *tauIDSF_VSjetVVVLoose = new TauIDSFTool(year_label, "DeepTau2017v2p1VSjet", "VVVLoose", false);
    TauIDSFTool *tauIDSF_VSeTight = new TauIDSFTool(year_label, "DeepTau2017v2p1VSe", "Tight", false);
    TauIDSFTool *tauIDSF_VSeVLoose = new TauIDSFTool(year_label, "DeepTau2017v2p1VSe", "VLoose", false);
    TauIDSFTool *tauIDSF_VSmuTight = new TauIDSFTool(year_label, "DeepTau2017v2p1VSmu", "Tight", false);
    TauIDSFTool *tauIDSF_VSmuVLoose = new TauIDSFTool(year_label, "DeepTau2017v2p1VSmu", "VLoose", false);

    //Debug switch
    bool isDebug = false;
    
    //Switch for synch exercise!
    bool SynchExer = false;
    ofstream myfile;
    if(SynchExer==true)
    {
        if(HiggsFinalState == "MT")
            myfile.open ("DESYSyn_mmmt.txt");
        if(HiggsFinalState == "ET")
            myfile.open ("DESYSyn_mmet.txt");
        if(HiggsFinalState == "EM")
            myfile.open ("DESYSyn_mmem.txt");
        if(HiggsFinalState == "TT")
            myfile.open ("DESYSyn_mmtt.txt");
        myfile << "run,lumi,evtid,category" << std::endl;    }
    
  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  int nTotalFiles = 0;
    
    //SV fit
  TH1::AddDirectory(false);
  TFile *inputFile_visPtResolution = new TFile(TString(cmsswBase)+"/src/"+svFitPtResFile.data());
    
    // initialize BTag scale factors =================================================================================================================================
    BTagCalibration calib(bTagAlgorithm, cmsswBase+"/src/"+bTagFile);
    BTagCalibrationReader reader_BTAG(BTagEntry::OP_MEDIUM,"central",{"up","down"});
    TFile * fileTagging = new TFile(TString(cmsswBase)+TString("/src/"+bTagEffFile));
       
    reader_BTAG.load(calib,BTagEntry::FLAV_B,"comb");
    reader_BTAG.load(calib,BTagEntry::FLAV_C,"comb");
    reader_BTAG.load(calib,BTagEntry::FLAV_UDSG,"incl");
    
    float etaBTAG[2] = {0.5,2.1};
    float ptBTAG[5] = {25.,35.,50.,100.,200.};
    
    for (int iEta=0; iEta<2; ++iEta) {
       for (int iPt=0; iPt<5; ++iPt) {
          float sfB = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_B, etaBTAG[iEta], ptBTAG[iPt]);
          float sfC = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_C, etaBTAG[iEta], ptBTAG[iPt]);
          float sfLight = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, etaBTAG[iEta], ptBTAG[iPt]);
          printf("pT = %3.0f   eta = %3.1f  ->  SFb = %5.3f   SFc = %5.3f   SFl = %5.3f\n",ptBTAG[iPt],etaBTAG[iEta],sfB,sfC,sfLight);
       }
    }
    
    TH1F *tagEff_B;
    TH1F *tagEff_C;
    TH1F *tagEff_Light;
    TRandom3 r;
    
    tagEff_B = (TH1F*)fileTagging->Get("btag_eff_b");
    tagEff_C = (TH1F*)fileTagging->Get("btag_eff_c");
    tagEff_Light = (TH1F*)fileTagging->Get("btag_eff_oth");
    
    TString lep1;
    TString lep2;
    if(HiggsFinalState == "EM")
    {
        lep1 ="e";
        lep2 ="m";
    }
    if(HiggsFinalState == "ET")
    {
        lep1 ="e";
        lep2 ="t";
    }
    if(HiggsFinalState == "MT")
    {
        lep1 ="m";
        lep2 ="t";
    }
    if(HiggsFinalState == "TT")
    {
        lep1 ="t";
        lep2 ="t";
    }
    if(HiggsFinalState == "E")//just to let the marco run
    {
        lep1 ="e";
        lep2 ="e";
    }
    if(HiggsFinalState == "M")
    {
        lep1 ="m";
        lep2 ="m";
    }
    if(HiggsFinalState == "T")
    {
        lep1 ="t";
        lep2 ="t";
    }
    FakeRateTools fakerate_lep1(FakeRateFiles1,true,lep1);
    FakeRateTools fakerate_lep2(FakeRateFiles2,true,lep2);
    FakeRateTools fakerate_Zlep1(FakeRateFilesZ,true,"m");
    FakeRateTools fakerate_Zlep2(FakeRateFilesZ,true,"m");

  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  unsigned int RunMin = 9999999;
  unsigned int RunMax = 0;

  std::vector<unsigned int> allRuns; allRuns.clear();

  for (int iF=0; iF<nTotalFiles; ++iF) {
  
    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
      
      TTree * _inittree = NULL;
      _inittree = (TTree*)file_->Get(TString(initNtupleName));
      if (_inittree!=NULL) {
          Float_t genweight;
          if (!isData)
              _inittree->SetBranchAddress("genweight",&genweight);
          Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
          std::cout << "      number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
          for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
              _inittree->GetEntry(iEntry);
              if (isData)
                  histWeightsH->Fill(0.,1.);
              else
                  histWeightsH->Fill(0.,genweight);
          }
      }
      
      TTree * _tree = NULL;
      _tree = (TTree*)file_->Get(TString(ntupleName));
      if (_tree==NULL) continue;
      Long64_t numberOfEntries = _tree->GetEntries();
      std::cout << "      number of entries in Tree      = " << numberOfEntries << std::endl;
      AC1B analysisTree(_tree);
    
    TH1D * histoInputEvents = NULL;
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    if (histoInputEvents==NULL) continue;
      
    int NE = int(histoInputEvents->GetEntries());
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);
    std::cout << "      number of input events         = " << NE << std::endl;

    // EVENT LOOP //
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
      analysisTree.GetEntry(iEntry);
      nEvents++;
      
      if (nEvents%10000==0) 
          cout << "      processed " << nEvents << " events" << endl;

        if(isDebug)
        {
            if(!((analysisTree.event_nr==2439767)))//input here the event(s) you want to check
            {
                continue;
            }
            
            std::cout << "Starting listing check event(s)" << std::endl;
        }
        
        float weight = 1;
        puweight = 1;
        mcweight = 1;
        TrigSF_1 = 1;
        TrigSF_2 = 1;
        prefiringweight = 1;
        prefiringweightUp = 1;
        prefiringweightDown = 1;
        mcweight = analysisTree.genweight;
        prefiringweight = analysisTree.prefiringweight;
        prefiringweightUp = analysisTree.prefiringweightup;
        prefiringweightDown = analysisTree.prefiringweightdown;
        
        fakeRateWeight_1 = 1;
        fakeRateWeightUp_1 = 1;
        fakeRateWeightDown_1 = 1;
        
        fakeRateWeight_2 = 1;
        fakeRateWeightUp_2 = 1;
        fakeRateWeightDown_2 = 1;

        // Pileup variable
        npv = analysisTree.primvertex_count;
        npu = analysisTree.numtruepileupinteractions;
        rho = analysisTree.rho;
        
        if (!isData)
          weight *=analysisTree.genweight;
        npartons = analysisTree.genparticles_noutgoing;

        
        if (!isData) {
            //	cout << analysisTree.numtruepileupinteractions << endl;
            if (applyPUreweighting) {
                nTruePUInteractionsH->Fill(analysisTree.numtruepileupinteractions,weight);
                double Ninteractions = analysisTree.numtruepileupinteractions;
                double PUweight = PUofficial->get_PUweight(Ninteractions);
                weight *= float(PUweight);
                PUweightsOfficialH->Fill(PUweight);
                //	  cout << PUweight << endl;
            }
        }

      if (isData && applyGoodRunSelection){
          bool lumi = false;
          int n=analysisTree.event_run;
          int lum = analysisTree.event_luminosityblock;
	
          std::string num = std::to_string(n);
          std::string lnum = std::to_string(lum);
          for(const auto& a : periods)
          {
	    
              if ( num.c_str() ==  a.name ) {
                  //std::cout<< " Eureka "<<num<<"  "<<a.name<<" ";
                  //     std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      
                  for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
                      //	cout<<b->lower<<"  "<<b->bigger<<endl;
                      if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
                  }
                  auto last = std::prev(a.ranges.end());
                  //    std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
                  if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;
              }
          }
    
          if (!lumi) continue;
          //if (lumi ) cout<<"  =============  Found good run"<<"  "<<n<<"  "<<lum<<endl;
          //std::remove("myinputfile");
      }     

      if (analysisTree.event_run<RunMin)
          RunMin = analysisTree.event_run;
      
      if (analysisTree.event_run>RunMax)
          RunMax = analysisTree.event_run;
      
      //std::cout << " Run : " << analysisTree.event_run << std::endl;
      
      bool isNewRun = true;
      if (allRuns.size()>0) {
          for (unsigned int iR=0; iR<allRuns.size(); ++iR) {
              if (analysisTree.event_run==allRuns.at(iR)) {
                  isNewRun = false;
                  break;
              }
          }
      }
        
        run = analysisTree.event_run;
        lumi = analysisTree.event_luminosityblock;
        evt = analysisTree.event_nr;
        
      if (isNewRun) 
          allRuns.push_back(analysisTree.event_run);
      
      //Trigger information
      bool isTriggerSingleMuon_1 = false; //HLT_IsoMu24_v
      bool isTriggerSingleMuon_2 = false; //HLT_IsoTkMu24_v
      bool isTriggerDoubleMuon_3 = false; //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v
      bool isTriggerDoubleMuon_4 = false; //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
          TString trigName_1(it->first);
          if (trigName_1.Contains(MuonTriggerName_1)) {
              //  std::cout << it->first << " : " << it->second << std::endl;
              if (it->second==1)
                  isTriggerSingleMuon_1 = true;
          }
          TString trigName_2(it->first);
          if (trigName_2.Contains(MuonTriggerName_2))
          {
              //  std::cout << it->first << " : " << it->second << std::endl;
              if (it->second==1)
                  isTriggerSingleMuon_2 = true;
          }
          TString trigName_3(it->first);
          if (trigName_3.Contains(MuonTriggerName_3))
          {
              //  std::cout << it->first << " : " << it->second << std::endl;
              if (it->second==1)
                  isTriggerDoubleMuon_3 = true;
          }
          TString trigName_4(it->first);
          if (trigName_4.Contains(MuonTriggerName_4))
          {
              //  std::cout << it->first << " : " << it->second << std::endl;
              if (it->second==1)
                  isTriggerDoubleMuon_4 = true;
          }
      }

      if (applyTrigger && !isTriggerSingleMuon_1 && !isTriggerSingleMuon_2 && !isTriggerDoubleMuon_3 && !isTriggerDoubleMuon_4) continue;
    
        //Here: count TriggerFound
        if (isData)
          histWeightsTriggerH->Fill(0.,1.);
        else
          histWeightsTriggerH->Fill(0.,analysisTree.genweight);
          
      unsigned int nSingleMuon_1_Filter_a = 0;
      unsigned int nSingleMuon_1_Filter_b = 0;
      unsigned int nSingleMuon_2_Filter = 0;
      unsigned int nDoubleMuon_3_FilterLeg1 = 0;
      unsigned int nDoubleMuon_3_FilterLeg2_a = 0;
      unsigned int nDoubleMuon_3_FilterLeg2_b = 0;
      unsigned int nDoubleMuon_3_FilterDz = 0;
      unsigned int nDoubleMuon_3_FilterMass = 0;
      unsigned int nDoubleMuon_4_FilterLeg1 = 0;
      unsigned int nDoubleMuon_4_FilterLeg2 = 0;
      unsigned int nDoubleMuon_4_FilterDz = 0;
      unsigned int nfilters = analysisTree.run_hltfilters->size();
        
      for (unsigned int i=0; i<nfilters; ++i) {
          TString HLTFilter(analysisTree.run_hltfilters->at(i));
          if (HLTFilter==SingleMuonFilterName_1_a)
          {
              nSingleMuon_1_Filter_a = i;
          }
          if (HLTFilter==SingleMuonFilterName_1_b)
          {
              nSingleMuon_1_Filter_b = i;
          }
          if (HLTFilter==SingleMuonFilterName_2)
          {
              nSingleMuon_2_Filter = i;
          }
          if (HLTFilter==DoubleMuonLeg1FilterName_3)
          {
              nDoubleMuon_3_FilterLeg1 = i;
          }
          if (HLTFilter==DoubleMuonLeg2FilterName_3_a)
          {
              nDoubleMuon_3_FilterLeg2_a = i;
          }
          if (HLTFilter==DoubleMuonLeg2FilterName_3_b)
          {
              nDoubleMuon_3_FilterLeg2_b = i;
          }
          if (HLTFilter==DoubleMuonDzFilterName_3)
          {
              nDoubleMuon_3_FilterDz = i;
          }
          if (HLTFilter==DoubleMuonMassFilterName_3)
          {
              nDoubleMuon_3_FilterMass = i;
          }
          if (HLTFilter==DoubleMuonLeg1FilterName_4)
          {
              nDoubleMuon_4_FilterLeg1 = i;
          }
          if (HLTFilter==DoubleMuonLeg2FilterName_4)
          {
              nDoubleMuon_4_FilterLeg2 = i;
          }
          if (HLTFilter==DoubleMuonDzFilterName_4)
          {
              nDoubleMuon_4_FilterDz = i;
          }
      }
      
        if(!((HiggsFinalState== "MT")||(HiggsFinalState== "ET")||(HiggsFinalState== "EM")||(HiggsFinalState== "EE")||(HiggsFinalState== "MM")||(HiggsFinalState== "TT")||(HiggsFinalState == "T")||(HiggsFinalState == "E")||(HiggsFinalState == "M")))
        {
            std::cout << "you are not considering one of the final states of Higgs decay"<< std::endl;
            exit(-1);
        }
        
        // vertex cuts
        if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
        if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
        float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+analysisTree.primvertex_y*analysisTree.primvertex_y);
        if (dVertex>dVertexCut) continue;
        
        //Here: count vertex cuts
        if (isData)
            histWeightsVertexH->Fill(0.,1.);
        else
            histWeightsVertexH->Fill(0.,analysisTree.genweight);
        
        // searching for btagging discriminant ========================================================================================================================
        unsigned int nBTagDiscriminant1 = 0;
        unsigned int nBTagDiscriminant2 = 0;
        unsigned int nBTagDiscriminant3 = 0;
        SearchForBtagDiscriminant(analysisTree, BTagDiscriminator1, BTagDiscriminator2, BTagDiscriminator3, nBTagDiscriminant1, nBTagDiscriminant2, nBTagDiscriminant3, era);
        
        // electron selection
        vector<unsigned int> idisoEles; idisoEles.clear(); //identified electrons
        vector<unsigned int> goodEles; goodEles.clear(); //good electrons
        for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
            if (analysisTree.electron_pt[ie]<ptEleCut) continue;
            if (fabs(analysisTree.electron_eta[ie])>=etaEleCut) continue; //the "=" is needed because scale factor tools may return error
            if (fabs(analysisTree.electron_dxy[ie])>dxyEleCut) continue;
            if (fabs(analysisTree.electron_dz[ie])>dzEleCut) continue;
            if (!analysisTree.electron_pass_conversion[ie]) continue;
            if (analysisTree.electron_nmissinghits[ie]>=2) continue;
            goodEles.push_back(ie);
            if (!analysisTree.electron_mva_wp90_noIso_Fall17_v2[ie]) continue;
            //float absIso = analysisTree.electron_r03_sumChargedHadronPt[ie];
            //float neutralIso = analysisTree.electron_r03_sumNeutralHadronEt[ie] + analysisTree.electron_r03_sumPhotonEt[ie] - 0.5*analysisTree.electron_r03_sumPUPt[ie];
            //neutralIso = TMath::Max(float(0),neutralIso);
            //absIso += neutralIso;
            //if((absIso/analysisTree.electron_pt[ie])>0.15)continue;
            if((analysisTree.electron_eaIsolation[ie]/analysisTree.electron_pt[ie])>0.15)continue;
            idisoEles.push_back(ie);
            //    cout << "pt:" << analysisTree.electron_pt[ie] << "  passed:" << elePassed << endl;
        }
        
        // muon selection
        vector<unsigned int> idisoMuons; idisoMuons.clear(); //identified muons
        vector<unsigned int> goodMuons; goodMuons.clear(); //good muons
        for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
            if (analysisTree.muon_pt[im]<ptMuonCut) continue;
            //std::cout << "every muon's pt: " << analysisTree.muon_pt[im] << std::endl;
            if (fabs(analysisTree.muon_eta[im])>=etaMuonCut) continue; //the "=" is needed because scale factor tools may return error
            if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
            if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
            if (!(analysisTree.muon_isGlobal[im]||analysisTree.muon_isTracker[im])) continue;
            goodMuons.push_back(im);
            if (!(analysisTree.muon_isLoose[im])) continue;
            float absIso = analysisTree.muon_r04_sumChargedHadronPt[im];
            float neutralIso = analysisTree.muon_r04_sumNeutralHadronEt[im] + analysisTree.muon_r04_sumPhotonEt[im] - 0.5*analysisTree.muon_r04_sumPUPt[im];
            neutralIso = TMath::Max(float(0),neutralIso);
            absIso += neutralIso;
            if(absIso/analysisTree.muon_pt[im]>0.15) continue;
            idisoMuons.push_back(im);
            //	cout << "pt:" << analysisTree.muon_pt[im] << "  passed:" << muPassed << endl;
            //cout << "eta:" << analysisTree.muon_eta[im] << endl;

        }
        
        //tau selection
        vector<unsigned int> goodTaus; goodTaus.clear(); //good taus
        vector<unsigned int> idisoTaus; idisoTaus.clear(); //identified taus
        for(unsigned int it=0;it<analysisTree.tau_count;++it)
        {
            if ((1 + tau_Shift(era, &analysisTree, it)) * analysisTree.tau_pt[it]<ptTauCut) continue;
            if (fabs(analysisTree.tau_eta[it])>etaTauCut) continue;
            if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it])>dzTauCut) continue;
            if (analysisTree.tau_decayModeFindingNewDMs[it]<=0.5) continue;
            if (analysisTree.tau_decayMode[it] == 5 || analysisTree.tau_decayMode[it] == 6) continue;
            //(no needed for Synchronization)
            //if (analysisTree.tau_byMediumDeepTau2017v2p1VSjet[it] <= 0.5) continue;
            if (analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[it] <= 0.5) continue; //For synch! Because nanoAOD has it!
            if (analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[it] <= 0.5) continue;
            if (analysisTree.tau_byVLooseDeepTau2017v2p1VSe[it] <= 0.5) continue;
            goodTaus.push_back(it);
            if (analysisTree.tau_byMediumDeepTau2017v2p1VSjet[it] <= 0.5) continue;
            idisoTaus.push_back(it);
        }
        
        //Final state selection: MuMu
        if(goodMuons.size() < 2) continue;
        
        if(isDebug)
        {
            std::cout << "Passing [goodMuons >= 2]" << std::endl;
        }
        
        //Here: count ZMM Found
        if (isData)
            histWeightsZFoundH->Fill(0.,1.);
        else
            histWeightsZFoundH->Fill(0.,analysisTree.genweight);
        
        //Extra lepton vetoes: Higgs candidates
        if(HiggsFinalState == "EE")
        {
            if(goodEles.size() < 2) continue;
            if(goodMuons.size() < 2) continue;
            if(idisoEles.size() > 2) continue;
            if(idisoMuons.size() > 2) continue;
        }
        if(HiggsFinalState == "EM")
        {
            if(goodEles.size() < 1) continue;
            if(goodMuons.size() < 3) continue;
            if(idisoEles.size() > 1) continue;
            if(idisoMuons.size() > 3) continue;
            if(!isData)
            {
                if(idisoEles.size() < 1) continue;
                if(idisoMuons.size() < 3) continue;
            }
        }
        if(HiggsFinalState == "ET")
        {
            if(goodEles.size() < 1) continue;
            if(goodMuons.size() < 2) continue;
            if(idisoEles.size() > 1) continue;
            if(idisoMuons.size() > 2) continue;
            if(goodTaus.size() < 1) continue;
            if(!isData)
            {
                if(idisoEles.size() < 1) continue;
                if(idisoMuons.size() < 2) continue;
                if(idisoTaus.size() < 1) continue;
            }
        }
        if(HiggsFinalState == "MT")
        {
            if(goodMuons.size() < 3) continue;
            if(idisoEles.size() > 0) continue;
            if(idisoMuons.size() > 3) continue;
            if(goodTaus.size() < 1) continue;
            if(!isData)
            {
                if(idisoMuons.size() < 3) continue;
                if(idisoTaus.size() < 1) continue;
            }
        }
        if(HiggsFinalState == "MM")
        {
            if(goodMuons.size() < 4) continue;
            if(idisoEles.size() > 0) continue;
            if(idisoMuons.size() > 4) continue;
        }
        if(HiggsFinalState == "TT")
        {
            if(goodMuons.size() < 2) continue;
            if(idisoEles.size() > 0) continue;
            if(idisoMuons.size() > 2) continue;
            if(goodTaus.size() < 2) continue;
            if(!isData)
            {
                if(idisoMuons.size() < 2) continue;
                if(idisoTaus.size() < 2) continue;
            }
        }
        if(HiggsFinalState == "T")
        {
            if(goodMuons.size() < 2) continue;
            if(idisoEles.size() > 0) continue;
            if(idisoMuons.size() > 2) continue;
            if(goodTaus.size() < 1) continue;
            if(idisoTaus.size() > 1) continue;
        }
        if(HiggsFinalState == "E")
        {
            if(goodMuons.size() < 2) continue;
            if(goodEles.size() < 1) continue;
            if(idisoEles.size() > 1) continue;
            if(idisoMuons.size() > 2) continue;
            //if(goodTaus.size() != 0) continue;
            if(idisoTaus.size() > 0) continue;
        }
        if(HiggsFinalState == "M")
        {
            if(goodMuons.size() < 3) continue;
            if(idisoEles.size() > 0) continue;
            if(idisoMuons.size() > 3) continue;
            //if(goodTaus.size() != 0) continue;
            if(idisoTaus.size() > 0) continue;
        }
        
        if(isDebug)
        {
            std::cout << "Passing Minimum baseline leptons and maximum id leptons" << std::endl;
        }
        
        if (isData)
            histWeightsExtraVetoH->Fill(0.,1.);
        else
            histWeightsExtraVetoH->Fill(0.,analysisTree.genweight);
        
        
        //std::cout << "id electrons:" << goodEles.size() << std::endl;
        //std::cout << "id muons:" << goodMuons.size() << std::endl;
        //std::cout << "id tau:" << idTaus.size() << std::endl;
        
        //Here: count HiggsCandidate Found
        if (isData)
            histWeightsHiggsFoundH->Fill(0.,1.);
        else
            histWeightsHiggsFoundH->Fill(0.,analysisTree.genweight);

        //Check sepration among the leptons----------

        //Good Z candidates selection
        int indexMuon1 = -1;
        int indexMuon2 = -1;
        bool foundZOS = false;
        TLorentzVector Mu_1, Mu_2;
        float m_MM = 0;
        float m_Dist = 9999;
        for(unsigned int im1=0;im1<goodMuons.size();++im1)
        {
            int indexTemp1 = goodMuons[im1];
            float q1 = analysisTree.muon_charge[indexTemp1];
            for(unsigned int im2=im1+1;im2<goodMuons.size();++im2)
            {
                int indexTemp2 = goodMuons[im2];
                if (indexTemp1 == indexTemp2) continue;
                float q2 = analysisTree.muon_charge[indexTemp2];
                if (q1*q2 > 0) continue;
                float DR_MM = deltaR(analysisTree.muon_eta[indexTemp1],analysisTree.muon_phi[indexTemp1],analysisTree.muon_eta[indexTemp2],analysisTree.muon_phi[indexTemp2]);
                if(DR_MM < 0.3) continue;
                foundZOS = true;
                Mu_1.SetXYZM(analysisTree.muon_px[indexTemp1],analysisTree.muon_py[indexTemp1],analysisTree.muon_pz[indexTemp1],muonMass);
                Mu_2.SetXYZM(analysisTree.muon_px[indexTemp2],analysisTree.muon_py[indexTemp2],analysisTree.muon_pz[indexTemp2],muonMass);
                if (fabs((Mu_1+Mu_2).M() - 91.1876) < m_Dist)
                {
                    m_MM = (Mu_1+Mu_2).M();
                    m_Dist = fabs(m_MM - 91.1876);
                    if (analysisTree.muon_pt[indexTemp1]>analysisTree.muon_pt[indexTemp2])
                    {
                        indexMuon1 = indexTemp1;
                        indexMuon2 = indexTemp2;
                    }
                    else
                    {
                        indexMuon1 = indexTemp2;
                        indexMuon2 = indexTemp1;
                    }
                }
            }
        }
        
        if(foundZOS==false) continue;
        
        if(isDebug)
        {
            std::cout << "Passing Z candidate" << std::endl;
        }
        
        //Here: count OS ZMM Found
        if (isData)
            histWeightsOSZH->Fill(0.,1.);
        else
            histWeightsOSZH->Fill(0.,analysisTree.genweight);
        
        //std::cout << "index muon 1: " << indexMuon1 << std::endl;
        //std::cout << "index muon 2: " << indexMuon2 << std::endl;
        //std::cout << "distance between true Z mass and reco mass: " << m_Dist << std::endl;
        //std::cout << "reco Z mass " << m_MM << std::endl;
        //std::cout << "__________________________________________________________ " << std::endl;

        //Trigger matching--underconstruction
        bool isSingleMuon_1_TriggerFilterMatched = false;
        bool isSingleMuon_2_TriggerFilterMatched = false;
        bool isDoubleMuon_3_TriggerFilterMatched = false;
        bool isDoubleMuon_4_TriggerFilterMatched = false;
        
        bool isMuon1SingleMuon_1_a_Matched = false;
        bool isMuon1SingleMuon_1_b_Matched = false;
        bool isMuon2SingleMuon_1_a_Matched = false;
        bool isMuon2SingleMuon_1_b_Matched = false;

        bool isMuon1SingleMuon_2_Matched = false;
        bool isMuon2SingleMuon_2_Matched = false;
        
        bool isMuon1SingleMuonMatched = false;
        bool isMuon2SingleMuonMatched = false;
        
        bool isMuon1DoubleMuon_3_Leg1Matched = false;
        bool isMuon1DoubleMuon_3_Leg2_a_Matched = false;
        bool isMuon1DoubleMuon_3_Leg2_b_Matched = false;
        bool isMuon2DoubleMuon_3_Leg1Matched = false;
        bool isMuon2DoubleMuon_3_Leg2_a_Matched = false;
        bool isMuon2DoubleMuon_3_Leg2_b_Matched = false;
        bool isMuon1DoubleMuon_3_dzMatched = false;
        bool isMuon2DoubleMuon_3_dzMatched = false;
        bool isMuon1DoubleMuon_3_massMatched = false;
        bool isMuon2DoubleMuon_3_massMatched = false;

        bool isMuon1DoubleMuon_4_Leg1Matched = false;
        bool isMuon1DoubleMuon_4_Leg2Matched = false;
        bool isMuon2DoubleMuon_4_Leg1Matched = false;
        bool isMuon2DoubleMuon_4_Leg2Matched = false;
        bool isMuon1DoubleMuon_4_dzMatched = false;
        bool isMuon2DoubleMuon_4_dzMatched = false;
        
        bool isSingleMuonTriggerFilterMatched = false;
        bool isDoubleMuonTriggerFilterMatched = false;

        if(applyTrigger)
        {
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nSingleMuon_1_Filter_a] && analysisTree.muon_pt[indexMuon1] > singleMuonTriggerPtCut_1)
                    isMuon1SingleMuon_1_a_Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nSingleMuon_1_Filter_b] && analysisTree.muon_pt[indexMuon1] > singleMuonTriggerPtCut_1)
                    isMuon1SingleMuon_1_b_Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nSingleMuon_1_Filter_a] && analysisTree.muon_pt[indexMuon2] > singleMuonTriggerPtCut_1)
                    isMuon2SingleMuon_1_a_Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nSingleMuon_1_Filter_b] && analysisTree.muon_pt[indexMuon2] > singleMuonTriggerPtCut_1)
                    isMuon2SingleMuon_1_b_Matched = true;
            }
            
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nSingleMuon_2_Filter] && analysisTree.muon_pt[indexMuon1] > singleMuonTriggerPtCut_1)
                    isMuon1SingleMuon_2_Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nSingleMuon_2_Filter] && analysisTree.muon_pt[indexMuon2] > singleMuonTriggerPtCut_1)
                    isMuon2SingleMuon_2_Matched = true;
            }
            
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_3_FilterLeg1] && analysisTree.muon_pt[indexMuon1] > doubleMuonLeg1TriggerPtCut_3)
                    isMuon1DoubleMuon_3_Leg1Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_3_FilterLeg2_a] && analysisTree.muon_pt[indexMuon1] > doubleMuonLeg2TriggerPtCut_3)
                    isMuon1DoubleMuon_3_Leg2_a_Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_3_FilterLeg2_b] && analysisTree.muon_pt[indexMuon1] > doubleMuonLeg2TriggerPtCut_3)
                    isMuon1DoubleMuon_3_Leg2_b_Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_3_FilterLeg1] && analysisTree.muon_pt[indexMuon2] > doubleMuonLeg1TriggerPtCut_3)
                    isMuon2DoubleMuon_3_Leg1Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_3_FilterLeg2_a] && analysisTree.muon_pt[indexMuon2] > doubleMuonLeg2TriggerPtCut_3)
                    isMuon2DoubleMuon_3_Leg2_a_Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_3_FilterLeg2_b] && analysisTree.muon_pt[indexMuon2] > doubleMuonLeg2TriggerPtCut_3)
                    isMuon2DoubleMuon_3_Leg2_b_Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_3_FilterDz])
                    isMuon1DoubleMuon_3_dzMatched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_3_FilterDz])
                    isMuon2DoubleMuon_3_dzMatched = true;
            }
            
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_3_FilterMass])
                    isMuon1DoubleMuon_3_massMatched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_3_FilterMass])
                    isMuon2DoubleMuon_3_massMatched = true;
            }
            
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_4_FilterLeg1] && analysisTree.muon_pt[indexMuon1] > doubleMuonLeg1TriggerPtCut_4)
                    isMuon1DoubleMuon_4_Leg1Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_4_FilterLeg2] && analysisTree.muon_pt[indexMuon1] > doubleMuonLeg2TriggerPtCut_4)
                    isMuon1DoubleMuon_4_Leg2Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_4_FilterLeg1] && analysisTree.muon_pt[indexMuon2] > doubleMuonLeg1TriggerPtCut_4)
                    isMuon2DoubleMuon_4_Leg1Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_4_FilterLeg2] && analysisTree.muon_pt[indexMuon2] > doubleMuonLeg2TriggerPtCut_4)
                    isMuon2DoubleMuon_4_Leg2Matched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_4_FilterDz])
                    isMuon1DoubleMuon_4_dzMatched = true;
            }
            for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
            {
                float dRtrig = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
                if (dRtrig>DRTrigMatch) continue;
                if (analysisTree.trigobject_filters[iT][nDoubleMuon_4_FilterDz])
                    isMuon2DoubleMuon_4_dzMatched = true;
            }
            
            //check which muon matching which trigger
            if((isMuon1SingleMuon_1_a_Matched || isMuon1SingleMuon_1_b_Matched) || (isMuon2SingleMuon_1_a_Matched || isMuon2SingleMuon_1_b_Matched))
                isSingleMuon_1_TriggerFilterMatched = true;
            if(isMuon1SingleMuon_2_Matched || isMuon2SingleMuon_2_Matched)
                isSingleMuon_2_TriggerFilterMatched = true;
            
            if((isMuon1SingleMuon_1_a_Matched || isMuon1SingleMuon_1_b_Matched) || isMuon1SingleMuon_2_Matched)
                isMuon1SingleMuonMatched = true;
            if((isMuon2SingleMuon_1_a_Matched || isMuon2SingleMuon_1_b_Matched) || isMuon2SingleMuon_2_Matched)
                isMuon2SingleMuonMatched = true;
            
            if((era == "2016") || (era == "UL2016preVFP") || (era == "UL2016postVFP"))
            {
                if(isMuon1DoubleMuon_3_Leg1Matched && (isMuon2DoubleMuon_3_Leg2_a_Matched || isMuon2DoubleMuon_3_Leg2_b_Matched) && isMuon1DoubleMuon_3_dzMatched && isMuon2DoubleMuon_3_dzMatched)
                    isDoubleMuon_3_TriggerFilterMatched = true;
            }
            else
            {
                if(isMuon1DoubleMuon_3_Leg1Matched && (isMuon2DoubleMuon_3_Leg2_a_Matched || isMuon2DoubleMuon_3_Leg2_b_Matched) && isMuon1DoubleMuon_3_dzMatched && isMuon2DoubleMuon_3_dzMatched && isMuon1DoubleMuon_3_massMatched && isMuon2DoubleMuon_3_massMatched)
                    isDoubleMuon_3_TriggerFilterMatched = true;
            }
            if((era == "2016") || (era == "UL2016preVFP") || (era == "UL2016postVFP"))
            {
                if(isMuon2DoubleMuon_3_Leg1Matched && (isMuon1DoubleMuon_3_Leg2_a_Matched || isMuon1DoubleMuon_3_Leg2_b_Matched) && isMuon1DoubleMuon_3_dzMatched && isMuon2DoubleMuon_3_dzMatched)
                    isDoubleMuon_3_TriggerFilterMatched = true;
            }
            else
            {
                if(isMuon2DoubleMuon_3_Leg1Matched && (isMuon1DoubleMuon_3_Leg2_a_Matched || isMuon1DoubleMuon_3_Leg2_b_Matched) && isMuon1DoubleMuon_3_dzMatched && isMuon2DoubleMuon_3_dzMatched && isMuon1DoubleMuon_3_massMatched && isMuon2DoubleMuon_3_massMatched)
                    isDoubleMuon_3_TriggerFilterMatched = true;
            }
            if(isMuon1DoubleMuon_4_Leg1Matched && isMuon2DoubleMuon_4_Leg2Matched && isMuon1DoubleMuon_4_dzMatched && isMuon2DoubleMuon_4_dzMatched)
                isDoubleMuon_4_TriggerFilterMatched = true;
                
            if(isMuon2DoubleMuon_4_Leg1Matched && isMuon1DoubleMuon_4_Leg2Matched && isMuon1DoubleMuon_4_dzMatched && isMuon2DoubleMuon_4_dzMatched)
                isDoubleMuon_4_TriggerFilterMatched = true;
                
            //
            if (applySingleMuTriggerOnly == true)//this setting is for switching on and off single electron trigger
            {
                if((era == "2018") || (era == "2017")||(era == "UL2018") || (era == "UL2017"))
                {
                    if(isSingleMuon_1_TriggerFilterMatched == true)
                        isSingleMuonTriggerFilterMatched = true;
                }
                if((era == "2016") || (era == "UL2016preVFP") || (era == "UL2016postVFP") )
                {
                    if((isSingleMuon_1_TriggerFilterMatched == true) || (isSingleMuon_2_TriggerFilterMatched == true))
                        isSingleMuonTriggerFilterMatched = true;
                }
                if(isSingleMuonTriggerFilterMatched == false)
                    continue;
            }
            else//"OR" double trigger apply
            {
                if((era == "2018") || (era == "2017")||(era == "UL2018") || (era == "UL2017"))
                {
                    if(isSingleMuon_1_TriggerFilterMatched == true)
                        isSingleMuonTriggerFilterMatched = true;
                        
                    if (isDoubleMuon_3_TriggerFilterMatched == true)
                        isDoubleMuonTriggerFilterMatched = true;
                    
                    if((isSingleMuonTriggerFilterMatched==false) && (isDoubleMuonTriggerFilterMatched==false))
                        continue;
                }
                if((era == "2016") || (era == "UL2016preVFP") || (era == "UL2016postVFP") )
                {
                    if((isSingleMuon_1_TriggerFilterMatched==true) || (isSingleMuon_2_TriggerFilterMatched==true))
                        isSingleMuonTriggerFilterMatched = true;
                    
                    if((isDoubleMuon_3_TriggerFilterMatched==true) || (isDoubleMuon_4_TriggerFilterMatched==true))
                        isDoubleMuonTriggerFilterMatched = true;
                    
                    if ((isSingleMuonTriggerFilterMatched==false) && (isDoubleMuonTriggerFilterMatched==false) )
                    continue;
                }
            }
        }
        //Here: count trigger objects matched
        if (isData)
            histWeightsTriggerMatchedH->Fill(0.,1.);
        else
            histWeightsTriggerMatchedH->Fill(0.,analysisTree.genweight);
        
        if(isDebug)
        {
            std::cout << "Passing trigger filter matched" << std::endl;
        }
          
        //Higgs boson candidates selection
        int indexLeg3 = -1;
        int indexLeg4 = -1;
        float Higgs_LT = -9999;
        if(HiggsFinalState == "EM")
        {
            for(unsigned int ie=0;ie<goodEles.size();++ie)
            {
                int indexTemp1 = goodEles[ie];
                float DR_1E = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.electron_eta[indexTemp1],analysisTree.electron_phi[indexTemp1]);
                if(DR_1E < 0.3) continue;
                float DR_2E = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.electron_eta[indexTemp1],analysisTree.electron_phi[indexTemp1]);
                if(DR_2E < 0.3) continue;
                for(unsigned int im=0;im<goodMuons.size();++im)
                {
                    int indexTemp2 = goodMuons[im];
                    if ((indexTemp2 == indexMuon1)||(indexTemp2 == indexMuon2))
                        continue;
                    float DR_1M = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.muon_eta[indexTemp2],analysisTree.muon_phi[indexTemp2]);
                    if(DR_1M < 0.3) continue;
                    float DR_2M = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.muon_eta[indexTemp2],analysisTree.muon_phi[indexTemp2]);
                    if(DR_2M < 0.3) continue;
                    float DR_EM = deltaR(analysisTree.electron_eta[indexTemp1],analysisTree.electron_phi[indexTemp1],analysisTree.muon_eta[indexTemp2],analysisTree.muon_phi[indexTemp2]);
                    if(DR_EM < 0.3) continue;
                    if (analysisTree.electron_pt[indexTemp1]+analysisTree.muon_pt[indexTemp2]>Higgs_LT)
                    {
                        Higgs_LT = analysisTree.electron_pt[indexTemp1] + analysisTree.muon_pt[indexTemp2];
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp2;
                    }
                }
            }
            if((indexLeg3 == -1)||(indexLeg4 == -1))
                continue;
        }

        if(HiggsFinalState =="ET")
        {
            for(unsigned int ie=0;ie<goodEles.size();++ie)
            {
                int indexTemp1 = goodEles[ie];
                float DR_1E = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.electron_eta[indexTemp1],analysisTree.electron_phi[indexTemp1]);
                if(DR_1E < 0.3) continue;
                float DR_2E = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.electron_eta[indexTemp1],analysisTree.electron_phi[indexTemp1]);
                if(DR_2E < 0.3) continue;
                for(unsigned int it=0;it<goodTaus.size();++it)
                {
                    int indexTemp2 = goodTaus[it];
                    //if (analysisTree.tau_byTightDeepTau2017v2p1VSe[indexTemp2]<0.5) continue;
                    //(no needed for Synchronization)
                    //if (analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexTemp2]<0.5) continue;
                    float DR_1Tau = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    float DR_2Tau = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    float DR_EleTau = deltaR(analysisTree.electron_eta[indexTemp1],analysisTree.electron_phi[indexTemp1],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    if((DR_1Tau < 0.5)||(DR_2Tau < 0.5)||(DR_EleTau < 0.5))//Select tau that is not overlapped with others first, then LT...
                    {
                        continue;
                    }

                    if (analysisTree.electron_pt[indexTemp1] + (1 + tau_Shift(era, &analysisTree, indexTemp2))*analysisTree.tau_pt[indexTemp2]>Higgs_LT)
                    {
                        Higgs_LT = analysisTree.electron_pt[indexTemp1] + (1 + tau_Shift(era, &analysisTree, indexTemp2))*analysisTree.tau_pt[indexTemp2];
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp2;
                    }
                }
            }
            if((indexLeg3 == -1)||(indexLeg4 == -1))
                continue;
        }
        
        if(HiggsFinalState =="E")
        {
            for(unsigned int ie=0;ie<goodEles.size();++ie)
            {
                int indexTemp1 = goodEles[ie];
                float DR_1M1 = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.electron_eta[indexTemp1],analysisTree.electron_phi[indexTemp1]);
                float DR_2M1 = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.electron_eta[indexTemp1],analysisTree.electron_phi[indexTemp1]);
                if((DR_1M1 < 0.3)||(DR_2M1 < 0.3))
                {
                    continue;
                }
                if(indexLeg3!= -1)
                {
                    if(analysisTree.electron_pt[indexTemp1]>analysisTree.electron_pt[indexLeg3])
                    {
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp1; //setting the 4th leg equal to 3rd leg in case
                    }
                }
                else
                {
                    indexLeg3 = indexTemp1;
                    indexLeg4 = indexTemp1;
                }
            }
            if(indexLeg3 == -1)
                continue;
        }
        
        if(HiggsFinalState == "MT")
        {
            //std::cout << "Muon size: " << goodMuons.size() << std::endl;
            for(unsigned int im=0;im<goodMuons.size();++im)
            {
                int indexTemp1 = goodMuons[im];
                if ((indexTemp1 == indexMuon1)||(indexTemp1 == indexMuon2))
                    continue;
                float DR_1M = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.muon_eta[indexTemp1],analysisTree.muon_phi[indexTemp1]);
                if(DR_1M < 0.3) continue;
                float DR_2M = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.muon_eta[indexTemp1],analysisTree.muon_phi[indexTemp1]);
                if(DR_2M < 0.3) continue;
                for(unsigned int it=0;it<goodTaus.size();++it)
                {
                    int indexTemp2 = goodTaus[it];
                    //indexLeg4 = indexTemp2;
                    //if (analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexTemp2]<0.5) continue;
                    float DR_1Tau = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    float DR_2Tau = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    float DR_MuTau = deltaR(analysisTree.muon_eta[indexTemp1],analysisTree.muon_phi[indexTemp1],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    if((DR_1Tau < 0.5)||(DR_2Tau < 0.5)||(DR_MuTau < 0.5))//Select tau that is not overlapped with others first, then LT...
                    {
                        continue;
                    }
                    if (analysisTree.electron_pt[indexTemp1] + (1 + tau_Shift(era, &analysisTree, indexTemp2))*analysisTree.tau_pt[indexTemp2]>Higgs_LT)
                    {
                        Higgs_LT = analysisTree.electron_pt[indexTemp1] + (1 + tau_Shift(era, &analysisTree, indexTemp2))*analysisTree.tau_pt[indexTemp2];
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp2;
                    }
                }
            }
            if((indexLeg3 == -1)||(indexLeg4 == -1))
                continue;
        }
        
        if(HiggsFinalState =="M")
        {
            for(unsigned int im=0;im<goodMuons.size();++im)
            {
                int indexTemp1 = goodMuons[im];
                if ((indexTemp1 == indexMuon1)||(indexTemp1 == indexMuon2))
                    continue;
                float DR_1M1 = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.muon_eta[indexTemp1],analysisTree.muon_phi[indexTemp1]);
                float DR_2M1 = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.muon_eta[indexTemp1],analysisTree.muon_phi[indexTemp1]);
                if((DR_1M1 < 0.3)||(DR_2M1 < 0.3))
                {
                    continue;
                }
                if(indexLeg3!= -1)
                {
                    if(analysisTree.muon_pt[indexTemp1]>analysisTree.muon_pt[indexLeg3])
                    {
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp1; //setting the 4th leg equal to 3rd leg in case
                    }
                }
                else
                {
                    indexLeg3 = indexTemp1;
                    indexLeg4 = indexTemp1;
                }
            }
            if(indexLeg3 == -1)
                continue;
        }
        
        if(HiggsFinalState =="TT")
        {
            for(unsigned int it=0;it<goodTaus.size();++it)
            {
                int indexTemp1 = goodTaus[it];
                float DR_1Tau1 = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.tau_eta[indexTemp1],analysisTree.tau_phi[indexTemp1]);
                float DR_2Tau1 = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.tau_eta[indexTemp1],analysisTree.tau_phi[indexTemp1]);
                if((DR_1Tau1 < 0.5)||(DR_2Tau1 < 0.5))//Select 1st tau that is not overlapped with electrons first, then next tau...
                {
                    continue;
                }
                for(unsigned int it2=it+1;it2<goodTaus.size();++it2)
                {
                    int indexTemp2 = goodTaus[it2];
                    float DR_1Tau2 = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    float DR_2Tau2 = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    float DR_TauTau = deltaR(analysisTree.tau_eta[indexTemp1],analysisTree.tau_phi[indexTemp1],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    if((DR_1Tau2 < 0.5)||(DR_2Tau2 < 0.5)||(DR_TauTau < 0.5))//Select tau that is not overlapped with others first, then LT...
                    {
                        continue;
                    }

                    if((1 + tau_Shift(era, &analysisTree, indexTemp1))*analysisTree.tau_pt[indexTemp1] + (1 + tau_Shift(era, &analysisTree, indexTemp2))*analysisTree.tau_pt[indexTemp2]>Higgs_LT)
                    {
                        Higgs_LT = (1 + tau_Shift(era, &analysisTree, indexTemp1))*analysisTree.tau_pt[indexTemp1] + (1 + tau_Shift(era, &analysisTree, indexTemp2))*analysisTree.tau_pt[indexTemp2];
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp2;
                    }
                }
            }
            if((indexLeg3 == -1)||(indexLeg4 == -1))
                continue;
        }
        
        if(HiggsFinalState == "T")
        {
            for(unsigned int it=0;it<goodTaus.size();++it)
            {
                int indexTemp1 = goodTaus[it];
                float DR_1Tau1 = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.tau_eta[indexTemp1],analysisTree.tau_phi[indexTemp1]);
                float DR_2Tau1 = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.tau_eta[indexTemp1],analysisTree.tau_phi[indexTemp1]);
                if((DR_1Tau1 < 0.5)||(DR_2Tau1 < 0.5))//Select 1st tau that is not overlapped with electrons first, then next tau...
                {
                    continue;
                }
                
                if(indexLeg3!= -1)
                {
                    if(analysisTree.tau_pt[indexTemp1]>analysisTree.tau_pt[indexLeg3])
                    {
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp1; //setting the 4th leg equal to 3rd leg in case
                    }
                }
                else
                {
                    indexLeg3 = indexTemp1;
                    indexLeg4 = indexTemp1;
                }
            }
            if(indexLeg3 == -1)
                continue;
        }

        if(isDebug)
        {
            std::cout << "Passing Higgs candiate" << std::endl;
        }
        //(No need for Synchronization)
        //if(Higgs_LT<=40) continue;
        //if(HiggsFinalState =="TT")
        //{
        //    if(Higgs_LT < 60) continue;
        //}
                
        //std::cout << "index muon 1: " << indexLeg3 << std::endl;
        //std::cout << "index tau 1: " << indexLeg4 << std::endl;
        //std::cout << "Higgs LT (pT_muon + pT_tau): " << Higgs_LT << std::endl;
        //std::cout << "__________________________________________________________ " << std::endl;
        //Filling up lepton kinematics
        if((!isData) || (HiggsFinalState == "E") || (HiggsFinalState == "M") || (HiggsFinalState == "T"))//filter out by signal region ID selection for MC
        {
            if(!analysisTree.muon_isLoose[indexMuon1]) continue;
            if(!analysisTree.muon_isLoose[indexMuon2]) continue;
            if(HiggsFinalState == "EM")
            {
                if(!analysisTree.electron_mva_wp90_noIso_Fall17_v2[indexLeg3]) continue;
                if(!analysisTree.muon_isLoose[indexLeg4]) continue;
            }
            if(HiggsFinalState == "ET")
            {
                if(!analysisTree.electron_mva_wp90_noIso_Fall17_v2[indexLeg3]) continue;
                if(!analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexLeg4]) continue;
                if(!analysisTree.tau_byTightDeepTau2017v2p1VSe[indexLeg4]) continue;
            }
            if(HiggsFinalState == "MT")
            {
                if(!analysisTree.muon_isLoose[indexLeg3]) continue;
                if(isDebug)
                {
                    std::cout << "Passing muon id leg3 imposed" << std::endl;
                    std::cout << "tau pt leg4:" << analysisTree.tau_pt[indexLeg4] << std::endl;
                }
                if(!analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexLeg4]) continue;
                if(isDebug)
                {
                    std::cout << "Passing Medium VSjet leg4 imposed" << std::endl;
                }
                if(!analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexLeg4]) continue;
                if(isDebug)
                {
                    std::cout << "Passing Tight VSmu leg4 imposed" << std::endl;
                }
            }
            if(HiggsFinalState == "TT")
            {
                if(!analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexLeg3]) continue;
                if(!analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexLeg4]) continue;
            }
        }
        
        if(isDebug)
        {
            std::cout << "Passing tight id imposed" << std::endl;
        }
        
        pt_1 = analysisTree.muon_pt[indexMuon1];
        phi_1 = analysisTree.muon_phi[indexMuon1];
        eta_1 = analysisTree.muon_eta[indexMuon1];
        
        float absIso_1 = analysisTree.muon_r04_sumChargedHadronPt[indexMuon1];
        float neutralIso_1 = analysisTree.muon_r04_sumNeutralHadronEt[indexMuon1] + analysisTree.muon_r04_sumPhotonEt[indexMuon1] - 0.5*analysisTree.muon_r04_sumPUPt[indexMuon1];
        neutralIso_1 = TMath::Max(float(0),neutralIso_1);
        absIso_1 += neutralIso_1;
        iso_1 = absIso_1/analysisTree.muon_pt[indexMuon1];
        
        if(!isData)
        {
            IdSF_1 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[indexMuon1]),double(analysisTree.muon_eta[indexMuon1]));
            //std::cout << "Muon 1 SFs, pt of input muon: " << analysisTree.muon_pt[indexMuon1] << std::endl;
            if(isSingleMuonTriggerFilterMatched && isMuon1SingleMuonMatched)
                TrigSF_1 = (float)SF_muonTrig->get_ScaleFactor(double(analysisTree.muon_pt[indexMuon1]),double(analysisTree.muon_eta[indexMuon1]));
        }
        muonLoose_1 = analysisTree.muon_isLoose[indexMuon1];
        muonMedium_1 = analysisTree.muon_isMedium[indexMuon1];
        muonTight_1 = analysisTree.muon_isTight[indexMuon1];
        
        pt_2 = analysisTree.muon_pt[indexMuon2];
        phi_2 = analysisTree.muon_phi[indexMuon2];
        eta_2 = analysisTree.muon_eta[indexMuon2];
        
        float absIso_2 = analysisTree.muon_r04_sumChargedHadronPt[indexMuon2];
        float neutralIso_2 = analysisTree.muon_r04_sumNeutralHadronEt[indexMuon2] + analysisTree.muon_r04_sumPhotonEt[indexMuon2] - 0.5*analysisTree.muon_r04_sumPUPt[indexMuon2];
        neutralIso_2 = TMath::Max(float(0),neutralIso_2);
        absIso_2 += neutralIso_2;
        iso_2 = absIso_2/analysisTree.muon_pt[indexMuon2];
        
        if(!isData)
        {
            IdSF_2 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[indexMuon2]),double(analysisTree.muon_eta[indexMuon2]));
            //std::cout << "Muon 2 SFs, pt of input muon: " << analysisTree.muon_pt[indexMuon2] << std::endl;
            if(isSingleMuonTriggerFilterMatched && isMuon2SingleMuonMatched)
                TrigSF_2 = (float)SF_muonTrig->get_ScaleFactor(double(analysisTree.muon_pt[indexMuon2]),double(analysisTree.muon_eta[indexMuon2]));
        }
        muonLoose_2 = analysisTree.muon_isLoose[indexMuon2];
        muonMedium_2 = analysisTree.muon_isMedium[indexMuon2];
        muonTight_2 = analysisTree.muon_isTight[indexMuon2];
        
        Z_DR = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2]);
        Z_SS = false;
        if (analysisTree.muon_charge[indexMuon1]*analysisTree.muon_charge[indexMuon2]>0)
            Z_SS = true;
        
        TLorentzVector Leg3;
        TLorentzVector Leg4;
        TLorentzVector Leg3_unCorrected_1;
        TLorentzVector Leg4_unCorrected_1;
        TLorentzVector Leg3_unCorrected_2;
        TLorentzVector Leg4_unCorrected_2;
        TLorentzVector Leg3_unCorrected_3;
        TLorentzVector Leg4_unCorrected_3;
        TLorentzVector Leg3_unCorrected_4;
        TLorentzVector Leg4_unCorrected_4;
        TLorentzVector Leg3_unCorrected_5;
        TLorentzVector Leg4_unCorrected_5;
        TLorentzVector Leg3_unCorrected_6;
        TLorentzVector Leg4_unCorrected_6;
        if(HiggsFinalState == "EM")
        {
            pt_3 = analysisTree.electron_pt[indexLeg3];
            phi_3 = analysisTree.electron_phi[indexLeg3];
            eta_3 = analysisTree.electron_eta[indexLeg3];
            m_3 = classic_svFit::electronMass;
            //float absIso_3 = analysisTree.electron_r03_sumChargedHadronPt[indexLeg3];
            //float neutralIso_3 = analysisTree.electron_r03_sumNeutralHadronEt[indexLeg3] + analysisTree.electron_r03_sumPhotonEt[indexLeg3] - 0.5*analysisTree.electron_r03_sumPUPt[indexLeg3];
            //neutralIso_3 = TMath::Max(float(0),neutralIso_3);
            //absIso_3 += neutralIso_3;
            //iso_3 = absIso_3/analysisTree.electron_pt[indexLeg3];
            iso_3 = analysisTree.electron_eaIsolation[indexLeg3]/analysisTree.electron_pt[indexLeg3];
            if(!isData)
            {
                IdSF_3 = (float)SF_eleIdIso->get_ScaleFactor(double(analysisTree.electron_pt[indexLeg3]),double(analysisTree.electron_eta[indexLeg3]));
            }
            Leg3.SetXYZM(analysisTree.electron_px[indexLeg3],analysisTree.electron_py[indexLeg3],analysisTree.electron_pz[indexLeg3],electronMass);
            electronIDWP80v2_3 = analysisTree.electron_mva_wp80_noIso_Fall17_v2[indexLeg3];
            electronIDWP90v2_3 = analysisTree.electron_mva_wp90_noIso_Fall17_v2[indexLeg3];

            pt_4 = analysisTree.muon_pt[indexLeg4];
            phi_4 = analysisTree.muon_phi[indexLeg4];
            eta_4 = analysisTree.muon_eta[indexLeg4];
            m_4 = classic_svFit::muonMass;
            float absIso_4 = analysisTree.muon_r04_sumChargedHadronPt[indexLeg4];
            float neutralIso_4 = analysisTree.muon_r04_sumNeutralHadronEt[indexLeg4] + analysisTree.muon_r04_sumPhotonEt[indexLeg4] - 0.5*analysisTree.muon_r04_sumPUPt[indexLeg4];
            neutralIso_4 = TMath::Max(float(0),neutralIso_4);
            absIso_4 += neutralIso_4;
            iso_4 = absIso_4/analysisTree.muon_pt[indexLeg4];
            if(!isData)
            {
                IdSF_4 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[indexLeg4]),double(analysisTree.muon_eta[indexLeg4]));
            }
            Leg4.SetXYZM(analysisTree.muon_px[indexLeg4],analysisTree.muon_py[indexLeg4],analysisTree.muon_pz[indexLeg4],muonMass);
            muonLoose_4 = analysisTree.muon_isLoose[indexLeg4];
            muonMedium_4 = analysisTree.muon_isMedium[indexLeg4];
            muonTight_4 = analysisTree.muon_isTight[indexLeg4];
            H_SS = false;
            if (analysisTree.electron_charge[indexLeg3]*analysisTree.muon_charge[indexLeg4]>0)
                H_SS = true;
        }
        if(HiggsFinalState == "ET")
        {
            pt_3 = analysisTree.electron_pt[indexLeg3];
            phi_3 = analysisTree.electron_phi[indexLeg3];
            eta_3 = analysisTree.electron_eta[indexLeg3];
            m_3 = classic_svFit::electronMass;
            //float absIso_3 = analysisTree.electron_r03_sumChargedHadronPt[indexLeg3];
            //float neutralIso_3 = analysisTree.electron_r03_sumNeutralHadronEt[indexLeg3] + analysisTree.electron_r03_sumPhotonEt[indexLeg3] - 0.5*analysisTree.electron_r03_sumPUPt[indexLeg3];
            //neutralIso_3 = TMath::Max(float(0),neutralIso_3);
            //absIso_3 += neutralIso_3;
            //iso_3 = absIso_3/analysisTree.electron_pt[indexLeg3];
            iso_3 = analysisTree.electron_eaIsolation[indexLeg3]/analysisTree.electron_pt[indexLeg3];
            if(!isData)
            {
                IdSF_3 = (float)SF_eleIdIso->get_ScaleFactor(double(analysisTree.electron_pt[indexLeg3]),double(analysisTree.electron_eta[indexLeg3]));
            }
            Leg3.SetXYZM(analysisTree.electron_px[indexLeg3],analysisTree.electron_py[indexLeg3],analysisTree.electron_pz[indexLeg3],electronMass);
            electronIDWP80v2_3 = analysisTree.electron_mva_wp80_noIso_Fall17_v2[indexLeg3];
            electronIDWP90v2_3 = analysisTree.electron_mva_wp90_noIso_Fall17_v2[indexLeg3];

            pt_4 = (1 + tau_Shift(era, &analysisTree, indexLeg4))*analysisTree.tau_pt[indexLeg4];
            phi_4 = analysisTree.tau_phi[indexLeg4];
            eta_4 = analysisTree.tau_eta[indexLeg4];
            dm_4 = analysisTree.tau_decayMode[indexLeg4];
            if(dm_4 == 0)
            {
                m_4 = analysisTree.tau_mass[indexLeg4];
            }
            else
            {
                m_4 = (1 + tau_Shift(era, &analysisTree, indexLeg4))*analysisTree.tau_mass[indexLeg4];
            }
            iso_4 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexLeg4];
            Leg4.SetPtEtaPhiM(analysisTree.tau_pt[indexLeg4],analysisTree.tau_eta[indexLeg4],analysisTree.tau_phi[indexLeg4],analysisTree.tau_mass[indexLeg4]);

            VSjetVVTight_4 = analysisTree.tau_byVVTightDeepTau2017v2p1VSjet[indexLeg4];
            VSjetVTight_4 = analysisTree.tau_byVTightDeepTau2017v2p1VSjet[indexLeg4];
            VSjetTight_4 = analysisTree.tau_byTightDeepTau2017v2p1VSjet[indexLeg4];
            VSjetMedium_4 = analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexLeg4];
            VSjetLoose_4 = analysisTree.tau_byLooseDeepTau2017v2p1VSjet[indexLeg4];
            VSjetVLoose_4 = analysisTree.tau_byVLooseDeepTau2017v2p1VSjet[indexLeg4];
            VSjetVVLoose_4 = analysisTree.tau_byVVLooseDeepTau2017v2p1VSjet[indexLeg4];
            VSjetVVVLoose_4 = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[indexLeg4];
            VSeVVTight_4 = analysisTree.tau_byVVTightDeepTau2017v2p1VSe[indexLeg4];
            VSeVTight_4 = analysisTree.tau_byVTightDeepTau2017v2p1VSe[indexLeg4];
            VSeTight_4 = analysisTree.tau_byTightDeepTau2017v2p1VSe[indexLeg4];
            VSeMedium_4 = analysisTree.tau_byMediumDeepTau2017v2p1VSe[indexLeg4];
            VSeLoose_4 = analysisTree.tau_byLooseDeepTau2017v2p1VSe[indexLeg4];
            VSeVLoose_4 = analysisTree.tau_byVLooseDeepTau2017v2p1VSe[indexLeg4];
            VSeVVLoose_4 = analysisTree.tau_byVVLooseDeepTau2017v2p1VSe[indexLeg4];
            VSeVVVLoose_4 = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSe[indexLeg4];
            VSmuTight_4 = analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexLeg4];
            VSmuMedium_4 = analysisTree.tau_byMediumDeepTau2017v2p1VSmu[indexLeg4];
            VSmuLoose_4 = analysisTree.tau_byLooseDeepTau2017v2p1VSmu[indexLeg4];
            VSmuVLoose_4 = analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[indexLeg4];
            
            H_SS = false;
            if (analysisTree.electron_charge[indexLeg3]*analysisTree.tau_charge[indexLeg4]>0)
                H_SS = true;
        }
        
        if(HiggsFinalState == "E")
        {
            pt_3 = analysisTree.electron_pt[indexLeg3];
            phi_3 = analysisTree.electron_phi[indexLeg3];
            eta_3 = analysisTree.electron_eta[indexLeg3];
            m_3 = classic_svFit::electronMass;

            //float absIso_3 = analysisTree.electron_r03_sumChargedHadronPt[indexLeg3];
            //float neutralIso_3 = analysisTree.electron_r03_sumNeutralHadronEt[indexLeg3] + analysisTree.electron_r03_sumPhotonEt[indexLeg3] - 0.5*analysisTree.electron_r03_sumPUPt[indexLeg3];
            //neutralIso_3 = TMath::Max(float(0),neutralIso_3);
            //absIso_3 += neutralIso_3;
            iso_3 = analysisTree.electron_eaIsolation[indexLeg3]/analysisTree.electron_pt[indexLeg3];
            if(!isData)
            {
                IdSF_3 = (float)SF_eleIdIso->get_ScaleFactor(double(analysisTree.electron_pt[indexLeg3]),double(analysisTree.electron_eta[indexLeg3]));
            }
            Leg3.SetXYZM(analysisTree.electron_px[indexLeg3],analysisTree.electron_py[indexLeg3],analysisTree.electron_pz[indexLeg3],electronMass);
            electronIDWP80v2_3 = analysisTree.electron_mva_wp80_noIso_Fall17_v2[indexLeg3];
            electronIDWP90v2_3 = analysisTree.electron_mva_wp90_noIso_Fall17_v2[indexLeg3];
            
            pt_4 = analysisTree.electron_pt[indexLeg4];
            phi_4 = analysisTree.electron_phi[indexLeg4];
            eta_4 = analysisTree.electron_eta[indexLeg4];
            m_4 = classic_svFit::electronMass;

            //float absIso_4 = analysisTree.electron_r03_sumChargedHadronPt[indexLeg4];
            //float neutralIso_4 = analysisTree.electron_r03_sumNeutralHadronEt[indexLeg4] + analysisTree.electron_r03_sumPhotonEt[indexLeg4] - 0.5*analysisTree.electron_r03_sumPUPt[indexLeg4];
            //neutralIso_4 = TMath::Max(float(0),neutralIso_4);
            //absIso_4 += neutralIso_4;
            iso_4 = analysisTree.electron_eaIsolation[indexLeg4]/analysisTree.electron_pt[indexLeg4];
            if(!isData)
            {
                IdSF_4 = (float)SF_eleIdIso->get_ScaleFactor(double(analysisTree.electron_pt[indexLeg4]),double(analysisTree.electron_eta[indexLeg4]));
            }
            Leg4.SetXYZM(analysisTree.electron_px[indexLeg4],analysisTree.electron_py[indexLeg4],analysisTree.electron_pz[indexLeg4],electronMass);
            electronIDWP80v2_4 = analysisTree.electron_mva_wp80_noIso_Fall17_v2[indexLeg4];
            electronIDWP90v2_4 = analysisTree.electron_mva_wp90_noIso_Fall17_v2[indexLeg4];
        }
        
        if(HiggsFinalState == "MT")
        {
            pt_3 = analysisTree.muon_pt[indexLeg3];
            phi_3 = analysisTree.muon_phi[indexLeg3];
            eta_3 = analysisTree.muon_eta[indexLeg3];
            m_3 = classic_svFit::muonMass;

            float absIso_3 = analysisTree.muon_r04_sumChargedHadronPt[indexLeg3];
            float neutralIso_3 = analysisTree.muon_r04_sumNeutralHadronEt[indexLeg3] + analysisTree.muon_r04_sumPhotonEt[indexLeg3] - 0.5*analysisTree.muon_r04_sumPUPt[indexLeg3];
            neutralIso_3 = TMath::Max(float(0),neutralIso_3);
            absIso_3 += neutralIso_3;
            iso_3 = absIso_3/analysisTree.muon_pt[indexLeg3];
            if(!isData)
            {
                IdSF_3 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[indexLeg3]),double(analysisTree.muon_eta[indexLeg3]));
            }
            Leg3.SetXYZM(analysisTree.muon_px[indexLeg3],analysisTree.muon_py[indexLeg3],analysisTree.muon_pz[indexLeg3],muonMass);
            muonLoose_3 = analysisTree.muon_isLoose[indexLeg3];
            muonMedium_3 = analysisTree.muon_isMedium[indexLeg3];
            muonTight_3 = analysisTree.muon_isTight[indexLeg3];
            
            pt_4 = (1 + tau_Shift(era, &analysisTree, indexLeg4))*analysisTree.tau_pt[indexLeg4];
            phi_4 = analysisTree.tau_phi[indexLeg4];
            eta_4 = analysisTree.tau_eta[indexLeg4];
            dm_4 = analysisTree.tau_decayMode[indexLeg4];
            if(dm_4 == 0)
            {
                m_4 = analysisTree.tau_mass[indexLeg4];
            }
            else
            {
                m_4 = (1 + tau_Shift(era, &analysisTree, indexLeg4))*analysisTree.tau_mass[indexLeg4];
            }
            iso_4 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexLeg4];
            Leg4.SetPtEtaPhiM(analysisTree.tau_pt[indexLeg4],analysisTree.tau_eta[indexLeg4],analysisTree.tau_phi[indexLeg4],analysisTree.tau_mass[indexLeg4]);

            VSjetVVTight_4 = analysisTree.tau_byVVTightDeepTau2017v2p1VSjet[indexLeg4];
            VSjetVTight_4 = analysisTree.tau_byVTightDeepTau2017v2p1VSjet[indexLeg4];
            VSjetTight_4 = analysisTree.tau_byTightDeepTau2017v2p1VSjet[indexLeg4];
            VSjetMedium_4 = analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexLeg4];
            VSjetLoose_4 = analysisTree.tau_byLooseDeepTau2017v2p1VSjet[indexLeg4];
            VSjetVLoose_4 = analysisTree.tau_byVLooseDeepTau2017v2p1VSjet[indexLeg4];
            VSjetVVLoose_4 = analysisTree.tau_byVVLooseDeepTau2017v2p1VSjet[indexLeg4];
            VSjetVVVLoose_4 = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[indexLeg4];
            VSeVVTight_4 = analysisTree.tau_byVVTightDeepTau2017v2p1VSe[indexLeg4];
            VSeVTight_4 = analysisTree.tau_byVTightDeepTau2017v2p1VSe[indexLeg4];
            VSeTight_4 = analysisTree.tau_byTightDeepTau2017v2p1VSe[indexLeg4];
            VSeMedium_4 = analysisTree.tau_byMediumDeepTau2017v2p1VSe[indexLeg4];
            VSeLoose_4 = analysisTree.tau_byLooseDeepTau2017v2p1VSe[indexLeg4];
            VSeVLoose_4 = analysisTree.tau_byVLooseDeepTau2017v2p1VSe[indexLeg4];
            VSeVVLoose_4 = analysisTree.tau_byVVLooseDeepTau2017v2p1VSe[indexLeg4];
            VSeVVVLoose_4 = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSe[indexLeg4];
            VSmuTight_4 = analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexLeg4];
            VSmuMedium_4 = analysisTree.tau_byMediumDeepTau2017v2p1VSmu[indexLeg4];
            VSmuLoose_4 = analysisTree.tau_byLooseDeepTau2017v2p1VSmu[indexLeg4];
            VSmuVLoose_4 = analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[indexLeg4];
            H_SS = false;
            if (analysisTree.muon_charge[indexLeg3]*analysisTree.tau_charge[indexLeg4]>0)
                H_SS = true;
            //std::cout << "H_SS:  " << H_SS << std::endl;
        }
        if(HiggsFinalState == "M")
        {
            pt_3 = analysisTree.muon_pt[indexLeg3];
            phi_3 = analysisTree.muon_phi[indexLeg3];
            eta_3 = analysisTree.muon_eta[indexLeg3];
            m_3 = classic_svFit::muonMass;

            float absIso_3 = analysisTree.muon_r04_sumChargedHadronPt[indexLeg3];
            float neutralIso_3 = analysisTree.muon_r04_sumNeutralHadronEt[indexLeg3] + analysisTree.muon_r04_sumPhotonEt[indexLeg3] - 0.5*analysisTree.muon_r04_sumPUPt[indexLeg3];
            neutralIso_3 = TMath::Max(float(0),neutralIso_3);
            absIso_3 += neutralIso_3;
            iso_3 = absIso_3/analysisTree.muon_pt[indexLeg3];
            if(!isData)
            {
                IdSF_3 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[indexLeg3]),double(analysisTree.muon_eta[indexLeg3]));
            }
            Leg3.SetXYZM(analysisTree.muon_px[indexLeg3],analysisTree.muon_py[indexLeg3],analysisTree.muon_pz[indexLeg3],muonMass);
            muonLoose_3 = analysisTree.muon_isLoose[indexLeg3];
            muonMedium_3 = analysisTree.muon_isMedium[indexLeg3];
            muonTight_3 = analysisTree.muon_isTight[indexLeg3];
            
            pt_4 = analysisTree.muon_pt[indexLeg4];
            phi_4 = analysisTree.muon_phi[indexLeg4];
            eta_4 = analysisTree.muon_eta[indexLeg4];
            m_4 = classic_svFit::muonMass;

            float absIso_4 = analysisTree.muon_r04_sumChargedHadronPt[indexLeg4];
            float neutralIso_4 = analysisTree.muon_r04_sumNeutralHadronEt[indexLeg4] + analysisTree.muon_r04_sumPhotonEt[indexLeg4] - 0.5*analysisTree.muon_r04_sumPUPt[indexLeg4];
            neutralIso_4 = TMath::Max(float(0),neutralIso_4);
            absIso_4 += neutralIso_4;
            iso_4 = absIso_4/analysisTree.muon_pt[indexLeg4];
            if(!isData)
            {
                IdSF_4 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[indexLeg4]),double(analysisTree.muon_eta[indexLeg4]));
            }
            Leg4.SetXYZM(analysisTree.muon_px[indexLeg4],analysisTree.muon_py[indexLeg4],analysisTree.muon_pz[indexLeg4],muonMass);
            muonLoose_4 = analysisTree.muon_isLoose[indexLeg4];
            muonMedium_4 = analysisTree.muon_isMedium[indexLeg4];
            muonTight_4 = analysisTree.muon_isTight[indexLeg4];
            
            H_SS = false;
            if (analysisTree.muon_charge[indexLeg3]*analysisTree.muon_charge[indexLeg4]>0)
                H_SS = true;
            //std::cout << "H_SS:  " << H_SS << std::endl;
        }
        
        if((HiggsFinalState == "TT") ||(HiggsFinalState == "T"))
        {
            pt_3 = (1 + tau_Shift(era, &analysisTree, indexLeg3))*analysisTree.tau_pt[indexLeg3];
            phi_3 = analysisTree.tau_phi[indexLeg3];
            eta_3 = analysisTree.tau_eta[indexLeg3];
            dm_3 = analysisTree.tau_decayMode[indexLeg3];
            if(dm_3 == 0)
            {
                m_3 = analysisTree.tau_mass[indexLeg3];
            }
            else
            {
                m_3 = (1 + tau_Shift(era, &analysisTree, indexLeg3))*analysisTree.tau_mass[indexLeg3];
            }
            iso_3 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexLeg3];
            Leg3.SetPtEtaPhiM(analysisTree.tau_pt[indexLeg3],analysisTree.tau_eta[indexLeg3],analysisTree.tau_phi[indexLeg3],analysisTree.tau_mass[indexLeg3]);

            VSjetVVTight_3 = analysisTree.tau_byVVTightDeepTau2017v2p1VSjet[indexLeg3];
            VSjetVTight_3 = analysisTree.tau_byVTightDeepTau2017v2p1VSjet[indexLeg3];
            VSjetTight_3 = analysisTree.tau_byTightDeepTau2017v2p1VSjet[indexLeg3];
            VSjetMedium_3 = analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexLeg3];
            VSjetLoose_3 = analysisTree.tau_byLooseDeepTau2017v2p1VSjet[indexLeg3];
            VSjetVLoose_3 = analysisTree.tau_byVLooseDeepTau2017v2p1VSjet[indexLeg3];
            VSjetVVLoose_3 = analysisTree.tau_byVVLooseDeepTau2017v2p1VSjet[indexLeg3];
            VSjetVVVLoose_3 = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[indexLeg3];
            VSeVVTight_3 = analysisTree.tau_byVVTightDeepTau2017v2p1VSe[indexLeg3];
            VSeVTight_3 = analysisTree.tau_byVTightDeepTau2017v2p1VSe[indexLeg3];
            VSeTight_3 = analysisTree.tau_byTightDeepTau2017v2p1VSe[indexLeg3];
            VSeMedium_3 = analysisTree.tau_byMediumDeepTau2017v2p1VSe[indexLeg3];
            VSeLoose_3 = analysisTree.tau_byLooseDeepTau2017v2p1VSe[indexLeg3];
            VSeVLoose_3 = analysisTree.tau_byVLooseDeepTau2017v2p1VSe[indexLeg3];
            VSeVVLoose_3 = analysisTree.tau_byVVLooseDeepTau2017v2p1VSe[indexLeg3];
            VSeVVVLoose_3 = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSe[indexLeg3];
            VSmuTight_3 = analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexLeg3];
            VSmuMedium_3 = analysisTree.tau_byMediumDeepTau2017v2p1VSmu[indexLeg3];
            VSmuLoose_3 = analysisTree.tau_byLooseDeepTau2017v2p1VSmu[indexLeg3];
            VSmuVLoose_3 = analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[indexLeg3];
            
            pt_4 = (1 + tau_Shift(era, &analysisTree, indexLeg4))*analysisTree.tau_pt[indexLeg4];
            phi_4 = analysisTree.tau_phi[indexLeg4];
            eta_4 = analysisTree.tau_eta[indexLeg4];
            dm_4 = analysisTree.tau_decayMode[indexLeg4];
            if(dm_4 == 0)
            {
                m_4 = analysisTree.tau_mass[indexLeg4];
            }
            else
            {
                m_4 = (1 + tau_Shift(era, &analysisTree, indexLeg4))*analysisTree.tau_mass[indexLeg4];
            }
            iso_4 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexLeg4];
            Leg4.SetPtEtaPhiM(analysisTree.tau_pt[indexLeg4],analysisTree.tau_eta[indexLeg4],analysisTree.tau_phi[indexLeg4],analysisTree.tau_mass[indexLeg4]);

            VSjetVVTight_4 = analysisTree.tau_byVVTightDeepTau2017v2p1VSjet[indexLeg4];
            VSjetVTight_4 = analysisTree.tau_byVTightDeepTau2017v2p1VSjet[indexLeg4];
            VSjetTight_4 = analysisTree.tau_byTightDeepTau2017v2p1VSjet[indexLeg4];
            VSjetMedium_4 = analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexLeg4];
            VSjetLoose_4 = analysisTree.tau_byLooseDeepTau2017v2p1VSjet[indexLeg4];
            VSjetVLoose_4 = analysisTree.tau_byVLooseDeepTau2017v2p1VSjet[indexLeg4];
            VSjetVVLoose_4 = analysisTree.tau_byVVLooseDeepTau2017v2p1VSjet[indexLeg4];
            VSjetVVVLoose_4 = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[indexLeg4];
            VSeVVTight_4 = analysisTree.tau_byVVTightDeepTau2017v2p1VSe[indexLeg4];
            VSeVTight_4 = analysisTree.tau_byVTightDeepTau2017v2p1VSe[indexLeg4];
            VSeTight_4 = analysisTree.tau_byTightDeepTau2017v2p1VSe[indexLeg4];
            VSeMedium_4 = analysisTree.tau_byMediumDeepTau2017v2p1VSe[indexLeg4];
            VSeLoose_4 = analysisTree.tau_byLooseDeepTau2017v2p1VSe[indexLeg4];
            VSeVLoose_4 = analysisTree.tau_byVLooseDeepTau2017v2p1VSe[indexLeg4];
            VSeVVLoose_4 = analysisTree.tau_byVVLooseDeepTau2017v2p1VSe[indexLeg4];
            VSeVVVLoose_4 = analysisTree.tau_byVVVLooseDeepTau2017v2p1VSe[indexLeg4];
            VSmuTight_4 = analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexLeg4];
            VSmuMedium_4 = analysisTree.tau_byMediumDeepTau2017v2p1VSmu[indexLeg4];
            VSmuLoose_4 = analysisTree.tau_byLooseDeepTau2017v2p1VSmu[indexLeg4];
            VSmuVLoose_4 = analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[indexLeg4];
            
            H_SS = false;
            if (analysisTree.tau_charge[indexLeg3]*analysisTree.tau_charge[indexLeg4]>0)
                H_SS = true;
            //std::cout << "H_SS:  " << H_SS << std::endl;
        }
        
        if((!isData) || (HiggsFinalState == "E") || (HiggsFinalState == "M") || (HiggsFinalState == "T"))//filter out by signal region isolation selection for MC
        {
            if(iso_1>0.15) continue;
            if(iso_2>0.15) continue;
            if(HiggsFinalState == "EM")
            {
                if(iso_3>0.15) continue;
                //if(iso_ea_3>0.15) continue;
                if(iso_4>0.15) continue;

            }
            if(HiggsFinalState == "ET")
            {
                if(iso_3>0.15) continue;
                //if(iso_ea_3>0.15) continue;
            }
            if(HiggsFinalState == "MT")
            {
                if(iso_3>0.15) continue;
            }
        }
        if(isDebug)
        {
            std::cout << "Passing isolation imposed" << std::endl;
        }
        
        
        Leg3_unCorrected_1 = Leg3; // for uncluster MET Up
        Leg4_unCorrected_1 = Leg4;
        Leg3_unCorrected_2 = Leg3; // for uncluster MET Down
        Leg4_unCorrected_2 = Leg4;
        Leg3_unCorrected_3 = Leg3; // for JEC MET Up
        Leg4_unCorrected_3 = Leg4;
        Leg3_unCorrected_4 = Leg3; // for JEC MET Down
        Leg4_unCorrected_4 = Leg4;
        Leg3_unCorrected_5 = Leg3; // for HEM 2018 JEC MET Up
        Leg4_unCorrected_5 = Leg4;
        Leg3_unCorrected_6 = Leg3; // for HEM 2018 JEC MET Down
        Leg4_unCorrected_6 = Leg4;
        //new MC matching piece of code
        gen_match_1 = 6;
        gen_match_2 = 6;
        gen_match_3 = 6;
        gen_match_4 = 6;
        float minDR_1 = 0.2;
        float minDR_2 = 0.2;
        float minDR_3 = 0.2;
        float minDR_4 = 0.2;
        if (!isData)
        {
            for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen)
            {
                float ptGen = PtoPt(analysisTree.genparticles_px[igen], analysisTree.genparticles_py[igen]);
                bool type1 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
                bool type2 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] && ptGen>8;
                bool type3 = abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
                bool type4 = abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] && ptGen>8;
                bool isAnyType = type1 || type2 || type3 || type4;
                if (isAnyType)
                {
                    float etaGen = PtoEta(analysisTree.genparticles_px[igen],analysisTree.genparticles_py[igen], analysisTree.genparticles_pz[igen]);
                    float phiGen = PtoPhi(analysisTree.genparticles_px[igen],analysisTree.genparticles_py[igen]);
                    float deltaR_1 = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],etaGen,phiGen);
                    if (deltaR_1<minDR_1)
                    {
                        minDR_1 = deltaR_1;
                        if (type1) gen_match_1 = 1;
                        else if (type2) gen_match_1 = 2;
                        else if (type3) gen_match_1 = 3;
                        else if (type4) gen_match_1 = 4;
                    }
                    float deltaR_2 = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],etaGen,phiGen);
                    if (deltaR_2<minDR_2)
                    {
                        minDR_2 = deltaR_2;
                        if (type1) gen_match_2 = 1;
                        else if (type2) gen_match_2 = 2;
                        else if (type3) gen_match_2 = 3;
                        else if (type4) gen_match_2 = 4;
                    }

                    float deltaR_3 = deltaR(Leg3.Eta(),Leg3.Phi(),etaGen,phiGen);
                    if (deltaR_3<minDR_3)
                    {
                        minDR_3 = deltaR_3;
                        if (type1) gen_match_3 = 1;
                        else if (type2) gen_match_3 = 2;
                        else if (type3) gen_match_3 = 3;
                        else if (type4) gen_match_3 = 4;
                    }
                    float deltaR_4 = deltaR(Leg4.Eta(),Leg4.Phi(),etaGen,phiGen);
                    if (deltaR_4<minDR_4)
                    {
                        minDR_4 = deltaR_4;
                        if (type1) gen_match_4 = 1;
                        else if (type2) gen_match_4 = 2;
                        else if (type3) gen_match_4 = 3;
                        else if (type4) gen_match_4 = 4;
                    }
                }
            }
            for (unsigned int igen=0; igen < analysisTree.gentau_count; ++igen)
            {
                if (analysisTree.gentau_visibleNoLep_pt[igen]>15.)
                {
                    float deltaR_1 = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.gentau_visibleNoLep_eta[igen],analysisTree.gentau_visibleNoLep_phi[igen]);
                    if (deltaR_1<minDR_1)
                    {
                        minDR_1 = deltaR_1;
                        gen_match_1 = 5;
                    }
                    float deltaR_2 = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.gentau_visibleNoLep_eta[igen],analysisTree.gentau_visibleNoLep_phi[igen]);
                    if (deltaR_2<minDR_2)
                    {
                        minDR_2 = deltaR_2;
                        gen_match_2 = 5;
                    }

                    float deltaR_3 = deltaR(Leg3.Eta(),Leg3.Phi(),analysisTree.gentau_visibleNoLep_eta[igen],analysisTree.gentau_visibleNoLep_phi[igen]);
                    if (deltaR_3<minDR_3)
                    {
                        minDR_3 = deltaR_3;
                        gen_match_3 = 5;
                    }
                    float deltaR_4 = deltaR(Leg4.Eta(),Leg4.Phi(),analysisTree.gentau_visibleNoLep_eta[igen],analysisTree.gentau_visibleNoLep_phi[igen]);
                    if (deltaR_4<minDR_4)
                    {
                        minDR_4 = deltaR_4;
                        gen_match_4 = 5;
                    }
                }
            }
        }
        tauVSjetSF_3 = 1;
        tauVSjetSFUp_3 = 1;
        tauVSjetSFDown_3 = 1;
        tauVSeSF_3 = 1;
        tauVSeSFUp_3 = 1;
        tauVSeSFDown_3 = 1;
        tauVSmuSF_3 = 1;
        tauVSmuSFUp_3 = 1;
        tauVSmuSFDown_3 = 1;
        tauVSjetSF_4 = 1;
        tauVSjetSFUp_4 = 1;
        tauVSjetSFDown_4 = 1;
        tauVSeSF_4 = 1;
        tauVSeSFUp_4 = 1;
        tauVSeSFDown_4 = 1;
        tauVSmuSF_4 = 1;
        tauVSmuSFUp_4 = 1;
        tauVSmuSFDown_4 = 1;
        tauVSjetVVVLooseSF_3 = 1;
        
        if(HiggsFinalState == "MT")
        {
            if (gen_match_4 == 5)
            {
                tauVSjetSF_4 = tauIDSF_VSjetMedium->getSFvsPT(pt_4);
                tauVSjetSFUp_4 = tauIDSF_VSjetMedium->getSFvsPT(pt_4,"Up");
                tauVSjetSFDown_4 = tauIDSF_VSjetMedium->getSFvsPT(pt_4,"Down");
            }
            if ((gen_match_4 == 1)||(gen_match_4 == 3))
            {
                tauVSeSF_4 = tauIDSF_VSeVLoose->getSFvsEta(eta_4,gen_match_4);
                tauVSeSFUp_4 = tauIDSF_VSeVLoose->getSFvsEta(eta_4,gen_match_4,"Up");
                tauVSeSFDown_4 = tauIDSF_VSeVLoose->getSFvsEta(eta_4,gen_match_4,"Down");
            }
            if ((gen_match_4 == 2)||(gen_match_4 == 4))
            {
                tauVSmuSF_4 = tauIDSF_VSmuTight->getSFvsEta(eta_4,gen_match_4);
                tauVSmuSFUp_4 = tauIDSF_VSmuTight->getSFvsEta(eta_4,gen_match_4,"Up");
                tauVSmuSFDown_4 = tauIDSF_VSmuTight->getSFvsEta(eta_4,gen_match_4,"Down");
            }
        }
        if(HiggsFinalState == "ET")
        {
            if (gen_match_4 == 5)
            {
                tauVSjetSF_4 = tauIDSF_VSjetMedium->getSFvsPT(pt_4);
                tauVSjetSFUp_4 = tauIDSF_VSjetMedium->getSFvsPT(pt_4,"Up");
                tauVSjetSFDown_4 = tauIDSF_VSjetMedium->getSFvsPT(pt_4,"Down");
            }
            if ((gen_match_4 == 1)||(gen_match_4 == 3))
            {
                tauVSeSF_4 = tauIDSF_VSeTight->getSFvsEta(eta_4,gen_match_4);
                tauVSeSFUp_4 = tauIDSF_VSeTight->getSFvsEta(eta_4,gen_match_4,"Up");
                tauVSeSFDown_4 = tauIDSF_VSeTight->getSFvsEta(eta_4,gen_match_4,"Down");
            }
            if ((gen_match_4 == 2)||(gen_match_4 == 4))
            {
                tauVSmuSF_4 = tauIDSF_VSmuVLoose->getSFvsEta(eta_4,gen_match_4);
                tauVSmuSFUp_4 = tauIDSF_VSmuVLoose->getSFvsEta(eta_4,gen_match_4,"Up");
                tauVSmuSFDown_4 = tauIDSF_VSmuVLoose->getSFvsEta(eta_4,gen_match_4,"Down");
            }
        }
        
        if(HiggsFinalState == "TT")
        {
            if (gen_match_3 == 5)
            {
                tauVSjetSF_3 = tauIDSF_VSjetMedium->getSFvsPT(pt_3);
                tauVSjetSFUp_3 = tauIDSF_VSjetMedium->getSFvsPT(pt_3,"Up");
                tauVSjetSFDown_3 = tauIDSF_VSjetMedium->getSFvsPT(pt_3,"Down");
            }
            if ((gen_match_3 == 1)||(gen_match_3 == 3))
            {
                tauVSeSF_3 = tauIDSF_VSeVLoose->getSFvsEta(eta_3,gen_match_3);
                tauVSeSFUp_3 = tauIDSF_VSeVLoose->getSFvsEta(eta_3,gen_match_3,"Up");
                tauVSeSFDown_3 = tauIDSF_VSeVLoose->getSFvsEta(eta_3,gen_match_3,"Down");
            }
            if ((gen_match_3 == 2)||(gen_match_3 == 4))
            {
                tauVSmuSF_3 = tauIDSF_VSmuVLoose->getSFvsEta(eta_3,gen_match_3);
                tauVSmuSFUp_3 = tauIDSF_VSmuVLoose->getSFvsEta(eta_3,gen_match_3,"Up");
                tauVSmuSFDown_3 = tauIDSF_VSmuVLoose->getSFvsEta(eta_3,gen_match_3,"Down");
            }
            if (gen_match_4 == 5)
            {
                tauVSjetSF_4 = tauIDSF_VSjetMedium->getSFvsPT(pt_4);
                tauVSjetSFUp_4 = tauIDSF_VSjetMedium->getSFvsPT(pt_4,"Up");
                tauVSjetSFDown_4 = tauIDSF_VSjetMedium->getSFvsPT(pt_4,"Down");
            }
            if ((gen_match_4 == 1)||(gen_match_4 == 3))
            {
                tauVSeSF_4 = tauIDSF_VSeVLoose->getSFvsEta(eta_4,gen_match_4);
                tauVSeSFUp_4 = tauIDSF_VSeVLoose->getSFvsEta(eta_4,gen_match_4,"Up");
                tauVSeSFDown_4 = tauIDSF_VSeVLoose->getSFvsEta(eta_4,gen_match_4,"Down");
            }
            if ((gen_match_4 == 2)||(gen_match_4 == 4))
            {
                tauVSmuSF_4 = tauIDSF_VSmuVLoose->getSFvsEta(eta_4,gen_match_4);
                tauVSmuSFUp_4 = tauIDSF_VSmuVLoose->getSFvsEta(eta_4,gen_match_4,"Up");
                tauVSmuSFDown_4 = tauIDSF_VSmuVLoose->getSFvsEta(eta_4,gen_match_4,"Down");
            }
        }
        
        if(HiggsFinalState == "T")
        {
            if (gen_match_3 == 5)
            {
                tauVSjetSF_3 = tauIDSF_VSjetMedium->getSFvsPT(pt_3);
                tauVSjetSFUp_3 = tauIDSF_VSjetMedium->getSFvsPT(pt_3,"Up");
                tauVSjetSFDown_3 = tauIDSF_VSjetMedium->getSFvsPT(pt_3,"Down");
            }
            if ((gen_match_3 == 1)||(gen_match_3 == 3))
                tauVSeSF_3 = tauIDSF_VSeVLoose->getSFvsEta(eta_3,gen_match_3);
            if ((gen_match_3 == 2)||(gen_match_3 == 4))
                tauVSmuSF_3 = tauIDSF_VSmuVLoose->getSFvsEta(eta_3,gen_match_3);
            if (gen_match_3 == 5)
                tauVSjetVVVLooseSF_3 = tauIDSF_VSjetVVVLoose->getSFvsPT(pt_3);
        }

        //uncorrected met
        met = TMath::Sqrt(analysisTree.pfmetcorr_ex*analysisTree.pfmetcorr_ex + analysisTree.pfmetcorr_ey*analysisTree.pfmetcorr_ey);
        metphi = analysisTree.pfmetcorr_phi;
        
        // MET Filters ================================================================================================================================================
        passMETFilters = metFiltersPasses(analysisTree,metFlags,isData);
        
        ////////////////////////////////////////////////////////////
        // Zpt reweighting
        ////////////////////////////////////////////////////////////
        TLorentzVector genV( 0., 0., 0., 0.);
        TLorentzVector genL( 0., 0., 0., 0.);
        UInt_t njetshad = analysisTree.pfjet_count;
        zptweight = 1.;
        
        if (!isData && isDY ){
            genV = genTools::genV(analysisTree);
            genL = genTools::genL(analysisTree);
            
            if(applyRecoilCorrections)
            {
                float met_rcmr = met;
                float metphi_rcmr = metphi;
                // PF MET
                genTools::RecoilCorrections( *recoilPFMetCorrector,
                                            (!isData && applyRecoilCorrections && isDY ) * genTools::MeanResolution,
                                            met, metphi,
                                            genV.Px(), genV.Py(),
                                            genL.Px(), genL.Py(),
                                            njetshad,
                                            met_rcmr, metphi_rcmr
                                            );
                met = met_rcmr;
                metphi = metphi_rcmr;
            }
            
            float bosonMass = genV.M();
            float bosonPt = genV.Pt();
            Float_t zptmassweight = 1;
            if (bosonMass>50.0) {
                float maxmass = h_zptweight->GetXaxis()->GetXmax();
                float maxpt = h_zptweight->GetYaxis()->GetXmax();
                float minmass = h_zptweight->GetXaxis()->GetXmin();
                float minpt = h_zptweight->GetYaxis()->GetXmin();
                float _bosonMass = bosonMass;
                float _bosonPt = bosonPt;
                if(bosonMass>maxmass) _bosonMass=maxmass;
                else if(bosonMass<minmass) _bosonMass=minmass;
                if(bosonPt>maxpt) _bosonPt = maxpt;
                else if(bosonPt<minpt) _bosonPt = minpt;
                zptmassweight = h_zptweight->GetBinContent(h_zptweight->GetXaxis()->FindBin(_bosonMass),
                                                           h_zptweight->GetYaxis()->FindBin(_bosonPt));
            }
            zptweight =zptmassweight;
        }
        
        metcov00 = analysisTree.pfmetcorr_sigxx;
        metcov01 = analysisTree.pfmetcorr_sigxy;
        metcov10 = analysisTree.pfmetcorr_sigyx;
        metcov11 = analysisTree.pfmetcorr_sigyy;
        
        //puppi met
        if(ApplyPuppiMet)
        {
            met = TMath::Sqrt(analysisTree.puppimet_ex*analysisTree.puppimet_ex + analysisTree.puppimet_ey*analysisTree.puppimet_ey);
            metphi = analysisTree.puppimet_phi;
            metcov00 = analysisTree.puppimet_sigxx;
            metcov01 = analysisTree.puppimet_sigxy;
            metcov10 = analysisTree.puppimet_sigyx;
            metcov11 = analysisTree.puppimet_sigyy;
        }

        TLorentzVector metLV;
        metLV.SetXYZT(met*TMath::Cos(metphi),met*TMath::Sin(metphi),0,TMath::Sqrt(met*TMath::Sin(metphi)*met*TMath::Sin(metphi) + met*TMath::Cos(metphi)*met*TMath::Cos(metphi)));
        
        //Corrected MET after TES shift
        if(!isData)
        {
            if((HiggsFinalState=="ET")||(HiggsFinalState=="MT"))
            {
                correctTauES(Leg4, metLV, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
            }
            if(HiggsFinalState=="TT")
            {
                correctTauES(Leg3, metLV, tau_Shift(era, &analysisTree, indexLeg3), dm_3);
                correctTauES(Leg4, metLV, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
            }
        }
        
        met = TMath::Sqrt(metLV.X() * metLV.X() + metLV.Y() * metLV.Y());
        metphi = metLV.Phi();
        
        TLorentzVector muon1;
        muon1.SetXYZM(analysisTree.muon_px[indexMuon1],analysisTree.muon_py[indexMuon1],analysisTree.muon_pz[indexMuon1],muonMass);
        pfmt_1=mT(muon1,metLV);
        
        TLorentzVector muon2;
        muon2.SetXYZM(analysisTree.muon_px[indexMuon2],analysisTree.muon_py[indexMuon2],analysisTree.muon_pz[indexMuon2],muonMass);
        pfmt_2=mT(muon2,metLV);
        
        pfmt_3=mT(Leg3,metLV);
        pfmt_4=mT(Leg4,metLV);

        m_vis = (Leg3+Leg4).M();
        pt_tt = (Leg3+Leg4).Pt();

        mll = (muon1+muon2).M();
        Z_Pt = (muon1+muon2).Pt();
        
        if ( mll < 60 || mll > 120)
            continue;
        
        if(isDebug)
        {
            std::cout << "Passing mll cut on Z candidat" << std::endl;
        }
        
        //MET systematics shifts
        float metx_mscaleUp = metLV.X();
        float mety_mscaleUp = metLV.Y();
        float metx_mscaleDown = metLV.X();
        float mety_mscaleDown = metLV.Y();
        float metx_escaleUp = metLV.X();
        float mety_escaleUp = metLV.Y();
        float metx_escaleDown = metLV.X();
        float mety_escaleDown = metLV.Y();
        float metx_tscaleUp = metLV.X();
        float mety_tscaleUp = metLV.Y();
        float metx_tscaleDown = metLV.X();
        float mety_tscaleDown = metLV.Y();
        float metx_eftscaleUp = metLV.X();
        float mety_eftscaleUp = metLV.Y();
        float metx_eftscaleDown = metLV.X();
        float mety_eftscaleDown = metLV.Y();
        float metx_mftscaleUp = metLV.X();
        float mety_mftscaleUp = metLV.Y();
        float metx_mftscaleDown = metLV.X();
        float mety_mftscaleDown = metLV.Y();
        float metx_unclmetscaleUp = metLV.X();
        float mety_unclmetscaleUp = metLV.Y();
        float metx_unclmetscaleDown = metLV.X();
        float mety_unclmetscaleDown = metLV.Y();
        float metx_jecscaleUp = metLV.X();
        float mety_jecscaleUp = metLV.Y();
        float metx_jecscaleDown = metLV.X();
        float mety_jecscaleDown = metLV.Y();
        float metx_HEM2018scaleUp = metLV.X();
        float mety_HEM2018scaleUp = metLV.Y();
        float metx_HEM2018scaleDown = metLV.X();
        float mety_HEM2018scaleDown = metLV.Y();
        
        if(applySystShift && !isData)
        {
            if(!ApplyPuppiMet)
            {
                metx_unclmetscaleUp = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
                mety_unclmetscaleUp = analysisTree.pfmetcorr_ey_UnclusteredEnUp;
                metx_unclmetscaleDown = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
                mety_unclmetscaleDown = analysisTree.pfmetcorr_ey_UnclusteredEnDown;
                metx_jecscaleUp = analysisTree.pfmetcorr_ex_JetEnUp;
                mety_jecscaleUp = analysisTree.pfmetcorr_ey_JetEnUp;
                metx_jecscaleDown = analysisTree.pfmetcorr_ex_JetEnDown;
                mety_jecscaleDown = analysisTree.pfmetcorr_ey_JetEnDown;
                
                metx_HEM2018scaleUp = analysisTree.pfmetcorr_ex;
                mety_HEM2018scaleUp = analysisTree.pfmetcorr_ey;
                metx_HEM2018scaleDown = analysisTree.pfmetcorr_ex;
                mety_HEM2018scaleDown = analysisTree.pfmetcorr_ey;
                
                //HEM issue of JEC propagating to MET
                if(era == "UL2018")
                {
                    for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet)
                    {
                        float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
                        float jetEta = analysisTree.pfjet_eta[jet];
                        
                        bool isPFJetId =false;
                        isPFJetId = tightJetiD_UL2018UL2017(analysisTree,int(jet));
                        if (!isPFJetId) continue;
                        if (analysisTree.pfjet_pt[jet] < 15) continue;
                        
                        bool cleanedJet = true;
                        float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                           eta_1,phi_1);
                        if (dR1<dRJetLeptonCut) cleanedJet = false;
                        float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                           eta_2,phi_2);
                        if (dR2<dRJetLeptonCut) cleanedJet = false;
                        if(HiggsFinalState == "EM")
                        {
                            float dR3 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                               eta_3,phi_3);
                            if (dR3<dRJetLeptonCut) cleanedJet = false;
                            float dR4 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                               eta_4,phi_4);
                            if (dR4<dRJetLeptonCut) cleanedJet = false;
                        }
                        if((HiggsFinalState == "ET") ||(HiggsFinalState == "MT"))
                        {
                            float dR3 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                               eta_3,phi_3);
                            if (dR3<dRJetLeptonCut) cleanedJet = false;
                        }
                        if (!cleanedJet) continue;

                        if (analysisTree.pfjet_phi[jet] > (-1.57) && analysisTree.pfjet_phi[jet] < (-0.87) && jetEta > (-2.5) && jetEta < (-1.3))
                        {
                            
                            metx_HEM2018scaleUp = metx_HEM2018scaleUp - analysisTree.pfjet_px[jet]*(0.2);
                            mety_HEM2018scaleUp = mety_HEM2018scaleUp - analysisTree.pfjet_py[jet]*(0.2);

                            metx_HEM2018scaleDown = metx_HEM2018scaleDown + analysisTree.pfjet_px[jet]*(0.2);
                            mety_HEM2018scaleDown = mety_HEM2018scaleDown + analysisTree.pfjet_px[jet]*(0.2);
                        }
                        
                        if (analysisTree.pfjet_phi[jet] > (-1.57) && analysisTree.pfjet_phi[jet] < (-0.87) && jetEta > (-3.0) && jetEta < (-2.5))
                        {
                            
                            metx_HEM2018scaleUp = metx_HEM2018scaleUp - analysisTree.pfjet_px[jet]*(0.35);
                            mety_HEM2018scaleUp = mety_HEM2018scaleUp - analysisTree.pfjet_py[jet]*(0.35);

                            metx_HEM2018scaleDown = metx_HEM2018scaleDown + analysisTree.pfjet_px[jet]*(0.35);
                            mety_HEM2018scaleDown = mety_HEM2018scaleDown + analysisTree.pfjet_px[jet]*(0.35);
                        }
                    }
                }
            }
            else
            {
                metx_unclmetscaleUp = analysisTree.puppimet_ex_UnclusteredEnUp;
                mety_unclmetscaleUp = analysisTree.puppimet_ey_UnclusteredEnUp;
                metx_unclmetscaleDown = analysisTree.puppimet_ex_UnclusteredEnDown;
                mety_unclmetscaleDown = analysisTree.puppimet_ey_UnclusteredEnDown;
                metx_jecscaleUp = analysisTree.puppimet_ex_JetEnUp;
                mety_jecscaleUp = analysisTree.puppimet_ey_JetEnUp;
                metx_jecscaleDown = analysisTree.puppimet_ex_JetEnDown;
                mety_jecscaleDown = analysisTree.puppimet_ey_JetEnDown;
                metx_HEM2018scaleUp = analysisTree.puppimet_ex;//Not switch-on yet
                mety_HEM2018scaleUp = analysisTree.puppimet_ey;
                metx_HEM2018scaleDown = analysisTree.puppimet_ex;
                mety_HEM2018scaleDown = analysisTree.puppimet_ey;
            }
            if(HiggsFinalState == "EM")
            {
                metx_mscaleUp = metLV.X() - (muon1.X() + muon2.X() + Leg4.X()) * MuonScale;
                mety_mscaleUp = metLV.Y() - (muon1.Y() + muon2.Y() + Leg4.Y()) * MuonScale;
                metx_mscaleDown = metLV.X() + (muon1.X() + muon2.X() + Leg4.X()) * MuonScale;
                mety_mscaleDown = metLV.Y() + (muon1.Y() + muon2.Y() + Leg4.Y()) * MuonScale;
                
                metx_escaleUp = metLV.X() - (analysisTree.electron_px_energyscale_up[indexLeg3] - analysisTree.electron_px[indexLeg3]);
                mety_escaleUp = metLV.Y() - (analysisTree.electron_py_energyscale_up[indexLeg3] - analysisTree.electron_py[indexLeg3]);
                metx_escaleDown = metLV.X() + (analysisTree.electron_px[indexLeg3] - analysisTree.electron_px_energyscale_down[indexLeg3]);
                mety_escaleDown = metLV.Y() + (analysisTree.electron_py[indexLeg3] - analysisTree.electron_py_energyscale_down[indexLeg3]);
            }
            
            if(HiggsFinalState == "ET")
            {
                metx_mscaleUp = metLV.X() - (muon1.X() + muon2.X()) * MuonScale;
                mety_mscaleUp = metLV.Y() - (muon1.Y() + muon2.Y()) * MuonScale;
                metx_mscaleDown = metLV.X() + (muon1.X() + muon2.X()) * MuonScale;
                mety_mscaleDown = metLV.Y() + (muon1.Y() + muon2.Y()) * MuonScale;
                
                metx_escaleUp = metLV.X() - (analysisTree.electron_px_energyscale_up[indexLeg3] - analysisTree.electron_px[indexLeg3]);
                mety_escaleUp = metLV.Y() - (analysisTree.electron_py_energyscale_up[indexLeg3] - analysisTree.electron_py[indexLeg3]);
                metx_escaleDown = metLV.X() + (analysisTree.electron_px[indexLeg3] - analysisTree.electron_px_energyscale_down[indexLeg3]);
                mety_escaleDown = metLV.Y() + (analysisTree.electron_py[indexLeg3] - analysisTree.electron_py_energyscale_down[indexLeg3]);
                
                metx_tscaleUp = metLV.X() - Leg4.X() * tau_ES_Syst(era,gen_match_4,dm_4);
                mety_tscaleUp = metLV.Y() - Leg4.Y() * tau_ES_Syst(era,gen_match_4,dm_4);
                metx_tscaleDown = metLV.X() + Leg4.X() * tau_ES_Syst(era,gen_match_4,dm_4);
                mety_tscaleDown = metLV.Y() + Leg4.Y() * tau_ES_Syst(era,gen_match_4,dm_4);
                
                float efes_error_up = 0.0;
                float efes_error_down = 0.0;
                tau_EFES_Syst(era,gen_match_4,dm_4,eta_4, efes_error_up, efes_error_down);
                metx_eftscaleUp = metLV.X() - Leg4.X() * efes_error_up;
                mety_eftscaleUp = metLV.Y() - Leg4.Y() * efes_error_up;
                metx_eftscaleDown = metLV.X() + Leg4.X() * efes_error_down;
                mety_eftscaleDown = metLV.Y() + Leg4.Y() * efes_error_down;
                
                metx_mftscaleUp = metLV.X() - Leg4.X() * tau_MFES_Syst(gen_match_4);
                mety_mftscaleUp = metLV.Y() - Leg4.Y() * tau_MFES_Syst(gen_match_4);
                metx_mftscaleDown = metLV.X() + Leg4.X() * tau_MFES_Syst(gen_match_4);
                mety_mftscaleDown = metLV.Y() + Leg4.Y() * tau_MFES_Syst(gen_match_4);
            }
            
            if(HiggsFinalState == "MT")
            {
                metx_mscaleUp = metLV.X() - (muon1.X() + muon2.X() + Leg3.X()) * MuonScale;
                mety_mscaleUp = metLV.Y() - (muon1.Y() + muon2.Y() + Leg3.Y()) * MuonScale;
                metx_mscaleDown = metLV.X() + (muon1.X() + muon2.X() + Leg3.X()) * MuonScale;
                mety_mscaleDown = metLV.Y() + (muon1.Y() + muon2.Y() + Leg3.Y()) * MuonScale;
                
                metx_tscaleUp = metLV.X() - Leg4.X() * tau_ES_Syst(era,gen_match_4,dm_4);
                mety_tscaleUp = metLV.Y() - Leg4.Y() * tau_ES_Syst(era,gen_match_4,dm_4);
                metx_tscaleDown = metLV.X() + Leg4.X() * tau_ES_Syst(era,gen_match_4,dm_4);
                mety_tscaleDown = metLV.Y() + Leg4.Y() * tau_ES_Syst(era,gen_match_4,dm_4);
                
                float efes_error_up = 0.0;
                float efes_error_down = 0.0;
                tau_EFES_Syst(era,gen_match_4,dm_4,eta_4, efes_error_up, efes_error_down);
                metx_eftscaleUp = metLV.X() - Leg4.X() * efes_error_up;
                mety_eftscaleUp = metLV.Y() - Leg4.Y() * efes_error_up;
                metx_eftscaleDown = metLV.X() + Leg4.X() * efes_error_down;
                mety_eftscaleDown = metLV.Y() + Leg4.Y() * efes_error_down;
                
                metx_mftscaleUp = metLV.X() - Leg4.X() * tau_MFES_Syst(gen_match_4);
                mety_mftscaleUp = metLV.Y() - Leg4.Y() * tau_MFES_Syst(gen_match_4);
                metx_mftscaleDown = metLV.X() + Leg4.X() * tau_MFES_Syst(gen_match_4);
                mety_mftscaleDown = metLV.Y() + Leg4.Y() * tau_MFES_Syst(gen_match_4);
            }
            
            if(HiggsFinalState == "TT")
            {
                metx_mscaleUp = metLV.X() - (muon1.X() + muon2.X()) * MuonScale;
                mety_mscaleUp = metLV.Y() - (muon1.Y() + muon2.Y()) * MuonScale;
                metx_mscaleDown = metLV.X() + (muon1.X() + muon2.X()) * MuonScale;
                mety_mscaleDown = metLV.Y() + (muon1.Y() + muon2.Y()) * MuonScale;
                
                metx_tscaleUp = metLV.X() - Leg3.X() * tau_ES_Syst(era,gen_match_3,dm_3) - Leg4.X() * tau_ES_Syst(era,gen_match_4,dm_4);
                mety_tscaleUp = metLV.Y() - Leg3.Y() * tau_ES_Syst(era,gen_match_3,dm_3) - Leg4.Y() * tau_ES_Syst(era,gen_match_4,dm_4);
                metx_tscaleDown = metLV.X() + Leg3.X() * tau_ES_Syst(era,gen_match_3,dm_3) + Leg4.X() * tau_ES_Syst(era,gen_match_4,dm_4);
                mety_tscaleDown = metLV.Y() + Leg3.Y() * tau_ES_Syst(era,gen_match_3,dm_3) + Leg4.Y() * tau_ES_Syst(era,gen_match_4,dm_4);
                
                float efes_error_up_3 = 0.0;
                float efes_error_down_3 = 0.0;
                tau_EFES_Syst(era,gen_match_3,dm_3,eta_3, efes_error_up_3, efes_error_down_3);
                float efes_error_up_4 = 0.0;
                float efes_error_down_4 = 0.0;
                tau_EFES_Syst(era,gen_match_4,dm_4,eta_4, efes_error_up_4, efes_error_down_4);
                
                metx_eftscaleUp = metLV.X() - Leg3.X() * efes_error_up_3 - Leg4.X() * efes_error_up_4;
                mety_eftscaleUp = metLV.Y() - Leg3.Y() * efes_error_up_3 - Leg4.Y() * efes_error_up_4;
                metx_eftscaleDown = metLV.X() + Leg3.X() * efes_error_up_3 + Leg4.X() * efes_error_down_4;
                mety_eftscaleDown = metLV.Y() + Leg3.Y() * efes_error_up_3 + Leg4.Y() * efes_error_down_4;
                
                metx_mftscaleUp = metLV.X() - Leg3.X() * tau_MFES_Syst(gen_match_3) - Leg4.X() * tau_MFES_Syst(gen_match_4);
                mety_mftscaleUp = metLV.Y() - Leg3.Y() * tau_MFES_Syst(gen_match_3) - Leg4.Y() * tau_MFES_Syst(gen_match_4);
                metx_mftscaleDown = metLV.X() + Leg3.X() * tau_MFES_Syst(gen_match_3) + Leg4.X() * tau_MFES_Syst(gen_match_4);
                mety_mftscaleDown = metLV.Y() + Leg3.Y() * tau_MFES_Syst(gen_match_3) + Leg4.Y() * tau_MFES_Syst(gen_match_4);
            }
        }
        float met_mscaleUp = TMath::Sqrt(metx_mscaleUp*metx_mscaleUp+mety_mscaleUp*mety_mscaleUp);
        float met_mscaleDown = TMath::Sqrt(metx_mscaleDown*metx_mscaleDown+mety_mscaleDown*mety_mscaleDown);
        float met_escaleUp = TMath::Sqrt(metx_escaleUp*metx_escaleUp+mety_escaleUp*mety_escaleUp);
        float met_escaleDown = TMath::Sqrt(metx_escaleDown*metx_escaleDown+mety_escaleDown*mety_escaleDown);
        float met_tscaleUp = TMath::Sqrt(metx_tscaleUp*metx_tscaleUp+mety_tscaleUp*mety_tscaleUp);
        float met_tscaleDown = TMath::Sqrt(metx_tscaleDown*metx_tscaleDown+mety_tscaleDown*mety_tscaleDown);
        float met_eftscaleUp = TMath::Sqrt(metx_eftscaleUp*metx_eftscaleUp+mety_eftscaleUp*mety_eftscaleUp);
        float met_eftscaleDown = TMath::Sqrt(metx_eftscaleDown*metx_eftscaleDown+mety_eftscaleDown*mety_eftscaleDown);
        float met_mftscaleUp = TMath::Sqrt(metx_mftscaleUp*metx_mftscaleUp+mety_mftscaleUp*mety_mftscaleUp);
        float met_mftscaleDown = TMath::Sqrt(metx_mftscaleDown*metx_mftscaleDown+mety_mftscaleDown*mety_mftscaleDown);
        float met_unclmetscaleUp = TMath::Sqrt(metx_unclmetscaleUp*metx_unclmetscaleUp+mety_unclmetscaleUp*mety_unclmetscaleUp);
        float met_unclmetscaleDown = TMath::Sqrt(metx_unclmetscaleDown*metx_unclmetscaleDown+mety_unclmetscaleDown*mety_unclmetscaleDown);
        float met_jecscaleUp = TMath::Sqrt(metx_jecscaleUp*metx_jecscaleUp+mety_jecscaleUp*mety_jecscaleUp);
        float met_jecscaleDown = TMath::Sqrt(metx_jecscaleDown*metx_jecscaleDown+mety_jecscaleDown*mety_jecscaleDown);
        float met_HEM2018scaleUp = TMath::Sqrt(metx_HEM2018scaleUp*metx_HEM2018scaleUp+mety_HEM2018scaleUp*mety_HEM2018scaleUp);
        float met_HEM2018scaleDown = TMath::Sqrt(metx_HEM2018scaleDown*metx_HEM2018scaleDown+mety_HEM2018scaleDown*mety_HEM2018scaleDown);
        
        TLorentzVector metLV_muonScaleUp;
        TLorentzVector metLV_muonScaleDown;
        TLorentzVector metLV_electronScaleUp;
        TLorentzVector metLV_electronScaleDown;
        TLorentzVector metLV_tauScaleUp;
        TLorentzVector metLV_tauScaleDown;
        TLorentzVector metLV_electronfakingtauScaleUp;
        TLorentzVector metLV_electronfakingtauScaleDown;
        TLorentzVector metLV_muonfakingtauScaleUp;
        TLorentzVector metLV_muonfakingtauScaleDown;
        TLorentzVector metLV_unclmetScaleUp;
        TLorentzVector metLV_unclmetScaleDown;
        TLorentzVector metLV_jecScaleUp;
        TLorentzVector metLV_jecScaleDown;
        TLorentzVector metLV_HEM2018ScaleUp;
        TLorentzVector metLV_HEM2018ScaleDown;
        
        metLV_muonScaleUp.SetXYZT(metx_mscaleUp,mety_mscaleUp,0.,met_mscaleUp);
        metLV_muonScaleDown.SetXYZT(metx_mscaleDown,mety_mscaleDown,0.,met_mscaleDown);
        metLV_electronScaleUp.SetXYZT(metx_escaleUp,mety_escaleUp,0.,met_escaleUp);
        metLV_electronScaleDown.SetXYZT(metx_escaleDown,mety_escaleDown,0.,met_escaleDown);
        metLV_tauScaleUp.SetXYZT(metx_tscaleUp,mety_tscaleUp,0.,met_tscaleUp);
        metLV_tauScaleDown.SetXYZT(metx_tscaleDown,mety_tscaleDown,0.,met_tscaleDown);
        metLV_electronfakingtauScaleUp.SetXYZT(metx_eftscaleUp,mety_eftscaleUp,0.,met_eftscaleUp);
        metLV_electronfakingtauScaleDown.SetXYZT(metx_eftscaleDown,mety_eftscaleDown,0.,met_eftscaleDown);
        metLV_muonfakingtauScaleUp.SetXYZT(metx_mftscaleUp,mety_mftscaleUp,0.,met_mftscaleUp);
        metLV_muonfakingtauScaleDown.SetXYZT(metx_mftscaleDown,mety_mftscaleDown,0.,met_mftscaleDown);
        metLV_unclmetScaleUp.SetXYZT(metx_unclmetscaleUp,mety_unclmetscaleUp,0.,met_unclmetscaleUp);
        metLV_unclmetScaleDown.SetXYZT(metx_unclmetscaleDown,mety_unclmetscaleDown,0.,met_unclmetscaleDown);
        metLV_jecScaleUp.SetXYZT(metx_jecscaleUp,mety_jecscaleUp,0.,met_jecscaleUp);
        metLV_jecScaleDown.SetXYZT(metx_jecscaleDown,mety_jecscaleDown,0.,met_jecscaleDown);
        metLV_HEM2018ScaleUp.SetXYZT(metx_HEM2018scaleUp,mety_HEM2018scaleUp,0.,met_HEM2018scaleUp);
        metLV_HEM2018ScaleDown.SetXYZT(metx_HEM2018scaleDown,mety_HEM2018scaleDown,0.,met_HEM2018scaleDown);

        //Corrected uncluster MET Up and Down template with TES shift
        if(!isData)
        {
            if((HiggsFinalState=="ET")||(HiggsFinalState=="MT"))
            {
                correctTauES(Leg4_unCorrected_1, metLV_unclmetScaleUp, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
                correctTauES(Leg4_unCorrected_2, metLV_unclmetScaleDown, tau_Shift(era, &analysisTree, indexLeg4), dm_4);

                correctTauES(Leg4_unCorrected_3, metLV_jecScaleUp, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
                correctTauES(Leg4_unCorrected_4, metLV_jecScaleDown, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
                if(era == "UL2018")
                {
                    correctTauES(Leg4_unCorrected_5, metLV_HEM2018ScaleUp, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
                    correctTauES(Leg4_unCorrected_6, metLV_HEM2018ScaleDown, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
                }
            }
            if(HiggsFinalState=="TT")
            {
                correctTauES(Leg3_unCorrected_1, metLV_unclmetScaleUp, tau_Shift(era, &analysisTree, indexLeg3), dm_3);
                correctTauES(Leg4_unCorrected_1, metLV_unclmetScaleUp, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
                correctTauES(Leg3_unCorrected_2, metLV_unclmetScaleDown, tau_Shift(era, &analysisTree, indexLeg3), dm_3);
                correctTauES(Leg4_unCorrected_2, metLV_unclmetScaleDown, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
                
                correctTauES(Leg3_unCorrected_3, metLV_jecScaleUp, tau_Shift(era, &analysisTree, indexLeg3), dm_3);
                correctTauES(Leg4_unCorrected_3, metLV_jecScaleUp, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
                correctTauES(Leg3_unCorrected_4, metLV_jecScaleDown, tau_Shift(era, &analysisTree, indexLeg3), dm_3);
                correctTauES(Leg4_unCorrected_4, metLV_jecScaleDown, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
                
                if(era == "UL2018")
                {
                    correctTauES(Leg3_unCorrected_5, metLV_HEM2018ScaleUp, tau_Shift(era, &analysisTree, indexLeg3), dm_3);
                    correctTauES(Leg4_unCorrected_5, metLV_HEM2018ScaleUp, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
                    correctTauES(Leg3_unCorrected_6, metLV_HEM2018ScaleDown, tau_Shift(era, &analysisTree, indexLeg3), dm_3);
                    correctTauES(Leg4_unCorrected_6, metLV_HEM2018ScaleDown, tau_Shift(era, &analysisTree, indexLeg4), dm_4);
                }
            }
        }
        //Lepton systematic shifts
        TLorentzVector muon1UpLV;
        TLorentzVector muon1DownLV;
        TLorentzVector muon2UpLV;
        TLorentzVector muon2DownLV;
        
        TLorentzVector Leg3UpLV;
        TLorentzVector Leg3DownLV;
        TLorentzVector Leg3electronfakingtauUpLV;
        TLorentzVector Leg3electronfakingtauDownLV;
        TLorentzVector Leg3muonfakingtauUpLV;
        TLorentzVector Leg3muonfakingtauDownLV;
        
        TLorentzVector Leg4UpLV;
        TLorentzVector Leg4DownLV;
        TLorentzVector Leg4electronfakingtauUpLV;
        TLorentzVector Leg4electronfakingtauDownLV;
        TLorentzVector Leg4muonfakingtauUpLV;
        TLorentzVector Leg4muonfakingtauDownLV;
        
        if(applySystShift && !isData)
        {
            muon1UpLV.SetXYZM((1.0+MuonScale)*analysisTree.muon_px[indexMuon1],
                              (1.0+MuonScale)*analysisTree.muon_py[indexMuon1],
                              (1.0+MuonScale)*analysisTree.muon_pz[indexMuon1],
                              classic_svFit::muonMass);
            muon1DownLV.SetXYZM((1.0-MuonScale)*analysisTree.muon_px[indexMuon1],
                                  (1.0-MuonScale)*analysisTree.muon_py[indexMuon1],
                                  (1.0-MuonScale)*analysisTree.muon_pz[indexMuon1],
                                  classic_svFit::muonMass);
            
            muon2UpLV.SetXYZM((1.0+MuonScale)*analysisTree.muon_px[indexMuon2],
                              (1.0+MuonScale)*analysisTree.muon_py[indexMuon2],
                              (1.0+MuonScale)*analysisTree.muon_pz[indexMuon2],
                              classic_svFit::muonMass);
            muon2DownLV.SetXYZM((1.0-MuonScale)*analysisTree.muon_px[indexMuon2],
                                  (1.0-MuonScale)*analysisTree.muon_py[indexMuon2],
                                  (1.0-MuonScale)*analysisTree.muon_pz[indexMuon2],
                                  classic_svFit::muonMass);
            if(HiggsFinalState == "EM")
            {
                Leg3UpLV.SetXYZM(analysisTree.electron_px_energyscale_up[indexLeg3],
                                 analysisTree.electron_py_energyscale_up[indexLeg3],
                                 analysisTree.electron_pz_energyscale_up[indexLeg3],
                                 classic_svFit::electronMass);
                Leg3DownLV.SetXYZM(analysisTree.electron_px_energyscale_down[indexLeg3],
                                      analysisTree.electron_py_energyscale_down[indexLeg3],
                                      analysisTree.electron_pz_energyscale_down[indexLeg3],
                                      classic_svFit::electronMass);
                Leg4UpLV.SetXYZM((1.0+MuonScale)*analysisTree.muon_px[indexLeg4],
                                      (1.0+MuonScale)*analysisTree.muon_py[indexLeg4],
                                      (1.0+MuonScale)*analysisTree.muon_pz[indexLeg4],
                                      classic_svFit::muonMass);
                Leg4DownLV.SetXYZM((1.0-MuonScale)*analysisTree.muon_px[indexLeg4],
                                      (1.0-MuonScale)*analysisTree.muon_py[indexLeg4],
                                      (1.0-MuonScale)*analysisTree.muon_pz[indexLeg4],
                                      classic_svFit::muonMass);
            }
            
            if(HiggsFinalState == "ET")
            {
                Leg3UpLV.SetXYZM(analysisTree.electron_px_energyscale_up[indexLeg3],
                                 analysisTree.electron_py_energyscale_up[indexLeg3],
                                 analysisTree.electron_pz_energyscale_up[indexLeg3],
                                 classic_svFit::electronMass);
                Leg3DownLV.SetXYZM(analysisTree.electron_px_energyscale_down[indexLeg3],
                                      analysisTree.electron_py_energyscale_down[indexLeg3],
                                      analysisTree.electron_pz_energyscale_down[indexLeg3],
                                      classic_svFit::electronMass);
                float efes_error_up = 0.0;
                float efes_error_down = 0.0;
                tau_EFES_Syst(era,gen_match_4,dm_4,eta_4,efes_error_up,efes_error_down);
                if(dm_4==0)
                {
                    Leg4UpLV.SetXYZM((1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.X(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Y(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Z(),
                                        m_4);
                    Leg4DownLV.SetXYZM((1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.X(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Y(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Z(),
                                          m_4);
                    
                    Leg4electronfakingtauUpLV.SetXYZM((1.0+efes_error_up)*Leg4.X(),
                                                      (1.0+efes_error_up)*Leg4.Y(),
                                                      (1.0+efes_error_up)*Leg4.Z(),
                                                         m_4);
                    Leg4electronfakingtauDownLV.SetXYZM((1.0-efes_error_down)*Leg4.X(),
                                                        (1.0-efes_error_down)*Leg4.Y(),
                                                        (1.0-efes_error_down)*Leg4.Z(),
                                                           m_4);
                    Leg4muonfakingtauUpLV.SetXYZM((1.0+tau_MFES_Syst(gen_match_4))*Leg4.X(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*Leg4.Y(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*Leg4.Z(),
                                                             m_4);
                    Leg4muonfakingtauDownLV.SetXYZM((1.0-tau_MFES_Syst(gen_match_4))*Leg4.X(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*Leg4.Y(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*Leg4.Z(),
                                                               m_4);
                }
                else
                {
                    Leg4UpLV.SetXYZM((1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.X(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Y(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Z(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*m_4);
                    Leg4DownLV.SetXYZM((1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.X(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Y(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Z(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*m_4);
                    
                    Leg4electronfakingtauUpLV.SetXYZM((1.0+efes_error_up)*Leg4.X(),
                                                      (1.0+efes_error_up)*Leg4.Y(),
                                                      (1.0+efes_error_up)*Leg4.Z(),
                                                      (1.0+efes_error_up)*m_4);
                    Leg4electronfakingtauDownLV.SetXYZM((1.0-efes_error_down)*Leg4.X(),
                                                        (1.0-efes_error_down)*Leg4.Y(),
                                                        (1.0-efes_error_down)*Leg4.Z(),
                                                        (1.0-efes_error_down)*m_4);
                    Leg4muonfakingtauUpLV.SetXYZM((1.0+tau_MFES_Syst(gen_match_4))*Leg4.X(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*Leg4.Y(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*Leg4.Z(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*m_4);
                    Leg4muonfakingtauDownLV.SetXYZM((1.0-tau_MFES_Syst(gen_match_4))*Leg4.X(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*Leg4.Y(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*Leg4.Z(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*m_4);
                }
            }
            
            if(HiggsFinalState == "MT")
            {
                Leg3UpLV.SetXYZM((1.0+MuonScale)*analysisTree.muon_px[indexLeg3],
                                      (1.0+MuonScale)*analysisTree.muon_py[indexLeg3],
                                      (1.0+MuonScale)*analysisTree.muon_pz[indexLeg3],
                                      classic_svFit::muonMass);
                Leg3DownLV.SetXYZM((1.0-MuonScale)*analysisTree.muon_px[indexLeg3],
                                      (1.0-MuonScale)*analysisTree.muon_py[indexLeg3],
                                      (1.0-MuonScale)*analysisTree.muon_pz[indexLeg3],
                                      classic_svFit::muonMass);
                float efes_error_up = 0.0;
                float efes_error_down = 0.0;
                tau_EFES_Syst(era,gen_match_4,dm_4,eta_4,efes_error_up,efes_error_down);
                if(dm_4==0)
                {
                    Leg4UpLV.SetXYZM((1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.X(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Y(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Z(),
                                        m_4);
                    Leg4DownLV.SetXYZM((1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.X(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Y(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Z(),
                                          m_4);
                    
                    Leg4electronfakingtauUpLV.SetXYZM((1.0+efes_error_up)*Leg4.X(),
                                                      (1.0+efes_error_up)*Leg4.Y(),
                                                      (1.0+efes_error_up)*Leg4.Z(),
                                                         m_4);
                    Leg4electronfakingtauDownLV.SetXYZM((1.0-efes_error_down)*Leg4.X(),
                                                        (1.0-efes_error_down)*Leg4.Y(),
                                                        (1.0-efes_error_down)*Leg4.Z(),
                                                           m_4);
                    Leg4muonfakingtauUpLV.SetXYZM((1.0+tau_MFES_Syst(gen_match_4))*Leg4.X(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*Leg4.Y(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*Leg4.Z(),
                                                             m_4);
                    Leg4muonfakingtauDownLV.SetXYZM((1.0-tau_MFES_Syst(gen_match_4))*Leg4.X(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*Leg4.Y(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*Leg4.Z(),
                                                               m_4);
                }
                else
                {
                    Leg4UpLV.SetXYZM((1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.X(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Y(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Z(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*m_4);
                    Leg4DownLV.SetXYZM((1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.X(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Y(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Z(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*m_4);
                    
                    Leg4electronfakingtauUpLV.SetXYZM((1.0+efes_error_up)*Leg4.X(),
                                                      (1.0+efes_error_up)*Leg4.Y(),
                                                      (1.0+efes_error_up)*Leg4.Z(),
                                                      (1.0+efes_error_up)*m_4);
                    Leg4electronfakingtauDownLV.SetXYZM((1.0-efes_error_down)*Leg4.X(),
                                                        (1.0-efes_error_down)*Leg4.Y(),
                                                        (1.0-efes_error_down)*Leg4.Z(),
                                                        (1.0-efes_error_down)*m_4);
                    Leg4muonfakingtauUpLV.SetXYZM((1.0+tau_MFES_Syst(gen_match_4))*Leg4.X(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*Leg4.Y(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*Leg4.Z(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*m_4);
                    Leg4muonfakingtauDownLV.SetXYZM((1.0-tau_MFES_Syst(gen_match_4))*Leg4.X(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*Leg4.Y(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*Leg4.Z(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*m_4);
                }
            }
            if(HiggsFinalState == "TT")
            {
                float efes_error_up_3 = 0.0;
                float efes_error_down_3 = 0.0;
                tau_EFES_Syst(era,gen_match_3,dm_3,eta_3,efes_error_up_3,efes_error_down_3);
                if(dm_3==0)
                {
                    Leg3UpLV.SetXYZM((1.0+tau_ES_Syst(era,gen_match_3,dm_3))*Leg3.X(),
                                     (1.0+tau_ES_Syst(era,gen_match_3,dm_3))*Leg3.Y(),
                                     (1.0+tau_ES_Syst(era,gen_match_3,dm_3))*Leg3.Z(),
                                        m_3);
                    Leg3DownLV.SetXYZM((1.0-tau_ES_Syst(era,gen_match_3,dm_3))*Leg3.X(),
                                       (1.0-tau_ES_Syst(era,gen_match_3,dm_3))*Leg3.Y(),
                                       (1.0-tau_ES_Syst(era,gen_match_3,dm_3))*Leg3.Z(),
                                          m_3);
                    Leg3electronfakingtauUpLV.SetXYZM((1.0+efes_error_up_3)*Leg3.X(),
                                                      (1.0+efes_error_up_3)*Leg3.Y(),
                                                      (1.0+efes_error_up_3)*Leg3.Z(),
                                                         m_3);
                    Leg3electronfakingtauDownLV.SetXYZM((1.0-efes_error_down_3)*Leg3.X(),
                                                        (1.0-efes_error_down_3)*Leg3.Y(),
                                                        (1.0-efes_error_down_3)*Leg3.Z(),
                                                           m_3);
                    Leg3muonfakingtauUpLV.SetXYZM((1.0+tau_MFES_Syst(gen_match_3))*Leg3.X(),
                                                      (1.0+tau_MFES_Syst(gen_match_3))*Leg3.Y(),
                                                      (1.0+tau_MFES_Syst(gen_match_3))*Leg3.Z(),
                                                             m_3);
                    Leg3muonfakingtauDownLV.SetXYZM((1.0-tau_MFES_Syst(gen_match_3))*Leg3.X(),
                                                            (1.0-tau_MFES_Syst(gen_match_3))*Leg3.Y(),
                                                            (1.0-tau_MFES_Syst(gen_match_3))*Leg3.Z(),
                                                               m_3);
                }
                else
                {
                    Leg3UpLV.SetXYZM((1.0+tau_ES_Syst(era,gen_match_3,dm_3))*Leg3.X(),
                                     (1.0+tau_ES_Syst(era,gen_match_3,dm_3))*Leg3.Y(),
                                     (1.0+tau_ES_Syst(era,gen_match_3,dm_3))*Leg3.Z(),
                                     (1.0+tau_ES_Syst(era,gen_match_3,dm_3))*m_3);
                    Leg3DownLV.SetXYZM((1.0-tau_ES_Syst(era,gen_match_3,dm_3))*Leg3.X(),
                                       (1.0-tau_ES_Syst(era,gen_match_3,dm_3))*Leg3.Y(),
                                       (1.0-tau_ES_Syst(era,gen_match_3,dm_3))*Leg3.Z(),
                                       (1.0-tau_ES_Syst(era,gen_match_3,dm_3))*m_3);
                    Leg3electronfakingtauUpLV.SetXYZM((1.0+efes_error_up_3)*Leg3.X(),
                                                      (1.0+efes_error_up_3)*Leg3.Y(),
                                                      (1.0+efes_error_up_3)*Leg3.Z(),
                                                      (1.0+efes_error_up_3)*m_3);
                    Leg3electronfakingtauDownLV.SetXYZM((1.0-efes_error_down_3)*Leg3.X(),
                                                        (1.0-efes_error_down_3)*Leg3.Y(),
                                                        (1.0-efes_error_down_3)*Leg3.Z(),
                                                        (1.0-efes_error_down_3)*m_3);
                    Leg3muonfakingtauUpLV.SetXYZM((1.0+tau_MFES_Syst(gen_match_3))*Leg3.X(),
                                                      (1.0+tau_MFES_Syst(gen_match_3))*Leg3.Y(),
                                                      (1.0+tau_MFES_Syst(gen_match_3))*Leg3.Z(),
                                                      (1.0+tau_MFES_Syst(gen_match_3))*m_3);
                    Leg3muonfakingtauDownLV.SetXYZM((1.0-tau_MFES_Syst(gen_match_3))*Leg3.X(),
                                                            (1.0-tau_MFES_Syst(gen_match_3))*Leg3.Y(),
                                                            (1.0-tau_MFES_Syst(gen_match_3))*Leg3.Z(),
                                                    (1.0-tau_MFES_Syst(gen_match_3))*m_3);
                }
                float efes_error_up_4 = 0.0;
                float efes_error_down_4 = 0.0;
                tau_EFES_Syst(era,gen_match_4,dm_4,eta_4,efes_error_up_4,efes_error_down_4);
                if(dm_4==0)
                {
                    Leg4UpLV.SetXYZM((1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.X(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Y(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Z(),
                                        m_4);
                    Leg4DownLV.SetXYZM((1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.X(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Y(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Z(),
                                          m_4);
                    Leg4electronfakingtauUpLV.SetXYZM((1.0+efes_error_up_4)*Leg4.X(),
                                                      (1.0+efes_error_up_4)*Leg4.Y(),
                                                      (1.0+efes_error_up_4)*Leg4.Z(),
                                                         m_4);
                    Leg4electronfakingtauDownLV.SetXYZM((1.0-efes_error_down_4)*Leg4.X(),
                                                        (1.0-efes_error_down_4)*Leg4.Y(),
                                                        (1.0-efes_error_down_4)*Leg4.Z(),
                                                           m_4);
                    Leg4muonfakingtauUpLV.SetXYZM((1.0+tau_MFES_Syst(gen_match_4))*Leg4.X(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*Leg4.Y(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*Leg4.Z(),
                                                             m_4);
                    Leg4muonfakingtauDownLV.SetXYZM((1.0-tau_MFES_Syst(gen_match_4))*Leg4.X(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*Leg4.Y(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*Leg4.Z(),
                                                               m_4);
                }
                else
                {
                    Leg4UpLV.SetXYZM((1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.X(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Y(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Z(),
                                     (1.0+tau_ES_Syst(era,gen_match_4,dm_4))*m_4);
                    Leg4DownLV.SetXYZM((1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.X(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Y(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*Leg4.Z(),
                                       (1.0-tau_ES_Syst(era,gen_match_4,dm_4))*m_4);
                    Leg4electronfakingtauUpLV.SetXYZM((1.0+efes_error_up_4)*Leg4.X(),
                                                      (1.0+efes_error_up_4)*Leg4.Y(),
                                                      (1.0+efes_error_up_4)*Leg4.Z(),
                                                      (1.0+efes_error_up_4)*m_4);
                    Leg4electronfakingtauDownLV.SetXYZM((1.0-efes_error_down_4)*Leg4.X(),
                                                        (1.0-efes_error_down_4)*Leg4.Y(),
                                                        (1.0-efes_error_down_4)*Leg4.Z(),
                                                        (1.0-efes_error_down_4)*m_4);
                    Leg4muonfakingtauUpLV.SetXYZM((1.0+tau_MFES_Syst(gen_match_4))*Leg4.X(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*Leg4.Y(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*Leg4.Z(),
                                                      (1.0+tau_MFES_Syst(gen_match_4))*m_4);
                    Leg4muonfakingtauDownLV.SetXYZM((1.0-tau_MFES_Syst(gen_match_4))*Leg4.X(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*Leg4.Y(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*Leg4.Z(),
                                                            (1.0-tau_MFES_Syst(gen_match_4))*m_4);
                }
            }
        }
        
        //Fake rate measurement selections
        float m2l_1, m2l_2, m2l_3, m3l;
        
        m2l_1 = (muon1+muon2).M();
        m2l_2 = (muon1+Leg3).M();
        m2l_3 = (muon2+Leg3).M();
        
        m3l = (muon1+muon2+Leg3).M();
        
        if((HiggsFinalState == "E"))
        {
            if(mll < 81.2)
                continue;
            
            if((m2l_1<5)||(m2l_2<5)||(m2l_3<5))
                continue;
            if(TMath::Max(iso_1,iso_2)>0.14)
                continue;
            if(abs(Leg3.Eta())<1.5)
            {
                if(pfmt_3>55)
                    continue;
            }
        }
        
        if((HiggsFinalState == "M"))//Fake rate measurement selections
        {
            if(mll < 81.2)
                continue;
            if(m3l > 250)
                continue;
            if(abs(Leg3.Eta())<1.2)
            {
                if(pfmt_3>55)
                    continue;
            }
        }
        
        if (isData)
            histWeightsmllH->Fill(0.,1.);
        else
            histWeightsmllH->Fill(0.,analysisTree.genweight);
        
        Mass = (muon1+muon2+Leg3+Leg4).M();
        H_DR = deltaR(Leg3.Eta(),Leg3.Phi(),Leg4.Eta(),Leg4.Phi());
        DR_12 = Z_DR;
        DR_13 = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],Leg3.Eta(),Leg3.Phi());
        DR_14 = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],Leg4.Eta(),Leg4.Phi());
        DR_23 = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],Leg3.Eta(),Leg3.Phi());
        DR_24 = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],Leg4.Eta(),Leg4.Phi());
        DR_34 = H_DR;
        
        Mass_gen = 0;
        //Gen mass
        for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {
            TLorentzVector genALV; genALV.SetXYZT(analysisTree.genparticles_px[igen],
                                                  analysisTree.genparticles_py[igen],
                                                  analysisTree.genparticles_pz[igen],
                                                  analysisTree.genparticles_e[igen]);
            
            if (fabs(analysisTree.genparticles_pdgid[igen])==36) {
                Mass_gen = genALV.M();
            }
        }

        //SV fit und FastMTT
        if(ApplySVFit||ApplyFastMTT)
        {
            double measuredMETx = met*TMath::Cos(metphi);
            double measuredMETy = met*TMath::Sin(metphi);
            
            TMatrixD covMET(2, 2);
            // using PF MET
            covMET[0][0] = metcov00;
            covMET[1][0] = metcov10;
            covMET[0][1] = metcov01;
            covMET[1][1] = metcov11;
            
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_MuonES_Up;
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_MuonES_Down;
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_ElectronES_Up;
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_ElectronES_Down;
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_TauES_Up;
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_TauES_Down;
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_ElectronFakingTauES_Up;
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_ElectronFakingTauES_Down;
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_MuonFakingTauES_Up;
            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_MuonFakingTauES_Down;

            if(HiggsFinalState == "MT")
            {
                classic_svFit::MeasuredTauLepton svFitLep1(classic_svFit::MeasuredTauLepton::kTauToMuDecay, pt_3, eta_3, phi_3, m_3);
                classic_svFit::MeasuredTauLepton svFitLep2(classic_svFit::MeasuredTauLepton::kTauToHadDecay, pt_4, eta_4, phi_4, m_4, dm_4);
                measuredTauLeptons.push_back(svFitLep1);
                measuredTauLeptons.push_back(svFitLep2);
                if(applySystShift && !isData)
                {
                    classic_svFit::MeasuredTauLepton svFitMuUP(classic_svFit::MeasuredTauLepton::kTauToMuDecay, Leg3UpLV.Pt(), Leg3UpLV.Eta(), Leg3UpLV.Phi(), Leg3UpLV.M());
                    measuredTauLeptons_MuonES_Up.push_back(svFitMuUP);
                    measuredTauLeptons_MuonES_Up.push_back(svFitLep2);
                    
                    classic_svFit::MeasuredTauLepton svFitMuDown(classic_svFit::MeasuredTauLepton::kTauToMuDecay, Leg3DownLV.Pt(), Leg3DownLV.Eta(), Leg3DownLV.Phi(), Leg3DownLV.M());
                    measuredTauLeptons_MuonES_Down.push_back(svFitMuDown);
                    measuredTauLeptons_MuonES_Down.push_back(svFitLep2);
                    
                    classic_svFit::MeasuredTauLepton svFitTauUP(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4UpLV.Pt(), Leg4UpLV.Eta(), Leg4UpLV.Phi(), Leg4UpLV.M());
                    measuredTauLeptons_TauES_Up.push_back(svFitLep1);
                    measuredTauLeptons_TauES_Up.push_back(svFitTauUP);
              
                    classic_svFit::MeasuredTauLepton svFitTauDown(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4DownLV.Pt(), Leg4DownLV.Eta(), Leg4DownLV.Phi(), Leg4DownLV.M());
                    measuredTauLeptons_TauES_Down.push_back(svFitLep1);
                    measuredTauLeptons_TauES_Down.push_back(svFitTauDown);
                    
                    classic_svFit::MeasuredTauLepton svFitElectronFakingTauUp(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4electronfakingtauUpLV.Pt(), Leg4electronfakingtauUpLV.Eta(), Leg4electronfakingtauUpLV.Phi(), Leg4electronfakingtauUpLV.M());
                    measuredTauLeptons_ElectronFakingTauES_Up.push_back(svFitLep1);
                    measuredTauLeptons_ElectronFakingTauES_Up.push_back(svFitElectronFakingTauUp);
              
                    classic_svFit::MeasuredTauLepton svFitElectronFakingTauDown(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4electronfakingtauDownLV.Pt(), Leg4electronfakingtauDownLV.Eta(), Leg4electronfakingtauDownLV.Phi(), Leg4electronfakingtauDownLV.M());
                    measuredTauLeptons_ElectronFakingTauES_Down.push_back(svFitLep1);
                    measuredTauLeptons_ElectronFakingTauES_Down.push_back(svFitElectronFakingTauDown);
                    
                    classic_svFit::MeasuredTauLepton svFitMuonFakingTauUp(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4muonfakingtauUpLV.Pt(), Leg4muonfakingtauUpLV.Eta(), Leg4muonfakingtauUpLV.Phi(), Leg4muonfakingtauUpLV.M());
                    measuredTauLeptons_MuonFakingTauES_Up.push_back(svFitLep1);
                    measuredTauLeptons_MuonFakingTauES_Up.push_back(svFitMuonFakingTauUp);
              
                    classic_svFit::MeasuredTauLepton svFitMuonFakingTauDown(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4muonfakingtauDownLV.Pt(), Leg4muonfakingtauDownLV.Eta(), Leg4muonfakingtauDownLV.Phi(), Leg4muonfakingtauDownLV.M());
                    measuredTauLeptons_MuonFakingTauES_Down.push_back(svFitLep1);
                    measuredTauLeptons_MuonFakingTauES_Down.push_back(svFitMuonFakingTauDown);
                }
            }
            
            if(HiggsFinalState == "ET")
            {
                classic_svFit::MeasuredTauLepton svFitLep1(classic_svFit::MeasuredTauLepton::kTauToElecDecay, pt_3, eta_3, phi_3, m_3);
                classic_svFit::MeasuredTauLepton svFitLep2(classic_svFit::MeasuredTauLepton::kTauToHadDecay, pt_4, eta_4, phi_4, m_4, dm_4);
                measuredTauLeptons.push_back(svFitLep1);
                measuredTauLeptons.push_back(svFitLep2);
                if(applySystShift && !isData)
                {
                    classic_svFit::MeasuredTauLepton svFitEleUP(classic_svFit::MeasuredTauLepton::kTauToElecDecay, Leg3UpLV.Pt(), Leg3UpLV.Eta(), Leg3UpLV.Phi(), Leg3UpLV.M());
                    measuredTauLeptons_ElectronES_Up.push_back(svFitEleUP);
                    measuredTauLeptons_ElectronES_Up.push_back(svFitLep2);
              
                    classic_svFit::MeasuredTauLepton svFitEleDown(classic_svFit::MeasuredTauLepton::kTauToElecDecay, Leg3DownLV.Pt(), Leg3DownLV.Eta(), Leg3DownLV.Phi(), Leg3DownLV.M());
                    measuredTauLeptons_ElectronES_Down.push_back(svFitEleDown);
                    measuredTauLeptons_ElectronES_Down.push_back(svFitLep2);
                    
                    classic_svFit::MeasuredTauLepton svFitTauUP(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4UpLV.Pt(), Leg4UpLV.Eta(), Leg4UpLV.Phi(), Leg4UpLV.M());
                    measuredTauLeptons_TauES_Up.push_back(svFitLep1);
                    measuredTauLeptons_TauES_Up.push_back(svFitTauUP);
              
                    classic_svFit::MeasuredTauLepton svFitTauDown(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4DownLV.Pt(), Leg4DownLV.Eta(), Leg4DownLV.Phi(), Leg4DownLV.M());
                    measuredTauLeptons_TauES_Down.push_back(svFitLep1);
                    measuredTauLeptons_TauES_Down.push_back(svFitTauDown);
                    
                    classic_svFit::MeasuredTauLepton svFitElectronFakingTauUp(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4electronfakingtauUpLV.Pt(), Leg4electronfakingtauUpLV.Eta(), Leg4electronfakingtauUpLV.Phi(), Leg4electronfakingtauUpLV.M());
                    measuredTauLeptons_ElectronFakingTauES_Up.push_back(svFitLep1);
                    measuredTauLeptons_ElectronFakingTauES_Up.push_back(svFitElectronFakingTauUp);
              
                    classic_svFit::MeasuredTauLepton svFitElectronFakingTauDown(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4electronfakingtauDownLV.Pt(), Leg4electronfakingtauDownLV.Eta(), Leg4electronfakingtauDownLV.Phi(), Leg4electronfakingtauDownLV.M());
                    measuredTauLeptons_ElectronFakingTauES_Down.push_back(svFitLep1);
                    measuredTauLeptons_ElectronFakingTauES_Down.push_back(svFitElectronFakingTauDown);
                    
                    classic_svFit::MeasuredTauLepton svFitMuonFakingTauUp(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4muonfakingtauUpLV.Pt(), Leg4muonfakingtauUpLV.Eta(), Leg4muonfakingtauUpLV.Phi(), Leg4muonfakingtauUpLV.M());
                    measuredTauLeptons_MuonFakingTauES_Up.push_back(svFitLep1);
                    measuredTauLeptons_MuonFakingTauES_Up.push_back(svFitMuonFakingTauUp);
              
                    classic_svFit::MeasuredTauLepton svFitMuonFakingTauDown(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4muonfakingtauDownLV.Pt(), Leg4muonfakingtauDownLV.Eta(), Leg4muonfakingtauDownLV.Phi(), Leg4muonfakingtauDownLV.M());
                    measuredTauLeptons_MuonFakingTauES_Down.push_back(svFitLep1);
                    measuredTauLeptons_MuonFakingTauES_Down.push_back(svFitMuonFakingTauDown);
                }
            }
            if(HiggsFinalState == "EM")
            {
                classic_svFit::MeasuredTauLepton svFitLep1(classic_svFit::MeasuredTauLepton::kTauToElecDecay, pt_3, eta_3, phi_3, m_3);
                classic_svFit::MeasuredTauLepton svFitLep2(classic_svFit::MeasuredTauLepton::kTauToMuDecay, pt_4, eta_4, phi_4, m_4);
                measuredTauLeptons.push_back(svFitLep1);
                measuredTauLeptons.push_back(svFitLep2);
                if(applySystShift && !isData)
                {
                    classic_svFit::MeasuredTauLepton svFitMuUP(classic_svFit::MeasuredTauLepton::kTauToMuDecay, Leg4UpLV.Pt(), Leg4UpLV.Eta(), Leg4UpLV.Phi(), Leg4UpLV.M());
                    measuredTauLeptons_MuonES_Up.push_back(svFitLep1);
                    measuredTauLeptons_MuonES_Up.push_back(svFitMuUP);
                    
                    classic_svFit::MeasuredTauLepton svFitMuDown(classic_svFit::MeasuredTauLepton::kTauToMuDecay, Leg4DownLV.Pt(), Leg4DownLV.Eta(), Leg4DownLV.Phi(), Leg4DownLV.M());
                    measuredTauLeptons_MuonES_Down.push_back(svFitLep1);
                    measuredTauLeptons_MuonES_Down.push_back(svFitMuDown);
                    
                    classic_svFit::MeasuredTauLepton svFitEleUP(classic_svFit::MeasuredTauLepton::kTauToElecDecay, Leg3UpLV.Pt(), Leg3UpLV.Eta(), Leg3UpLV.Phi(), Leg3UpLV.M());
                    measuredTauLeptons_ElectronES_Up.push_back(svFitEleUP);
                    measuredTauLeptons_ElectronES_Up.push_back(svFitLep2);
              
                    classic_svFit::MeasuredTauLepton svFitEleDown(classic_svFit::MeasuredTauLepton::kTauToElecDecay, Leg3DownLV.Pt(), Leg3DownLV.Eta(), Leg3DownLV.Phi(), Leg3DownLV.M());
                    measuredTauLeptons_ElectronES_Down.push_back(svFitEleDown);
                    measuredTauLeptons_ElectronES_Down.push_back(svFitLep2);
                }
            }
            if(HiggsFinalState == "TT")
            {
                classic_svFit::MeasuredTauLepton svFitLep1(classic_svFit::MeasuredTauLepton::kTauToHadDecay, pt_3, eta_3, phi_3, m_3, dm_3);
                classic_svFit::MeasuredTauLepton svFitLep2(classic_svFit::MeasuredTauLepton::kTauToHadDecay, pt_4, eta_4, phi_4, m_4, dm_4);
                measuredTauLeptons.push_back(svFitLep1);
                measuredTauLeptons.push_back(svFitLep2);
                if(applySystShift && !isData)
                {
                    classic_svFit::MeasuredTauLepton svFitTau1Up(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg3UpLV.Pt(), Leg3UpLV.Eta(), Leg3UpLV.Phi(), Leg3UpLV.M());
                    classic_svFit::MeasuredTauLepton svFitTau2Up(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4UpLV.Pt(), Leg4UpLV.Eta(), Leg4UpLV.Phi(), Leg4UpLV.M());
                    measuredTauLeptons_TauES_Up.push_back(svFitTau1Up);
                    measuredTauLeptons_TauES_Up.push_back(svFitTau2Up);
              
                    classic_svFit::MeasuredTauLepton svFitTau1Down(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg3DownLV.Pt(), Leg3DownLV.Eta(), Leg3DownLV.Phi(), Leg3DownLV.M());
                    classic_svFit::MeasuredTauLepton svFitTau2Down(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4DownLV.Pt(), Leg4DownLV.Eta(), Leg4DownLV.Phi(), Leg4DownLV.M());
                    measuredTauLeptons_TauES_Down.push_back(svFitTau1Down);
                    measuredTauLeptons_TauES_Down.push_back(svFitTau2Down);
                    
                    classic_svFit::MeasuredTauLepton svFitElectronFakingTau1Up(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg3electronfakingtauUpLV.Pt(), Leg3electronfakingtauUpLV.Eta(), Leg3electronfakingtauUpLV.Phi(), Leg3electronfakingtauUpLV.M());
                    classic_svFit::MeasuredTauLepton svFitElectronFakingTau2Up(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4electronfakingtauUpLV.Pt(), Leg4electronfakingtauUpLV.Eta(), Leg4electronfakingtauUpLV.Phi(), Leg4electronfakingtauUpLV.M());
                    measuredTauLeptons_ElectronFakingTauES_Up.push_back(svFitElectronFakingTau1Up);
                    measuredTauLeptons_ElectronFakingTauES_Up.push_back(svFitElectronFakingTau2Up);
              
                    classic_svFit::MeasuredTauLepton svFitElectronFakingTau1Down(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg3electronfakingtauDownLV.Pt(), Leg3electronfakingtauDownLV.Eta(), Leg3electronfakingtauDownLV.Phi(), Leg3electronfakingtauDownLV.M());
                    classic_svFit::MeasuredTauLepton svFitElectronFakingTau2Down(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4electronfakingtauDownLV.Pt(), Leg4electronfakingtauDownLV.Eta(), Leg4electronfakingtauDownLV.Phi(), Leg4electronfakingtauDownLV.M());
                    measuredTauLeptons_ElectronFakingTauES_Down.push_back(svFitElectronFakingTau1Down);
                    measuredTauLeptons_ElectronFakingTauES_Down.push_back(svFitElectronFakingTau2Down);
                    
                    classic_svFit::MeasuredTauLepton svFitMuonFakingTau1Up(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg3muonfakingtauUpLV.Pt(), Leg3muonfakingtauUpLV.Eta(), Leg3muonfakingtauUpLV.Phi(), Leg3muonfakingtauUpLV.M());
                    classic_svFit::MeasuredTauLepton svFitMuonFakingTau2Up(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4muonfakingtauUpLV.Pt(), Leg4muonfakingtauUpLV.Eta(), Leg4muonfakingtauUpLV.Phi(), Leg4muonfakingtauUpLV.M());
                    measuredTauLeptons_MuonFakingTauES_Up.push_back(svFitMuonFakingTau1Up);
                    measuredTauLeptons_MuonFakingTauES_Up.push_back(svFitMuonFakingTau2Up);
              
                    classic_svFit::MeasuredTauLepton svFitMuonFakingTau1Down(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg3muonfakingtauDownLV.Pt(), Leg3muonfakingtauDownLV.Eta(), Leg3muonfakingtauDownLV.Phi(), Leg3muonfakingtauDownLV.M());
                    classic_svFit::MeasuredTauLepton svFitMuonFakingTau2Down(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Leg4muonfakingtauDownLV.Pt(), Leg4muonfakingtauDownLV.Eta(), Leg4muonfakingtauDownLV.Phi(), Leg4muonfakingtauDownLV.M());
                    measuredTauLeptons_MuonFakingTauES_Down.push_back(svFitMuonFakingTau1Down);
                    measuredTauLeptons_MuonFakingTauES_Down.push_back(svFitMuonFakingTau2Down);
                }
            }
            if (ApplySVFit)
            {

                ClassicSVfit svFitAlgo = SVFitMassComputation(HiggsFinalState, measuredTauLeptons, measuredMETx, measuredMETy, covMET, inputFile_visPtResolution);
                m_sv   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getMass();
                pt_sv  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getPt();
                eta_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getEta();
                phi_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getPhi();
                //met_sv = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getfittedMET().Rho();
                mt_sv  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
                
                TLorentzVector ditau_sv;
                ditau_sv.SetPtEtaPhiM(pt_sv, eta_sv, phi_sv, m_sv);
                Mass_sv = (muon1+muon2+ditau_sv).M();
                res_sv = (Mass_sv - Mass_gen)/Mass_sv;

                ClassicSVfit svFitAlgo_C = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, measuredMETx, measuredMETy, covMET, inputFile_visPtResolution);
                m_svc   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C.getHistogramAdapter())->getMass();
                pt_svc  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C.getHistogramAdapter())->getPt();
                eta_svc = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C.getHistogramAdapter())->getEta();
                phi_svc = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C.getHistogramAdapter())->getPhi();
                //met_svc = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo_C.getHistogramAdapter())->getfittedMET().Rho();
                mt_svc  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C.getHistogramAdapter())->getTransverseMass();
                
                TLorentzVector ditau_svc;
                ditau_svc.SetPtEtaPhiM(pt_svc, eta_svc, phi_svc, m_svc);
                Mass_svc = (muon1+muon2+ditau_svc).M();
                res_svc = (Mass_svc - Mass_gen)/Mass_svc;
                
                if(applySystShift && !isData)
                {
                    if(HiggsFinalState == "EM")
                    {
                        ClassicSVfit svFitAlgo_C_MuonESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_MuonES_Up, metLV_muonScaleUp.X(), metLV_muonScaleUp.Y(), covMET, inputFile_visPtResolution);
                        ClassicSVfit svFitAlgo_C_MuonESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_MuonES_Down, metLV_muonScaleDown.X(), metLV_muonScaleDown.Y(), covMET, inputFile_visPtResolution);
                        
                        m_svc_mscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESUp.getHistogramAdapter())->getMass();
                        pt_svc_mscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESUp.getHistogramAdapter())->getPt();
                        eta_svc_mscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESUp.getHistogramAdapter())->getEta();
                        phi_svc_mscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESUp.getHistogramAdapter())->getPhi();
                        mt_svc_mscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_mscaleUp;
                        ditau_svc_mscaleUp.SetPtEtaPhiM(pt_svc_mscaleUp, eta_svc_mscaleUp, phi_svc_mscaleUp, m_svc_mscaleUp);
                        Mass_svc_mscaleUp = (muon1UpLV+muon2UpLV+ditau_svc_mscaleUp).M();
                        res_svc_mscaleUp = (Mass_svc_mscaleUp - Mass_gen)/Mass_svc_mscaleUp;
                        
                        m_svc_mscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESDown.getHistogramAdapter())->getMass();
                        pt_svc_mscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESDown.getHistogramAdapter())->getPt();
                        eta_svc_mscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESDown.getHistogramAdapter())->getEta();
                        phi_svc_mscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESDown.getHistogramAdapter())->getPhi();
                        mt_svc_mscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_mscaleDown;
                        ditau_svc_mscaleDown.SetPtEtaPhiM(pt_svc_mscaleDown, eta_svc_mscaleDown, phi_svc_mscaleDown, m_svc_mscaleDown);
                        Mass_svc_mscaleDown = (muon1DownLV+muon2DownLV+ditau_svc_mscaleDown).M();
                        res_svc_mscaleDown = (Mass_svc_mscaleDown - Mass_gen)/Mass_svc_mscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_ElectronESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_ElectronES_Up, metLV_electronScaleUp.X(), metLV_electronScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_escaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESUp.getHistogramAdapter())->getMass();
                        pt_svc_escaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESUp.getHistogramAdapter())->getPt();
                        eta_svc_escaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESUp.getHistogramAdapter())->getEta();
                        phi_svc_escaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESUp.getHistogramAdapter())->getPhi();
                        mt_svc_escaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_escaleUp;
                        ditau_svc_escaleUp.SetPtEtaPhiM(pt_svc_escaleUp, eta_svc_escaleUp, phi_svc_escaleUp, m_svc_escaleUp);
                        Mass_svc_escaleUp = (muon1+muon2+ditau_svc_escaleUp).M();
                        res_svc_escaleUp = (Mass_svc_escaleUp - Mass_gen)/Mass_svc_escaleUp;
                        
                        ClassicSVfit svFitAlgo_C_ElectronESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_ElectronES_Down, metLV_electronScaleDown.X(), metLV_electronScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_escaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESDown.getHistogramAdapter())->getMass();
                        pt_svc_escaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESDown.getHistogramAdapter())->getPt();
                        eta_svc_escaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESDown.getHistogramAdapter())->getEta();
                        phi_svc_escaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESDown.getHistogramAdapter())->getPhi();
                        mt_svc_escaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_escaleDown;
                        ditau_svc_escaleDown.SetPtEtaPhiM(pt_svc_escaleDown, eta_svc_escaleDown, phi_svc_escaleDown, m_svc_escaleDown);
                        Mass_svc_escaleDown = (muon1+muon2+ditau_svc_escaleDown).M();
                        res_svc_escaleDown = (Mass_svc_escaleDown - Mass_gen)/Mass_svc_escaleDown;
                        
                        ClassicSVfit svFitAlgo_C_UnclMETESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_unclmetScaleUp.X(), metLV_unclmetScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_unclmetscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getMass();
                        pt_svc_unclmetscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getPt();
                        eta_svc_unclmetscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getEta();
                        phi_svc_unclmetscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getPhi();
                        mt_svc_unclmetscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_unclmetscaleUp;
                        ditau_svc_unclmetscaleUp.SetPtEtaPhiM(pt_svc_unclmetscaleUp, eta_svc_unclmetscaleUp, phi_svc_unclmetscaleUp, m_svc_unclmetscaleUp);
                        Mass_svc_unclmetscaleUp = (muon1+muon2+ditau_svc_unclmetscaleUp).M();
                        res_svc_unclmetscaleUp = (Mass_svc_unclmetscaleUp - Mass_gen)/Mass_svc_unclmetscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_UnclMETESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_unclmetScaleDown.X(), metLV_unclmetScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_unclmetscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getMass();
                        pt_svc_unclmetscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getPt();
                        eta_svc_unclmetscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getEta();
                        phi_svc_unclmetscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getPhi();
                        mt_svc_unclmetscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_unclmetscaleDown;
                        ditau_svc_unclmetscaleDown.SetPtEtaPhiM(pt_svc_unclmetscaleDown, eta_svc_unclmetscaleDown, phi_svc_unclmetscaleDown, m_svc_unclmetscaleDown);
                        Mass_svc_unclmetscaleDown = (muon1+muon2+ditau_svc_unclmetscaleDown).M();
                        res_svc_unclmetscaleDown = (Mass_svc_unclmetscaleDown - Mass_gen)/Mass_svc_unclmetscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_JESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_jecScaleUp.X(), metLV_jecScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_jecscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getMass();
                        pt_svc_jecscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getPt();
                        eta_svc_jecscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getEta();
                        phi_svc_jecscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getPhi();
                        mt_svc_jecscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_jecscaleUp;
                        ditau_svc_jecscaleUp.SetPtEtaPhiM(pt_svc_jecscaleUp, eta_svc_jecscaleUp, phi_svc_jecscaleUp, m_svc_jecscaleUp);
                        Mass_svc_jecscaleUp = (muon1+muon2+ditau_svc_jecscaleUp).M();
                        res_svc_jecscaleUp = (Mass_svc_jecscaleUp - Mass_gen)/Mass_svc_jecscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_JESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_jecScaleDown.X(), metLV_jecScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_jecscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getMass();
                        pt_svc_jecscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getPt();
                        eta_svc_jecscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getEta();
                        phi_svc_jecscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getPhi();
                        mt_svc_jecscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_jecscaleDown;
                        ditau_svc_jecscaleDown.SetPtEtaPhiM(pt_svc_jecscaleDown, eta_svc_jecscaleDown, phi_svc_jecscaleDown, m_svc_jecscaleDown);
                        Mass_svc_jecscaleDown = (muon1+muon2+ditau_svc_jecscaleDown).M();
                        res_svc_jecscaleDown = (Mass_svc_jecscaleDown - Mass_gen)/Mass_svc_jecscaleDown;
                        
                        if(era == "UL2018")
                        {
                            ClassicSVfit svFitAlgo_C_HEM2018Up = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_HEM2018ScaleUp.X(), metLV_HEM2018ScaleUp.Y(), covMET, inputFile_visPtResolution);
                            m_svc_HEM2018scaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getMass();
                            pt_svc_HEM2018scaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getPt();
                            eta_svc_HEM2018scaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getEta();
                            phi_svc_HEM2018scaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getPhi();
                            mt_svc_HEM2018scaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getTransverseMass();

                            TLorentzVector ditau_svc_HEM2018scaleUp;
                            ditau_svc_HEM2018scaleUp.SetPtEtaPhiM(pt_svc_HEM2018scaleUp, eta_svc_HEM2018scaleUp, phi_svc_HEM2018scaleUp, m_svc_HEM2018scaleUp);
                            Mass_svc_HEM2018scaleUp = (muon1+muon2+ditau_svc_HEM2018scaleUp).M();
                            res_svc_HEM2018scaleUp = (Mass_svc_HEM2018scaleUp - Mass_gen)/Mass_svc_HEM2018scaleUp;
                             
                            ClassicSVfit svFitAlgo_C_HEM2018Down = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_HEM2018ScaleDown.X(), metLV_HEM2018ScaleDown.Y(), covMET, inputFile_visPtResolution);
                            m_svc_HEM2018scaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getMass();
                            pt_svc_HEM2018scaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getPt();
                            eta_svc_HEM2018scaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getEta();
                            phi_svc_HEM2018scaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getPhi();
                            mt_svc_HEM2018scaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getTransverseMass();

                            TLorentzVector ditau_svc_HEM2018scaleDown;
                            ditau_svc_HEM2018scaleDown.SetPtEtaPhiM(pt_svc_HEM2018scaleDown, eta_svc_HEM2018scaleDown, phi_svc_HEM2018scaleDown, m_svc_HEM2018scaleDown);
                            Mass_svc_HEM2018scaleDown = (muon1+muon2+ditau_svc_HEM2018scaleDown).M();
                            res_svc_HEM2018scaleDown = (Mass_svc_HEM2018scaleDown - Mass_gen)/Mass_svc_HEM2018scaleDown;
                        }
                    }
                    
                    if(HiggsFinalState == "ET")
                    {
                        ClassicSVfit svFitAlgo_C_ElectronESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_ElectronES_Up, metLV_electronScaleUp.X(), metLV_electronScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_escaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESUp.getHistogramAdapter())->getMass();
                        pt_svc_escaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESUp.getHistogramAdapter())->getPt();
                        eta_svc_escaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESUp.getHistogramAdapter())->getEta();
                        phi_svc_escaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESUp.getHistogramAdapter())->getPhi();
                        mt_svc_escaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_escaleUp;
                        ditau_svc_escaleUp.SetPtEtaPhiM(pt_svc_escaleUp, eta_svc_escaleUp, phi_svc_escaleUp, m_svc_escaleUp);
                        Mass_svc_escaleUp = (muon1+muon2+ditau_svc_escaleUp).M();
                        res_svc_escaleUp = (Mass_svc_escaleUp - Mass_gen)/Mass_svc_escaleUp;
                        
                        ClassicSVfit svFitAlgo_C_ElectronESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_ElectronES_Down, metLV_electronScaleDown.X(), metLV_electronScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_escaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESDown.getHistogramAdapter())->getMass();
                        pt_svc_escaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESDown.getHistogramAdapter())->getPt();
                        eta_svc_escaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESDown.getHistogramAdapter())->getEta();
                        phi_svc_escaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESDown.getHistogramAdapter())->getPhi();
                        mt_svc_escaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_escaleDown;
                        ditau_svc_escaleDown.SetPtEtaPhiM(pt_svc_escaleDown, eta_svc_escaleDown, phi_svc_escaleDown, m_svc_escaleDown);
                        Mass_svc_escaleDown = (muon1+muon2+ditau_svc_escaleDown).M();
                        res_svc_escaleDown = (Mass_svc_escaleDown - Mass_gen)/Mass_svc_escaleDown;
                        
                        Mass_svc_mscaleUp = (muon1UpLV+muon2UpLV+ditau_svc).M();
                        res_svc_mscaleUp = (Mass_svc_mscaleUp - Mass_gen)/Mass_svc_mscaleUp;
                        Mass_svc_mscaleDown = (muon1DownLV+muon2DownLV+ditau_svc).M();
                        res_svc_mscaleDown = (Mass_svc_mscaleDown - Mass_gen)/Mass_svc_mscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_TauESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_TauES_Up, metLV_tauScaleUp.X(), metLV_tauScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_tscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getMass();
                        pt_svc_tscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getPt();
                        eta_svc_tscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getEta();
                        phi_svc_tscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getPhi();
                        mt_svc_tscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_tscaleUp;
                        ditau_svc_tscaleUp.SetPtEtaPhiM(pt_svc_tscaleUp, eta_svc_tscaleUp, phi_svc_tscaleUp, m_svc_tscaleUp);
                        Mass_svc_tscaleUp = (muon1+muon2+ditau_svc_tscaleUp).M();
                        res_svc_tscaleUp = (Mass_svc_tscaleUp - Mass_gen)/Mass_svc_tscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_TauESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_TauES_Down, metLV_tauScaleDown.X(), metLV_tauScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_tscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getMass();
                        pt_svc_tscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getPt();
                        eta_svc_tscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getEta();
                        phi_svc_tscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getPhi();
                        mt_svc_tscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_tscaleDown;
                        ditau_svc_tscaleDown.SetPtEtaPhiM(pt_svc_tscaleDown, eta_svc_tscaleDown, phi_svc_tscaleDown, m_svc_tscaleDown);
                        Mass_svc_tscaleDown = (muon1+muon2+ditau_svc_tscaleDown).M();
                        res_svc_tscaleDown = (Mass_svc_tscaleDown - Mass_gen)/Mass_svc_tscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_ElectronFakingTauESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_ElectronFakingTauES_Up, metLV_electronfakingtauScaleUp.X(), metLV_electronfakingtauScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_eftscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getMass();
                        pt_svc_eftscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getPt();
                        eta_svc_eftscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getEta();
                        phi_svc_eftscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getPhi();
                        mt_svc_eftscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_eftscaleUp;
                        ditau_svc_eftscaleUp.SetPtEtaPhiM(pt_svc_eftscaleUp, eta_svc_eftscaleUp, phi_svc_eftscaleUp, m_svc_eftscaleUp);
                        Mass_svc_eftscaleUp = (muon1+muon2+ditau_svc_eftscaleUp).M();
                        res_svc_eftscaleUp = (Mass_svc_eftscaleUp - Mass_gen)/Mass_svc_eftscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_ElectronFakingTauESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_ElectronFakingTauES_Down, metLV_electronfakingtauScaleDown.X(), metLV_electronfakingtauScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_eftscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getMass();
                        pt_svc_eftscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getPt();
                        eta_svc_eftscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getEta();
                        phi_svc_eftscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getPhi();
                        mt_svc_eftscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_eftscaleDown;
                        ditau_svc_eftscaleDown.SetPtEtaPhiM(pt_svc_eftscaleDown, eta_svc_eftscaleDown, phi_svc_eftscaleDown, m_svc_eftscaleDown);
                        Mass_svc_eftscaleDown = (muon1+muon2+ditau_svc_eftscaleDown).M();
                        res_svc_eftscaleDown = (Mass_svc_eftscaleDown - Mass_gen)/Mass_svc_eftscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_MuonFakingTauESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_MuonFakingTauES_Up, metLV_muonfakingtauScaleUp.X(), metLV_muonfakingtauScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_mftscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getMass();
                        pt_svc_mftscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getPt();
                        eta_svc_mftscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getEta();
                        phi_svc_mftscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getPhi();
                        mt_svc_mftscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_mftscaleUp;
                        ditau_svc_mftscaleUp.SetPtEtaPhiM(pt_svc_mftscaleUp, eta_svc_mftscaleUp, phi_svc_mftscaleUp, m_svc_mftscaleUp);
                        Mass_svc_mftscaleUp = (muon1+muon2+ditau_svc_mftscaleUp).M();
                        res_svc_mftscaleUp = (Mass_svc_mftscaleUp - Mass_gen)/Mass_svc_mftscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_MuonFakingTauESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_MuonFakingTauES_Down, metLV_muonfakingtauScaleDown.X(), metLV_muonfakingtauScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_mftscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getMass();
                        pt_svc_mftscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getPt();
                        eta_svc_mftscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getEta();
                        phi_svc_mftscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getPhi();
                        mt_svc_mftscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_mftscaleDown;
                        ditau_svc_mftscaleDown.SetPtEtaPhiM(pt_svc_mftscaleDown, eta_svc_mftscaleDown, phi_svc_mftscaleDown, m_svc_mftscaleDown);
                        Mass_svc_mftscaleDown = (muon1+muon2+ditau_svc_mftscaleDown).M();
                        res_svc_mftscaleDown = (Mass_svc_mftscaleDown - Mass_gen)/Mass_svc_mftscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_UnclMETESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_unclmetScaleUp.X(), metLV_unclmetScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_unclmetscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getMass();
                        pt_svc_unclmetscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getPt();
                        eta_svc_unclmetscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getEta();
                        phi_svc_unclmetscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getPhi();
                        mt_svc_unclmetscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_unclmetscaleUp;
                        ditau_svc_unclmetscaleUp.SetPtEtaPhiM(pt_svc_unclmetscaleUp, eta_svc_unclmetscaleUp, phi_svc_unclmetscaleUp, m_svc_unclmetscaleUp);
                        Mass_svc_unclmetscaleUp = (muon1+muon2+ditau_svc_unclmetscaleUp).M();
                        res_svc_unclmetscaleUp = (Mass_svc_unclmetscaleUp - Mass_gen)/Mass_svc_unclmetscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_UnclMETESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_unclmetScaleDown.X(), metLV_unclmetScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_unclmetscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getMass();
                        pt_svc_unclmetscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getPt();
                        eta_svc_unclmetscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getEta();
                        phi_svc_unclmetscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getPhi();
                        mt_svc_unclmetscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_unclmetscaleDown;
                        ditau_svc_unclmetscaleDown.SetPtEtaPhiM(pt_svc_unclmetscaleDown, eta_svc_unclmetscaleDown, phi_svc_unclmetscaleDown, m_svc_unclmetscaleDown);
                        Mass_svc_unclmetscaleDown = (muon1+muon2+ditau_svc_unclmetscaleDown).M();
                        res_svc_unclmetscaleDown = (Mass_svc_unclmetscaleDown - Mass_gen)/Mass_svc_unclmetscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_JESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_jecScaleUp.X(), metLV_jecScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_jecscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getMass();
                        pt_svc_jecscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getPt();
                        eta_svc_jecscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getEta();
                        phi_svc_jecscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getPhi();
                        mt_svc_jecscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_jecscaleUp;
                        ditau_svc_jecscaleUp.SetPtEtaPhiM(pt_svc_jecscaleUp, eta_svc_jecscaleUp, phi_svc_jecscaleUp, m_svc_jecscaleUp);
                        Mass_svc_jecscaleUp = (muon1+muon2+ditau_svc_jecscaleUp).M();
                        res_svc_jecscaleUp = (Mass_svc_jecscaleUp - Mass_gen)/Mass_svc_jecscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_JESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_jecScaleDown.X(), metLV_jecScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_jecscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getMass();
                        pt_svc_jecscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getPt();
                        eta_svc_jecscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getEta();
                        phi_svc_jecscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getPhi();
                        mt_svc_jecscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_jecscaleDown;
                        ditau_svc_jecscaleDown.SetPtEtaPhiM(pt_svc_jecscaleDown, eta_svc_jecscaleDown, phi_svc_jecscaleDown, m_svc_jecscaleDown);
                        Mass_svc_jecscaleDown = (muon1+muon2+ditau_svc_jecscaleDown).M();
                        res_svc_jecscaleDown = (Mass_svc_jecscaleDown - Mass_gen)/Mass_svc_jecscaleDown;
                        
                        if(era == "UL2018")
                        {
                            ClassicSVfit svFitAlgo_C_HEM2018Up = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_HEM2018ScaleUp.X(), metLV_HEM2018ScaleUp.Y(), covMET, inputFile_visPtResolution);
                            m_svc_HEM2018scaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getMass();
                            pt_svc_HEM2018scaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getPt();
                            eta_svc_HEM2018scaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getEta();
                            phi_svc_HEM2018scaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getPhi();
                            mt_svc_HEM2018scaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getTransverseMass();

                            TLorentzVector ditau_svc_HEM2018scaleUp;
                            ditau_svc_HEM2018scaleUp.SetPtEtaPhiM(pt_svc_HEM2018scaleUp, eta_svc_HEM2018scaleUp, phi_svc_HEM2018scaleUp, m_svc_HEM2018scaleUp);
                            Mass_svc_HEM2018scaleUp = (muon1+muon2+ditau_svc_HEM2018scaleUp).M();
                            res_svc_HEM2018scaleUp = (Mass_svc_HEM2018scaleUp - Mass_gen)/Mass_svc_HEM2018scaleUp;
                             
                            ClassicSVfit svFitAlgo_C_HEM2018Down = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_HEM2018ScaleDown.X(), metLV_HEM2018ScaleDown.Y(), covMET, inputFile_visPtResolution);
                            m_svc_HEM2018scaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getMass();
                            pt_svc_HEM2018scaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getPt();
                            eta_svc_HEM2018scaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getEta();
                            phi_svc_HEM2018scaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getPhi();
                            mt_svc_HEM2018scaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getTransverseMass();

                            TLorentzVector ditau_svc_HEM2018scaleDown;
                            ditau_svc_HEM2018scaleDown.SetPtEtaPhiM(pt_svc_HEM2018scaleDown, eta_svc_HEM2018scaleDown, phi_svc_HEM2018scaleDown, m_svc_HEM2018scaleDown);
                            Mass_svc_HEM2018scaleDown = (muon1+muon2+ditau_svc_HEM2018scaleDown).M();
                            res_svc_HEM2018scaleDown = (Mass_svc_HEM2018scaleDown - Mass_gen)/Mass_svc_HEM2018scaleDown;
                        }
                    }
                    
                    if(HiggsFinalState == "MT")
                    {
                        ClassicSVfit svFitAlgo_C_MuonESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_MuonES_Up, metLV_muonScaleUp.X(), metLV_muonScaleUp.Y(), covMET, inputFile_visPtResolution);
                        ClassicSVfit svFitAlgo_C_MuonESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_MuonES_Down, metLV_muonScaleDown.X(), metLV_muonScaleDown.Y(), covMET, inputFile_visPtResolution);
                        
                        m_svc_mscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESUp.getHistogramAdapter())->getMass();
                        pt_svc_mscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESUp.getHistogramAdapter())->getPt();
                        eta_svc_mscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESUp.getHistogramAdapter())->getEta();
                        phi_svc_mscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESUp.getHistogramAdapter())->getPhi();
                        mt_svc_mscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_mscaleUp;
                        ditau_svc_mscaleUp.SetPtEtaPhiM(pt_svc_mscaleUp, eta_svc_mscaleUp, phi_svc_mscaleUp, m_svc_mscaleUp);
                        Mass_svc_mscaleUp = (muon1UpLV+muon2UpLV+ditau_svc_mscaleUp).M();
                        res_svc_mscaleUp = (Mass_svc_mscaleUp - Mass_gen)/Mass_svc_mscaleUp;
                        
                        m_svc_mscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESDown.getHistogramAdapter())->getMass();
                        pt_svc_mscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESDown.getHistogramAdapter())->getPt();
                        eta_svc_mscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESDown.getHistogramAdapter())->getEta();
                        phi_svc_mscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESDown.getHistogramAdapter())->getPhi();
                        mt_svc_mscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_mscaleDown;
                        ditau_svc_mscaleDown.SetPtEtaPhiM(pt_svc_mscaleDown, eta_svc_mscaleDown, phi_svc_mscaleDown, m_svc_mscaleDown);
                        Mass_svc_mscaleDown = (muon1DownLV+muon2DownLV+ditau_svc_mscaleDown).M();
                        res_svc_mscaleDown = (Mass_svc_mscaleDown - Mass_gen)/Mass_svc_mscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_TauESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_TauES_Up, metLV_tauScaleUp.X(), metLV_tauScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_tscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getMass();
                        pt_svc_tscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getPt();
                        eta_svc_tscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getEta();
                        phi_svc_tscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getPhi();
                        mt_svc_tscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_tscaleUp;
                        ditau_svc_tscaleUp.SetPtEtaPhiM(pt_svc_tscaleUp, eta_svc_tscaleUp, phi_svc_tscaleUp, m_svc_tscaleUp);
                        Mass_svc_tscaleUp = (muon1+muon2+ditau_svc_tscaleUp).M();
                        res_svc_tscaleUp = (Mass_svc_tscaleUp - Mass_gen)/Mass_svc_tscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_TauESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_TauES_Down, metLV_tauScaleDown.X(), metLV_tauScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_tscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getMass();
                        pt_svc_tscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getPt();
                        eta_svc_tscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getEta();
                        phi_svc_tscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getPhi();
                        mt_svc_tscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_tscaleDown;
                        ditau_svc_tscaleDown.SetPtEtaPhiM(pt_svc_tscaleDown, eta_svc_tscaleDown, phi_svc_tscaleDown, m_svc_tscaleDown);
                        Mass_svc_tscaleDown = (muon1+muon2+ditau_svc_tscaleDown).M();
                        res_svc_tscaleDown = (Mass_svc_tscaleDown - Mass_gen)/Mass_svc_tscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_ElectronFakingTauESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_ElectronFakingTauES_Up, metLV_electronfakingtauScaleUp.X(), metLV_electronfakingtauScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_eftscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getMass();
                        pt_svc_eftscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getPt();
                        eta_svc_eftscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getEta();
                        phi_svc_eftscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getPhi();
                        mt_svc_eftscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_eftscaleUp;
                        ditau_svc_eftscaleUp.SetPtEtaPhiM(pt_svc_eftscaleUp, eta_svc_eftscaleUp, phi_svc_eftscaleUp, m_svc_eftscaleUp);
                        Mass_svc_eftscaleUp = (muon1+muon2+ditau_svc_eftscaleUp).M();
                        res_svc_eftscaleUp = (Mass_svc_eftscaleUp - Mass_gen)/Mass_svc_eftscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_ElectronFakingTauESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_ElectronFakingTauES_Down, metLV_electronfakingtauScaleDown.X(), metLV_electronfakingtauScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_eftscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getMass();
                        pt_svc_eftscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getPt();
                        eta_svc_eftscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getEta();
                        phi_svc_eftscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getPhi();
                        mt_svc_eftscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_eftscaleDown;
                        ditau_svc_eftscaleDown.SetPtEtaPhiM(pt_svc_eftscaleDown, eta_svc_eftscaleDown, phi_svc_eftscaleDown, m_svc_eftscaleDown);
                        Mass_svc_eftscaleDown = (muon1+muon2+ditau_svc_eftscaleDown).M();
                        res_svc_eftscaleDown = (Mass_svc_eftscaleDown - Mass_gen)/Mass_svc_eftscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_MuonFakingTauESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_MuonFakingTauES_Up, metLV_muonfakingtauScaleUp.X(), metLV_muonfakingtauScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_mftscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getMass();
                        pt_svc_mftscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getPt();
                        eta_svc_mftscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getEta();
                        phi_svc_mftscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getPhi();
                        mt_svc_mftscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_mftscaleUp;
                        ditau_svc_mftscaleUp.SetPtEtaPhiM(pt_svc_mftscaleUp, eta_svc_mftscaleUp, phi_svc_mftscaleUp, m_svc_mftscaleUp);
                        Mass_svc_mftscaleUp = (muon1+muon2+ditau_svc_mftscaleUp).M();
                        res_svc_mftscaleUp = (Mass_svc_mftscaleUp - Mass_gen)/Mass_svc_mftscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_MuonFakingTauESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_MuonFakingTauES_Down, metLV_muonfakingtauScaleDown.X(), metLV_muonfakingtauScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_mftscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getMass();
                        pt_svc_mftscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getPt();
                        eta_svc_mftscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getEta();
                        phi_svc_mftscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getPhi();
                        mt_svc_mftscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_mftscaleDown;
                        ditau_svc_mftscaleDown.SetPtEtaPhiM(pt_svc_mftscaleDown, eta_svc_mftscaleDown, phi_svc_mftscaleDown, m_svc_mftscaleDown);
                        Mass_svc_mftscaleDown = (muon1+muon2+ditau_svc_mftscaleDown).M();
                        res_svc_mftscaleDown = (Mass_svc_mftscaleDown - Mass_gen)/Mass_svc_mftscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_UnclMETESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_unclmetScaleUp.X(), metLV_unclmetScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_unclmetscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getMass();
                        pt_svc_unclmetscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getPt();
                        eta_svc_unclmetscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getEta();
                        phi_svc_unclmetscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getPhi();
                        mt_svc_unclmetscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_unclmetscaleUp;
                        ditau_svc_unclmetscaleUp.SetPtEtaPhiM(pt_svc_unclmetscaleUp, eta_svc_unclmetscaleUp, phi_svc_unclmetscaleUp, m_svc_unclmetscaleUp);
                        Mass_svc_unclmetscaleUp = (muon1+muon2+ditau_svc_unclmetscaleUp).M();
                        res_svc_unclmetscaleUp = (Mass_svc_unclmetscaleUp - Mass_gen)/Mass_svc_unclmetscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_UnclMETESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_unclmetScaleDown.X(), metLV_unclmetScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_unclmetscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getMass();
                        pt_svc_unclmetscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getPt();
                        eta_svc_unclmetscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getEta();
                        phi_svc_unclmetscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getPhi();
                        mt_svc_unclmetscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_unclmetscaleDown;
                        ditau_svc_unclmetscaleDown.SetPtEtaPhiM(pt_svc_unclmetscaleDown, eta_svc_unclmetscaleDown, phi_svc_unclmetscaleDown, m_svc_unclmetscaleDown);
                        Mass_svc_unclmetscaleDown = (muon1+muon2+ditau_svc_unclmetscaleDown).M();
                        res_svc_unclmetscaleDown = (Mass_svc_unclmetscaleDown - Mass_gen)/Mass_svc_unclmetscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_JESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_jecScaleUp.X(), metLV_jecScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_jecscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getMass();
                        pt_svc_jecscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getPt();
                        eta_svc_jecscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getEta();
                        phi_svc_jecscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getPhi();
                        mt_svc_jecscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_jecscaleUp;
                        ditau_svc_jecscaleUp.SetPtEtaPhiM(pt_svc_jecscaleUp, eta_svc_jecscaleUp, phi_svc_jecscaleUp, m_svc_jecscaleUp);
                        Mass_svc_jecscaleUp = (muon1+muon2+ditau_svc_jecscaleUp).M();
                        res_svc_jecscaleUp = (Mass_svc_jecscaleUp - Mass_gen)/Mass_svc_jecscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_JESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_jecScaleDown.X(), metLV_jecScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_jecscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getMass();
                        pt_svc_jecscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getPt();
                        eta_svc_jecscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getEta();
                        phi_svc_jecscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getPhi();
                        mt_svc_jecscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_jecscaleDown;
                        ditau_svc_jecscaleDown.SetPtEtaPhiM(pt_svc_jecscaleDown, eta_svc_jecscaleDown, phi_svc_jecscaleDown, m_svc_jecscaleDown);
                        Mass_svc_jecscaleDown = (muon1+muon2+ditau_svc_jecscaleDown).M();
                        res_svc_jecscaleDown = (Mass_svc_jecscaleDown - Mass_gen)/Mass_svc_jecscaleDown;
                        
                        if(era == "UL2018")
                        {
                            ClassicSVfit svFitAlgo_C_HEM2018Up = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_HEM2018ScaleUp.X(), metLV_HEM2018ScaleUp.Y(), covMET, inputFile_visPtResolution);
                            m_svc_HEM2018scaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getMass();
                            pt_svc_HEM2018scaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getPt();
                            eta_svc_HEM2018scaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getEta();
                            phi_svc_HEM2018scaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getPhi();
                            mt_svc_HEM2018scaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getTransverseMass();

                            TLorentzVector ditau_svc_HEM2018scaleUp;
                            ditau_svc_HEM2018scaleUp.SetPtEtaPhiM(pt_svc_HEM2018scaleUp, eta_svc_HEM2018scaleUp, phi_svc_HEM2018scaleUp, m_svc_HEM2018scaleUp);
                            Mass_svc_HEM2018scaleUp = (muon1+muon2+ditau_svc_HEM2018scaleUp).M();
                            res_svc_HEM2018scaleUp = (Mass_svc_HEM2018scaleUp - Mass_gen)/Mass_svc_HEM2018scaleUp;
                             
                            ClassicSVfit svFitAlgo_C_HEM2018Down = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_HEM2018ScaleDown.X(), metLV_HEM2018ScaleDown.Y(), covMET, inputFile_visPtResolution);
                            m_svc_HEM2018scaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getMass();
                            pt_svc_HEM2018scaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getPt();
                            eta_svc_HEM2018scaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getEta();
                            phi_svc_HEM2018scaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getPhi();
                            mt_svc_HEM2018scaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getTransverseMass();

                            TLorentzVector ditau_svc_HEM2018scaleDown;
                            ditau_svc_HEM2018scaleDown.SetPtEtaPhiM(pt_svc_HEM2018scaleDown, eta_svc_HEM2018scaleDown, phi_svc_HEM2018scaleDown, m_svc_HEM2018scaleDown);
                            Mass_svc_HEM2018scaleDown = (muon1+muon2+ditau_svc_HEM2018scaleDown).M();
                            res_svc_HEM2018scaleDown = (Mass_svc_HEM2018scaleDown - Mass_gen)/Mass_svc_HEM2018scaleDown;
                        }
                    }
                    
                    if(HiggsFinalState == "TT")
                    {
                        Mass_svc_mscaleUp = (muon1UpLV+muon2UpLV+ditau_svc).M();
                        res_svc_mscaleUp = (Mass_svc_mscaleUp - Mass_gen)/Mass_svc_mscaleUp;
                        Mass_svc_mscaleDown = (muon1DownLV+muon2DownLV+ditau_svc).M();
                        res_svc_mscaleDown = (Mass_svc_mscaleDown - Mass_gen)/Mass_svc_mscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_TauESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_TauES_Up, metLV_tauScaleUp.X(), metLV_tauScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_tscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getMass();
                        pt_svc_tscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getPt();
                        eta_svc_tscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getEta();
                        phi_svc_tscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getPhi();
                        mt_svc_tscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_tscaleUp;
                        ditau_svc_tscaleUp.SetPtEtaPhiM(pt_svc_tscaleUp, eta_svc_tscaleUp, phi_svc_tscaleUp, m_svc_tscaleUp);
                        Mass_svc_tscaleUp = (muon1+muon2+ditau_svc_tscaleUp).M();
                        res_svc_tscaleUp = (Mass_svc_tscaleUp - Mass_gen)/Mass_svc_tscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_TauESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_TauES_Down, metLV_tauScaleDown.X(), metLV_tauScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_tscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getMass();
                        pt_svc_tscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getPt();
                        eta_svc_tscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getEta();
                        phi_svc_tscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getPhi();
                        mt_svc_tscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_TauESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_tscaleDown;
                        ditau_svc_tscaleDown.SetPtEtaPhiM(pt_svc_tscaleDown, eta_svc_tscaleDown, phi_svc_tscaleDown, m_svc_tscaleDown);
                        Mass_svc_tscaleDown = (muon1+muon2+ditau_svc_tscaleDown).M();
                        res_svc_tscaleDown = (Mass_svc_tscaleDown - Mass_gen)/Mass_svc_tscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_ElectronFakingTauESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_ElectronFakingTauES_Up, metLV_electronfakingtauScaleUp.X(), metLV_electronfakingtauScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_eftscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getMass();
                        pt_svc_eftscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getPt();
                        eta_svc_eftscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getEta();
                        phi_svc_eftscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getPhi();
                        mt_svc_eftscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_eftscaleUp;
                        ditau_svc_eftscaleUp.SetPtEtaPhiM(pt_svc_eftscaleUp, eta_svc_eftscaleUp, phi_svc_eftscaleUp, m_svc_eftscaleUp);
                        Mass_svc_eftscaleUp = (muon1+muon2+ditau_svc_eftscaleUp).M();
                        res_svc_eftscaleUp = (Mass_svc_eftscaleUp - Mass_gen)/Mass_svc_eftscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_ElectronFakingTauESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_ElectronFakingTauES_Down, metLV_electronfakingtauScaleDown.X(), metLV_electronfakingtauScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_eftscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getMass();
                        pt_svc_eftscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getPt();
                        eta_svc_eftscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getEta();
                        phi_svc_eftscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getPhi();
                        mt_svc_eftscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_ElectronFakingTauESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_eftscaleDown;
                        ditau_svc_eftscaleDown.SetPtEtaPhiM(pt_svc_eftscaleDown, eta_svc_eftscaleDown, phi_svc_eftscaleDown, m_svc_eftscaleDown);
                        Mass_svc_eftscaleDown = (muon1+muon2+ditau_svc_eftscaleDown).M();
                        res_svc_eftscaleDown = (Mass_svc_eftscaleDown - Mass_gen)/Mass_svc_eftscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_MuonFakingTauESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_MuonFakingTauES_Up, metLV_muonfakingtauScaleUp.X(), metLV_muonfakingtauScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_mftscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getMass();
                        pt_svc_mftscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getPt();
                        eta_svc_mftscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getEta();
                        phi_svc_mftscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getPhi();
                        mt_svc_mftscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_mftscaleUp;
                        ditau_svc_mftscaleUp.SetPtEtaPhiM(pt_svc_mftscaleUp, eta_svc_mftscaleUp, phi_svc_mftscaleUp, m_svc_mftscaleUp);
                        Mass_svc_mftscaleUp = (muon1+muon2+ditau_svc_mftscaleUp).M();
                        res_svc_mftscaleUp = (Mass_svc_mftscaleUp - Mass_gen)/Mass_svc_mftscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_MuonFakingTauESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons_MuonFakingTauES_Down, metLV_muonfakingtauScaleDown.X(), metLV_muonfakingtauScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_mftscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getMass();
                        pt_svc_mftscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getPt();
                        eta_svc_mftscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getEta();
                        phi_svc_mftscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getPhi();
                        mt_svc_mftscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_MuonFakingTauESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_mftscaleDown;
                        ditau_svc_mftscaleDown.SetPtEtaPhiM(pt_svc_mftscaleDown, eta_svc_mftscaleDown, phi_svc_mftscaleDown, m_svc_mftscaleDown);
                        Mass_svc_mftscaleDown = (muon1+muon2+ditau_svc_mftscaleDown).M();
                        res_svc_mftscaleDown = (Mass_svc_mftscaleDown - Mass_gen)/Mass_svc_mftscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_UnclMETESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_unclmetScaleUp.X(), metLV_unclmetScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_unclmetscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getMass();
                        pt_svc_unclmetscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getPt();
                        eta_svc_unclmetscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getEta();
                        phi_svc_unclmetscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getPhi();
                        mt_svc_unclmetscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_unclmetscaleUp;
                        ditau_svc_unclmetscaleUp.SetPtEtaPhiM(pt_svc_unclmetscaleUp, eta_svc_unclmetscaleUp, phi_svc_unclmetscaleUp, m_svc_unclmetscaleUp);
                        Mass_svc_unclmetscaleUp = (muon1+muon2+ditau_svc_unclmetscaleUp).M();
                        res_svc_unclmetscaleUp = (Mass_svc_unclmetscaleUp - Mass_gen)/Mass_svc_unclmetscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_UnclMETESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_unclmetScaleDown.X(), metLV_unclmetScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_unclmetscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getMass();
                        pt_svc_unclmetscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getPt();
                        eta_svc_unclmetscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getEta();
                        phi_svc_unclmetscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getPhi();
                        mt_svc_unclmetscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_unclmetscaleDown;
                        ditau_svc_unclmetscaleDown.SetPtEtaPhiM(pt_svc_unclmetscaleDown, eta_svc_unclmetscaleDown, phi_svc_unclmetscaleDown, m_svc_unclmetscaleDown);
                        Mass_svc_unclmetscaleDown = (muon1+muon2+ditau_svc_unclmetscaleDown).M();
                        res_svc_unclmetscaleDown = (Mass_svc_unclmetscaleDown - Mass_gen)/Mass_svc_unclmetscaleDown;
                        
                        ClassicSVfit svFitAlgo_C_JESUp = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_jecScaleUp.X(), metLV_jecScaleUp.Y(), covMET, inputFile_visPtResolution);
                        m_svc_jecscaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getMass();
                        pt_svc_jecscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getPt();
                        eta_svc_jecscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getEta();
                        phi_svc_jecscaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getPhi();
                        mt_svc_jecscaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESUp.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_jecscaleUp;
                        ditau_svc_jecscaleUp.SetPtEtaPhiM(pt_svc_jecscaleUp, eta_svc_jecscaleUp, phi_svc_jecscaleUp, m_svc_jecscaleUp);
                        Mass_svc_jecscaleUp = (muon1+muon2+ditau_svc_jecscaleUp).M();
                        res_svc_jecscaleUp = (Mass_svc_jecscaleUp - Mass_gen)/Mass_svc_jecscaleUp;
                        
                        ClassicSVfit svFitAlgo_C_JESDown = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_jecScaleDown.X(), metLV_jecScaleDown.Y(), covMET, inputFile_visPtResolution);
                        m_svc_jecscaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getMass();
                        pt_svc_jecscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getPt();
                        eta_svc_jecscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getEta();
                        phi_svc_jecscaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_JESDown.getHistogramAdapter())->getPhi();
                        mt_svc_jecscaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_UnclMETESDown.getHistogramAdapter())->getTransverseMass();

                        TLorentzVector ditau_svc_jecscaleDown;
                        ditau_svc_jecscaleDown.SetPtEtaPhiM(pt_svc_jecscaleDown, eta_svc_jecscaleDown, phi_svc_jecscaleDown, m_svc_jecscaleDown);
                        Mass_svc_jecscaleDown = (muon1+muon2+ditau_svc_jecscaleDown).M();
                        res_svc_jecscaleDown = (Mass_svc_jecscaleDown - Mass_gen)/Mass_svc_jecscaleDown;
                        
                        if(era == "UL2018")
                        {
                            ClassicSVfit svFitAlgo_C_HEM2018Up = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_HEM2018ScaleUp.X(), metLV_HEM2018ScaleUp.Y(), covMET, inputFile_visPtResolution);
                            m_svc_HEM2018scaleUp   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getMass();
                            pt_svc_HEM2018scaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getPt();
                            eta_svc_HEM2018scaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getEta();
                            phi_svc_HEM2018scaleUp = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getPhi();
                            mt_svc_HEM2018scaleUp  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Up.getHistogramAdapter())->getTransverseMass();

                            TLorentzVector ditau_svc_HEM2018scaleUp;
                            ditau_svc_HEM2018scaleUp.SetPtEtaPhiM(pt_svc_HEM2018scaleUp, eta_svc_HEM2018scaleUp, phi_svc_HEM2018scaleUp, m_svc_HEM2018scaleUp);
                            Mass_svc_HEM2018scaleUp = (muon1+muon2+ditau_svc_HEM2018scaleUp).M();
                            res_svc_HEM2018scaleUp = (Mass_svc_HEM2018scaleUp - Mass_gen)/Mass_svc_HEM2018scaleUp;
                             
                            ClassicSVfit svFitAlgo_C_HEM2018Down = SVFitMassConstraintComputation(HiggsFinalState, measuredTauLeptons, metLV_HEM2018ScaleDown.X(), metLV_HEM2018ScaleDown.Y(), covMET, inputFile_visPtResolution);
                            m_svc_HEM2018scaleDown   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getMass();
                            pt_svc_HEM2018scaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getPt();
                            eta_svc_HEM2018scaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getEta();
                            phi_svc_HEM2018scaleDown = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getPhi();
                            mt_svc_HEM2018scaleDown  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C_HEM2018Down.getHistogramAdapter())->getTransverseMass();

                            TLorentzVector ditau_svc_HEM2018scaleDown;
                            ditau_svc_HEM2018scaleDown.SetPtEtaPhiM(pt_svc_HEM2018scaleDown, eta_svc_HEM2018scaleDown, phi_svc_HEM2018scaleDown, m_svc_HEM2018scaleDown);
                            Mass_svc_HEM2018scaleDown = (muon1+muon2+ditau_svc_HEM2018scaleDown).M();
                            res_svc_HEM2018scaleDown = (Mass_svc_HEM2018scaleDown - Mass_gen)/Mass_svc_HEM2018scaleDown;
                        }
                    }
                }
            }
            if (ApplyFastMTT)
            {
                // FasMTT
                LorentzVector tau1P4;
                LorentzVector tau2P4;
                LorentzVector ttP4 = FastMTTComputation(measuredTauLeptons, measuredMETx, measuredMETy, covMET, tau1P4, tau2P4);

                double dPhiTT = dPhiFrom2P( tau1P4.Px(), tau1P4.Py(), tau2P4.Px(), tau2P4.Py() );
                mt_fast = TMath::Sqrt(2*tau1P4.Pt()*tau2P4.Pt()*(1 - TMath::Cos(dPhiTT)));
                m_fast = ttP4.M();
                pt_fast = ttP4.Pt();
                eta_fast = ttP4.Eta();
                phi_fast = ttP4.Phi();
                
                TLorentzVector ditau_fast;
                ditau_fast.SetPtEtaPhiM(pt_fast, eta_fast, phi_fast, m_fast);
                Mass_fast = (muon1+muon2+ditau_fast).M();
                res_fast = (Mass_fast - Mass_gen)/Mass_fast;

                // FasMTT mass constraint
                LorentzVector tau1P4_c;
                LorentzVector tau2P4_c;
                LorentzVector ttP4_c = FastMTTConstraintComputation(measuredTauLeptons, measuredMETx, measuredMETy, covMET, tau1P4_c, tau2P4_c);

                double dPhiTT_c = dPhiFrom2P( tau1P4_c.Px(), tau1P4_c.Py(), tau2P4_c.Px(), tau2P4_c.Py() );
                mt_fastc = TMath::Sqrt(2*tau1P4_c.Pt()*tau2P4_c.Pt()*(1 - TMath::Cos(dPhiTT_c)));
                m_fastc = ttP4_c.M();
                pt_fastc = ttP4_c.Pt();
                eta_fastc = ttP4_c.Eta();
                phi_fastc = ttP4_c.Phi();
                
                TLorentzVector ditau_fastc;
                ditau_fastc.SetPtEtaPhiM(pt_fastc, eta_fastc, phi_fastc, m_fastc);
                Mass_fastc = (muon1+muon2+ditau_fastc).M();
                res_fastc = (Mass_fastc - Mass_gen)/Mass_fastc;
            }
        }
        
        //fake rate weight calculation
        if(isData && (HiggsFinalState != "E") && (HiggsFinalState != "M") && (HiggsFinalState != "T"))
        {
            float fakerate_Z1 = fakerate_Zlep1.GetMuFakeRate(pt_1,eta_1,true);
            float fakerateUp_Z1 = fakerate_Zlep1.GetMuFakeRateUp(pt_1,eta_1,true,"68");
            float fakerateDown_Z1 = fakerate_Zlep1.GetMuFakeRateDown(pt_1,eta_1,true,"68");
            
            float fakerate_Z2 = fakerate_Zlep2.GetMuFakeRate(pt_2,eta_2,true);
            float fakerateUp_Z2 = fakerate_Zlep2.GetMuFakeRateUp(pt_2,eta_2,true,"68");
            float fakerateDown_Z2 = fakerate_Zlep2.GetMuFakeRateDown(pt_2,eta_2,true,"68");
            
            fakeRateWeight_Z1 = fakerate_Z1/(1 - fakerate_Z1);
            fakeRateWeightUp_Z1 = fakerateUp_Z1/(1 - fakerateUp_Z1);
            fakeRateWeightDown_Z1 = fakerateDown_Z1/(1 - fakerateDown_Z1);

            fakeRateWeight_Z2 = fakerate_Z2/(1 - fakerate_Z2);
            fakeRateWeightUp_Z2 = fakerateUp_Z2/(1 - fakerateUp_Z2);
            fakeRateWeightDown_Z2 = fakerateDown_Z2/(1 - fakerateDown_Z2);
        }
        
        if(isData && (HiggsFinalState == "EM"))
        {
            float fakerate_1 = fakerate_lep1.GetEleFakeRate(pt_3,eta_3,true);
            float fakerateUp_1 = fakerate_lep1.GetEleFakeRateUp(pt_3,eta_3,true,"68");
            float fakerateDown_1 = fakerate_lep1.GetEleFakeRateDown(pt_3,eta_3,true,"68");
            
            float fakerate_2 = fakerate_lep2.GetMuFakeRate(pt_4,eta_4,true);
            float fakerateUp_2 = fakerate_lep2.GetMuFakeRateUp(pt_4,eta_4,true,"68");
            float fakerateDown_2 = fakerate_lep2.GetMuFakeRateDown(pt_4,eta_4,true,"68");

            fakeRateWeight_1 = fakerate_1/(1 - fakerate_1);
            fakeRateWeightUp_1 = fakerateUp_1/(1 - fakerateUp_1);
            fakeRateWeightDown_1 = fakerateDown_1/(1 - fakerateDown_1);

            fakeRateWeight_2 = fakerate_2/(1 - fakerate_2);
            fakeRateWeightUp_2 = fakerateUp_2/(1 - fakerateUp_2);
            fakeRateWeightDown_2 = fakerateDown_2/(1 - fakerateDown_2);
        }
        if(isData && (HiggsFinalState == "ET"))
        {
            float fakerate_1 = fakerate_lep1.GetEleFakeRate(pt_3,eta_3,true);
            float fakerateUp_1 = fakerate_lep1.GetEleFakeRateUp(pt_3,eta_3,true,"68");
            float fakerateDown_1 = fakerate_lep1.GetEleFakeRateDown(pt_3,eta_3,true,"68");
            
            float fakerate_2 = fakerate_lep2.GetTauFakeRate(pt_4,dm_4,true);
            float fakerateUp_2 = fakerate_lep2.GetTauFakeRateUp(pt_4,dm_4,true,"68");
            float fakerateDown_2 = fakerate_lep2.GetTauFakeRateDown(pt_4,dm_4,true,"68");

            fakeRateWeight_1 = fakerate_1/(1 - fakerate_1);
            fakeRateWeightUp_1 = fakerateUp_1/(1 - fakerateUp_1);
            fakeRateWeightDown_1 = fakerateDown_1/(1 - fakerateDown_1);

            fakeRateWeight_2 = fakerate_2/(1 - fakerate_2);
            fakeRateWeightUp_2 = fakerateUp_2/(1 - fakerateUp_2);
            fakeRateWeightDown_2 = fakerateDown_2/(1 - fakerateDown_2);
        }
        if(isData && (HiggsFinalState == "MT"))
        {
            float fakerate_1 = fakerate_lep1.GetMuFakeRate(pt_3,eta_3,true);
            float fakerateUp_1 = fakerate_lep1.GetMuFakeRateUp(pt_3,eta_3,true,"68");
            float fakerateDown_1 = fakerate_lep1.GetMuFakeRateDown(pt_3,eta_3,true,"68");
            
            float fakerate_2 = fakerate_lep2.GetTauFakeRate(pt_4,dm_4,true);
            float fakerateUp_2 = fakerate_lep2.GetTauFakeRateUp(pt_4,dm_4,true,"68");
            float fakerateDown_2 = fakerate_lep2.GetTauFakeRateDown(pt_4,dm_4,true,"68");

            fakeRateWeight_1 = fakerate_1/(1 - fakerate_1);
            fakeRateWeightUp_1 = fakerateUp_1/(1 - fakerateUp_1);
            fakeRateWeightDown_1 = fakerateDown_1/(1 - fakerateDown_1);

            fakeRateWeight_2 = fakerate_2/(1 - fakerate_2);
            fakeRateWeightUp_2 = fakerateUp_2/(1 - fakerateUp_2);
            fakeRateWeightDown_2 = fakerateDown_2/(1 - fakerateDown_2);
        }
        if(isData && (HiggsFinalState == "TT"))
        {
            float fakerate_1 = fakerate_lep1.GetTauFakeRate(pt_3,dm_3,true);
            float fakerateUp_1 = fakerate_lep1.GetTauFakeRateUp(pt_3,dm_3,true,"68");
            float fakerateDown_1 = fakerate_lep1.GetTauFakeRateDown(pt_3,dm_3,true,"68");
            
            float fakerate_2 = fakerate_lep2.GetTauFakeRate(pt_4,dm_4,true);
            float fakerateUp_2 = fakerate_lep2.GetTauFakeRateUp(pt_4,dm_4,true,"68");
            float fakerateDown_2 = fakerate_lep2.GetTauFakeRateDown(pt_4,dm_4,true,"68");

            fakeRateWeight_1 = fakerate_1/(1 - fakerate_1);
            fakeRateWeightUp_1 = fakerateUp_1/(1 - fakerateUp_1);
            fakeRateWeightDown_1 = fakerateDown_1/(1 - fakerateDown_1);

            fakeRateWeight_2 = fakerate_2/(1 - fakerate_2);
            fakeRateWeightUp_2 = fakerateUp_2/(1 - fakerateUp_2);
            fakeRateWeightDown_2 = fakerateDown_2/(1 - fakerateDown_2);
            /*std::cout << "Fake Rate weight 1:  " << fakeRateWeight_1 << std::endl;
            std::cout << "Fake Rate weight Up 1:  " << fakeRateWeightUp_1 << std::endl;
            std::cout << "Fake Rate weight Down 1:  " << fakeRateWeightDown_1 << std::endl;
            std::cout << "pt 3:  " << pt_3 << std::endl;
            std::cout << "dm 3:  " << dm_3 << std::endl;
            std::cout << "Fake Rate weight 2:  " << fakeRateWeight_2 << std::endl;
            std::cout << "Fake Rate weight Up 2:  " << fakeRateWeightUp_2 << std::endl;
            std::cout << "Fake Rate weight Down 2:  " << fakeRateWeightDown_2 << std::endl;
            std::cout << "pt 4:  " << pt_4 << std::endl;
            std::cout << "dm 4:  " << dm_4 << std::endl;
            std::cout << "-------" << std::endl;
            std::cout << "-------" << std::endl;*/
        }
        
        //pile up weight variable
        if (!isData) {
            puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
        }
        
        // Jets selection
        vector<unsigned int> jets; jets.clear();
        vector<unsigned int> jetspt20; jetspt20.clear();
        vector<unsigned int> jetsUp; jetsUp.clear();
        vector<unsigned int> jetsDown; jetsDown.clear();
        vector<unsigned int> jetspt20Up; jetspt20Up.clear();
        vector<unsigned int> jetspt20Down; jetspt20Down.clear();
        vector<unsigned int> jetsHEM2018Up; jetsHEM2018Up.clear();
        vector<unsigned int> jetsHEM2018Down; jetsHEM2018Down.clear();
        vector<unsigned int> jetspt20HEM2018Up; jetspt20HEM2018Up.clear();
        vector<unsigned int> jetspt20HEM2018Down; jetspt20HEM2018Down.clear();
        

        vector<unsigned int> bjets; bjets.clear();
        vector<unsigned int> bjetsUp; bjetsUp.clear();
        vector<unsigned int> bjetsDown; bjetsDown.clear();
        vector<unsigned int> bjetsHEM2018Up; bjetsHEM2018Up.clear();
        vector<unsigned int> bjetsHEM2018Down; bjetsHEM2018Down.clear();
        vector<unsigned int> bjets_nocleaned; bjets_nocleaned.clear();
        vector<unsigned int> bjetsRaw; bjetsRaw.clear();
        vector<unsigned int> bjetsRawUp; bjetsRawUp.clear();
        vector<unsigned int> bjetsRawDown; bjetsRawDown.clear();
        vector<unsigned int> bjets_mistagUp;bjets_mistagUp.clear();
        vector<unsigned int> bjets_mistagDown;bjets_mistagDown.clear();
        vector<unsigned int> bjets_btagUp;bjets_btagUp.clear();
        vector<unsigned int> bjets_btagDown;bjets_btagDown.clear();
        
        float MaxBJetPt = 1000.;
        float MaxLJetPt = 1000.;
        float MinLJetPt = 20.;
        float MinBJetPt = 20.; // !!!!!
        
        int indexLeadingJet = -1;
        float ptLeadingJet = -1;
        int indexLeadingJetUp = -1;
        float ptLeadingJetUp = -1;
        int indexLeadingJetDown = -1;
        float ptLeadingJetDown = -1;
        
        int indexSubLeadingJet = -1;
        float ptSubLeadingJet = -1;
        int indexSubLeadingJetUp = -1;
        float ptSubLeadingJetUp = -1;
        int indexSubLeadingJetDown = -1;
        float ptSubLeadingJetDown = -1;
        
        int indexLeadingBJet = -1;
        float ptLeadingBJet = -1;
        int indexLeadingBJetUp = -1;
        float ptLeadingBJetUp = -1;
        int indexLeadingBJetDown = -1;
        float ptLeadingBJetDown = -1;
        
        int indexSubLeadingBJet = -1;
        float ptSubLeadingBJet = -1;
        int indexSubLeadingBJetUp = -1;
        float ptSubLeadingBJetUp = -1;
        int indexSubLeadingBJetDown = -1;
        float ptSubLeadingBJetDown = -1;
        
        for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet)
        {
            float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
            float jetEta = analysisTree.pfjet_eta[jet];
            if (absJetEta>jetEtaCut) continue;
                
            float jetPt = analysisTree.pfjet_pt[jet];

            //std::cout << jet << " : uncertainty = " << analysisTree.pfjet_jecUncertainty[jet] << std::endl;
            if (jetPt<jetPtLowCut) continue;
                
            bool isPFJetId =false;
            if (era=="2016") isPFJetId= looseJetiD_2016(analysisTree,int(jet));
            else if (era=="2017") isPFJetId = tightJetiD_2017(analysisTree,int(jet));
            else if (era=="2018") isPFJetId = tightJetiD_2018(analysisTree,int(jet));
            else if ((era=="UL2018")||(era=="UL2017")) isPFJetId = tightJetiD_UL2018UL2017(analysisTree,int(jet));
            else if ((era=="UL2016postVFP")||(era=="UL2016preVFP")) isPFJetId = tightJetiD_UL2016(analysisTree,int(jet));
            if (!isPFJetId) continue;
            
            //bool isPFJetId = tightJetiD_2018(analysisTree,int(jet));
            //if (!isPFJetId) continue;
            // cleaned jet at the moment!
            
            bool cleanedJet = true;
            float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                               eta_1,phi_1);
            if (dR1<dRJetLeptonCut) cleanedJet = false;
            float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                               eta_2,phi_2);
            if (dR2<dRJetLeptonCut) cleanedJet = false;
            if(HiggsFinalState == "EM")
            {
                float dR3 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                   eta_3,phi_3);
                if (dR3<dRJetLeptonCut) cleanedJet = false;
                float dR4 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                   eta_4,phi_4);
                if (dR4<dRJetLeptonCut) cleanedJet = false;
            }
            if((HiggsFinalState == "ET") ||(HiggsFinalState == "MT"))
            {
                float dR3 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                   eta_3,phi_3);
                if (dR3<dRJetLeptonCut) cleanedJet = false;
            }
            if (!cleanedJet) continue;
            
            // jetId
            if (jetPt>jetPtLowCut)
                jetspt20.push_back(jet);
            
            if (absJetEta<bJetEtaCut)
            {   // jet within b-tagging acceptance
                bool tagged = false;
                bool tagged_mistagUp = false;
                bool tagged_mistagDown = false;
                bool tagged_btagUp = false;
                bool tagged_btagDown = false;
                //std::cout << "D1 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant1] << std::endl;
                //std::cout << "D2 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant2] << std::endl;
                //std::cout << "D3 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant3] << std::endl;

                tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2] + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] > btagCut;
                tagged_mistagUp = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2]  + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] > btagCut;
                tagged_mistagDown = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2]  + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] > btagCut;
                tagged_btagUp = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2]  + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] > btagCut;
                tagged_btagDown = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2]  + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] > btagCut;
                bool taggedRaw = tagged;
                
                if (!isData)
                {
                    int flavor = abs(analysisTree.pfjet_flavour[jet]);
                    
                    double jet_scalefactor      = 1;
                    double jet_scalefactor_up   = 1;
                    double jet_scalefactor_down = 1;
                    double JetPtForBTag = jetPt;
                    double tageff = 1;
                    if (flavor==5)
                    {
                        if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                        if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                        jet_scalefactor      = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                        jet_scalefactor_up   = reader_BTAG.eval_auto_bounds("up"     ,BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                        jet_scalefactor_down = reader_BTAG.eval_auto_bounds("down"   ,BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                        tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
                    }
                    else if (flavor==4)
                    {
                        if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                        if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                        jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                        jet_scalefactor_up   = reader_BTAG.eval_auto_bounds("up"     , BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                        jet_scalefactor_down = reader_BTAG.eval_auto_bounds("down"   , BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                        tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
                    }
                    else
                    {
                        if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
                        if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
                        jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                        jet_scalefactor_up   = reader_BTAG.eval_auto_bounds("up"     , BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                        jet_scalefactor_down = reader_BTAG.eval_auto_bounds("down"   , BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                        tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
                    }
                    if (tageff<1e-5)      tageff = 1e-5;
                    if (tageff>0.99999)   tageff = 0.99999;
                    r.SetSeed((int)((jetEta+5)*100000));
                    double rannum = r.Rndm();
                    if (tagged)
                    { // demote
                        if(jet_scalefactor<1)
                        {
                            double fraction = 1-jet_scalefactor;
                            if (rannum<fraction) tagged = false;
                        }
                        if(jet_scalefactor_up<1)
                        {
                            double fraction_up = 1-jet_scalefactor_up;
                            if (rannum<fraction_up) tagged_mistagUp = false;
                        }
                        if(jet_scalefactor_down<1)
                        {
                            double fraction_down = 1-jet_scalefactor_down;
                            if (rannum<fraction_down) tagged_mistagDown = false;
                        }
                        tagged_btagUp   = tagged;
                        tagged_btagDown = tagged;
                    }
                    else if (!tagged)
                    { // promote
                        if(jet_scalefactor>1)
                        {
                            double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
                            if (rannum<fraction) tagged = true;
                        }
                        if(jet_scalefactor_up>1)
                        {
                            double fraction_up = (jet_scalefactor_up-1.0)/(1.0/tageff-1.0);
                            if (rannum<fraction_up) tagged_btagUp = true;
                        }
                        if(jet_scalefactor_down>1)
                        {
                            double fraction_down = (jet_scalefactor_down-1.0)/(1.0/tageff-1.0);
                            if (rannum<fraction_down) tagged_btagDown = true;
                        }
                            tagged_mistagUp   = tagged;
                            tagged_mistagDown = tagged;
                    }
                }
                //bool tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet
                //missing b-tag SFs
                if (taggedRaw) bjetsRaw.push_back(jet);

                if (tagged)
                {
                    bjets.push_back(jet);
                    if (jetPt<ptLeadingBJet&&jetPt>ptSubLeadingBJet)
                    {
                        indexSubLeadingBJet = jet;
                        ptSubLeadingBJet = jetPt;
                    }
                    if (jetPt>ptLeadingBJet)
                    {
                        ptLeadingBJet = jetPt;
                        indexLeadingBJet = jet;
                    }
                }
                if(tagged_mistagUp)   bjets_mistagUp.push_back(jet);
                if(tagged_mistagDown) bjets_mistagDown.push_back(jet);
                if(tagged_btagUp)     bjets_btagUp.push_back(jet);
                if(tagged_btagDown)   bjets_btagDown.push_back(jet);
            }

            if (jetPt>jetPtHighCut)
                jets.push_back(jet);
                
            if (indexLeadingJet>=0)
            {
                if (jetPt<ptLeadingJet&&jetPt>ptSubLeadingJet)
                {
                    indexSubLeadingJet = jet;
                    ptSubLeadingJet = jetPt;
                }
            }
                
            if (jetPt>ptLeadingJet)
            {
                indexSubLeadingJet = indexLeadingJet;
                ptSubLeadingJet = ptLeadingJet;
                indexLeadingJet = jet;
                ptLeadingJet = jetPt;
            }
        }//end jet selection
        
        if(!isData)
        {
            //start jet selection again with UP JEC
            for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet)
            {
                float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
                float jetEta = analysisTree.pfjet_eta[jet];
                if (absJetEta>jetEtaCut) continue;
                    
                float jetPtUp = analysisTree.pfjet_pt[jet]*(1.0+analysisTree.pfjet_jecUncertainty[jet]);

                //std::cout << jet << " : uncertainty = " << analysisTree.pfjet_jecUncertainty[jet] << std::endl;
                if (jetPtUp<jetPtLowCut) continue;
                    
                bool isPFJetId =false;
                if (era=="2016") isPFJetId= looseJetiD_2016(analysisTree,int(jet));
                else if (era=="2017") isPFJetId = tightJetiD_2017(analysisTree,int(jet));
                else if (era=="2018") isPFJetId = tightJetiD_2018(analysisTree,int(jet));
                else if ((era=="UL2018")||(era=="UL2017")) isPFJetId = tightJetiD_UL2018UL2017(analysisTree,int(jet));
                else if ((era=="UL2016postVFP")||(era=="UL2016preVFP")) isPFJetId = tightJetiD_UL2016(analysisTree,int(jet));
                if (!isPFJetId) continue;
                
                //bool isPFJetId = tightJetiD_2018(analysisTree,int(jet));
                //if (!isPFJetId) continue;
                // cleaned jet at the moment!
                
                bool cleanedJet = true;
                float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                   eta_1,phi_1);
                if (dR1<dRJetLeptonCut) cleanedJet = false;
                float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                   eta_2,phi_2);
                if (dR2<dRJetLeptonCut) cleanedJet = false;
                if(HiggsFinalState == "EM")
                {
                    float dR3 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                       eta_3,phi_3);
                    if (dR3<dRJetLeptonCut) cleanedJet = false;
                    float dR4 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                       eta_4,phi_4);
                    if (dR4<dRJetLeptonCut) cleanedJet = false;
                }
                if((HiggsFinalState == "ET") ||(HiggsFinalState == "MT"))
                {
                    float dR3 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                       eta_3,phi_3);
                    if (dR3<dRJetLeptonCut) cleanedJet = false;
                }
                if (!cleanedJet) continue;
                
                // jetId
                if (jetPtUp>jetPtLowCut)
                    jetspt20Up.push_back(jet);
                
                if (absJetEta<bJetEtaCut)
                {   // jet within b-tagging acceptance
                    bool tagged = false;
                    bool tagged_mistagUp = false;
                    bool tagged_mistagDown = false;
                    bool tagged_btagUp = false;
                    bool tagged_btagDown = false;
                    //std::cout << "D1 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant1] << std::endl;
                    //std::cout << "D2 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant2] << std::endl;
                    //std::cout << "D3 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant3] << std::endl;

                    tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2] + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] > btagCut;
                    bool taggedRaw = tagged;
                    
                    if (!isData)
                    {
                        int flavor = abs(analysisTree.pfjet_flavour[jet]);
                        
                        double jet_scalefactor      = 1;
                        double JetPtForBTag = jetPtUp;
                        double tageff = 1;
                        if (flavor==5)
                        {
                            if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                            if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                            jet_scalefactor      = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                            tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
                        }
                        else if (flavor==4)
                        {
                            if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                            if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                            jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                            tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
                        }
                        else
                        {
                            if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
                            if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
                            jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                            tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
                        }
                        if (tageff<1e-5)      tageff = 1e-5;
                        if (tageff>0.99999)   tageff = 0.99999;
                        r.SetSeed((int)((jetEta+5)*100000));
                        double rannum = r.Rndm();
                        if (tagged)
                        { // demote
                            if(jet_scalefactor<1)
                            {
                                double fraction = 1-jet_scalefactor;
                                if (rannum<fraction) tagged = false;
                            }
                        }
                        else if (!tagged)
                        { // promote
                            if(jet_scalefactor>1)
                            {
                                double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
                                if (rannum<fraction) tagged = true;
                            }
                        }
                    }
                    //missing b-tag SFs
                    if (taggedRaw) bjetsRawUp.push_back(jet);

                    if (tagged)
                    {
                        bjetsUp.push_back(jet);
                        if (jetPtUp<ptLeadingBJetUp&&jetPtUp>ptSubLeadingBJetUp)
                        {
                            indexSubLeadingBJetUp = jet;
                            ptSubLeadingBJetUp = jetPtUp;
                        }
                        if (jetPtUp>ptLeadingBJetUp)
                        {
                            ptLeadingBJetUp = jetPtUp;
                            indexLeadingBJetUp = jet;
                        }
                    }
                }

                if (jetPtUp>jetPtHighCut)
                    jetsUp.push_back(jet);
                    
                if (indexLeadingJetUp>=0)
                {
                    if (jetPtUp<ptLeadingJetUp&&jetPtUp>ptSubLeadingJetUp)
                    {
                        indexSubLeadingJetUp = jet;
                        ptSubLeadingJetUp = jetPtUp;
                    }
                }
                    
                if (jetPtUp>ptLeadingJetUp)
                {
                    indexSubLeadingJetUp = indexLeadingJetUp;
                    ptSubLeadingJetUp = ptLeadingJetUp;
                    indexLeadingJetUp = jet;
                    ptLeadingJetUp = jetPtUp;
                }
            }//end jet selection for UP JEC
            
            //start jet selection again with Down JEC
            for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet)
            {
                float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
                float jetEta = analysisTree.pfjet_eta[jet];
                if (absJetEta>jetEtaCut) continue;
                    
                float jetPtDown = analysisTree.pfjet_pt[jet]*(1.0-analysisTree.pfjet_jecUncertainty[jet]);

                //std::cout << jet << " : uncertainty = " << analysisTree.pfjet_jecUncertainty[jet] << std::endl;
                if (jetPtDown<jetPtLowCut) continue;
                    
                bool isPFJetId =false;
                if (era=="2016") isPFJetId= looseJetiD_2016(analysisTree,int(jet));
                else if (era=="2017") isPFJetId = tightJetiD_2017(analysisTree,int(jet));
                else if (era=="2018") isPFJetId = tightJetiD_2018(analysisTree,int(jet));
                else if ((era=="UL2018")||(era=="UL2017")) isPFJetId = tightJetiD_UL2018UL2017(analysisTree,int(jet));
                else if ((era=="UL2016postVFP")||(era=="UL2016preVFP")) isPFJetId = tightJetiD_UL2016(analysisTree,int(jet));
                if (!isPFJetId) continue;
                
                //bool isPFJetId = tightJetiD_2018(analysisTree,int(jet));
                //if (!isPFJetId) continue;
                // cleaned jet at the moment!
                
                bool cleanedJet = true;
                float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                   eta_1,phi_1);
                if (dR1<dRJetLeptonCut) cleanedJet = false;
                float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                   eta_2,phi_2);
                if (dR2<dRJetLeptonCut) cleanedJet = false;
                if(HiggsFinalState == "EM")
                {
                    float dR3 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                       eta_3,phi_3);
                    if (dR3<dRJetLeptonCut) cleanedJet = false;
                    float dR4 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                       eta_4,phi_4);
                    if (dR4<dRJetLeptonCut) cleanedJet = false;
                }
                if((HiggsFinalState == "ET") ||(HiggsFinalState == "MT"))
                {
                    float dR3 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                       eta_3,phi_3);
                    if (dR3<dRJetLeptonCut) cleanedJet = false;
                }
                if (!cleanedJet) continue;
                
                // jetId
                if (jetPtDown>jetPtLowCut)
                    jetspt20Down.push_back(jet);
                
                if (absJetEta<bJetEtaCut)
                {   // jet within b-tagging acceptance
                    bool tagged = false;
                    //std::cout << "D1 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant1] << std::endl;
                    //std::cout << "D2 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant2] << std::endl;
                    //std::cout << "D3 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant3] << std::endl;

                    tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2] + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] > btagCut;
                    bool taggedRaw = tagged;
                    
                    if (!isData)
                    {
                        int flavor = abs(analysisTree.pfjet_flavour[jet]);
                        
                        double jet_scalefactor      = 1;
                        double JetPtForBTag = jetPtDown;
                        double tageff = 1;
                        if (flavor==5)
                        {
                            if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                            if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                            jet_scalefactor      = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                            tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
                        }
                        else if (flavor==4)
                        {
                            if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                            if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                            jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                            tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
                        }
                        else
                        {
                            if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
                            if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
                            jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                            tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
                        }
                        if (tageff<1e-5)      tageff = 1e-5;
                        if (tageff>0.99999)   tageff = 0.99999;
                        r.SetSeed((int)((jetEta+5)*100000));
                        double rannum = r.Rndm();
                        if (tagged)
                        { // demote
                            if(jet_scalefactor<1)
                            {
                                double fraction = 1-jet_scalefactor;
                                if (rannum<fraction) tagged = false;
                            }
                        }
                        else if (!tagged)
                        { // promote
                            if(jet_scalefactor>1)
                            {
                                double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
                                if (rannum<fraction) tagged = true;
                            }
                        }
                    }
                    //missing b-tag SFs
                    if (taggedRaw) bjetsRawDown.push_back(jet);

                    if (tagged)
                    {
                        bjetsDown.push_back(jet);
                        if (jetPtDown<ptLeadingBJetDown&&jetPtDown>ptSubLeadingBJetDown)
                        {
                            indexSubLeadingBJetDown = jet;
                            ptSubLeadingBJetDown = jetPtDown;
                        }
                        if (jetPtDown>ptLeadingBJetDown)
                        {
                            ptLeadingBJetDown = jetPtDown;
                            indexLeadingBJetDown = jet;
                        }
                    }
                }

                if (jetPtDown>jetPtHighCut)
                    jetsDown.push_back(jet);
                    
                if (indexLeadingJetDown>=0)
                {
                    if (jetPtDown<ptLeadingJetDown&&jetPtDown>ptSubLeadingJetDown)
                    {
                        indexSubLeadingJetDown = jet;
                        ptSubLeadingJetDown = jetPtDown;
                    }
                }
                    
                if (jetPtDown>ptLeadingJetDown)
                {
                    indexSubLeadingJetDown = indexLeadingJetDown;
                    ptSubLeadingJetDown = ptLeadingJetDown;
                    indexLeadingJetDown = jet;
                    ptLeadingJetDown = jetPtDown;
                }
            }//end jet selection for Down JEC
            
            //Start JEC HEM 2018 issue
            if(era == "UL2018")
            {
                //Up HEM 2018 JEC
                for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet)
                {
                    float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
                    float jetEta = analysisTree.pfjet_eta[jet];
                    if (absJetEta>jetEtaCut) continue;

                    bool isPFJetId =false;
                    isPFJetId = tightJetiD_UL2018UL2017(analysisTree,int(jet));
                    if (!isPFJetId) continue;
                    if (analysisTree.pfjet_pt[jet] < 15) continue;
                    
                    float jetPtUp = analysisTree.pfjet_pt[jet];

                    if (analysisTree.pfjet_phi[jet] > (-1.57) && analysisTree.pfjet_phi[jet] < (-0.87) && jetEta > (-2.5) && jetEta < (-1.3))
                    {
                        jetPtUp = analysisTree.pfjet_pt[jet]*(1.2);
                    }

                    if (analysisTree.pfjet_phi[jet] > (-1.57) && analysisTree.pfjet_phi[jet] < (-0.87) && jetEta > (-3.0) && jetEta < (-2.5))
                    {
                        jetPtUp = analysisTree.pfjet_pt[jet]*(1.35);
                    }
                    
                    //std::cout << jet << " : uncertainty = " << analysisTree.pfjet_jecUncertainty[jet] << std::endl;
                    if (jetPtUp<jetPtLowCut) continue;
                    
                    bool cleanedJet = true;
                    float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                       eta_1,phi_1);
                    if (dR1<dRJetLeptonCut) cleanedJet = false;
                    float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                       eta_2,phi_2);
                    if (dR2<dRJetLeptonCut) cleanedJet = false;
                    if(HiggsFinalState == "EM")
                    {
                        float dR3 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                           eta_3,phi_3);
                        if (dR3<dRJetLeptonCut) cleanedJet = false;
                        float dR4 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                           eta_4,phi_4);
                        if (dR4<dRJetLeptonCut) cleanedJet = false;
                    }
                    if((HiggsFinalState == "ET") ||(HiggsFinalState == "MT"))
                    {
                        float dR3 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                           eta_3,phi_3);
                        if (dR3<dRJetLeptonCut) cleanedJet = false;
                    }
                    if (!cleanedJet) continue;
                    
                    // jetId
                    if (jetPtUp>jetPtLowCut)
                        jetspt20HEM2018Up.push_back(jet);
                    
                    if (absJetEta<bJetEtaCut)
                    {   // jet within b-tagging acceptance
                        bool tagged = false;
                        //std::cout << "D1 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant1] << std::endl;
                        //std::cout << "D2 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant2] << std::endl;
                        //std::cout << "D3 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant3] << std::endl;

                        tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2] + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] > btagCut;
                        bool taggedRaw = tagged;
                        
                        if (!isData)
                        {
                            int flavor = abs(analysisTree.pfjet_flavour[jet]);
                            
                            double jet_scalefactor      = 1;
                            double JetPtForBTag = jetPtUp;
                            double tageff = 1;
                            if (flavor==5)
                            {
                                if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                                if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                                jet_scalefactor      = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                                tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
                            }
                            else if (flavor==4)
                            {
                                if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                                if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                                jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                                tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
                            }
                            else
                            {
                                if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
                                if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
                                jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                                tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
                            }
                            if (tageff<1e-5)      tageff = 1e-5;
                            if (tageff>0.99999)   tageff = 0.99999;
                            r.SetSeed((int)((jetEta+5)*100000));
                            double rannum = r.Rndm();
                            if (tagged)
                            { // demote
                                if(jet_scalefactor<1)
                                {
                                    double fraction = 1-jet_scalefactor;
                                    if (rannum<fraction) tagged = false;
                                }
                            }
                            else if (!tagged)
                            { // promote
                                if(jet_scalefactor>1)
                                {
                                    double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
                                    if (rannum<fraction) tagged = true;
                                }
                            }
                        }
                        if (tagged)
                        {
                            bjetsHEM2018Up.push_back(jet);
                        }
                    }
                    if (jetPtUp>jetPtHighCut)
                        jetsHEM2018Up.push_back(jet);
                }//end jet selection for HEM 2018 Up JEC
                
                //Down HEM 2018 JEC
                for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet)
                {
                    float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
                    float jetEta = analysisTree.pfjet_eta[jet];
                    if (absJetEta>jetEtaCut) continue;

                    bool isPFJetId =false;
                    isPFJetId = tightJetiD_UL2018UL2017(analysisTree,int(jet));
                    if (!isPFJetId) continue;
                    if (analysisTree.pfjet_pt[jet] < 15) continue;
                    float jetPtDown = analysisTree.pfjet_pt[jet];

                    if (analysisTree.pfjet_phi[jet] > (-1.57) && analysisTree.pfjet_phi[jet] < (-0.87) && jetEta > (-2.5) && jetEta < (-1.3))
                    {
                        jetPtDown = analysisTree.pfjet_pt[jet]*(0.8);
                    }
                    if (analysisTree.pfjet_phi[jet] > (-1.57) && analysisTree.pfjet_phi[jet] < (-0.87) && jetEta > (-3.0) && jetEta < (-2.5))
                    {
                        jetPtDown = analysisTree.pfjet_pt[jet]*(0.65);
                    }
                    //std::cout << jet << " : uncertainty = " << analysisTree.pfjet_jecUncertainty[jet] << std::endl;
                    if (jetPtDown<jetPtLowCut) continue;
                        
                    bool cleanedJet = true;
                    float dR1 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                       eta_1,phi_1);
                    if (dR1<dRJetLeptonCut) cleanedJet = false;
                    float dR2 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                       eta_2,phi_2);
                    if (dR2<dRJetLeptonCut) cleanedJet = false;
                    if(HiggsFinalState == "EM")
                    {
                        float dR3 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                           eta_3,phi_3);
                        if (dR3<dRJetLeptonCut) cleanedJet = false;
                        float dR4 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                           eta_4,phi_4);
                        if (dR4<dRJetLeptonCut) cleanedJet = false;
                    }
                    if((HiggsFinalState == "ET") ||(HiggsFinalState == "MT"))
                    {
                        float dR3 = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
                                           eta_3,phi_3);
                        if (dR3<dRJetLeptonCut) cleanedJet = false;
                    }
                    if (!cleanedJet) continue;
                    
                    // jetId
                    if (jetPtDown>jetPtLowCut)
                        jetspt20HEM2018Down.push_back(jet);
                    
                    if (absJetEta<bJetEtaCut)
                    {   // jet within b-tagging acceptance
                        bool tagged = false;
                        //std::cout << "D1 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant1] << std::endl;
                        //std::cout << "D2 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant2] << std::endl;
                        //std::cout << "D3 = " << analysisTree.pfjet_btag[jet][nBTagDiscriminant3] << std::endl;

                        tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant1] + analysisTree.pfjet_btag[jet][nBTagDiscriminant2] + analysisTree.pfjet_btag[jet][nBTagDiscriminant3] > btagCut;
                        bool taggedRaw = tagged;
                        
                        if (!isData)
                        {
                            int flavor = abs(analysisTree.pfjet_flavour[jet]);
                            
                            double jet_scalefactor      = 1;
                            double JetPtForBTag = jetPtDown;
                            double tageff = 1;
                            if (flavor==5)
                            {
                                if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                                if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                                jet_scalefactor      = reader_BTAG.eval_auto_bounds("central",BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
                                tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
                            }
                            else if (flavor==4)
                            {
                                if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
                                if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
                                jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
                                tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
                            }
                            else
                            {
                                if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
                                if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
                                jet_scalefactor      = reader_BTAG.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
                                tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
                            }
                            if (tageff<1e-5)      tageff = 1e-5;
                            if (tageff>0.99999)   tageff = 0.99999;
                            r.SetSeed((int)((jetEta+5)*100000));
                            double rannum = r.Rndm();
                            if (tagged)
                            { // demote
                                if(jet_scalefactor<1)
                                {
                                    double fraction = 1-jet_scalefactor;
                                    if (rannum<fraction) tagged = false;
                                }
                            }
                            else if (!tagged)
                            { // promote
                                if(jet_scalefactor>1)
                                {
                                    double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
                                    if (rannum<fraction) tagged = true;
                                }
                            }
                        }
                        if (tagged)
                        {
                            bjetsHEM2018Down.push_back(jet);
                        }
                    }

                    if (jetPtDown>jetPtHighCut)
                        jetsHEM2018Down.push_back(jet);
                }//end jet selection for HEM 2018 Down JEC
            }
        }
        
        njets = jets.size();
        njetsUp = jetsUp.size();
        njetsDown = jetsDown.size();
        njetsHEM2018Up = jetsHEM2018Up.size();
        njetsHEM2018Down = jetsHEM2018Down.size();

        int njetsMax = njets;

        njetspt20 = jetspt20.size();
        njetspt20Up = jetspt20Up.size();
        njetspt20Down = jetspt20Down.size();
        njetspt20HEM2018Up = jetspt20HEM2018Up.size();
        njetspt20HEM2018Down = jetspt20HEM2018Down.size();

        nbtag = bjets.size();
        nbtagUp = bjetsUp.size();
        nbtagDown = bjetsDown.size();
        nbtagHEM2018Up = bjetsHEM2018Up.size();
        nbtagHEM2018Down = bjetsHEM2018Down.size();

        nbtag_mistagUp   = bjets_mistagUp.size();
        nbtag_mistagDown = bjets_mistagDown.size();
        nbtag_btagUp   = bjets_btagUp.size();
        nbtag_btagDown = bjets_btagDown.size();
        nbtag_noSF = bjetsRaw.size();
        
        bpt_1 = -9999;
        beta_1 = -9999;
        bphi_1 = -9999;
        bhadronFlavor_1 = -1;
        bpt_1_Up = -9999;
        beta_1_Up = -9999;
        bphi_1_Up = -9999;
        bpt_1_Down = -9999;
        beta_1_Down = -9999;
        bphi_1_Down = -9999;
        
        bpt_2 = -9999;
        beta_2 = -9999;
        bphi_2 = -9999;
        bhadronFlavor_2 = -1;
        bpt_2_Up = -9999;
        beta_2_Up = -9999;
        bphi_2_Up = -9999;
        bpt_2_Down = -9999;
        beta_2_Down = -9999;
        bphi_2_Down = -9999;

        if (indexLeadingBJet>=0)
        {
            bpt_1 = analysisTree.pfjet_pt[indexLeadingBJet];
            beta_1 = analysisTree.pfjet_eta[indexLeadingBJet];
            bphi_1 = analysisTree.pfjet_phi[indexLeadingBJet];
        }

        if(!isData)
        {
            bhadronFlavor_1 = abs(analysisTree.pfjet_flavour[indexLeadingBJet]);
            if (indexLeadingBJetUp>=0)
            {
                bpt_1_Up = analysisTree.pfjet_pt[indexLeadingBJetUp]*(1.0+analysisTree.pfjet_jecUncertainty[indexLeadingBJetUp]);
                beta_1_Up = analysisTree.pfjet_eta[indexLeadingBJetUp];
                bphi_1_Up = analysisTree.pfjet_phi[indexLeadingBJetUp];
            }
            if (indexLeadingBJetDown>=0)
            {
                bpt_1_Down = analysisTree.pfjet_pt[indexLeadingBJetDown]*(1.0-analysisTree.pfjet_jecUncertainty[indexLeadingBJetDown]);
                beta_1_Down = analysisTree.pfjet_eta[indexLeadingBJetDown];
                bphi_1_Down = analysisTree.pfjet_phi[indexLeadingBJetDown];
            }
        }
        if (indexSubLeadingBJet>=0)
        {
            bpt_2 = analysisTree.pfjet_pt[indexSubLeadingBJet];
            beta_2 = analysisTree.pfjet_eta[indexSubLeadingBJet];
            bphi_2 = analysisTree.pfjet_phi[indexSubLeadingBJet];
        }
        
        if(!isData)
        {
            bhadronFlavor_2 = abs(analysisTree.pfjet_flavour[indexSubLeadingJet]);
            if (indexSubLeadingBJetUp>=0)
            {
                bpt_2_Up = analysisTree.pfjet_pt[indexSubLeadingBJetUp]*(1.0+analysisTree.pfjet_jecUncertainty[indexSubLeadingBJetUp]);
                beta_2_Up = analysisTree.pfjet_eta[indexSubLeadingBJetUp];
                bphi_2_Up = analysisTree.pfjet_phi[indexSubLeadingBJetUp];
            }
            if (indexSubLeadingBJetDown>=0)
            {
                bpt_2_Down = analysisTree.pfjet_pt[indexSubLeadingBJetDown]*(1.0-analysisTree.pfjet_jecUncertainty[indexSubLeadingBJetDown]);
                beta_2_Down = analysisTree.pfjet_eta[indexSubLeadingBJetDown];
                bphi_2_Down = analysisTree.pfjet_phi[indexSubLeadingBJetDown];
            }
        }
        jpt_1 = -9999;
        jeta_1 = -9999;
        jphi_1 = -9999;
        jpt_1_Up = -9999;
        jeta_1_Up = -9999;
        jphi_1_Up = -9999;
        jpt_1_Down = -9999;
        jeta_1_Down = -9999;
        jphi_1_Down = -9999;
        if (indexLeadingJet>=0&&indexSubLeadingJet>=0&&indexLeadingJet==indexSubLeadingJet)
            cout << "warning : indexLeadingJet ==indexSubLeadingJet = " << indexSubLeadingJet << endl;
            
        if (indexLeadingJet>=0)
        {
            jpt_1 = analysisTree.pfjet_pt[indexLeadingJet];
            jeta_1 = analysisTree.pfjet_eta[indexLeadingJet];
            jphi_1 = analysisTree.pfjet_phi[indexLeadingJet];
        //        cout << "Leading jet pt = " << jpt_1 << "   eta = " << jeta_1 << endl;
        //        for (auto const& Name : uncertNames) {
        //          cout << "    " << Name << " : " << jecUncertainties->getUncertainty(Name,jpt_1,jeta_1) << endl;
        //        }
        }
        if(!isData)
        {
            if (indexLeadingJetUp>=0)
            {
                jpt_1_Up =analysisTree.pfjet_pt[indexLeadingJetUp]*(1.0+analysisTree.pfjet_jecUncertainty[indexLeadingJetUp]);
                jeta_1_Up = analysisTree.pfjet_eta[indexLeadingJetUp];
                jphi_1_Up = analysisTree.pfjet_phi[indexLeadingJetUp];
            }
            if (indexLeadingJetDown>=0)
            {
                jpt_1_Down =analysisTree.pfjet_pt[indexLeadingJetDown]*(1.0-analysisTree.pfjet_jecUncertainty[indexLeadingJetDown]);
                jeta_1_Down = analysisTree.pfjet_eta[indexLeadingJetDown];
                jphi_1_Down = analysisTree.pfjet_phi[indexLeadingJetDown];
            }
        }

        jpt_2 = -9999;
        jeta_2 = -9999;
        jphi_2 = -9999;
        jpt_2_Up = -9999;
        jeta_2_Up = -9999;
        jphi_2_Up = -9999;
        jpt_2_Down = -9999;
        jeta_2_Down = -9999;
        jphi_2_Down = -9999;
        
        if (indexSubLeadingJet>=0)
        {
            jpt_2 = analysisTree.pfjet_pt[indexSubLeadingJet];
            jeta_2 = analysisTree.pfjet_eta[indexSubLeadingJet];
            jphi_2 = analysisTree.pfjet_phi[indexSubLeadingJet];
        }
        
        if(!isData)
        {
            if (indexSubLeadingJetUp>=0)
            {
                jpt_2_Up =analysisTree.pfjet_pt[indexSubLeadingJetUp]*(1.0+analysisTree.pfjet_jecUncertainty[indexSubLeadingJetUp]);
                jeta_2_Up = analysisTree.pfjet_eta[indexSubLeadingJetUp];
                jphi_2_Up = analysisTree.pfjet_phi[indexSubLeadingJetUp];
            }
            if (indexSubLeadingJetDown>=0)
            {
                jpt_2_Down =analysisTree.pfjet_pt[indexSubLeadingJetDown]*(1.0-analysisTree.pfjet_jecUncertainty[indexSubLeadingJetDown]);
                jeta_2_Down = analysisTree.pfjet_eta[indexSubLeadingJetDown];
                jphi_2_Down = analysisTree.pfjet_phi[indexSubLeadingJetDown];
            }
        }
        
        TString SynCategory("");
        
        if(HiggsFinalState == "ET")
            SynCategory = "mmet";
        if(HiggsFinalState == "MT")
            SynCategory = "mmmt";
        if(HiggsFinalState == "EM")
            SynCategory = "mmem";
        if(HiggsFinalState == "TT")
            SynCategory = "mmtt";
        
        if(SynchExer==true)
        {
            myfile << run << "," << lumi << "," << evt << "," << SynCategory << std::endl;
        }
        
        SynTree->Fill();
        selEvents++;
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
    myfile.close();
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events                     = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree                   = " << nEvents << std::endl;
  std::cout << "Total number of selected events                  = " << selEvents << std::endl;
  std::cout << std::endl;
  std::cout << "RunMin = " << RunMin << std::endl;
  std::cout << "RunMax = " << RunMax << std::endl;

  //cout << "weight used:" << weight << std::endl;

  // using object as comp
  std::sort (allRuns.begin(), allRuns.end(), myobject);
  std::cout << "Runs   :";
  for (unsigned int iR=0; iR<allRuns.size(); ++iR)
    std::cout << " " << allRuns.at(iR);
  std::cout << std::endl;
  
  file->Write();
  file->Close();
  delete file;
  
  }

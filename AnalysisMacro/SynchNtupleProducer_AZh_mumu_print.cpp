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
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TMath.h"

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
    const float btagCut = cfg.get<float>("btagCut");
    const string bTagDiscriminator = cfg.get<string>("BTagAlgorithm");
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
    
    // vertex cuts
    const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");
    const float zVertexCut     = cfg.get<float>("ZVertexCut");
    const float dVertexCut     = cfg.get<float>("DVertexCut");
    
    //scale factor
    const string MuonIdIsoFile = cfg.get<string>("MuonIdIsoFile");
    const string MuonTriggerFile = cfg.get<string>("MuonTriggerFile");
    
    //svfit
    const string svFitPtResFile = cfg.get<string>("svFitPtResFile").data();

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
    Int_t   njets;
    Int_t   njetspt20;
    Float_t jpt_1;
    Float_t jeta_1;
    Float_t jphi_1;
    Float_t jpt_2;
    Float_t jeta_2;
    Float_t jphi_2;
    Float_t bpt_1;
    Float_t beta_1;
    Float_t bphi_1;
    Float_t bhadronFlavor_1;
    Float_t bpt_2;
    Float_t beta_2;
    Float_t bphi_2;
    Float_t bhadronFlavor_2;
    
    //Leg 1 (highest pt elec/mu from Z candidate)
    Float_t pt_1;
    Float_t phi_1;
    Float_t eta_1;
    Float_t iso_1;
    Float_t iso_ea_1;
    Int_t gen_match_1;
    Float_t IdRawMva_1;
    Float_t pfmt_1;
    Float_t muonSF_1;
    Float_t muonIDSF_1;
    Float_t muonTrigSF_1;
    Float_t electronSF_1;
    Float_t electronIDSF_1;
    Bool_t muonLoose_1;
    Bool_t muonMedium_1;
    Bool_t muonTight_1;

    //Leg 2 (subleading pt elec/mu from Z candidate)
    Float_t pt_2;
    Float_t phi_2;
    Float_t eta_2;
    Float_t iso_2;
    Float_t iso_ea_2;
    Int_t gen_match_2;
    Float_t IdRawMva_2;
    Float_t pfmt_2;
    Float_t muonSF_2;
    Float_t muonIDSF_2;
    Float_t muonTrigSF_2;
    Float_t electronSF_2;
    Float_t electronIDSF_2;
    Bool_t muonLoose_2;
    Bool_t muonMedium_2;
    Bool_t muonTight_2;

    //Leg 3 (Based off of Higgs legs: ET -> E, MT -> M, EM -> E, TT -> highest pt tau)
    Float_t pt_3;
    Float_t phi_3;
    Float_t eta_3;
    Float_t iso_3; //if tau, then it is VSjet raw value
    Float_t iso_ea_3;
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
    Float_t muonIDSF_3;
    Float_t electronSF_3;
    Float_t electronIDSF_3;
    Float_t tauSF_3;
    Bool_t electronIDWP90v2_3;
    Bool_t electronIDWP80v2_3;
    Bool_t muonLoose_3;
    Bool_t muonMedium_3;
    Bool_t muonTight_3;
    
    //Leg 4 (Based off of Higgs legs: ET -> T, MT -> T, EM -> M, TT -> subleading pt tau)
    Float_t pt_4;
    Float_t phi_4;
    Float_t eta_4;
    Float_t iso_4; //if tau, then it is VSjet raw value
    Float_t iso_ea_4;
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
    Float_t muonIDSF_4;
    Float_t electronSF_4;
    Float_t electronIDSF_4;
    Float_t tauSF_4;
    Bool_t electronIDWP90v2_4;
    Bool_t electronIDWP80v2_4;
    Bool_t muonLoose_4;
    Bool_t muonMedium_4;
    Bool_t muonTight_4;
    
    Float_t met;
    Float_t metphi;
    Float_t metcov00;
    Float_t metcov01;
    Float_t metcov10;
    Float_t metcov11;
    
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
    
    Float_t Mass;
    Float_t Mass_sv;
    Float_t Mass_svc;
    Float_t Mass_fast;

    UInt_t  npartons;

    SynTree->Branch("run",&run,"run/l");
    SynTree->Branch("lumi",&lumi,"lumi/i");
    SynTree->Branch("evt",&evt,"evt/l");
    SynTree->Branch("npv",&npv,"npv/I");
    SynTree->Branch("npu",&npu,"npu/I");
    SynTree->Branch("rho",&rho,"rho/F");
    SynTree->Branch("mcweight",&mcweight,"mcweight/F");
    SynTree->Branch("puweight", &puweight, "puweight/F");
    SynTree->Branch("nbtag",&nbtag,"nbtag/I");
    SynTree->Branch("njets",&njets,"njets/I");
    SynTree->Branch("njetspt20",&njetspt20,"njetspt20/I");
    SynTree->Branch("jpt_1",&jpt_1,"jpt_1/F");
    SynTree->Branch("jeta_1",&jeta_1,"jeta_1/F");
    SynTree->Branch("jphi_1",&jphi_1,"jphi_1/F");
    SynTree->Branch("jpt_2",&jpt_2,"jpt_2/F");
    SynTree->Branch("jeta_2",&jeta_2,"jeta_2/F");
    SynTree->Branch("jphi_2",&jphi_2,"jphi_2/F");
    SynTree->Branch("bpt_1",&bpt_1,"bpt_1/F");
    SynTree->Branch("beta_1",&beta_1,"beta_1/F");
    SynTree->Branch("bphi_1",&bphi_1,"bphi_1/F");
    SynTree->Branch("bhadronFlavor_1",&bhadronFlavor_1,"bhadronFlavor_1/F");
    SynTree->Branch("bpt_2",&bpt_2,"bpt_2/F");
    SynTree->Branch("beta_2",&beta_2,"beta_2/F");
    SynTree->Branch("bphi_2",&bphi_2,"bphi_2/F");
    SynTree->Branch("bhadronFlavor_2",&bhadronFlavor_2,"bhadronFlavor_2/F");
    
    SynTree->Branch("pt_1",&pt_1,"pt_1/F");
    SynTree->Branch("phi_1",&phi_1,"phi_1/F");
    SynTree->Branch("eta_1",&eta_1,"eta_1/F");
    SynTree->Branch("iso_1",&iso_1,"iso_1/F");
    SynTree->Branch("iso_ea_1",&iso_ea_1,"iso_ea_1/F");
    SynTree->Branch("gen_match_1",&gen_match_1,"gen_match_1/I");
    SynTree->Branch("IdRawMva_1",&IdRawMva_1,"IdRawMva_1/F");
    SynTree->Branch("pfmt_1",&pfmt_1,"pfmt_1/F");
    SynTree->Branch("muonSF_1",&muonSF_1,"muonSF_1/F");
    SynTree->Branch("muonIDSF_1",&muonIDSF_1,"muonIDSF_1/F");
    SynTree->Branch("muonTrigSF_1",&muonTrigSF_1,"muonTrigSF_1/F");
    SynTree->Branch("electronSF_1",&electronSF_1,"electronSF_1/F");
    SynTree->Branch("electronIDSF_1",&electronIDSF_1,"electronIDSF_1/F");
    SynTree->Branch("muonLoose_1",&muonLoose_1,"muonLoose_1/O");
    SynTree->Branch("muonMedium_1",&muonMedium_1,"muonMedium_1/O");
    SynTree->Branch("muonTight_1",&muonTight_1,"muonTight_1/O");

    SynTree->Branch("pt_2",&pt_2,"pt_2/F");
    SynTree->Branch("phi_2",&phi_2,"phi_2/F");
    SynTree->Branch("eta_2",&eta_2,"eta_2/F");
    SynTree->Branch("iso_2",&iso_2,"iso_2/F");
    SynTree->Branch("iso_ea_2",&iso_ea_2,"iso_ea_2/F");
    SynTree->Branch("gen_match_2",&gen_match_2,"gen_match_2/I");
    SynTree->Branch("IdRawMva_2",&IdRawMva_2,"IdRawMva_2/F");
    SynTree->Branch("pfmt_2",&pfmt_2,"pfmt_2/F");
    SynTree->Branch("muonSF_2",&muonSF_2,"muonSF_2/F");
    SynTree->Branch("muonIDSF_2",&muonIDSF_2,"muonIDSF_2/F");
    SynTree->Branch("muonTrigSF_2",&muonTrigSF_2,"muonTrigSF_2/F");
    SynTree->Branch("electronSF_2",&electronSF_2,"electronSF_2/F");
    SynTree->Branch("electronIDSF_2",&electronIDSF_2,"electronIDSF_2/F");
    SynTree->Branch("muonLoose_2",&muonLoose_2,"muonLoose_2/O");
    SynTree->Branch("muonMedium_2",&muonMedium_2,"muonMedium_2/O");
    SynTree->Branch("muonTight_2",&muonTight_2,"muonTight_2/O");

    SynTree->Branch("pt_3",&pt_3,"pt_3/F");
    SynTree->Branch("phi_3",&phi_3,"phi_3/F");
    SynTree->Branch("eta_3",&eta_3,"eta_3/F");
    SynTree->Branch("iso_3",&iso_3,"iso_3/F");
    SynTree->Branch("iso_ea_3",&iso_ea_3,"iso_ea_3/F");
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
    SynTree->Branch("muonIDSF_3",&muonIDSF_3,"muonIDSF_3/F");
    SynTree->Branch("electronSF_3",&electronSF_3,"electronSF_3/F");
    SynTree->Branch("electronIDSF_3",&electronIDSF_3,"electronIDSF_3/F");
    SynTree->Branch("tauSF_3",&tauSF_3,"tauSF_3/F");
    SynTree->Branch("muonLoose_3",&muonLoose_3,"muonLoose_3/O");
    SynTree->Branch("muonMedium_3",&muonMedium_3,"muonMedium_3/O");
    SynTree->Branch("muonTight_3",&muonTight_3,"muonTight_3/O");
    SynTree->Branch("electronSF_3",&electronSF_3,"electronSF_3/F");
    SynTree->Branch("electronIDSF_3",&electronIDSF_3,"electronIDSF_3/F");
    SynTree->Branch("electronIDWP90v2_3",&electronIDWP90v2_3,"electronIDWP90v2_3/O");
    SynTree->Branch("electronIDWP80v2_3",&electronIDWP80v2_3,"electronIDWP80v2_3/O");
    
    SynTree->Branch("pt_4",&pt_4,"pt_4/F");
    SynTree->Branch("phi_4",&phi_4,"phi_4/F");
    SynTree->Branch("eta_4",&eta_4,"eta_4/F");
    SynTree->Branch("iso_4",&iso_4,"iso_4/F");
    SynTree->Branch("iso_ea_4",&iso_ea_3,"iso_ea_4/F");
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
    SynTree->Branch("muonIDSF_4",&muonIDSF_4,"muonIDSF_4/F");
    SynTree->Branch("electronSF_4",&electronSF_4,"electronSF_4/F");
    SynTree->Branch("electronIDSF_4",&electronIDSF_4,"electronIDSF_4/F");
    SynTree->Branch("tauSF_4",&tauSF_4,"tauSF_4/F");
    SynTree->Branch("muonLoose_4",&muonLoose_4,"muonLoose_4/O");
    SynTree->Branch("muonMedium_4",&muonMedium_4,"muonMedium_4/O");
    SynTree->Branch("muonTight_4",&muonTight_4,"muonTight_4/O");
    SynTree->Branch("electronSF_4",&electronSF_4,"electronSF_4/F");
    SynTree->Branch("electronIDSF_4",&electronIDSF_4,"electronIDSF_4/F");
    SynTree->Branch("electronIDWP90v2_4",&electronIDWP90v2_4,"electronIDWP90v2_4/O");
    SynTree->Branch("electronIDWP80v2_4",&electronIDWP80v2_4,"electronIDWP80v2_4/O");
    
    SynTree->Branch("met",&met,"met/F");
    SynTree->Branch("metphi",&metphi,"metphi/F");
    SynTree->Branch("metcov00",&metcov00,"metcov00/F");
    SynTree->Branch("metcov01",&metcov01,"metcov01/F");
    SynTree->Branch("metcov10",&metcov10,"metcov10/F");
    SynTree->Branch("metcov11",&metcov11,"metcov11/F");
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

    SynTree->Branch("Mass",&Mass,"Mass/F");
    SynTree->Branch("Mass_sv",&Mass_sv,"Mass_sv/F");
    SynTree->Branch("Mass_svc",&Mass_svc,"Mass_svc/F");
    SynTree->Branch("Mass_fast",&Mass_fast,"Mass_fast/F");
    SynTree->Branch("npartons",&npartons,"npartons/i");

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
    // Muon scale factors
    ScaleFactor * SF_muonIdIso = new ScaleFactor();
    SF_muonIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonIdIsoFile));
    ScaleFactor * SF_muonTrig = new ScaleFactor();
    SF_muonTrig->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(MuonTriggerFile));
    
    ofstream myfile;
    if(HiggsFinalState == "MT")
        myfile.open ("DESYPrint_mmmt.txt");
    if(HiggsFinalState == "ET")
        myfile.open ("DESYPrint_mmet.txt");
    if(HiggsFinalState == "EM")
        myfile.open ("DESYPrint_mmem.txt");
    if(HiggsFinalState == "TT")
        myfile.open ("DESYPrint_mmtt.txt");
    
    
    // Zpt reweighting for LO DY samples
    TFile * f_zptweight = new TFile(TString(cmsswBase)+"/src/"+ZptweightFile,"read");
    TH2D * h_zptweight = (TH2D*)f_zptweight->Get("zptmass_histo");

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  int nTotalFiles = 0;
    
    //SV fit
  TH1::AddDirectory(false);
  TFile *inputFile_visPtResolution = new TFile(TString(cmsswBase)+"/src/"+svFitPtResFile.data());

  std::cout << "inputFile_visPtResolution : " << std::endl;

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

        float weight = 1;
        puweight = 1;
        mcweight = 1;
        muonTrigSF_1 = 1;
        muonTrigSF_2 = 1;
        mcweight = analysisTree.genweight;

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
        
        //pile up weight variable
        if (!isData) {
            puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
        }
        
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
      
        if(!((HiggsFinalState== "MT")||(HiggsFinalState== "ET")||(HiggsFinalState== "EM")||(HiggsFinalState== "EE")||(HiggsFinalState== "MM")||(HiggsFinalState== "TT")))
        {
            std::cout << "you are not considering one of the final states of Higgs decay"<< std::endl;
            exit(-1);
        }
        
        //select events to be printed
        //if(!((evt==)||(evt==)||(evt==)||(evt==)||(evt==)||(evt==)||(evt==)))//PU-only template
        //if(!((evt==2435542)||(evt==2450214)||(evt==2496060)||(evt==2478721)||(evt==1631840)||(evt==1700143)||(evt==1709409)))//PU-only mmem 2022.05.23
        //if(!((evt==1549901)||(evt==2741245)||(evt==2464017)||(evt==2527463)||(evt==2526128)||(evt==2470164)||(evt==2558360)))//PU-only mmet 2022.05.23
        if(!((evt==2473159)||(evt==2293503)||(evt==2503104)||(evt==2565017)||(evt==2559449)||(evt==2482196)||(evt==2439767)))//PU-only mmmt 2022.05.23
        //if(!((evt==1629981)||(evt==1806621)||(evt==194847)||(evt==162484)||(evt==1816150)||(evt==1988055)||(evt==1721662)))//PU-only mmtt 2022.05.23
        //if(!((evt==2457597)))//DESY-only mmet 2022.05.23
        //if(!((evt==2174715)))//DESY-only mmmt 2022.05.24
        //if(!((evt==2532089)))//DESY-only mmtt 2022.05.24
        {
            continue;
        }
        
        // vertex cuts
        float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+analysisTree.primvertex_y*analysisTree.primvertex_y);

        myfile << "---- Event Number: " << evt <<std::endl;
        myfile << "-------------------------" <<std::endl;
        myfile << "--PV (Primary Vertex) " <<std::endl;
        myfile << "--PV z   : " << fabs(analysisTree.primvertex_z) << std::endl;
        myfile << "--PV ndof: " << analysisTree.primvertex_ndof << std::endl;
        myfile << "--PV dist: " << dVertex << std::endl;
        if((fabs(analysisTree.primvertex_z)<zVertexCut)&&(analysisTree.primvertex_ndof>ndofVertexCut)&&(dVertex<dVertexCut))
        {
            myfile << "--PV is good! " << std::endl;
        }
        else
        {
            myfile << "--PV is bad! " << std::endl;
        }
        //initiation of btag
        unsigned int nBTagDiscriminant = 0;
        for (unsigned int iBTag=0; iBTag < analysisTree.run_btagdiscriminators->size(); ++iBTag)
        {
            TString discr(analysisTree.run_btagdiscriminators->at(iBTag));
            if (discr.Contains(bTagDiscriminator))
                nBTagDiscriminant = iBTag;
        }
        
        // electron selection
        vector<unsigned int> idEles; idEles.clear(); //identified electrons
        vector<unsigned int> allEles; allEles.clear(); //all electrons
        for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
            allEles.push_back(ie);
            if (analysisTree.electron_pt[ie]<ptEleCut) continue;
            if (fabs(analysisTree.electron_eta[ie])>=etaEleCut) continue; //the "=" is needed because scale factor tools may return error
            if (fabs(analysisTree.electron_dxy[ie])>dxyEleCut) continue;
            if (fabs(analysisTree.electron_dz[ie])>dzEleCut) continue;
            //if (!analysisTree.electron_mva_wp90_noIso_Fall17_v2[ie]) continue;
            if (!analysisTree.electron_pass_conversion[ie]) continue;
            if (analysisTree.electron_nmissinghits[ie]>=2) continue;
            idEles.push_back(ie);
            //    cout << "pt:" << analysisTree.electron_pt[ie] << "  passed:" << elePassed << endl;
        }
        
        // muon selection
        vector<unsigned int> idMuons; idMuons.clear(); //identified muons
        vector<unsigned int> allMuons; allMuons.clear(); //all muons
        for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
            allMuons.push_back(im);
            if (analysisTree.muon_pt[im]<ptMuonCut) continue;
            if (fabs(analysisTree.muon_eta[im])>=etaMuonCut) continue; //the "=" is needed because scale factor tools may return error
            if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
            if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
            //if (!(analysisTree.muon_isLoose[im])) continue;
            if (!(analysisTree.muon_isGlobal[im]||analysisTree.muon_isTracker[im])) continue;
            idMuons.push_back(im);
            //	cout << "pt:" << analysisTree.muon_pt[im] << "  passed:" << muPassed << endl;
            //cout << "eta:" << analysisTree.muon_eta[im] << endl;
        }
        
        //tau selection
        vector<unsigned int> idTaus; idTaus.clear(); //identified taus
        vector<unsigned int> allTaus; allTaus.clear(); //all taus
        for(unsigned int it=0;it<analysisTree.tau_count;++it)
        {
            allTaus.push_back(it);
            if (analysisTree.tau_pt[it]<ptTauCut) continue;
            if (fabs(analysisTree.tau_eta[it])>etaTauCut) continue;
            if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it])>dzTauCut) continue;
            if (analysisTree.tau_decayModeFindingNewDMs[it]<=0.5) continue;
            if (analysisTree.tau_decayMode[it] == 5 || analysisTree.tau_decayMode[it] == 6) continue;
            //(no needed for Synchronization)
            //if (analysisTree.tau_byMediumDeepTau2017v2p1VSjet[it] <= 0.5) continue;
            if (analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[it] <= 0.5) continue; //For synch! Because nanoAOD has it!
            if (analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[it] <= 0.5) continue;
            if (analysisTree.tau_byVLooseDeepTau2017v2p1VSe[it] <= 0.5) continue;
            idTaus.push_back(it);
        }
        


        vector<unsigned int> goodEles; goodEles.clear(); //good electrons
        vector<unsigned int> goodMuons; goodMuons.clear(); //good muons
        vector<unsigned int> goodTaus; goodTaus.clear(); //good taus
        goodEles = idEles;
        goodMuons = idMuons;
        goodTaus = idTaus;

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
            
            if(era == "2016"||era == "UL2016postVFP"||era == "UL2016preVFP")
            {
                if(isMuon1DoubleMuon_3_Leg1Matched && (isMuon2DoubleMuon_3_Leg2_a_Matched || isMuon2DoubleMuon_3_Leg2_b_Matched) && isMuon1DoubleMuon_3_dzMatched && isMuon2DoubleMuon_3_dzMatched)
                    isDoubleMuon_3_TriggerFilterMatched = true;
            }
            else
            {
                if(isMuon1DoubleMuon_3_Leg1Matched && (isMuon2DoubleMuon_3_Leg2_a_Matched || isMuon2DoubleMuon_3_Leg2_b_Matched) && isMuon1DoubleMuon_3_dzMatched && isMuon2DoubleMuon_3_dzMatched && isMuon1DoubleMuon_3_massMatched && isMuon2DoubleMuon_3_massMatched)
                    isDoubleMuon_3_TriggerFilterMatched = true;
            }
            if(era == "2016"||era == "UL2016postVFP"||era == "UL2016preVFP")
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
                if(era == "2016"|| era == "UL2016postVFP"|| era == "UL2016preVFP")
                {
                    if((isSingleMuon_1_TriggerFilterMatched == true) || (isSingleMuon_2_TriggerFilterMatched == true))
                        isSingleMuonTriggerFilterMatched = true;
                }
            }
            else//"OR" double trigger apply
            {
                if((era == "2018") || (era == "2017")||(era == "UL2018") || (era == "UL2017"))
                {
                    if(isSingleMuon_1_TriggerFilterMatched == true)
                        isSingleMuonTriggerFilterMatched = true;
                        
                    if (isDoubleMuon_3_TriggerFilterMatched == true)
                        isDoubleMuonTriggerFilterMatched = true;
                }
                if(era == "2016"|| era == "UL2016postVFP"|| era == "UL2016preVFP")
                {
                    if((isSingleMuon_1_TriggerFilterMatched==true) || (isSingleMuon_2_TriggerFilterMatched==true))
                        isSingleMuonTriggerFilterMatched = true;
                    
                    if((isDoubleMuon_3_TriggerFilterMatched==true) || (isDoubleMuon_4_TriggerFilterMatched==true))
                        isDoubleMuonTriggerFilterMatched = true;
                }
            }
        }
 
        //Higgs boson candidates selection
        int indexLeg3 = -1;
        int indexLeg4 = -1;
        float Higgs_LT = -9999;
        if(HiggsFinalState =="EM")
        for(unsigned int ie=0;ie<goodEles.size();++ie)
        {
                int indexTemp1 = goodEles[ie];
                for(unsigned int im=0;im<goodMuons.size();++im)
                {
                    int indexTemp2 = goodMuons[im];
                    if ((indexTemp2 == indexMuon1)||(indexTemp2 == indexMuon2))
                        continue;
                    if (analysisTree.electron_pt[indexTemp1]+analysisTree.muon_pt[indexTemp2]>Higgs_LT)
                    {
                        Higgs_LT = analysisTree.electron_pt[indexTemp1] + analysisTree.muon_pt[indexTemp2];
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp2;
                    }
                }
        }
        if(HiggsFinalState =="ET")
        {
            myfile << "-------------------------" <<std::endl;
            myfile << "--Checking DR of taus" <<std::endl;
            myfile << "-------------------------" <<std::endl;
            for(unsigned int ie=0;ie<goodEles.size();++ie)
            {
                int indexTemp1 = goodEles[ie];
                for(unsigned int it=0;it<goodTaus.size();++it)
                {
                    int indexTemp2 = goodTaus[it];
                    if (analysisTree.tau_byTightDeepTau2017v2p1VSe[indexTemp2]<0.5) continue;
                    //(no needed for Synchronization)
                    //if (analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexTemp2]<0.5) continue;
                    float DR_1Tau = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    float DR_2Tau = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    float DR_EleTau = deltaR(analysisTree.electron_eta[indexTemp1],analysisTree.electron_phi[indexTemp1],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    myfile << "DR between Z leg 1 and tau #" << indexTemp2 << ": " << DR_1Tau << std::endl;
                    myfile << "DR between Z leg 2 and tau #" << indexTemp2 << ": " << DR_2Tau << std::endl;
                    myfile << "DR between electron #" << indexTemp1 << " and tau #" << indexTemp2 << ": " << DR_EleTau << std::endl;
                    if((DR_1Tau < 0.5)||(DR_2Tau < 0.5)||(DR_EleTau < 0.5))//Select tau that is not overlapped with others first, then LT...
                    {
                        continue;
                    }
                    if (analysisTree.electron_pt[indexTemp1]+analysisTree.tau_pt[indexTemp2]>Higgs_LT)
                    {
                        Higgs_LT = analysisTree.electron_pt[indexTemp1] + analysisTree.tau_pt[indexTemp2];
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp2;
                    }
                }
            }
        }
        
        if(HiggsFinalState =="MT")
        {
            myfile << "-------------------------" <<std::endl;
            myfile << "--Checking DR of taus" <<std::endl;
            myfile << "-------------------------" <<std::endl;
            for(unsigned int im=0;im<goodMuons.size();++im)
            {
                int indexTemp1 = goodMuons[im];
                if ((indexTemp1 == indexMuon1)||(indexTemp1 == indexMuon2))
                    continue;
                for(unsigned int it=0;it<goodTaus.size();++it)
                {
                    int indexTemp2 = goodTaus[it];
                    if (analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexTemp2]<0.5) continue;
                    float DR_1Tau = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    float DR_2Tau = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    float DR_MuTau = deltaR(analysisTree.muon_eta[indexTemp1],analysisTree.muon_phi[indexTemp1],analysisTree.tau_eta[indexTemp2],analysisTree.tau_phi[indexTemp2]);
                    myfile << "DR between Z leg 1 and tau #" << indexTemp2 << ": " << DR_1Tau << std::endl;
                    myfile << "DR between Z leg 2 and tau #" << indexTemp2 << ": " << DR_2Tau << std::endl;
                    myfile << "DR between muon #" << indexTemp1 << " and tau #" << indexTemp2 << ": " << DR_MuTau << std::endl;
                    if((DR_1Tau < 0.5)||(DR_2Tau < 0.5)||(DR_MuTau < 0.5))//Select tau that is not overlapped with others first, then LT...
                    {
                        continue;
                    }
                    if (analysisTree.muon_pt[indexTemp1]+analysisTree.tau_pt[indexTemp2]>Higgs_LT)
                    {
                        Higgs_LT = analysisTree.muon_pt[indexTemp1] + analysisTree.tau_pt[indexTemp2];
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp2;
                    }
                }
            }
        }
        
        if(HiggsFinalState =="TT")
        {
            myfile << "-------------------------" <<std::endl;
            myfile << "--Checking DR of taus" <<std::endl;
            myfile << "-------------------------" <<std::endl;
            for(unsigned int it=0;it<goodTaus.size();++it)
            {
                int indexTemp1 = goodTaus[it];
                float q1 = analysisTree.tau_charge[indexTemp1];
                float DR_1Tau1 = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.tau_eta[indexTemp1],analysisTree.tau_phi[indexTemp1]);
                float DR_2Tau1 = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.tau_eta[indexTemp1],analysisTree.tau_phi[indexTemp1]);
                myfile << "DR between Z leg 1 and tau #" << indexTemp1 << ": " << DR_1Tau1 << std::endl;
                myfile << "DR between Z leg 2 and tau #" << indexTemp1 << ": " << DR_2Tau1 << std::endl;
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
                    myfile << "DR between Z leg 1 and tau #" << indexTemp2 << ": " << DR_1Tau2 << std::endl;
                    myfile << "DR between Z leg 2 and tau #" << indexTemp2 << ": " << DR_2Tau2 << std::endl;
                    myfile << "DR between tau #"<< indexTemp1 << " and tau #" << indexTemp2 << ": " << DR_TauTau << std::endl;
                    if((DR_1Tau2 < 0.5)||(DR_2Tau2 < 0.5)||(DR_TauTau < 0.5))//Select tau that is not overlapped with others first, then LT...
                    {
                        continue;
                    }
                    if (analysisTree.tau_pt[indexTemp1]+analysisTree.tau_pt[indexTemp2]>Higgs_LT)
                    {
                        Higgs_LT = analysisTree.tau_pt[indexTemp1] + analysisTree.tau_pt[indexTemp2];
                        indexLeg3 = indexTemp1;
                        indexLeg4 = indexTemp2;
                    }
                }
            }
        }

        //Filling up lepton kinematics
        pt_1 = analysisTree.muon_pt[indexMuon1];
        phi_1 = analysisTree.muon_phi[indexMuon1];
        eta_1 = analysisTree.muon_eta[indexMuon1];
        
        float absIso_1 = analysisTree.muon_chargedHadIso[indexMuon1];
        float neutralIso_1 = analysisTree.muon_neutralHadIso[indexMuon1] + analysisTree.muon_photonIso[indexMuon1] - 0.5*analysisTree.muon_puIso[indexMuon1];
        neutralIso_1 = TMath::Max(float(0),neutralIso_1);
        absIso_1 += neutralIso_1;
        iso_1 = absIso_1/analysisTree.muon_pt[indexMuon1];
        if(!isData)
        {
            muonIDSF_1 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[indexMuon1]),double(analysisTree.muon_eta[indexMuon1]));
            if(isSingleMuonTriggerFilterMatched && isMuon1SingleMuonMatched)
                muonTrigSF_1 = (float)SF_muonTrig->get_ScaleFactor(double(analysisTree.muon_pt[indexMuon1]),double(analysisTree.muon_eta[indexMuon1]));
        }
        muonLoose_1 = analysisTree.muon_isLoose[indexMuon1];
        muonMedium_1 = analysisTree.muon_isMedium[indexMuon1];
        muonTight_1 = analysisTree.muon_isTight[indexMuon1];
        
        pt_2 = analysisTree.muon_pt[indexMuon2];
        phi_2 = analysisTree.muon_phi[indexMuon2];
        eta_2 = analysisTree.muon_eta[indexMuon2];
        
        float absIso_2 = analysisTree.muon_chargedHadIso[indexMuon2];
        float neutralIso_2 = analysisTree.muon_neutralHadIso[indexMuon2] + analysisTree.muon_photonIso[indexMuon2] - 0.5*analysisTree.muon_puIso[indexMuon2];
        neutralIso_2 = TMath::Max(float(0),neutralIso_2);
        absIso_2 += neutralIso_2;
        iso_2 = absIso_2/analysisTree.muon_pt[indexMuon2];
        if(!isData)
        {
            muonIDSF_2 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[indexMuon2]),double(analysisTree.muon_eta[indexMuon2]));
            if(isSingleMuonTriggerFilterMatched && isMuon2SingleMuonMatched)
                muonTrigSF_2 = (float)SF_muonTrig->get_ScaleFactor(double(analysisTree.muon_pt[indexMuon2]),double(analysisTree.muon_eta[indexMuon2]));
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
        if(HiggsFinalState == "EM")
        {
            pt_3 = analysisTree.electron_pt[indexLeg3];
            phi_3 = analysisTree.electron_phi[indexLeg3];
            eta_3 = analysisTree.electron_eta[indexLeg3];
            m_3 = classic_svFit::electronMass;
            float absIso_3 = analysisTree.electron_r03_sumChargedHadronPt[indexLeg3];
            float neutralIso_3 = analysisTree.electron_r03_sumNeutralHadronEt[indexLeg3] + analysisTree.electron_r03_sumPhotonEt[indexLeg3] - 0.5*analysisTree.electron_r03_sumPUPt[indexLeg3];
            neutralIso_3 = TMath::Max(float(0),neutralIso_3);
            absIso_3 += neutralIso_3;
            iso_3 = absIso_3/analysisTree.electron_pt[indexLeg3];
            iso_ea_3 = analysisTree.electron_eaIsolation[indexLeg3]/analysisTree.electron_pt[indexLeg3];

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
                muonIDSF_4 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[indexLeg4]),double(analysisTree.muon_eta[indexLeg4]));
            }
            Leg4.SetXYZM(analysisTree.muon_px[indexLeg4],analysisTree.muon_py[indexLeg4],analysisTree.muon_pz[indexLeg4],muonMass);
            muonLoose_4 = analysisTree.muon_isLoose[indexLeg4];
            muonMedium_4 = analysisTree.muon_isMedium[indexLeg4];
            muonTight_4 = analysisTree.muon_isTight[indexLeg4];
        }
        if(HiggsFinalState == "ET")
        {
            pt_3 = analysisTree.electron_pt[indexLeg3];
            phi_3 = analysisTree.electron_phi[indexLeg3];
            eta_3 = analysisTree.electron_eta[indexLeg3];
            m_3 = classic_svFit::electronMass;
            float absIso_3 = analysisTree.electron_r03_sumChargedHadronPt[indexLeg3];
            float neutralIso_3 = analysisTree.electron_r03_sumNeutralHadronEt[indexLeg3] + analysisTree.electron_r03_sumPhotonEt[indexLeg3] - 0.5*analysisTree.electron_r03_sumPUPt[indexLeg3];
            neutralIso_3 = TMath::Max(float(0),neutralIso_3);
            absIso_3 += neutralIso_3;
            iso_3 = absIso_3/analysisTree.electron_pt[indexLeg3];
            iso_ea_3 = analysisTree.electron_eaIsolation[indexLeg3]/analysisTree.electron_pt[indexLeg3];
            Leg3.SetXYZM(analysisTree.electron_px[indexLeg3],analysisTree.electron_py[indexLeg3],analysisTree.electron_pz[indexLeg3],electronMass);
            electronIDWP80v2_3 = analysisTree.electron_mva_wp80_noIso_Fall17_v2[indexLeg3];
            electronIDWP90v2_3 = analysisTree.electron_mva_wp90_noIso_Fall17_v2[indexLeg3];

            pt_4 = analysisTree.tau_pt[indexLeg4];
            phi_4 = analysisTree.tau_phi[indexLeg4];
            eta_4 = analysisTree.tau_eta[indexLeg4];
            m_4 = analysisTree.tau_mass[indexLeg4];
            dm_4 = analysisTree.tau_decayMode[indexLeg4];
            iso_4 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexLeg4];
            Leg4.SetXYZM(analysisTree.tau_px[indexLeg4],analysisTree.tau_py[indexLeg4],analysisTree.tau_pz[indexLeg4],analysisTree.tau_mass[indexLeg4]);
            
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
                muonIDSF_3 = (float)SF_muonIdIso->get_ScaleFactor(double(analysisTree.muon_pt[indexLeg3]),double(analysisTree.muon_eta[indexLeg3]));
            }
            Leg3.SetXYZM(analysisTree.muon_px[indexLeg3],analysisTree.muon_py[indexLeg3],analysisTree.muon_pz[indexLeg3],muonMass);
            muonLoose_3 = analysisTree.muon_isLoose[indexLeg3];
            muonMedium_3 = analysisTree.muon_isMedium[indexLeg3];
            muonTight_3 = analysisTree.muon_isTight[indexLeg3];
            
            pt_4 = analysisTree.tau_pt[indexLeg4];
            phi_4 = analysisTree.tau_phi[indexLeg4];
            eta_4 = analysisTree.tau_eta[indexLeg4];
            m_4 = analysisTree.tau_mass[indexLeg4];
            dm_4 = analysisTree.tau_decayMode[indexLeg4];
            iso_4 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexLeg4];
            Leg4.SetXYZM(analysisTree.tau_px[indexLeg4],analysisTree.tau_py[indexLeg4],analysisTree.tau_pz[indexLeg4],analysisTree.tau_mass[indexLeg4]);
            
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
        if(HiggsFinalState == "TT")
        {
            pt_3 = analysisTree.tau_pt[indexLeg3];
            phi_3 = analysisTree.tau_phi[indexLeg3];
            eta_3 = analysisTree.tau_eta[indexLeg3];
            m_3 = analysisTree.tau_mass[indexLeg3];
            dm_3 = analysisTree.tau_decayMode[indexLeg3];
            iso_3 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexLeg3];
            Leg3.SetXYZM(analysisTree.tau_px[indexLeg3],analysisTree.tau_py[indexLeg3],analysisTree.tau_pz[indexLeg3],analysisTree.tau_mass[indexLeg3]);

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
            
            pt_4 = analysisTree.tau_pt[indexLeg4];
            phi_4 = analysisTree.tau_phi[indexLeg4];
            eta_4 = analysisTree.tau_eta[indexLeg4];
            m_4 = analysisTree.tau_mass[indexLeg4];
            dm_4 = analysisTree.tau_decayMode[indexLeg4];
            iso_4 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexLeg4];
            Leg4.SetXYZM(analysisTree.tau_px[indexLeg4],analysisTree.tau_py[indexLeg4],analysisTree.tau_pz[indexLeg4],analysisTree.tau_mass[indexLeg4]);
            
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
        
        myfile << "-------------------------" <<std::endl;
        myfile << "--Electron information " <<std::endl;
        myfile << "-------------------------" <<std::endl;

        for(unsigned int ie=0;ie<allEles.size();++ie)
        {
            int indexPrintEle = allEles[ie];
            myfile << "--Electron number: " << indexPrintEle << std::endl;
            if((HiggsFinalState == "ET") || (HiggsFinalState == "EM"))
            {
                if(indexPrintEle == indexLeg3)
                    myfile << "--It is selected as Higgs leg1" << std::endl;
            }
            myfile << "--pT                 : " << analysisTree.electron_pt[indexPrintEle] << std::endl;
            myfile << "--eta                : " << analysisTree.electron_eta[indexPrintEle] << std::endl;
            myfile << "--phi                : " << analysisTree.electron_phi[indexPrintEle] << std::endl;
            myfile << "--dxy                : " << analysisTree.electron_dxy[indexPrintEle] << std::endl;
            myfile << "--dz                 : " << analysisTree.electron_dz[indexPrintEle] << std::endl;
            myfile << "--MVA value noIso v2 : " << analysisTree.electron_mva_value_noIso_Fall17_v2[indexPrintEle] << std::endl;
            myfile << "--MVA WP90 noIso v2  : " << analysisTree.electron_mva_wp90_noIso_Fall17_v2[indexPrintEle] << std::endl;
            myfile << "--MVA WP80 noIso v2  : " << analysisTree.electron_mva_wp80_noIso_Fall17_v2[indexPrintEle] << std::endl;
            myfile << "--Pass Coversion     : " << analysisTree.electron_pass_conversion[indexPrintEle] << std::endl;
            myfile << "--N Missing Hits     : " << (int) analysisTree.electron_nmissinghits[indexPrintEle] << std::endl;
            myfile << "--Charge             : " << analysisTree.electron_charge[indexPrintEle] << std::endl;
            myfile << "-------------------------" << std::endl;
        }
        
        myfile << "-------------------------" <<std::endl;
        myfile << "--Muon information " <<std::endl;
        myfile << "-------------------------" <<std::endl;

        for(unsigned int im=0;im<allMuons.size();++im)
        {
            int indexPrintMu = allMuons[im];
            myfile << "--Muon number    : " << indexPrintMu << std::endl;
            if(indexPrintMu == indexMuon1)
                myfile << "--It is selected as Z leg1" << std::endl;
            if(indexPrintMu == indexMuon2)
                myfile << "--It is selected as Z leg2" << std::endl;
            if(HiggsFinalState == "MT")
            {
                if(indexPrintMu == indexLeg3)
                    myfile << "--It is selected as Higgs leg1" << std::endl;
            }
            if(HiggsFinalState == "EM")
            {
                if(indexPrintMu == indexLeg4)
                    myfile << "--It is selected as Higgs leg2" << std::endl;
            }
            myfile << "--pT                 : " << analysisTree.muon_pt[indexPrintMu] << std::endl;
            myfile << "--eta                : " << analysisTree.muon_eta[indexPrintMu] << std::endl;
            myfile << "--phi                : " << analysisTree.muon_phi[indexPrintMu] << std::endl;
            myfile << "--dxy                : " << analysisTree.muon_dxy[indexPrintMu] << std::endl;
            myfile << "--dz                 : " << analysisTree.muon_dz[indexPrintMu] << std::endl;
            myfile << "--Loose id           : " << analysisTree.muon_isLoose[indexPrintMu] << std::endl;
            myfile << "--Global Muon        : " << analysisTree.muon_isGlobal[indexPrintMu] << std::endl;
            myfile << "--Tracker Muon       : " << analysisTree.muon_isTracker[indexPrintMu] << std::endl;
            myfile << "--Charge             : " << analysisTree.muon_charge[indexPrintMu] << std::endl;
            myfile << "-------------------------" << std::endl;
        }
        myfile << "-------------------------" <<std::endl;
        myfile << "--Tau information " <<std::endl;
        myfile << "-------------------------" <<std::endl;
        for(unsigned int it=0;it<allTaus.size();++it)
        {
            int indexPrintTau = allTaus[it];
            myfile << "--Tau number     : " << indexPrintTau << std::endl;
            if(HiggsFinalState == "TT")
            {
                if(indexPrintTau == indexLeg3)
                    myfile << "--It is selected as Higgs leg1" << std::endl;
            }
            if((HiggsFinalState == "MT")||(HiggsFinalState == "ET")||(HiggsFinalState == "TT"))
            {
                if(indexPrintTau == indexLeg4)
                    myfile << "--It is selected as Higgs leg2" << std::endl;
            }
            myfile << "--pT                 : " << analysisTree.tau_pt[indexPrintTau] << std::endl;
            myfile << "--eta                : " << analysisTree.tau_eta[indexPrintTau] << std::endl;
            myfile << "--phi                : " << analysisTree.tau_phi[indexPrintTau] << std::endl;
            myfile << "--dz(leading Track)  : " << analysisTree.tau_leadchargedhadrcand_dz[indexPrintTau] << std::endl;
            myfile << "--NewDM              : " << analysisTree.tau_decayModeFindingNewDMs[indexPrintTau] << std::endl;
            myfile << "--Decay mode         : " << analysisTree.tau_decayMode[indexPrintTau] << std::endl;
            myfile << "--VSjet score        : " << analysisTree.tau_byDeepTau2017v2p1VSjetraw[indexPrintTau] << std::endl;
            myfile << "--VSe score          : " << analysisTree.tau_byDeepTau2017v2p1VSeraw[indexPrintTau] << std::endl;
            myfile << "--VSmu score         : " << analysisTree.tau_byDeepTau2017v2p1VSmuraw[indexPrintTau] << std::endl;
            myfile << "--VVVLoose VSjet     : " << analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[indexPrintTau] << std::endl;
            myfile << "--Medium VSjet     : " << analysisTree.tau_byMediumDeepTau2017v2p1VSjet[indexPrintTau] << std::endl;
            myfile << "--Tight VSe         : " << analysisTree.tau_byTightDeepTau2017v2p1VSe[indexPrintTau] << std::endl;
            myfile << "--VLoose VSe         : " << analysisTree.tau_byVLooseDeepTau2017v2p1VSe[indexPrintTau] << std::endl;
            myfile << "--Tight VSmu         : " << analysisTree.tau_byTightDeepTau2017v2p1VSmu[indexPrintTau] << std::endl;
            myfile << "--VLoose VSmu         : " << analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[indexPrintTau] << std::endl;
            myfile << "--Charge             : " << analysisTree.tau_charge[indexPrintTau] << std::endl;
            myfile << "-------------------------" << std::endl;
        }
        
        myfile << "-------------------------" <<std::endl;
        myfile << "--Other information " <<std::endl;
        myfile << "-------------------------" <<std::endl;
        myfile << "--SingleMuon Trigger:" << std::endl;
        myfile << "--The trigger filter    :" << SingleMuonFilterName_1_a << " in HLT filter #" << nSingleMuon_1_Filter_a << std::endl;
        for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
        {
            float dRtrig = deltaR(analysisTree.muon_eta[indexMuon1],analysisTree.muon_phi[indexMuon1],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
            if (dRtrig<DRTrigMatch)
            {
                myfile << "Trigger object #" << iT << " is matched to Z leg1 with DR " << dRtrig << std::endl;
                myfile << "Filter is :" << analysisTree.trigobject_filters[iT][nSingleMuon_1_Filter_a] << std::endl;
            }
        }
        myfile << "Z leg1 pT :"<< analysisTree.muon_pt[indexMuon1] << std::endl;
        myfile << "Z leg1 eta :" << analysisTree.muon_eta[indexMuon1] << std::endl;
        myfile << "-------------------------" <<std::endl;
        for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT)
        {
            float dRtrig = deltaR(analysisTree.muon_eta[indexMuon2],analysisTree.muon_phi[indexMuon2],analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
            if (dRtrig<DRTrigMatch)
            {
                myfile << "Trigger object #" << iT << " is matched to Z leg2 with DR " << dRtrig << std::endl;
                myfile << "Filter is :" << analysisTree.trigobject_filters[iT][nSingleMuon_1_Filter_a] << std::endl;
            }
        }
        myfile << "Z leg2 pT :"<< analysisTree.muon_pt[indexMuon2] << std::endl;
        myfile << "Z leg2 eta :" << analysisTree.muon_eta[indexMuon2] << std::endl;
        myfile << "-------------------------" <<std::endl;

        if(!((isMuon1SingleMuon_1_a_Matched || isMuon1SingleMuon_1_b_Matched)||(isMuon2SingleMuon_1_a_Matched || isMuon2SingleMuon_1_b_Matched)))
        {
            myfile << "--Either no SingleMuon trigger filter matched to any Z leg or passed pT and eta cuts!" << endl;
        }
        if((isMuon1SingleMuon_1_a_Matched || isMuon1SingleMuon_1_b_Matched))
        {
            myfile << "--Muon " << indexMuon1 << " matched SingleMuon filter and passed pT and eta cuts" << endl;
        }
        if((isMuon2SingleMuon_1_a_Matched || isMuon2SingleMuon_1_b_Matched))
        {
            myfile << "--Muon " << indexMuon2 << " matched SingleMuon filter and passed pT and eta cuts" << endl;

        }
        myfile << "-------------------------" <<std::endl;


        //uncorrected met
        met = TMath::Sqrt(analysisTree.pfmetcorr_ex*analysisTree.pfmetcorr_ex + analysisTree.pfmetcorr_ey*analysisTree.pfmetcorr_ey);
        metphi = analysisTree.pfmetcorr_phi;
        metcov00 = analysisTree.pfmetcorr_sigxx;
        metcov01 = analysisTree.pfmetcorr_sigxy;
        metcov10 = analysisTree.pfmetcorr_sigyx;
        metcov11 = analysisTree.pfmetcorr_sigyy;

        TLorentzVector metLV;
        metLV.SetXYZT(met*TMath::Cos(metphi),met*TMath::Sin(metphi),0,TMath::Sqrt(met*TMath::Sin(metphi)*met*TMath::Sin(metphi) + met*TMath::Cos(metphi)*met*TMath::Cos(metphi)));
        
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
        
        myfile << "--mll                   :" << mll <<std::endl;
        
        myfile << "--# of baseline Electron:" << goodEles.size() << std::endl;
        myfile << "--# of baseline Muon    :" << goodMuons.size() << std::endl;
        myfile << "--# of baseline Tau     :" << goodTaus.size() << std::endl;

        if(foundZOS)
        {
            myfile << "--Z candidate found!" << std::endl;
            myfile << "--DR of leg1 and leg2   :" << Z_DR <<std::endl;
        }
        else
        {
            myfile << "--Z candidate NOT found!" << std::endl;

        }
        if(indexLeg3!= (-1) && indexLeg4!=(-1))
        {
            myfile << "--Higgs candidate found!" << std::endl;
        }
        else
        {
            myfile << "--Higgs candidate NOT found!" << std::endl;

        }
        if(indexMuon1!= (-1) && indexMuon2!=(-1)&&indexLeg3!= (-1) && indexLeg4!=(-1))
        {
            myfile << "--A candidate found!" << std::endl;
            myfile << "--DR of leg1 and leg3   :" << DR_13 <<std::endl;
            myfile << "--DR of leg1 and leg4   :" << DR_14 <<std::endl;
            myfile << "--DR of leg2 and leg3   :" << DR_23 <<std::endl;
            myfile << "--DR of leg2 and leg4   :" << DR_24 <<std::endl;
            myfile << "--DR of leg3 and leg4   :" << H_DR <<std::endl;
            if(HiggsFinalState=="MT")
            {
                myfile << "--Extra # of baseline Electron   :" << goodEles.size() - 0 <<std::endl;
                myfile << "--Extra # of baseline Muon       :" << goodMuons.size() - 3 <<std::endl;
                myfile << "--Extra # of baseline Tau        :" << goodTaus.size() - 1 <<std::endl;
            }
            if(HiggsFinalState=="ET")
            {
                myfile << "--Extra # of baseline Electron   :" << goodEles.size() - 1 <<std::endl;
                myfile << "--Extra # of baseline Muon       :" << goodMuons.size() - 2 <<std::endl;
                myfile << "--Extra # of baseline Tau        :" << goodTaus.size() - 1 <<std::endl;
            }
            if(HiggsFinalState=="EM")
            {
                myfile << "--Extra # of baseline Electron   :" << goodEles.size() - 1 <<std::endl;
                myfile << "--Extra # of baseline Muon       :" << goodMuons.size() - 3 <<std::endl;
                myfile << "--Extra # of baseline Tau        :" << goodTaus.size() - 0 <<std::endl;
            }
            if(HiggsFinalState=="TT")
            {
                myfile << "--Extra # of baseline Electron   :" << goodEles.size() - 0 <<std::endl;
                myfile << "--Extra # of baseline Muon       :" << goodMuons.size() - 2 <<std::endl;
                myfile << "--Extra # of baseline Tau        :" << goodTaus.size() - 2 <<std::endl;
            }
            myfile << "--iso 1   :" << iso_1 <<std::endl;
            myfile << "--iso 2   :" << iso_2 <<std::endl;
            myfile << "--iso 3   :" << iso_3 <<std::endl;
            myfile << "--iso 4   :" << iso_4 <<std::endl;
            
            myfile << "--iso ea 1   :" << iso_ea_1 <<std::endl;
            myfile << "--iso ea 2   :" << iso_ea_2 <<std::endl;
            myfile << "--iso ea 3   :" << iso_ea_3 <<std::endl;
            myfile << "--iso ea 4   :" << iso_ea_4 <<std::endl;
        }
        else
        {
            myfile << "--A candidate NOT found!" << std::endl;
        }
        myfile << "-------------------------" <<std::endl;
        myfile << "" <<std::endl;
        myfile << "" <<std::endl;
        myfile << "" <<std::endl;

        

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
            
            if(HiggsFinalState == "MT")
            {
                classic_svFit::MeasuredTauLepton svFitMu(classic_svFit::MeasuredTauLepton::kTauToMuDecay, pt_3, eta_3, phi_3, m_3);
                classic_svFit::MeasuredTauLepton svFitTau(classic_svFit::MeasuredTauLepton::kTauToHadDecay, pt_4, eta_4, phi_4, m_4, dm_4);
                measuredTauLeptons.push_back(svFitMu);
                measuredTauLeptons.push_back(svFitTau);
            }
            
            if(HiggsFinalState == "ET")
            {
                classic_svFit::MeasuredTauLepton svFitEle(classic_svFit::MeasuredTauLepton::kTauToElecDecay, pt_3, eta_3, phi_3, m_3);
                classic_svFit::MeasuredTauLepton svFitTau(classic_svFit::MeasuredTauLepton::kTauToHadDecay, pt_4, eta_4, phi_4, m_4, dm_4);
                measuredTauLeptons.push_back(svFitEle);
                measuredTauLeptons.push_back(svFitTau);
            }
            if(HiggsFinalState == "EM")
            {
                classic_svFit::MeasuredTauLepton svFitEle(classic_svFit::MeasuredTauLepton::kTauToElecDecay, pt_3, eta_3, phi_3, m_3);
                classic_svFit::MeasuredTauLepton svFitMu(classic_svFit::MeasuredTauLepton::kTauToMuDecay, pt_4, eta_4, phi_4, m_4);
                measuredTauLeptons.push_back(svFitEle);
                measuredTauLeptons.push_back(svFitMu);
            }
            if(HiggsFinalState == "TT")
            {
                classic_svFit::MeasuredTauLepton svFitTau1(classic_svFit::MeasuredTauLepton::kTauToHadDecay, pt_3, eta_3, phi_3, m_3, dm_3);
                classic_svFit::MeasuredTauLepton svFitTau2(classic_svFit::MeasuredTauLepton::kTauToHadDecay, pt_4, eta_4, phi_4, m_4, dm_4);
                measuredTauLeptons.push_back(svFitTau1);
                measuredTauLeptons.push_back(svFitTau2);
            }
            if (ApplySVFit)
            {
                int verbosity = 0;
                ClassicSVfit svFitAlgo(verbosity);
                double kappa = 0.; // use 3 for emu, 4 for etau and mutau, 5 for tautau channel
                if((HiggsFinalState == "MT")||(HiggsFinalState == "ET"))
                    kappa = 4.;
                if(HiggsFinalState == "EM")
                    kappa = 3.;
                if(HiggsFinalState == "TT")
                    kappa = 5.;
                
                svFitAlgo.addLogM_fixed(true, kappa);
                //svFitAlgo.setIntegrationParams(true);
                svFitAlgo.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
                //bool isValidSolution = svFitAlgo.isValidSolution();
                m_sv   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getMass();
                pt_sv  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getPt();
                eta_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getEta();
                phi_sv = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getPhi();
                //met_sv = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getfittedMET().Rho();
                mt_sv  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
                
                TLorentzVector ditau_sv;
                ditau_sv.SetPtEtaPhiM(pt_sv, eta_sv, phi_sv, m_sv);
                Mass_sv = (muon1+muon2+ditau_sv).M();
                
                ClassicSVfit svFitAlgo_C(verbosity);//define mass constraint SV fit
                svFitAlgo_C.addLogM_fixed(true, kappa);
                //the most accurate number 125.38 +/- 0.14 from Physics Letters B Volume 805, 10 June 2020, 135425
                svFitAlgo_C.setDiTauMassConstraint(125.38);
                svFitAlgo_C.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
                m_svc   = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C.getHistogramAdapter())->getMass();
                pt_svc  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C.getHistogramAdapter())->getPt();
                eta_svc = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C.getHistogramAdapter())->getEta();
                phi_svc = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C.getHistogramAdapter())->getPhi();
                //met_svc = static_cast<DiTauSystemHistogramAdapter*>(svFitAlgo_C.getHistogramAdapter())->getfittedMET().Rho();
                mt_svc  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_C.getHistogramAdapter())->getTransverseMass();
                
                TLorentzVector ditau_svc;
                ditau_svc.SetPtEtaPhiM(pt_svc, eta_svc, phi_svc, m_svc);
                Mass_svc = (muon1+muon2+ditau_svc).M();
                
            }
            if (ApplyFastMTT)
            {
                // FasMTT
                LorentzVector tau1P4;
                LorentzVector tau2P4;
                FastMTT aFastMTTAlgo;
                  
                aFastMTTAlgo.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
                LorentzVector ttP4 = aFastMTTAlgo.getBestP4();
                tau1P4 = aFastMTTAlgo.getTau1P4();
                tau2P4 = aFastMTTAlgo.getTau2P4();

                double dPhiTT = dPhiFrom2P( tau1P4.Px(), tau1P4.Py(), tau2P4.Px(), tau2P4.Py() );
                mt_fast = TMath::Sqrt(2*tau1P4.Pt()*tau2P4.Pt()*(1 - TMath::Cos(dPhiTT)));
                m_fast = ttP4.M();
                pt_fast = ttP4.Pt();
                eta_fast = ttP4.Eta();
                phi_fast = ttP4.Phi();
                
                TLorentzVector ditau_fast;
                ditau_fast.SetPtEtaPhiM(pt_fast, eta_fast, phi_fast, m_fast);
                Mass_fast = (muon1+muon2+ditau_fast).M();
            }
        }
        
        // Jets selection
        vector<unsigned int> jets; jets.clear();
        vector<unsigned int> jetspt20; jetspt20.clear();
        vector<unsigned int> bjets; bjets.clear();
        vector<unsigned int> bjets_nocleaned; bjets_nocleaned.clear();
        vector<unsigned int> bjetsRaw; bjetsRaw.clear();
            
        int indexLeadingJet = -1;
        float ptLeadingJet = -1;
            
        int indexSubLeadingJet = -1;
        float ptSubLeadingJet = -1;
            
        int indexLeadingBJet = -1;
        float ptLeadingBJet = -1;
        
        int indexSubLeadingBJet = -1;
        float ptSubLeadingBJet = -1;
        for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet)
        {
            float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
            float jetEta = analysisTree.pfjet_eta[jet];
            if (absJetEta>jetEtaCut) continue;
                
            float jetPt = analysisTree.pfjet_pt[jet];
            //std::cout << jet << " : uncertainty = " << analysisTree.pfjet_jecUncertainty[jet] << std::endl;
            if (jetPt<jetPtLowCut) continue;
                
            bool isPFJetId = looseJetiD(analysisTree,int(jet));
            if (!isPFJetId) continue;
            // non-cleaned jet at the moment
            
            // jetId
            if (jetPt>jetPtLowCut)
                jetspt20.push_back(jet);
            
            if (absJetEta<bJetEtaCut)
            {   // jet within b-tagging acceptance
                    
                bool tagged = analysisTree.pfjet_btag[jet][nBTagDiscriminant]>btagCut; // b-jet
                //missing b-tag SFs

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
        
        njets = jets.size();
        
        int njetsMax = njets;

        njetspt20 = jetspt20.size();
        nbtag = bjets.size();

        //missing b-tag weight
        
        bpt_1 = -9999;
        beta_1 = -9999;
        bphi_1 = -9999;
        bhadronFlavor_1 = -1;
        
        bpt_2 = -9999;
        beta_2 = -9999;
        bphi_2 = -9999;
        bhadronFlavor_2 = -1;

        if (indexLeadingBJet>=0)
        {
            bpt_1 = analysisTree.pfjet_pt[indexLeadingBJet];
            beta_1 = analysisTree.pfjet_eta[indexLeadingBJet];
            bphi_1 = analysisTree.pfjet_phi[indexLeadingBJet];
        }
        if(!isData)
            bhadronFlavor_1 = abs(analysisTree.pfjet_flavour[indexLeadingBJet]);

        if (indexSubLeadingBJet>=0)
        {
            bpt_2 = analysisTree.pfjet_pt[indexSubLeadingBJet];
            beta_2 = analysisTree.pfjet_eta[indexSubLeadingBJet];
            bphi_2 = analysisTree.pfjet_phi[indexSubLeadingBJet];
        }
        if(!isData)
            bhadronFlavor_2 = abs(analysisTree.pfjet_flavour[indexSubLeadingJet]);

        jpt_1 = -9999;
        jeta_1 = -9999;
        jphi_1 = -9999;
            
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
            
        jpt_2 = -9999;
        jeta_2 = -9999;
        jphi_2 = -9999;
            
        if (indexSubLeadingJet>=0)
        {
            jpt_2 = analysisTree.pfjet_pt[indexSubLeadingJet];
            jeta_2 = analysisTree.pfjet_eta[indexSubLeadingJet];
            jphi_2 = analysisTree.pfjet_phi[indexSubLeadingJet];
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

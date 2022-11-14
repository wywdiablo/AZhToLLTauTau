#include "HttStylesNew.cc"
#include "CMS_lumi.C"
Double_t pol2(Double_t * x, Double_t * par) {

  double r = par[0]+par[1]*x[0]+par[2]*x[0]*x[0];//+par[3]*x[0]*x[0]*x[0];

  return r;

}


void FitFR(TString LepFlav = "Ele", //Options: Ele, Mu and Tau
           TString ID = "Fall17MVAv2WP80_noIso_Iso0p15",
           TString DMEtaBin = "barrel", //Options: barrel and endcap for Ele/Mu; DM0, DM1, DM10 and DM11 for Tau
	       TString era = "2018") {

  SetStyle();
  bool logX = false;
  bool logY = false;
  TString LepLatex = "";
  if(LepFlav == "Ele")
      LepLatex = "e";
  if(LepFlav == "Mu")
      LepLatex = "#mu";
  if(LepFlav == "Tau")
      LepLatex = "#tau";
 
  TFile * fileFR = new TFile("FakeRateFiles/Jet"+LepFlav+"FakeRate_"+ID+"_"+era+".root");
  TH1D * histFR = (TH1D*)fileFR->Get("PostfitFR_"+DMEtaBin);

  TH1D * FRH = (TH1D*)histFR->Clone("FRH");
  FRH->GetYaxis()->SetTitle("Jet to " + LepLatex + " Fake Rate");
  FRH->GetXaxis()->SetTitle(LepLatex +" P_{T} (GeV)");
  FRH->SetTitle("");
  InitData(FRH);
    
  float yUpper = FRH->GetMaximum();

  FRH->GetYaxis()->SetRangeUser(0.,2.5*yUpper);
  TCanvas * dummy1 = new TCanvas("dummy","",400,400);
    int ptLowThreshold = 10;
    int ptHighThreshold = 80;
    if((era == "UL2016postVFP") && (LepFlav == "Mu"))
    {
        ptLowThreshold = 10;
        ptHighThreshold = 40;
    }
    if(LepFlav == "Tau")
    {
        ptLowThreshold = 20;
        ptHighThreshold = 100;
    }
  //TF1 * fitFunc = new TF1("fitFunc",pol2,0.5,5.3,3);
    TF1 * fitFunc = new TF1("fitFunc",pol2,ptLowThreshold,ptHighThreshold,3);
  fitFunc->SetLineColor(2);
  fitFunc->SetParameter(0,1.);
  fitFunc->SetParameter(1,0.);
  fitFunc->SetParameter(2,0.);
  fitFunc->SetParameter(3,0.);
  
  FRH->Fit("fitFunc","Q");
  delete dummy1;
    
  TH1D * hint = new TH1D("ff_"+era,"",200,ptLowThreshold,ptHighThreshold);
  TH1D * hint_68 = new TH1D("ff_"+era+"_68","",200,ptLowThreshold,ptHighThreshold);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint,0.95);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_68,0.68);
  TH1D * hint_central = (TH1D*)hint->Clone("hint_central");
  for (int iB=1; iB<=200; ++iB)
    hint_central->SetBinError(iB,0.);
  InitHist(hint,"","",kCyan,1001);
  InitHist(hint_68,"","",kYellow,1001);
  hint->SetMarkerSize(0);
  hint->SetMarkerStyle(0);
  hint_68->SetMarkerSize(0);
  hint_68->SetMarkerStyle(0);
  hint_central->SetLineColor(2);
  hint_central->SetLineWidth(2);
    
    TFile * outputfile = new TFile("JetTo"+LepFlav+"FR_POL2Fit_"+ID+"_"+DMEtaBin+"_"+era+".root","recreate");
    TH1D * Uncert95 = (TH1D*)hint->Clone("POL2FitFR_Uncert95_"+DMEtaBin);
    TH1D * Uncert68 = (TH1D*)hint_68->Clone("POL2FitFR_Uncert68_"+DMEtaBin);
    TH1D * Central = (TH1D*)hint_central->Clone("POL2FitFR_Central_"+DMEtaBin);
    outputfile->Write();
    outputfile->Close();
  
  TCanvas * canv = MakeCanvas("canv1", "", 500, 500);

  FRH->Draw("e1");
  hint->Draw("e3same");
  hint_68->Draw("e3same");
  //  hint_central->Draw("hsame");
  FRH->Draw("e1same");

  TLegend * leg= new TLegend(0.20,0.85,0.50,0.90);
  SetLegendStyle(leg);
  leg->SetTextSize(0.03);
  leg->SetHeader(ID + " " + DMEtaBin);
  leg->Draw();

  if (logX) canv->SetLogx(true);
  if (logY) canv->SetLogy(true);
    if(era == "2016")
        lumi_13TeV = "2016, 35.9 fb^{-1}";
    if(era == "2017")
        lumi_13TeV = "2017, 41.5 fb^{-1}";
    if(era == "2018")
        lumi_13TeV = "2018, 59.7 fb^{-1}";
    if(era == "UL2016postVFP")
        lumi_13TeV = "UL2016postVFP, 16.8 fb^{-1}";
    if(era == "UL2016preVFP")
        lumi_13TeV = "UL2016preVFP, 19.5 fb^{-1}";
  writeExtraText = true;
  extraText = "Work in progress";
  CMS_lumi(canv,4,33); 

  canv->Draw("SAME");
  canv->RedrawAxis();
  canv->Modified();
  canv->Update();

  canv->Print("FakeRateFit_JetTo"+LepFlav+"_"+ID+"_"+DMEtaBin+"_"+era+".pdf");

}

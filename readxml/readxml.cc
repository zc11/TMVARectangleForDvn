#include "readxml.h"
#include "Tools.h"

void readxml(Float_t ptMin=5.5, Float_t ptMax=7., Float_t RAA=1.)
{
  void calRatio(Float_t* results, Bool_t verbose=false);
  ptmin = ptMin;
  ptmax = ptMax;
  raa = RAA;

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  //read weight file
  const char* filename = "../myTMVA/weights/TMVAClassification_CutsSA.weights.xml";
  void *doc = TMVA::gTools().xmlengine().ParseFile(filename,TMVA::gTools().xmlenginebuffersize());
  void* rootnode = TMVA::gTools().xmlengine().DocGetRootElement(doc); // node "MethodSetup"
  TString fullMethodName("");  
  TMVA::gTools().ReadAttr(rootnode, "Method", fullMethodName);

  cout<<endl;
  cout<<" ╒══════════════════════════════════════════════════╕"<<endl;
  cout<<" |               Cut Opt Configuration              |"<<endl;
  cout<<" ├────────────┬────────────────────────────┬────────┤"<<endl;
  cout<<" | "<<setiosflags(ios::left)<<setw(10)<<"Method"<<" | "<<setiosflags(ios::left)<<setw(26)<<fullMethodName<<" | "<<setiosflags(ios::left)<<setw(6)<<" "<<" |"<<endl;

  void *opts = TMVA::gTools().GetChild(rootnode,"Options");
  void* opt = TMVA::gTools().GetChild(opts,"Option");

  TString varProp("");
  while (opt)
    {
      TString optname("");
      TMVA::gTools().ReadAttr(opt, "name", optname);
      if (optname=="VarProp") varProp = TMVA::gTools().GetContent(opt);
      opt = TMVA::gTools().GetNextChild(opt);
    }

  TObjArray *marginclass = varProp.Tokenize(" ");
  std::vector<TString> margins;//avoid objarrays
  for(int i=0;i<marginclass->GetEntries();i++)
    {
      margins.push_back(((TObjString *)(marginclass->At(i)))->String());
    }
  void* variables = TMVA::gTools().GetChild(rootnode,"Variables");
  UInt_t nVar=0;
  std::vector<TString> varnames;
  TMVA::gTools().ReadAttr(variables, "NVar", nVar);

  void* var = TMVA::gTools().GetChild(variables,"Variable");
  for(unsigned int k=0;k<nVar;k++)
    {
      TString varname("");
      TMVA::gTools().ReadAttr(var, "Expression", varname);
      TString tem = Form("Variable%i",k);
      varval[k] = varname;
      cout<<" ├────────────┼────────────────────────────┼────────┤"<<endl;
      cout<<" | "<<setiosflags(ios::left)<<setw(10)<<tem<<" | "<<setiosflags(ios::left)<<setw(26)<<varname<<" | "<<setiosflags(ios::left)<<setw(6)<<margins[k]<<" |"<<endl;
      var = TMVA::gTools().GetNextChild(var);
      varnames.push_back(varname);
    }
  cout<<" ╞════════════╪════════════════════════════╪════════╡"<<endl;
    
  void* weight = TMVA::gTools().GetChild(rootnode,"Weights");
  void* eff = TMVA::gTools().GetChild(weight,"Bin");
  int n=0;
  while(eff)
    {
      TMVA::gTools().ReadAttr(eff, "effS", effS[n]);
      TMVA::gTools().ReadAttr(eff, "effB", effB[n]);
      void* cutsnode = TMVA::gTools().GetChild(eff,"Cuts");

      TString cut;
      for(ULong_t l=0;l<varnames.size();l++)
	{
	  Double_t min,max;
	  TMVA::gTools().ReadAttr(cutsnode, TString("cutMin_")+l, min);
	  TMVA::gTools().ReadAttr(cutsnode, TString("cutMax_")+l, max);
	  TString lessmax = "<"; lessmax+=max;
	  TString moremin = ">"; moremin+=min;
	  if(margins[l]=="FMin")
	    {
	      cut+=" && "+varnames[l]+lessmax;
	      cutval[l].push_back(max);
	    }
	  if(margins[l]=="FMax")
	    {
	      cut+=" && "+varnames[l]+moremin;
	      cutval[l].push_back(min);
	    }
	}
      cuts.push_back(cut);
      eff = TMVA::gTools().GetNextChild(eff);
      n++;
    }
  TMVA::gTools().xmlengine().FreeDoc(doc);

  Float_t wSignal=0;
  Float_t wBackground=0;
  Float_t* weights = new Float_t[2];
  //
  calRatio(weights);//weight signal and background
  //

  wSignal = weights[1];
  wBackground = weights[0];

  cout<<endl;
  cout<<"Looking for max significance ..."<<endl;

  Double_t max = wSignal*effS[1]/sqrt(wSignal*effS[1]+wBackground*effB[1]);
  int maxindex = 1;
  effS[0]=0;
  for(int i=1;i<100;i++)
    {
      effSig[i] = wSignal*effS[i]/sqrt(wSignal*effS[i]+wBackground*effB[i]);
      if(effSig[i]>max)
	{
	  max=effSig[i];
	  maxindex=i;
	}
    }
  cout<<endl;
  cout<<" ╒══════════════════════════════════════════════════╕"<<endl;
  cout<<" |                     Opt Result                   |"<<endl;
  cout<<" ├────────────┬────────────┬───────────────┬────────┤"<<endl;
  cout<<" | "<<setiosflags(ios::left)<<setw(10)<<"Sig eff"<<" | "<<setiosflags(ios::left)<<setw(10)<<effS[maxindex]<<" | "<<setiosflags(ios::left)<<setw(13)<<"S/sqrt(S+B)"<<" | "<<setiosflags(ios::left)<<setw(6)<<max<<" |"<<endl;
  cout<<" ├────────────┴────────────┴───┬───────────┴────────┤"<<endl;

  for(int m=0;m<nVar;m++)
    {
      if(m) cout<<" ├─────────────────────────────┼────────────────────┤"<<endl;
      cout<<" | "<<setiosflags(ios::left)<<setw(27)<<varval[m]<<" | "<<setiosflags(ios::left)<<setw(18)<<cutval[m].at(maxindex)<<" |"<<endl;
    }
  cout<<" ╘═════════════════════════════╧════════════════════╛"<<endl;
  cout<<endl;

  TH2F* hempty = new TH2F("hempty","",50,0,1.,10,0.,max*1.2);  
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetYaxis()->SetTitle("Signal efficiency");
  hempty->GetXaxis()->SetTitle("S/sqrt(S+B)");
  hempty->GetXaxis()->SetTitleOffset(0.9);
  hempty->GetYaxis()->SetTitleOffset(1.0);
  hempty->GetXaxis()->SetTitleSize(0.05);
  hempty->GetYaxis()->SetTitleSize(0.05);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.035);
  hempty->GetYaxis()->SetLabelSize(0.035);
  TLatex* texPar = new TLatex(0.18,0.93, "PbPb 2.76 TeV D^{0}");
  texPar->SetNDC();
  texPar->SetTextAlign(12);
  texPar->SetTextSize(0.04);
  texPar->SetTextFont(42);
  TLatex* texPtY = new TLatex(0.96,0.93, Form("|y|<1, %.1f<p_{T}<%.1f GeV/c",ptmin,ptmax));
  texPtY->SetNDC();
  texPtY->SetTextAlign(32);
  texPtY->SetTextSize(0.04);
  texPtY->SetTextFont(42);

  TGraph* gsig = new TGraph(100,effS,effSig);
  TCanvas* csig = new TCanvas("csig","",600,600);
  hempty->Draw();
  texPar->Draw();
  texPtY->Draw();
  gsig->Draw("same*");
  csig->SaveAs("plots/Significance.pdf");
}

void calRatio(Float_t* results, Bool_t verbose=false)
{
  TFile *inputS = TFile::Open(inputSname);
  TFile *inputB = TFile::Open(inputBname);

  TTree *signal = (TTree*)inputS->Get("recodmesontree");
  TTree *background = (TTree*)inputB->Get("recodmesontree");
  TTree *generated = (TTree*)inputS->Get("gendmesontree");

  TString sels = Form("%s&&dcandpt>%f&&dcandpt<%f",mycuts.Data(),ptmin,ptmax);
  TString selb = Form("%s&&dcandpt>%f&&dcandpt<%f",mycutb.Data(),ptmin,ptmax);
  TString selg = Form("%s&&dpt>%f&&dpt<%f",mycutg.Data(),ptmin,ptmax);

  TString ptstring = Form("(%.1f,%.1f)",ptmin,ptmax);
  cout<<" | "<<setiosflags(ios::left)<<setw(10)<<"Pt"<<" | "<<setiosflags(ios::left)<<setw(26)<<ptstring<<" | "<<setiosflags(ios::left)<<setw(6)<<" "<<" |"<<endl;
  cout<<" ├────────────┼────────────────────────────┼────────┤"<<endl;
  cout<<" | "<<setiosflags(ios::left)<<setw(10)<<"Raa"<<" | "<<setiosflags(ios::left)<<setw(26)<<raa<<" | "<<setiosflags(ios::left)<<setw(6)<<" "<<" |"<<endl;
  cout<<" ╘════════════╧════════════════════════════╧════════╛"<<endl;
  cout<<endl;

  //Get signal peak sigma
  TH1D* hmassS = new TH1D("hmassS",";D mass (GeV/c^{2});Signal Entries",50,1.7,2.1);
  signal->Project("hmassS","dcandmass",sels);
  hmassS->Sumw2();
  TCanvas* cmassS = new TCanvas("cmassS","",600,600);
  hmassS->Draw();
  TF1* fmass = new TF1("fmass","[0]*Gaus(x,[1],[2])");
  fmass->SetParLimits(2,0.005,0.05);
  Float_t setparam1 = 1.86;
  Float_t setparam2 = 0.01;
  fmass->SetParameter(1,setparam1);
  fmass->SetParameter(2,setparam2);
  if(verbose) hmassS->Fit("fmass","L","",1.7,2.1);
  else hmassS->Fit("fmass","L q","",1.7,2.1);
  cmassS->SaveAs("plots/Signal.pdf");
  Float_t sigma = fmass->GetParameter(2);

  //Background candidate number
  TH1D* hmassB = new TH1D("hmassB",";D mass (GeV/c^{2});Background Entries",50,0,10);
  background->Project("hmassB","dcandmass",selb);
  TCanvas* cmassB = new TCanvas("cmassB","",600,600);
  hmassB->Draw();
  cmassB->SaveAs("plots/Background.pdf");
  Int_t nentriesB = hmassB->Integral();

  //FONLL
  ifstream getdata("fonlls/fo_Dzero_pp_2p76_y1.dat");
  if(!getdata.is_open()) cout<<"Opening the file fails"<<endl;
  Float_t tem;
  Int_t nbin=0;
  while (!getdata.eof())
    {
      getdata>>pt[nbin]>>central[nbin]>>tem>>tem>>tem>>tem>>tem>>tem>>tem>>tem>>tem>>tem>>tem>>tem;
      if(pt[nbin]>=ptmin&&pt[nbin]<=ptmax) nbin++;
    }
  TH1D* hfonll = new TH1D("hfonll",";D p_{T} (GeV/c);FONLL differential xsection",nbin-1,pt);
  for(int i=0;i<nbin;i++)
    {
      hfonll->SetBinContent(i,central[i]);
    }
  TCanvas* cfonll = new TCanvas("cfonll","",600,600);
  hfonll->Draw();
  cfonll->SaveAs("plots/Fonll.pdf");

  TH1D* hrec = new TH1D("hrec",";D p_{T} (GeV/c);Signal reco entries",nbin-1,pt);
  TH1D* hgen = new TH1D("hgen",";D p_{T} (GeV/c);Generated entries",nbin-1,pt);
  TH1D* heff = new TH1D("heff",";D p_{T} (GeV/c);Efficiency",nbin-1,pt);
  signal->Project("hrec","dcandpt",mycuts);
  generated->Project("hgen","dpt",mycutg);
  heff->Divide(hrec,hgen,1.,1.,"B");
  TCanvas* ceff = new TCanvas("ceff","",600,600);
  heff->Draw();

  TH1D* htheoryreco = new TH1D("htheoryreco","",nbin-1,pt);
  htheoryreco->Multiply(heff,hfonll,1,1,"B");

  Double_t nevent = background->GetEntries();
  Double_t Taa = 5.65; //mb^-1
  Double_t BR = 0.0387;
  Double_t deltapt = 0.25;
  //central[i] - in pb/GeV/c

  Double_t yieldDzero = htheoryreco->Integral();
  yieldDzero*=(1.e-9)*BR*deltapt*Taa*nevent*raa;

  results[0] = nentriesB*Nsigma*sigma/0.05;//0.05: half of sideband width
  results[1] = yieldDzero;
  cout<<endl;
  cout<<" ╒══════════════════════════════════════════════════╕"<<endl;
  cout<<" |                   Weight Result                  |"<<endl;
  cout<<" ├────────────┬────────────┬────────────┬───────────┤"<<endl;
  cout<<" | "<<setiosflags(ios::left)<<setw(10)<<"Bkg #"<<" | "<<setiosflags(ios::left)<<setw(10)<<nentriesB<<" | "<<setiosflags(ios::left)<<setw(10)<<"Sig reg"<<" | "<<setiosflags(ios::left)<<setw(9)<<setprecision(3)<<sigma*Nsigma*2<<" |"<<endl;
  cout<<" ├────────────┼────────────┼────────────┼───────────┤"<<endl;
  cout<<" | "<<setiosflags(ios::left)<<setw(10)<<"SigWeight"<<" | "<<setiosflags(ios::left)<<setw(10)<<yieldDzero<<" | "<<setiosflags(ios::left)<<setw(10)<<"BkgWeight"<<" | "<<setiosflags(ios::left)<<setw(9)<<nentriesB*Nsigma*sigma/0.05<<" |"<<endl;
  cout<<" ╘════════════╧════════════╧════════════╧═══════════╛"<<endl;
}

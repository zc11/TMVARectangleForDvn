#include "readxml.h"
#include "Tools.h"

#include "TGraph.h"
#include <iostream>
#include <vector>
#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "TTree.h"

#include <iostream>
#include <fstream>

using namespace RooFit;
using namespace std;
//cut order overall: pTmax,pTmin,abs(eta),abs(daueta),VtxProb,3DPointingAngle,3DDecayLengthSignificance,pTD,NHitD
//cut order TMVA: VtxProb,3DPointingAngle,pTD1,pTD2

void readxml(Float_t ptMin=2.5, Float_t ptMax=3.0)
{
  ptmin = ptMin;
  ptmax = ptMax;

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
    TString ptstring = Form("(%.1f,%.1f)",ptmin,ptmax);
    cout<<" | "<<setiosflags(ios::left)<<setw(10)<<"Pt"<<" | "<<setiosflags(ios::left)<<setw(26)<<ptstring<<" | "<<setiosflags(ios::left)<<setw(6)<<" "<<" |"<<endl;
    cout<<" ╘════════════╧════════════════════════════╧════════╛"<<endl;
    cout<<endl;
    
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
    
    cout<<"finish reading cuts"<<endl;
    //construct histos with TMVA cuts
    TFile *inputB = TFile::Open(inputBname);
TFile *inputS = TFile::Open(inputSname);
    if(!inputB || !inputS)
    {
        cout<<"file not found"<<endl;
        return;
    }
    
    TTree *background = (TTree*) inputB->Get("demo/D0para");
    TTree *signal = (TTree*) inputS->Get("D0para/D0para");

    TH1D* Dmass[100];
    TH1D* DmassS = new TH1D("DmassS","DmassS",60,1.7,2.0);
    TH1D* DmassB = new TH1D("DmassB","DmassB",60,1.7,2.0);

    for(int icut=0;icut<100;icut++)
    {
        Dmass[icut] = new TH1D(Form("mass_cut%d",icut),Form("mass_cut%d",icut),60,1.7,2.0);
    }
    
    cout<<"finish histo construction"<<endl;

    int nMult_ass_goodS;
    float ptS;
    float etaS;
    float massS;
    float VtxProbS;
    int nhit1S;
    int nhit2S;
    float dlosS;
    float pt1S;
    float pt2S;
    float ptErr1S;
    float ptErr2S;
    float eta1S;
    float eta2S;
    float agl_absS;
    float dlS;
    float dlerrS;

    //signal->SetBranchAddress("Ntrkoffline",&nMult_ass_goodS);
    signal->SetBranchAddress("pT",&ptS);
    //signal->SetBranchAddress("eta",&etaS);
    signal->SetBranchAddress("mass",&massS);
    signal->SetBranchAddress("VtxProb",&VtxProbS);
    signal->SetBranchAddress("3DPointingAngle",&agl_absS);
    signal->SetBranchAddress("3DDecayLengthSignificance",&dlosS);
    //signal->SetBranchAddress("NHitD1",&nhit1S);
    //signal->SetBranchAddress("NHitD2",&nhit2S);
    //signal->SetBranchAddress("pTD1",&pt1S);
    //signal->SetBranchAddress("pTD2",&pt2S);
    //signal->SetBranchAddress("pTerrD1",&ptErr1S);
    //signal->SetBranchAddress("pTerrD2",&ptErr2S);
    //signal->SetBranchAddress("EtaD1",&eta1S);
    //signal->SetBranchAddress("EtaD2",&eta2S);
    //signal->SetBranchAddress("3DDecayLength",&dlS);
    //signal->SetBranchAddress("3DDecayLengthError", &dlerrS);

    Int_t nentriesS = signal->GetEntries();
    for (Int_t i=0;i<nentriesS;i++)
    {
        signal->GetEntry(i);
        
        //first implement non-tuning cut; must be sychronized with the setting in TMVA tuning (mycutb)
        //if(nMult_ass_goodS<185 || nMult_ass_goodS>=250) continue;
        if(ptS>3.0 || ptS<2.5) continue;
        if(dlosS<=3.5) continue;
        if(VtxProbS<0.05) continue;
        if(agl_absS>0.20) continue;
        
            DmassS->Fill(massS);
        
    }
    
    int nMult_ass_good;
    float pt;
    float eta;
    float mass;
    float VtxProb;
    int nhit1;
    int nhit2;
    float dlos;
    float pt1;
    float pt2;
    float ptErr1;
    float ptErr2;
    float eta1;
    float eta2;
    float agl_abs;
    float dl;
    float dlerr;
    
    //background->SetBranchAddress("Ntrkoffline",&nMult_ass_good);
    background->SetBranchAddress("pT",&pt);
    //background->SetBranchAddress("eta",&eta);
    background->SetBranchAddress("mass",&mass);
    background->SetBranchAddress("VtxProb",&VtxProb);
    background->SetBranchAddress("3DPointingAngle",&agl_abs);
    background->SetBranchAddress("3DDecayLengthSignificance",&dlos);
    //background->SetBranchAddress("NHitD1",&nhit1);
    //background->SetBranchAddress("NHitD2",&nhit2);
    //background->SetBranchAddress("pTD1",&pt1);
    //background->SetBranchAddress("pTD2",&pt2);
    //background->SetBranchAddress("pTerrD1",&ptErr1);
    //background->SetBranchAddress("pTerrD2",&ptErr2);
    //background->SetBranchAddress("EtaD1",&eta1);
    //background->SetBranchAddress("EtaD2",&eta2);
    //background->SetBranchAddress("3DDecayLength",&dl);
    //background->SetBranchAddress("3DDecayLengthError",&dlerr);


    
    Int_t nentries = background->GetEntries();
    for (Int_t i=0;i<nentries;i++)
    {
        background->GetEntry(i);
        
        //first implement non-tuning cut; must be sychronized with the setting in TMVA tuning (mycutb)
        if(pt>3.0 || pt<2.5) continue;
        if(dlos<=3.5) continue;
        if(VtxProb<0.05) continue;
        if(agl_abs>0.20) continue;

        DmassB->Fill(mass);
        
        for(int icut=0;icut<100;icut++)
        {
            if(VtxProb<cutval[0].at(icut)) continue;
            if(agl_abs>cutval[1].at(icut)) continue;
            if(dlos<cutval[2].at(icut)) continue;
          //  if(dlerr>cutval[3].at(icut)) continue;
          //  if(pt1<cutval[0].at(icut)) continue;
          //  if(pt2<cutval[1].at(icut)) continue;
            Dmass[icut]->Fill(mass);
        }
    }

    cout<<"finish fill histos"<<endl;
    
    //fit the histos to get sigsig
    double sigsig[100];
    double sigeff[100];
    
    for(int icut=0;icut<100;icut++)      
    {
        RooMsgService::instance().setStreamStatus(0,kFALSE);
        RooMsgService::instance().setStreamStatus(1,kFALSE);

        TCanvas *c = new TCanvas("c","c",1600,1600);
        c->cd();
        TH1D* h;
        gDirectory->GetObject(Form("mass_cut%d",icut),h);
        
        RooRealVar x("x","mass",1.70, 2.00);
        RooDataHist data("data","dataset", x, h);
        RooPlot* xframe = x.frame(60);
        data.plotOn(xframe,Name("data"));
        RooRealVar mean("mean","mean",1.86, 1.82, 1.90);
        RooRealVar sigma1("sigma1","sigma1",0.015,0.005,0.03);
        RooRealVar sigma2("sigma2","sigma2",0.015,0.005,0.03);
        RooRealVar sig1("sig1","signal1",10,0,10000000);
        RooRealVar sig2("sig2","signal2",10,0,10000000);
        RooRealVar a("a","a",0,-100000,100000);
        RooRealVar b("b","b",0,-100000,100000);
        RooRealVar cp("cp","cp",0,-100000,100000);
        RooRealVar d("d","d",0,-100000,100000);
        RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
        RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
        RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp,d));
        RooRealVar polysig("polysig","polysig",10,0,10000000);
        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));
        
        RooFitResult *r = sum.fitTo(data,Save(),Minos(kTRUE),PrintLevel(-1));
        r = sum.fitTo(data,Save(),Minos(kTRUE),PrintLevel(-1));
        r = sum.fitTo(data,Save(),Minos(kTRUE),PrintLevel(-1));
        r = sum.fitTo(data,Save(),Minos(kTRUE),PrintLevel(-1));
        r = sum.fitTo(data,Save(),Minos(kTRUE),PrintLevel(-1));
        
        sum.plotOn(xframe,Name("sum"));
        sum.plotOn(xframe,Components(poly),LineStyle(kDashed));

        // r->Print();
        
        double chi2 = xframe->chiSquare("sum","data");
        double meanf = mean.getVal();
        double meanfe = mean.getError();
        double sigmaf1 = sigma1.getVal();
        double sigmaf2 = sigma2.getVal();
        double bkgf = polysig.getVal();
        double sigf1 = sig1.getVal();
        double sigf2 = sig2.getVal();
        double sigwf1 = sigf1/(sigf1 + sigf2);
        double sigwf2 = sigf2/(sigf1 + sigf2);
        double c1 = a.getVal();
        double c2 = b.getVal();
        
        double sigmaf = sqrt(sigmaf1 * sigmaf1 * sigwf1 + sigmaf2 * sigmaf2 * sigwf2);
        double massmin = meanf - 2.0*sigmaf;
        double massmax = meanf + 2.0*sigmaf;

        //TLine* l1 = new TLine(massmin,0,massmin,500);
        //TLine* l2 = new TLine(massmax,0,massmax,500);
        //l1->Draw("same");
        //l2->Draw("same");

        int nmin = h->GetXaxis()->FindBin(massmin);
        int nmax = h->GetXaxis()->FindBin(massmax);
        int anmin = h->GetXaxis()->FindBin(1.70);
        int anmax = h->GetXaxis()->FindBin(2.00);
        
        double awyh1 = h->Integral(anmin,nmin);
        double awyh2 = h->Integral(nmax,anmax);
        double awyh = awyh1 + awyh2;
        double totyh = h->Integral(nmin,nmax);
        
        x.setRange("cut",massmin,massmax);
        RooAbsReal* ibkg = poly.createIntegral(x,NormSet(x),Range("cut"));
        RooAbsReal* isig1 = gaus1.createIntegral(x,NormSet(x),Range("cut"));
        RooAbsReal* isig2 = gaus2.createIntegral(x,NormSet(x),Range("cut"));
        double ibkgf = ibkg->getVal();
        double bkgfe = polysig.getError();
        double isig1f = isig1->getVal();
        double isig2f = isig2->getVal();
        
        double bkgy = ibkgf*bkgf;
        double bkgye = ibkgf*bkgfe;
        double sigy1 = isig1f*sigf1;
        double sigy2 = isig2f*sigf2;
        double sigy = sigy1 + sigy2;
        double toty = bkgy + sigy;
        
        double abkgy = (1-ibkgf)*bkgf;
        double asigy1 = (1-isig1f)*sigf1;
        double asigy2 = (1-isig2f)*sigf2;
        double asigy = asigy1 + asigy2;
        double awy = abkgy + asigy;
        
        double sigfrac = sigy/toty;
        double bkgfrac = bkgy/toty;
        
        double sigyh = totyh - bkgy;
        double sigfrach = sigy/totyh;
        double bkgfrach = bkgy/totyh;
        
        double signif = sigy/sqrt(totyh);
        sigsig[icut] = signif;
        sigeff[icut] = icut;
        xframe->SetXTitle("M(K#pi)[GeV]");
        xframe->GetXaxis()->CenterTitle();
        xframe->GetXaxis()->SetTitleOffset(1.2);
        xframe->GetXaxis()->SetLabelSize(0.02);
        xframe->SetYTitle("Events(per 0.005GeV)");
        xframe->GetYaxis()->CenterTitle();
        xframe->GetYaxis()->SetTitleOffset(1.2);
        xframe->GetYaxis()->SetLabelSize(0.02);
        xframe->SetTitle("");
        xframe->Draw();
        h->Draw("PEsame");

        TLatex latex;
        latex.SetTextAlign(12);
        latex.SetTextSize(0.03);
        latex.SetNDC();

        latex.DrawLatex(0.70,0.85,Form("#frac{s}{#sqrt{s+b}} = %.3f", signif));
        latex.DrawLatex(0.70,0.80,Form("#frac{s}{s+b} = %.3f", sigfrach));
        latex.DrawLatex(0.70,0.75,Form("s+b = %.1f", totyh));
        latex.DrawLatex(0.70,0.70,Form("mean = %.3f", meanf));
        latex.DrawLatex(0.70,0.65,Form("#sigma = %.4f", sigmaf));

        c->Print(Form("./plots/fits/cut%d.pdf", icut));

    }

    for(int icut=0;icut<100;icut++)
    {
        cout<<"icut: "<<icut<<" ,sigsig: "<<sigsig[icut];
      //  cout<<"VtxProb > "<<cutval[0].at(icut)<<", angle < "<<cutval[1].at(icut);
        cout<<", VtxProb > "<<cutval[0].at(icut)<<", angle < "<<cutval[1].at(icut)<<", dlos > "<<cutval[2].at(icut)<<endl;
      //  cout<<", dl > "<<cutval[2].at(icut)<<", dlerr <"<<cutval[3].at(icut)<<endl;       
    }

    TCanvas* cg = new TCanvas("cg","cg",800,800);
    cg->cd();
    TGraph* g = new TGraph (100 ,sigeff, sigsig);
    g->Draw("A");
    cg->Print("./plots/sigeff.pdf");

    TCanvas *cS = new TCanvas("cS","cS",1600,1600);
    cS->cd();
    DmassS->Draw("E");
    cS->Print("./plots/signal.pdf");

    TCanvas *cB = new TCanvas("cB","cB",1600,1600);
    cB->cd();
    DmassB->Draw("E");
    cB->Print("./plots/background.pdf");    


}





#include "uti.h"

const int NmaxFonll = 401; //fonll data points number
Float_t central[NmaxFonll],pt[NmaxFonll];
const int NEff = 100;
Double_t effS[NEff], effB[NEff],effSig[NEff],effBac[NEff];
const int NmaxVar = 5;
std::vector<TString> cuts;
std::vector<Double_t> cutval[NmaxVar];
TString varval[NmaxVar];

Int_t Nsigma = 2;
Float_t ptmin,ptmax,raa;

TString inputSname = "/uscms_data/d3/azhang/MYWORK/CMSSW_8_0_24/src/newall/D0tree_Pythia8_pp502_Pthat0_prompt_loosecut_matched_deltaR0p5.root";
TString inputBname = "/uscms_data/d3/azhang/MYWORK/CMSSW_8_0_24/src/newall/a.root";
//TString mycuts = "pT<3.0&&pT>2.5&&VtxProb>0.05&&3DPointingAngle<0.30&&TMath::Abs(3DDecayLength/3DDecayLengthError)>1.5";
//TString mycutb = "pT<3.0&&pT>2.5&&VtxProb>0.05&&3DPointingAngle<0.30&&TMath::Abs(3DDecayLength/3DDecayLengthError)>1.5";
//TString mycutg = "eta>-1.&&eta<1.";

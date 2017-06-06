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

TString inputSname = "/data/aozhang/TMVA/D0tree_Pythia8_pp502_Pthat0_5_10_15_prompt_loosecut_matched_deltaR0p5.root";
TString inputBname = "/data/aozhang/TMVA/a.root";
TString mycuts = "pT<2&&pT>1&&eta<1.0&&eta>-1.0&&EtaD1<1.5&&EtaD1>-1.5&&EtaD2<1.5&&EtaD2>-1.5&&VtxProb>0.05&&3DPointingAngle<0.30&&pTD1>0.7&&pTD2>0.7&&NHitD1>=11&&NHitD2>=11";
TString mycutb = "Ntrkoffline>=185&&Ntrkoffline<250&&pT<2&&pT>1&&eta<1.0&&eta>-1.0&&EtaD1<1.5&&EtaD1>-1.5&&EtaD2<1.5&&EtaD2>-1.5&&VtxProb>0.05&&3DPointingAngle<0.30&&pTD1>0.7&&pTD2>0.7&&NHitD1>=11&&NHitD2>=11";
TString mycutg = "eta>-1.&&eta<1.";

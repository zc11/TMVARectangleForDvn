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

TString inputSname = "/data/zhchen/D0tree_Pythia8_pp502_Pthat0_5_10_15_prompt_loosecut_matched_deltaR0p5.root";
TString inputBname = "/data/zhchen/D0tree_pPb_8160_N185_250_defaultcut_testsample.root";
TString mycuts = "pT<2&&pT>1&&eta<1&&eta>-1&&EtaD1<1.5&&EtaD1>-1.5&&EtaD2<1.5&&EtaD2>-1.5&&VtxProb>0.01&&3DPointingAngle<1.5&&3DDecayLengthSignificance>0&&pTD1>0.5&&pTD2>0.5&&NHitD1>0&&NHitD2>0";
TString mycutb = "Ntrkoffline>=185&&Ntrkoffline<250&&pT<2&&pT>1&&eta<1&&eta>-1&&EtaD1<1.5&&EtaD1>-1.5&&EtaD2<1.5&&EtaD2>-1.5&&VtxProb>0.01&&3DPointingAngle<1.5&&3DDecayLengthSignificance>0&&pTD1>0.5&&pTD2>0.5&&NHitD1>0&&NHitD2>0";
TString mycutg = "eta>-1.&&eta<1.";

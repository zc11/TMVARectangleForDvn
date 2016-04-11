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

TString inputSname = "/data/wangj/TutorialsSamples/Dmesonana_hiforest_official_PbPbD0tokaonpion_Pt0153050_tkpt1p0eta1p1_2760GeV_0613_genmatched.root";
TString inputBname = "/data/wangj/TutorialsSamples/Dmesonana_Rereco_MBtrig_d0pt0_y1p2_tk1p0_eta1p1_d2p0_alpha0p2_tight_0710_6lumi.root";
TString mycuts = "dcandy>-1.&&dcandy<1.&&dcanddau1pt>1.0&&dcanddau2pt>1.0&&(matchedtogen&&nongendoublecounted)&&dcandffls3d>2.0&&TMath::ACos(dcandcosalpha)<0.12&&dcandfprob>0.05";
TString mycutb = "MinBias&&dcandy>-1.&&dcandy<1.&&dcanddau1pt>1.0&&dcanddau2pt>1.0&&(TMath::Abs(dcandmass-1.865)>0.10&&TMath::Abs(dcandmass-1.865)<0.15)&&dcandffls3d>2.0&&TMath::ACos(dcandcosalpha)<0.12&&dcandfprob>0.05";
TString mycutg = "dy>-1.&&dy<1.";

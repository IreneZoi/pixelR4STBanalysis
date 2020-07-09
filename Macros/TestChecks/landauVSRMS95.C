#include <stdio.h>
#include "TCanvas.h"
#include "TF1.h"
#include <TTreeFormula.h>
#include <TMath.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TVectorT.h>
#include <TMap.h>
#include <TAxis.h>
#include <sstream>

#include "fileHandler.h"
#define Cuts 5
#define conversionFactor 1.
#define comparisons 1
#define DreiMasterPlanes 3
bool print=true;
using namespace std;


void landauVSRMS95(TString function = "RMS95"){

  TString detectorA="146";
  TString detectorB="148";
  TString detectorC="163";
  TString labelA="FDB150Y_2_R4S100x25-Y2_2, thr 22 ADC";
  TString labelC="FDB150Y_2_R4S100x25-Y6_1, thr 22 ADC";
  TString labelB="Pstop_RD53Apads_FDB, thr 22 ADC";
  
  Int_t run=2743;        //{2735,2743,2758};
  TString Run;
  Run.Form("%d",run);
  if(print) cout << "Run: " << run << endl;
  
  //  TString Angle[comparisons] = {"0","8.75","27.5"};
  TString info = "beam energy 5.6 GeV, angle 8.75 deg";
  
  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";
  TString label = "1clusterABC";
  TString ss_dphcut = "22";


  TString  inputfile = inputDir+"drei-r"+Run+"_irene_dphcutB"+ss_dphcut+"_"+label+".root";
  cout << inputfile << endl;
  TFile * file = new TFile(inputfile);
  TTree * tree = (TTree*)file->Get("charge_res");
  tree->Print();
    
  TH1I * hdx3tree = new TH1I("hdx3tree", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  tree->Draw("dx3tree>>hdx3tree","","goff");
  hdx3tree = (TH1I*)gDirectory->Get("hdx3tree");
  
  double tolerance = (hdx3tree->GetBinContent(0)+hdx3tree->GetBinContent(hdx3tree->GetNbinsX()+1))/hdx3tree->GetEntries(); // since we are not working with continuous quantities but with binned hists, the difference between the integ ral and the 95% it may not be zero, so we ask it to be lower than 15% (before I was using 1.1% but since I changed to use the full range of hists, the overflow bin plays a too bigger role when stats is low and the 1.1% is not reached.

  cout << " RMS method charge based with tolerance " << tolerance << endl;
  double maximum = hdx3tree->GetMaximum();
  int maxbin = hdx3tree->GetMaximumBin();
  cout <<  " maxbin " << maxbin << " at " << hdx3tree->GetBinCenter(maxbin)<<endl;
  cout << "intital sigma = " << hdx3tree->GetRMS() * 1000 << " sigmaerr = " << hdx3tree->GetRMSError() * 1000 << endl;
  double Integral = hdx3tree->Integral(0,hdx3tree->GetNbinsX()+1);
  cout << " integral " << Integral << " entries " << hdx3tree->GetEntries() << endl;
  double integral95 = 0.95*Integral;
  cout << " integral95 " << integral95 << endl;
  int i = 0;
  
  double low = 0;
  double high = 0;
  for(int i =0; i<hdx3tree->GetNbinsX()/2; i++)  {
    low = maxbin-i;
    high = maxbin+i;
    
    Integral = hdx3tree->Integral(low,high);
    cout << " integral " << Integral << " low " << low << " high " << high << endl;
    cout << " while "<< i << " fabs(integral-integral95)/integral95 " << fabs(Integral-integral95)/integral95 << endl;
    
    //      if(fabs(Integral-integral95)/integral95 < 0.011 || Integral>integral95)//integral>integral95)
    if(fabs(Integral-integral95)/integral95 < tolerance || Integral>integral95)//integral>integral95)
      break;
    
    
  }
  cout << "final integral " << Integral << " low " << low <<" high " << high << endl;
  //  hdx3tree->GetXaxis()->SetRange(low,high);

  TString ss_low,ss_high;
  ss_low.Form("%f",hdx3tree->GetBinCenter(low));
  ss_high.Form("%f",hdx3tree->GetBinCenter(high));
  
  cout << ss_low << " " << ss_high << endl;

  TH1I * hclphB95 = new TH1I( "clphB", "B isolated cluster charge;cluster charge [ADC];B isolatewd clusters", 200, 0., 1000. );
  TH1I * hclphB5 = new TH1I( "clphB", "B isolated cluster charge;cluster charge [ADC];B isolatewd clusters", 200, 0., 1000. );

  tree->Draw("clphBiiitree>>hclphB95","dx3tree>"+ss_low+"&&dx3tree<"+ss_high,"goff");
  hclphB95 = (TH1I*)gDirectory->Get("hclphB95");

  cout  << "95 " << hclphB95->GetEntries() << " perc " << hclphB95->GetEntries()/hdx3tree->Integral(0,hdx3tree->GetNbinsX()+1) <<endl;

  tree->Draw("clphBiiitree>>hclphB5","dx3tree<"+ss_low+"||dx3tree>"+ss_high,"goff");
  hclphB5 = (TH1I*)gDirectory->Get("hclphB5");

  cout  << "5 " << hclphB5->GetEntries() <<" perc " << hclphB5->GetEntries()/hdx3tree->Integral(0,hdx3tree->GetNbinsX()+1) << endl;

 
  TCanvas *c2 = new TCanvas("c2", "resolution vs dphcut", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  hclphB95->SetTitle(" ");
  hclphB95->GetYaxis()->SetTitle("Entries");
  hclphB95->GetXaxis()->SetTitle("Cluster charge [ADC]");
  hclphB95->SetLineColor(kRed);
  hclphB95->SetMarkerStyle(20);
  hclphB95->SetLineWidth(2);
  hclphB95->SetLineStyle(2);
  hclphB95->GetXaxis()->SetRangeUser(0.,1000.);
  hclphB95->GetYaxis()->SetRangeUser(0.1,50000.);
  hclphB95->Draw("hist");
  hclphB5->Draw("histsame");
  c2->SetLogy();  
  TLegend* leg2 = new TLegend(0.6,0.6,0.8,0.8);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.04);
  leg2->AddEntry(hclphB95,"95%" , "l");
  leg2->AddEntry(hclphB5,"5%" ,"l");
  leg2->Draw();

  TString name = outputDir+"Residual95vs5_vs_clusterchergeABC_"+function+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+"_run"+Run+"_"+label;
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".pdf");
  
}//resolution 


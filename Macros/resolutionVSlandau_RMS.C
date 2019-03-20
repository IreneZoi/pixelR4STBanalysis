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
#define Cuts 9
#define conversionFactor 1.
#define comparisons 1
bool print=true;
using namespace std;


void resolutionVSlandau_RMS(TString func = "RMS")
{

  TString detectorA="146";
  TString detectorB="148";
  TString detectorC="163";
  TString labelA="FDB150Y_2_R4S100x25-Y2_2, thr 12 ADC";
  
  TString labelC="FDB150Y_2_R4S100x25-Y6_1, thr 12 ADC";
  
  
  TString labelB="Pstop_RD53Apads_FDB, thr 12 ADC";
  
  Int_t run=2743;
  TString Run[comparisons];
  Run[0].Form("%d",run);
  if(print) cout << "Run: " << run << endl;
  TString info = "beam energy 5.6 GeV, angle 8.75 deg";
  
  TString Landau[Cuts];
  Double_t Resolution[Cuts];
  Double_t ResolutionError[Cuts];
  
  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";
  TString inputfile;
  TFile * file[Cuts];
  TH1F * h_res[Cuts];

  TString Path;
  TString hist[Cuts];
  ostringstream strs[Cuts];
  hist[0] = "dx3_clchargeB2e";
  hist[1] = "dx3_clchargeB2e5e";
  hist[2] = "dx3_clchargeB5e8e";
  hist[3] = "dx3_clchargeB8e10e";
  hist[4] = "dx3_clchargeB10e13e";
  hist[5] = "dx3_clchargeB13e16e";
  hist[6] = "dx3_clchargeB16e22e";
  hist[7] = "dx3_clchargeB22e40e";
  hist[8] = "dx3_clchargeB40e";
																																	    
  Landau[0] = "0_2";
  Landau[1] = "2_5";
  Landau[2] = "5_8";
  Landau[3] = "8_10";
  Landau[4] = "10_13";
  Landau[5] = "13_16";
  Landau[6] = "16_22";
  Landau[7] = "22_40";
  Landau[8] = "40_inf";

  TString dirInFile[comparisons] = {"straightTracksY_isoAandCandB_straightTracksX"};
  int dph[comparisons] = {12};
  const Int_t NBINS = 9;
  Double_t edges[NBINS + 1] = {0.0, 2., 5., 8., 10., 13., 16., 22., 40., 60.};
  // Bin 1 corresponds to range [0.0, 0.2]
  // Bin 2 corresponds to range [0.2, 0.3] etc...

  TH1* h_res_vs_charge = new TH1D("h_res_vs_charge","Resolution vs cluster charge",NBINS,edges);
  MapTH1 m_hists;
  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  
  if(print) cout << "Getting files "  << endl;
  for(int i=0; i<Cuts; i++)
    {
      GetHists(&m_hists,1, 1, comparisons,  true,Run, dph, dirInFile, hist[i], h_res[i]);

      auto it2 = m_hists.find(std::make_pair(Run[0]+"_12_"+dirInFile[0],hist[i]));
      if(it2  != m_hists.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	  
	  FitTH1(it2->second, &(Resolution[i]), &(ResolutionError[i]), Landau[i], detectorA, detectorB, detectorC,func );
	  
	}
      
      if(print) cout << "bin " << i<< ": " << Landau[i] << " [ke] -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      
      h_res_vs_charge->SetBinContent(i+1,Resolution[i]); 
      h_res_vs_charge->SetBinError(i+1,ResolutionError[i]); 
    }
  
  
  if(print) cout << "Plotting Res vs dhcut"  << endl;
 
  TCanvas *c2 = new TCanvas("c2", "resolution vs dphcut", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  h_res_vs_charge->SetTitle(" ");
  h_res_vs_charge->GetYaxis()->SetTitle("RMS 95% [#mum]");
  h_res_vs_charge->GetXaxis()->SetTitle("Cluster charge [ke]");
  h_res_vs_charge->SetMarkerSize(2.5);
  h_res_vs_charge->SetLineColor(2);
  h_res_vs_charge->SetMarkerColor(2);
  h_res_vs_charge->SetMarkerStyle(20);
  h_res_vs_charge->SetLineWidth(2);
  //  h_res_vs_charge->SetLineStyle(2);
  h_res_vs_charge->GetXaxis()->SetRangeUser(0.,60.);
  h_res_vs_charge->GetYaxis()->SetRangeUser(0.,60.);
  //  h_res_vs_charge->Draw("L");
  h_res_vs_charge->Draw("E1X0");

  TLegend* leg2 = new TLegend(0.15,0.6,0.5,0.8);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(h_res_vs_charge,"c"+detectorB+" "+labelB , "lp");
  leg2->AddEntry(h_res_vs_charge,"A: c"+detectorA+" "+labelA ,"");
  leg2->AddEntry(h_res_vs_charge,"C: c"+detectorC+" "+labelC ,"");
  leg2->AddEntry(h_res_vs_charge,info ,"");
  leg2->Draw();
  

  TString name = outputDir+"Res_vs_clustercherge_"+func+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+"_run"+Run[0];
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".pdf");

}//resolution 


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
#include "TH1.h"
#include "TArrayD.h"
#include <sstream>
#include "fileHandler.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"


bool print=true;
using namespace std;
#define runNONirr "2743"
#define runPirr2 "2801"
#define runNirr4 "3839"
#define irradiations 3

void TDR();
void TDR2(TCanvas * c_all, int period=0, int pos= 11);
Double_t ScaleX(Double_t x);
void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t));

void controlPlots_forPaper(TString name = "preliminary_TreeCorr")
{

  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  TString irr[irradiations];
  irr[0] = "Non-irradiated, 120 V"; // "no irr, 5.6 GeV";
  irr[1] = "#phi_{eq} = 2.1 #times 10^{15} cm^{-2}, proton, 800 V";
  irr[2] = "#phi_{eq} = 3.6 #times 10^{15} cm^{-2}, neutron, 800 V";


  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  TH1F * h_res;  
  TString Label[1] = {"straightTracksY_isoAandCandB_straightTracksX"};
  TString Hist = "clphB";
  bool dphcut = true;
  //TString pitch[1] = {"25"}; //,"25","25"};
  TString pitch = "25"; //,"25","25"};

  // Non-irr 
  MapTH1 landau_map;
  TString Run[1] = {runNONirr};
  int i_dphcut[1] = {12}; 
  Double_t d_dphcut[1] = {12};
  TString ss_dphcut[1] = {"12"};
  TString inputfileNONirr= "closest_A13C14_bestnonirr";
  GetHists(&landau_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res,true,inputfileNONirr);
  
  //proton

  Run[0] = runPirr2;
  i_dphcut[0] = 15; 
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  TString  inputfilePirr2="closest_A12C15_bestproton";
  GetHists(&landau_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res,true,inputfilePirr2);
    
  

  //neutron

  Run[0] = runNirr4;
  i_dphcut[0] = 15; 
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  TString inputfileNirr4="closest_A12C13";
  GetHists(&landau_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res,true,inputfileNirr4);
    
  TH1F  h_landau[irradiations];
  TH1F  h_landau_orig[irradiations];
  
  float integral[irradiations];
  float integral90[irradiations];
  float bin90[irradiations];
  TString ss_dphcuts[irradiations] = {"12","15","15"};
  TString runs[irradiations] = {runNONirr,runPirr2,runNirr4};
  for(int i=0; i<irradiations; i++)
    {
      auto it2 = landau_map.find(std::make_pair(runs[i]+"_"+ss_dphcuts[i]+"_"+Label[0],Hist));
      if(it2  != landau_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	  h_landau_orig[i] = *it2->second;
          h_landau[i] = *it2->second;
	  integral[i] = h_landau[i].Integral(0, h_landau[i].GetNbinsX()+1);
	  cout << "integral " << integral[i] << endl;
	  
      	  h_landau_orig[i].Scale(1./integral[i]);
	  integral90[i]= 0.9*h_landau[i].Integral();

          int j = 0;
	  while(integral[i]>integral90[i])
	    {
	      cout << " while "<< j << endl;
	      integral[i] =  h_landau[i].Integral(0, h_landau[i].GetNbinsX()-j);
	      bin90[i] =  h_landau[i].GetBinCenter( h_landau[i].GetNbinsX()-j);
	      cout << " integral " << integral[i] << " high " << bin90[i] << endl;
	      j++;
	    }
	  cout << " integral90 " << integral90[i] << endl;
	}
      else cout << "map key " << it2->first.first << " not found " << endl;
    }

  /////// plots!

  
  TCanvas *cFDB2 = new TCanvas("cFDB2", "FDB resolution", 600, 600);
  cFDB2->SetLeftMargin(0.15);  
  cFDB2->SetRightMargin(1.5);  
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(15);
  cout << "starting making hists pretty" << endl;
  h_landau_orig[0].SetTitle(" ");
  h_landau_orig[0].GetYaxis()->SetTitle("Normalized number of events");
  h_landau_orig[0].GetXaxis()->SetTitle("Cluster charge [ADC]");
  //  h_landau_orig[0].SetMarkerSize(2.5);
  h_landau_orig[0].SetLineColor(kBlack);
  h_landau_orig[0].SetLineStyle(1);
  h_landau_orig[0].SetLineWidth(2);
  h_landau_orig[0].Rebin(2);
  h_landau_orig[0].GetXaxis()->SetLimits(0,1000.);
  h_landau_orig[0].GetXaxis()->SetMaxDigits(3); //SetNoExponent(true);
  h_landau_orig[0].GetYaxis()->SetRangeUser(0.0001,1.);
  h_landau_orig[0].GetXaxis()->SetRange(0.,h_landau_orig[0].GetNbinsX()+1);
  h_landau_orig[0].Draw("histe");
  //  TGaxis::SetMaxDigits(3);
  //  h_landau_orig[1].SetMarkerSize(2.5);
  h_landau_orig[1].GetXaxis()->SetLimits(0,1000.);
  //h_landau_orig[1].GetYaxis()->SetRangeUser(0.,0.08);
  h_landau_orig[1].Rebin(2);
  h_landau_orig[1].GetYaxis()->SetRangeUser(0.0001,1.);
  h_landau_orig[1].SetLineColor(kGreen+1);
  h_landau_orig[1].SetLineStyle(2);
  h_landau_orig[1].SetLineWidth(2);
  h_landau_orig[1].GetXaxis()->SetRange(0.,h_landau_orig[1].GetNbinsX()+1);
  //  h_landau_orig[1].GetXaxis()->SetNoExponent(true);
  h_landau_orig[1].Draw("histesame");
 				      
  // h_landau_orig[2].SetMarkerSize(2.5);
  h_landau_orig[2].GetXaxis()->SetLimits(0,1000.);
  h_landau_orig[2].GetYaxis()->SetRangeUser(0.0001,1.);
  h_landau_orig[2].GetXaxis()->SetRange(0.,h_landau_orig[2].GetNbinsX()+1);
  //h_landau_orig[2].GetYaxis()->SetRangeUser(0.,0.08);
  h_landau_orig[2].SetLineColor(kMagenta);
  h_landau_orig[2].Rebin(2);
  h_landau_orig[2].SetLineWidth(2);
  h_landau_orig[2].SetLineStyle(3);
  //  h_landau_orig[2].GetXaxis()->SetNoExponent(true);
  h_landau_orig[2].Draw("histesame");

  TLegend* legFDB2 = new TLegend(0.4,0.7,0.8,0.85);
  legFDB2->SetLineColor(0);
  //legFDB2->SetLegendFont(43);
  for(int i =0; i < irradiations; i++)
    legFDB2->AddEntry(&(h_landau_orig[i]),irr[i],"le");

  legFDB2->Draw();

  //landau90[i]
  int color[irradiations] = {1,417,616};
  TLine *  line2[irradiations];
  TArrow *ar[irradiations];
  //  float height[irradiations] = {;
  
  for(int i =0; i < irradiations; i++){
    //    cout << landau90[i] << endl;
    line2[i]    = new TLine( bin90[i],0.,bin90[i],0.015);
    line2[i]->SetLineColor(color[i]);
    line2[i]->SetLineWidth(2);
    line2[i]->SetLineStyle(1);
    line2[i]->Draw("same");
    ar[i]  = new TArrow(bin90[i],0.001,bin90[i]-50,0.001,0.03,"|>"); //,"<|");
    ar[i]->SetLineColor(color[i]);
    ar[i]->SetFillColor(color[i]);
    ar[i]->SetAngle(30);
    
    ar[i]->Draw();

  }
  cFDB2->SetLogy();
  
  TString  outname = outputDir+"Landau_3irr_bestAngle_"+name;
  cFDB2->SaveAs(outname+".eps");
  cFDB2->SaveAs(outname+".png");
  cFDB2->SaveAs(outname+".pdf");
  cFDB2->SaveAs(outname+".root");
  cFDB2->SaveAs(outname+".C");


  /*
  // ####################    residual   ################
  Hist = "dx3_clchargeABC90evR95";
  
  // Non-irr 
  MapTH1 rmsq_map;
  Run[0] = runNONirr;
  i_dphcut[0] = 22; 
  d_dphcut[0] = 22;
  ss_dphcut[0] = "22";
  //  GetHists(&rmsq_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res);
  TH1F * hdx3_clchargeABC90evR; // = new TH1I("dx3_clchargeABC90evR ", "triplet dx_clchargeABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
  hdx3_clchargeABC90evR =    GetCheckHists(&rmsq_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,"testEntries");
  
  //proton

  Run[0] = runPirr2;
  i_dphcut[0] = 15; 
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  //  GetHists(&rmsq_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res);
  hdx3_clchargeABC90evR =    GetCheckHists(&rmsq_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,"testEntries");
  
  

  //neutron

  Run[0] = runNirr4;
  i_dphcut[0] = 18; 
  d_dphcut[0] = 18;
  ss_dphcut[0] = "18";
  //GetHists(&rmsq_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res);
  hdx3_clchargeABC90evR =    GetCheckHists(&rmsq_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,"testEntries");
  
  TH1F  h_resq[irradiations];
  for(int i=0; i<irradiations; i++)
    {
      auto it2 = rmsq_map.find(std::make_pair(runs[i]+"_"+ss_dphcuts[i]+"_"+pitch[0],Hist));
      if(it2  != rmsq_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	  h_resq[i] = *it2->second;
	  cout << "integral " << h_resq[i].Integral() << endl;
	  //	  landau90[i]= 0.9*h_resq[i].Integral();
	  h_resq[i].Scale(1./h_resq[i].Integral());
	}
      else cout << "map key " << it2->first.first << " not found " << endl;
    }


  /////// plots!

  
  TCanvas *cresq = new TCanvas("cresq", "FDB resolution", 600, 600);
  cresq->SetLeftMargin(0.12);  
  cresq->SetRightMargin(-0.1);  
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  cout << "starting making hists pretty" << endl;
  h_resq[0].SetTitle(" ");
  h_resq[0].GetYaxis()->SetTitle("Normalized number of events");
  h_resq[0].GetXaxis()->SetTitle("residual [mm]");
  //  h_resq[0].SetMarkerSize(2.5);
  h_resq[0].SetLineColor(kRed);
  h_resq[0].SetLineStyle(1);
  h_resq[0].SetLineWidth(2);
  h_resq[0].GetXaxis()->SetLimits(-0.5,0.5);
  h_resq[0].GetXaxis()->SetMaxDigits(3); //SetNoExponent(true);
  h_resq[0].GetYaxis()->SetRangeUser(0.00001,10.);
  h_resq[0].GetXaxis()->SetRange(0.,h_resq[0].GetNbinsX() + 1);
  h_resq[0].Draw("histe");
  //  TGaxis::SetMaxDigits(3);
  //  h_resq[1].SetMarkerSize(2.5);
  h_resq[1].GetXaxis()->SetLimits(-0.5,0.5);
  h_resq[1].GetXaxis()->SetRange(0.,h_resq[1].GetNbinsX() + 1);
  h_resq[1].GetYaxis()->SetRangeUser(0.00001,10.);
  h_resq[1].SetLineColor(kBlack);
  h_resq[1].SetLineStyle(2);
  h_resq[1].SetLineWidth(2);
  //  h_resq[1].GetXaxis()->SetNoExponent(true);
  h_resq[1].Draw("histesame");
 				      
  // h_resq[2].SetMarkerSize(2.5);
  h_resq[2].GetXaxis()->SetLimits(-0.5,0.5);
  h_resq[2].GetXaxis()->SetRange(0, h_resq[2].GetNbinsX() + 1);
  h_resq[2].GetYaxis()->SetRangeUser(0.00001,10.);
  //h_resq[2].GetYaxis()->SetRangeUser(0.,0.08);
  h_resq[2].SetLineColor(kBlue);
  h_resq[2].SetLineWidth(2);
  h_resq[2].SetLineStyle(3);
  //  h_resq[2].GetXaxis()->SetNoExponent(true);
  h_resq[2].Draw("histesame");

  TLegend* legresq = new TLegend(0.45,0.72,0.8,0.87);
  legresq->SetLineColor(0);
  legresq->SetTextSize(0.02);
  for(int i =0; i < irradiations; i++)
    legresq->AddEntry(&(h_resq[i]),irr[i],"le");

  legresq->Draw();

  /*
  int color[irradiations] = {632,1,600};
  TLine *  line2[irradiations];
  TArrow *ar[irradiations];
  
  for(int i =0; i < irradiations; i++){
    cout << landau90[i] << endl;
    line2[i]    = new TLine( landau90[i],0.,landau90[i],0.007);
    line2[i]->SetLineColor(color[i]);
    line2[i]->SetLineWidth(2);
    line2[i]->SetLineStyle(1);
    line2[i]->Draw("same");
    ar[i]  = new TArrow(landau90[i],0.006,landau90[i]-6000,0.006,0.03,"|>"); //,"<|");
    ar[i]->SetLineColor(color[i]);
    ar[i]->SetFillColor(color[i]);
    ar[i]->SetAngle(30);
    
    ar[i]->Draw();

  }
  */
  /*
  cresq->SetLogy();
  
  outname = outputDir+"ResQ95_3irr_bestAngle_"+name;
  cresq->SaveAs(outname+".eps");
  cresq->SaveAs(outname+".png");
  cresq->SaveAs(outname+".pdf");
  cresq->SaveAs(outname+".root");
  cresq->SaveAs(outname+".C");
  */

  // ####################    residual   ################
  Hist = "dx3_clphABC90evR";
  TH1F * hdx3_clchargeABC90evR;  
  // Non-irr 
  MapTH1 rmsph_map;
  Run[0] = runNONirr;
  i_dphcut[0] = 12; 
  d_dphcut[0] = 12;
  ss_dphcut[0] = "12";
  //  GetHists(&rmsph_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res);
  //hdx3_clchargeABC90evR =
  rmsph_map = GetCheckHists(&rmsph_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfileNONirr);
  
  //proton

  Run[0] = runPirr2;
  i_dphcut[0] = 15; 
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  //  GetHists(&rmsph_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res);
  //hdx3_clchargeABC90evR =
  rmsph_map = GetCheckHists(&rmsph_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfilePirr2);
  
  

  //neutron

  Run[0] = runNirr4;
  i_dphcut[0] = 15; 
  d_dphcut[0] = 15;
  ss_dphcut[0] = "15";
  //GetHists(&rmsph_map,1, 1, 1,dphcut, Run, i_dphcut, Label, Hist, h_res);
  //hdx3_clchargeABC90evR =
  rmsph_map = GetCheckHists(&rmsph_map, 1, 1,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR,true,inputfileNirr4);



  TH1F  h_resq[irradiations];   
  //  TH1F  h_resq[irradiations];
  for(int i=0; i<irradiations; i++)
    {
      auto it2 = rmsph_map.find(std::make_pair(runs[i]+"_"+ss_dphcuts[i]+"_"+pitch,Hist));
      if(it2  != rmsph_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	  h_resq[i] = *it2->second;
	  cout << "integral " << h_resq[i].Integral() << endl;
	  //	  landau90[i]= 0.9*h_resq[i].Integral();
	  h_resq[i].Scale(1./h_resq[i].Integral());
	  h_resq[i].Draw();
	  ScaleAxis(h_resq[i].GetXaxis(), ScaleX);
	  
	}
      else cout << "map key " << it2->first.first << " not found " << endl;
    }


  /////// plots!

  
  TCanvas *cresph = new TCanvas("cresph", "FDB resolution", 600, 600);
  cresph->SetLeftMargin(0.12);  
  cresph->SetRightMargin(-0.1);  
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  
  cout << "starting making hists pretty" << endl;
  h_resq[0].SetTitle(" ");
  h_resq[0].GetYaxis()->SetTitle("Normalized number of events");
  h_resq[0].GetXaxis()->SetTitle("#Deltax [#mum]");
  //  h_resq[0].SetMarkerSize(2.5);
  h_resq[0].SetLineColor(kBlack);
  h_resq[0].SetLineStyle(1);
  h_resq[0].SetLineWidth(2);
  h_resq[0].GetXaxis()->SetNdivisions(20,5,3);
  h_resq[0].GetXaxis()->SetRangeUser(-0.1,0.1);
  h_resq[0].GetXaxis()->SetMaxDigits(3); //SetNoExponent(true);
  h_resq[0].GetYaxis()->SetRangeUser(0.00001,10.);
  h_resq[0].GetXaxis()->SetRange(0.,h_resq[0].GetNbinsX() + 1);
  h_resq[0].Draw("histe");
  cout << " Xaxis " << h_resq[0].GetXaxis()->GetNdivisions() << endl;
  cout << " RMS " << h_resq[0].GetRMS() << endl;

  //  TGaxis::SetMaxDigits(3);
  //  h_resq[1].SetMarkerSize(2.5);
  h_resq[1].GetXaxis()->SetRangeUser(-0.1,0.1);
  h_resq[1].GetXaxis()->SetRange(0.,h_resq[1].GetNbinsX() + 1);
  h_resq[1].GetYaxis()->SetRangeUser(0.00001,10.);
  h_resq[1].SetLineColor(kGreen+1);
  h_resq[1].SetLineStyle(2);
  h_resq[1].SetLineWidth(2);
  //  h_resq[1].GetXaxis()->SetNoExponent(true);
  h_resq[1].Draw("histesame");
 				      
  // h_resq[2].SetMarkerSize(2.5);
  h_resq[2].GetXaxis()->SetRangeUser(-0.1,0.1);
  h_resq[2].GetYaxis()->SetRangeUser(0.00001,10.);
  h_resq[2].GetXaxis()->SetRange(0.,h_resq[2].GetNbinsX() + 1);
  h_resq[2].SetLineColor(kMagenta);
  h_resq[2].SetLineWidth(2);
  h_resq[2].SetLineStyle(3);
  //  h_resq[2].GetXaxis()->SetNoExponent(true);
  h_resq[2].Draw("histesame");
  
  TLegend* legresq = new TLegend(0.35,0.72,0.8,0.87);
  legresq->SetLineColor(0);
  //  legresq->SetTextSize(0.03);
  for(int i =0; i < irradiations; i++)
    legresq->AddEntry(&(h_resq[i]),irr[i],"le");
  
  legresq->Draw();

  /*
  int color[irradiations] = {632,1,600};
  TLine *  line2[irradiations];
  TArrow *ar[irradiations];
  
  for(int i =0; i < irradiations; i++){
    cout << landau90[i] << endl;
    line2[i]    = new TLine( landau90[i],0.,landau90[i],0.007);
    line2[i]->SetLineColor(color[i]);
    line2[i]->SetLineWidth(2);
    line2[i]->SetLineStyle(1);
    line2[i]->Draw("same");
    ar[i]  = new TArrow(landau90[i],0.006,landau90[i]-6000,0.006,0.03,"|>"); //,"<|");
    ar[i]->SetLineColor(color[i]);
    ar[i]->SetFillColor(color[i]);
    ar[i]->SetAngle(30);
    
    ar[i]->Draw();

  }
  */

  cresph->SetLogy();
  
  outname = outputDir+"ResPHAll_3irr_bestAngle_"+name;
  cresph->SaveAs(outname+".eps");
  cresph->SaveAs(outname+".png");
  cresph->SaveAs(outname+".pdf");
  cresph->SaveAs(outname+".root");
  cresph->SaveAs(outname+".C");











  
}//resolution 



Double_t ScaleX(Double_t x)
{
  Double_t v;
  v = 1000 * x; // "linear scaling" function example
  return v;
}

void ScaleAxis(TAxis *a, Double_t (*Scale)(Double_t))
{
  if (!a) return; // just a precaution
  if (a->GetXbins()->GetSize())
    {
      // an axis with variable bins
      // note: bins must remain in increasing order, hence the "Scale"
      // function must be strictly (monotonically) increasing
      TArrayD X(*(a->GetXbins()));
      for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i]);
      a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
    }
  else
    {
      // an axis with fix bins
      // note: we modify Xmin and Xmax only, hence the "Scale" function
      // must be linear (and Xmax must remain greater than Xmin)
      a->Set( a->GetNbins(),
	      Scale(a->GetXmin()), // new Xmin
	      Scale(a->GetXmax()) ); // new Xmax
    }
  return;
}



void TDR()
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_sqrtS = "";
  bool drawLogo = true;
  int H_ref = 600;
  int W_ref = 600;
  float T = 0.07*H_ref;//0.07
  float B = 0.11*H_ref;//0.12
  float L = 0.12*W_ref;
  float R = 0.01*W_ref;


}

void TDR2(TCanvas * c_all, int period = 0, int pos = 11)
{
  CMS_lumi( c_all, period,pos);
  //  CMS_lumi( c_all, 0, 11);
  c_all->Update();
  c_all->RedrawAxis();
}

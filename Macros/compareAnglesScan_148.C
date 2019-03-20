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
#include "resolution.h"

#define AnglesExtra 3
#define Angles 29
bool print=true;
using namespace std;


void compareAnglesScan_148()
{

  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  TString detectorA="146";
  TString detectorB="148";
  TString detectorC="163";
  TString labelA="FDB150Y_2_R4S100x25-Y2_2, thr 12 ADC";
  TString labelB="FDB150P_12_R4S100x25-P4_1, 120V, thr 12 ADC";
  TString labelC="FDB150Y_2_R4S100x25-Y6_1, thr 12 ADC";
  TString info = "beam energy 5.6 GeV";


  Float_t AngleExtra[AnglesExtra]={9.5,10.5,11.5};
  Float_t Angle[Angles]={-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
  TString ss_Angle[Angles];
  TString ss_File[Angles];

  float angleunit = 1.25;
  
  Int_t run[Angles]={2731,2732,2733,2734,2735,2736,2737,2738,2739,2740,2741,2743,2744,2745,2746,2747,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2760};//,2764,2766,2767};
  Int_t runExtra[AnglesExtra]={2767,2766,2764};
  


  /*  TString sensor[measurements];
  TString irr[measurements];
  Double_t Resolution[measurements];
  TString pitch[measurements];
  TString bias[measurements];
  TString angle[measurements];
  TString beam[measurements];

  Int_t runs[measurements];
  TString specific[measurements];
  TString Short[measurements];
  TString thr[measurements];
  TString chargecut[measurements];
  TString gain[measurements];
  TString typeunf[measurements];
  TString note[measurements];

  TString filename = inputDir+"/TextFiles/resolution_RMS.txt";
  cout << filename << endl;
  ifstream stream(filename);
  int i=0;
  std::string line;
  if(!stream.is_open())
    {
      cout << " File " << filename << " not opened" << endl;
    }
  else
    {
      std::getline(stream,line);
      if(print) cout << " first line " << line << endl;
      //      while(!stream.eof())
      for(i= 0; i < measurements; i++)
	{
	  stream  >> sensor[i] >> irr[i] >> pitch[i] >> bias[i] >> angle[i] >> beam[i] >> Resolution[i] >> runs[i] >> specific[i] >> Short[i] >> thr[i] >> chargecut[i] >> gain[i] >> typeunf[i] >> note[i];
	  if(print)       cout << "line " << i << endl;
	  if(print)       cout  << sensor[i] << " " << irr[i] << " " << pitch[i] << " " << bias[i] << " " << angle[i] << " " << beam[i] << " " << Resolution[i] << " " << runs[i] << " " << specific[i] << " " << Short[i] << " " << thr[i] << " " << chargecut[i] << " " << gain[i] << " " << typeunf[i] << " " << note[i] <<endl;
	  //  i++;

	  cout << " resolution for sensor " << sensor[i] << " is " << Resolution[i] << endl;
	}
    }
  */


  TString inputfile;
  TFile * file[Angles+AnglesExtra];
  TH1F * h_ncol[Angles+AnglesExtra];
  TH1F * h_madx[Angles+AnglesExtra];
  TH1F * h_madxvsq[Angles+AnglesExtra];

  TString Run, Path;
  ostringstream strs[Angles+AnglesExtra];
  
  if(print) cout << "Getting files and hists"  << endl;
  for(int i=0; i<Angles+AnglesExtra; i++)
    {
      ss_File[i].Form("%d",i);
      if(i<Angles)
	{
	  if(print) cout << "Run: " << run[i] << endl;
	  Run.Form("%d",run[i]);
	  Angle[i]=Angle[i]*angleunit; //grad
	  stringstream stream;
	  stream << fixed << setprecision(2) << Angle[i];
	  string s = stream.str();
	  ss_Angle[i] = s;
	  //	  ss_Angle[i].Form("%f",Angle[i]);
      
	  if(print) cout << ss_Angle[i] << " : " << Angle[i] << endl;
	}
      if(i>=Angles)
	{
	  if(print) cout << "Run: " << runExtra[i-Angles] << endl;
	  Run.Form("%d",runExtra[i-Angles]);
	  AngleExtra[i-Angles]=AngleExtra[i-Angles]*angleunit; //grad
	  stringstream stream;
	  stream << fixed << setprecision(2) << AngleExtra[i-Angles];
	  string s = stream.str();
	  ss_Angle[i] = s;
	 
	  //	  ss_Angle[i].Form("%f",AngleExtra[i-Angles]);
	  if(print) cout << ss_Angle[i] << " : " << AngleExtra[i-Angles] << endl;
	}	      
	  if(print) cout << "Run: " << Run << endl;
	  inputfile = "drei-r"+Run+"_irene.root";
	  if(print) cout << "File Name: " << inputfile << endl;
	  Path=inputDir+inputfile;
	  file[i] = new TFile(Path);
	  if(print) cout << "File Path: " << Path << endl;
	  TString hist = "nrowvsxmB3";
	  h_ncol[i] = (TH1F*)file[i]->Get(hist);
	  if(print) cout << hist << " " << h_ncol[i]->GetEntries() <<endl;
	  hist = "madx3vsxm";
	  h_madx[i] = (TH1F*)file[i]->Get(hist);
	  if(print) cout << hist << " " << h_madx[i]->GetEntries() <<endl;
	  hist = "madx3vsq";
	  h_madxvsq[i] = (TH1F*)file[i]->Get(hist);
	  if(print) cout << hist << " " << h_madxvsq[i]->GetEntries() <<endl;
    }
 
if(print) cout << "Plotting nrow big scan "  << endl;
  //i da 0 a 6 

  TCanvas *c2 = new TCanvas("c2", "nrow vs xm - all scan", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend* leg2 = new TLegend(0.15,0.15,0.75,0.3);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.018);


  for(int i = 5; i<Angles+AnglesExtra ; i++)
    {
      h_ncol[i]->SetLineWidth(2);
      if(i>=Angles) h_ncol[i]->SetLineStyle(2);

      int color = i-4;//i+1;
      if(color == 5) color = 38;
      if(color == 10) color = 46;
      if(color == 19) color = 39;
      if(i>=Angles) color = i+11;
      h_ncol[i]->SetLineColor(color);
      h_ncol[i]->GetYaxis()->SetRangeUser(0.,4.5);
      h_ncol[i]->Draw("histsame");
      //      leg2->AddEntry(h_ncol[i],pitch[i]+" #mum pitch, "+Short[i]+" "+bias[i]+"V, irr "+irr[i]+"x10^{15} neq" , "l");
      leg2->AddEntry(h_ncol[i],ss_Angle[i], "l");
    }
  leg2->SetNColumns(5);
  leg2->Draw();
  TString name = outputDir+"Nrow_vs_xmod_148_allScan";
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".pdf");


  if(print) cout << "Plotting nrow scan "  << endl;
  //i da 0 a 6 

  TCanvas *c3 = new TCanvas("c3", "nrow vs xm - investigate scan", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend* leg3 = new TLegend(0.15,0.15,0.75,0.3);
  leg3->SetLineColor(0);
  leg3->SetTextSize(0.018);


  for(int i = 13; i<Angles+AnglesExtra ; i++)
    {
      h_ncol[i]->SetLineWidth(2);
      if(i>=Angles) h_ncol[i]->SetLineStyle(2);

      int color = i-4;//i+1;
      if(color == 5) color = 38;
      if(color == 10) color = 46;
      if(color == 19) color = 39;
      if(i>=Angles) color = i+11;
      h_ncol[i]->SetLineColor(color);
      h_ncol[i]->GetYaxis()->SetRangeUser(0.,4.5);
      if(i < 17 || i >=Angles) h_ncol[i]->Draw("histsame");
      //      leg2->AddEntry(h_ncol[i],pitch[i]+" #mum pitch, "+Short[i]+" "+bias[i]+"V, irr "+irr[i]+"x10^{15} neq" , "l");
      if(i < 17 || i >=Angles)      leg3->AddEntry(h_ncol[i],ss_Angle[i], "l");
    }
  leg3->SetNColumns(4);
  leg3->Draw();
  name = outputDir+"Nrow_vs_xmod_148_investigateScan";
  c3->SaveAs(name+".eps");
  c3->SaveAs(name+".png");
  c3->SaveAs(name+".pdf");

  
  
    if(print) cout << "Plotting madx "  << endl;
  //i da 0 a 6 

  TCanvas *c4 = new TCanvas("c4", "madx vs xm ", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend* leg4 = new TLegend(0.25,0.7,0.75,0.8);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.03);


  for(int i = 5; i<Angles+AnglesExtra ; i++)
    {
      h_madx[i]->SetLineWidth(2);
      if(i>=Angles) h_madx[i]->SetLineStyle(2);

      int color = i-4;
      if(color == 5) color = 38;
      if(color == 10) color = 46;
      if(color == 19) color = 39;
      if(i>=Angles) color = i+11;
      h_madx[i]->SetLineColor(color);
      
      h_madx[i]->GetYaxis()->SetRangeUser(0.001,0.009);
  
      if(i<14) h_madx[i]->Draw("histsame");
      //      leg2->AddEntry(h_madx[i],pitch[i]+" #mum pitch, "+Short[i]+" "+bias[i]+"V, irr "+irr[i]+"x10^{15} neq" , "l");
      if(i<14) leg4->AddEntry(h_madx[i],ss_Angle[i], "l");
    }
  
  leg4->SetNColumns(4);
  leg4->Draw();
  name = outputDir+"Madx_vs_xmod_148_allScan";
  c4->SaveAs(name+".eps");
  c4->SaveAs(name+".png");
  c4->SaveAs(name+".pdf");


  
  if(print) cout << "Plotting madx vs q "  << endl;
  //i da 0 a 6 

  TCanvas *c6 = new TCanvas("c6", "madx vs q ", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend* leg6 = new TLegend(0.15,0.75,0.8,0.85);
  leg6->SetLineColor(0);
  leg6->SetTextSize(0.025);

  for(int i = 5; i<Angles ; i++)
    {
      h_madxvsq[i]->SetLineWidth(2);
      if(i>=Angles) h_madxvsq[i]->SetLineStyle(2);

      int color = i-4;
      if(color == 5) color = 38;
      if(color == 10) color = 46;
      if(color == 19) color = 39;
      if(i>=Angles) color = i+11;
      h_madxvsq[i]->SetLineColor(color);


      h_madxvsq[i]->GetYaxis()->SetRangeUser(0.,0.012);
      h_madxvsq[i]->GetXaxis()->SetRangeUser(0.,40.);
      h_madxvsq[i]->Draw("histsame");
      //      leg2->AddEntry(h_madxvsq[i],pitch[i]+" #mum pitch, "+Short[i]+" "+bias[i]+"V, irr "+irr[i]+"x10^{15} neq" , "l");
      leg6->AddEntry(h_madxvsq[i],ss_Angle[i], "l");
    }

  leg6->SetNColumns(6);
  
  
  leg6->Draw();
  name = outputDir+"Madx_vs_q_148_allScan";
  c6->SaveAs(name+".eps");
  c6->SaveAs(name+".png");
  c6->SaveAs(name+".pdf");

}//resolution 

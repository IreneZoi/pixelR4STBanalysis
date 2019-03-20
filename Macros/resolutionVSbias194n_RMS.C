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

#define BiasVoltages 8
#define Angles 8
bool print=true;
using namespace std;
#define dphcuts  1
#define comparisons 1


void resolutionVSbias194n_RMS(TString func = "RMS")
{

  TString detectorA="148";
  TString detectorB="194i";
  TString detectorC="150";
  TString labelA="FDB150P_12_R4S100x25-P4_1, thr 12 ADC";
  TString labelB="FDB150P_12_R4S100x25-P1_1, 120 V, thr 15 ADC";
  TString labelC="FDB150P_12_R4S100x25-P1_3, thr 12 ADC";

  
  
  double BiasVoltage[BiasVoltages];
  BiasVoltage[0]=200;//V  
  BiasVoltage[1]=250;//V
  BiasVoltage[2]=300;//V
  BiasVoltage[3]=400;//V
  BiasVoltage[4]=500;//V
  BiasVoltage[5]=600;//V
  BiasVoltage[6]=700;//V
  BiasVoltage[7]=800;//V  

  double BiasVoltageError[BiasVoltages];

  double Angle[Angles]={8,10,12,14,16,18,20,22};
  TString ss_Angle[Angles];//={10};//,12};//,14,16,18,20};
  
  Double_t Resolution[Angles][BiasVoltages];
  Double_t ResolutionError[Angles][BiasVoltages];
  TString ss_BiasVoltage[BiasVoltages];
  
  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";
  TString inputfile;

  TH1F * h_res;

  // Int_t run[BiasVoltages]={3763,3764,3765,3766,3767,3768,3769,3770};
  Int_t run[Angles][BiasVoltages];

  //8
  run[0][7] = 3823;
  run[0][6] = 3824;
  run[0][5] = 3825;
  run[0][4] = 3826;
  run[0][3] = 3827;
  run[0][2] = 3828;
  run[0][1] = 3830;
  run[0][0] = 3829;
  //10
  run[1][0] = 3770;
  run[1][1] = 3769;
  run[1][2] = 3768;
  run[1][3] = 3767;
  run[1][4] = 3766;
  run[1][5] = 3765;
  run[1][6] = 3764;
  run[1][7] = 3763;
  
  //12
  run[2][0] = 3783;
  run[2][1] = 3782;
  run[2][2] = 3781;
  run[2][3] = 3780;
  run[2][4] = 3779;
  run[2][5] = 3777;
  run[2][6] = 3776;
  run[2][7] = 3778;
  
  //14
  run[3][0] = 3791;
  run[3][1] = 3790;
  run[3][2] = 3789;
  run[3][3] = 3788;
  run[3][4] = 3787;
  run[3][5] = 3786;
  run[3][6] = 3785;
  run[3][7] = 3784;
  //16
  run[4][0] = 3799;
  run[4][1] = 3798;
  run[4][2] = 3797;
  run[4][3] = 3796;
  run[4][4] = 3795;
  run[4][5] = 3794;
  run[4][6] = 3793;
  run[4][7] = 3792;
  
  //18
  run[5][0] = 3813;
  run[5][1] = 3812;
  run[5][2] = 3811;
  run[5][3] = 3810;
  run[5][4] = 3809;
  run[5][5] = 3808;
  run[5][6] = 3807;
  run[5][7] = 3806;

  //20
  run[6][0] = 3821;
  run[6][1] = 3820;
  run[6][2] = 3819;
  run[6][3] = 3818;
  run[6][4] = 3817;
  run[6][5] = 3816;
  run[6][6] = 3815;
  run[6][7] = 3814;
  
//22
  run[7][7] = 3831;
  run[7][6] = 3832;
  run[7][5] = 3833;
  run[7][4] = 3834;
  run[7][3] = 3835;
  run[7][2] = 3836;
  run[7][1] = 3837;
  run[7][0] = 3838;
  

  TString Run[Angles][BiasVoltages];
  ostringstream strs[BiasVoltages];
  ostringstream strsa[Angles];
  
  if(print) cout << "Iniazializating bias voltages and runs "  << endl;

  MapTH1 res_map[Angles];
  int i_dphcut[dphcuts] = {12};
  TString ss_dphcut[dphcuts] = {"12"};

  TString Label[comparisons] = {"straightTracksY_isoAandCandB_straightTracksX"};
  TString Hist = "dx3_clchargeABC90evR";
  bool dphcut = false;

  for(int j = 0; j < Angles; j++)
    {
      strsa[j] << Angle[j];
      ss_Angle[j]=strsa[j].str();

      for(int i=0; i<BiasVoltages; i++)
	{ 
	  
	  if(print) cout << "Run: "<< i <<" " << run[j][i] << endl;
	  Run[j][i].Form("%d",run[j][i]);
	  if(print) cout << "Run: " << Run[j][i] << endl;
	  // ss_BiasVoltage[i].Form("%d",BiasVoltage[i]);	  
	  if(j==0)
	    {
	      BiasVoltageError[i]=0; //GeV
	      if(print) cout << "Beam energy " << i<< ": " << BiasVoltage[i] << " GeV" << endl;
	      strs[i] << BiasVoltage[i];
	      ss_BiasVoltage[i]=strs[i].str();
	      if(print) cout << "Beam energy " << i<< ": " << ss_BiasVoltage[i] << " GeV" << endl;
	    }
	}
      GetHists(&(res_map[j]),BiasVoltages, dphcuts, comparisons,dphcut, Run[j], i_dphcut, Label, Hist, h_res);
    }
  


  
  
  

  if(print) cout << "Fit:"  << endl;
  

  ofstream myfile[Angles];
  for(int j = 0; j < Angles; j++)
    {
      
      myfile[j].open ("/home/zoiirene/Output/TextFiles/Bscan_"+detectorB+"_"+func+"_"+ss_Angle[j]+"deg_rms.txt");
      myfile[j] << "A " << detectorA << "\n";
      myfile[j] << labelA << "\n";
      myfile[j] << "B " << detectorB<< "\n";
      myfile[j] << labelB << "\n";
      myfile[j] << "C " << detectorC<< "\n";
      myfile[j] << labelC << "\n";

      myfile[j] << "Momentum(GeV) RMS(um) Error\n";
      
      for(int i=0; i<BiasVoltages; i++)
	{
	  auto it2 = res_map[j].find(std::make_pair(Run[j][i]+"_"+Label[0],Hist));
	  if(it2  != res_map[j].end())
	    {
	      if(print)               cout << " found map " << endl;
	      if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	      
	      FitTH1(it2->second, &(Resolution[j][i]), &(ResolutionError[j][i]), ss_Angle[j]+"_"+ss_BiasVoltage[i], detectorA, detectorB, detectorC,func );
	      if(print) cout << "Beam energy " << i<< ": " << ss_BiasVoltage[i] << " GeV -> Resolution: " << Resolution[j][i] << " and res err: " << ResolutionError[j][i] << endl;
	      myfile[j] << ss_BiasVoltage[i] << " " << Resolution[j][i] << " " << ResolutionError[j][i] << "\n";
	    }
	}
      
      myfile[j].close();
      DrawTGraphWithErrorDouble(BiasVoltages, BiasVoltage, Resolution[j], ResolutionError[j], Hist,Run[j][0], Label[0], "biasScan_"+ss_Angle[j]+"_rms", "RMS 95% [#mum]",0.,7.,0.,8., "biasScan_"+ss_Angle[j], "Bias Voltage [V]"  );
    }
  

  //////////////          RESOLUTION ///////////////////////
  double freshres = 2.03;
  double freshres_err = 0.01;

  for(int j = 0; j < Angles; j++)
    {
      
      for(int i=0; i<BiasVoltages; i++)
	{
	  
	  if(print) cout << "Energy " << i<< ": " << ss_BiasVoltage[i] << " GeV -> Resolution: " << Resolution[j][i] << " and res err: " << ResolutionError[j][i] << endl;
	  ExtractRes(&(Resolution[j][i]),&(ResolutionError[j][i]),true, freshres, freshres_err);
	  if(print) cout << "Energy " << i<< ": " << ss_BiasVoltage[i] << " GeV -> Resolution: " << Resolution[j][i] << " and res err: " << ResolutionError[j][i] << endl;
	  
	  
	}
      DrawTGraphWithErrorDouble(BiasVoltages, BiasVoltage, Resolution[j], ResolutionError[j], Hist,Run[j][0], Label[0], "biasScan_"+ss_Angle[j]+"_res", "Resolution [#mum]",0.,7.,0.,8., "biasScan_"+ss_Angle[j], "Bias Voltage [V]"  );
    }
  

  Double_t resolution[BiasVoltages][Angles];
  Double_t resolutionError[BiasVoltages][Angles];
  


  double errx[Angles];
  for(int i =0; i< Angles;i++)
    {
      errx[i]=0;
      for(int j=0; j<BiasVoltages; j++)
	{
	  cout << Angle[i] << " " << BiasVoltage[j] << " " << Resolution[i][j] << endl;
	  resolution[j][i] = Resolution[i][j];
	  resolutionError[j][i] = ResolutionError[i][j];
	  cout << Angle[i] << " " << BiasVoltage[j] << " " << resolution[j][i] << endl;
	  
	}
    }

  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //  gStyle->SetOptStat();
  //gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  int color;

  TLegend* leg2 = new TLegend(0.3,0.75,0.85,0.85);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.04);
  

  

  TGraphErrors * gr[BiasVoltages];
  for(int j=0; j<BiasVoltages; j++)
    {
      gr[j]  = new TGraphErrors(Angles, Angle, resolution[j],errx,resolutionError[j]);//, errx, rms);
      color = j+1;
      if(j+1 ==5) color = 95;
      if(j+1 ==10) color = 29;
      
      gr[j]->SetLineColor(color);
      gr[j]->SetLineWidth(1);
      gr[j]->SetLineStyle(2);
      gr[j]->SetMarkerColor(color);
      gr[j]->SetMarkerSize(1.5);
      gr[j]->SetMarkerStyle(21);
      gr[j]->GetXaxis()->SetRangeUser(8.,22.);
      gr[j]->SetMinimum(0.);
      gr[j]->SetMaximum(10.);
      
      //  gr[j]->SetTitle("Option ACP example");
      gr[j]->GetXaxis()->SetTitle("Angle [deg]");
      gr[j]->GetYaxis()->SetTitle("Resolution [#mum]");
      //gr[j]->SetFillStyle(3003);
      //gr[j]->SetFillColor(kRed-8);
            leg2->AddEntry(gr[j],ss_BiasVoltage[j], "lp");      
	    cout  << "try  " << BiasVoltage[j] << " " << resolution[j][0] <<" " << resolution[j][1] << endl;

	    if(j==0)   gr[j]->Draw("ACPE");
	    else gr[j]->Draw("CPEsame");      
    }

  leg2->SetNColumns(4);

  leg2->Draw("same");
  TString name;
  name = outputDir+"compare_angle_bias_Scan_"+detectorB+"_"+Label[0]+"_"+Hist;
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");


  ////// zoom in 
  
  TCanvas *c2z = new TCanvas("c2z", "c2z", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //  gStyle->SetOptStat();
  //gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  
  for(int j=0; j<BiasVoltages; j++)
    {
      gr[j]  = new TGraphErrors(Angles, Angle, resolution[j],errx,resolutionError[j]);//, errx, rms);
      color = j+1;
      if(j+1 ==5) color = 95;
      if(j+1 ==10) color = 29;

      gr[j]->SetLineColor(color);
      gr[j]->SetLineWidth(1);
      gr[j]->SetLineStyle(2);
      gr[j]->SetMarkerColor(color);
      gr[j]->SetMarkerSize(1.5);
      gr[j]->SetMarkerStyle(21);
      gr[j]->GetXaxis()->SetRangeUser(8.,22.);
      gr[j]->SetMinimum(4.);
      gr[j]->SetMaximum(7.);

      //  gr[j]->SetTitle("Option ACP example");
      gr[j]->GetXaxis()->SetTitle("Angle [deg]");
      gr[j]->GetYaxis()->SetTitle("Resolution [#mum]");
      //gr[j]->SetFillStyle(3003);
      //gr[j]->SetFillColor(kRed-8);
      //      leg2->AddEntry(gr[j],ss_BiasVoltage[j], "lp");
      cout  << "try  " << BiasVoltage[j] << " " << resolution[j][0] <<" " << resolution[j][1] << endl;

      if(j==0)   gr[j]->Draw("ACPE");
      else gr[j]->Draw("CPEsame");
    }

  
  leg2->SetTextSize(0.03);
  leg2->Draw("same");
  
  name = outputDir+"compare_angle_bias_Scan_zoom_"+detectorB+"_"+Label[0]+"_"+Hist;
  c2z->SaveAs(name+".eps");
  c2z->SaveAs(name+".pdf");
  c2z->SaveAs(name+".png");
  c2z->SaveAs(name+".root");


  ////// cluster size
  Hist = "nrowB";

  for(int j = 0; j < Angles; j++)
    GetHists(&(res_map[j]),BiasVoltages, dphcuts, comparisons,dphcut, Run[j], i_dphcut, Label, Hist, h_res);



  TCanvas *csza[Angles];/// = new TCanvas("csza", "csza", 1500, 900);

  double NrowB[Angles][BiasVoltages];



  TLegend* leg2a = new TLegend(0.7,0.35,0.85,0.85);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.04);

  
  for(int j = 0; j < Angles; j++)
    {
      csza[j] = new TCanvas("csza", "csza", 1500, 900);  
      gPad->SetTicks(1,1);
      gROOT->SetStyle("Plain");
      gStyle->SetPadGridX(0);
      gStyle->SetPadGridY(0);
      gStyle->SetPalette(1);

      //gStyle->SetOptStat(1110);
      gStyle->SetOptTitle(0);
      
      for(int i=0; i<BiasVoltages; i++)
	{
	  gStyle->SetOptStat(0);
	  auto it2 = res_map[j].find(std::make_pair(Run[j][i]+"_"+Label[0],Hist));
	  if(it2  != res_map[j].end())
	    {
	      if(print)  cout << " found map " << endl;
	      if(print)  cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetMean() << endl;
	      NrowB[j][i] = it2->second->GetMean();
	    }
	  color = i+1;
	  if(i+1 ==5) color = 95;
	  if(i+1 ==10) color = 29;

	  it2->second->SetLineColor(color);
	  it2->second->SetLineWidth(1);
	  it2->second->SetLineStyle(1);
	  //	  it2->second->SetMarkerColor(color);
	  //it2->second->SetMarkerSize(1.5);
	  //it2->second->SetMarkerStyle(21);
	  it2->second->GetXaxis()->SetRangeUser(0.,5.);
	  it2->second->SetMinimum(0.);
	  it2->second->SetMaximum(40000.);

	  //  it2->second->SetTitle("Option ACP example");
	  it2->second->GetXaxis()->SetTitle("cluster size B");
	  it2->second->GetYaxis()->SetTitle("Entries");
	  //it2->second->SetFillStyle(3003);
	  //it2->second->SetFillColor(kRed-8);
	  //      leg2->AddEntry(it2->second,ss_BiasVoltage[j], "lp");
	  if(j==0)          leg2a->AddEntry(it2->second,ss_BiasVoltage[i], "l");

	  it2->second->Draw("histsame");
	  //	  it2->second->Draw("CPEsame");

	  
	}
      leg2a->SetNColumns(2);
      leg2a->SetTextSize(0.03);
      leg2a->Draw("same");

      name = outputDir+"compare_angle"+ss_Angle[j]+"_bias_Scan_clustersize_"+detectorB+"_"+Label[0]+"_"+Hist;
      csza[j]->SaveAs(name+".eps");
      csza[j]->SaveAs(name+".pdf");
      csza[j]->SaveAs(name+".png");
      csza[j]->SaveAs(name+".root");
    }


  
  TCanvas *csz = new TCanvas("csz", "csz", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //  gStyle->SetOptStat();
  //gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);


  TGraph * grB[BiasVoltages];
  for(int j=0; j<BiasVoltages; j++)
    {
       
      grB[j]  = new TGraph(Angles, Angle, NrowB[j]);
      color = j+1;
      if(j+1 ==5) color = 95;
      if(j+1 ==10) color = 29;

      grB[j]->SetLineColor(color);
      grB[j]->SetLineWidth(1);
      grB[j]->SetLineStyle(2);
      grB[j]->SetMarkerColor(color);
      grB[j]->SetMarkerSize(1.5);
      grB[j]->SetMarkerStyle(21);
      grB[j]->GetXaxis()->SetRangeUser(8.,22.);
      grB[j]->SetMinimum(0.);
      grB[j]->SetMaximum(4.);

      //  grB[j]->SetTitle("Option ACP example");
      grB[j]->GetXaxis()->SetTitle("Angle [deg]");
      grB[j]->GetYaxis()->SetTitle("Mean cluster size B");
      //grB[j]->SetFillStyle(3003);
      //grB[j]->SetFillColor(kRed-8);
      //      leg2->AddEntry(grB[j],ss_BiasVoltage[j], "lp");


      if(j==0)   grB[j]->Draw("ACPE");
      else grB[j]->Draw("CPEsame");
	
    }

  
  leg2->SetTextSize(0.03);
  leg2->Draw("same");
  
  name = outputDir+"compare_angle_bias_Scan_clustersize_"+detectorB+"_"+Label[0]+"_"+Hist;
  csz->SaveAs(name+".eps");
  csz->SaveAs(name+".pdf");
  csz->SaveAs(name+".png");
  csz->SaveAs(name+".root");










  
}//resolution 


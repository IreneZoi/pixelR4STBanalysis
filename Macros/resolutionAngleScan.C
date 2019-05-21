
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

bool print=true;
using namespace std;
#define anglesNONirr 29
#define anglesPirr2 28
#define anglesNirr4 9
#define IEEE 2
#define irradiations 3

void resolutionAngleScan(TString name = "preliminary", TString func = "RMS", bool unfolding = false)
{



  // 194i  4 25   800   12   5.6   3.86636   0.0322234   3839   FDB150P_12_R4S100x25-P1_1   Pstop_default_FDB   15 

  
  TString inputDir="/home/zoiirene/Output/TextFiles/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  int measurements[irradiations+2];
  measurements[0] = anglesNONirr;
  measurements[1] = anglesPirr2;
  measurements[2] = anglesNirr4;
  measurements[3] = 1;
  measurements[4] = anglesNONirr;
  
  TString filenames[irradiations+1];
  filenames[0] = "Ascan_res_148_RMS.txt";
  filenames[1] = "Ascan_res_120i_RMS.txt";
  filenames[2] = "Ascan_194i_RMS_800_res.txt";
  filenames[3] = "Ascan_pixelav_RMS.txt";

  TString filenames2[IEEE];
  filenames2[0] = "Ascan_clsizeB_148_RMS.txt";
  filenames2[1] = "Ascan_clsizeB_120i_RMS.txt";

  
  TString irr[irradiations+2];
  irr[0] = "no irr, 5.6 GeV";
  irr[1] = "proton irr at #phi_{eq}=2#times10^{15} cm^{-2}, 5.6 GeV";
  irr[2] = "neutron irr at #phi_{eq}=4#times10^{15} cm^{-2}, 5.2 GeV";
  irr[3] = "neutron irr at #phi_{eq}=4#times10^{15} cm^{-2}, 5.6 GeV";
  irr[4] = "no irr, pixelav";
  
  double angle_0[anglesNONirr];
  double angle_1[anglesPirr2];
  double angle_2[anglesNirr4];
  double angle_3[] = {12};
  double angle_4[anglesNONirr];

  double angleerr_0[anglesNONirr];
  double angleerr_1[anglesPirr2];
  double angleerr_2[anglesNirr4];
  double angleerr_3[]= {1};
  double angleerr_4[anglesNONirr];

  double res_0[anglesNONirr];
  double res_1[anglesPirr2];
  double res_2[anglesNirr4];
  double res_3[] = {3.87};
  double res_4[anglesNONirr];

  double reserr_0[anglesNONirr];
  double reserr_1[anglesPirr2];
  double reserr_2[anglesNirr4];
  double reserr_3[] = {0.03};
  double reserr_4[anglesNONirr];
  
  double angle;
  double res,reserror,angleerr;
  
  TString filename[irradiations];
  for(int i =0; i< irradiations; i++)
    {
      filename[i]      = inputDir+filenames[i];
      cout << filename[i] << endl;
      ifstream stream(filename[i]);
      
      std::string line;
      if(!stream.is_open())
	{
	  cout << " File " << filename[i] << " not opened" << endl;
	}
      else
	{
	  for(int k =0; k<8;k++)
	    {
	      std::getline(stream,line);
	      if(print) cout << k << " line " << line << endl;
	    }
	  //      while(!stream.eof())
	  for(int j= 0; j < measurements[i]; j++)
	    {
	      stream  >> angle >> res >> reserror ;
	      if(print)	 cout << "line " << j+7 << endl;
	      if(i==0)
		{
		  angle_0[j] = angle;
		  angleerr_0[j] = 1;
		  res_0[j] = res;
		  reserr_0[j] = reserror;
		  if(print)	 cout  << angle_0[j] << " " << res_0[j] << " " << reserr_0[j] << endl;
		}
	      if(i==1)
		{
		  angle_1[j] = angle;
		  angleerr_1[j] = 1;
		  res_1[j] = res;
		  reserr_1[j] = reserror;
		  if(print)	 cout  << angle_1[j] << " " << res_1[j] << " " << reserr_1[j] << endl;
		}
	      if(i==2)
		{
		  angle_2[j] = angle;
		  angleerr_2[j] = 1;
		  res_2[j] = res;
		  reserr_2[j] = reserror;
		  if(print)	 cout  << angle_2[j] << " " << res_2[j] << " " << reserr_2[j] << endl;
		}

	    }//measurements
	}//file scanning
    }//irradiations

  //clsize file for IEEE
  double clsize_0[anglesNONirr];
  double clsize_1[anglesPirr2];
  double clsizeerr_0[anglesNONirr];
  double clsizeerr_1[anglesPirr2];
  double clsize,clsizeerr;
  TString filename2[IEEE];
  
  for(int i =0; i< IEEE; i++)
    {
      filename2[i]      = inputDir+filenames2[i];
      cout << filename2[i] << endl;
      ifstream stream(filename2[i]);
      
      std::string line;
      if(!stream.is_open())
	{
	  cout << " File " << filename2[i] << " not opened" << endl;
	}
      else
	{
	  for(int k =0; k<8;k++)
	    {
	      std::getline(stream,line);
	      if(print) cout << k << " line " << line << endl;
	    }
	  //      while(!stream.eof())
	  for(int j= 0; j < measurements[i]; j++)
	    {
	      stream  >> angle >> clsize >> clsizeerr ;
	      if(print)	 cout << "line " << j+7 << endl;
	      if(i==0)
		{
		  clsize_0[j] = clsize;
		  clsizeerr_0[j] = clsizeerr;
		  if(print)	 cout  << angle_0[j] << " " << clsize_0[j] << " " << clsizeerr_0[j] << endl;
		}
	      if(i==1)
		{
		  clsize_1[j] = clsize;
		  clsizeerr_1[j] = clsizeerr;
		  if(print)	 cout  << angle_1[j] << " " << clsize_1[j] << " " << clsizeerr_1[j] << endl;
		}

	    }//measurements
	}//file scanning
    }//IEEE



  
  if(print) cout << " opening sim file "  << endl;
  
  TString pixelavFile =inputDir+filenames[3];
  cout << pixelavFile << endl;
  ifstream stream(pixelavFile);
      
  std::string line;
  if(!stream.is_open())
    {
      cout << " File " << pixelavFile << " not opened" << endl;
    }
  else
    {
      for(int k =0; k<2;k++)
	{
	  std::getline(stream,line);
	  if(print) cout << k << " line " << line << endl;
	}
      //      while(!stream.eof())
      for(int j= 0; j < measurements[4]; j++)
	{
	  stream  >> angle >> res >> reserror ;
	  if(print)	 cout << "line " << j+2 << endl;

	  angle_4[j] = angle;
	  angleerr_4[j] = 1;
	  res_4[j] = res;
	  reserr_4[j] = reserror;
	  if(print)	 cout  << angle_4[j] << " " << res_4[j] << " " << reserr_4[j] << endl;
	}//measurements
    }//file scanning

  TCanvas *cFDB2 = new TCanvas("cFDB2", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TGraphErrors* resolutionPlot[irradiations+2];
  
  resolutionPlot[0] = new TGraphErrors(anglesNONirr,angle_0,res_0,angleerr_0,reserr_0);
  resolutionPlot[1] = new TGraphErrors(anglesPirr2,angle_1,res_1,angleerr_1,reserr_1);
  resolutionPlot[2] = new TGraphErrors(anglesNirr4,angle_2,res_2,angleerr_2,reserr_2);
  resolutionPlot[3] = new TGraphErrors(1,angle_3,res_3,angleerr_3,reserr_3);
  resolutionPlot[4] = new TGraphErrors(anglesNONirr,angle_4,res_4,angleerr_4,reserr_4);

  resolutionPlot[0]->SetTitle(" ");
  resolutionPlot[0]->GetYaxis()->SetTitle("Resolution [#mum]");
  resolutionPlot[0]->GetXaxis()->SetTitle("Angle [deg]");
  resolutionPlot[0]->SetMarkerSize(2.5);
  resolutionPlot[0]->SetMarkerColor(kRed);
  resolutionPlot[0]->SetLineColor(kRed);
  resolutionPlot[0]->SetMarkerStyle(20);
  resolutionPlot[0]->GetXaxis()->SetLimits(-6.,31.);
  resolutionPlot[0]->GetYaxis()->SetRangeUser(0.,10.);
  resolutionPlot[0]->Draw("AEP");

  resolutionPlot[1]->SetMarkerSize(2.5);
  resolutionPlot[1]->SetMarkerColor(kBlack);
  resolutionPlot[1]->SetMarkerStyle(21);
  resolutionPlot[1]->Draw("EPsame");
				      
  resolutionPlot[2]->SetMarkerSize(2.5);
  resolutionPlot[2]->SetMarkerColor(kBlue);
  resolutionPlot[2]->SetLineColor(kBlue);
  resolutionPlot[2]->SetMarkerStyle(22);
  resolutionPlot[2]->Draw("EPsame");

  resolutionPlot[3]->SetMarkerSize(2.5);
  resolutionPlot[3]->SetMarkerColor(kBlue);
  resolutionPlot[3]->SetLineColor(kBlue);
  resolutionPlot[3]->SetMarkerStyle(23);
  resolutionPlot[3]->Draw("EPsame");

  resolutionPlot[4]->SetMarkerSize(2.5);
  resolutionPlot[4]->SetMarkerColor(kRed);
  resolutionPlot[4]->SetLineColor(kRed);
  resolutionPlot[4]->SetMarkerStyle(24);
  //  resolutionPlot[4]->Draw("EPsame");


  TLegend* legFDB2 = new TLegend(0.65,0.15,0.85,0.35);
  legFDB2->SetLineColor(0);
  legFDB2->SetTextSize(0.027);
  for(int i =0; i < irradiations+1; i++)
    legFDB2->AddEntry(resolutionPlot[i],irr[i],"lp");

  legFDB2->Draw();

  
  TString  outname = outputDir+"ResolutionSummary_angleScans25_"+name;
  cFDB2->SaveAs(outname+".eps");
  cFDB2->SaveAs(outname+".png");
  cFDB2->SaveAs(outname+".pdf");
  cFDB2->SaveAs(outname+".root");
  cFDB2->SaveAs(outname+".C");

  ////// for finn
  TCanvas *c2 = new TCanvas("c2", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  resolutionPlot[0]->Draw("AEP");

  resolutionPlot[1]->Draw("EPsame");
				      
  TLegend* leg2 = new TLegend(0.65,0.15,0.85,0.25);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.027);
  leg2->AddEntry(resolutionPlot[0],irr[0],"lp");
  leg2->AddEntry(resolutionPlot[1],irr[1],"lp");

  leg2->Draw();




  
  outname = outputDir+"ResolutionSummary_angleScans25_freshVSprot_"+name;
  c2->SaveAs(outname+".eps");
  c2->SaveAs(outname+".png");
  c2->SaveAs(outname+".pdf");
  c2->SaveAs(outname+".root");
  c2->SaveAs(outname+".C");



  ////// for IEEE
  TCanvas *c3 = new TCanvas("c3", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillStyle(4000); //will be transparent
  pad2->SetFrameFillStyle(0);
  
  pad1->Draw();
  pad1->cd();
  resolutionPlot[0]->Draw("AEP");

  resolutionPlot[1]->Draw("EPsame");
				      
  pad2->Draw();
  pad2->cd();
  TGraphErrors* clsizePlot[IEEE];
  
  clsizePlot[0] = new TGraphErrors(anglesNONirr,angle_0,clsize_0,angleerr_0,clsizeerr_0);
  clsizePlot[1] = new TGraphErrors(anglesPirr2,angle_1,clsize_1,angleerr_1,clsizeerr_1);

  clsizePlot[0]->SetTitle(" ");
  clsizePlot[0]->GetYaxis()->SetTitle("Average cluster size");
  clsizePlot[0]->GetXaxis()->SetTitle("Angle [deg]");
  clsizePlot[0]->SetMarkerSize(2.5);
  clsizePlot[0]->SetMarkerColor(kRed);
  clsizePlot[0]->SetLineColor(kRed);
  clsizePlot[0]->SetMarkerStyle(24);
  clsizePlot[0]->GetXaxis()->SetLimits(-6.,31.);
  clsizePlot[0]->GetYaxis()->SetRangeUser(0.,10.);
  clsizePlot[0]->Draw("AEPY+");

  clsizePlot[1]->SetMarkerSize(2.5);
  clsizePlot[1]->SetMarkerColor(kBlack);
  clsizePlot[1]->SetMarkerStyle(25);
  clsizePlot[1]->Draw("EPsame");




  TLegend* leg3 = new TLegend(0.5,0.15,0.85,0.25);
  leg3->SetLineColor(0);
  leg3->SetTextSize(0.027);
  leg3->AddEntry(resolutionPlot[0],irr[0],"lp");
  leg3->AddEntry(resolutionPlot[1],irr[1],"lp");

  leg3->Draw();




  
  outname = outputDir+"ResolutionSummary_angleScans25_freshVSprot_withclsize"+name;
  c3->SaveAs(outname+".eps");
  c3->SaveAs(outname+".png");
  c3->SaveAs(outname+".pdf");
  c3->SaveAs(outname+".root");
  c3->SaveAs(outname+".C");


  


  //data vs simulation
  TCanvas *c1 = new TCanvas("c1","3 Graphs",700,900);

  auto *p2 = new TPad("p2","p3",0.,0.,1.,0.3); p2->Draw();
  p2->SetTopMargin(0.001);
  p2->SetBottomMargin(0.3);
  //p2->SetLogx ();
  //p2->SetGrid();

  auto *p1 = new TPad("p1","p1",0.,0.3,1.,1.);  p1->Draw();
  p1->SetBottomMargin(0.001);
  p1->cd();
  //p1->SetGrid();
  //  p1->SetLogx ();
  //p1->SetLogy();


  resolutionPlot[0]->Draw("AEP");
  resolutionPlot[4]->Draw("EPsame");
  
  TLegend* leg1 = new TLegend(0.65,0.1,0.85,0.25);
  leg1->SetLineColor(0);
  leg1->SetTextSize(0.027);
  leg1->AddEntry(resolutionPlot[0],irr[0],"lp");
  leg1->AddEntry(resolutionPlot[4],irr[4],"lp");

  leg1->Draw();
  

  // ratio
  p2->cd();
  TGraph*r = new TGraph(anglesNONirr); r->SetTitle("");
  for (int i=0; i<anglesNONirr; i++) r->SetPoint(i, angle_0[i], res_0[i]/res_4[i]);
  r->GetXaxis()->SetLabelSize(0.075);
  r->GetYaxis()->SetLabelSize(0.075);
  r->GetYaxis()->SetTitle("data/MC");
  r->Draw("AL");
  
  outname = outputDir+"ResolutionSummary_angleScans25_dataVSpixelav_"+name;
  c1->SaveAs(outname+".eps");
  c1->SaveAs(outname+".png");
  c1->SaveAs(outname+".pdf");
  c1->SaveAs(outname+".root");

  

  
}//resolution 



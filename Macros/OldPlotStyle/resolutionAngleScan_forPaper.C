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

#include "tdrstyle.C"
#include "CMS_lumi.C"


bool print=true;
using namespace std;
#define anglesNONirr 29
#define anglesPirr2 28
#define anglesNirr4 9
#define irradiations 3
#define biasScan 4


void TDR();
void TDR2(TCanvas * c_all);

void resolutionAngleScan_forPaper(TString thr="500",TString name = "preliminary", TString func = "RMS", bool unfolding = false)
{



  // 194i  4 25   800   12   5.6   3.86636   0.0322234   3839   FDB150P_12_R4S100x25-P1_1   Pstop_default_FDB   15 

  
  TString inputDir="/home/zoiirene/Output/TextFiles/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  int measurements[irradiations+1];
  measurements[0] = anglesNONirr;
  measurements[1] = anglesPirr2;
  measurements[2] = anglesNirr4;
  measurements[3] = 1;


  //  Ascan_resTree_148_dphcutB22_RMS.txt
  TString filenames[irradiations];
  filenames[0] = "Ascan_res_148_RMS.txt";
  filenames[1] = "Ascan_res_120i_RMS.txt";
  filenames[2] = "Ascan_194i_RMS_800_res.txt";


  TString filenames2[irradiations+1];
  filenames2[0] = "Ascan_clsizeB_148_RMS.txt";
  filenames2[1] = "Ascan_clsizeB_120i_RMS.txt";
  filenames2[2] = "Ascan_clsizeB_194i_800.txt"; //5.2 GeV

  TString ss_bias[biasScan] = {"200","400","600","800"}; 
  TString filenamesNbias[biasScan];
  for (int i = 0; i< biasScan; i++){
    filenamesNbias[i] = "Ascan_194i_RMS_"+ss_bias[i]+"_res.txt";
    if(print) cout << "files for neutron bias scan " << filenamesNbias[i] << endl;
  }
  
  TString irr[irradiations+1];
  irr[0] = " Non-irradiated, 5.6 GeV"; // "no irr, 5.6 GeV";
  //  irr[1] = "proton irr at #phi_{eq}=2#times10^{15} cm^{-2}, 5.6 GeV";
  irr[1] = "#phi_{eq} = 2.1#times10^{15} cm^{-2}, proton, 5.6 GeV";
  irr[2] = "#phi_{eq} = 4#times10^{15} cm^{-2}, neutron, 5.2 GeV";
  irr[3] = "#phi_{eq} = 4#times10^{15} cm^{-2}, neutron, 5.6 GeV";
  
  double angle_0[anglesNONirr];
  double angle_1[anglesPirr2];
  double angle_2[anglesNirr4];
  double angle_3[] = {12};

  double angleerr_0[anglesNONirr];
  double angleerr_1[anglesPirr2];
  double angleerr_2[anglesNirr4];
  double angleerr_3[]= {1};

  double res_0[anglesNONirr];
  double res_1[anglesPirr2];
  double res_2[anglesNirr4];
  double res_3[] = {3.87};

  double reserr_0[anglesNONirr];
  double reserr_1[anglesPirr2];
  double reserr_2[anglesNirr4];
  double reserr_3[] = {0.03};

  double res_N[biasScan][anglesNirr4];
  double reserr_N[biasScan][anglesNirr4];
  
  double angle;
  double res,reserror,angleerr;
  int lines[irradiations] = {8,8,7};
  
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
	  for(int k =0; k<lines[i];k++)
	    {
	      std::getline(stream,line);
	      if(print) cout << k << " line " << line << endl;
	    }
	  //      while(!stream.eof())
	  for(int j= 0; j < measurements[i]; j++)
	    {
	      stream  >> angle >> res >> reserror ;
	      if(print)	 cout << "line " << j+lines[i]-1 << endl;
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

  // neutron angle & bias scan
  TString filenameN[biasScan];
  for(int i =0; i< biasScan; i++)
    {
      filenameN[i]      = inputDir+filenamesNbias[i];
      cout << filenameN[i] << endl;
      ifstream stream(filenameN[i]);
      
      std::string line;
      if(!stream.is_open())
	{
	  cout << " File " << filenameN[i] << " not opened" << endl;
	}
      else
	{
	  for(int k =0; k<lines[2];k++)
	    {
	      std::getline(stream,line);
	      if(print) cout << k << " line " << line << endl;
	    }
	  //      while(!stream.eof())
	  for(int j= 0; j < anglesNirr4; j++)
	    {
	      stream  >> angle >> res >> reserror ;
	      if(print)	 cout << "line " << j+7 << endl;
	      angle_2[j] = angle;
	      angleerr_2[j] = 1;
	      res_N[i][j] = res;
	      reserr_N[i][j] = reserror;
	      if(print)	 cout  << angle_2[j] << " " << res_N[i][j] << " " << reserr_N[i][j] << " at " << ss_bias[i] << endl;
	      
	    }//angle measurements
	}//file scanning
    }//bias scan






   //clsize file for IEEE
  double clsize_0[anglesNONirr];
  double clsize_2[anglesNONirr];
  double clsize_1[anglesPirr2];
  double clsizeerr_0[anglesNONirr];
  double clsizeerr_2[anglesNONirr];
  double clsizeerr_1[anglesPirr2];

  double clsize_3[] = {2.002}; //mean of nrowB of straightTracksY_isoAandCandB_straightTracksX (same as after cut on charge), drei-r3839_irene_dphcutB15.root neutron irr, 5.6 GeV
  double clsizeerr_3[] = {0.003}; //mean error
  double clsize_4[anglesNirr4]; // neutron irr 5.2 GeV
  double clsizeerr_4[anglesNirr4]; //mean error
  
  double clsize,clsizeerr;
  TString filename2[irradiations];
  int lines2[irradiations] = {8,8,9};
  
  for(int i =0; i< irradiations; i++)
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
	  for(int k =0; k<lines2[i];k++)
	    {
	      std::getline(stream,line);
	      if(print) cout << k << " line " << line << endl;
	    }
	  //      while(!stream.eof())
	  for(int j= 0; j < measurements[i]; j++)
	    {
	      stream  >> angle >> clsize >> clsizeerr ;
	      if(print)	 cout << "line " << j+lines2[i]-1 << endl;
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
	      if(i==2)
		{
		  clsize_4[j] = clsize;
		  clsizeerr_4[j] = clsizeerr;
		  if(print)	 cout  << angle_2[j] << " " << clsize_4[j] << " " << clsizeerr_4[j] << endl;
		}

	    }//measurements
	}//file scanning
    }//IEEE



  /////// plots!

  
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


  resolutionPlot[0]->SetTitle(" ");
  resolutionPlot[0]->GetYaxis()->SetTitle("Single hit resolution [#mum]");
  resolutionPlot[0]->GetXaxis()->SetTitle("#theta [deg]");
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



  TLegend* legFDB2 = new TLegend(0.65,0.15,0.85,0.35);
  legFDB2->SetLineColor(0);
  legFDB2->SetTextSize(0.035);
  for(int i =0; i < irradiations+1; i++)
    legFDB2->AddEntry(resolutionPlot[i],irr[i],"lp");

  legFDB2->Draw();

  
  TString  outname = outputDir+"ResolutionSummaryPaper_angleScans25_"+name;
  cFDB2->SaveAs(outname+".eps");
  cFDB2->SaveAs(outname+".png");
  cFDB2->SaveAs(outname+".pdf");
  cFDB2->SaveAs(outname+".root");
  cFDB2->SaveAs(outname+".C");


  // neutron angle & bias scan
  
  TCanvas *cN = new TCanvas("cN", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TGraphErrors* resolutionPlotN[biasScan];
  for (int i = 0; i< biasScan; i++)
    {
      resolutionPlotN[i] = new TGraphErrors(anglesNirr4,angle_2,res_N[i],angleerr_2,reserr_N[i]);
      resolutionPlotN[i]->SetMarkerSize(2.);
      //      resolutionPlotN[i]->SetLineStyle(2);
    }


  resolutionPlotN[0]->SetTitle(" ");
  resolutionPlotN[0]->GetYaxis()->SetTitle("Single hit resolution [#mum]");
  resolutionPlotN[0]->GetXaxis()->SetTitle("#theta [deg]");

  resolutionPlotN[0]->SetMarkerColor(kCyan-3);
  resolutionPlotN[0]->SetLineColor(kCyan-3);
  resolutionPlotN[0]->SetMarkerStyle(20);
  resolutionPlotN[0]->GetXaxis()->SetLimits(0.,25.);
  resolutionPlotN[0]->GetYaxis()->SetRangeUser(0.,7.);
  //resolutionPlotN[0]->GetYaxis()->SetRangeUser(4.,7.);
  resolutionPlotN[0]->Draw("AEPC");

  resolutionPlotN[1]->SetMarkerColor(kCyan+2);
  resolutionPlotN[1]->SetLineColor(kCyan+2);
  resolutionPlotN[1]->SetMarkerStyle(21);
  resolutionPlotN[1]->Draw("EPCsame");
				      
  resolutionPlotN[2]->SetMarkerColor(kBlue-9);
  resolutionPlotN[2]->SetLineColor(kBlue-9);
  resolutionPlotN[2]->SetMarkerStyle(23);
  resolutionPlotN[2]->Draw("EPCsame");

  resolutionPlotN[3]->SetMarkerColor(kBlue);
  resolutionPlotN[3]->SetLineColor(kBlue);
  resolutionPlotN[3]->SetMarkerStyle(22);
  resolutionPlotN[3]->Draw("EPCsame");



  TLegend* legN = new TLegend(0.12,0.25,0.27,0.4);
  legN->SetLineColor(0);
  legN->SetTextSize(0.035);
  for(int i =0; i < irradiations+1; i++)
    legN->AddEntry(resolutionPlotN[i],ss_bias[i]+" V","lp");

  legN->Draw();

  
  outname = outputDir+"ResolutionSummaryPaper_angleBiasScans25Neutron_"+name;
  cN->SaveAs(outname+".eps");
  cN->SaveAs(outname+".png");
  cN->SaveAs(outname+".pdf");
  cN->SaveAs(outname+".root");
  cN->SaveAs(outname+".C");


  ////// for finn
  //  setTDRStyle();
  TDR();

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
  resolutionPlot[2]->Draw("EPsame");
  resolutionPlot[3]->Draw("EPsame");
			      
  TLegend* leg2 = new TLegend(0.65,0.15,0.85,0.25);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.035);
  leg2->AddEntry(resolutionPlot[0],irr[0],"lp");
  leg2->AddEntry(resolutionPlot[1],irr[1],"lp");

  leg2->Draw();




  TDR2(c2);
  outname = outputDir+"ResolutionSummaryPaper_angleScans25_freshVSprotVSneut_"+name;
  c2->SaveAs(outname+".eps");
  c2->SaveAs(outname+".png");
  c2->SaveAs(outname+".pdf");
  c2->SaveAs(outname+".root");
  c2->SaveAs(outname+".C");



  ////// for IEEE
  //  TDR();
  TCanvas *c3 = new TCanvas("c3", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  c3->SetLeftMargin(1.);
  c3->SetTopMargin(1.);
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
  resolutionPlot[2]->Draw("EPsame");
  resolutionPlot[3]->Draw("EPsame");
				      
  pad2->Draw();
  pad2->cd();
  TGraphErrors* clsizePlot[irradiations+1];
  
  clsizePlot[0] = new TGraphErrors(anglesNONirr,angle_0,clsize_0,angleerr_0,clsizeerr_0);
  clsizePlot[1] = new TGraphErrors(anglesPirr2,angle_1,clsize_1,angleerr_1,clsizeerr_1);
  clsizePlot[2] = new TGraphErrors(anglesNirr4,angle_2,clsize_4,angleerr_2,clsizeerr_4);
  clsizePlot[3] = new TGraphErrors(1,angle_3,clsize_3,angleerr_3,clsizeerr_3);


  
  clsizePlot[0]->SetTitle(" ");
  clsizePlot[0]->GetYaxis()->SetTitle("Average cluster size");
  clsizePlot[0]->GetXaxis()->SetTitle("#theta [deg]");
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

  clsizePlot[2]->SetMarkerSize(2.5);
  clsizePlot[2]->SetMarkerColor(kBlue);
  clsizePlot[2]->SetMarkerStyle(26);
  clsizePlot[2]->Draw("EPsame");

  clsizePlot[3]->SetMarkerSize(2.5);
  clsizePlot[3]->SetMarkerColor(kBlue);
  clsizePlot[3]->SetMarkerStyle(32);
  clsizePlot[3]->Draw("EPsame");

  TLine *  linea = new TLine( -6.,2.,31.,2.);
  linea->SetLineColor(kGray);
  linea->SetLineWidth(2);
  linea->SetLineStyle(2);
  linea->Draw("same");

  TLine *  lineb = new TLine( 9.5,0.,9.5,6.5);
  lineb->SetLineColor(kGray);
  lineb->SetLineWidth(2);
  lineb->SetLineStyle(2);
  lineb->Draw("same");


  TLegend* leg3 = new TLegend(0.3,0.65,0.65,0.85);
  leg3->SetLineColor(0);
  leg3->SetTextSize(0.035);
  //  leg3->AddEntry(resolutionPlot[0],"beam momentum 5.6 GeV","");
  leg3->AddEntry(resolutionPlot[0],irr[0],"lp");
  leg3->AddEntry(resolutionPlot[1],irr[1],"lp");
  leg3->AddEntry(resolutionPlot[2],irr[2],"lp");
  leg3->AddEntry(resolutionPlot[3],irr[3],"lp");

  leg3->Draw();



  pad1->cd(); 
  //  TDR2(c3);
  outname = outputDir+"ResolutionSummaryPaper_angleScans25_freshVSprotVSneutr_withclsize_line2_newlabel"+name;
  c3->SaveAs(outname+".eps");
  c3->SaveAs(outname+".png");
  c3->SaveAs(outname+".pdf");
  c3->SaveAs(outname+".root");
  c3->SaveAs(outname+".C");


  
}//resolution 


void TDR()
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  extraText  = "Work in progress";  // default extra text is "Preliminary"
  lumi_sqrtS = "";
  bool drawLogo = true;
  int H_ref = 600;
  int W_ref = 600;
  float T = 0.07*H_ref;//0.07
  float B = 0.11*H_ref;//0.12
  float L = 0.12*W_ref;
  float R = 0.01*W_ref;


}

void TDR2(TCanvas * c_all)
{
  CMS_lumi( c_all, 0, 11);
  //  CMS_lumi( c_all, 0, 11);
  c_all->Update();
  c_all->RedrawAxis();
}

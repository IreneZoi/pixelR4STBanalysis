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
#define dphcuts 7
#define comparison 2

void TDR();
void TDR2(TCanvas * c_all, int period=0, int pos= 11);

void resolutionAngleScanNonIrrThrRatio(TString name = "preliminary_ThrScanProtIrrABCratio85", TString func = "RMS95", bool unfolding = false){
  
  TString inputDir="/home/zoiirene/Output/TextFiles/";
  TString outputDir="/home/zoiirene/Output/Plots/";
  int i_dphcut[dphcuts] = {10, 12, 15, 18, 22, 30, 40};
  TString ss_dphcut[dphcuts] = {"10", "12", "15", "18", "22", "30", "40"};
  
  TString filenames[dphcuts][comparison];
  TString filenames2[dphcuts][comparison];
  for(int i = 0; i < dphcuts; i++){
    filenames[i][0] = "Ascan_resTreeCorr_120i_dphcutB"+ss_dphcut[i]+"_RMS95_testThr.txt";
    filenames2[i][0] = "Ascan_clsizeB_120i_dphcutB"+ss_dphcut[i]+"_RMS95_testThr.txt";
    filenames[i][1] = "Ascan_resTreeCorr_120i_dphcutB"+ss_dphcut[i]+"_RMS95_testThr85.txt";
    filenames2[i][1] = "Ascan_clsizeB_120i_dphcutB"+ss_dphcut[i]+"_RMS95_testThr85.txt";
  }
  
  TString irr =  " Non-irradiated, 5.6 GeV"; // "no irr, 5.6 GeV";
  double angle_0[dphcuts][anglesNONirr];
  double angleerr_0[dphcuts][anglesNONirr];
  double res_0_ABC[dphcuts][anglesNONirr];
  double reserr_0_ABC[dphcuts][anglesNONirr];
  double res_0_B[dphcuts][anglesNONirr];
  double reserr_0_B[dphcuts][anglesNONirr];
  double angle;
  double res,reserror,angleerr;
  int lines = 8;
  
  TString filename;
  for(int i =0; i< dphcuts; i++){

    filename      = inputDir+filenames[i][0];
    cout << filename << endl;
    ifstream stream(filename);
      
    std::string line;
    if(!stream.is_open()){
      cout << " File " << filename[i] << " not opened" << endl;
    }else{
      for(int k =0; k<lines;k++)  {
	std::getline(stream,line);
	if(print) cout << k << " line " << line << endl;
      }

      for(int j= 0; j < anglesNONirr; j++){
	stream  >> angle >> res >> reserror ;
	if(print)	 cout << "line " << j+lines-1 << endl;
	angle_0[i][j] = angle;
	angleerr_0[i][j] = 1;
	res_0_ABC[i][j] = res;
	reserr_0_ABC[i][j] = reserror;
	if(print)	 cout  << angle_0[i][j] << " " << res_0_ABC[i][j] << " " << reserr_0_ABC[i][j] << endl;
      }//measurements
    }//file scanning
  }//irradiations

  
  for(int i =0; i< dphcuts; i++){

    filename      = inputDir+filenames[i][1];
    cout << filename << endl;
    ifstream stream(filename);
      
    std::string line;
    if(!stream.is_open()){
      cout << " File " << filename[i] << " not opened" << endl;
    }else{
      for(int k =0; k<lines;k++)  {
	std::getline(stream,line);
	if(print) cout << k << " line " << line << endl;
      }

      for(int j= 0; j < anglesNONirr; j++){
	stream  >> angle >> res >> reserror ;
	if(print)	 cout << "line " << j+lines-1 << endl;
	res_0_B[i][j] = res;
	reserr_0_B[i][j] = reserror;
	if(print)	 cout  << angle<< " " << res_0_B[i][j] << " " << reserr_0_B[i][j] << endl;
      }//measurements
    }//file scanning
  }//irradiations


  double ratio[dphcuts][anglesNONirr];
  double ratioerr[dphcuts][anglesNONirr];

  for(int i =0; i< dphcuts; i++){
    cout << " dphcut " <<    ss_dphcut[i]  << endl;
    for(int j= 0; j < anglesNONirr; j++){
      cout << " Angle " <<    angle_0[i][j]  << endl;
      ratio[i][j] = res_0_ABC[i][j] / res_0_B[i][j];
      cout << ratio[i][j] << endl;
    }
  }
  /*
   //clsize file for IEEE
  double clsize_0[dphcuts][anglesNONirr];
  double clsizeerr_0[dphcuts][anglesNONirr];
  double clsize,clsizeerr;
  
  for(int i =0; i< dphcuts; i++){
      filename  = inputDir+filenames2[i];
      cout << filename << endl;
      ifstream stream(filename);
      
      std::string line;
      if(!stream.is_open())	{
	cout << " File " << filename << " not opened" << endl;
      }else{
	for(int k =0; k<lines;k++){
	  std::getline(stream,line);
	  if(print) cout << k << " line " << line << endl;
	}
	
	for(int j= 0; j <anglesNONirr; j++){
	  stream  >> angle >> clsize >> clsizeerr ;
	  if(print)	 cout << "line " << j+lines-1 << endl;
	  clsize_0[i][j] = clsize;
	  clsizeerr_0[i][j] = clsizeerr;
	  if(print)	 cout  << angle_0[i][j] << " " << clsize_0[i][j] << " " << clsizeerr_0[i][j] << endl;
	  
	}//measurements
      }//file scanning
  }//IEEE

  */


  

  /////// plots!

  
  TCanvas *cFDB2 = new TCanvas("cFDB2", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  bool logar = false;
  TGraph* resolutionPlot[dphcuts];
  for(int i =0; i< dphcuts; i++){
    resolutionPlot[i] = new TGraph(anglesNONirr,angle_0[i],ratio[i]); 
    resolutionPlot[i]->SetTitle(" ");
    resolutionPlot[i]->GetYaxis()->SetTitle("ratio 90/85");
    resolutionPlot[i]->GetXaxis()->SetTitle("#theta [deg]");
    resolutionPlot[i]->SetMarkerSize(2.5);
    /*
    resolutionPlot[i]->SetMarkerColor(kRed+i-10);
    resolutionPlot[i]->SetLineColor(kRed+i-10);
    if(i>15){
      resolutionPlot[i]->SetMarkerColor(kMagenta+4+i-dphcuts);
      resolutionPlot[i]->SetLineColor(kMagenta+4+i-dphcuts);
    }
    */
      resolutionPlot[i]->SetMarkerColor(616+i);
      resolutionPlot[i]->SetLineColor(616+i);

    if(i==5){
      resolutionPlot[i]->SetMarkerColor(615);
      resolutionPlot[i]->SetLineColor(615);
      
    }
    resolutionPlot[i]->SetMarkerStyle(20+i);
    resolutionPlot[i]->GetXaxis()->SetLimits(-6.,31.);
    resolutionPlot[i]->GetYaxis()->SetRangeUser(0.8,1.2);
    if(logar) resolutionPlot[i]->GetYaxis()->SetRangeUser(0.05,50.);
    
    if(i==0)    resolutionPlot[i]->Draw("AEP");
    else resolutionPlot[i]->Draw("EPsame");
  }


  if(logar) cFDB2->SetLogy();   
  TLegend* legFDB2 = new TLegend(0.5,0.7,0.85,0.85);
  legFDB2->SetLineColor(0);
  legFDB2->SetTextSize(0.035);
  legFDB2->SetNColumns(5);
  for(int i =0; i < dphcuts; i++)
    legFDB2->AddEntry(resolutionPlot[i],ss_dphcut[i],"lp");

  legFDB2->Draw();

  TLine *  linea = new TLine( -6.,1.,31.,1.);
  linea->SetLineColor(kGray);
  linea->SetLineWidth(2);
  linea->SetLineStyle(2);
  linea->Draw("same");
  
  
  TString  outname = outputDir+"Ratio_angleScans25_"+name;
  if(logar) outname+="_log";
  cFDB2->SaveAs(outname+".eps");
  cFDB2->SaveAs(outname+".png");
  cFDB2->SaveAs(outname+".pdf");
  cFDB2->SaveAs(outname+".root");
  cFDB2->SaveAs(outname+".C");


  /*
  //cluster size

  TCanvas *cCL = new TCanvas("cCL", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  logar = false;
  TGraphErrors* clsizeolutionPlot[dphcuts];
  for(int i =0; i< dphcuts; i++){
    clsizeolutionPlot[i] = new TGraphErrors(anglesNONirr,angle_0[i],clsize_0[i],angleerr_0[i],clsizeerr_0[i]);
    clsizeolutionPlot[i]->SetTitle(" ");
    clsizeolutionPlot[i]->GetYaxis()->SetTitle("Average cluster size");
    clsizeolutionPlot[i]->GetXaxis()->SetTitle("#theta [deg]");
    clsizeolutionPlot[i]->SetMarkerSize(2.5);
      clsizeolutionPlot[i]->SetMarkerColor(616+i);
      clsizeolutionPlot[i]->SetLineColor(616+i);

    if(i==5){
      clsizeolutionPlot[i]->SetMarkerColor(615);
      clsizeolutionPlot[i]->SetLineColor(615);
      
    }
    clsizeolutionPlot[i]->SetMarkerStyle(20+i);
    clsizeolutionPlot[i]->GetXaxis()->SetLimits(-6.,31.);
    clsizeolutionPlot[i]->GetYaxis()->SetRangeUser(0.,6.);
    if(logar) clsizeolutionPlot[i]->GetYaxis()->SetRangeUser(0.05,50.);
    
    if(i==0)    clsizeolutionPlot[i]->Draw("AEP");
    else clsizeolutionPlot[i]->Draw("EPsame");
  }


  if(logar) cCL->SetLogy();   
  TLegend* legCL = new TLegend(0.15,0.7,0.5,0.85);
  legCL->SetLineColor(0);
  legCL->SetTextSize(0.035);
  legCL->SetNColumns(5);
  for(int i =0; i < dphcuts; i++)
    legCL->AddEntry(clsizeolutionPlot[i],ss_dphcut[i],"lp");

  legCL->Draw();

  TLine *  linea = new TLine( 0.,2.,31.,2.);
  linea->SetLineColor(kGray);
  linea->SetLineWidth(2);
  linea->SetLineStyle(2);
  linea->Draw("same");

  
  outname = outputDir+"Clsize_angleScans25_"+name;
  if(logar) outname+="_log";
  cCL->SaveAs(outname+".eps");
  cCL->SaveAs(outname+".png");
  cCL->SaveAs(outname+".pdf");
  cCL->SaveAs(outname+".root");
  cCL->SaveAs(outname+".C");
  */
  

  /*
  ////// for paper
  //   TDR();
  TCanvas *c3 = new TCanvas("c3", "FDB resolution", 600, 600);
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
  resolutionPlot[0]->GetXaxis()->SetLimits(0.,31.);
  resolutionPlot[0]->SetMarkerSize(1.);
  resolutionPlot[0]->Draw("AEP");
  resolutionPlot[1]->SetMarkerSize(1.);
  resolutionPlot[1]->Draw("EPsame");
  resolutionPlot[2]->SetMarkerSize(1.);
  resolutionPlot[2]->Draw("EPsame");
  resolutionPlot[3]->SetMarkerSize(1.);
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
  clsizePlot[0]->SetMarkerSize(1.);
  clsizePlot[0]->GetXaxis()->SetLimits(0.,31.);
  clsizePlot[0]->GetYaxis()->SetRangeUser(0.,10.);
  clsizePlot[0]->Draw("AEPY+");

  clsizePlot[1]->SetMarkerSize(2.5);
  clsizePlot[1]->SetMarkerColor(kBlack);
  clsizePlot[1]->SetMarkerStyle(25);
  clsizePlot[1]->SetMarkerSize(1.);
  clsizePlot[1]->Draw("EPsame");

  clsizePlot[2]->SetMarkerSize(2.5);
  clsizePlot[2]->SetMarkerColor(kBlue);
  clsizePlot[2]->SetMarkerStyle(26);
  clsizePlot[2]->SetMarkerSize(1.);  
  clsizePlot[2]->Draw("EPsame");

  clsizePlot[3]->SetMarkerSize(2.5);
  clsizePlot[3]->SetMarkerColor(kBlue);
  clsizePlot[3]->SetMarkerStyle(32);
  clsizePlot[3]->SetMarkerSize(1.);
  clsizePlot[3]->Draw("EPsame");

  TLine *  linea = new TLine( 0.,2.,31.,2.);
  linea->SetLineColor(kGray);
  linea->SetLineWidth(2);
  linea->SetLineStyle(2);
  linea->Draw("same");

  TLine *  lineb = new TLine( 9.5,0.,9.5,6.5);
  lineb->SetLineColor(kGray);
  lineb->SetLineWidth(2);
  lineb->SetLineStyle(2);
  lineb->Draw("same");


  TLegend* leg3 = new TLegend(0.2,0.67,0.6,0.87);
  leg3->SetLineColor(0);
  leg3->SetTextSize(0.035);
  //  leg3->AddEntry(resolutionPlot[0],"beam momentum 5.6 GeV","");
  leg3->AddEntry(resolutionPlot[0],irr[0],"lp");
  leg3->AddEntry(resolutionPlot[1],irr[1],"lp");
  leg3->AddEntry(resolutionPlot[2],irr[2],"lp");
  leg3->AddEntry(resolutionPlot[3],irr[3],"lp");

  leg3->Draw();



  pad1->cd(); 
  //  TDR2(c3,0,0);
  outname = outputDir+"ResolutionSummaryPaper_angleScans25_noNeg_freshVSprotVSneutr_withclsize_line2_newlabel"+name;
  c3->SaveAs(outname+".eps");
  c3->SaveAs(outname+".png");
  c3->SaveAs(outname+".pdf");
  c3->SaveAs(outname+".root");
  c3->SaveAs(outname+".C");

  */




  
  /*
  ////// for paper
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
*/














  
}//resolution 


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

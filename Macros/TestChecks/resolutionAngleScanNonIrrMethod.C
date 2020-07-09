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
#define Cuts 7

void TDR();
void TDR2(TCanvas * c_all, int period=0, int pos= 11);

void resolutionAngleScanNonIrrMethod(TString name = "preliminary_MethodScanNonIrrABC1clusterClosestBugFix", TString func = "RMS95", bool unfolding = false){
  
  TString inputDir="/home/zoiirene/Output/TextFiles/";
  TString outputDir="/home/zoiirene/Output/Plots/";
  
  TString filenames[Cuts];

  //  TString interval[Cuts]={"landau90RMS","landau90RMS95","landauAllRMS","landauAllRMS95"};
  TString interval[Cuts]={"RMS99_1clusterABC","RMS99_closest","RMSself6sig_closest","RMSself5sig_closest","RMSself4sig_closest","RMS99_closest99p9","RMS99_closest99"}; 
  for(int i = 0; i < Cuts; i++){
    
    //    filenames[i] = "Ascan_resTree_148_dphcutB22_RMS95_testThrABC_"+interval[i]+".txt";
    filenames[i] = "Ascan_resTree_148_dphcutB22_"+interval[i]+".txt";
  }
  
  TString irr =  " Non-irradiated, 5.6 GeV"; // "no irr, 5.6 GeV";
  double angle_0[Cuts][anglesNONirr];
  double angleerr_0[Cuts][anglesNONirr];
  double res_0[Cuts][anglesNONirr];
  double reserr_0[Cuts][anglesNONirr];
  double angle;
  double res,reserror,angleerr;
  int lines = 8;
  
  TString filename;
  for(int i =0; i< Cuts; i++){

    filename      = inputDir+filenames[i];
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
	res_0[i][j] = res;
	reserr_0[i][j] = reserror;
	if(print)	 cout  << angle_0[i][j] << " " << res_0[i][j] << " " << reserr_0[i][j] << endl;
      }//measurements
    }//file scanning
  }//irradiations



  for(int i = 0; i < Cuts; i++){
    
    //    filenames[i] = "Ascan_resTree_148_dphcutB22_RMS95_testThrABC_"+interval[i]+".txt";
    filenames[i] = "Ascan_Entries_148_dphcutB22_"+interval[i]+".txt";
  }
  

  double entries_0[Cuts][anglesNONirr];
  double entries;
  
  for(int i =0; i< Cuts; i++){

    filename      = inputDir+filenames[i];
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
	stream  >> angle >> entries;
	if(print)	 cout << "line " << j+lines-1 << endl;
	entries_0[i][j] = entries;
	if(print)	 cout  << angle << " " << entries_0[i][j] << endl;
      }//measurements
    }//file scanning
  }//irradiations


  for(int i = 2; i < 5; i++){

    //    filenames[i] = "Ascan_resTree_148_dphcutB22_RMS95_testThrABC_"+interval[i]+".txt";
    filenames[i] = "Ascan_Percentage_148_dphcutB22_"+interval[i]+".txt";
  }


  double percent_0[3][anglesNONirr];
  double percent;

  for(int i =2; i< 5; i++){

    filename      = inputDir+filenames[i];
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
	stream  >> angle >> percent;
	if(print)        cout << "line " << j+lines-1 << endl;
	percent_0[i][j] = percent;
	if(print)        cout  << angle << " " << percent_0[i][j] << endl;
      }//measurements
    }//file scanning
  }//irradiations
  
  

  /////// plots!

  
  TCanvas *cFDB2 = new TCanvas("cFDB2", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TGraphErrors* resolutionPlot[Cuts];
  for(int i =0; i< Cuts; i++){
    resolutionPlot[i] = new TGraphErrors(anglesNONirr,angle_0[i],res_0[i],angleerr_0[i],reserr_0[i]);
    resolutionPlot[i]->SetTitle(" ");
    resolutionPlot[i]->GetYaxis()->SetTitle("Single hit resolution [#mum]");
    resolutionPlot[i]->GetXaxis()->SetTitle("#theta [deg]");
    resolutionPlot[i]->SetMarkerSize(2.5);
    /*
    resolutionPlot[i]->SetMarkerColor(kRed+i-10);
    resolutionPlot[i]->SetLineColor(kRed+i-10);
    if(i>15){
      resolutionPlot[i]->SetMarkerColor(kMagenta+4+i-Cuts);
      resolutionPlot[i]->SetLineColor(kMagenta+4+i-Cuts);
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
    resolutionPlot[i]->GetYaxis()->SetRangeUser(0.,10.);
    
    if(i==0)    resolutionPlot[i]->Draw("AEP");
    else resolutionPlot[i]->Draw("EPsame");
  }


  //  if(logar) cFDB2->SetLogy();   
  TLegend* legFDB2 = new TLegend(0.2,0.7,0.7,0.85);
  legFDB2->SetLineColor(0);
  legFDB2->SetTextSize(0.035);
  legFDB2->SetNColumns(2);
  for(int i =0; i < Cuts; i++){
    //if(i==0) legFDB2->AddEntry(resolutionPlot[i],"all clusters RMS95","lp");
    //else
    legFDB2->AddEntry(resolutionPlot[i],interval[i],"lp");
  }
  legFDB2->Draw();

  
  TString  outname = outputDir+"Resolution_angleScans25_"+name;
  cFDB2->SaveAs(outname+".eps");
  cFDB2->SaveAs(outname+".png");
  cFDB2->SaveAs(outname+".pdf");
  cFDB2->SaveAs(outname+".root");
  cFDB2->SaveAs(outname+".C");



  // entries


  TCanvas *cEntries2 = new TCanvas("cEntries2", "Entries resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TGraph* EntriesPlot[Cuts];
  for(int i =0; i< Cuts; i++){
    EntriesPlot[i] = new TGraph(anglesNONirr,angle_0[i],entries_0[i]);
    EntriesPlot[i]->SetTitle(" ");
    EntriesPlot[i]->GetYaxis()->SetTitle("Entries");
    EntriesPlot[i]->GetXaxis()->SetTitle("#theta [deg]");
    EntriesPlot[i]->SetMarkerSize(2.5);
    /*
    EntriesPlot[i]->SetMarkerColor(kRed+i-10);
    EntriesPlot[i]->SetLineColor(kRed+i-10);
    if(i>15){
      EntriesPlot[i]->SetMarkerColor(kMagenta+4+i-Cuts);
      EntriesPlot[i]->SetLineColor(kMagenta+4+i-Cuts);
    }
    */
      EntriesPlot[i]->SetMarkerColor(616+i);
      EntriesPlot[i]->SetLineColor(616+i);

    if(i==5){
      EntriesPlot[i]->SetMarkerColor(615);
      EntriesPlot[i]->SetLineColor(615);
      
    }
    EntriesPlot[i]->SetMarkerStyle(20+i);
    EntriesPlot[i]->GetXaxis()->SetLimits(-6.,31.);
    EntriesPlot[i]->GetYaxis()->SetRangeUser(0.,36000.);
    
    if(i==0)    EntriesPlot[i]->Draw("AP");
    else EntriesPlot[i]->Draw("Psame");
  }


  //  if(logar) cEntries2->SetLogy();   
  TLegend* legEntries2 = new TLegend(0.2,0.3,0.7,0.4);
  legEntries2->SetLineColor(0);
  legEntries2->SetTextSize(0.035);
  legEntries2->SetNColumns(2);
  for(int i =0; i < Cuts; i++){
    //if(i==0) legEntries2->AddEntry(EntriesPlot[i],"all clusters RMS95","lp");
    //else
    legEntries2->AddEntry(EntriesPlot[i],interval[i],"lp");
  }
  legEntries2->Draw();

  
  outname = outputDir+"Entries_angleScans25_"+name;
  cEntries2->SaveAs(outname+".eps");
  cEntries2->SaveAs(outname+".png");
  cEntries2->SaveAs(outname+".pdf");
  cEntries2->SaveAs(outname+".root");
  cEntries2->SaveAs(outname+".C");

  // entries


  TCanvas *cPercent = new TCanvas("cPercent", "Percent resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TGraph* PercentPlot[3];
  for(int i =2; i< 5; i++){
    PercentPlot[i-2] = new TGraph(anglesNONirr,angle_0[i],percent_0[i]);
    PercentPlot[i-2]->SetTitle(" ");
    PercentPlot[i-2]->GetYaxis()->SetTitle("Percent");
    PercentPlot[i-2]->GetXaxis()->SetTitle("#theta [deg]");
    PercentPlot[i-2]->SetMarkerSize(2.5);
    /*
    PercentPlot[i-2]->SetMarkerColor(kRed+i-10);
    PercentPlot[i-2]->SetLineColor(kRed+i-10);
    if(i>15){
      PercentPlot[i-2]->SetMarkerColor(kMagenta+4+i-Cuts);
      PercentPlot[i-2]->SetLineColor(kMagenta+4+i-Cuts);
    }
    */
      PercentPlot[i-2]->SetMarkerColor(616+i);
      PercentPlot[i-2]->SetLineColor(616+i);

    if(i==5){
      PercentPlot[i-2]->SetMarkerColor(615);
      PercentPlot[i-2]->SetLineColor(615);
      
    }
    PercentPlot[i-2]->SetMarkerStyle(20+i);
    PercentPlot[i-2]->GetXaxis()->SetLimits(-6.,31.);
    PercentPlot[i-2]->GetYaxis()->SetRangeUser(0.985,1.);
    
    if(i==2)    PercentPlot[i-2]->Draw("AP");
    else PercentPlot[i-2]->Draw("Psame");
  }


  //  if(logar) cPercent->SetLogy();   
  TLegend* legPercent = new TLegend(0.2,0.3,0.7,0.4);
  legPercent->SetLineColor(0);
  legPercent->SetTextSize(0.035);
  legPercent->SetNColumns(2);
  for(int i =2; i < 5; i++){
    //if(i==0) legPercent->AddEntry(PercentPlot[i-2],"all clusters RMS95","lp");
    //else
    legPercent->AddEntry(PercentPlot[i-2],interval[i],"lp");
  }
  legPercent->Draw();

  
  outname = outputDir+"Percent_angleScans25_"+name;
  cPercent->SaveAs(outname+".eps");
  cPercent->SaveAs(outname+".png");
  cPercent->SaveAs(outname+".pdf");
  cPercent->SaveAs(outname+".root");
  cPercent->SaveAs(outname+".C");


  

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


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
#define measurements 32


void resolutionLiteratureSummary()
{

  Double_t Resolution[measurements];
  Double_t ResolutionError[measurements];
  Double_t Pitch[measurements];
  TString Legend[measurements];
  TString Label[measurements] ={
    "n-in-n, 285 #mum, psi46dig, fit Gen. Err.",
    "n-in-n, 285 #mum, psi46dig, fit Gen. Err.",
    "n-in-n, 285 #mum, psi46dig, RMS #pm pitch/2",
    "n-in-p, 200 #mum, timepix3, fit Gauss",
    "n-in-p, 150 #mum, RD53A (CMS), fit Gen. Err.",
    "n-in-n, 280 #mum, FE-A/FE-B, fit Gauss",
    "3D, n-in-p, 130 #mum, RD53A (CMS), fit Student-t",
    "n-in-n, 280 #mum, FE-A/FE-B, fit Gauss",
    "n-in-n, 280 #mum, FE-A/FE-B, fit Gauss",
    "n-in-n, 280 #mum, FE-A/FE-B, fit Gauss",
    "n-in-n, 280 #mum, FE-A/FE-B, fit Gauss",
    "n-in-p, 100 #mum, RD53A (ATLAS), RMS",
    "n-in-p, 150 #mum, R4S, Recursive RMS",
    "n-in-n, 285 #mum, psi46dig, fit Gauss",
    "ATLASpix Simple, 50 #mum, (?)",
    "DEPFET, 450 #mum, RMS",
    "SOI, 500 #mum, fit Gauss",
    "High-Resistivity CMOS, 25 #mum, (?)",
    "High-Voltage CMOS C3PD, 50 #mum, (?)",
    "n-in-n, 285 #mum, psi46dig, fit Gauss",
    "n-in-p, 100 #mum, RD53A (ATLAS), RMS",
    "3D, n-in-p, 130 #mum, RD53A (CMS), fit Student-t",
    "n-in-p, 150 #mum, RD53A (CMS), fit Gen. Err.",
    "n-in-p, 150 #mum, R4S, Recursive RMS",
    "CLICpix2 prototype, 130 #mum, (?)",
    "DEPFET, 450 #mum, RMS",
    "n in p, High Voltage CMOS, 15 #mum, (?)",
    "SOI 0.2 #mum, 200 #mum, (?)",
    "DEPFET, 450 #mum, RMS",
    "Mimosa26 sensors, 50 #mum, fit Gauss",
    "Mimosa18 sensors, 14 #mum, fit Gauss",
    "FPIX SOI 0.2 #mum, 400 #mum, (?)"
  };


    
  // change color and markers: full = normal sensors, empty = cmos
  // colors dark = 120 GeV, medium = ~ GeV, light = unknown
  int colors[measurements];
  int markers[measurements];
  
  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  TString filename = inputDir+"/TextFiles/LiteratureResolution.txt";
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
	  stream  >> Legend[i] >> Pitch[i] >> Resolution[i] >> ResolutionError[i] >> colors[i] >> markers[i]; 
	  if(print)	 cout << "line " << i << endl;
	  if(print)	 cout  << Legend[i] << " " << Pitch[i] << " " << Resolution[i] << " " << ResolutionError[i] << endl;

	}
    }




  

  TCanvas *cFDB2 = new TCanvas("cFDB2", "FDB resolution", 600, 600);
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
  
  TLegend* legFDB2 = new TLegend(0.13,0.37,0.5,0.88);
  legFDB2->SetLineColor(0);
  legFDB2->SetTextSize(10);
  //legFDB2->SetNColumns(2);

  TLegend* legFDB4 = new TLegend(0.5,0.6,0.8,0.88);
  legFDB4->SetLineColor(0);
  legFDB4->SetTextSize(10);
  //legFDB2->SetNColumns(2);

  TLegend* legFDB3 = new TLegend(0.55,0.12,0.88,0.35);
  legFDB3->SetLineColor(0);
  legFDB3->SetTextSize(10);
  //  legFDB3->SetNColumns(2);

  TGraphErrors* resolutionPlot[measurements];
  for (int i = 0 ; i< measurements ; i++){
    Double_t pitch[]={Pitch[i]};
    Double_t pitchE[]={0.};
    Double_t res[]={Resolution[i]};
    Double_t resE[]={ResolutionError[i]};
    
    resolutionPlot[i] = new TGraphErrors(1,pitch,res,pitchE,resE);

    resolutionPlot[i]->SetMarkerSize(1.5);
    resolutionPlot[i]->SetMarkerColor(colors[i]);
    resolutionPlot[i]->SetMarkerStyle(markers[i]);


    if(i==0){
      resolutionPlot[i]->SetTitle(" ");
      resolutionPlot[i]->GetYaxis()->SetTitle("Resolution [#mum]");
      resolutionPlot[i]->GetXaxis()->SetTitle("Pitch [#mum]");
      resolutionPlot[0]->GetXaxis()->SetTitleFont(43);
      resolutionPlot[0]->GetXaxis()->SetTitleOffset(1.5);
      resolutionPlot[0]->GetXaxis()->SetTitleSize(15); // labels will be 14 pixels

      resolutionPlot[0]->GetYaxis()->SetTitleFont(43);
      resolutionPlot[0]->GetYaxis()->SetTitleSize(15); // labels will be 14 pixels
      resolutionPlot[0]->GetYaxis()->SetTitleOffset(1.8);


      resolutionPlot[0]->GetXaxis()->SetLabelFont(43);
      resolutionPlot[0]->GetXaxis()->SetLabelSize(15); // labels will be 14 pixels
      resolutionPlot[0]->GetYaxis()->SetLabelFont(43);
      resolutionPlot[0]->GetYaxis()->SetLabelSize(15); // labels will be 14 pixels
      
      
      resolutionPlot[i]->GetXaxis()->SetLimits(0.,160.);
      resolutionPlot[i]->GetYaxis()->SetRangeUser(0.,55.);
      resolutionPlot[i]->Draw("AEP");
    }
    else {
      /*
      resolutionPlot[i]->SetMarkerSize(2.5);
      resolutionPlot[i]->SetMarkerColor(kBlue+i);
      resolutionPlot[i]->SetMarkerStyle(20+i);
      if(i==5) resolutionPlot[i]->SetMarkerColor(kGreen+1);
      if(i==21) resolutionPlot[i]->SetMarkerColor(kGreen+2);
      */
      resolutionPlot[i]->Draw("EPsame");      
    }
    if (i <= 17) 
      legFDB2->AddEntry(resolutionPlot[i],Label[i],"p");
    else if(i > 17 && i < 25)
      legFDB4->AddEntry(resolutionPlot[i],Label[i],"p");
    else 
      legFDB3->AddEntry(resolutionPlot[i],Label[i],"p");

      //legFDB2->AddEntry(resolutionPlot[i],Legend[i],"p"); 

  }

  //cFDB2->SetLogy();
  legFDB2->Draw();
  legFDB4->Draw();
  legFDB3->Draw();
  TF1* pi12 = new TF1("pi12","x/TMath::Sqrt(12)",0,160);
  pi12->Draw("lsame");
  
  TString  outname = outputDir+"ResolutionLiterature"; //+name;
  cFDB2->SaveAs(outname+".eps");
  cFDB2->SaveAs(outname+".png");
  cFDB2->SaveAs(outname+".pdf");
  cFDB2->SaveAs(outname+".root");


  /////// now log version! 
  TCanvas *cLog = new TCanvas("cLog", "FDB resolution", 600, 600);
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
  
  for (int i = 0 ; i< measurements ; i++){
    if(i==0){
      resolutionPlot[i]->Draw("AEP");
      resolutionPlot[i]->GetXaxis()->SetLimits(1.,200.);
      resolutionPlot[i]->GetYaxis()->SetRangeUser(0.1,100.);

    }
    else {
      resolutionPlot[i]->Draw("EPsame");      
    }
  }

  pi12->Draw("lsame");
  cLog->SetLogy();
  cLog->SetLogx();
  
  outname = outputDir+"ResolutionLiteratureLog"; //+name;
  cLog->SaveAs(outname+".eps");
  cLog->SaveAs(outname+".png");
  cLog->SaveAs(outname+".pdf");
  cLog->SaveAs(outname+".root");



  
  
}//resolution 



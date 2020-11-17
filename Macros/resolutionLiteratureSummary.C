
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
#define measurements 36


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
    "3D, 285 #mum, timepix3, fit Gauss",
    "n-in-p, 200 #mum, timepix3, fit Gauss",
    "p-in-n, 300 #mum, timepix3, fit Gauss",
    "p-in-n, 300 #mum, timepix3, fit Gauss",
    "n-in-p, 150 #mum, RD53A (CMS), fit Gen. Err.",
    "n-in-n, 280 #mum, FE-A/FE-B, fit Gauss",
    "3D, n-in-p, 130 #mum, RD53A (CMS), Student's t", //Student-t",
    "n-in-n, 280 #mum, FE-A/FE-B, fit Gauss",
    "n-in-n, 280 #mum, FE-A/FE-B, fit Gauss",
    "n-in-n, 280 #mum, FE-A/FE-B, fit Gauss",
    "n-in-n, 280 #mum, FE-A/FE-B, fit Gauss",
    "n-in-p, 100 #mum, RD53A (ATLAS), RMS",
    "n-in-p, 150 #mum, R4S, Recursive RMS",
    "n-in-n, 285 #mum, psi46dig, fit Gauss",
    "ATLASpix Simple, 50 #mum", //, (?)",
    "DEPFET, 450 #mum, RMS",
    "SOI, 500 #mum, fit Gauss",
    "High-Resistivity CMOS, 25 #mum", //, (?)",
    "High-Voltage CMOS C3PD, 50 #mum", //, (?)",
    "n-in-n, 285 #mum, psi46dig, fit Gauss",
    "n-in-p, 100 #mum, RD53A (ATLAS), RMS",
    "3D, n-in-p, 130 #mum, RD53A (CMS), Student's t",// Student-t",
    "n-in-p, 150 #mum, RD53A (CMS), fit Gen. Err.",
    "n-in-p, 150 #mum, R4S, Recursive RMS",
    "CLICpix2 prototype, 130 #mum", //, (?)",
    "DEPFET, 450 #mum, RMS",
    "n in p, High Voltage CMOS, 15 #mum", //, (?)",
    "SOI 0.2 #mum, 200 #mum", //, (?)",
    "DEPFET, 450 #mum, RMS",
    "Mimosa26 sensors, 50 #mum, fit Gauss",
    "SOI 0.2 #mum, 114 #mum, Gauss",
    "Mimosa18 sensors, 14 #mum, fit Gauss",
    "FPIX SOI 0.2 #mum, 400 #mum" //, (?)"
  };
  bool draw[measurements] ={
    true,
    true,
    true,
    true,
    true,
    true,
    true,
    true,
    true,
    true,
    true,
    true,
    true,
    true,
    false,
    true,
    true,
    false,
    false,
    true,
    true,
    true,
    true,
    true,
    false,
    true,
    false,
    false,
    true,
    true,
    true,
    false};


    
  // change color and markers: full = normal sensors, empty = cmos
  // colors dark = 120 GeV, medium = ~ GeV, light = unknown
  int colors[measurements];
  int markers[measurements];
  
  TString inputDir="/home/zoiirene/Programming/r4s_clientsw/Macros/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  TString filename = inputDir+"LiteratureResolution.txt";
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
  legFDB2->SetBorderSize(0);
  legFDB2->SetFillStyle(0);
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

  TF1* pi12 = new TF1("pi12","x/TMath::Sqrt(12)",0,160);
  pi12->Draw("lsame");

  legFDB3->AddEntry(pi12,"pitch/#sqrt{12}","l");
  //cFDB2->SetLogy();
  legFDB2->Draw();
  legFDB4->Draw();
  legFDB3->Draw();
  
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

  TLegend* legFDB2log = new TLegend(0.11,0.37,0.47,0.88);
  legFDB2log->SetLineColor(0);
  legFDB2log->SetBorderSize(0);
  legFDB2log->SetFillStyle(0);
  legFDB2log->SetTextSize(10);
  //legFDB2->SetNColumns(2);

  TLegend* legFDB4log = new TLegend(0.48,0.6,0.78,0.88);
  legFDB4log->SetLineColor(0);
  legFDB4log->SetTextSize(10);
  //legFDB2->SetNColumns(2);

  TLegend* legFDB3log = new TLegend(0.11,0.12,0.88,0.23);
  legFDB3log->SetLineColor(0);
  legFDB3log->SetTextSize(10);
  legFDB3log->SetNColumns(2);
  

  
  for (int i = 0 ; i< measurements ; i++){
    if(i==0){
      resolutionPlot[i]->Draw("AEP");
      resolutionPlot[i]->GetXaxis()->SetLimits(1.,200.);
      resolutionPlot[i]->GetYaxis()->SetRangeUser(0.1,1000.);

    }
    else {
      resolutionPlot[i]->Draw("EPsame");      
    }
    if (i <= 19)
      legFDB2log->AddEntry(resolutionPlot[i],Label[i],"p");
    else if(i > 19 && i < 28)
      legFDB4log->AddEntry(resolutionPlot[i],Label[i],"p");
    else
      legFDB3log->AddEntry(resolutionPlot[i],Label[i],"p");
    
  }

  pi12->Draw("lsame");
  cLog->SetLogy();
  cLog->SetLogx();
  legFDB3log->AddEntry(pi12,"pitch/#sqrt{12}","l");
  
  legFDB2log->Draw();
  legFDB4log->Draw();
  legFDB3log->Draw();
    
  outname = outputDir+"ResolutionLiteratureLog"; //+name;
  cLog->SaveAs(outname+".eps");
  cLog->SaveAs(outname+".png");
  cLog->SaveAs(outname+".pdf");
  cLog->SaveAs(outname+".root");



  


  //remove measurement with not enough info
  TCanvas *cFDB2sel = new TCanvas("cFDB2sel", "FDB resolution", 600, 600);
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
  
  TLegend* legFDB2sel = new TLegend(0.13,0.37,0.5,0.88);
  legFDB2sel->SetLineColor(0);
  legFDB2sel->SetBorderSize(0);
  legFDB2sel->SetFillStyle(0);
  legFDB2sel->SetTextSize(10);
  //legFDB2sel->SetNColumns(2);

  TLegend* legFDB4sel = new TLegend(0.5,0.7,0.8,0.88);
  legFDB4sel->SetLineColor(0);
  legFDB4sel->SetTextSize(10);
  //legFDB2sel->SetNColumns(2);

  TLegend* legFDB3sel = new TLegend(0.55,0.12,0.88,0.35);
  legFDB3sel->SetLineColor(0);
  legFDB3sel->SetTextSize(10);
  //  legFDB3sel->SetNColumns(2);

  //  TGraphErrors* resolutionPlot[measurements];
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
      if (draw[i] == true) resolutionPlot[i]->Draw("EPsame");      
    }
    if (i <= 19 && draw[i] == true) 
      legFDB2sel->AddEntry(resolutionPlot[i],Label[i],"p");
    else if(i > 19 && i < 27 && draw[i] == true)
      legFDB4sel->AddEntry(resolutionPlot[i],Label[i],"p");
    else 
      if(draw[i] == true)  legFDB3sel->AddEntry(resolutionPlot[i],Label[i],"p");

      //legFDB2sel->AddEntry(resolutionPlot[i],Legend[i],"p"); 

  }

  //  TF1* pi12 = new TF1("pi12","x/TMath::Sqrt(12)",0,160);
  pi12->Draw("lsame");

  legFDB3sel->AddEntry(pi12,"pitch/#sqrt{12}","l");
  //cFDB2sel->SetLogy();
  legFDB2sel->Draw();
  legFDB4sel->Draw();
  legFDB3sel->Draw();
  
  outname = outputDir+"ResolutionLiteratureSelected"; //+name;
  cFDB2sel->SaveAs(outname+".eps");
  cFDB2sel->SaveAs(outname+".png");
  cFDB2sel->SaveAs(outname+".pdf");
  cFDB2sel->SaveAs(outname+".root");


  //divided by pitch / sqrt 12
  TCanvas *cdiv = new TCanvas("cdiv", "FDB resolution", 600, 600);
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
  
  TLegend* legFDB2div = new TLegend(0.12,0.37,0.9,0.88);
  legFDB2div->SetLineColor(0);
  legFDB2div->SetBorderSize(0);
  legFDB2div->SetFillStyle(0);
  legFDB2div->SetTextSize(10);
  legFDB2div->SetNColumns(2);
  /*
  TLegend* legFDB4 = new TLegend(0.5,0.6,0.8,0.88);
  legFDB4->SetLineColor(0);
  legFDB4->SetTextSize(10);
  //legFDB2div->SetNColumns(2);

  TLegend* legFDB3 = new TLegend(0.55,0.12,0.88,0.35);
  legFDB3->SetLineColor(0);
  legFDB3->SetTextSize(10);
  //  legFDB3->SetNColumns(2);
  */


  TGraphErrors* resolutionPlotDiv[measurements];
  for (int i = 0 ; i< measurements ; i++){
    Double_t pitch[]={Pitch[i]};
    Double_t pitchE[]={0.};
    Double_t res[]={Resolution[i]/(Pitch[i]/TMath::Sqrt(12.))};
    Double_t resE[]={ResolutionError[i]/(Pitch[i]/TMath::Sqrt(12.))};
    
    resolutionPlotDiv[i] = new TGraphErrors(1,pitch,res,pitchE,resE);

    resolutionPlotDiv[i]->SetMarkerSize(1.5);
    resolutionPlotDiv[i]->SetMarkerColor(colors[i]);
    resolutionPlotDiv[i]->SetMarkerStyle(markers[i]);


    if(i==0){
      resolutionPlotDiv[i]->SetTitle(" ");
      resolutionPlotDiv[i]->GetYaxis()->SetTitle("Resolution/(Pitch/#sqrt{12})");
      resolutionPlotDiv[i]->GetXaxis()->SetTitle("Pitch [#mum]");
      resolutionPlotDiv[0]->GetXaxis()->SetTitleFont(43);
      resolutionPlotDiv[0]->GetXaxis()->SetTitleOffset(1.5);
      resolutionPlotDiv[0]->GetXaxis()->SetTitleSize(15); // labels will be 14 pixels

      resolutionPlotDiv[0]->GetYaxis()->SetTitleFont(43);
      resolutionPlotDiv[0]->GetYaxis()->SetTitleSize(15); // labels will be 14 pixels
      resolutionPlotDiv[0]->GetYaxis()->SetTitleOffset(1.8);


      resolutionPlotDiv[0]->GetXaxis()->SetLabelFont(43);
      resolutionPlotDiv[0]->GetXaxis()->SetLabelSize(15); // labels will be 14 pixels
      resolutionPlotDiv[0]->GetYaxis()->SetLabelFont(43);
      resolutionPlotDiv[0]->GetYaxis()->SetLabelSize(15); // labels will be 14 pixels
      
      
      resolutionPlotDiv[i]->GetXaxis()->SetLimits(0.,160.);
      resolutionPlotDiv[i]->GetYaxis()->SetRangeUser(0.,1.5);
      resolutionPlotDiv[i]->Draw("AEP");
    }
    else {
      /*
      resolutionPlotDiv[i]->SetMarkerSize(2.5);
      resolutionPlotDiv[i]->SetMarkerColor(kBlue+i);
      resolutionPlotDiv[i]->SetMarkerStyle(20+i);
      if(i==5) resolutionPlotDiv[i]->SetMarkerColor(kGreen+1);
      if(i==21) resolutionPlotDiv[i]->SetMarkerColor(kGreen+2);
      */
      resolutionPlotDiv[i]->Draw("EPsame");      
    }
    
    //if (i <= 17) 
    legFDB2div->AddEntry(resolutionPlotDiv[i],Label[i],"p");
      /*
      else if(i > 17 && i < 25)
      legFDB4->AddEntry(resolutionPlotDiv[i],Label[i],"p");
    else 
      legFDB3->AddEntry(resolutionPlotDiv[i],Label[i],"p");
    */
      //legFDB2div->AddEntry(resolutionPlotDiv[i],Legend[i],"p"); 

  }

  //  TF1* pi12 = new TF1("pi12","x/TMath::Sqrt(12)",0,160);
  //  pi12->Draw("lsame");
  TF1* div12 = new TF1("div12","1.",0,160);
  div12->SetLineColor(kRed+2);
  div12->SetLineStyle(2); 
  div12->Draw("lsame");


  
  // legFDB2div->AddEntry(pi12,"pitch/#sqrt{12}","l");
  //cdiv->SetLogy();
  //legFDB2div->Draw();
  //legFDB4->Draw();
  //legFDB3->Draw();
  
  outname = outputDir+"ResolutionLiterature_div_pdivsqrt12_noleg"; //+name;
  cdiv->SaveAs(outname+".eps");
  cdiv->SaveAs(outname+".png");
  cdiv->SaveAs(outname+".pdf");
  cdiv->SaveAs(outname+".root");



  


  
}//resolution 



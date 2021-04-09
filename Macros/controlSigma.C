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
#define measurements 9
#define irradiation 3

void controlSigma(){
  TString outputDir="/home/zoiirene/Output/Plots/";
Double_t Resolution[irradiation][measurements];
Double_t Tracks[irradiation][measurements];
 Double_t sigma[measurements] = {2,3,4,5,6,7,9,11,12};;
 Int_t sigma_i[measurements] = {2,3,4,5,6,7,9,11,12};;
 int irr [irradiation];
 TString ss_irr[irradiation];
 TString filename = "res_tracks_RMS";
 for(int k = 0;k <measurements; k++){
   cout << 0 << endl;
   std::string ns = std::to_string(sigma_i[k]);

   
   TString allname = filename+ns+".txt";

   cout << allname << endl;
   ifstream stream(allname);
   std::string line;
   if(!stream.is_open())
     {
       cout << " File " << allname << " not opened" << endl;
     }
   else
     {
       std::getline(stream,line);
       if(print) cout << " first line " << line << endl;
       //      while(!stream.eof())
       for(int i= 0; i < irradiation; i++)
	 {
	   stream  >> irr[i] >> Resolution[i][k] >> Tracks[i][k];
	   if(print)	 cout << "line " << i << endl;
	   if(print)	 cout  <<  " " << irr[i] << " " << Resolution[i][k] << " " << Tracks[i][k] << endl;
	   //  i++;
	   if(irr[i] == 0) ss_irr[i] = "Non-irradiated, 120 V";
	   else if (irr[i] == 1) ss_irr[i] = "2.1#times10^{15}, proton, 800 V";
	   else if (irr[i] == 2) ss_irr[i] = "3.6#times10^{15}, neutron, 800 V";
	   //	  ss_irr[i].Form("%f",irr[i]);
	 }
     }
   
 } // k measurements




TCanvas *c3 = new TCanvas("c3", "FDB resolution", 600, 600);
  gROOT->SetStyle("Plain");
  c3->SetLeftMargin(0.12);
  c3->SetRightMargin(0.02);
  //   c3->SetBottomMargin(.9);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //  gEnv->Print();

  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(24);

  
  
  TPad *pad2 = new TPad("pad2","",0,0.,1,0.4);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.25);
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.02);
  
  pad2->Draw();

  //pad2->SetFillStyle(4000); //will be transparent
  //#pad2->SetFrameFillStyle(0);
  TPad *pad1 = new TPad("pad1","",0,0.4,1,1);  
  pad1->Draw();
  pad1->SetBottomMargin(0.03);
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.02);
  
  pad1->cd();
  gPad->SetTicks(1,1);

  TGraph* resolutionPlot[irradiation];
  
  resolutionPlot[0] = new TGraph(measurements,sigma,Resolution[0]);
  resolutionPlot[1] = new TGraph(measurements,sigma,Resolution[1]);
  resolutionPlot[2] = new TGraph(measurements,sigma,Resolution[2]);


  resolutionPlot[0]->SetTitle(" ");
  resolutionPlot[0]->GetYaxis()->SetTitle("Reduced RMS #delta_{#Delta_{x}} [#mum]"); //Single hit resolution [#mum]");
  //resolutionPlot[0]->GetXaxis()->SetTitle("#theta [deg]");

  resolutionPlot[0]->GetYaxis()->SetTitleFont(43);
  resolutionPlot[0]->GetYaxis()->SetTitleSize(30); // labels will be 14 pixels

  resolutionPlot[0]->GetYaxis()->SetLabelFont(43);
  resolutionPlot[0]->GetYaxis()->SetLabelSize(30); // labels will be 14 pixels
  
  resolutionPlot[0]->GetXaxis()->SetLabelFont(43);
  resolutionPlot[0]->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels

  resolutionPlot[0]->SetMarkerColor(kBlack);
  resolutionPlot[0]->SetLineColor(kBlack);
  resolutionPlot[0]->SetMarkerStyle(20);

  resolutionPlot[0]->GetYaxis()->SetRangeUser(0.,10.);



  resolutionPlot[1]->SetMarkerColor(kGreen+1);
  resolutionPlot[1]->SetLineColor(kGreen+1);
  resolutionPlot[1]->SetMarkerStyle(21);

				      

  resolutionPlot[2]->SetMarkerColor(kMagenta);
  resolutionPlot[2]->SetLineColor(kMagenta);
  resolutionPlot[2]->SetMarkerStyle(26);



  resolutionPlot[0]->GetXaxis()->SetLimits(0.,13.);
  resolutionPlot[0]->SetMarkerSize(1.);
  resolutionPlot[0]->Draw("AEP");
  resolutionPlot[1]->SetMarkerSize(1.);
  resolutionPlot[1]->Draw("EPsame");
  resolutionPlot[2]->SetMarkerSize(1.);
  resolutionPlot[2]->Draw("EPsame");


  TLine *  linemax2 = new TLine( 6.,0.,6.,10.);
  linemax2->SetLineColor(kGray);
  linemax2->SetLineWidth(2);
  linemax2->SetLineStyle(2);
  linemax2->Draw("same");


  TLegend* leg3 = new TLegend(0.15,0.67,0.55,0.87);
  leg3->SetLineColor(0);
  leg3->AddEntry(resolutionPlot[0],ss_irr[0],"lp");
  leg3->AddEntry(resolutionPlot[1],ss_irr[1],"lp");
  leg3->AddEntry(resolutionPlot[2],ss_irr[2],"lp");


  leg3->Draw();

 
  pad2->cd();
  //  pad2->SetLogy();
  //  gROOT->SetStyle("Plain");
  TGraph* clsizePlot[irradiation];
  gPad->SetTicks(1,1);

  clsizePlot[0] = new TGraph(measurements,sigma,Tracks[0]);
  clsizePlot[1] = new TGraph(measurements,sigma,Tracks[1]);
  clsizePlot[2] = new TGraph(measurements,sigma,Tracks[2]);

  
  clsizePlot[0]->SetTitle(" ");
  clsizePlot[0]->GetYaxis()->SetTitle("Tracks [%]");
  clsizePlot[0]->GetYaxis()->SetNdivisions(5);
  clsizePlot[0]->GetXaxis()->SetTitle("N RMS");


  clsizePlot[0]->GetXaxis()->SetTitleFont(43);
  clsizePlot[0]->GetXaxis()->SetTitleSize(30); // labels will be 14 pixels
  clsizePlot[0]->GetXaxis()->SetTitleOffset(1.8); // labels will be 14 pixels
  clsizePlot[0]->GetXaxis()->SetLabelFont(43);
  clsizePlot[0]->GetXaxis()->SetLabelSize(30); // labels will be 14 pixels
  
  clsizePlot[0]->GetYaxis()->SetTitleFont(43);
  clsizePlot[0]->GetYaxis()->SetTitleSize(30); // labels will be 14 pixels
  clsizePlot[0]->GetYaxis()->SetLabelFont(43);
  clsizePlot[0]->GetYaxis()->SetLabelSize(30); // labels will be 14 pixels
  clsizePlot[0]->GetYaxis()->SetTitleOffset(1.1); // labels will be 14 pixels
  
  clsizePlot[0]->SetMarkerColor(kBlack);
  clsizePlot[0]->SetLineColor(kBlack);
  clsizePlot[0]->SetMarkerStyle(20);
  clsizePlot[0]->SetMarkerSize(1.);
  clsizePlot[0]->GetXaxis()->SetLimits(0.,13.);
  clsizePlot[0]->GetYaxis()->SetRangeUser(92.,100.);
  clsizePlot[0]->Draw("AEP");
  //  clsizePlot[0]->Draw("AEPY+");


  clsizePlot[1]->SetMarkerColor(kGreen+1);
  clsizePlot[1]->SetLineColor(kGreen+1);
  clsizePlot[1]->SetMarkerStyle(21);
  clsizePlot[1]->SetMarkerSize(1.);
  clsizePlot[1]->Draw("EPsame");


  clsizePlot[2]->SetMarkerColor(kMagenta);
  clsizePlot[2]->SetMarkerStyle(26);
  clsizePlot[2]->SetMarkerSize(1.);  
  clsizePlot[2]->SetLineColor(kMagenta);
  clsizePlot[2]->Draw("EPsame");

  
  TLine *  linemax = new TLine( 6.,92.,6.,100.);
  linemax->SetLineColor(kGray);
  linemax->SetLineWidth(2);
  linemax->SetLineStyle(2);
  linemax->Draw("same");




  pad1->cd(); 
  //  TDR2(c3,0,0);
  TString outname = outputDir+"Resolution_Tracks_RMS";
  c3->SaveAs(outname+".eps");
  c3->SaveAs(outname+".png");
  c3->SaveAs(outname+".pdf");
  c3->SaveAs(outname+".root");
  c3->SaveAs(outname+".C");
}

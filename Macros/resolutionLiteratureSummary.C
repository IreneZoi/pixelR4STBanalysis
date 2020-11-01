
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
    "n-in-n, 285 #mum, psi46dig",
    "n-in-n, 285 #mum, psi46dig",
    "n-in-n, 285 #mum, psi46dig",
    "n-in-p, 200 #mum, timepix3",
    "n-in-p, 150 #mum, RD53A (CMS)",
    "n-in-n, 280 #mum, FE-A/FE-B",
    "3D, n-in-p, 130 #mum, RD53A (CMS)",
    "n-in-n, 280 #mum, FE-A/FE-B",
    "n-in-n, 280 #mum, FE-A/FE-B",
    "n-in-n, 280 #mum, FE-A/FE-B",
    "n-in-n, 280 #mum, FE-A/FE-B",
    "n-in-p, 100 #mum, RD53A (ATLAS)",
    "n-in-p, 150 #mum, R4S",
    "n-in-n, 285 #mum, psi46dig",
    "ATLASpix Simple, 50 #mum",
    "DEPFET, 450 #mum",
    "SOI, 500 #mum",
    "High-Resistivity CMOS, 25 #mum",
    "High-Voltage CMOS C3PD, 50 #mum",
    "n-in-n, 285 #mum, psi46dig",
    "n-in-p, 100 #mum, RD53A (ATLAS)",
    "3D, n-in-p, 130 #mum, RD53A (CMS)",
    "n-in-p, 150 #mum, RD53A (CMS)",
    "n-in-p, 150 #mum, R4S",
    "CLICpix2 prototype, 130 #mum",
    "DEPFET, 450 #mum",
    "n in p, High Voltage CMOS, 15 #mum",
    "SOI 0.2 #mum, 200 #mum",
    "DEPFET, 450 #mum",
    "Mimosa26 sensors, 50 #mum",
    "Mimosa18 sensors, 14 #mum",
    "FPIX SOI 0.2 #mum, 400 #mum"
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
      resolutionPlot[i]->GetYaxis()->SetRangeUser(0.1,55.);
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

  /*
  TCanvas *c = new TCanvas("c", "resolution vs inverse beam energy", 1500, 900);
  gPad->SetTicks(1,1);

  TGraph* resolutionPlotInv = new TGraphErrors(measurements,BeamEnergyInverse,Resolution,BeamEnergyError,ResolutionError);

  resolutionPlotInv->SetTitle(" ");
  resolutionPlotInv->GetYaxis()->SetTitle("Resolution [#mum]");
  resolutionPlotInv->GetXaxis()->SetTitle("1/Beam Energy [GeV^{-1}]");
  resolutionPlotInv->SetMarkerSize(2.5);
  resolutionPlotInv->SetLineColor(2);
  resolutionPlotInv->SetMarkerColor(2);
  resolutionPlotInv->SetMarkerStyle(20);
  resolutionPlotInv->SetLineWidth(2);
  resolutionPlotInv->SetLineStyle(2);
  resolutionPlotInv->GetYaxis()->SetRangeUser(0.,20.);
  resolutionPlotInv->GetXaxis()->SetLimits(0.,1.);

  //  TF1 *fit1 = new TF1("fit1","pol1", 0.15, 1.255); 
  TF1 *fit1 = new TF1("fit1","sqrt([0]*[0]+([1]*x)*([1]*x))", 0., 0.8); 
  fit1->SetLineColor(kBlue);
  fit1->SetParameter(0,3);
  fit1->SetParameter(1,15);
  fit1->SetParName(0,"#sigma_{hit}");
  fit1->SetParName(1,"#sigma_{MS}");
  
  resolutionPlotInv->Fit("fit1","R");
  resolutionPlotInv->Draw("AEP");

  // TPaveStats * ps = (TPaveStats *)fit1->FindObject("stats");
  // ps->SetX1NDC(0.2);
  // ps->SetX2NDC(0.2);
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.85);
  
  TLegend* leg = new TLegend(0.65,0.4,0.75,0.6);
  leg->SetLineColor(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(resolutionPlotInv,"DUT "+detectorB , "ep");
  //  leg->AddEntry(fit1,"p_{0}+x*p_{1}" , "l");
  leg->AddEntry(fit1,"#sqrt{#sigma_{hit}^{2}+(#sigma_{MS}/p)^{2}}" , "l");
  leg->Draw();

  name = outputDir+"Res_vs_InvEnergy_newFit_updatedcode_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC;
  c->SaveAs(name+".eps");
  c->SaveAs(name+".pdf");
  c->SaveAs(name+".png");


  TCanvas *c4 = new TCanvas("c4", "resolution vs inverse beam energy squared", 1500, 900);
  gPad->SetTicks(1,1);

  TGraph* resolutionPlotInvGeorg = new TGraphErrors(measurements,BeamEnergyInverseSquare,ResolutionSquare,BeamEnergyError,ResolutionErrorSquare);

  resolutionPlotInvGeorg->SetTitle(" ");
  resolutionPlotInvGeorg->GetYaxis()->SetTitle("Resolution^{2} [#mum^{2}]");
  resolutionPlotInvGeorg->GetXaxis()->SetTitle("1/(Beam Energy)^{2} [GeV^{-2}]");
  resolutionPlotInvGeorg->SetMarkerSize(2.5);
  resolutionPlotInvGeorg->SetLineColor(2);
  resolutionPlotInvGeorg->SetMarkerColor(2);
  resolutionPlotInvGeorg->SetMarkerStyle(20);
  resolutionPlotInvGeorg->SetLineWidth(2);
  resolutionPlotInvGeorg->SetLineStyle(2);
  resolutionPlotInvGeorg->GetYaxis()->SetRangeUser(0.,300.);
  resolutionPlotInvGeorg->GetXaxis()->SetLimits(0.,1.);

  //  TF1 *fit1 = new TF1("fit1","pol1", 0.15, 1.255); 
  TF1 *fit3 = new TF1("fit3","pol1", 0., 0.8); 
  fit3->SetLineColor(kBlue);
  //  fit1->SetParameter(0,3);
  //fit1->SetParameter(1,15);
  fit3->SetParName(0,"#sigma_{hit}^{2}");
  fit3->SetParName(1,"#sigma_{MS}^{2}");
  
  resolutionPlotInvGeorg->Fit("fit3","R");
  resolutionPlotInvGeorg->Draw("AEP");

  // TPaveStats * ps = (TPaveStats *)fit1->FindObject("stats");
  // ps->SetX1NDC(0.2);
  // ps->SetX2NDC(0.2);
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.85);

  ostringstream strfit[2];
  TString ss_fit[2];
  ostringstream strfit_err[2];
  TString ss_fit_err[2];
  for(int i=0; i<2;i++)
    {
      strfit[i] << setprecision(3) << sqrt(fit3->GetParameter(i));
      ss_fit[i]=strfit[i].str();
      strfit_err[i] << setprecision(1) << 0.5*fit3->GetParError(i)/sqrt(fit3->GetParameter(i));
      ss_fit_err[i]=strfit_err[i].str();
    }
  
  TLatex Tl_2;
  Tl_2.SetTextAlign(12);
  Tl_2.SetTextSize(0.03);
  Tl_2.DrawLatexNDC(0.55,0.8,"#sigma_{hit} = ("+ss_fit[0]+" #pm "+ss_fit_err[0]+") #mum");
  Tl_2.DrawLatexNDC(0.55,0.72,"#sigma_{MS} = ("+ss_fit[1]+" #pm "+ss_fit_err[1]+") #mum*GeV");
  
  
  TLegend* leg4 = new TLegend(0.45,0.15,0.8,0.35);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.025);
  leg4->AddEntry(resolutionPlot,"c"+detectorB+" "+labelB , "ep");
  leg4->AddEntry(resolutionPlot,"A: c"+detectorA+" "+labelA ,"");
  leg4->AddEntry(resolutionPlot,"C: c"+detectorC+" "+labelC ,"");
  leg4->AddEntry(resolutionPlot,info ,"");
  leg4->AddEntry(fit3,"#sigma_{hit}^{2}+(#sigma_{MS}/p)^{2}" , "l");
  leg4->Draw();

  name = outputDir+"Res_vs_InvEnergySquared_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC;  
  c4->SaveAs(name+".eps");
  c4->SaveAs(name+".pdf");
  c4->SaveAs(name+".png");

  

  // TCanvas *c3 = new TCanvas("c3", "resolution vs inverse beam energy old", 1500, 900);
  // gPad->SetTicks(1,1);

  // TGraph* resolutionPlotInvOLD = new TGraphErrors(measurements,BeamEnergyInverse,Resolution,BeamEnergyError,ResolutionError);

  // resolutionPlotInvOLD->SetTitle(" ");
  // resolutionPlotInvOLD->GetYaxis()->SetTitle("Resolution [#mum]");
  // resolutionPlotInvOLD->GetXaxis()->SetTitle("1/Beam Energy [GeV^{-1}]");
  // resolutionPlotInvOLD->SetMarkerSize(2.5);
  // resolutionPlotInvOLD->SetLineColor(2);
  // resolutionPlotInvOLD->SetMarkerColor(2);
  // resolutionPlotInvOLD->SetMarkerStyle(20);
  // resolutionPlotInvOLD->SetLineWidth(2);
  // resolutionPlotInvOLD->SetLineStyle(2);
  // resolutionPlotInvOLD->GetYaxis()->SetRangeUser(0.,20.);
  // resolutionPlotInvOLD->GetXaxis()->SetLimits(0.,1.6);

  // TF1 *fit2 = new TF1("fit2","pol1", 0.15, 1.255); 
  // //  TF1 *fit1 = new TF1("fit1","sqrt([0]*[0]+([1]*x)*([1]*x))", 0., 0.8); 
  // fit2->SetLineColor(kBlue);
  // // fit1->SetParameter(0,3);
  // // fit1->SetParameter(1,15);
  // // fit1->SetParameterName(0,"#sigma_{hit}");
  // // fit1->SetParameterName(1,"#sigma_{MS}");
  
  // resolutionPlotInvOLD->Fit("fit2");
  // resolutionPlotInvOLD->Draw("AEP");

  // // TPaveStats * ps = (TPaveStats *)fit1->FindObject("stats");
  // // ps->SetX1NDC(0.2);
  // // ps->SetX2NDC(0.2);
  // gStyle->SetStatX(0.5);
  // gStyle->SetStatY(0.85);
  
  // TLegend* leg3 = new TLegend(0.65,0.4,0.75,0.6);
  // leg3->SetLineColor(0);
  // leg3->SetTextSize(0.03);
  // leg3->AddEntry(resolutionPlotInvOLD,"DUT "+detectorB , "ep");
  // leg3->AddEntry(fit2,"p_{0}+x*p_{1}" , "l");
  // //  leg->AddEntry(fit1,"#sqrt{#sigma_{hit}^{2}+(#sigma_{MS}/p)^{2}}" , "l");
  // leg3->Draw();

  // name = outputDir+"Res_vs_InvEnergy_oldFit_updatedcode_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC;
  // c3->SaveAs(name+".png"); 
  // c3->SaveAs(name+".pdf"); 
  // c3->SaveAs(name+".eps"); 
  */

  
}//resolution 



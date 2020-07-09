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

#define measurements 12
bool print=true;
using namespace std;


void compareAngles()
{

  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  TString sensor[measurements];
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
  


  TString inputfile;
  TFile * file[measurements];
  TH1F * h_ncol[measurements];
  TH1F * h_madx[measurements];
  TH1F * h_madxvsq[measurements];

  TString Run, Path;
  ostringstream strs[measurements];
  
  if(print) cout << "Getting files and hists"  << endl;
  for(int i=0; i<measurements; i++)
    {
      if(print) cout << "Run: " << runs[i] << endl;
      Run.Form("%d",runs[i]);
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


  if(print) cout << "Plotting nrow for 25 um "  << endl;
  //i da 0 a 6 

  TCanvas *c2 = new TCanvas("c2", "nrow vs xm - 25 um", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend* leg2 = new TLegend(0.15,0.15,0.75,0.35);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.024);


  for(int i = 0; i<7 ; i++)
    {
      h_ncol[i]->SetLineWidth(2);
      int color = i+1;
      if(i+1 == 5) color = 25;
      h_ncol[i]->SetLineColor(color);
      h_ncol[i]->GetYaxis()->SetRangeUser(1.,2.5);
      h_ncol[i]->Draw("histsame");
      leg2->AddEntry(h_ncol[i],pitch[i]+" #mum pitch, "+Short[i]+" "+bias[i]+"V, irr "+irr[i]+"x10^{15} neq" , "l");

    }
  leg2->Draw();
  TString name = outputDir+"Nrow_vs_xmod_25_forResolution";
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".pdf");

  if(print) cout << "Plotting nrow for 50 um "  << endl;

  TCanvas *c3 = new TCanvas("c3", "nrow vs xm - 50 um", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend* leg3 = new TLegend(0.15,0.15,0.75,0.35);
  leg3->SetLineColor(0);
  leg3->SetTextSize(0.024);

  for(int i = 7; i<measurements ; i++)
    {
      h_ncol[i]->SetLineWidth(2);
      int color = i+1;
      if(i+1 == 10) color = 46;
      h_ncol[i]->SetLineColor(color);
      h_ncol[i]->GetYaxis()->SetRangeUser(1.,2.5);
      h_ncol[i]->Draw("histsame");
      leg3->AddEntry(h_ncol[i],pitch[i]+" #mum pitch, "+Short[i]+" "+bias[i]+"V, irr "+irr[i]+"x10^{15} neq" , "l");

    }
  leg3->Draw();
  name = outputDir+"Nrow_vs_xmod_50_forResolution";
  c3->SaveAs(name+".eps");
  c3->SaveAs(name+".png");
  c3->SaveAs(name+".pdf");




    if(print) cout << "Plotting madx for 25 um "  << endl;
  //i da 0 a 6 

  TCanvas *c4 = new TCanvas("c4", "madx vs xm - 25 um", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend* leg4 = new TLegend(0.15,0.2,0.75,0.4);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.03);


  for(int i = 0; i<7 ; i++)
    {
      h_madx[i]->SetLineWidth(2);
      int color = i+1;
      if(i+1 == 5) color = 25;
      h_madx[i]->SetLineColor(color);
      h_madx[i]->GetYaxis()->SetRangeUser(0.,0.006);
      h_madx[i]->Draw("histsame");
      leg4->AddEntry(h_madx[i],pitch[i]+" #mum pitch, "+Short[i]+" "+bias[i]+"V, irr "+irr[i]+"x10^{15} neq" , "l");

    }
  leg4->Draw();
  name = outputDir+"Madx_vs_xmod_25_forResolution";
  c4->SaveAs(name+".eps");
  c4->SaveAs(name+".png");
  c4->SaveAs(name+".pdf");

  if(print) cout << "Plotting madx for 50 um "  << endl;

  TCanvas *c5 = new TCanvas("c5", "madx vs xm - 50 um", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend* leg5 = new TLegend(0.15,0.2,0.75,0.4);
  leg5->SetLineColor(0);
  leg5->SetTextSize(0.03);

  for(int i = 7; i<measurements ; i++)
    {
      h_madx[i]->SetLineWidth(2);
      int color = i+1;
      if(i+1 == 10) color = 46;
      h_madx[i]->SetLineColor(color);
      h_madx[i]->GetYaxis()->SetRangeUser(0.,0.008);
      h_madx[i]->Draw("histsame");
      leg5->AddEntry(h_madx[i],pitch[i]+" #mum pitch, "+Short[i]+" "+bias[i]+"V, irr "+irr[i]+"x10^{15} neq" , "l");

    }
  leg5->Draw();
  name = outputDir+"Madx_vs_xmod_50_forResolution";
  c5->SaveAs(name+".eps");
  c5->SaveAs(name+".png");
  c5->SaveAs(name+".pdf");


  if(print) cout << "Plotting madx vs qfor 25 um "  << endl;
  //i da 0 a 6 

  TCanvas *c6 = new TCanvas("c6", "madx vs q - 25 um", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend* leg6 = new TLegend(0.15,0.7,0.75,0.85);
  leg6->SetLineColor(0);
  leg6->SetTextSize(0.025);


  for(int i = 0; i<7 ; i++)
    {
      h_madxvsq[i]->SetLineWidth(2);
      int color = i+1;
      if(i+1 == 5) color = 25;
      h_madxvsq[i]->SetLineColor(color);
      h_madxvsq[i]->GetYaxis()->SetRangeUser(0.,0.1);
      h_madxvsq[i]->Draw("histsame");
      leg6->AddEntry(h_madxvsq[i],pitch[i]+" #mum pitch, "+Short[i]+" "+bias[i]+"V, irr "+irr[i]+"x10^{15} neq" , "l");

    }
  leg6->Draw();
  name = outputDir+"Madx_vs_q_25_forResolution";
  c6->SaveAs(name+".eps");
  c6->SaveAs(name+".png");
  c6->SaveAs(name+".pdf");

  if(print) cout << "Plotting madxvsq for 50 um "  << endl;

  TCanvas *c7 = new TCanvas("c7", "madxvsq vs xm - 50 um", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend* leg7 = new TLegend(0.15,0.7,0.75,0.85);
  leg7->SetLineColor(0);
  leg7->SetTextSize(0.025);

  for(int i = 7; i<measurements ; i++)
    {
      h_madxvsq[i]->SetLineWidth(2);
      int color = i+1;
      if(i+1 == 10) color = 46;
      h_madxvsq[i]->SetLineColor(color);
      h_madxvsq[i]->GetYaxis()->SetRangeUser(0.,0.1);
      h_madxvsq[i]->Draw("histsame");
      leg7->AddEntry(h_madxvsq[i],pitch[i]+" #mum pitch, "+Short[i]+" "+bias[i]+"V, irr "+irr[i]+"x10^{15} neq" , "l");

    }
  leg7->Draw();
  name = outputDir+"Madx_vs_q_50_forResolution";
  c7->SaveAs(name+".eps");
  c7->SaveAs(name+".png");
  c7->SaveAs(name+".pdf");

  
  /*

  

  TCanvas *c = new TCanvas("c", "resolution vs inverse beam energy", 1500, 900);
  gPad->SetTicks(1,1);

  TGraph* resolutionPlotInv = new TGraphErrors(BeamEnergies,BeamEnergyInverse,Resolution,BeamEnergyError,ResolutionError);

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

  name = outputDir+"Res_vs_InvEnergy_"+func+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC;
  c->SaveAs(name+".eps");
  c->SaveAs(name+".pdf");
  c->SaveAs(name+".png");


  TCanvas *c4 = new TCanvas("c4", "resolution vs inverse beam energy squared", 1500, 900);
  gPad->SetTicks(1,1);

  TGraph* resolutionPlotInvGeorg = new TGraphErrors(BeamEnergies,BeamEnergyInverseSquare,ResolutionSquare,BeamEnergyError,ResolutionErrorSquare);

  resolutionPlotInvGeorg->SetTitle(" ");
  resolutionPlotInvGeorg->GetYaxis()->SetTitle("Resolution^{2} [#mum^{2}]");
  resolutionPlotInvGeorg->GetXaxis()->SetTitle("1/(Beam Energy)^{2} [GeV^{-2}]");
  resolutionPlotInvGeorg->SetMarkerSize(2.5);
  resolutionPlotInvGeorg->SetLineColor(2);
  resolutionPlotInvGeorg->SetMarkerColor(2);
  resolutionPlotInvGeorg->SetMarkerStyle(20);
  resolutionPlotInvGeorg->SetLineWidth(2);
  resolutionPlotInvGeorg->SetLineStyle(2);
  resolutionPlotInvGeorg->GetYaxis()->SetRangeUser(0.,200.);
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
  Tl_2.DrawLatexNDC(0.55,0.6,"#sigma_{hit} = ("+ss_fit[0]+" #pm "+ss_fit_err[0]+") #mum");
  Tl_2.DrawLatexNDC(0.55,0.52,"#sigma_{MS} = ("+ss_fit[1]+" #pm "+ss_fit_err[1]+") #mum*GeV");
  
  
  TLegend* leg4 = new TLegend(0.45,0.2,0.8,0.4);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.025);
  leg4->AddEntry(resolutionPlot,"c"+detectorB+" "+labelB , "ep");
  leg4->AddEntry(resolutionPlot,"A: c"+detectorA+" "+labelA ,"");
  leg4->AddEntry(resolutionPlot,"C: c"+detectorC+" "+labelC ,"");
  leg4->AddEntry(resolutionPlot,info ,"");
  leg4->AddEntry(fit3,"#sigma_{hit}^{2}+(#sigma_{MS}/p)^{2}" , "l");
  leg4->Draw();

  name = outputDir+"Res_vs_InvEnergySquared_"+func+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC;  
  c4->SaveAs(name+".eps");
  c4->SaveAs(name+".pdf");
  c4->SaveAs(name+".png");

  

  // TCanvas *c3 = new TCanvas("c3", "resolution vs inverse beam energy old", 1500, 900);
  // gPad->SetTicks(1,1);

  // TGraph* resolutionPlotInvOLD = new TGraphErrors(BeamEnergies,BeamEnergyInverse,Resolution,BeamEnergyError,ResolutionError);

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

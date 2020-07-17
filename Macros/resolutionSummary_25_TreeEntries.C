
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
#define measurements 7
#define irradiation 3

void resolutionSummary_25_TreeEntries(TString name = "preliminary_RMSself", TString func = "RMSself", bool unfolding = false)
{

  Double_t Resolution[measurements];
  Double_t ResolutionUnfolded[measurements];
  Double_t ResolutionError[measurements];
  Double_t ResolutionUnfoldedError[measurements];
  
  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  TString sensor[measurements];
  Double_t irr[measurements];
  TString ss_irr[measurements];
  TString pitch[measurements];
  TString bias[measurements];
  TString angle[measurements];
  TString beam[measurements];

  TString runs[measurements];
  TString specific[measurements];
  TString Short[measurements];
  TString thr[measurements];

  TString filename = inputDir+"/TextFiles/resolution_25gain1_RMSself.txt";
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
	  stream  >> sensor[i] >> irr[i] >> pitch[i] >> bias[i] >> angle[i] >> beam[i] >> Resolution[i] >> ResolutionError[i] >>runs[i] >> specific[i] >> Short[i] >> thr[i];
	  if(print)	 cout << "line " << i << endl;
	  if(print)	 cout  << sensor[i] << " " << irr[i] << " " << pitch[i] << " " << bias[i] << " " << angle[i] << " " << beam[i] << " " << Resolution[i] << " " << ResolutionError[i] <<" " << runs[i] << " " << specific[i] << " " << Short[i] << " " << thr[i] << endl;
	  //  i++;
	  if(irr[i] == 0) ss_irr[i] = "no irr";
	  else if (irr[i] == 2) ss_irr[i] = "2.1 #times 10^{15}";
	  else if (irr[i] == 4) ss_irr[i] = "3.6 #times 10^{15}";
	  //	  ss_irr[i].Form("%f",irr[i]);
	  cout << " resolution for sensor " << sensor[i] << " is " << Resolution[i] << endl;
	}
    }


  Double_t x1[irradiation] ={irr[0],irr[1],irr[2]};
  Double_t Pstop_default_FTH_25_irr[] = {Resolution[0]};
  Double_t Pstop_default_FTH_25_irr_err[] = {ResolutionError[0]};
  Double_t Pstop_default_FDB_25_preirr[] = {Resolution[1]};
  Double_t Pstop_default_FDB_25_preirr_err[] = {ResolutionError[1]};
  Double_t x2[] = {irr[0]};
  Double_t x3[] = {irr[1]};
  Double_t x4[] = {irr[2]};
  Double_t Pstop_RD53Apads_FDB_25[] = {Resolution[3]};
  Double_t Pstop_RD53Apads_FDB_25_err[] = {ResolutionError[3]};
  Double_t Pspray_default_FDB_25[] = {Resolution[4]};
  Double_t Pspray_default_FDB_25_err[] = {ResolutionError[4]};
  Double_t Pspray_RD53Apads_FDB_25[] = {Resolution[5]};
  Double_t Pspray_RD53Apads_FDB_25_err[] = {ResolutionError[5]};

  Double_t Pstop_default_FDB_25_irr[] = {Resolution[6]};
  Double_t Pstop_default_FDB_25_irr_err[] = {ResolutionError[6]};
  Double_t Pstop_default_FDB_25_gain2_irr[] = {Resolution[2]};
  Double_t Pstop_default_FDB_25_gain2_irr_err[] = {ResolutionError[2]};
  

  TCanvas *cFDB2 = new TCanvas("cFDB2", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  Double_t FDB1[] = {0};
  Double_t FDB2[] = {1};
  Double_t FDB3[] = {2};
  Double_t FDB4[] = {3};
  Double_t FTH[] = {4};
  Double_t FDB5[] = {5};
  Double_t FDB6[] = {6};

  int gain1_25 = measurements;
  TGraphErrors* resolutionPlot[gain1_25];
  
  resolutionPlot[1] = new TGraphErrors(1,FDB1,Pstop_default_FDB_25_preirr,0,Pstop_default_FDB_25_preirr_err);
  resolutionPlot[3] = new TGraphErrors(1,FDB2,Pstop_RD53Apads_FDB_25,0,Pstop_RD53Apads_FDB_25_err);
  resolutionPlot[4] = new TGraphErrors(1,FDB3,Pspray_default_FDB_25,0,Pspray_default_FDB_25_err);
  resolutionPlot[5] = new TGraphErrors(1,FDB4,Pspray_RD53Apads_FDB_25,0,Pspray_RD53Apads_FDB_25_err);//,0,0_err);
  resolutionPlot[0] = new TGraphErrors(1,FTH,Pstop_default_FTH_25_irr,0,Pstop_default_FTH_25_irr_err);//Pspray_RD53Apads_FDB_25);//,0,0);
  resolutionPlot[6] = new TGraphErrors(1,FDB5,Pstop_default_FDB_25_irr,0,Pstop_default_FDB_25_irr_err);
  resolutionPlot[2] = new TGraphErrors(1,FDB6,Pstop_default_FDB_25_gain2_irr,0,Pstop_default_FDB_25_gain2_irr_err);

  resolutionPlot[1]->SetTitle(" ");
  resolutionPlot[1]->GetYaxis()->SetTitle("Resolution [#mum]");
  resolutionPlot[1]->GetXaxis()->SetTitle("Tested samples");
  resolutionPlot[1]->SetMarkerSize(2.5);
  resolutionPlot[1]->SetMarkerColor(kBlue);
  resolutionPlot[1]->SetMarkerStyle(20);
  resolutionPlot[1]->GetXaxis()->SetLimits(-1.,7.);
  resolutionPlot[1]->GetYaxis()->SetRangeUser(0.,6.);
  resolutionPlot[1]->Draw("AEP");

  resolutionPlot[3]->SetMarkerSize(2.5);
  resolutionPlot[3]->SetMarkerColor(kBlue);
  resolutionPlot[3]->SetMarkerStyle(21);
  resolutionPlot[3]->Draw("EPsame");
				      
  resolutionPlot[4]->SetMarkerSize(2.5);
  resolutionPlot[4]->SetMarkerColor(kBlue);
  resolutionPlot[4]->SetMarkerStyle(24);
  resolutionPlot[4]->Draw("EPsame");

  resolutionPlot[5]->SetMarkerSize(2.5);
  resolutionPlot[5]->SetMarkerColor(kBlue);
  resolutionPlot[5]->SetMarkerStyle(25);
  resolutionPlot[5]->Draw("EPsame");

  resolutionPlot[0]->SetMarkerSize(2.5);
  resolutionPlot[0]->SetMarkerColor(kBlue+2);
  resolutionPlot[0]->SetMarkerStyle(20);
  resolutionPlot[0]->Draw("EPsame");


  resolutionPlot[6]->SetMarkerSize(2.5);
  resolutionPlot[6]->SetMarkerColor(kBlue);
  resolutionPlot[6]->SetMarkerStyle(20);
  resolutionPlot[6]->Draw("EPsame");

  resolutionPlot[2]->SetMarkerSize(2.5);
  resolutionPlot[2]->SetMarkerColor(kBlue);
  resolutionPlot[2]->SetMarkerStyle(3);
  resolutionPlot[2]->Draw("EPsame");
  

  TLegend* legFDB2 = new TLegend(0.15,0.2,0.35,0.4);
  legFDB2->SetLineColor(0);
  legFDB2->SetTextSize(0.027);

  int sens[] = {1,3,4,5,0,6,2};
  for(int i =0; i<gain1_25; i++){
    if( i != 6) legFDB2->AddEntry(resolutionPlot[sens[i]],pitch[sens[i]]+"#mum, "+Short[sens[i]]+", "+bias[sens[i]]+"V ("+ss_irr[sens[i]]+"), thr"+thr[sens[i]], "p");
    if( i == 6) legFDB2->AddEntry(resolutionPlot[sens[i]],pitch[sens[i]]+"#mum, "+Short[sens[i]]+", "+bias[sens[i]]+"V ("+ss_irr[sens[i]]+"), thr"+thr[sens[i]]+" (gain2)", "p");

  }


  TLine *  linea    = new TLine( 3.5,0.,3.5,6.);
  linea->SetLineColor(kGray);
  linea->SetLineWidth(2);
  linea->SetLineStyle(2);
  linea->Draw("same");

  /*
  legFDB2->AddEntry(resolutionPlotFDB2,"FDB, default, p-stop", "p"); //Pstop_default_FDB_25
  legFDB2->AddEntry(resolutionPlotFDB3,"FDB, RD53A pads, p-stop", "p"); //Pstop_RD53Apads_FDB_25
  legFDB2->AddEntry(resolutionPlotFDB4,"FDB, default, p-spray", "p");//Pspray_default_FDB_25
  legFDB2->AddEntry(resolutionPlotFDB5,"FDB, RD53A pads, p-spray", "p");//Pspray_RD53Apads_FDB_25
  */

  //  legFDB2->AddEntry(resolutionPlot25,"25 #mum pitch", "p");
  //  legFDB2->AddEntry(resolutionPlot50,"50 #mum pitch", "p");
  /*
  legFDB2->AddEntry(resolutionPlotPstop,"P stop", "p");
  legFDB2->AddEntry(resolutionPlotPspray,"P spray", "p");
  */

  //legFDB2->AddEntry(resolutionPlotFTH25,"25 #mum pitch, bulk FTH", "p");

  /*
  legFDB2->AddEntry(resolutionPlotFDB25,"25 #mum pitch, bulk FDB", "p");
  */

  //legFDB2->AddEntry(resolutionPlotFDD25,"25 #mum pitch, bulk FDD", "p");
  //legFDB2->AddEntry(resolutionPlotFTH50,"50 #mum pitch, bulk FTH", "p");
  //legFDB2->AddEntry(resolutionPlotFDB50,"50 #mum pitch, bulk FDB", "p");
  //legFDB2->AddEntry(resolutionPlotFDD50,"50 #mum pitch, bulk FDD", "p");

/*
  legFDB2->AddEntry(resolutionPlotDefault,"Default", "p");
  legFDB2->AddEntry(resolutionPlotRD53,"RD53A pads", "p");
*/ 
 //legFDB2->AddEntry(resolutionPlotBdot,"Bias dot large", "p");
  //legFDB2->AddEntry(resolutionPlotBdotW,"Bias dot wiggle", "p");
  legFDB2->Draw();

  
  TString  outname = outputDir+"ResolutionSummary_25gain1_"+name;
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



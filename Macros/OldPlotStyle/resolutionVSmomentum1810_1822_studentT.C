#include <stdio.h>
#include "TCanvas.h"
#include "TF1.h"
#include <TTreeFormula.h>
#include <TMath.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TVectorT.h>
#include <TMap.h>
#include <TMath.h>
#include <TAxis.h>
#include <sstream>
#include "fileHandler.h"

#define BeamEnergies 13
bool print=true;
using namespace std;
#define dphcuts  1
#define comparisons 1


void resolutionVSmomentum1810_1822_studentT(TString func = "RMS")
{

  TString detectorA="148";
  TString detectorB="130i";
  TString detectorC="150";
  TString labelA="FDB150P_12_R4S100x25-P4_1, thr 12 ADC";
  TString labelB="FDB150P_12_R4S100x25-P1_4, 600 V, thr";
  TString labelC="FDB150P_12_R4S100x25-P1_3, thr 12 ADC";
  TString info = "angle 11.25 deg";
  
  
  double BeamEnergy[BeamEnergies];
  double BeamEnergyInverse[BeamEnergies];
  double BeamEnergyInverseSquare[BeamEnergies];
  double BeamEnergyError[BeamEnergies];
  Double_t Resolution[BeamEnergies];
  Double_t ResolutionSquare[BeamEnergies];
  Double_t ResolutionError[BeamEnergies];
  Double_t ResolutionErrorSquare[BeamEnergies];
  TString ss_BeamEnergy[BeamEnergies];
  
  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";
  TString inputfile;

  TH1F * h_res;

  Int_t run[BeamEnergies]={1819,1818,1817,1816,1815,1814,1813,1812,1811,1810,1800,1820,1822};
  TString Run[BeamEnergies];
  ostringstream strs[BeamEnergies];
  
  if(print) cout << "Iniazializating beam energies "  << endl;
  for(int i=0; i<BeamEnergies; i++)
    { 

      if(print) cout << "Run: "<< i <<" " << run[i] << endl;
      Run[i].Form("%d",run[i]);
      if(print) cout << "Run: " << Run << endl;
      

      BeamEnergy[i]=1.2+0.4*i;//GeV

      BeamEnergyInverse[i]=1/BeamEnergy[i];
      BeamEnergyInverseSquare[i]=BeamEnergyInverse[i]*BeamEnergyInverse[i];
      BeamEnergyError[i]=0; //GeV
      if(print) cout << "Beam energy " << i<< ": " << BeamEnergy[i] << " GeV" << endl;
      strs[i] << BeamEnergy[i];
      ss_BeamEnergy[i]=strs[i].str();
      if(print) cout << "Beam energy " << i<< ": " << ss_BeamEnergy[i] << " GeV" << endl;
      
    }

  MapTH1 res_map;
  int i_dphcut[dphcuts] = {12};
  TString ss_dphcut[dphcuts] = {"12"};

  TString Label[comparisons] = {"straightTracksY_isoAandCandB_straightTracksX"};
  TString Hist = "dx3_clchargeABC90evR";
  bool dphcut = false;


  
  GetHists(&res_map,BeamEnergies, dphcuts, comparisons,dphcut, Run, i_dphcut, Label, Hist, h_res);
  

  if(print) cout << "Fit:"  << endl;
  

  ofstream myfile;
  myfile.open ("/home/zoiirene/Output/TextFiles/Mscan_"+detectorB+"_"+func+"_rms.txt");
  myfile << "A " << detectorA << "\n";
  myfile << labelA << "\n";
  myfile << "B " << detectorB<< "\n";
  myfile << labelB << "\n";
  myfile << "C " << detectorC<< "\n";
  myfile << labelC << "\n";
  myfile << info << "\n";
  myfile << "Momentum(GeV) RMS(um) Error\n";
  
  for(int i=0; i<BeamEnergies; i++)
    {
      auto it2 = res_map.find(std::make_pair(Run[i]+"_"+Label[0],Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	  
	  FitTH1(it2->second, &(Resolution[i]), &(ResolutionError[i]), ss_BeamEnergy[i], detectorA, detectorB, detectorC,func );
	  if(print) cout << "Beam energy " << i<< ": " << ss_BeamEnergy[i] << " GeV -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
	  ResolutionSquare[i]=Resolution[i]*Resolution[i];
	  ResolutionErrorSquare[i]=2*ResolutionError[i]*Resolution[i];
	  myfile << ss_BeamEnergy[i] << " " << Resolution[i] << " " << ResolutionError[i] << "\n";
	}
    }

  myfile.close();


  myfile.close();

  if(print) cout << "Plotting Res vs angle"  << endl;

  DrawTGraphWithErrorDouble(BeamEnergies, BeamEnergy, Resolution, ResolutionError, Hist,"2832", Label[0], "momentumScan_rms", "RMS 95% [#mum]",0.,7.,0.,8., "momentumScan", "Beam momentum [GeV]"  );

  TCanvas *c4 = new TCanvas("c4", "resolution vs inverse beam energy squared", 1500, 900);
  gPad->SetTicks(1,1);

  TGraph* resolutionPlotInvGeorg = new TGraphErrors(BeamEnergies,BeamEnergyInverseSquare,ResolutionSquare,BeamEnergyError,ResolutionErrorSquare);

  resolutionPlotInvGeorg->SetTitle(" ");
  resolutionPlotInvGeorg->GetYaxis()->SetTitle("Resolution^{2} [#mum^{2}]");
  resolutionPlotInvGeorg->GetXaxis()->SetTitle("1/(Beam Momentum)^{2} [GeV^{-2}]");
  resolutionPlotInvGeorg->SetMarkerSize(2.5);
  resolutionPlotInvGeorg->SetLineColor(2);
  resolutionPlotInvGeorg->SetMarkerColor(2);
  resolutionPlotInvGeorg->SetMarkerStyle(20);
  resolutionPlotInvGeorg->SetLineWidth(2);
  //resolutionPlotInvGeorg->SetLineStyle(2);
  resolutionPlotInvGeorg->GetYaxis()->SetRangeUser(0.,200.);
  resolutionPlotInvGeorg->GetXaxis()->SetLimits(0.,1.);

  //  TF1 *fit1 = new TF1("fit1","pol1", 0.15, 1.255); 
  TF1 *fit3 = new TF1("fit3","pol1", 0., 0.45); 
  fit3->SetLineColor(kBlue);
  //  fit1->SetParameter(0,3);
  //fit1->SetParameter(1,15);
  fit3->SetParName(0,"#sigma_{hit}^{2}");
  fit3->SetParName(1,"#sigma_{MS}^{2}");
  
  resolutionPlotInvGeorg->Fit("fit3","R");
  resolutionPlotInvGeorg->Draw("AEP");

  //  gStyle->SetStatX(0.5);
  gStyle->SetStatX(0.85);
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
  // Tl_2.DrawLatexNDC(0.55,0.8,"#sigma_{hit} = ("+ss_fit[0]+" #pm "+ss_fit_err[0]+") #mum");
  // Tl_2.DrawLatexNDC(0.55,0.72,"#sigma_{MS} = ("+ss_fit[1]+" #pm "+ss_fit_err[1]+") #mum*GeV");
  Tl_2.DrawLatexNDC(0.55,0.6,"#sigma_{hit} = ("+ss_fit[0]+" #pm "+ss_fit_err[0]+") #mum");
  Tl_2.DrawLatexNDC(0.55,0.52,"#sigma_{MS} = ("+ss_fit[1]+" #pm "+ss_fit_err[1]+") #mum*GeV");
  
  
  TLegend* leg4 = new TLegend(0.45,0.15,0.8,0.35);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.025);
  leg4->AddEntry(resolutionPlotInvGeorg,"c"+detectorB+" "+labelB , "ep");
  leg4->AddEntry(resolutionPlotInvGeorg,"A: c"+detectorA+" "+labelA ,"");
  leg4->AddEntry(resolutionPlotInvGeorg,"C: c"+detectorC+" "+labelC ,"");
  leg4->AddEntry(resolutionPlotInvGeorg,info ,"");
  leg4->AddEntry(fit3,"#sigma_{hit}^{2}+(#sigma_{MS}/p)^{2}" , "l");
  leg4->Draw();

  TString  name = outputDir+"RMS_vs_InvEnergySquared_"+func+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC;  
  c4->SaveAs(name+".eps");
  c4->SaveAs(name+".pdf");
  c4->SaveAs(name+".png");

  
  //////////////          RESOLUTION ///////////////////////
  double freshres = 2.03;
  double freshres_err = 0.01;

  
  for(int i=0; i<BeamEnergies; i++)
    {

      if(print) cout << "Energy " << i<< ": " << ss_BeamEnergy[i] << " GeV -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      ExtractRes(&(Resolution[i]),&(ResolutionError[i]),true, freshres, freshres_err);

      if(print) cout << "Energy " << i<< ": " << ss_BeamEnergy[i] << " GeV -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      ResolutionSquare[i]=Resolution[i]*Resolution[i];
      ResolutionErrorSquare[i]=2*ResolutionError[i]*Resolution[i];
      if(print) cout << "Energy^2 " << i<< ": " << BeamEnergyInverseSquare[i] << " GeV^2 -> Resolution: " << ResolutionSquare[i] << " and res err: " << ResolutionErrorSquare[i] << endl;
      
    }

  
  DrawTGraphWithErrorDouble(BeamEnergies, BeamEnergy, Resolution, ResolutionError, Hist,"2832", Label[0], "momentumScan_res", "Resolution [#mum]",0.,7.,0.,8., "momentumScan", "Beam Momentum [GeV]"  );

  TCanvas *c42 = new TCanvas("c42", "resolution vs inverse beam energy squared", 1500, 900);
  c42->SetLeftMargin(0.1);
  c42->SetRightMargin(0.05);
  c42->SetTopMargin(0.05);
  c42->SetBottomMargin(0.2);
  
  gPad->SetTicks(1,1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TGraph* resolutionPlotInvSquare = new TGraphErrors(BeamEnergies,BeamEnergyInverseSquare,ResolutionSquare,BeamEnergyError,ResolutionErrorSquare);

  resolutionPlotInvSquare->SetTitle(" ");
  resolutionPlotInvSquare->GetYaxis()->SetTitle("Single hit resolution^{2} [#mum^{2}]");
  resolutionPlotInvSquare->GetXaxis()->SetTitle("1/(Beam Momentum)^{2} [GeV^{-2}]");
  resolutionPlotInvSquare->SetMarkerSize(2.5);
  resolutionPlotInvSquare->SetLineColor(2);
  resolutionPlotInvSquare->SetMarkerColor(2);
  resolutionPlotInvSquare->SetMarkerStyle(20);
  resolutionPlotInvSquare->SetLineWidth(2);
  //  resolutionPlotInvSquare->SetLineStyle(2);


 resolutionPlotInvSquare->GetXaxis()->SetLabelFont(42);
 resolutionPlotInvSquare->GetXaxis()->SetLabelSize(0.06);
 resolutionPlotInvSquare->GetXaxis()->SetTitleSize(0.06);
 resolutionPlotInvSquare->GetXaxis()->SetTitleOffset(1.);
 resolutionPlotInvSquare->GetXaxis()->SetTitleFont(42);

 resolutionPlotInvSquare->GetYaxis()->SetLabelFont(42);
 resolutionPlotInvSquare->GetYaxis()->SetLabelSize(0.06);
 resolutionPlotInvSquare->GetYaxis()->SetTitleSize(0.06);
 resolutionPlotInvSquare->GetYaxis()->SetTitleOffset(0.8);
 resolutionPlotInvSquare->GetYaxis()->SetTitleFont(42);
 resolutionPlotInvSquare->GetYaxis()->SetNoExponent(3);
  


  resolutionPlotInvSquare->GetYaxis()->SetRangeUser(0.,200.);
  resolutionPlotInvSquare->GetXaxis()->SetLimits(0.,0.6);

  //  TF1 *fit1 = new TF1("fit1","pol1", 0.15, 1.255); 
  TF1 *fit32 = new TF1("fit32","pol1", 0., 0.45); 
  fit32->SetLineColor(kBlue);
  //  fit1->SetParameter(0,3);
  //fit1->SetParameter(1,15);
  fit32->SetParName(0,"#sigma_{intr}^{2}");
  fit32->SetParName(1,"#sigma_{MS}^{2}");
  
  resolutionPlotInvSquare->Fit("fit32","R");
  resolutionPlotInvSquare->Draw("AEP");

  // TPaveStats * ps = (TPaveStats *)fit1->FindObject("stats");
  // ps->SetX1NDC(0.2);
  // ps->SetX2NDC(0.2);

  //gStyle->SetStatX(0.85);
  //gStyle->SetStatY(0.85);

  ostringstream strfit2[2];
  TString ss_fit2[2];
  ostringstream strfit_err2[2];
  TString ss_fit_err2[2];

  for(int i=0; i<2;i++)
    {
      strfit2[i] << setprecision(3) << sqrt(fit32->GetParameter(i));
      ss_fit2[i]=strfit2[i].str();
      strfit_err2[i] << setprecision(1) << 0.5*fit32->GetParError(i)/sqrt(fit32->GetParameter(i));
      ss_fit_err2[i]=strfit_err2[i].str();
    }
  
  TLatex Tl_22;
  Tl_22.SetTextAlign(12);
  Tl_22.SetTextSize(0.07);
  Tl_22.DrawLatexNDC(0.2,0.8,"#sigma_{intr} = ("+ss_fit2[0]+" #pm "+ss_fit_err2[0]+") #mum");

  Tl_22.SetTextSize(0.05);
  Tl_22.DrawLatexNDC(0.2,0.72,"#sigma_{MS} = ("+ss_fit2[1]+" #pm "+ss_fit_err2[1]+") #mum*GeV");
  
  
  TLegend* leg42 = new TLegend(0.35,0.25,0.8,0.45);
  leg42->SetLineColor(0);
  leg42->SetTextSize(0.04);
  // leg42->AddEntry(resolutionPlotInvSquare,"c"+detectorB+" "+labelB , "ep");
  // leg42->AddEntry(resolutionPlotInvSquare,"A: c"+detectorA+" "+labelA ,"");
  // leg42->AddEntry(resolutionPlotInvSquare,"C: c"+detectorC+" "+labelC ,"");
  //  leg42->AddEntry(resolutionPlotInvSquare,info ,"");
  leg42->AddEntry(resolutionPlotInvSquare,"#phi_{eq} = 2.1#times10^{15} cm^{-2}, PS, 600V, optimal angle" ,"ep");
  leg42->AddEntry(fit32,"#sigma_{intr}^{2}+(#sigma_{MS}/momentum)^{2}" , "l");
  leg42->Draw();
  c42->Update();
  name = outputDir+"Res_vs_InvEnergySquared_"+func+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC;  
  c42->SaveAs(name+".eps");
  c42->SaveAs(name+".pdf");
  c42->SaveAs(name+".png");
  c42->SaveAs(name+".root");


  TCanvas *c4b = new TCanvas("c4b", "resolution vs beam momentum", 1500, 900);
  c4b->SetLeftMargin(0.1);
  c4b->SetRightMargin(0.05);
  c4b->SetTopMargin(0.05);
  c4b->SetBottomMargin(0.2);
  
  gPad->SetTicks(1,1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TGraph* resolutionPlotMyParam = new TGraphErrors(BeamEnergies,BeamEnergy,Resolution,BeamEnergyError,ResolutionError);

  resolutionPlotMyParam->SetTitle(" ");
  resolutionPlotMyParam->GetYaxis()->SetTitle("Resolution [#mum]");
  resolutionPlotMyParam->GetXaxis()->SetTitle("Beam Momentum [GeV]");
  resolutionPlotMyParam->SetMarkerSize(2.5);
  resolutionPlotMyParam->SetLineColor(2);
  resolutionPlotMyParam->SetMarkerColor(2);
  resolutionPlotMyParam->SetMarkerStyle(20);
  resolutionPlotMyParam->SetLineWidth(2);
  //  resolutionPlotMyParam->SetLineStyle(2);


 resolutionPlotMyParam->GetXaxis()->SetLabelFont(42);
 resolutionPlotMyParam->GetXaxis()->SetLabelSize(0.07);
 resolutionPlotMyParam->GetXaxis()->SetTitleSize(0.07);
 resolutionPlotMyParam->GetXaxis()->SetTitleOffset(0.9);
 resolutionPlotMyParam->GetXaxis()->SetTitleFont(42);

 resolutionPlotMyParam->GetYaxis()->SetLabelFont(42);
 resolutionPlotMyParam->GetYaxis()->SetLabelSize(0.07);
 resolutionPlotMyParam->GetYaxis()->SetTitleSize(0.07);
 resolutionPlotMyParam->GetYaxis()->SetTitleOffset(0.7);
 resolutionPlotMyParam->GetYaxis()->SetTitleFont(42);
 resolutionPlotMyParam->GetYaxis()->SetNoExponent(3);
  


  resolutionPlotMyParam->GetYaxis()->SetRangeUser(0.,8.);
  resolutionPlotMyParam->GetXaxis()->SetLimits(0.,7.);

  //  TF1 *fit1 = new TF1("fit1","pol1", 0.15, 1.255); 
  TF1 *fnotsq = new TF1("fnotsq","TMath::Sqrt([0]+[1]/(x*x))", 0., 7.); 
  fnotsq->SetLineColor(kBlue);
  //  fit1->SetParameter(0,3);
  //fit1->SetParameter(1,15);
  fnotsq->SetParName(0,"#sigma_{hit}^{2}");
  fnotsq->SetParName(1,"#sigma_{MS}^{2}");
  fnotsq->SetParameter(0,fit32->GetParameter(0));
  fnotsq->SetParameter(1,fit32->GetParameter(1));
  

  resolutionPlotMyParam->Draw("AEP");
  fnotsq->Draw("same");

  c4b->Update();
  name = outputDir+"Res_vs_Momentum_myparam_"+func+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC;  
  c4b->SaveAs(name+".eps");
  c4b->SaveAs(name+".pdf");
  c4b->SaveAs(name+".png");
  c4b->SaveAs(name+".root");


  


  
}//resolution 


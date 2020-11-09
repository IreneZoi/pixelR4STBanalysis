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

#include "tdrstyle.C"
#include "CMS_lumi.C"


#define BeamEnergies 8
bool print=true;
using namespace std;
#define dphcuts  1
#define comparisons 1

void TDR();
void TDR2(TCanvas * c_all, int period=0, int pos= 11);


void resolutionVSmomentum2832_2840_Tree(TString func = "RMSself")
{

  TString detectorA="148";
  TString detectorB="163";
  TString detectorC="150";
  TString labelA="FDB150P_12_R4S100x25-P4_1, thr 12 ADC";
  TString labelB="FDB150Y_2_R4S100x25-Y6_1, 120 V, thr 29 ADC";
  TString labelC="FDB150P_12_R4S100x25-P1_3, thr 12 ADC";
  TString info = "angle 11.25 deg";
  
  
  double BeamEnergy[BeamEnergies];
  double BeamEnergyInverse[BeamEnergies];
  double BeamEnergyInverseSquare[BeamEnergies];
  double BeamEnergyError[BeamEnergies];
  double BeamEnergyInverseSquareError[BeamEnergies];
  
  Double_t Resolution[BeamEnergies];
  Double_t ResolutionSquare[BeamEnergies];
  Double_t ResolutionError[BeamEnergies];
  Double_t ResolutionErrorSquare[BeamEnergies];
  Double_t Percentage[BeamEnergies];
  TString ss_BeamEnergy[BeamEnergies];
  
  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";
  TString inputfile;
  TString label = "dycut_A11C12";//"closest2_A11C12";
  TH1F * h_res;

  //Int_t run[BeamEnergies]={2837,2836,2835,2834,2833,2840,2838};
  Int_t run[BeamEnergies]={2837,2836,2835,2834,2833,2840,2832,2838};
  TString Run[BeamEnergies];
  ostringstream strs[BeamEnergies];
  TString pitch = "25";
  if(print) cout << "Iniazializating beam energies "  << endl;
  for(int i=0; i<BeamEnergies; i++)
    { 

      if(print) cout << "Run: "<< i <<" " << run[i] << endl;
      Run[i].Form("%d",run[i]);
      if(print) cout << "Run: " << Run[i] << endl;


      BeamEnergy[i]=1.6+0.8*i;//GeV
      if(i==5) BeamEnergy[i]=5.4;//GeV
  
      if(i==6) BeamEnergy[i]=5.6;//GeV
      if(i==7) BeamEnergy[i]=6.0;//GeV      
  
      //if(i==6) BeamEnergy[i]=6.;//GeV
      BeamEnergyInverse[i]=1/BeamEnergy[i];
      BeamEnergyInverseSquare[i]=BeamEnergyInverse[i]*BeamEnergyInverse[i];
      BeamEnergyError[i]=0.158; //GeV https://www.sciencedirect.com/science/article/pii/S0168900218317868?via%3Dihub#sec7.3
      BeamEnergyInverseSquareError[i]=BeamEnergyError[i]*2*TMath::Power(BeamEnergyInverse[i],3);
      

      if(print) cout << "Beam energy " << i<< ": " << BeamEnergy[i] << " GeV" << endl;
      strs[i] << BeamEnergy[i];
      ss_BeamEnergy[i]=strs[i].str();
      if(print) cout << "Beam energy " << i<< ": " << ss_BeamEnergy[i] << " GeV" << endl;
      
    }

  MapTH1 res_map;
  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  int i_dphcut[dphcuts] = {14}; // {29};
  TString ss_dphcut[dphcuts] = {"14"}; //{"29"};
  
  
  TString Label[comparisons] = {"straightTracksY_isoAandCandB_straightTracksX"};
  TString Hist = "dx3_clphABC90evR";
  bool dphcut = true;

  TH1F * hdx3_clchargeABC90evR; // = new TH1I("dx3_clchargeABC90evR ", "triplet dx_clchargeABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
  res_map =    GetCheckHists(&res_map, BeamEnergies, dphcuts,dphcut, Run,ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR, true,label);

  for ( it = res_map.begin(); it != res_map.end(); it++ )       {
    cout << it->first.first << " " << it->first.second << " " <<  it->second->GetEntries() << endl;
  }
  
  
   if(print) cout << "Fit:"  << endl;
  

  ofstream myfile;
  myfile.open ("/home/zoiirene/Output/TextFiles/Mscan_"+detectorB+"_"+func+"_rmsTree.txt");
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
      auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[0]+"_"+pitch,Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	  
	  FitTH1(it2->second, &(Resolution[i]), &(ResolutionError[i]), ss_BeamEnergy[i], detectorA, detectorB, detectorC,func,&(Percentage[i]) );
	  if(print) cout << "Beam energy " << i<< ": " << ss_BeamEnergy[i] << " GeV -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
	  ResolutionSquare[i]=Resolution[i]*Resolution[i];
	  ResolutionErrorSquare[i]=2*ResolutionError[i]*Resolution[i];
	  myfile << ss_BeamEnergy[i] << " " << Resolution[i] << " " << ResolutionError[i] << " perc "<< Percentage[i] << "\n";
	}
    }

  myfile.close();


  myfile.close();

  if(print) cout << "Plotting Res vs angle"  << endl;

  DrawTGraphWithErrorDouble(BeamEnergies, BeamEnergy, Resolution, ResolutionError, Hist,"2832", Label[0], "momentumScan_rmsTree", "RMS 95% [#mum]",0.,7.,0.,8., "momentumScan", "Beam momentum [GeV]"  );

  TCanvas *c4 = new TCanvas("c4", "resolution vs inverse beam energy squared", 1500, 900);
  gPad->SetTicks(1,1);

  TGraph* resolutionPlotInvGeorg = new TGraphErrors(BeamEnergies,BeamEnergyInverseSquare,ResolutionSquare,BeamEnergyInverseSquareError,ResolutionErrorSquare);

  resolutionPlotInvGeorg->SetTitle(" ");
  resolutionPlotInvGeorg->GetYaxis()->SetTitle("Single hit resolution^{2} [#mum^{2}]");
  resolutionPlotInvGeorg->GetXaxis()->SetTitle("1/(Beam Momentum)^{2} [GeV^{-2}]");
  resolutionPlotInvGeorg->SetMarkerSize(2.5);
  resolutionPlotInvGeorg->SetLineColor(2);
  resolutionPlotInvGeorg->SetMarkerColor(2);
  resolutionPlotInvGeorg->SetMarkerStyle(20);
  resolutionPlotInvGeorg->SetLineWidth(2);
  //resolutionPlotInvGeorg->SetLineStyle(2);
  resolutionPlotInvGeorg->GetYaxis()->SetRangeUser(0.,100.);
  resolutionPlotInvGeorg->GetXaxis()->SetLimits(0.,1.);

  //  TF1 *fit1 = new TF1("fit1","pol1", 0.15, 1.255); 
  TF1 *fit3 = new TF1("fit3","pol1", 0., 0.45); 
  fit3->SetLineColor(kBlue);
  //  fit1->SetParameter(0,3);
  //fit1->SetParameter(1,15);
  fit3->SetParName(0,"#sigma_{int}^{2}");
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
  Tl_2.DrawLatexNDC(0.55,0.6,"#sigma_{intr} = ("+ss_fit[0]+" #pm "+ss_fit_err[0]+") #mum");
  Tl_2.DrawLatexNDC(0.55,0.52,"#sigma_{MS} = ("+ss_fit[1]+" #pm "+ss_fit_err[1]+") #mum*GeV");
  
  
  TLegend* leg4 = new TLegend(0.45,0.15,0.8,0.35);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.025);
  leg4->AddEntry(resolutionPlotInvGeorg,"c"+detectorB+" "+labelB , "ep");
  leg4->AddEntry(resolutionPlotInvGeorg,"A: c"+detectorA+" "+labelA ,"");
  leg4->AddEntry(resolutionPlotInvGeorg,"C: c"+detectorC+" "+labelC ,"");
  leg4->AddEntry(resolutionPlotInvGeorg,info ,"");
  leg4->AddEntry(fit3,"#sigma_{intr}^{2}+(#sigma_{MS}/p)^{2}" , "l");
  leg4->Draw();

  TString  name = outputDir+"RMS_vs_InvEnergySquared_"+func+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC;  
  c4->SaveAs(name+".eps");
  c4->SaveAs(name+".pdf");
  c4->SaveAs(name+".png");




  ofstream myfile2;
  myfile2.open ("/home/zoiirene/Output/TextFiles/Mscan_"+detectorB+"_"+func+"_resTree.txt");
  myfile2 << "A " << detectorA << "\n";
  myfile2 << labelA << "\n";
  myfile2 << "B " << detectorB<< "\n";
  myfile2 << labelB << "\n";
  myfile2 << "C " << detectorC<< "\n";
  myfile2 << labelC << "\n";
  myfile2 << info << "\n";
  myfile2 << "Momentum(GeV) RES(um) Error\n";
  
  //////////////          RESOLUTION ///////////////////////
  for(int i=0; i<BeamEnergies; i++)
    {

      if(print) cout << "Energy " << i<< ": " << ss_BeamEnergy[i] << " GeV -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      ExtractRes(&(Resolution[i]),&(ResolutionError[i]));
      if(print) cout << "Energy " << i<< ": " << ss_BeamEnergy[i] << " GeV -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      myfile2 << ss_BeamEnergy[i] << " " << Resolution[i] << " " << ResolutionError[i] << "\n";

      ResolutionSquare[i]=Resolution[i]*Resolution[i];
      ResolutionErrorSquare[i]=2*ResolutionError[i]*Resolution[i];

      if(print) cout << "Energy^2 " << i<< ": " << BeamEnergyInverseSquare[i] << " GeV^2 -> Resolution: " << ResolutionSquare[i] << " and res err: " << ResolutionErrorSquare[i] << endl;
      
    }
  myfile2.close();  
  myfile2.close();

  
  DrawTGraphWithErrorDouble(BeamEnergies, BeamEnergy, Resolution, ResolutionError, Hist,"2832", Label[0], "momentumScan_resTree", "Resolution [#mum]",0.,7.,0.,8., "momentumScan", "Beam Momentum [GeV]"  );

  ///////////////////////////////////////////////////
  
  //  TDR();
  TCanvas *c42 = new TCanvas("c42", "resolution vs inverse beam energy squared", 1000, 700); //1500,900
  c42->SetLeftMargin(0.1);
  c42->SetRightMargin(0.05);
  c42->SetTopMargin(0.05);
  c42->SetBottomMargin(0.2);
  
  gPad->SetTicks(1,1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TGraph* resolutionPlotInvSquare = new TGraphErrors(BeamEnergies,BeamEnergyInverseSquare,ResolutionSquare,BeamEnergyInverseSquareError,ResolutionErrorSquare);

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
  


  resolutionPlotInvSquare->GetYaxis()->SetRangeUser(0.,60.);
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
      strfit2[i] << fixed <<setprecision(2) << sqrt(fit32->GetParameter(i));
      ss_fit2[i]=strfit2[i].str();
      strfit_err2[i] << setprecision(1) << 0.5*fit32->GetParError(i)/sqrt(fit32->GetParameter(i));
      ss_fit_err2[i]=strfit_err2[i].str();
    }
  
  TLatex Tl_22;
  Tl_22.SetTextAlign(12);
  Tl_22.SetTextSize(0.05);
  Tl_22.DrawLatexNDC(0.2,0.85,"#sigma_{intr} = ("+ss_fit2[0]+" #pm "+ss_fit_err2[0]+") #mum");

  Tl_22.SetTextSize(0.05);
  Tl_22.DrawLatexNDC(0.2,0.77,"#sigma_{MS} = ("+ss_fit2[1]+" #pm "+ss_fit_err2[1]+") #mum*GeV");
  
  
  TLegend* leg42 = new TLegend(0.38,0.25,0.8,0.45);
  leg42->SetLineColor(0);
  leg42->SetTextSize(0.04);
  leg42->SetBorderSize(0);
  // leg42->AddEntry(resolutionPlotInvSquare,"c"+detectorB+" "+labelB , "ep");
  // leg42->AddEntry(resolutionPlotInvSquare,"A: c"+detectorA+" "+labelA ,"");
  // leg42->AddEntry(resolutionPlotInvSquare,"C: c"+detectorC+" "+labelC ,"");
  //  leg42->AddEntry(resolutionPlotInvSquare,info ,"");
  leg42->AddEntry(resolutionPlotInvSquare,"Non-irradiated, 120V, optimal angle" ,"ep");
  leg42->AddEntry(fit32,"#sigma_{intr}^{2}+(#sigma_{MS}/momentum)^{2}" , "l");
  leg42->Draw();
  c42->Update();
  //TDR2(c42);
  name = outputDir+"ResTree_vs_InvEnergySquared_"+func+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC;  
  c42->SaveAs(name+".eps");
  c42->SaveAs(name+".pdf");
  c42->SaveAs(name+".png");
  c42->SaveAs(name+".root");


  /*
  //par from lin fit
  //TDR();
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
  


  resolutionPlotMyParam->GetYaxis()->SetRangeUser(0.,10.);
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
  leg42->Draw();
  c4b->Update();
  //  TDR2(c4b);
  name = outputDir+"ResTree_vs_Momentum_myparam_"+func+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC;  
  c4b->SaveAs(name+".eps");
  c4b->SaveAs(name+".pdf");
  c4b->SaveAs(name+".png");
  c4b->SaveAs(name+".root");
  */

  
  //hyperbolic (?) fit
  //TDR();
  TCanvas *cf = new TCanvas("cf", "resolution vs beam momentum", 1500, 900);
  cf->SetLeftMargin(0.1);
  cf->SetRightMargin(0.05);
  cf->SetTopMargin(0.05);
  cf->SetBottomMargin(0.2);
  
  gPad->SetTicks(1,1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TGraph* resolutionPlotMyParam = new TGraphErrors(BeamEnergies,BeamEnergy,Resolution,BeamEnergyError,ResolutionError);
  
  resolutionPlotMyParam->SetTitle(" ");
  resolutionPlotMyParam->GetYaxis()->SetTitle("Single hit resolution [#mum]");
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
  


  resolutionPlotMyParam->GetYaxis()->SetRangeUser(0.,10.);
  resolutionPlotMyParam->GetXaxis()->SetLimits(0.,7.);

  cout << " new fit! " << endl;
  //  TF1 *fit1 = new TF1("fit1","pol1", 0.15, 1.255); 
  TF1 *fnotsq = new TF1("fnotsq","TMath::Sqrt([0]+[1]/(x*x))", 0., 7.); 
  fnotsq->SetLineColor(kBlue);
  fnotsq->SetParName(0,"#sigma_{hit}^{2}");
  fnotsq->SetParName(1,"#sigma_{MS}^{2}");
  fnotsq->SetParLimits(0,0,20);
  fnotsq->SetParLimits(1,0,200);
  
  //fnotsq->SetParameter(0,fit32->GetParameter(0));
  //fnotsq->SetParameter(1,fit32->GetParameter(1));
  resolutionPlotMyParam->Fit("fnotsq");

  resolutionPlotMyParam->Draw("AEP");
  fnotsq->Draw("same");
  ostringstream strfit4[2];
  TString ss_fit4[2];
  ostringstream strfit_err4[2];
  TString ss_fit_err4[2];
  
  TF1 *fit = resolutionPlotMyParam->GetFunction("fnotsq");
  Double_t chi2 = fit->GetChisquare();
  Double_t ndof = fit->GetNDF();
  cout << " chi2 "<< chi2 << " / " << " ndof "<< ndof << " = "<< chi2/ndof << endl;

  

  for(int i=0; i<2;i++)
    {
      strfit4[i] << fixed <<setprecision(2) << sqrt(fnotsq->GetParameter(i));
      ss_fit4[i]=strfit4[i].str();
      strfit_err4[i] << setprecision(1) << 0.5*fnotsq->GetParError(i)/sqrt(fnotsq->GetParameter(i));
      ss_fit_err4[i]=strfit_err4[i].str();
    }

  TLatex Tl_44;
  Tl_44.SetTextAlign(12);
  Tl_44.SetTextSize(0.05);
  Tl_44.DrawLatexNDC(0.15,0.35,"#sigma_{intr} = ("+ss_fit4[0]+" #pm "+ss_fit_err4[0]+") #mum");

  Tl_44.SetTextSize(0.05);
  Tl_44.DrawLatexNDC(0.15,0.27,"#sigma_{MS} = ("+ss_fit4[1]+" #pm "+ss_fit_err4[1]+") #mum*GeV");

  TLegend* leg42f = new TLegend(0.45,0.65,0.85,0.8);
  leg42f->SetLineColor(0);
  leg42f->SetTextSize(0.04);
  leg42f->SetBorderSize(0);
  // leg42f->AddEntry(resolutionPlotInvSquare,"c"+detectorB+" "+labelB , "ep");
  // leg42f->AddEntry(resolutionPlotInvSquare,"A: c"+detectorA+" "+labelA ,"");
  // leg42f->AddEntry(resolutionPlotInvSquare,"C: c"+detectorC+" "+labelC ,"");
  //  leg42f->AddEntry(resolutionPlotInvSquare,info ,"");
  leg42f->AddEntry(resolutionPlotInvSquare,"Non-irradiated, 120V, optimal angle" ,"ep");
  leg42f->AddEntry(fit32,"#sqrt{#sigma_{intr}^{2}+(#sigma_{MS}/momentum)^{2}}" , "l");
  leg42f->Draw();
  

  cf->Update();
  //  TDR2(cf);
  name = outputDir+"ResTree_vs_Momentum_fit_"+func+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC;  
  cf->SaveAs(name+".eps");
  cf->SaveAs(name+".pdf");
  cf->SaveAs(name+".png");
  cf->SaveAs(name+".root");



  

  


  
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


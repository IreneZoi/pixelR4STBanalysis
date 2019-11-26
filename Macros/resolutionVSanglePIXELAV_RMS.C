#include <stdio.h>
#include "TCanvas.h"
#include "TF1.h"
#include <TTreeFormula.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TVectorT.h>
#include <TMap.h>
#include <TAxis.h>
#include <sstream>
#include "fileHandler.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"


#define Angles 29
bool print=true;
#define dphcuts  1
#define comparisons 1

using namespace std;
TString outputDir = "/home/zoiirene/Output/Plots/";

void TDR();
void TDR2(TCanvas * c_all);



void resolutionVSanglePIXELAV_RMS(TString thr ="700", TString function = "RMS")
{

  Double_t Angle[Angles] = {-5,-3.75,-2.5,-1.25,0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,13.75,15,16.25,17.5,18.75,20,21.25,22.5,23.75,25,26.25,27.5,28.75,30 };
  TString Angle0[Angles] = {"-005","-003","-002","-001","0000","0001","0002","0003","0005","0006","0007","0008","0010","0011","0012","0013","0015","0016","0017","0018","0020","0021","0022","0023","0025","0026","0027","0028","0030"};
  Double_t AngleError[Angles];
  Double_t TanAngle[Angles];
  Double_t TanAngleError[Angles];
  Double_t Resolution[Angles];
  Double_t RMS[Angles];
  //  Double_t Ncol[Angles];
  //Double_t NcolError[Angles];
  Double_t ResolutionError[Angles];
  Double_t RMSError[Angles];
  TString ss_Angle[Angles];
  
  TH1F * h_res;//[Angles];
  TString pitch[Angles];
  

  
  if(print) cout << "Getting files "  << endl;

  for(int i=0; i<Angles; i++)
    {
      pitch[i] = "25";
      ss_Angle[i].Form("%d",i);
      TanAngle[i]=TMath::Tan(Angle[i]*TMath::Pi()/((Double_t)180.));
      cout << " Angle " <<  Angle[i] << " TanAngle " <<  TanAngle[i] << endl;
      AngleError[i]=1; //grad
      TanAngleError[i]=(1+TMath::Power(TMath::Tan(AngleError[i]*TMath::Pi()/((Double_t) 180.)),2))*AngleError[i]*TMath::Pi()/((Double_t) 180.); //grad
      
      if(print) cout << "Angle " << i<< ": " << Angle[i] << " grads" << endl;
      
   
}



  MapTH1 res_map;
  //hdycq0
  TString Hist = "landau90/hdy";

  GetPixelavHists(&res_map,Angles, Angle0, pitch,Hist, h_res);



  
  

  if(print) cout << "Fit:"  << endl;
  
  ofstream myfile;
  myfile.open ("/home/zoiirene/Output/TextFiles/Ascan_pixelav_"+function+"_"+thr+".txt");
  myfile << "Non irradiated \n";
  myfile << "Angle Resolution(um) Error\n";

  for(int i=0; i<Angles; i++)
    {


      auto it2 = res_map.find(std::make_pair(Angle0[i]+"_"+pitch[i],Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
	  
	  FitTH1(it2->second, &(RMS[i]), &(RMSError[i]), ss_Angle[i], "pixelav" , "nonirr", "2um", function );
	  RMS[i]/=1000.;
	  RMSError[i]/=1000.;
	  if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> Resolution: " << RMS[i] << " and res err: " << RMSError[i] << endl;
	  myfile << Angle[i] << " "  << RMS[i] << " " << RMSError[i] << "\n";
	}
    }

 

  myfile.close();

  if(print) cout << "Plotting Res vs angle"  << endl;

  Hist = "landau90_hdy";
  

  DrawTGraphWithErrorDouble(Angles, Angle, RMS, RMSError, Hist,"pixelav_"+thr, "nonirr", "angleScan_rms", "Resolution (RMS 95%) [#mum]",-7.,29.,0.,8., "angleScan", "Angle [degrees]"  );


  ///////
  Hist = "landau90/hnrow";
  
  GetPixelavHists(&res_map,Angles, Angle0, pitch,Hist, h_res);
  
  ofstream myfile3;
  myfile3.open ("/home/zoiirene/Output/TextFiles/Ascan_clsizeB_pixelav_"+function+"_"+thr+".txt");
  myfile3 << "Non irradiated \n";
  myfile3 << "Angle MeanNrowB Error\n";

  Double_t NrowB[Angles];
  Double_t NrowBerror[Angles];
  for(int i=0; i<Angles; i++)
    {
      auto it2 = res_map.find(std::make_pair(Angle0[i]+"_"+pitch[i],Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;



	  NrowB[i] = it2->second->GetMean();
	  NrowBerror[i] = it2->second->GetMeanError();

	  if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> nrowB: " << NrowB[i] << " and err: " << NrowBerror[i] << endl;
	  myfile3 << " " << Angle[i] << " " << NrowB[i]<< " " << NrowBerror[i] << "\n";
	}

    }
  myfile3.close();

  Hist = "landau90_hnrow";
  
  //  DrawTGraphWithErrorDouble(Angles, TanAngle, NrowB, NrowBerror, Hist,"pixelav", Label[0], "angleScan_nrowB_vs_tanangle", "nrowB",-7.,29.,0.,5., "angleScan", "Tan angle"  );
  DrawTGraphWithErrorDouble(Angles, Angle, NrowB, NrowBerror, Hist,"pixelav_"+thr, "nonirr", "angleScan_nrowB_vs_angle", "nrow",-7.,29.,0.,5., "angleScan", "Angle [degrees]"  );
  //DrawTGraphWithErrorDouble(Angles, Angle, RMS, RMSError, Hist,"pixelav", "nonirr", "angleScan_rms", "Resolution (RMS 95%) [#mum]",-7.,29.,0.,8., "angleScan", "Angle [degrees]"  );
  
  

  /*
  Hist = "nrowB";
  GetHists(&res_map,Angles, dphcuts, comparisons,dphcut, Run, i_dphcut, Label, Hist, h_res);

  Double_t NrowB[Angles];
  Double_t NrowBerror[Angles];
  for(int i=0; i<Angles; i++)
    {
      auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[0]+"_"+Label[0],Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;



	  NrowB[i] = it2->second->GetMean();
	  NrowBerror[i] = it2->second->GetMeanError();

	  if(print) cout << "Angle " << i<< ": " << ss_Angle[i] << " degrees -> nrowB: " << NrowB[i] << " and err: " << NrowBerror[i] << endl;

	}

    }
  DrawTGraphWithErrorDouble(Angles, TanAngle, NrowB, NrowBerror, Hist,"2789", Label[0], "angleScan_nrowB_vs_tanangle", "nrowB",-7.,29.,0.,5., "angleScan", "Tan angle"  );
  DrawTGraphWithErrorDouble(Angles, Angle, NrowB, NrowBerror, Hist,"2789", Label[0], "angleScan_nrowB_vs_angle", "nrowB",-7.,29.,0.,5., "angleScan", "Angle [degrees]"  );


  Hist = "clsizeB";
  GetHists(&res_map,Angles, dphcuts, comparisons,dphcut, Run, i_dphcut, Label, Hist, h_res);

  for(int i=0; i<Angles; i++)
    {
      auto it2 = res_map.find(std::make_pair(Run[i]+"_"+ss_dphcut[0]+"_"+Label[0],Hist));
      if(it2  != res_map.end())
	{
	  if(print)               cout << " found map " << endl;
	  if(print)               cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;

	  DrawHist(it2->second, Run[i]+"_"+ss_Angle[i], ss_dphcut[0], Label[0], Hist,0.,10.,0., 50000.)  ;

	}

    }
  */
  
  
  /*
  //*******************size plot
  TCanvas *c3 = new TCanvas("c3", "Mean Number of columns vs angle", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TGraph* NcolPlot = new TGraphErrors(Angles,Angle,Ncol,AngleError,NcolError);

  NcolPlot->SetTitle(" ");
  NcolPlot->GetYaxis()->SetTitle("Mean Number of columns");
  NcolPlot->GetXaxis()->SetTitle("Angle [degrees]");
  NcolPlot->SetMarkerSize(2.5);
  NcolPlot->SetLineColor(kBlue);
  NcolPlot->SetMarkerColor(kBlue);
  NcolPlot->SetMarkerStyle(20);
  NcolPlot->SetLineWidth(2);
  NcolPlot->SetLineStyle(2);
  NcolPlot->GetXaxis()->SetRangeUser(0.,29.);
  NcolPlot->GetYaxis()->SetRangeUser(0.,3.);
  NcolPlot->Draw("AEPL");

  TLegend* leg3 = new TLegend(0.65,0.4,0.75,0.5);
  leg3->SetLineColor(0);
  leg3->SetTextSize(0.03);
  leg3->AddEntry(NcolPlot,"DUT "+detectorB , "ep");
  leg3->Draw();

  
  c3->SaveAs(outputDir+"Ncol_vs_angle_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps");
  c3->SaveAs(outputDir+"Ncol_vs_angle_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".root");

  //*******************size plot vs tang

  for(int i=0; i<Angles; i++)
    cout << " Angle " <<  Angle[i] << " TanAngle " <<  TanAngle[i] << endl;
       

  TCanvas *c5 = new TCanvas("c5", "Number of columns vs tan angle", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TGraph* NcolTanPlot = new TGraphErrors(Angles,TanAngle,Ncol,TanAngleError,NcolError);

  NcolTanPlot->SetTitle(" ");
  NcolTanPlot->GetYaxis()->SetTitle("Mean Number of columns");
  NcolTanPlot->GetXaxis()->SetTitle("Tan(Angle)");
  NcolTanPlot->SetMarkerSize(2.5);
  NcolTanPlot->SetLineColor(kBlue);
  NcolTanPlot->SetMarkerColor(kBlue);
  NcolTanPlot->SetMarkerStyle(20);
  NcolTanPlot->SetLineWidth(2);
  NcolTanPlot->SetLineStyle(2);
  //  NcolTanPlot->GetXaxis()->SetRangeUser(0.017,0.018);
  NcolTanPlot->GetXaxis()->SetLimits(0.,0.6);
  NcolTanPlot->GetYaxis()->SetRangeUser(0.,3.);
  NcolTanPlot->Draw("AEPL");
  NcolTanPlot->Fit("pol1");
  NcolTanPlot->Draw("ALPEsame");

  TLegend* leg4 = new TLegend(0.65,0.4,0.75,0.5);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.03);
  leg4->AddEntry(NcolTanPlot,"DUT "+detectorB , "ep");
  leg4->Draw();

  
  c5->SaveAs(outputDir+"Ncol_vs_tanangle_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps");
  c5->SaveAs(outputDir+"Ncol_vs_tanangle_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".root");

  /*
  TCanvas *c = new TCanvas("c", "resolution vs inverse beam energy", 1500, 900);
  gPad->SetTicks(1,1);
  TGraph* resolutionPlotInv = new TGraphErrors(Angles,BeamEnergyInverse,Resolution,BeamEnergyError,ResolutionError);

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

  c->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergy_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".png");  
  c->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergy_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps");
  c->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergy_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".pdf");


  TCanvas *c4 = new TCanvas("c4", "resolution vs inverse beam energy squared", 1500, 900);
  gPad->SetTicks(1,1);

  TGraph* resolutionPlotInvGeorg = new TGraphErrors(Angles,BeamEnergyInverseSquare,ResolutionSquare,BeamEnergyError,ResolutionErrorSquare);

  resolutionPlotInvGeorg->SetTitle(" ");
  resolutionPlotInvGeorg->GetYaxis()->SetTitle("Resolution^{2} [#mum^{2}]");
  resolutionPlotInvGeorg->GetXaxis()->SetTitle("1/(Beam Energy)^{2} [GeV^{-2}]");
  resolutionPlotInvGeorg->SetMarkerSize(2.5);
  resolutionPlotInvGeorg->SetLineColor(2);
  resolutionPlotInvGeorg->SetMarkerColor(2);
  resolutionPlotInvGeorg->SetMarkerStyle(20);
  resolutionPlotInvGeorg->SetLineWidth(2);
  resolutionPlotInvGeorg->SetLineStyle(2);
  //  resolutionPlotInvGeorg->GetYaxis()->SetRangeUser(0.,40.);
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
  Tl_2.SetTextSize(0.04);
  Tl_2.DrawLatexNDC(0.15,0.6,"#sigma_{hit} = ("+ss_fit[0]+" #pm "+ss_fit_err[0]+") #mum");
  Tl_2.DrawLatexNDC(0.15,0.52,"#sigma_{MS} = ("+ss_fit[1]+" #pm "+ss_fit_err[1]+") #mum*GeV");
  
  
  TLegend* leg4 = new TLegend(0.75,0.4,0.85,0.6);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.03);
  leg4->AddEntry(resolutionPlotInvGeorg,"DUT "+detectorB , "ep");
  leg4->AddEntry(fit3,"#sigma_{hit}^{2}+(#sigma_{MS}/p)^{2}" , "l");
  leg4->Draw();

  
  c4->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergySquared_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".png");  
  c4->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergySquared_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps");
  c4->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergySquared_newFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".pdf");

  

  TCanvas *c3 = new TCanvas("c3", "resolution vs inverse beam energy old", 1500, 900);
  gPad->SetTicks(1,1);

  TGraph* resolutionPlotInvOLD = new TGraphErrors(Angles,BeamEnergyInverse,Resolution,BeamEnergyError,ResolutionError);

  resolutionPlotInvOLD->SetTitle(" ");
  resolutionPlotInvOLD->GetYaxis()->SetTitle("Resolution [#mum]");
  resolutionPlotInvOLD->GetXaxis()->SetTitle("1/Beam Energy [GeV^{-1}]");
  resolutionPlotInvOLD->SetMarkerSize(2.5);
  resolutionPlotInvOLD->SetLineColor(2);
  resolutionPlotInvOLD->SetMarkerColor(2);
  resolutionPlotInvOLD->SetMarkerStyle(20);
  resolutionPlotInvOLD->SetLineWidth(2);
  resolutionPlotInvOLD->SetLineStyle(2);
  resolutionPlotInvOLD->GetYaxis()->SetRangeUser(0.,20.);
  resolutionPlotInvOLD->GetXaxis()->SetLimits(0.,1.6);

  TF1 *fit2 = new TF1("fit2","pol1", 0.15, 1.255); 
  //  TF1 *fit1 = new TF1("fit1","sqrt([0]*[0]+([1]*x)*([1]*x))", 0., 0.8); 
  fit2->SetLineColor(kBlue);
  // fit1->SetParameter(0,3);
  // fit1->SetParameter(1,15);
  // fit1->SetParameterName(0,"#sigma_{hit}");
  // fit1->SetParameterName(1,"#sigma_{MS}");
  
  resolutionPlotInvOLD->Fit("fit2");
  resolutionPlotInvOLD->Draw("AEP");

  // TPaveStats * ps = (TPaveStats *)fit1->FindObject("stats");
  // ps->SetX1NDC(0.2);
  // ps->SetX2NDC(0.2);
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.85);
  
  TLegend* leg3 = new TLegend(0.65,0.4,0.75,0.6);
  leg3->SetLineColor(0);
  leg3->SetTextSize(0.03);
  leg3->AddEntry(resolutionPlotInvOLD,"DUT "+detectorB , "ep");
  leg3->AddEntry(fit2,"p_{0}+x*p_{1}" , "l");
  //  leg->AddEntry(fit1,"#sqrt{#sigma_{hit}^{2}+(#sigma_{MS}/p)^{2}}" , "l");
  leg3->Draw();

  
  c3->SaveAs("/home/irene/Documents/CMS_material/DESY_TB/Plot/Res_vs_InvEnergy_oldFit_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps");

  */

  
}//resolution 

void TDR()
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  extraText  = "Work in progress";  // default extra text is "Preliminary"
  lumi_sqrtS = "";
  bool drawLogo = true;
  int H_ref = 600;
  int W_ref = 600;
  float T = 0.07*H_ref;//0.07
  float B = 0.11*H_ref;//0.12
  float L = 0.12*W_ref;
  float R = 0.01*W_ref;


}

void TDR2(TCanvas * c_all)
{
  CMS_lumi( c_all, 0, 11);
  //  CMS_lumi( c_all, 0, 11);
  c_all->Update();
  c_all->RedrawAxis();
}



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
#define momsNONirr 8
#define momsPirr2 13
#define irradiations 2



void TDR();
void TDR2(TCanvas * c_all, int period=0, int pos= 11);

void resolutionMomScanTreeEntries_forPaper(TString name = "preliminary_TreeCorrEntries", TString func = "RMSself", bool unfolding = false)
{


  TString inputDir="/home/zoiirene/Output/TextFiles/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  int measurements[irradiations];
  measurements[0] = momsNONirr;
  measurements[1] = momsPirr2;

  TString filenames[irradiations];
  filenames[0] = "Mscan_163_RMSself_resTree.txt";
  filenames[1] = "Mscan_130i_RMSself_resTree_beamdiv_A32C42.txt"; 

  double mom_0[momsNONirr];
  double mom_1[momsPirr2];
  double momerr_0[momsNONirr];
  double momerr_1[momsPirr2];

  double res_0[momsNONirr];
  double res_1[momsPirr2];

  double reserr_0[momsNONirr];
  double reserr_1[momsPirr2];

  double mom;
  double res,reserror,momerr;
  int lines[irradiations] = {8,8};
  
  TString filename[irradiations];
  for(int i =0; i< irradiations; i++)
    {
      filename[i]      = inputDir+filenames[i];
      cout << filename[i] << endl;
      ifstream stream(filename[i]);
      
      std::string line;
      if(!stream.is_open())
	{
	  cout << " File " << filename[i] << " not opened" << endl;
	}
      else
	{
	  for(int k =0; k<lines[i];k++)
	    {
	      std::getline(stream,line);
	      if(print) cout << k << " line " << line << endl;
	    }
	  //      while(!stream.eof())
	  for(int j= 0; j < measurements[i]; j++)
	    {
	      stream  >> mom >> res >> reserror ;
	      if(print)	 cout << "line " << j+lines[i]-1 << endl;
	      if(i==0)
		{
		  mom_0[j] = mom;
		  momerr_0[j] = 0.158;//GeV https://www.sciencedirect.com/science/article/pii/S0168900218317868?via%3Dihub#sec7.3
		  res_0[j] = res;
		  reserr_0[j] = reserror;
		  if(print)	 cout  << mom_0[j] << " " << res_0[j] << " " << reserr_0[j] << endl;
		}
	      if(i==1)
		{
		  mom_1[j] = mom;
		  momerr_1[j] = 0.158;
		  res_1[j] = res;
		  reserr_1[j] = reserror;
		  if(print)	 cout  << mom_1[j] << " " << res_1[j] << " " << reserr_1[j] << endl;
		}

	    }//measurements
	}//file scanning
    }//irradiations

  if( print ) cout << " got all measurements from file" << endl;


  /////// plots!

  
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
  

  TGraphErrors* resolutionPlot[irradiations];
  
  resolutionPlot[0] = new TGraphErrors(momsNONirr,mom_0,res_0,momerr_0,reserr_0);
  resolutionPlot[1] = new TGraphErrors(momsPirr2,mom_1,res_1,momerr_1,reserr_1);

  resolutionPlot[0]->SetTitle(" ");
  resolutionPlot[0]->GetYaxis()->SetTitle("#sigma_{x} [#mum]");
  resolutionPlot[0]->GetXaxis()->SetTitle("p_{beam} [GeV]");

  resolutionPlot[0]->GetYaxis()->SetTitleFont(43);
  resolutionPlot[0]->GetYaxis()->SetTitleSize(20); // labels will be 14 pixels
  resolutionPlot[0]->GetYaxis()->SetTitleOffset(1.2); // labels will be 14 pixels

  resolutionPlot[0]->GetYaxis()->SetLabelFont(43);
  resolutionPlot[0]->GetYaxis()->SetLabelSize(20); // labels will be 14 pixels

  resolutionPlot[0]->GetXaxis()->SetLabelFont(43);
  resolutionPlot[0]->GetXaxis()->SetLabelSize(20); // labels will be 14 pixels
  resolutionPlot[0]->GetXaxis()->SetTitleOffset(1.2); // labels will be 14 pixels
  
  resolutionPlot[0]->GetXaxis()->SetTitleFont(43);
  resolutionPlot[0]->GetXaxis()->SetTitleSize(20); // labels will be 14 pixels
  
  resolutionPlot[0]->SetMarkerSize(1.);
  resolutionPlot[0]->SetMarkerColor(kBlack);
  resolutionPlot[0]->SetLineColor(kBlack);
  resolutionPlot[0]->SetMarkerStyle(20);
  resolutionPlot[0]->GetXaxis()->SetLimits(0.,7.);
  resolutionPlot[0]->GetYaxis()->SetRangeUser(0.,25.);
  resolutionPlot[0]->Draw("AEP");
  cout << "fit non irr " << endl;
  TF1 *fnotsq = new TF1("fnotsq","TMath::Sqrt([0]+[1]/(x*x))", 0.2, 7.);
  fnotsq->SetLineColor(kGray);
  fnotsq->SetParName(0,"#sigma_{extr}^{2}");
  fnotsq->SetParName(1,"#sigma_{MS}^{2}");
  fnotsq->SetParLimits(0,0,20);
  fnotsq->SetParLimits(1,0,200);
  resolutionPlot[0]->Fit("fnotsq","R");

  TF1 *fit = resolutionPlot[0]->GetFunction("fnotsq");
  Double_t chi2 = fit->GetChisquare();
  Double_t ndof = fit->GetNDF();
  cout << " chi2 "<< chi2 << " / " << " ndof "<< ndof << " = "<< chi2/ndof << endl;

  
  fnotsq->Draw("same");
  ostringstream strfit4[2];
  TString ss_fit4[2];
  ostringstream strfit_err4[2];
  TString ss_fit_err4[2];


  for(int i=0; i<2;i++)
    {
      strfit4[i] << fixed <<setprecision(1) << sqrt(fnotsq->GetParameter(i));
      ss_fit4[i]=strfit4[i].str();
      strfit_err4[i] << setprecision(1) << 0.5*fnotsq->GetParError(i)/sqrt(fnotsq->GetParameter(i));
      ss_fit_err4[i]=strfit_err4[i].str();
    }

  
  resolutionPlot[1]->GetYaxis()->SetRangeUser(0.,20.);
  resolutionPlot[1]->SetMarkerSize(1.);
  resolutionPlot[1]->SetMarkerColor(kGreen+1);
  resolutionPlot[1]->SetLineColor(kGreen+1);
  resolutionPlot[1]->SetMarkerStyle(21);
  resolutionPlot[1]->Draw("EPsame");
  cout << "fit prot irr " << endl;
  TF1 *fnotsqi = new TF1("fnotsqi","TMath::Sqrt([0]+[1]/(x*x))", 0., 7.);
  fnotsqi->SetLineColor(kOrange);
  fnotsqi->SetParName(0,"#sigma_{extr}^{2}");
  fnotsqi->SetParName(1,"#sigma_{MS}^{2}");
  fnotsqi->SetParLimits(0,0,30);
  fnotsqi->SetParLimits(1,100,1000);
  //fnotsqi->SetParameter(0,fit32->GetParameter(0));
  //fnotsqi->SetParameter(1,fit32->GetParameter(1));
  resolutionPlot[1]->Fit("fnotsqi","R");
  resolutionPlot[1]->Fit("fnotsqi","R");
  resolutionPlot[1]->Fit("fnotsqi","R");
  resolutionPlot[1]->Fit("fnotsqi","R");

  fnotsqi->Draw("same");
  ostringstream strfit4i[2];
  TString ss_fit4i[2];
  ostringstream strfit_err4i[2];
  TString ss_fit_err4i[2];

  TF1 *fiti = resolutionPlot[1]->GetFunction("fnotsqi");
  Double_t chi2i = fiti->GetChisquare();
  Double_t ndofi = fiti->GetNDF();
  cout << " chi2 "<< chi2i << " / " << " ndof "<< ndofi << " = "<< chi2i/ndofi << endl;

  

  for(int i=0; i<2;i++)
    {
      strfit4i[i] << fixed <<setprecision(1) << sqrt(fnotsqi->GetParameter(i));
      ss_fit4i[i]=strfit4i[i].str();
      strfit_err4i[i] << setprecision(1) << 0.5*fnotsqi->GetParError(i)/sqrt(fnotsqi->GetParameter(i));
      ss_fit_err4i[i]=strfit_err4i[i].str();
    }


  TLegend* legFDB2 = new TLegend(0.4,0.48,0.85,0.87);
  legFDB2->SetLineColor(0);
  //legFDB2->SetTextSize(0.035);
  legFDB2->AddEntry(resolutionPlot[0],"Optimal incidence angle","");
  legFDB2->AddEntry(resolutionPlot[0],"Non-irradiated, 120 V","lp");
  legFDB2->AddEntry(fnotsq,"#sqrt{#sigma_{extr}^{2}+(#sigma_{MS}/p_{beam})^{2}}","l");
  legFDB2->AddEntry(fnotsq,"#sigma_{extr} = ("+ss_fit4[0]+" #pm "+ss_fit_err4[0]+") #mum","");
  legFDB2->AddEntry(resolutionPlot[1],"#phi_{eq} = 2.1 #times 10^{15} cm^{-2}, proton, 600 V","lp");
  //  legFDB2->AddEntry(resolutionPlot[1],"600V, optimal angle","");
  legFDB2->AddEntry(fnotsqi,"#sqrt{#sigma_{extr}^{2}+(#sigma_{MS}/p_{beam})^{2}}","l");
  legFDB2->AddEntry(fnotsqi,"#sigma_{extr} = ("+ss_fit4i[0]+" #pm "+ss_fit_err4i[0]+") #mum","");

  legFDB2->Draw();

  
  TString  outname = outputDir+"ResolutionSummaryPaper_momScans25_"+name;
  cFDB2->SaveAs(outname+".eps");
  cFDB2->SaveAs(outname+".png");
  cFDB2->SaveAs(outname+".pdf");
  cFDB2->SaveAs(outname+".root");
  cFDB2->SaveAs(outname+".C");

  // -----------------------------------------------------------------------
  double mom1_0[momsNONirr];
  double mom1_1[momsPirr2];
  double mom2_0[momsNONirr];
  double mom2_1[momsPirr2];
  double momerr2_0[momsNONirr];
  double momerr2_1[momsPirr2];

  double res2_0[momsNONirr];
  double res2_1[momsPirr2];

  double reserr2_0[momsNONirr];
  double reserr2_1[momsPirr2];
  
  cout << " non irr " << endl;
  for(int i=0; i<momsNONirr; i++){
    if(print) cout << "Energy " << mom_0[i] << "+- " << momerr_0[i] << " GeV -> Resolution: " << res_0[i] << " and res err: " << reserr_0[i] << endl;
      mom1_0[i]=1./mom_0[i];
      mom2_0[i]=mom1_0[i]*mom1_0[i];
      momerr2_0[i]=momerr_0[i]*2*TMath::Power(mom1_0[i],3);
      res2_0[i]=res_0[i]*res_0[i];
      reserr2_0[i]=2*reserr_0[i]*res_0[i];
      if(print) cout << "Energy^2 " << mom2_0[i] << "+- " << momerr2_0[i] << " GeV^2 -> Resolution2: " << res2_0[i] << " and res err: " << reserr2_0[i] << endl;
  }

  cout << " prot irr " << endl;
  for(int i=0; i<momsPirr2; i++){
      mom1_1[i]=1./mom_1[i];
      mom2_1[i]=mom1_1[i]*mom1_1[i];
      momerr2_1[i]=momerr_1[i]*2*TMath::Power(mom1_1[i],3);
      res2_1[i]=res_1[i]*res_1[i];
      reserr2_1[i]=2*reserr_1[i]*res_1[i];
      if(print) cout << "Energy^2 " << mom_1[i] << ": " << mom2_1[i] << " GeV^2 -> Resolution: " << res2_1[i] << " and res err: " << reserr2_1[i] << endl;
  }

  

  TCanvas *c42 = new TCanvas("c42", "resolution vs inverse beam energy squared", 600, 600);
  c42->SetLeftMargin(0.12);
  c42->SetRightMargin(0.05);
  c42->SetTopMargin(0.05);
  c42->SetBottomMargin(0.1);

  gPad->SetTicks(1,1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(15);


  TGraphErrors* resolutionPlotInvSquare[irradiations];

  resolutionPlotInvSquare[0] = new TGraphErrors(momsNONirr,mom2_0,res2_0,momerr2_0,reserr2_0);
  resolutionPlotInvSquare[1] = new TGraphErrors(momsPirr2,mom2_1,res2_1,momerr2_1,reserr2_1);

  resolutionPlotInvSquare[0]->SetTitle(" ");
  resolutionPlotInvSquare[0]->GetYaxis()->SetTitle("#sigma_{x}^{2} [#mum^{2}]");
  resolutionPlotInvSquare[0]->GetXaxis()->SetTitle("p_{beam}^{2} [GeV^{2}]");

  resolutionPlotInvSquare[0]->GetYaxis()->SetTitleFont(43);
  resolutionPlotInvSquare[0]->GetYaxis()->SetTitleSize(20); // labels will be 14 pixels                                                                                                                                  
  resolutionPlotInvSquare[0]->GetYaxis()->SetTitleOffset(1.4); // labels will be 14 pixels                                                                                                                               

  resolutionPlotInvSquare[0]->GetYaxis()->SetLabelFont(43);
  resolutionPlotInvSquare[0]->GetYaxis()->SetLabelSize(20); // labels will be 14 pixels                                                                                                                                  

  resolutionPlotInvSquare[0]->GetXaxis()->SetLabelFont(43);
  resolutionPlotInvSquare[0]->GetXaxis()->SetLabelSize(20); // labels will be 14 pixels                                                                                                                                  
  resolutionPlotInvSquare[0]->GetXaxis()->SetTitleOffset(1.2); // labels will be 14 pixels                                                                                                                               

  resolutionPlotInvSquare[0]->GetXaxis()->SetTitleFont(43);
  resolutionPlotInvSquare[0]->GetXaxis()->SetTitleSize(20); // labels will be 14 pixels                                                                                                                                  

  resolutionPlotInvSquare[0]->SetMarkerSize(1.);
  resolutionPlotInvSquare[0]->SetMarkerColor(kBlack);
  resolutionPlotInvSquare[0]->SetLineColor(kBlack);
  resolutionPlotInvSquare[0]->SetMarkerStyle(20);
  resolutionPlotInvSquare[0]->GetXaxis()->SetLimits(0.,7.);
  resolutionPlotInvSquare[0]->GetYaxis()->SetRangeUser(0.,25.);
  resolutionPlotInvSquare[0]->Draw("AEP");
  resolutionPlotInvSquare[1]->GetYaxis()->SetRangeUser(0.,20.);
  resolutionPlotInvSquare[1]->SetMarkerSize(1.);
  resolutionPlotInvSquare[1]->SetMarkerColor(kGreen+1);
  resolutionPlotInvSquare[1]->SetLineColor(kGreen+1);
  resolutionPlotInvSquare[1]->SetMarkerStyle(21);
  resolutionPlotInvSquare[1]->Draw("EPsame");


  
  resolutionPlotInvSquare[0]->GetYaxis()->SetRangeUser(0.,300.);
  resolutionPlotInvSquare[0]->GetXaxis()->SetLimits(0.,1.);

  TF1 *fit32 = new TF1("fit32","pol1", 0., 0.45);
  fit32->SetLineColor(kGray);
  fit32->SetParName(0,"#sigma_{intr}^{2}");
  fit32->SetParName(1,"#sigma_{MS}^{2}");

  resolutionPlotInvSquare[0]->Fit("fit32","R");
  ostringstream strfit2[2];
  TString ss_fit2[2];
  ostringstream strfit_err2[2];
  TString ss_fit_err2[2];


  TF1 *fitI = new TF1("fitI","pol1", 0., 0.7);
  fitI->SetLineColor(kOrange);
  fitI->SetParName(0,"#sigma_{intr}^{2}");
  fitI->SetParName(1,"#sigma_{MS}^{2}");

  resolutionPlotInvSquare[1]->Fit("fitI","R");
  ostringstream strfitI[2];
  TString ss_fitI[2];
  ostringstream strfit_errI[2];
  TString ss_fit_errI[2];


  for(int i=0; i<2;i++)
    {
      strfit2[i] << fixed <<setprecision(1) << sqrt(fit32->GetParameter(i));
      ss_fit2[i]=strfit2[i].str();
      strfit_err2[i] << setprecision(1) << 0.5*fit32->GetParError(i)/sqrt(fit32->GetParameter(i));
      ss_fit_err2[i]=strfit_err2[i].str();

      strfitI[i] << fixed <<setprecision(1) << sqrt(fitI->GetParameter(i));
      ss_fitI[i]=strfitI[i].str();
      strfit_errI[i] << setprecision(1) << 0.5*fitI->GetParError(i)/sqrt(fitI->GetParameter(i));
      ss_fit_errI[i]=strfit_errI[i].str();

    }



  TLegend* lgI2 = new TLegend(0.45,0.3,0.88,0.6);
  lgI2->SetLineWidth(0);
  lgI2->SetFillStyle(0);
  lgI2->AddEntry(resolutionPlotInvSquare[0],"Optimal incidence angle","");
  lgI2->AddEntry(resolutionPlotInvSquare[0],"Non-irradiated, 120 V","lp");
  lgI2->AddEntry(fit32,"#sigma_{intr}^{2}+(#sigma_{MS}/p_{beam})^{2}","l");
  lgI2->AddEntry(fnotsq,"#sigma_{extr} = ("+ss_fit2[0]+" #pm "+ss_fit_err2[0]+") #mum","");
  lgI2->AddEntry(resolutionPlotInvSquare[1],"#phi_{eq} = 2.1 #times 10^{15} cm^{-2}, proton, 600 V","lp");
  //  lgI2->AddEntry(resolutionPlot[1],"600V, optimal angle","");                                                                                                                                            
  lgI2->AddEntry(fitI,"#sigma_{intr}^{2}+(#sigma_{MS}/p_{beam})^{2}","l");
  lgI2->AddEntry(fnotsqi,"#sigma_{extr} = ("+ss_fitI[0]+" #pm "+ss_fit_errI[0]+") #mum","");

  lgI2->Draw();

  /*
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
                               
  leg42->AddEntry(resolutionPlotInvSquare[0],"Non-irradiated, 120V, optimal angle" ,"ep");
  leg42->AddEntry(fit32,"#sigma_{intr}^{2}+(#sigma_{MS}/momentum)^{2}" , "l");
  leg42->Draw();
  */
  c42->Update();
  //TDR2(c42);                                                                                                                                                                                                                                
  name = outputDir+"ResolutionSummaryPaper_momScans25Squared_"+name;
  c42->SaveAs(name+".eps");
  c42->SaveAs(name+".pdf");
  c42->SaveAs(name+".png");
  c42->SaveAs(name+".root");





  
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

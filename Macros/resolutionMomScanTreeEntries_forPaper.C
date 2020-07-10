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
  filenames[1] = "Mscan_130i_RMSself_resTree_closest_A32C42.txt"; 

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
		  momerr_0[j] = 0.35;//GeV
		  res_0[j] = res;
		  reserr_0[j] = reserror;
		  if(print)	 cout  << mom_0[j] << " " << res_0[j] << " " << reserr_0[j] << endl;
		}
	      if(i==1)
		{
		  mom_1[j] = mom;
		  momerr_1[j] = 0.35;
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
  fnotsq->SetParName(0,"#sigma_{hit}^{2}");
  fnotsq->SetParName(1,"#sigma_{MS}^{2}");
  fnotsq->SetParLimits(0,0,20);
  fnotsq->SetParLimits(1,0,200);
  resolutionPlot[0]->Fit("fnotsq","R");


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
  fnotsqi->SetParName(0,"#sigma_{hit}^{2}");
  fnotsqi->SetParName(1,"#sigma_{MS}^{2}");
  fnotsqi->SetParLimits(0,0,50);
  fnotsqi->SetParLimits(1,0,1000);
  //fnotsqi->SetParameter(0,fit32->GetParameter(0));
  //fnotsqi->SetParameter(1,fit32->GetParameter(1));
  resolutionPlot[1]->Fit("fnotsqi","R");

  fnotsqi->Draw("same");
  ostringstream strfit4i[2];
  TString ss_fit4i[2];
  ostringstream strfit_err4i[2];
  TString ss_fit_err4i[2];


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
  legFDB2->AddEntry(resolutionPlot[0],"Optimal incident angle","");
  legFDB2->AddEntry(resolutionPlot[0],"Non-irradiated, 120 V","lp");
  legFDB2->AddEntry(fnotsq,"#sqrt{#sigma_{intr}^{2}+(#sigma_{MS}/p_{beam})^{2}}","l");
  legFDB2->AddEntry(fnotsq,"#sigma_{intr} = ("+ss_fit4[0]+" #pm "+ss_fit_err4[0]+") #mum","");
  legFDB2->AddEntry(resolutionPlot[1],"#phi_{eq} = 2.1 #times 10^{15} cm^{-2}, proton, 600 V","lp");
  //  legFDB2->AddEntry(resolutionPlot[1],"600V, optimal angle","");
  legFDB2->AddEntry(fnotsqi,"#sqrt{#sigma_{intr}^{2}+(#sigma_{MS}/p_{beam})^{2}}","l");
  legFDB2->AddEntry(fnotsqi,"#sigma_{intr} = ("+ss_fit4i[0]+" #pm "+ss_fit_err4i[0]+") #mum","");

  legFDB2->Draw();

  
  TString  outname = outputDir+"ResolutionSummaryPaper_momScans25_"+name;
  cFDB2->SaveAs(outname+".eps");
  cFDB2->SaveAs(outname+".png");
  cFDB2->SaveAs(outname+".pdf");
  cFDB2->SaveAs(outname+".root");
  cFDB2->SaveAs(outname+".C");







  
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

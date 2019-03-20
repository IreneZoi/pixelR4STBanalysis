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

#define BeamEnergies 9
#define planes 3
bool print=true;
using namespace std;

void FitTH1(TH1F* h1, Double_t *  sigma, Double_t *  sigmaerr, TString name, TString detectorA, TString detectorB, TString detectorC, TString func );
Double_t tp0Fit( Double_t *x, Double_t *par );
double G = 1E0;

void prepareRootFile()
{

  TString detectorA="158";
  TString detectorB="152";
  TString detectorC="159";
  TString labelA="FDD150P_22_R4S50x50-P9_2, thr 12 ADC";
  TString labelB="FDD150P_22_R4S50x50-P1_1, 120 V, thr 12 ADC";
  TString labelC="FDB150P_12_R4S50x50-P3_1, thr 12 ADC";
  TString info = "angle 17.5 deg";
  
  
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
  TFile * file[BeamEnergies];
  TH1F * h_res[BeamEnergies];
  TH2F * h_eff[BeamEnergies];
  TH1F * h_clsz[planes][BeamEnergies];
  TH1F * h_landau[planes][BeamEnergies];

  Int_t run=1027;
  TString Run, Path;
  ostringstream strs[BeamEnergies];
  
  if(print) cout << "Getting files "  << endl;
  for(int i=0; i<BeamEnergies; i++)
    {
      if(print) cout << "Run: " << run << endl;
      Run.Form("%d",run);

      if(print) cout << "Run: " << Run << endl;
      inputfile = "drei-r"+Run+"_irene.root";
      if(print) cout << "File Name: " << inputfile << endl;
      Path=inputDir+inputfile;
      file[i] = new TFile(Path);
      if(print) cout << "File Path: " << Path << endl;
      run +=1;
      
    }

  if(print) cout << "Iniazializating beam energies "  << endl;
  for(int i=0; i<BeamEnergies; i++)
    { 
      BeamEnergy[i]=1.2+0.4*i;//GeV
      if(i>3) BeamEnergy[i]=2.4+0.8*(i-3);//GeV
      if(i==8) BeamEnergy[i]=6.0;//GeV      
      strs[i] << BeamEnergy[i];
      ss_BeamEnergy[i]=strs[i].str();
      if(print) cout << "Beam energy " << i<< ": " << ss_BeamEnergy[i] << " GeV" << endl;
      
    }

  if(print) cout << "Ordering histos by beam energy and not run number:"  << endl;
  TString hist = "dx3ciiiqr3";
  h_res[0]=(TH1F*)file[7]->Get(hist); 
  h_res[1]=(TH1F*)file[6]->Get(hist);
  h_res[2]=(TH1F*)file[5]->Get(hist);
  h_res[3]=(TH1F*)file[4]->Get(hist);
  h_res[4]=(TH1F*)file[3]->Get(hist);
  h_res[5]=(TH1F*)file[2]->Get(hist);
  h_res[6]=(TH1F*)file[1]->Get(hist);
  h_res[7]=(TH1F*)file[0]->Get(hist);
  h_res[8]=(TH1F*)file[8]->Get(hist);
  hist = "effvsxy";
  h_eff[0]=(TH2F*)file[7]->Get(hist); 
  h_eff[1]=(TH2F*)file[6]->Get(hist);
  h_eff[2]=(TH2F*)file[5]->Get(hist);
  h_eff[3]=(TH2F*)file[4]->Get(hist);
  h_eff[4]=(TH2F*)file[3]->Get(hist);
  h_eff[5]=(TH2F*)file[2]->Get(hist);
  h_eff[6]=(TH2F*)file[1]->Get(hist);
  h_eff[7]=(TH2F*)file[0]->Get(hist);
  h_eff[8]=(TH2F*)file[8]->Get(hist);
  for(int i =0; i< planes; i++)
    {
      if(i==0)hist="clszAiii";
      if(i==1)hist="clszBiii";
      if(i==2)hist="clszCiii";
      h_clsz[i][0]=(TH1F*)file[7]->Get(hist); 
      h_clsz[i][1]=(TH1F*)file[6]->Get(hist);
      h_clsz[i][2]=(TH1F*)file[5]->Get(hist);
      h_clsz[i][3]=(TH1F*)file[4]->Get(hist);
      h_clsz[i][4]=(TH1F*)file[3]->Get(hist);
      h_clsz[i][5]=(TH1F*)file[2]->Get(hist);
      h_clsz[i][6]=(TH1F*)file[1]->Get(hist);
      h_clsz[i][7]=(TH1F*)file[0]->Get(hist);
      h_clsz[i][8]=(TH1F*)file[8]->Get(hist);
  
    }
for(int i =0; i< planes; i++)
    {
      if(i==0)hist="clqAiii";
      if(i==1)hist="clqBiii";
      if(i==2)hist="clqCiii";
      h_landau[i][0]=(TH1F*)file[7]->Get(hist); 
      h_landau[i][1]=(TH1F*)file[6]->Get(hist);
      h_landau[i][2]=(TH1F*)file[5]->Get(hist);
      h_landau[i][3]=(TH1F*)file[4]->Get(hist);
      h_landau[i][4]=(TH1F*)file[3]->Get(hist);
      h_landau[i][5]=(TH1F*)file[2]->Get(hist);
      h_landau[i][6]=(TH1F*)file[1]->Get(hist);
      h_landau[i][7]=(TH1F*)file[0]->Get(hist);
      h_landau[i][8]=(TH1F*)file[8]->Get(hist);
  
    }


 for(int i=0; i< BeamEnergies; i++)
   {
     TFile *MyFile = new TFile("hist_"+ss_BeamEnergy[i]+".root","RECREATE");
     h_eff[i]->Write("efficiency");
     h_clsz[0][i]->Write("cluster_size_A");
     h_clsz[1][i]->Write("cluster_size_B");
     h_clsz[2][i]->Write("cluster_size_C");
     h_landau[0][i]->Write("cluster_charge_A");
     h_landau[1][i]->Write("cluster_charge_B");
     h_landau[2][i]->Write("cluster_charge_C");
     h_res[i]->Write("residuals");
     
     MyFile->Close();
   }
  
}//resolution 


void FitTH1(TH1F* h1, Double_t *  sigma, Double_t *  sigmaerr, TString name, TString detectorA, TString detectorB, TString detectorC, TString func)
{
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  //  gPad->SetTickx();
  //gPad->SetTicky();
  
  TCanvas * c = new TCanvas("c","c",700,700);

  c->cd();
  c->SetFrameFillStyle(1000);
  c->SetFrameFillColor(0);
  gPad->SetTicks(1,1);


  h1->SetTitle("");
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetLabelSize(0.025);
  h1->GetXaxis()->SetTitleSize(0.035);
  h1->GetXaxis()->SetTitleOffset(0.8);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetRangeUser(-0.1,0.1);

  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelSize(0.025);
  h1->GetYaxis()->SetTitleSize(0.035);
  h1->GetYaxis()->SetTitleOffset(1.4);
  h1->GetYaxis()->SetTitleFont(42);

  h1->GetYaxis()->SetTitle("Entries");
  h1->GetXaxis()->SetTitle("residual [mm]");

  //TF1 * MyGaus = new TF1("MyGaus","gaus", -0.06,0.06);
  h1->SetMarkerStyle(20);
  h1->Draw("PZ");
  //h1->Fit("MyGaus","R");
if(func == "gaus")
    {
      TF1 * MyGaus = new TF1("MyGaus","gaus", -0.06,0.06);
      h1->Fit("MyGaus","RQ");
    
      double scale=1.5;

      double lower_bound;
      double upper_bound;

      for(int k=0;k<4;k++)
	{
	  lower_bound = MyGaus->GetParameter(1)-scale*MyGaus->GetParameter(2);
	  upper_bound = MyGaus->GetParameter(1)+scale*MyGaus->GetParameter(2);
	  MyGaus = new TF1("MyGaus","gaus", lower_bound,upper_bound);
	  
	  h1->Fit("MyGaus","R");
	}


      *sigma = MyGaus->GetParameter(2)*1000;
      *sigmaerr = MyGaus->GetParError(2)*1000;
      
      MyGaus->SetLineColor(kRed);  
      MyGaus->SetLineWidth(2);  
      MyGaus->Draw("same");
    }
  if(func == "studentT")
    {
      const int mpar = 5;
      double x1=1;
      double x9=0;
      double dx = h1->GetBinWidth(1);
      double nmax = h1->GetBinContent(h1->GetMaximumBin());
      double xmax = h1->GetBinCenter(h1->GetMaximumBin());
      int nb = h1->GetNbinsX();
      if( x9 < x1 ) {
	x1 = h1->GetBinCenter(1);
	x9 = h1->GetBinCenter(nb);
      }
      int i1 = h1->FindBin(x1);
      int i9 = h1->FindBin(x9);
      double n1 = h1->GetBinContent(i1);
      double n9 = h1->GetBinContent(i9);
      double bg = 0.5*(n1+n9);

      double nn = 7*(nmax-bg);
      double G = 1E0;
      
      TF1 *tp0Fcn = new TF1( "tp0Fcn", tp0Fit, x1 , x9, mpar );

      tp0Fcn->SetParName( 0, "mean" );
      tp0Fcn->SetParName( 1, "sigma" );
      tp0Fcn->SetParName( 2, "nu" );
      tp0Fcn->SetParName( 3, "area" );
      tp0Fcn->SetParName( 4, "BG" );

      // set start values for some parameters:

      cout << "start dx " << dx
	   << ", max " << nmax
	   << " at " << xmax
	   << ", bg " << bg
	   << ", area " << nn
	   << endl;

      tp0Fcn->SetParameter( 0, xmax ); // peak position
      tp0Fcn->SetParameter( 1, 4*dx ); // width
      tp0Fcn->SetParameter( 2, 2.2 ); // nu
      tp0Fcn->SetParameter( 3, nn/G ); // N
      tp0Fcn->SetParameter( 4, bg );

      tp0Fcn->SetNpx(500);

      cout << endl << "Minos:" << endl << endl;
      h1->Fit( "tp0Fcn", "RME", "ep" );
      tp0Fcn->SetLineColor(kRed);  
      tp0Fcn->SetLineWidth(2);  
      tp0Fcn->Draw("same");

      *sigma = tp0Fcn->GetParameter(1)*1000.;
      *sigmaerr = tp0Fcn->GetParError(1)*1000.;

    }


  

  c->Update();
  gStyle->SetOptFit(1111);
  TString outputDir="/home/zoiirene/Output/Plots/";

  //  TString outputDir = outputDir+"";
  TString outputFile = outputDir+"residual_"+func+"Fit_"+name+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps";
  c->SaveAs(outputFile);
}//



Double_t tp0Fit( Double_t *x, Double_t *par )
{
  static int nn = 0;
  nn++;
  static double dx = 0.1;
  static double b1 = 0;
  if( nn == 1 ) b1 = x[0];
  if( nn == 2 ) dx = x[0] - b1;

  // Mean and width:

  double xm = par[0];
  double t = ( x[0] - xm ) / par[1];
  double tt = t*t;

  // exponent:

  double rn = par[2];
  double xn = 0.5 * ( rn + 1.0 );

  // Normalization needs Gamma function:

  double pk = 0.0;

  if( rn > 0.0 && fabs( xn * log( 1.0 + tt/rn ) ) < 333 ) {

    double pi = 3.14159265358979323846;
    double aa = dx / par[1] / sqrt(rn*pi) * TMath::Gamma(xn) / TMath::Gamma(0.5*rn);

    pk = G * par[3] * aa * exp( -xn * log( 1.0 + tt/rn ) );

    // lim n->inf (1+a/n)^n = e^a

  }

  return pk + par[4];
}

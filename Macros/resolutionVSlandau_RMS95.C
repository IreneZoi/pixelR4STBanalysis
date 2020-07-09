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

#include "fileHandler.h"
#define Cuts 5
#define conversionFactor 1.
#define comparisons 1
#define DreiMasterPlanes 3
bool print=true;
using namespace std;


void resolutionVSlandau_RMS95(TString function = "RMS95"){

  TString detectorA="146";
  TString detectorB="148";
  TString detectorC="163";
  TString labelA="FDB150Y_2_R4S100x25-Y2_2, thr 22 ADC";
  TString labelC="FDB150Y_2_R4S100x25-Y6_1, thr 22 ADC";
  TString labelB="Pstop_RD53Apads_FDB, thr 22 ADC";
  
  Int_t run=2735;        //{2735,2743,2758};
  TString Run;
  Run.Form("%d",run);
  if(print) cout << "Run: " << run << endl;
  
  //  TString Angle[comparisons] = {"0","8.75","27.5"};
  TString info = "beam energy 5.6 GeV, angle 0 deg";
  
  TString Landau[Cuts];
  Double_t Resolution[Cuts];
  Double_t ResolutionError[Cuts];
  
  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";
  TH1F * h_res[Cuts];
  TString label = "testABC";
  TString ss_dphcut = "22";


  TString  inputfile = inputDir+"drei-r"+Run+"_irene_dphcutB"+ss_dphcut+"_"+label+".root";
  cout << inputfile << endl;
  TFile * file = new TFile(inputfile);
  TTree * tree = (TTree*)file->Get("charge_res");
  tree->Print();
    
  TH1I * hdx3tree = new TH1I("hdx3tree", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
  tree->Draw("dx3tree>>hdx3tree","","goff");
  hdx3tree = (TH1I*)gDirectory->Get("hdx3tree");
  
  double integral[DreiMasterPlanes];
  double integral_per[DreiMasterPlanes][Cuts];
  
  int high[DreiMasterPlanes][Cuts];
  
  
  TString dir = "straightTracksY_isoAandCandB_straightTracksX";
  TH1I * hclph[DreiMasterPlanes];
  TDirectory * histodir =(TDirectoryFile*)file->Get(dir);
  hclph[0]= (TH1I*)histodir->Get("clphA");
  cout << "hclphAiii "<< hclph[0]->GetEntries() << endl;
  hclph[1]= (TH1I*)histodir->Get("clphB");
  hclph[2]= (TH1I*)histodir->Get("clphC");
  
  for(int j =0; j < DreiMasterPlanes; j++)   {
    if(print)   cout << " plane " << j << endl;
    integral[j] = hclph[j]->Integral(0,hclph[j]->GetNbinsX()+1);
    if(print)      cout << " integral " << integral[j] << " entries " << hclph[j]->GetEntries()<<  endl;
    
    for(int k = 0; k< Cuts; k++){
      float perc = (float) (k+1) * 0.2 ;
      cout << " perc " << perc << endl;
      integral_per[j][k] = perc *integral[j]; //careful!! you should take into account the peak at low value! Hist is now filled only for isolated clusters and the effect is reduced
      cout << " integral " << integral_per[j][k] << endl;
      
      high[j][k] = 0;
	
      double integrating = 0;
      int i = 1;
      while(integrating<integral_per[j][k])    {
	
	cout << " while "<< i << endl;
	  integrating= hclph[j]->Integral(1,i); //hclph[j]->GetNbinsX()-i);
	  high[j][k] = hclph[j]->GetBinCenter(i); //hclph[j]->GetNbinsX()-i);
	  cout << " integral " << integrating << " high " << high[j][k] << endl;
	  i++;
      }
      cout << " integral "<< perc << " " << integral_per[j][k] << " obtained " << integrating << endl;
    }
  }//drei planes

  
  if(print){

    for(int j =0; j < DreiMasterPlanes; j++)   {
      cout << " plane " << j << " total " <<  integral[j] << endl ;
      for(int k = 0; k< Cuts; k++){
	cout << " perc " << (float) (k+1) * 0.1 << " integral " << integral_per[j][k] << " high bin " << high[j][k] << " integral high " << hclph[j]->Integral(1,hclph[j]->FindBin(high[j][k])) << endl;
      }
    }
  }

    
 

  
  
  TString ph[DreiMasterPlanes][Cuts];
  for(int j =0; j < DreiMasterPlanes; j++)   {
    for(int k = 0; k< Cuts; k++){    
      ph[j][k].Form("%d",high[j][k]);
    }
  }

  TString planes[DreiMasterPlanes] = {"A","B","C"};
  TFile *MyFile = new TFile("dx3_landau_"+Run+"_dphcut"+ss_dphcut+"_"+label+".root","RECREATE"); 
  //plot landaus
  TCanvas *cFDB2[DreiMasterPlanes];
  for(int j =0; j < DreiMasterPlanes; j++)   {
    cFDB2[j] = new TCanvas("cFDB2", "FDB resolution", 1500, 900);
    gPad->SetTicks(1,1);
    gROOT->SetStyle("Plain");
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetPalette(1);
    //    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    hclph[j]->Draw();
    TLine *  linea[Cuts];
    for(int k = 0; k< Cuts; k++){
      linea[k]    = new TLine( high[j][k],0.,high[j][k],2000.);
      linea[k]->SetLineColor(kRed);
      linea[k]->SetLineWidth(2);
      linea[k]->SetLineStyle(2);
      linea[k]->Draw("same");
    }
    cFDB2[j]->Write("landau_"+planes[j]);
  }







  
  TH1I * hdx3ph[Cuts];
  
  TString totbins;
  int maximumedge = hclph[0]->GetBinLowEdge(hclph[0]->GetNbinsX()+1)+hclph[0]->GetBinWidth(1);
  totbins.Form("%d",maximumedge);
  cout <<" totbins "<< totbins<< endl;
  TString ss_perc[Cuts];

  for(int k = 0; k< Cuts; k++){
    ss_perc[k].Form("%d",k);
    TH1I * hdx3treeph = new TH1I("hdx3treeph", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
    hdx3ph[k] = new TH1I("hdx3treeph", "triplet dx3  ; dx [mm];triplets", 500, -0.5, 0.5 );

    TString phA_low,phA_high, phB_low,phB_high, phC_low,phC_high;
    if(k==0){
      phA_low="0";
      phB_low="0";
      phC_low="0";
      phA_high=ph[0][0];     
      phB_high=ph[1][0];     
      phC_high=ph[2][0];     

    }
    else{
      phA_low=ph[0][k-1];
      phB_low=ph[1][k-1];
      phC_low=ph[2][k-1];
      phA_high=ph[0][k];
      phB_high=ph[1][k];
      phC_high=ph[2][k];
      
    }


    cout << " cut " << k << endl;
    cout << " A " << phA_low << " - " << phA_high << endl;
    cout << " B " << phB_low << " - " << phB_high << endl;
    cout << " C " << phC_low << " - " << phC_high << endl;
    
    tree->Draw("dx3tree>>hdx3treeph","clphAiiitree>"+phA_low+"&&clphBiiitree>"+phB_low+"&&clphCiiitree>"+phC_low+"&&clphAiiitree<"+phA_high+"&&clphBiiitree<"+phB_high+"&&clphCiiitree<"+phC_high,"goff");
    //tree->Draw("dx3tree>>hdx3treeph","clphBiiitree>"+phB_low+"&&clphBiiitree<"+phB_high) ; 
    //hdx3treeph =
    hdx3ph[k] = (TH1I*)gDirectory->Get("hdx3treeph");
    
    hdx3ph[k]->SetDirectory(0);
    cout << "filled" << endl;
    hdx3treeph->Delete();

    cout << "deleted" << endl;
    cout << hdx3ph[k]->GetTitle() << " entries " << hdx3ph[k]->GetEntries() << endl;
    hdx3ph[k]->GetXaxis()->SetRangeUser(-0.05,0.05);
    hdx3ph[k]->Write("hdx3ph"+ss_perc[k]);
    hdx3ph[k]->GetXaxis()->SetRangeUser(-0.5,0.5);

  } // cuts

  MyFile->Close();

  
  // now fill the 95 % rms

  double sigma[Cuts] ;
  double sigmaerr[Cuts] ;
  

  for(int k = 0; k< Cuts; k++){

  double tolerance = (hdx3ph[k]->GetBinContent(0)+hdx3ph[k]->GetBinContent(hdx3ph[k]->GetNbinsX()+1))/hdx3ph[k]->GetEntries(); // since we are not working with continuous quantities but with binned hists, the difference between the integ ral and the 95% it may not be zero, so we ask it to be lower than 15% (before I was using 1.1% but since I changed to use the full range of hists, the overflow bin plays a too bigger role when stats is low and the 1.1% is not reached.

  cout << " RMS method charge based with tolerance " << tolerance << endl;
  double maximum = hdx3ph[k]->GetMaximum();
  int maxbin = hdx3ph[k]->GetMaximumBin();
  cout <<  " maxbin " << maxbin << " at " << hdx3ph[k]->GetBinCenter(maxbin)<<endl;
  cout << "intital sigma = " << hdx3ph[k]->GetRMS() * 1000 << " sigmaerr = " << hdx3ph[k]->GetRMSError() * 1000 << endl;
  double Integral = hdx3ph[k]->Integral(0,hdx3ph[k]->GetNbinsX()+1);
  cout << " integral " << Integral << " entries " << hdx3ph[k]->GetEntries() << endl;
  double integral95 = 0.95*Integral;
  cout << " integral95 " << integral95 << endl;
  int i = 0;
  
  double low = 0;
  double high = 0;
  for(int i =0; i<hdx3ph[k]->GetNbinsX()/2; i++)  {
    low = maxbin-i;
    high = maxbin+i;
    
    Integral = hdx3ph[k]->Integral(low,high);
    cout << " integral " << Integral << " low " << low << " high " << high << endl;
    cout << " while "<< i << " fabs(integral-integral95)/integral95 " << fabs(Integral-integral95)/integral95 << endl;
    
    //      if(fabs(Integral-integral95)/integral95 < 0.011 || Integral>integral95)//integral>integral95)
    if(fabs(Integral-integral95)/integral95 < tolerance || Integral>integral95)//integral>integral95)
      break;
    
    
  }
  cout << "final integral " << Integral << " low " << low <<" high " << high << endl;
  hdx3ph[k]->GetXaxis()->SetRange(low,high);


  
  
  sigma[k] = hdx3ph[k]->GetRMS() * 1000;
  sigmaerr[k] = hdx3ph[k]->GetRMSError() * 1000;

  cout << " RMS " << sigma[k] << " ± " << sigmaerr[k] << endl;
  FitTH1(hdx3ph[k], &(sigma[k]), &(sigmaerr[k]), ss_perc[k], detectorA, detectorB, detectorC, function );
  ExtractRes(&(sigma[k]), &(sigmaerr[k]));
  cout << " RES " << sigma[k] << " ± " << sigmaerr[k] << endl;

  }




  /*
  
  Double_t edges[Cuts + 2] = {0.0, (Double_t) high[1][0], (Double_t)high[1][1],(Double_t)high[1][2],(Double_t)high[1][3],(Double_t)high[1][4],1000.};
  // Bin 1 corresponds to range [0.0, 0.2]
  // Bin 2 corresponds to range [0.2, 0.3] etc...

  //TH1* h_res_vs_charge = new TH1D("h_res_vs_charge","Resolution vs cluster charge",Cuts,edges);


  for(int i=0; i<Cuts; i++)    {
      h_res_vs_charge->SetBinContent(i+1,sigma[i]);
      h_res_vs_charge->SetBinError(i+1,sigmaerr[i]);
    }
  */
  
  Double_t charge[Cuts] = {
     (Double_t) high[1][0]/2.,
     ((Double_t)high[1][1]+(Double_t) high[1][0])/2.,
     ((Double_t)high[1][1]+(Double_t)high[1][2])/2.,
     ((Double_t)high[1][3]+(Double_t)high[1][2])/2.,
     ((Double_t)high[1][4]+(Double_t)high[1][3])/2.};

  TGraph * h_res_vs_charge = new TGraph(Cuts,charge,sigma);  
  
 
  TCanvas *c2 = new TCanvas("c2", "resolution vs dphcut", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  h_res_vs_charge->SetTitle(" ");
  h_res_vs_charge->GetYaxis()->SetTitle("RMS 95% [#mum]");
  h_res_vs_charge->GetXaxis()->SetTitle("Cluster charge [ADC]");
  h_res_vs_charge->SetMarkerSize(2.5);
  h_res_vs_charge->SetLineColor(2);
  h_res_vs_charge->SetMarkerColor(2);
  h_res_vs_charge->SetMarkerStyle(20);
  h_res_vs_charge->SetLineWidth(2);
  h_res_vs_charge->SetLineStyle(2);
  h_res_vs_charge->GetXaxis()->SetRangeUser(0.,1000.);
  h_res_vs_charge->GetYaxis()->SetRangeUser(0.,15.);
  //h_res_vs_charge->Draw("E1X0");
  
  //h_res_vs_charge->Draw("P");
  //  h_res_vs_charge->Draw("HITS L");

  h_res_vs_charge->Draw("APL");
    
  TLegend* leg2 = new TLegend(0.15,0.6,0.5,0.8);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.02);
  leg2->AddEntry(h_res_vs_charge,"c"+detectorB+" "+labelB , "lp");
  leg2->AddEntry(h_res_vs_charge,"A: c"+detectorA+" "+labelA ,"");
  leg2->AddEntry(h_res_vs_charge,"C: c"+detectorC+" "+labelC ,"");
  leg2->AddEntry(h_res_vs_charge,info ,"");
  leg2->Draw();
  TLine *  linea[Cuts];
  for(int k = 0; k< Cuts; k++){
    linea[k]    = new TLine( high[1][k],0.,high[1][k],6.);
    linea[k]->SetLineColor(kGray);
    linea[k]->SetLineWidth(2);
    linea[k]->SetLineStyle(2);
    linea[k]->Draw("same");
  }
  

  TString name = outputDir+"Res_vs_clusterchergeABC_"+function+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+"_run"+Run;
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".pdf");
  
}//resolution 


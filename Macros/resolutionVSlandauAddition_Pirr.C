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
#define Cuts 20
#define conversionFactor 1.
#define Angles 28 //32 
#define DreiMasterPlanes 3
bool print=true;
using namespace std;


void resolutionVSlandauAddition_Pirr(TString function = "RMSself"){



  TString detectorA="148";
  TString detectorB="120i";
  TString detectorC="163";
  TString labelA="Pstop_RD53Apads_FDB, thr 12 ADC";
  TString labelB="Pstop_default_FTH, 800V, thr 15 ADC";
  TString labelC="Pspray_RD53Apads_FDB, thr 12 ADC";
  TString info = "beam energy 5.6 GeV";
  TString testname="closest_A12C15_bestproton";
  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  TString testnameUnf="RMSself6sig_closest_A13C14_bestnonirr";
  TString  filenames = "Ascan_resTree_148_dphcutB12_"+testnameUnf+".txt";
  filenames      = inputDir+filenames;
  cout << filenames << endl;
  ifstream stream(filenames);
  double angle;
  double res,reserror;
  double freshresA[Angles];
  double freshresA_err[Angles];
  double average[Cuts];
  std::string line;
  if(!stream.is_open())
    {
      cout << " File " << filenames << " not opened" << endl;
    }
  else
    {
      for(int k =0; k<8;k++)
	{
	  std::getline(stream,line);
	  if(print) cout << k << " line " << line << endl;
	}
      //      while(!stream.eof())
      for(int j= 0; j < Angles; j++)
	{
	  stream  >> angle >> res >> reserror ;
	  if(print)  cout << "line " << j+7 << endl;

	  freshresA[j] = res;
	  freshresA_err[j] = reserror;
	  if(print)      cout  << angle << " " << freshresA[j] << " " << freshresA_err[j] << endl;
	}//angles
    }// file scanning

  
  Double_t Angle[Angles];
  Double_t AngleError[Angles];
  Double_t TanAngle[Angles];
  Double_t TanAngleError[Angles];
  Double_t Resolution[Angles];
  Double_t RMS[Angles];
  Double_t Ncol[Angles];
  Double_t NcolError[Angles];
  Double_t ResolutionError[Angles];
  Double_t RMSError[Angles];

  TString ss_Angle[Angles];

  TH1F * h_res;//[Angles];

  Int_t run=2789;
  TString Run[Angles];
  ostringstream strs[Angles];

  TString label = "closest_A12C15_bestproton";

  TString ss_dphcut = "15";
  TString ss_perc[Cuts];
  double sigma[Angles][Cuts] ;
  double sigmaerr[Angles][Cuts] ;
  double Percentage[Angles][Cuts] ;
  
  int high[Angles][Cuts];
  TString inputfile;
  
  float angleunit = 1.25;

  TString pitch = "25";
  float perc[Cuts];        
  for(int i=0; i<Angles; i++)  {
    if(print) cout << "Getting files "  << endl;
    if( i<11)      {
      Run[i].Form("%d",run+i);
      if(print) cout << "Run: " << Run << endl;
    }
    else      {
	Run[i].Form("%d",run+i+1);
	if(print) cout << "Run: " << Run << endl;
      }



	Angle[i] = -4+i;
	Angle[i]=Angle[i]*angleunit; //grad
	ss_Angle[i].Form("%f",Angle[i]);

	TanAngle[i]=TMath::Tan(Angle[i]*TMath::Pi()/((Double_t)180.));
	cout << " Angle " <<  Angle[i] << " TanAngle " <<  TanAngle[i] << endl;
	AngleError[i]=1; //grad
	TanAngleError[i]=(1+TMath::Power(TMath::Tan(AngleError[i]*TMath::Pi()/((Double_t) 180.)),2))*AngleError[i]*TMath::Pi()/((Double_t) 180.); //grad
	
    if(print) cout << "Angle " << i<< ": " << Angle[i] << " grads" << endl;
    

    inputfile = inputDir+"drei-r"+Run[i]+"_irene_dphcutB"+ss_dphcut+"_"+label+".root";
    cout << inputfile << endl;
    TFile * file = new TFile(inputfile);
    TTree * tree = (TTree*)file->Get("charge_res");
    tree->Print();
    
    TH1I * hdx3tree = new TH1I("hdx3tree", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
    tree->Draw("dx3tree>>hdx3tree","","goff");
    hdx3tree = (TH1I*)gDirectory->Get("hdx3tree");
    
    double integral;
    double integral_per[Cuts];
  
  
  
    TString dir = "straightTracksY_isoAandCandB_straightTracksX";
    TH1I * hclphB;
    TDirectory * histodir =(TDirectoryFile*)file->Get(dir);
    hclphB= (TH1I*)histodir->Get("clphB");
    

    
    integral = hclphB->Integral(0,hclphB->GetNbinsX()+1);
    if(print)      cout << " integral " << integral << " entries " << hclphB->GetEntries()<<  endl;

    for(int k = 0; k< Cuts; k++){
      perc[k] = (float) (k+1) * 0.05 ;
      cout << " perc " << perc[k] << endl;

      integral_per[k] = perc[k] *integral; //careful!! you should take into account the peak at low value! Hist is now filled only for isolated clusters and the effect is reduced
      cout << " integral " << integral_per[k] << endl;
	
      high[i][k] = 0;
	
      double integrating = 0;
      int m = 1;
      while(integrating<integral_per[k])    {
	  
	//if(print) cout << " while "<< i << endl;
	integrating= hclphB->Integral(1,m); //hclph[j]->GetNbinsX()-i);
	high[i][k] = hclphB->GetBinCenter(m); //hclph[j]->GetNbinsX()-i);
	  //if(print) cout << " integral " << integrating << " high " << high[i][j][k] << endl;
	m++;
      }
      cout << " integral "<< perc << " " << integral_per[k] << " obtained " << integrating << endl;
    }


    

    for(int k = 0; k< Cuts; k++){
      cout << " perc " << (float) (k+1) * 0.05 << " integral " << integral_per[k] << " high bin " << high[i][k] << " integral high " << hclphB->Integral(1,hclphB->FindBin(high[i][k])) << endl;
	  hclphB->GetXaxis()->SetRange(1,hclphB->FindBin(high[i][k]));
	  average[k] = hclphB->GetMean();
	  hclphB->GetXaxis()->SetRange(1,hclphB->GetNbinsX());
    }
  
    TString ph[Cuts]; 
    for(int k = 0; k< Cuts; k++){    
      ph[k].Form("%d",high[i][k]);
    }


  TFile *MyFile = new TFile(inputDir+"dx3_landauAddition_"+ss_Angle[i]+"_dphcut"+ss_dphcut+"_"+label+".root","RECREATE"); 
  //plot landaus
  TCanvas *cFDB2  = new TCanvas("cFDB2", "FDB resolution", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //    gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  hclphB->Draw();
  TLine *  linea[Cuts];
  for(int k = 0; k< Cuts; k++){
    linea[k]    = new TLine( high[i][k],0.,high[i][k],2000.);
    linea[k]->SetLineColor(kRed);
    linea[k]->SetLineWidth(2);
    linea[k]->SetLineStyle(2);
    linea[k]->Draw("same");
  }
  cFDB2->Write("landau_B");



  
  
  TH1I * hdx3ph[Cuts];
  
  TString totbins;
  int maximumedge = hclphB->GetBinLowEdge(hclphB->GetNbinsX()+1)+hclphB->GetBinWidth(1);
  totbins.Form("%d",maximumedge);
  cout <<" totbins "<< totbins<< endl;


  for(int k = 0; k< Cuts; k++){
    if(i==0) ss_perc[k].Form("%d",perc[k]);
    TH1I * hdx3treeph = new TH1I("hdx3treeph", "triplet dx3 ; dx [mm];triplets", 500, -0.5, 0.5 );
    hdx3ph[k] = new TH1I("hdx3treeph", "triplet dx3  ; dx [mm];triplets", 500, -0.5, 0.5 );

    TString phB_high;
    phB_high=ph[k];


    cout << " cut " << k << endl;
    cout << " B " << phB_high << endl;
    
    tree->Draw("dx3tree>>hdx3treeph","clphBiiitree<"+phB_high);
    //tree->Draw("dx3tree>>hdx3treeph","clphBiiitree>"+phB_low+"&&clphBiiitree<"+phB_high) ; 
    //hdx3treeph =
    hdx3ph[k] = (TH1I*)gDirectory->Get("hdx3treeph");
    
    hdx3ph[k]->SetDirectory(0);
    cout << "filled" << endl;
    hdx3treeph->Delete();

    cout << "deleted" << endl;
    cout << hdx3ph[k]->GetTitle() << " entries " << hdx3ph[k]->GetEntries() << endl;
    hdx3ph[k]->GetXaxis()->SetRangeUser(-0.05,0.05);
    hdx3ph[k]->Write("hdx3ph_"+ss_perc[k]);
    hdx3ph[k]->GetXaxis()->SetRangeUser(-0.5,0.5);

  } // cuts

  MyFile->Close();

  
  // now fill the rms & res
  
  for(int k = 0; k< Cuts; k++){

  sigma[i][k] = hdx3ph[k]->GetRMS() * 1000;
  sigmaerr[i][k] = hdx3ph[k]->GetRMSError() * 1000;


  cout << " RMS " << sigma[i][k] << " ± " << sigmaerr[i][k] << endl;
  FitTH1(hdx3ph[k], &(sigma[i][k]), &(sigmaerr[i][k]), ss_perc[k], detectorA, detectorB, detectorC, function,&(Percentage[i][k]) );
  ExtractRes(&(sigma[i][k]), &(sigmaerr[i][k]), true,freshresA[i],freshresA_err[i]);
  cout << " RES " << sigma[i][k] << " ± " << sigmaerr[i][k] << endl;

  }

  }//Angles


  /*
  
  Double_t edges[Cuts + 2] = {0.0, (Double_t) high[i][1][0], (Double_t)high[i][1][1],(Double_t)high[i][1][2],(Double_t)high[i][1][3],(Double_t)high[i][1][4],1000.};
  // Bin 1 corresponds to range [0.0, 0.2]
  // Bin 2 corresponds to range [0.2, 0.3] etc...

  //TH1* h_res_vs_charge = new TH1D("h_res_vs_charge","Resolution vs cluster charge",Cuts,edges);


  for(int i=0; i<Cuts; i++)    {
      h_res_vs_charge->SetBinContent(i+1,sigma[i]);
      h_res_vs_charge->SetBinError(i+1,sigmaerr[i]);
    }
  */
  


  for(int i=0; i<Angles; i++){
    
    ofstream myfile;
    myfile.open ("/home/zoiirene/Output/TextFiles/Ascan_resTree_"+detectorB+"_dphcutB"+ss_dphcut+"_"+function+"_"+label+"_"+ss_Angle[i]+".txt");
    myfile << "A " << detectorA << "\n";
    myfile << labelA << "\n";
    myfile << "B " << detectorB<< "\n";
    myfile << labelB << "\n";
    myfile << "C " << detectorC<< "\n";
    myfile << labelC << "\n";
    myfile << info << "\n";
    myfile << "phBcenter RES(um) Error High Percentage\n";

    Double_t charge[Cuts] = {
     (Double_t) high[i][0]/2.,
     ((Double_t)high[i][1]+(Double_t) high[i][0])/2.,
     ((Double_t)high[i][1]+(Double_t)high[i][2])/2.,
     ((Double_t)high[i][3]+(Double_t)high[i][2])/2.,
     ((Double_t)high[i][4]+(Double_t)high[i][3])/2.
    };
    TGraph * h_res_vs_charge = new TGraph(Cuts,charge,sigma[i]);    
 
    TCanvas *c2 = new TCanvas("c2", "resolution vs dphcut", 1500, 900);
    gPad->SetTicks(1,1);
    gROOT->SetStyle("Plain");
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    

    h_res_vs_charge->SetTitle(" ");
    h_res_vs_charge->GetYaxis()->SetTitle("Resolution [#mum]");
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
      linea[k]    = new TLine( high[i][k],0.,high[i][k],6.);
      linea[k]->SetLineColor(kGray);
      linea[k]->SetLineWidth(2);
      linea[k]->SetLineStyle(2);
      linea[k]->Draw("same");
      myfile << " " << average[k] << " " <<sigma[i][k] << " " << sigmaerr[i][k] << " " << high[i][k]  << " " << perc[k] <<"\n";
    }
    myfile.close();
    
    
    TString name = outputDir+"Res_vs_clusterchergeAddition_ABC_"+function+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+"_run"+Run[i];
    c2->SaveAs(name+".eps");
    c2->SaveAs(name+".png");
    c2->SaveAs(name+".pdf");
    c2->Delete();
  }//Angles




  
}//resolution 


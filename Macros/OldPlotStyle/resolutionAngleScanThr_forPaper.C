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
#define anglesNONirr 29
#define anglesPirr2 28
#define anglesNirr4 9
#define irradiations 3
#define dphcuts 6


void TDR();
void TDR2(TCanvas * c_all, int period=0, int pos= 11);

void resolutionAngleScanThr_forPaper(TString thr="500",TString name = "preliminary_SmallThrScan", TString func = "RMSself", bool unfolding = false){
  
  TString inputDir="/home/zoiirene/Output/TextFiles/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  int measurements[irradiations];
  measurements[0] = anglesNONirr;
  measurements[1] = anglesPirr2;
  measurements[2] = anglesNirr4;

  /*
  TString ss_dphcutPerc[dphcuts]={"4","6","9","12","15","18","20","25","30"};
  TString ss_dphcut_0[dphcuts]={"8","12","18","23","29","35", "38", "48", "57"};
  TString ss_dphcut_1[dphcuts]={"5", "8", "11", "15", "19", "23", "25", "31", "38"};
  TString ss_dphcut_2[dphcuts]={"6", "10", "15", "19", "24", "29", "32", "40", "48"};
  */
  
  TString ss_dphcutPerc[dphcuts]={"4","6","9","12","20","30"};

  TString ss_dphcut_0[dphcuts]={"8","12","18","23", "38", "57"};
  TString ss_dphcut_1[dphcuts]={"5", "8", "11", "15", "25", "38"};
  TString ss_dphcut_2[dphcuts]={"6", "10", "15", "19", "32", "48"};


  TString filenames_pt1[irradiations];
  filenames_pt1[0] = "Ascan_resTree_148_dphcutB"; 
  filenames_pt1[1] = "Ascan_resTreeCorr_120i_dphcutB";
  filenames_pt1[2] = "Ascan_194i_dphcutB"; 
  TString filenames_pt2[irradiations];
  filenames_pt2[0] = "_RMSselfthrScan_thrScan_A13C14.txt"; 
  filenames_pt2[1] = "_RMSself_RMSself6sig_closest_A13C14_bestnonirr.txt"; 
  filenames_pt2[2] = "_RMSself_800_resTreeCorr_thrScan_A12C13.txt"; //5.2 GeV


  TString filenames2_pt1[irradiations];
  filenames2_pt1[0] = "Ascan_clsizeB_148_dphcutB";
  filenames2_pt1[1] = "Ascan_clsizeB_120i_dphcutB";
  filenames2_pt1[2] = "Ascan_clsizeB_194i_dphcutB";
  TString filenames2_pt2[irradiations];
  filenames2_pt2[0] = "_RMSself_thrScan_A13C14.txt";
  filenames2_pt2[1] = "__RMSself_thrScan_A12C15.txt";
  filenames2_pt2[2] = "_800_thrScan_A12C13.txt";

  
  TString irr[irradiations];
  irr[0] = " Non-irradiated, 120 V, 5.6 GeV"; // "no irr, 5.6 GeV";
  //  irr[1] = "proton irr at #phi_{eq}=2#times10^{15} cm^{-2}, 5.6 GeV";
  irr[1] = "#phi_{eq} = 2.1#times10^{15} cm^{-2}, proton, 800 V, 5.6 GeV";
  irr[2] = "#phi_{eq} = 4#times10^{15} cm^{-2}, neutron, 800 V, 5.2 GeV";

  TString ss_irr[irradiations] = {"nonirr","2E15prot","4E15neut"};
  double angle_0[dphcuts][anglesNONirr];
  double angle_1[dphcuts][anglesPirr2];
  double angle_2[dphcuts][anglesNirr4];

  double angleerr_0[dphcuts][anglesNONirr];
  double angleerr_1[dphcuts][anglesPirr2];
  double angleerr_2[dphcuts][anglesNirr4];

  double res_0[dphcuts][anglesNONirr];
  double res_1[dphcuts][anglesPirr2];
  double res_2[dphcuts][anglesNirr4];

  double reserr_0[dphcuts][anglesNONirr];
  double reserr_1[dphcuts][anglesPirr2];
  double reserr_2[dphcuts][anglesNirr4];


  
  double angle;
  double res,reserror,angleerr;
  int lines[irradiations] = {8,8,7}; //lines to be skipped in input files before reading the measurements
  
  for(int i =0; i< irradiations; i++)    {
    cout << irr[i] << endl;
    for(int l = 0; l< dphcuts ; l++){
      cout << "dphcut "<< ss_dphcutPerc[l] << endl;
      TString  ss_dphcut=ss_dphcut_0[l];
      if(i==1) ss_dphcut=ss_dphcut_1[l];
      if(i==2) ss_dphcut=ss_dphcut_2[l];

      TString filename      = inputDir+filenames_pt1[i]+ss_dphcut+filenames_pt2[i];
      cout << filename << endl;
      ifstream stream(filename);
      
      std::string line;
      if(!stream.is_open())	{
	  cout << " File " << filename << " not opened" << endl;
      }
      else{
	for(int k =0; k<lines[i];k++) {
	  std::getline(stream,line);
	  if(print) cout << k << " line " << line << endl;
	}
	for(int j= 0; j < measurements[i]; j++){
	  stream  >> angle >> res >> reserror ;
	  if(print)	 cout << "line " << j+lines[i]-1 << endl;
	  if(i==0){
	    angle_0[l][j] = angle;
	    angleerr_0[l][j] = 1;
	    res_0[l][j] = res;
	    reserr_0[l][j] = reserror;
	    if(print)	 cout  << angle_0[l][j] << " " << res_0[l][j] << " " << reserr_0[l][j] << endl;
	  }
	  if(i==1){
	    angle_1[l][j] = angle;
	    angleerr_1[l][j] = 1;
	    res_1[l][j] = res;
	    reserr_1[l][j] = reserror;
	    if(print)	 cout  << angle_1[l][j] << " " << res_1[l][j] << " " << reserr_1[l][j] << endl;
	  }
	  if(i==2){
	    angle_2[l][j] = angle;
	    angleerr_2[l][j] = 1;
	    res_2[l][j] = res;
	    reserr_2[l][j] = reserror;
	    if(print)	 cout  << angle_2[l][j] << " " << res_2[l][j] << " " << reserr_2[l][j] << endl;
	  }

	}//measurements
      }//file scanning
    }//cuts
  }//irradiations



   //clsize file for IEEE
  double clsize_0[dphcuts][anglesNONirr];
  double clsize_1[dphcuts][anglesPirr2];
  double clsizeerr_0[dphcuts][anglesNONirr];
  double clsizeerr_1[dphcuts][anglesPirr2];
  double clsize_4[dphcuts][anglesNirr4]; // neutron irr 5.2 GeV
  double clsizeerr_4[dphcuts][anglesNirr4]; //mean error
  
  double clsize,clsizeerr;
  int lines2[irradiations] = {8,8,8};
  
  for(int i =0; i< irradiations; i++)  {
    cout << irr[i] << endl;
    for(int l = 0; l< dphcuts ; l++){
      cout << "dphcut "<< ss_dphcutPerc[l] << endl;
      TString  ss_dphcut=ss_dphcut_0[l];
      if(i==1) ss_dphcut=ss_dphcut_1[l];
      if(i==2) ss_dphcut=ss_dphcut_2[l];

      TString filename      = inputDir+filenames2_pt1[i]+ss_dphcut+filenames2_pt2[i];
      cout << filename << endl;
      ifstream stream(filename);
      
      
      std::string line;
      if(!stream.is_open()){
	cout << " File " << filename<< " not opened" << endl;
      }
      else{
	for(int k =0; k<lines2[i];k++){
	  std::getline(stream,line);
	  if(print) cout << k << " line " << line << endl;
	}
	cout << " number of measurements " << measurements[i] << endl;
	for(int j= 0; j < measurements[i]; j++){
	  stream  >> angle >> clsize >> clsizeerr ;
	  if(print)	 cout << "line " << j+lines2[i]-1 << endl;
	  if(i==0){
	    clsize_0[l][j] = clsize;
	    clsizeerr_0[l][j] = clsizeerr;
	    if(print)	 cout  << angle_0[l][j] << " " << clsize_0[l][j] << " " << clsizeerr_0[l][j] << endl;
	  }
	  if(i==1){
	    clsize_1[l][j] = clsize;
	    clsizeerr_1[l][j] = clsizeerr;
	    if(print)	 cout  << angle_1[l][j] << " " << clsize_1[l][j] << " " << clsizeerr_1[l][j] << endl;
	  }
	  if(i==2){
	    clsize_4[l][j] = clsize;
	    clsizeerr_4[l][j] = clsizeerr;
	    if(print)	 cout  << angle_2[l][j] << " " << clsize_4[l][j] << " " << clsizeerr_4[l][j] << endl;
	  }

	}//measurements
      }//file scanning
    }//dphcuts
  }//IEEE



  /////// plots!
  int symbols_res[dphcuts]={20,21,22,23,29,33}; //,34,39,43};
  int symbols_cls[dphcuts]={24,25,26,32,30,27};//,28,37,42};
  int colors_0[dphcuts]={625,632,626,634,630,627};//,635,631,636};
  int colors_1[dphcuts]={14,13,12,1,620,615};//,619,618,611};
  int colors_2[dphcuts]={591,593,600,602,598,603};//,599,604,1};
  
  TCanvas *cFDB2[irradiations];

  for(int i = 0; i<irradiations; i++){
    cFDB2[i]  = new TCanvas("c3", "FDB resolution", 600, 600);
    gPad->SetTicks(1,1);
    gROOT->SetStyle("Plain");
    cFDB2[i]->SetLeftMargin(1.);
    cFDB2[i]->SetTopMargin(1.);
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    pad2->SetFillStyle(4000); //will be transparent
    pad2->SetFrameFillStyle(0);

    pad1->Draw();
    pad1->cd();
    TLegend* legFDB2 = new TLegend(0.15,0.65,0.45,0.8);
    legFDB2->SetNColumns(3);
    legFDB2->SetLineColor(0);
    legFDB2->SetTextSize(0.035);
    TLegend* leg = new TLegend(0.05,0.8,0.45,0.85);
    leg->SetLineColor(0);
    leg->SetTextSize(0.035);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    //TPaveText *pt = new TPaveText(.15,.8,.65,.85,"NDC");
    //
    
    //pt->Draw();
    TGraphErrors* resolutionPlot[dphcuts];
    for (int l =0; l < dphcuts; l++){
      resolutionPlot[l] = new TGraphErrors(anglesNONirr,angle_0[l],res_0[l],angleerr_0[l],reserr_0[l]);
      int color = colors_0[l];
      if(i==1){
	resolutionPlot[l] = new TGraphErrors(anglesPirr2,angle_1[l],res_1[l],angleerr_1[l],reserr_1[l]);
	color = colors_1[l];
      }
      if(i==2){
	resolutionPlot[l] = new TGraphErrors(anglesNirr4,angle_2[l],res_2[l],angleerr_2[l],reserr_2[l]);
	color = colors_2[l];
      }   
      resolutionPlot[l]->SetTitle(" ");

      resolutionPlot[l]->SetMarkerSize(1.);
      resolutionPlot[l]->SetMarkerColor(color);
      resolutionPlot[l]->SetLineColor(color);
      resolutionPlot[l]->SetMarkerStyle(symbols_res[l]);
      if(l==0){
	leg->AddEntry(resolutionPlot[l],irr[i],"");

	resolutionPlot[0]->GetYaxis()->SetTitle("Single hit resolution [#mum]");
	resolutionPlot[0]->GetXaxis()->SetTitle("#theta [deg]");
	resolutionPlot[0]->GetXaxis()->SetLimits(0.,31.);
	resolutionPlot[0]->GetYaxis()->SetRangeUser(0.,15.);
	resolutionPlot[0]->Draw("AEP");
      }	  
      else
	resolutionPlot[l]->Draw("EPsame");
      

  
      legFDB2->AddEntry(resolutionPlot[l],ss_dphcutPerc[l]+" %","lp");
    }
    legFDB2->Draw();
    leg->Draw();
    pad2->Draw();
    pad2->cd();
    TGraphErrors* clsizePlot[dphcuts];
    for (int l =0; l < dphcuts; l++){
      clsizePlot[l] = new TGraphErrors(anglesNONirr,angle_0[l],clsize_0[l],angleerr_0[l],clsizeerr_0[l]);
      int color = colors_0[l];
      if(i==1){
	clsizePlot[l] = new TGraphErrors(anglesPirr2,angle_1[l],clsize_1[l],angleerr_1[l],clsizeerr_1[l]);
	color = colors_1[l];
      }
      if(i==2){
	clsizePlot[l] = new TGraphErrors(anglesNirr4,angle_2[l],clsize_4[l],angleerr_2[l],clsizeerr_4[l]);
	color = colors_2[l];
      }
      clsizePlot[l]->SetTitle(" ");
      clsizePlot[l]->SetMarkerColor(color);
      clsizePlot[l]->SetLineColor(color);
      clsizePlot[l]->SetMarkerStyle(symbols_cls[l]);
      clsizePlot[l]->SetMarkerSize(1.);
      
      if(l==0){
	clsizePlot[0]->GetYaxis()->SetTitle("Average cluster size");
	clsizePlot[0]->GetXaxis()->SetTitle("#theta [deg]");
	clsizePlot[0]->GetXaxis()->SetLimits(0.,31.);
	clsizePlot[0]->GetYaxis()->SetRangeUser(0.,15.);
	clsizePlot[0]->Draw("AEPY+");
      }
      else
	clsizePlot[l]->Draw("EPsame");
    }
    TString  outname = outputDir+"ResolutionSummaryPaper_angleScans25_"+name+"_"+ss_irr[i];
    cFDB2[i]->SaveAs(outname+".eps");
    cFDB2[i]->SaveAs(outname+".png");
    cFDB2[i]->SaveAs(outname+".pdf");
    cFDB2[i]->SaveAs(outname+".root");
    cFDB2[i]->SaveAs(outname+".C");
  }











  
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

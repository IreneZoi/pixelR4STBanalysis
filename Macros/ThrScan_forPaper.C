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
#include "resolution.h"

bool print=true;
using namespace std;
#define anglesNONirr 29
#define anglesNONirrD 5
#define anglesPirr2 28
#define anglesNirr4 9
#define anglesNirr4D 4
#define irradiations 3
#define dphcuts 35


void TDR();
void TDR2(TCanvas * c_all, int period=0, int pos= 11);

void ThrScan_forPaper(TString thr="500",TString name = "preliminary_ThrScan", TString func = "RMSself", bool unfolding = false){
  
  TString inputDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  int measurements[irradiations+1];
  measurements[0] = anglesNONirr;
  measurements[1] = anglesPirr2;
  measurements[2] = anglesNirr4;
  measurements[3] = 1; 

 
  double i_dphcutPerc[dphcuts]; //={0, 0.5,1,1.6,2,2.7,3.3,3.8,  4 , 6 , 9 , 12 , 15 , 18 , 20 , 25 , 30 };
  TString ss_dphcutPerc[dphcuts];//={"0","0.5","1","1.6","2","2.7","3.3","3.8","4","6","9","12","15","18","20","25","30"};
  int i_dphcut_0[dphcuts]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,23,26,29,35,38,48,57,60,70,80,90,100,110,120,130};
  TString ss_dphcut_0[dphcuts]; //={"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","20","23","26","29","35", "38", "48", "57","60","70","80","90","100","110","120","130"};
  TString ss_dphcut_1[dphcuts]={"5", "8", "11", "15", "19", "23", "25", "31", "38"};
  TString ss_dphcut_2[dphcuts]={"6", "10", "15", "19", "24", "29", "32", "40", "48"};
  double err_0[dphcuts];
  TString label = "beamdiv_A13C14";
  TString run="2743";
  double MPV=183;
  double events[dphcuts];
  double clsize[dphcuts];
  double events_err[dphcuts];
  double clsize_err[dphcuts];
  double Resolution[dphcuts];
  double ResolutionError[dphcuts];
  double Percentage[dphcuts];
  double Min[dphcuts];
  double Max[dphcuts];

  
  for(int i=0; i<dphcuts; i++)  {
    err_0[i]=0;
    i_dphcutPerc[i]=(double)i_dphcut_0[i]/MPV*100.;
    ss_dphcutPerc[i].Form("%f",i_dphcutPerc[i]);
    ss_dphcut_0[i].Form("%d",i_dphcut_0[i]);
    TString inputfile = inputDir+"drei-r"+run+"_irene_dphcutB"+ss_dphcut_0[i]+"_"+label+".root";
    cout << inputfile << endl;
    TFile * file = new TFile(inputfile);
    TString dir = "straightTracksY_isoAandCandB_straightTracksX";
    TH1I * hnrowB;
    TDirectory * histodir =(TDirectoryFile*)file->Get(dir);
    hnrowB = (TH1I*)histodir->Get("nrowB");
    clsize[i]=hnrowB->GetMean();
    clsize_err[i]=hnrowB->GetRMS();
    events[i]=hnrowB->GetEntries();
    events_err[i]=TMath::Sqrt(events[i]);
    cout << " thr "<< ss_dphcutPerc[i] << " entries " << events[i] << " clsize " << clsize[i] << endl;
    TH1I * hres;
    hres = (TH1I*)file->Get("dx3_clphABC90evR");
    FitTH1(hres, &(Resolution[i]), &(ResolutionError[i]), ss_dphcutPerc[i], "A", "148", "C", func, &(Percentage[i]),&(Min[i]),&(Max[i]));
    ExtractRes(&(Resolution[i]),&(ResolutionError[i]));



  }


  TString legend = " Non-irradiated, 120 V, 5.6 GeV, #theta= 8.8 deg";




  ///////////////////////////////
  //now compare threshold and cluster size at best angle


  TCanvas *cc  = new TCanvas("cc", "FDB resolution", 600, 600);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  cc->SetLeftMargin(-2);
  cc->SetRightMargin(-10);
  cc->SetTopMargin(1.);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(15);
  TLegend* legFDB2 = new TLegend(0.3,0.3,0.8,0.5);
  legFDB2->SetNColumns(4);
  legFDB2->SetLineColor(0);

  cout << " initialized new canvas" << endl;
  TPad *pad2 = new TPad("pad2","",0,0.,1,0.45);
  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.4);
  pad2->SetLeftMargin(0.15);
  pad2->SetRightMargin(0.1);
  pad2->Draw();

  TPad *pad1 = new TPad("pad1","",0,0.45,1,1);
  pad1->Draw();
  pad1->SetBottomMargin(0.03);
  pad1->SetRightMargin(0.1);
  pad1->SetLeftMargin(0.15);

  pad1->cd();
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(15);

  gPad->SetTicks(1,1);
  TGraphErrors* Evts = new TGraphErrors(dphcuts,i_dphcutPerc,events,err_0,events_err);
  Evts->GetYaxis()->SetTitle("Reconstructed clusters");  
  Evts->GetXaxis()->SetTitleSize(0); // labels will be 14 pixels
  Evts->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels
  Evts->GetYaxis()->SetTitleFont(43); // labels will be 14 pixels
  Evts->GetYaxis()->SetLabelFont(43); // labels will be 14 pixels
  Evts->GetYaxis()->SetTitleSize(15); // labels will be 14 pixels
  Evts->GetYaxis()->SetTitleOffset(2.5); // labels will be 14 pixels
  Evts->GetYaxis()->SetLabelSize(15); // labels will be 14 pixels
  Evts->GetXaxis()->SetLimits(-1,75); // labels will be 14 pixels
  Evts->GetYaxis()->SetRangeUser(0,50000.);
  //  gPad->SetLogy();
  Evts->Draw("AEP");
  legFDB2->AddEntry(Evts,legend,"");
  legFDB2->Draw();
  TLine *  lineb = new TLine( 12./MPV*100,0.,12./MPV*100,50000);
  lineb->SetLineColor(kGray);
  lineb->SetLineWidth(2);
  lineb->SetLineStyle(2);
  lineb->Draw("same");



  pad2->cd();
  TGraphErrors* clsizePlot = new TGraphErrors(dphcuts,i_dphcutPerc,clsize,err_0,clsize_err);

  gPad->SetTicks(1,1);


  clsizePlot->GetYaxis()->SetTitle("Average cluster size");
  clsizePlot->GetXaxis()->SetTitle("Threshold [%]");

  clsizePlot->GetYaxis()->SetNdivisions(5,0,5);
			       
  clsizePlot->GetXaxis()->SetTitleFont(43);
  clsizePlot->GetXaxis()->SetTitleSize(15); // labels will be 14 pixels
  clsizePlot->GetXaxis()->SetTitleOffset(3); // labels will be 14 pixels
  clsizePlot->GetXaxis()->SetLabelFont(43);
  clsizePlot->GetXaxis()->SetLabelSize(15); // labels will be 14 pixels

  clsizePlot->GetYaxis()->SetTitleFont(43);
  clsizePlot->GetYaxis()->SetTitleSize(15); // labels will be 14 pixels
  clsizePlot->GetYaxis()->SetLabelFont(43);
  clsizePlot->GetYaxis()->SetLabelSize(15); // labels will be 14 pixels
  clsizePlot->GetXaxis()->SetLimits(-1.,75.);
  clsizePlot->GetYaxis()->SetRangeUser(0.,12.);

  
    clsizePlot->SetTitle(" ");
    clsizePlot->SetMarkerSize(1.);
    //clsizePlot->SetMarkerColor(colors_irr[l]);
    //clsizePlot->SetLineColor(colors_irr[l]);
    //clsizePlot->SetMarkerStyle(marker_irr[l]);
    cout << " initialized styles" << endl;
    clsizePlot->Draw("AEP");
    
    TLine *  linea = new TLine( -1.,2.,75.,2.);
    linea->SetLineColor(kGray);
    linea->SetLineWidth(2);
    linea->SetLineStyle(2);
    linea->Draw("same");

    TLine *  linev = new TLine( 12./MPV*100,0.,12./MPV*100,12);
    linev->SetLineColor(kGray);
    linev->SetLineWidth(2);
    linev->SetLineStyle(2);
    linev->Draw("same");
  
    pad1->cd();
    cout << " back to pad 1" << endl;		  
    TString outname = outputDir+"ResolutionSummaryPaper_ThrScan_bestAngle_nonirr_"+name+"_clsize";
    cc->SaveAs(outname+".eps");
    cc->SaveAs(outname+".png");
    cc->SaveAs(outname+".pdf");
    cc->SaveAs(outname+".root");
    cc->SaveAs(outname+".C");

    
    ///////////////////////////////
  //now compare threshold and cluster size at best angle


  TCanvas *c2  = new TCanvas("c2", "FDB resolution", 600, 600);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  c2->SetLeftMargin(-2);
  c2->SetRightMargin(-10);
  c2->SetTopMargin(1.);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(15);
  TLegend* legFDB2r = new TLegend(0.3,0.15,0.8,0.25);
  legFDB2r->SetNColumns(4);
  legFDB2r->SetLineColor(0);

  cout << " initialized new canvas" << endl;
  TPad *pad2r = new TPad("pad2r","",0,0.,1,0.45);
  pad2r->SetTopMargin(0.03);
  pad2r->SetBottomMargin(0.4);
  pad2r->SetLeftMargin(0.15);
  pad2r->SetRightMargin(0.1);
  pad2r->Draw();

  TPad *pad1r = new TPad("pad1r","",0,0.45,1,1);
  pad1r->Draw();
  pad1r->SetBottomMargin(0.03);
  pad1r->SetRightMargin(0.1);
  pad1r->SetLeftMargin(0.15);

  pad1r->cd();
  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(15);

  gPad->SetTicks(1,1);
  TGraphErrors* Res = new TGraphErrors(dphcuts,i_dphcutPerc,Resolution,err_0,ResolutionError);
  Res->GetYaxis()->SetTitle("#sigma^{0}_{x}");  
  Res->GetXaxis()->SetTitleSize(0); // labels will be 14 pixels
  Res->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels
  Res->GetYaxis()->SetTitleFont(43); // labels will be 14 pixels
  Res->GetYaxis()->SetLabelFont(43); // labels will be 14 pixels
  Res->GetYaxis()->SetTitleSize(15); // labels will be 14 pixels
  Res->GetYaxis()->SetTitleOffset(2.5); // labels will be 14 pixels
  Res->GetYaxis()->SetLabelSize(15); // labels will be 14 pixels
  Res->GetXaxis()->SetLimits(-1,75); // labels will be 14 pixels
  Res->GetYaxis()->SetRangeUser(0,10.);
  Res->SetMarkerStyle(20);
  Res->SetMarkerSize(0.5);

  //  gPad->SetLogy();
  
  Res->Draw("AEP");
  legFDB2r->AddEntry(Res,legend,"");
  legFDB2r->Draw();
  TLine *  linebr = new TLine( 12./MPV*100,0.,12./MPV*100,10);
  linebr->SetLineColor(kGray);
  linebr->SetLineWidth(2);
  linebr->SetLineStyle(2);
  linebr->Draw("same");



  pad2r->cd();
  //TGraphErrors* clsizePlot = new TGraphErrors(dphcuts,i_dphcutPerc,clsize,err_0,clsize_err);

  gPad->SetTicks(1,1);


  clsizePlot->GetYaxis()->SetTitle("Average cluster size");
  clsizePlot->GetXaxis()->SetTitle("Threshold [%]");

  clsizePlot->GetYaxis()->SetNdivisions(5,0,5);
			       
  clsizePlot->GetXaxis()->SetTitleFont(43);
  clsizePlot->GetXaxis()->SetTitleSize(15); // labels will be 14 pixels
  clsizePlot->GetXaxis()->SetTitleOffset(3); // labels will be 14 pixels
  clsizePlot->GetXaxis()->SetLabelFont(43);
  clsizePlot->GetXaxis()->SetLabelSize(15); // labels will be 14 pixels

  clsizePlot->GetYaxis()->SetTitleFont(43);
  clsizePlot->GetYaxis()->SetTitleSize(15); // labels will be 14 pixels
  clsizePlot->GetYaxis()->SetLabelFont(43);
  clsizePlot->GetYaxis()->SetLabelSize(15); // labels will be 14 pixels
  clsizePlot->GetXaxis()->SetLimits(-1.,75.);
  clsizePlot->GetYaxis()->SetRangeUser(0.,5.);

  
    clsizePlot->SetTitle(" ");
    clsizePlot->SetMarkerStyle(20);
    clsizePlot->SetMarkerSize(0.5);
    //clsizePlot->SetMarkerColor(colors_irr[l]);
    //clsizePlot->SetLineColor(colors_irr[l]);
    //clsizePlot->SetMarkerStyle(marker_irr[l]);
    cout << " initialized styles" << endl;
    clsizePlot->Draw("AEP");
    
    //TLine *  linea = new TLine( -1.,2.,75.,2.);
    linea->SetLineColor(kRed);
    //linea->SetLineWidth(2);
    //linea->SetLineStyle(2);
    linea->Draw("same");

    //TLine *  linev = new TLine( 12./MPV*100,0.,12./MPV*100,12);
    // linev->SetLineColor(kGray);
    //linev->SetLineWidth(2);
    //linev->SetLineStyle(2);
    linev->Draw("same");
  
    pad1r->cd();
    cout << " back to pad 1" << endl;		  
    outname = outputDir+"ResolutionSummaryPaper_ResThrScan_bestAngle_nonirr_"+name+"_clsize";
    c2->SaveAs(outname+".eps");
    c2->SaveAs(outname+".png");
    c2->SaveAs(outname+".pdf");
    c2->SaveAs(outname+".root");
    c2->SaveAs(outname+".C");











  
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

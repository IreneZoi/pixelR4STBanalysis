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
#define angles 3
#define irradiations 2
#define perc 20
void TDR();
void TDR2(TCanvas * c_all, int period=0, int pos= 11);

void resVSchargePerc_forPaper(){

  TString inputDir="/home/zoiirene/Output/TextFiles/";
  TString rootDir="/home/zoiirene/Output/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  TString filenames[irradiations][angles];
  TString rootfilenames[irradiations][angles];
  TString Angles[angles]={"0.000000","8.750000","27.500000"};
  TString AnglesNice[angles]={"0 deg","8.8 deg","27.5 deg"};
  TString baseNONirr="Ascan_resTree_148_dphcutB12_RMSself_dycut_A13C14_angle";
  TString basePirr2="Ascan_resTree_120i_dphcutB15_RMSself_dycut_A12C15_";
  TString baseNONirrROOT="_dphcut12_dycut_A13C14.root";
  TString basePirr2ROOT="_dphcut15_dycut_A12C15.root";
  for(int i = 0; i< angles; i++){
    filenames[0][i]=inputDir+baseNONirr+Angles[i]+".txt";
    filenames[1][i]=inputDir+basePirr2+Angles[i]+".txt";
    rootfilenames[0][i]=rootDir+"dx3_landauAddition_"+Angles[i]+baseNONirrROOT;
    rootfilenames[1][i]=rootDir+"dx3_landauAddition_"+Angles[i]+basePirr2ROOT;
    
  }

  TString irr[irradiations];
  irr[0] = " Non-irradiated, 120 V"; //, 5.6 GeV"; // "no irr, 5.6 GeV";
  //  irr[1] = "proton irr at #phi_{eq}=2#times10^{15} cm^{-2}, 5.6 GeV";
  irr[1] = "#phi_{eq} = 2.1 #times 10^{15} cm^{-2}, proton, 800 V"; //, 5.6 GeV";
  
  double charge_0[angles][perc];
  double charge_1[angles][perc];
  double chargeerr[angles][perc];

  double high_0[angles][perc];
  double high_1[angles][perc];

  double  Percentage[perc];
  double res_0[angles][perc];
  double res_1[angles][perc];
  double reserr_0[angles][perc];
  double reserr_1[angles][perc];


  double charge,high;
  double res,reserror,percentage;
  int lines[irradiations] = {8,8};

  TString filename;
  for(int i =0; i< irradiations; i++)    {
    for(int l = 0; l< angles; l++){
    
      filename= filenames[i][l];
      cout << filename << endl;
      ifstream stream(filename);

      std::string line;
      if(!stream.is_open())	{
	  cout << " File " << filename << " not opened" << endl;
	}      else	{
	for(int k =0; k<lines[i];k++)	    {
	  std::getline(stream,line);
	  if(print) cout << k << " line " << line << endl;
	}
	for(int j= 0; j < perc; j++)	    {
	  stream  >> charge >> res >> reserror >> high >> percentage;
	  if(print)  cout << "line " << j+lines[i]-1 << endl;
	  cout << " charge " << charge << " res " << res << " reserror " << reserror << " high " << high << " percentage " << percentage << endl;
	  if(i==0){
	    Percentage[j] = percentage;
	    charge_0[l][j] = charge;
	    chargeerr[l][j] = 0;
	    res_0[l][j] = res;
	    reserr_0[l][j] = reserror;
	    high_0[l][j] = high;
	    if(print)      cout  << charge_0[l][j] << " " << res_0[l][j] << " " << reserr_0[l][j] << " " <<Percentage[j]<<endl;
	  }
	  if(i==1){
	    charge_1[l][j] = charge;
	    res_1[l][j] = res;
	    reserr_1[l][j] = reserror;
	    high_1[l][j] = high;
	    if(print)      cout  << charge_1[l][j] << " " << res_1[l][j] << " " << reserr_1[l][j] << endl;
	  }

	}//measurements
      }//file scanning
    }//angles
  }//irradiations
  
  cout << " plotting "<< endl;

  /////// plots!
  TH1I * Landau[irradiations][angles];
  TH1F * fLandau[irradiations][angles];
  //  TGraph * gLandau[irradiations][angles];
  for(int l = 0; l < irradiations; l++){
    cout <<"%%%%%%%% irr "<< irr[l] << endl;
    for(int i = 0; i< angles; i++){
      cout << AnglesNice[i] << endl;
      TFile *file = new TFile(rootfilenames[l][i]);
      TCanvas *canvas = (TCanvas*)file->Get("landau_B");
      Landau[l][i] = (TH1I*)canvas->GetPrimitive("clphB");
      //cout << Landau[l][i]->GetEntries() << endl;
      Landau[l][i]->Rebin(2);
      //cout << " bins "<< Landau[l][i]->GetNbinsX() << endl;
      //cout << "maximum " << Landau[l][i]->GetBinCenter(Landau[l][i]->GetMaximumBin()) << endl;
      fLandau[l][i] = new TH1F("f","f",Landau[l][i]->GetNbinsX(),0.,1000.);
      //Landau[l][i]->Copy(fLandau[l][i]);
      //cout << fLandau[l][i].Integral() << endl;
      //fLandau[l][i].Scale(1./fLandau[l][i].Integral());
      //cout << fLandau[l][i].Integral() << endl;
      //Landau[l][i].Rebin(2);
      //Double_t x[Landau[l][i].GetNbinsX()];
      //Double_t y[Landau[l][i].GetNbinsX()];
      
      for(int j = 1; j <= Landau[l][i]->GetNbinsX(); j++ ){
        //cout << "BIN "<< j << " Bin center " << Landau[l][i]->GetBinCenter(j) << " bin content " << Landau[l][i]->GetBinContent(j) << endl;
        //cout << ((float) Landau[l][i]->GetBinContent(j))/((float) Landau[l][i]->GetEntries()) << endl;
	fLandau[l][i]->SetBinContent(j,((float) Landau[l][i]->GetBinContent(j))/((float) Landau[l][i]->GetEntries()));
	//cout << fLandau[l][i]->GetBinCenter(j) << " " << fLandau[l][i]->GetBinContent(j) << endl;
	//gLandau[l][i].SetPoint(j-1, (Double_t) Landau[l][i].GetBinCenter(j),(Double_t) Landau[l][i].GetBinContent(j));
	//x[j] = Landau[l][i].GetBinCenter(j);
	//y[j] = (double) Landau[l][i].GetBinContent(j)/(double)  Landau[l][i].GetEntries() ; ///Landau[l][i].GetBinWidth(j);
     	//cout << x[j]  << " " << y[j]  << endl;
      }
      
      //gLandau[l][i] = new TGraph(Landau[l][i].GetNbinsX(),x,y);
      //cout << "got TGraph " << endl;
    }
  }
  
  TCanvas *cFDB2[angles];
  
  for(int i = 0; i< angles; i++){

    cFDB2[i] = new TCanvas("cFDB2", "FDB resolution", 600, 600);
    
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
    TPad *pad2 = new TPad("pad2","",0,0.,1,0.5);
    pad2->SetBottomMargin(0.18);
    pad2->SetTopMargin(0.05);
    pad2->Draw();

    //pad2->SetFillStyle(4000); //will be transparent
    //#pad2->SetFrameFillStyle(0);
    TPad *pad1 = new TPad("pad1","",0,0.5,1,1);
    pad1->Draw();

    pad1->SetTopMargin(0.05);
    pad1->SetBottomMargin(0.2);

    pad2->cd();
    gPad->SetTicks(1,1);
    /*
    gLandau[0][i]->SetLineColor(kBlack);
    gLandau[0][i]->GetXaxis()->SetLimits(0.,1000.);
    //    gLandau[0][i]->GetXaxis()->SetMaxDigits(3);
    gLandau[0][i]->GetYaxis()->SetMaxDigits(3);
    gLandau[0][i]->GetYaxis()->SetTitle("Number of clusters");
    gLandau[0][i]->GetXaxis()->SetTitle("Cluster charge [ADC]");
    //gLandau[0][i]->SetLineStyle(2);
    //gLandau[1][i]->SetLineStyle(2);
    gLandau[0][i]->SetMarkerStyle(20);
    gLandau[0][i]->SetMarkerColor(kGreen+1);
    gLandau[1][i]->SetMarkerStyle(21);
    /*
    gLandau[0][i]->GetYaxis()->SetTitleFont(43);
    cout << " ######### title size " <<     gLandau[0][i]->GetYaxis()->GetTitleSize() << endl;
    gLandau[0][i]->GetYaxis()->SetTitleSize(1000);
    cout << " ######### title size " <<     gLandau[0][i]->GetYaxis()->GetTitleSize() << endl;
    //    gLandau[0][i]->GetYaxis()->SetTitleOffset(0.07);
    //gLandau[0][i]->GetYaxis()->SetLabelSize(0.07);
    gLandau[0][i]->GetXaxis()->SetTitleSize(0.07);
    gLandau[0][i]->GetXaxis()->SetTitleOffset(0.8);
    gLandau[0][i]->GetXaxis()->SetLabelSize(0.07);
    */



    fLandau[0][i]->SetLineColor(kBlack);
    fLandau[0][i]->SetLineWidth(2);
    fLandau[0][i]->GetXaxis()->SetRangeUser(0.,500.);
    //    fLandau[0][i]->GetXaxis()->SetLimits(0.,1000.);
    //fLandau[1][i]->GetXaxis()->SetLimits(0.,1000.);
    fLandau[0][i]->GetYaxis()->SetRangeUser(0.,0.16);

    //    fLandau[0][i]->GetXaxis()->SetMaxDigits(3);
    fLandau[0][i]->GetYaxis()->SetMaxDigits(3);
    fLandau[0][i]->GetYaxis()->SetTitle("Normalized number of clusters");
    fLandau[0][i]->GetXaxis()->SetTitle("Cluster charge [ADC]");

    fLandau[0][i]->GetXaxis()->SetTitleFont(43);
    fLandau[0][i]->GetXaxis()->SetTitleOffset(2.5);
    fLandau[0][i]->GetXaxis()->SetTitleSize(15); // labels will be 14 pixels

    fLandau[0][i]->GetYaxis()->SetTitleFont(43);
    fLandau[0][i]->GetYaxis()->SetTitleSize(15); // labels will be 14 pixels
    fLandau[0][i]->GetYaxis()->SetTitleOffset(1.8); //1.5);
    
    fLandau[0][i]->GetXaxis()->SetLabelFont(43);
    fLandau[0][i]->GetXaxis()->SetLabelSize(15); // labels will be 14 pixels
    fLandau[0][i]->GetYaxis()->SetLabelFont(43);
    fLandau[0][i]->GetYaxis()->SetLabelSize(15); // labels will be 14 pixels

    



    //fLandau[0][i]->SetLineStyle(2);
    //fLandau[1][i]->SetLineStyle(2);
    fLandau[0][i]->SetMarkerStyle(20);

    fLandau[1][i]->SetMarkerColor(kGreen+1);
    fLandau[1][i]->SetLineColor(kGreen+1);

    fLandau[1][i]->SetLineWidth(2);
    fLandau[1][i]->SetLineStyle(2);

    fLandau[1][i]->SetMarkerStyle(21);
    //fLandau[0][i]->GetXaxis()->SetLimits(0.,1000.);
    //    fLandau[0][i]->GetXaxis()->SetMaxDigits(3);
    

    
    fLandau[0][i]->Draw("hist");
    fLandau[1][i]->Draw("histsame");

    TLine *  lineano = new TLine( high_0[i][17],0.,high_0[i][17],0.1);
    lineano->SetLineColor(kGray);
    cout << " high_0[i][17] " <<  high_0[i][17] << endl;
    lineano->Draw("same");
    TLine *  lineairr = new TLine( high_1[i][17],0.,high_1[i][17],0.1);
    lineairr->SetLineColor(kGray);
    lineairr->SetLineStyle(2);
    lineairr->Draw("same");
    


    TLegend* leg = new TLegend(0.35,0.7,0.85,0.88);
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    //    leg->SetTextSize(0.05);
    resolutionPlot[0] = new TGraphErrors(perc,Percentage,res_0[i],chargeerr[i],reserr_0[i]);
    resolutionPlot[1] = new TGraphErrors(perc,Percentage,res_1[i],chargeerr[i],reserr_1[i]);
    //resolutionPlot[0] = new TGraphErrors(perc,charge_0[i],res_0[i],chargeerr[i],reserr_0[i]);
    //resolutionPlot[1] = new TGraphErrors(perc,charge_1[i],res_1[i],chargeerr[i],reserr_1[i]);

    leg->AddEntry(resolutionPlot[0],AnglesNice[i]+", 5.6 GeV","");

    for(int l =0; l < irradiations; l++)
      leg->AddEntry(fLandau[l][i],irr[l],"lp");

    leg->Draw();

    /*
    for( int j = 0 ; j < perc ; j++){
      cout << Percentage[j] << " charge " << charge_0[i][j] << " res " << res_0[i][j] << endl;
      cout << Percentage[j] << " charge " << charge_1[i][j] << " res " << res_1[i][j] << endl;
      
    }
    */
    
    //    pad2->Draw();
    pad1->cd();

    gPad->SetTicks(1,1);
    
    resolutionPlot[0]->SetTitle(" ");
    resolutionPlot[0]->GetYaxis()->SetTitle("#sigma_{x} [#mum]");
    resolutionPlot[0]->GetXaxis()->SetTitle("Cluster charge [%]");
    resolutionPlot[0]->SetMarkerSize(1.);


    resolutionPlot[0]->GetXaxis()->SetTitleFont(43);
    resolutionPlot[0]->GetXaxis()->SetTitleOffset(2.5);
    resolutionPlot[0]->GetXaxis()->SetTitleSize(15); // labels will be 14 pixels
    
    resolutionPlot[0]->GetYaxis()->SetTitleFont(43);
    resolutionPlot[0]->GetYaxis()->SetTitleSize(15); // labels will be 14 pixels
    resolutionPlot[0]->GetYaxis()->SetTitleOffset(1.8);
    
    resolutionPlot[0]->GetXaxis()->SetLabelFont(43);
    resolutionPlot[0]->GetXaxis()->SetLabelSize(15); // labels will be 14 pixels
    resolutionPlot[0]->GetYaxis()->SetLabelFont(43);
    resolutionPlot[0]->GetYaxis()->SetLabelSize(15); // labels will be 14 pixels
    resolutionPlot[0]->SetMarkerColor(kBlack);
    resolutionPlot[0]->SetLineColor(kBlack);
    resolutionPlot[0]->SetMarkerStyle(20);
    resolutionPlot[0]->GetXaxis()->SetLimits(0.,1.);
    resolutionPlot[0]->GetYaxis()->SetRangeUser(0.,15.);
    resolutionPlot[0]->Draw("AEP");
    
    resolutionPlot[1]->SetMarkerSize(1.);
    resolutionPlot[1]->SetMarkerColor(kGreen+1);
    resolutionPlot[1]->SetLineColor(kGreen+1);
    resolutionPlot[1]->SetMarkerStyle(21);
    resolutionPlot[1]->Draw("EPsame");



    cout << " plotted "<< endl;
    
    TLine *  linea9;
    linea9    = new TLine( 0.9,0.,0.9,15.);
    linea9->SetLineColor(kGray);
    linea9->SetLineStyle(3);
    linea9->Draw("same");
      
    
    
    //TLegend* legFDB2 = new TLegend(0.35,0.75,0.85,0.89);
    //legFDB2->SetLineColor(0);
    //legFDB2->SetFillStyle(0);
    //legFDB2->SetTextSize(0.035);
    //legFDB2->AddEntry(resolutionPlot[0],AnglesNice[i],"");
    //for(int i =0; i < irradiations; i++)
    //legFDB2->AddEntry(resolutionPlot[i],irr[i],"lp");

    //legFDB2->Draw();
    

    TString  outname = outputDir+"ResolutionSummaryPaper_angleScans25_charge_"+Angles[i];
    cFDB2[i]->SaveAs(outname+".eps");
    cFDB2[i]->SaveAs(outname+".png");
    cFDB2[i]->SaveAs(outname+".pdf");
    cFDB2[i]->SaveAs(outname+".root");
    cFDB2[i]->SaveAs(outname+".C");
    cout << "saved" << endl;
  }//angles
  cout << "done"<< endl;

  //////// NOW ADD COMPARISON of 0 and 27 NON IRR starting aroun the Landau peak
  cout << " comparison of NON IRR " << endl;
  TCanvas *cnonirr = new TCanvas("cnonirr", "resolution", 600, 600);
    
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
  
  TGraphErrors* resolutionPlotA[angles];
  TPad *pad2 = new TPad("pad2","",0,0.,1,0.5);
  pad2->SetBottomMargin(0.18);
  pad2->SetTopMargin(0.05);
  pad2->Draw();

  TPad *pad1 = new TPad("pad1","",0,0.5,1,1);
  pad1->Draw();

  pad1->SetTopMargin(0.05);
  pad1->SetBottomMargin(0.2);
  
  pad2->cd();
  gPad->SetTicks(1,1);

  fLandau[0][0]->SetLineColor(kBlack);
  fLandau[0][0]->SetLineWidth(2);
  fLandau[0][0]->GetXaxis()->SetRangeUser(0.,700.);
  fLandau[0][0]->GetYaxis()->SetRangeUser(0.,0.16);

  fLandau[0][0]->GetYaxis()->SetMaxDigits(3);
  fLandau[0][0]->GetYaxis()->SetTitle("Normalized number of clusters");
  fLandau[0][0]->GetXaxis()->SetTitle("Cluster charge [ADC]");

  fLandau[0][0]->GetXaxis()->SetTitleFont(43);
  fLandau[0][0]->GetXaxis()->SetTitleOffset(2.5);
  fLandau[0][0]->GetXaxis()->SetTitleSize(15); // labels will be 14 pixels
  
  fLandau[0][0]->GetYaxis()->SetTitleFont(43);
  fLandau[0][0]->GetYaxis()->SetTitleSize(15); // labels will be 14 pixels
  fLandau[0][0]->GetYaxis()->SetTitleOffset(1.8); //1.5);
    
  fLandau[0][0]->GetXaxis()->SetLabelFont(43);
  fLandau[0][0]->GetXaxis()->SetLabelSize(15); // labels will be 14 pixels
  fLandau[0][0]->GetYaxis()->SetLabelFont(43);
  fLandau[0][0]->GetYaxis()->SetLabelSize(15); // labels will be 14 pixels
    
    
  fLandau[0][0]->SetMarkerStyle(20);

  fLandau[0][2]->SetMarkerColor(kBlack);
  fLandau[0][2]->SetLineColor(kBlack);

  fLandau[0][2]->SetLineWidth(2);
  fLandau[0][2]->SetLineStyle(2);

  fLandau[0][2]->SetMarkerStyle(24);
    
  fLandau[0][0]->Draw("hist");
  fLandau[0][2]->Draw("histsame");

  /*
  TLine *  linea[irradiations];
    for(int i = 0; i<irradiations; i++){

	
      linea[i]    = new TLine( high_0[i][17],0.,high_0[i][17],0.1);
      linea[i]->SetLineColor(kBlack);
      if(i==1){
	linea[i]    = new TLine( high_1[i][17],0.,high_1[i][17],0.1);
	linea[i]->SetLineColor(kGreen+1);
      }
      //linea[i][k]->SetLineWidth(2);
      linea[i]->SetLineStyle(3);
      linea[i]->Draw("same");
      
    }
  */
  TLine *  linea0 = new TLine( high_0[0][17],0.,high_0[0][17],0.1);
  linea0->SetLineColor(kGray);
  linea0->Draw("same");
  TLine *  linea27 = new TLine( high_0[2][17],0.,high_0[2][17],0.1);
  linea27->SetLineColor(kGray);
  linea27->SetLineStyle(2);
  linea27->Draw("same");
  TLine *  linea0l = new TLine( high_0[0][4],0.,high_0[0][4],0.1);
  linea0l->SetLineColor(kGray);
  linea0l->Draw("same");
  TLine *  linea27l = new TLine( high_0[2][4],0.,high_0[2][4],0.1);
  linea27l->SetLineColor(kGray);
  linea27l->SetLineStyle(2);
  linea27l->Draw("same");
  
  TArrow * ar0  = new TArrow( high_0[0][17],0.08, high_0[0][17]-20,0.08,0.02,"|>"); //,"<|");
  ar0->SetLineColor(kGray);
  ar0->SetFillColor(kGray);
  ar0->SetAngle(30);

  ar0->Draw();

  TArrow * ar27  = new TArrow( high_0[2][17],0.08, high_0[2][17]-20,0.08,0.02,"|>"); //,"<|");
  ar27->SetLineColor(kGray);
  ar27->SetFillColor(0);
  ar27->SetAngle(30);

  ar27->Draw();

  
  
    TLegend* leg = new TLegend(0.35,0.7,0.85,0.88);
    leg->SetLineColor(0);
    leg->SetFillStyle(0);
    //    leg->SetTextSize(0.05);
    resolutionPlotA[0] = new TGraphErrors(perc,Percentage,res_0[0],chargeerr[0],reserr_0[0]);
    resolutionPlotA[1] = new TGraphErrors(perc,Percentage,res_0[2],chargeerr[2],reserr_0[2]);
    //resolutionPlot[0] = new TGraphErrors(perc,charge_0[i],res_0[i],chargeerr[i],reserr_0[i]);
    //resolutionPlot[1] = new TGraphErrors(perc,charge_1[i],res_1[i],chargeerr[i],reserr_1[i]);

    leg->AddEntry(fLandau[0][0],irr[0],"");
    leg->AddEntry(resolutionPlotA[0],AnglesNice[0],"lp");
    leg->AddEntry(resolutionPlotA[1],AnglesNice[2],"lp");

    leg->Draw();

    /*
    for( int j = 0 ; j < perc ; j++){
      cout << Percentage[j] << " charge " << charge_0[i][j] << " res " << res_0[i][j] << endl;
      cout << Percentage[j] << " charge " << charge_1[i][j] << " res " << res_1[i][j] << endl;
      
    }
    */
    
    //    pad2->Draw();
    pad1->cd();

    gPad->SetTicks(1,1);
    
    resolutionPlotA[0]->SetTitle(" ");
    resolutionPlotA[0]->GetYaxis()->SetTitle("#sigma_{x} [#mum]");
    resolutionPlotA[0]->GetXaxis()->SetTitle("Cluster charge [%]");
    resolutionPlotA[0]->SetMarkerSize(1.);

    resolutionPlotA[0]->GetXaxis()->SetNdivisions(15, kTRUE);
    
    resolutionPlotA[0]->GetXaxis()->SetTitleFont(43);
    resolutionPlotA[0]->GetXaxis()->SetTitleOffset(2.5);
    resolutionPlotA[0]->GetXaxis()->SetTitleSize(15); // labels will be 14 pixels
    
    resolutionPlotA[0]->GetYaxis()->SetTitleFont(43);
    resolutionPlotA[0]->GetYaxis()->SetTitleSize(15); // labels will be 14 pixels
    resolutionPlotA[0]->GetYaxis()->SetTitleOffset(1.8);
    
    resolutionPlotA[0]->GetXaxis()->SetLabelFont(43);
    resolutionPlotA[0]->GetXaxis()->SetLabelSize(15); // labels will be 14 pixels
    resolutionPlotA[0]->GetYaxis()->SetLabelFont(43);
    resolutionPlotA[0]->GetYaxis()->SetLabelSize(15); // labels will be 14 pixels
    
    //    resolutionPlotA[0]->GetYaxis()->SetTitleSize(0.07);
    //resolutionPlotA[0]->GetYaxis()->SetTitleOffset(0.5);
    //    resolutionPlotA[0]->GetYaxis()->SetNdivisions(5);
    //resolutionPlotA[0]->GetYaxis()->SetLabelSize(0.07);

    //resolutionPlotA[0]->GetXaxis()->SetTitleSize(0.07);
    //resolutionPlotA[0]->GetXaxis()->SetLabelSize(0.07);
    
    resolutionPlotA[0]->SetMarkerColor(kBlack);
    resolutionPlotA[0]->SetLineColor(kBlack);
    resolutionPlotA[0]->SetMarkerStyle(20);
    resolutionPlotA[0]->GetXaxis()->SetLimits(0.,1.05);
    //    resolutionPlotA[0]->GetXaxis()->SetRangeUser(0.,1.05);
    resolutionPlotA[0]->GetYaxis()->SetRangeUser(0.,15.);
    resolutionPlotA[0]->Draw("AEP");
    
    resolutionPlotA[1]->SetMarkerSize(1.);
    resolutionPlotA[1]->SetMarkerColor(kBlack);
    resolutionPlotA[1]->SetLineColor(kBlack);
    resolutionPlotA[1]->SetLineStyle(2);
    resolutionPlotA[1]->SetMarkerStyle(24);
    resolutionPlotA[1]->Draw("EPsame");

    TBox * box = new TBox(0.04,7,0.22,13.5);
    box->SetFillColor(kWhite);
    box->Draw();

    cout << " plotted "<< endl;
    
      TLine *  linea9;
    linea9    = new TLine( 0.9,0.,0.9,15.);
    linea9->SetLineColor(kGray);
    linea9->SetLineStyle(3);
    linea9->Draw("same");
      
    TArrow * ar  = new TArrow(0.9,4,0.9-0.05,4,0.02,"|>"); //,"<|");
    ar->SetLineColor(kGray);
    ar->SetFillColor(kGray);
    ar->SetAngle(30);
    
    ar->Draw();
    /*  
    TLine *  linea5;
    linea5    = new TLine( 0.25,0.,0.25,15.);
    linea5->SetLineColor(kGray);
    linea5->SetLineStyle(3);
    linea5->Draw("same");
    */
    
    //TLegend* legFDB2 = new TLegend(0.35,0.75,0.85,0.89);
    //legFDB2->SetLineColor(0);
    //legFDB2->SetFillStyle(0);
    //legFDB2->SetTextSize(0.035);
    //legFDB2->AddEntry(resolutionPlotA[0],AnglesNice[i],"");
    //for(int i =0; i < irradiations; i++)
    //legFDB2->AddEntry(resolutionPlotA[i],irr[i],"lp");

    //legFDB2->Draw();
    

    TString  outname = outputDir+"ResolutionSummaryPaper_angleScans25NONIRR_charge_0andshallow";
    cnonirr->SaveAs(outname+".eps");
    cnonirr->SaveAs(outname+".png");
    cnonirr->SaveAs(outname+".pdf");
    cnonirr->SaveAs(outname+".root");
    cnonirr->SaveAs(outname+".C");
    cout << "saved" << endl;
    cout << "done"<< endl;




    /// second attempt with ratio plot


    TCanvas *cnonirr2 = new TCanvas("cnonirr", "resolution", 600, 600);
    
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
    
    TPad *pad2a = new TPad("pad2","",0,0.,1,0.4);
    pad2a->SetBottomMargin(0.18);
    pad2a->SetTopMargin(0.05);
    pad2a->Draw();

    TPad *pad1a = new TPad("pad1","",0,0.7,1,1);
    pad1a->Draw();
    TPad *pad1b = new TPad("pad1","",0,0.4,1,0.7);
    pad1b->Draw();

    pad1a->SetTopMargin(0.05);
    pad1b->SetBottomMargin(0.02);
    pad1b->SetTopMargin(0.02);
    pad1b->SetBottomMargin(0.2);
  
    pad2a->cd();
    gPad->SetTicks(1,1);
    fLandau[0][0]->GetXaxis()->SetRangeUser(0.,500.);
    fLandau[0][0]->Draw("hist");
    fLandau[0][2]->Draw("histsame");

  
  
    linea0->Draw("same");
    linea0l->Draw("same");
    linea27->Draw("same");
    linea27l->Draw("same");


    
    leg->Draw();

    pad1a->cd();

    gPad->SetTicks(1,1);
    
    resolutionPlotA[0]->GetXaxis()->SetTitleFont(0);
    resolutionPlotA[0]->GetXaxis()->SetTitleOffset(0);
    resolutionPlotA[0]->GetXaxis()->SetTitleSize(0); // labels will be 14 pixels
    resolutionPlotA[0]->GetXaxis()->SetLabelFont(0);
    resolutionPlotA[0]->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels
    resolutionPlotA[0]->GetXaxis()->SetLimits(0.,1.05);
    resolutionPlotA[0]->GetXaxis()->SetNdivisions(20, kTRUE);
    resolutionPlotA[0]->Draw("AEP");
    resolutionPlotA[1]->Draw("EPsame");



    cout << " plotted "<< endl;
    
    linea9->Draw("same");
      
        TLine *  linea5;
    linea5    = new TLine( 0.25,0.,0.25,15.);
    linea5->SetLineColor(kGray);
    linea5->SetLineStyle(3);
    linea5->Draw("same");
    
    // ratio
    pad1b->cd();
    gPad->SetTicks(1,1);
    TGraph*ratio = new TGraph(perc);
    ratio->SetTitle("");
    for (int i=0; i<perc; i++) ratio->SetPoint(i, Percentage[i], res_0[2][i]/res_0[0][i]);
    ratio->SetTitle(" ");
    ratio->SetMarkerSize(1.);


    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleOffset(4);
    ratio->GetXaxis()->SetTitleSize(15); // labels will be 14 pixels

    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(15); // labels will be 14 pixels
    ratio->GetYaxis()->SetTitleOffset(1.8);

    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(15); // labels will be 14 pixels
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(15); // labels will be 14 pixels

    //    ratio->GetYaxis()->SetTitleSize(0.07);
    //ratio->GetYaxis()->SetTitleOffset(0.5);
    //    ratio->GetYaxis()->SetNdivisions(5);
    //ratio->GetYaxis()->SetLabelSize(0.07);

    //ratio->GetXaxis()->SetTitleSize(0.07);
    //ratio->GetXaxis()->SetLabelSize(0.07);

    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerStyle(20);
    ratio->GetXaxis()->SetLimits(0.,1.05);
    ratio->GetXaxis()->SetNdivisions(20, kTRUE);
    ratio->GetYaxis()->SetRangeUser(1.,1.25);
        
    ratio->GetYaxis()->SetTitle("#sigma_{x}^{27.5 deg}/#sigma_{x}^{0 deg}");
    ratio->GetXaxis()->SetTitle("Cluster charge [%]");
    
    ratio->Draw("AP");
    linea9    = new TLine( 0.9,1.,0.9,1.25);
    linea9->SetLineColor(kGray);
    linea9->SetLineStyle(3);
    linea9->Draw("same");
    
    
    // TLine *
      linea5    = new TLine( 0.25,1.,0.25,1.25);
    linea5->SetLineColor(kGray);
    linea5->SetLineStyle(3);
    linea5->Draw("same");
    

    outname = outputDir+"ResolutionSummaryPaper_angleScans25NONIRRratio_charge_0andshallow";
    cnonirr2->SaveAs(outname+".eps");
    cnonirr2->SaveAs(outname+".png");
    cnonirr2->SaveAs(outname+".pdf");
    cnonirr2->SaveAs(outname+".root");
    cnonirr2->SaveAs(outname+".C");
    cout << "saved" << endl;
    cout << "done"<< endl;


}//resVScharge_forPaper

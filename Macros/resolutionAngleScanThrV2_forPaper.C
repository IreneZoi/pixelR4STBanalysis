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
#define anglesNONirrD 5
#define anglesPirr2 28
#define anglesNirr4 9
#define anglesNirr4D 4
#define irradiations 3
#define dphcuts 9


void TDR();
void TDR2(TCanvas * c_all, int period=0, int pos= 11);

void resolutionAngleScanThrV2_forPaper(TString thr="500",TString name = "preliminary_ThrScan", TString func = "RMSself", bool unfolding = false){
  
  TString inputDir="/home/zoiirene/Output/TextFiles/";
  TString outputDir="/home/zoiirene/Output/Plots/";

  int measurements[irradiations+1];
  measurements[0] = anglesNONirr;
  measurements[1] = anglesPirr2;
  measurements[2] = anglesNirr4;
  measurements[3] = 1; 

 
  double i_dphcutPerc[dphcuts]={ 4 , 6 , 9 , 12 , 15 , 18 , 20 , 25 , 30 };
  TString ss_dphcutPerc[dphcuts]={"4","6","9","12","15","18","20","25","30"};
  TString ss_dphcut_0[dphcuts]={"8","12","18","23","29","35", "38", "48", "57"};
  TString ss_dphcut_1[dphcuts]={"5", "8", "11", "15", "19", "23", "25", "31", "38"};
  TString ss_dphcut_2[dphcuts]={"6", "10", "15", "19", "24", "29", "32", "40", "48"};
 
  /* 
  TString ss_dphcutPerc[dphcuts]={"4","6","9","12","20","30"};

  TString ss_dphcut_0[dphcuts]={"8","12","18","23", "38", "57"};
  TString ss_dphcut_1[dphcuts]={"5", "8", "11", "15", "25", "38"};
  TString ss_dphcut_2[dphcuts]={"6", "10", "15", "19", "32", "48"};
  */

  TString filenames_pt1[irradiations+1];
  filenames_pt1[0] = "Ascan_resTree_148_dphcutB"; 
  filenames_pt1[1] = "Ascan_resTreeCorr_120i_dphcutB";
  filenames_pt1[2] = "Ascan_194i_dphcutB"; 
  filenames_pt1[3] = "AScanTreeCorr_194i_5p6_RMSself_800_bestAngle_res_dphcut"; // 6_thrScan_A12C13.txt
  TString filenames_pt2[irradiations+1];
  filenames_pt2[0] = "_RMSselfiso600_beamdiv_A13C14.txt"; //"_RMSselfthrScan_thrScan_A13C14.txt"; 
  filenames_pt2[1] = "_RMSself_beamdiv_A12C15_RMSselfiso600_beamdiv_A13C14.txt"; //"_RMSself_RMSself6sig_closest_A13C14_bestnonirr.txt"; 
  filenames_pt2[2] = "_RMSself_800_resTreeCorr_beamdiv_A12C13.txt"; // thrScan_A12C13.txt"; //5.2 GeV
  filenames_pt2[3] = "_beamdiv_A12C13.txt";

  TString filenames2_pt1[irradiations+1];
  filenames2_pt1[0] = "Ascan_clsizeB_148_dphcutB";
  filenames2_pt1[1] = "Ascan_clsizeB_120i_dphcutB";
  filenames2_pt1[2] = "Ascan_clsizeB_194i_dphcutB";
  filenames2_pt1[3] = "ThrScan_clsizeB_194i_800_dphcuts"; 
  
  TString filenames2_pt2[irradiations+1];
  filenames2_pt2[0] = "_RMSself_beamdiv_A13C14_ABC90.txt"; //"_RMSself_thrScan_A13C14.txt";
  filenames2_pt2[1] = "_RMSself_beamdiv_A12C15_ABC90.txt";//"__RMSself_thrScan_A12C15.txt";
  filenames2_pt2[2] = "_800_beamdiv_A12C13_ABC90.txt"; //"_800_thrScan_A12C13.txt";
  filenames2_pt2[3] = "_bestAngle_5p6_beamdiv_A12C13.txt"; //thrScan_A12C13.txt";
  
  TString irr[irradiations+1];
  irr[0] = " Non-irradiated, 120 V, 5.6 GeV"; // "no irr, 5.6 GeV";
  //  irr[1] = "proton irr at #phi_{eq}=2#times10^{15} cm^{-2}, 5.6 GeV";
  irr[1] = "#phi_{eq} = 2.1 #times 10^{15} cm^{-2}, proton, 800 V, 5.6 GeV";
  irr[2] = "#phi_{eq} = 3.6 #times 10^{15} cm^{-2}, neutron, 800 V, 5.2 GeV";
  irr[3] = "#phi_{eq} = 3.6 #times 10^{15} cm^{-2}, neutron, 800 V, 5.2 GeV";


  TString irrshort[irradiations+1];
  irrshort[0] = " Non-irradiated, 120 V"; // "no irr, 5.6 GeV";
  //  irr[1] = "proton irr at #phi_{eq}=2#times10^{15} cm^{-2}, 5.6 GeV";
  irrshort[1] = "#phi_{eq} = 2.1 #times 10^{15} cm^{-2}, proton, 800 V";
  irrshort[2] = "#phi_{eq} = 3.6 #times 10^{15} cm^{-2}, neutron, 800 V";
  irrshort[3] = "#phi_{eq} = 3.6 #times 10^{15} cm^{-2}, neutron, 800 V";
  
  
  
  TString ss_irr[irradiations] = {"nonirr","2E15prot","4E15neut"};
  TString SS_irr[irradiations] = {"NonIrr","ProtIrr2E15","NeutrIrr3p6E15"};
  double angle_0[anglesNONirr];
  double angle_1[anglesPirr2];
  double angle_2[anglesNirr4];

  double angleerr_0[dphcuts][anglesNONirr];
  double angleerr_1[dphcuts][anglesPirr2];
  double angleerr_2[dphcuts][anglesNirr4];

  double res_0[dphcuts][anglesNONirr];
  double res_1[dphcuts][anglesPirr2];
  double res_2[dphcuts][anglesNirr4];
  double res_3[dphcuts];

  double reserr_0[dphcuts][anglesNONirr];
  double reserr_1[dphcuts][anglesPirr2];
  double reserr_2[dphcuts][anglesNirr4];
  double reserr_3[dphcuts];


  
  double angle;
  double res,reserror,angleerr;
  int lines[irradiations+1] = {8,8,7,7}; //lines to be skipped in input files before reading the measurements
  
  for(int i =0; i< irradiations+1; i++)    {
    cout << irr[i] << endl;
    for(int l = 0; l< dphcuts ; l++){
      cout << "dphcut % "<< ss_dphcutPerc[l] << endl;
      TString  ss_dphcut=ss_dphcut_0[l];
      if(i==1) ss_dphcut=ss_dphcut_1[l];
      if(i==2 || i==3) {
	cout << "dphcut ADC "<< ss_dphcut_2[l] << endl;
	ss_dphcut=ss_dphcut_2[l];
      }
      cout << "dphcut ADC "<< ss_dphcut << endl;
       
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
	    if(l==0)	    angle_0[j] = angle;
	    angleerr_0[l][j] = 0;
	    res_0[l][j] = res;
	    reserr_0[l][j] = reserror;
	    if(print)	 cout  << angle_0[j] << " " << res_0[l][j] << " " << reserr_0[l][j] << endl;
	  }
	  if(i==1){
	    if(l==0) angle_1[j] = angle;
	    angleerr_1[l][j] = 0;
	    res_1[l][j] = res;
	    reserr_1[l][j] = reserror;
	    if(print)	 cout  << angle_1[j] << " " << res_1[l][j] << " " << reserr_1[l][j] << endl;
	  }
	  if(i==2){
	    if(l==0)    angle_2[j] = angle;
	    angleerr_2[l][j] = 0;
	    res_2[l][j] = res;
	    reserr_2[l][j] = reserror;
	    if(print)	 cout  << angle_2[j] << " " << res_2[l][j] << " " << reserr_2[l][j] << endl;
	  }
          if(i==3){
	    res_3[l] = res;
	    reserr_3[l] = reserror;
	    if(print)    cout  << angle << " " << res_3[l] << " " << reserr_3[l] << endl;
	  }
	  
	}//measurements
      }//file scanning
    }//cuts
  }//irradiations



  //clsize file for IEEE
  double clsize_0[anglesNONirr][dphcuts];
  double clsize_1[anglesPirr2][dphcuts];
  double clsizeerr_0[anglesNONirr][dphcuts];
  double clsizeerr_1[anglesPirr2][dphcuts];
  double clsizeerr_2[anglesNirr4][dphcuts];
  double clsize_2[anglesNirr4][dphcuts];
  double clsize_3[dphcuts]; //mean of nrowB of straightTracksY_isoAandCandB_straightTracksX (same as after cut on charge), drei-r3839_irene_dphcutB15_closest_A12C13.root neutron irr, 5.6 GeV
  double clsizeerr_3[dphcuts]; //mean error

  double clsize,clsizeerror;
  cout << " %%%%%%%%% now reading cluster size " << endl;
  int lines2[irradiations+1] = {8,8,8,8};
  for(int i =0; i< irradiations+1; ++i)    {
    TString pt1= filenames2_pt1[i];
    cout << irr[i] << " " << filenames2_pt1[i] << endl;
    for(int l = 0; l< dphcuts ; l++){
      cout << "dphcut % "<< ss_dphcutPerc[l] << endl;
      TString  ss_dphcut=ss_dphcut_0[l];
      if(i==1) ss_dphcut=ss_dphcut_1[l];
      if(i==2 || i==3) {
	ss_dphcut=ss_dphcut_2[l];
	cout << l << " dphcut ADC "<< ss_dphcut_2[l] << endl;
      }
      cout << "dphcut ADC "<< ss_dphcut << endl;

      cout << " inputDir "<< inputDir <<  " " << filenames2_pt1[i]<<  endl;
      cout << inputDir+filenames2_pt1[i]+ss_dphcut+filenames2_pt2[i] << endl;
      TString filename      = inputDir+pt1+ss_dphcut+filenames2_pt2[i];
      
      cout << filename << endl;
      ifstream stream(filename);
      
      std::string line;
      if(!stream.is_open())	{
	cout << " File " << filename << " not opened" << endl;
      }
      else{
	for(int k =0; k<lines2[i];k++) {
	  std::getline(stream,line);
	  if(print) cout << k << " line " << line << endl;
	}
        cout << " getting #  measurements "<< measurements[i] << endl;
	for(int j= 0; j < measurements[i]; j++){
	  stream  >> angle >> clsize >> clsizeerror ;
	  if(print)	 cout << "line " << j+lines[i]-1 << endl;
	  if(i==0){
	    clsize_0[j][l] = clsize;
	    clsizeerr_0[j][l] = clsizeerror;
	    if(print)	 cout  << angle_0[j] << " " << clsize_0[j][l] << " " << clsizeerr_0[j][l] << endl;
	  }
	  if(i==1){
	    clsize_1[j][l] = clsize;
	    clsizeerr_1[j][l] = clsizeerror;
	    if(print)	 cout  << angle_1[j] << " " << clsize_1[j][l] << " " << clsizeerr_1[j][l] << endl;
	  }
	  if(i==2){
	    clsize_2[j][l] = clsize;
	    clsizeerr_2[j][l] = clsizeerror;
	    if(print)	 cout  << angle_2[j] << " " << clsize_2[j][l] << " " << clsizeerr_2[j][l] << endl;
	  }
          if(i==3){
	    clsize_3[l] = clsize;
	    clsizeerr_3[l] = clsizeerror;
	    if(print)    cout  << angle << " " << clsize_3[l] << " " << clsizeerr_3[l] << endl;
	  }
	  cout << " cut "<< ss_dphcut <<  " " << filenames2_pt1[i]<<  endl;
	}//measurements
      }//file scanning
    }//cuts
  }//irradiations







  measurements[0] = anglesNONirrD;
  measurements[1] = anglesNONirrD;
  measurements[2] = anglesNirr4D;
  //  int index[anglesNONirrD]={4,8,10,11,12,13,14,17,20,23}; //,26};
  int index[anglesNONirrD]={4,11,13,17,23}; //,26};
  int index2[anglesNirr4D]={1,2,5,7}; //,26};
  double resd_0[anglesNONirrD][dphcuts];
  double err_0[anglesNONirrD][dphcuts];
  double reserrord_0[anglesNONirrD][dphcuts];
  double angled_0[anglesNONirrD]={angle_0[4],angle_0[11],angle_0[13],angle_0[17],angle_0[23]}; //,angle_0[26]};
  double resd_1[anglesNONirrD][dphcuts];
  double reserrord_1[anglesNONirrD][dphcuts];
  double angled_2[anglesNirr4D];
  double resd_2[anglesNirr4D][dphcuts];
  double reserrord_2[anglesNirr4D][dphcuts];


  for(int l = 0; l<dphcuts; l++){
    cout <<"dph "<< l << " " <<res_3[l]<<endl;
    for(int i = 0 ; i <anglesNONirrD; i++){
      cout <<" non irr " <<  angle_0[index[i]] << " "<< res_0[l][index[i]] << endl;
      resd_0[i][l]=res_0[l][index[i]]; //{res_0[l][4],res_0[l][8],res_0[l][10],res_0[l][11],res_0[l][12],res_0[l][13],res_0[l][14],res_0[l][17],res_0[l][20],res_0[l][23],res_0[l][26]};
      reserrord_0[i][l]=reserr_0[l][index[i]]; //{res_0[l][4],res_0[l][8],res_0[l][10],res_0[l][11],res_0[l][12],res_0[l][13],res_0[l][14],res_0[l][17],res_0[l][20],res_0[l][23],res_0[l][26]};
      err_0[i][l] = 0;
      cout << "err" << endl;
      cout << "prot "<< endl;
      resd_1[i][l]=res_1[l][index[i]]; //{res_0[l][4],res_0[l][8],res_0[l][10],res_0[l][11],res_0[l][12],res_0[l][13],res_0[l][14],res_0[l][17],res_0[l][20],res_0[l][23],res_0[l][26]};
      reserrord_1[i][l]=reserr_1[l][index[i]]; //{res_0[l][4],res_0[l][8],res_0[l][10],res_0[l][11],res_0[l][12],res_0[l][13],res_0[l][14],res_0[l][17],res_0[l][20],res_0[l][23],res_0[l][26]};
      cout << reserrord_1[i][l] << endl;
      if(i<anglesNirr4D){
	cout << "neut"<< endl;
	resd_2[i][l]=res_2[l][index2[i]]; //{res_0[l][4],res_0[l][8],res_0[l][10],res_0[l][11],res_0[l][12],res_0[l][13],res_0[l][14],res_0[l][17],res_0[l][20],res_0[l][23],res_0[l][26]};
	reserrord_2[i][l]=reserr_2[l][index2[i]]; //{res_0[l][4],res_0[l][8],res_0[l][10],res_0[l][11],res_0[l][12],res_0[l][13],res_0[l][14],res_0[l][17],res_0[l][20],res_0[l][23],res_0[l][26]};
	angled_2[i]=angle_2[index2[i]];
      }
      cout << resd_0[i][l] << endl;
      
  }
}


  
	        
  
  //  int start[irradiations] = {5,5,0};
  /////// plots!
  //  int symbols_res[dphcuts]={20,21,22,23,29,33,34,39,43};
  //  int colors_0[anglesNONirrD]={625,632,634,626,630};//,627};//,635,631,636};
  int colors_0[anglesNONirrD]={921,1,802,794,803};
  //int colors_1[dphcuts]={14,13,12,1,620,615};//,619,618,611};
  //int colors_2[dphcuts]={591,593,600,602,598,603};//,599,604,1};
  
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

    gStyle->SetTextFont(43);
    gStyle->SetTextSize(10);
    gStyle->SetLegendFont(43);
    gStyle->SetLegendTextSize(20);
    
    TLegend* legFDB2 = new TLegend(0.15,0.65,0.8,0.8);
    legFDB2->SetNColumns(4);
    legFDB2->SetLineColor(0);
    //legFDB2->SetTextSize(0.03);
    TLegend* leg = new TLegend(0.05,0.8,0.45,0.85);
    leg->SetLineColor(0);
    //leg->SetTextSize(0.035);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    //TPaveText *pt = new TPaveText(.15,.8,.65,.85,"NDC");
    //
    
    //pt->Draw();
    TGraphErrors* resolutionPlot[measurements[i]];
    for (int l =0; l < measurements[i]; l++){

      resolutionPlot[l] = new TGraphErrors(dphcuts,i_dphcutPerc,resd_0[l],err_0[l],reserrord_0[l]);
      double ANGLE = angled_0[l];
      if(i==1){
	resolutionPlot[l] = new TGraphErrors(dphcuts,i_dphcutPerc,resd_1[l],err_0[l],reserrord_1[l]);
	ANGLE = angled_0[l];
      }
      if(i==2){
	resolutionPlot[l] = new TGraphErrors(dphcuts,i_dphcutPerc,resd_2[l],err_0[l],reserrord_2[l]);
	ANGLE = angled_2[l];
      }   
      resolutionPlot[l]->SetTitle(" ");

      resolutionPlot[l]->SetMarkerSize(1.);
      resolutionPlot[l]->SetMarkerColor(colors_0[l]);
      resolutionPlot[l]->SetLineColor(colors_0[l]);
      resolutionPlot[l]->SetMarkerStyle(l+20);
      if(l==1) resolutionPlot[l]->SetMarkerStyle(20);
      if(l==0){
	leg->AddEntry(resolutionPlot[l],irr[i],"");
	resolutionPlot[l]->SetMarkerStyle(21);
	resolutionPlot[0]->GetYaxis()->SetTitle("#sigma_{x} [#mum]");
	resolutionPlot[0]->GetXaxis()->SetTitle("Threshold [%]");
	resolutionPlot[0]->GetYaxis()->SetTitleFont(43);
	resolutionPlot[0]->GetYaxis()->SetTitleSize(20); // labels will be 14 pixels

	resolutionPlot[0]->GetYaxis()->SetLabelFont(43);
	resolutionPlot[0]->GetYaxis()->SetLabelSize(20); // labels will be 14 pixels

	resolutionPlot[0]->GetXaxis()->SetLabelFont(43);
	resolutionPlot[0]->GetXaxis()->SetLabelSize(20); // labels will be 14 pixels
	
	resolutionPlot[0]->GetXaxis()->SetTitleFont(43);
	resolutionPlot[0]->GetXaxis()->SetTitleSize(20); // labels will be 14 pixels


	resolutionPlot[0]->GetXaxis()->SetLimits(0.,35.);
	resolutionPlot[0]->GetYaxis()->SetRangeUser(0.,12.);
	resolutionPlot[0]->Draw("AEP");
      }	  
      if(l !=3)
	resolutionPlot[l]->Draw("EPsame");
      

      TString ss_angle;
      ss_angle.Form("%.2g",ANGLE);
      if(l !=3) legFDB2->AddEntry(resolutionPlot[l],ss_angle+" deg","lp");
    }
    legFDB2->Draw();
    leg->Draw();
    TString  outname = outputDir+"ResolutionSummaryPaper_angleScans25_"+name+"_"+ss_irr[i];
    cFDB2[i]->SaveAs(outname+".eps");
    cFDB2[i]->SaveAs(outname+".png");
    cFDB2[i]->SaveAs(outname+".pdf");
    cFDB2[i]->SaveAs(outname+".root");
    cFDB2[i]->SaveAs(outname+".C");
  }

  ///////////////////////////////
  //now compare irradiations at best angle

  int colors_irr[irradiations]={1,417,616};
  int marker_irr[irradiations]={20,21,32};
  int bestThr[irradiations]={6,12,9};
  
  TCanvas *c  = new TCanvas("c", "FDB resolution", 600, 600);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  c->SetLeftMargin(1.);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.03);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(24);
  
  TLegend* legAng = new TLegend(0.2,0.86,0.45,0.9);
  //legFDB2->SetNColumns(5);
  legAng->SetLineColor(0);
  //legAng->SetTextSize(0.03);
  TLegend* legIrr = new TLegend(0.15,0.63,0.8,0.85);
  legIrr->SetLineColor(0);
		  //legIrr->SetTextSize(0.035);
  legIrr->SetFillStyle(0);
  legIrr->SetBorderSize(0);
    //TPaveText *pt = new TPaveText(.15,.8,.65,.85,"NDC");
    //
    
    //pt->Draw();

   TGraphErrors* resolutionPlotIrr[irradiations];
   for (int l =0; l < irradiations; l++){
    
    if(l==0) resolutionPlotIrr[0] = new TGraphErrors(dphcuts,i_dphcutPerc,resd_0[1],err_0[1],reserrord_0[1]);
    if(l==1) resolutionPlotIrr[1] = new TGraphErrors(dphcuts,i_dphcutPerc,resd_1[1],err_0[1],reserrord_1[1]);
    if(l==2) resolutionPlotIrr[2] = new TGraphErrors(dphcuts,i_dphcutPerc,res_3,err_0[2],reserr_3);

      resolutionPlotIrr[l]->SetTitle(" ");
      legIrr->AddEntry(resolutionPlotIrr[l],irrshort[l],"lp");
      resolutionPlotIrr[l]->SetMarkerSize(1.);
      resolutionPlotIrr[l]->SetMarkerColor(colors_irr[l]);
      resolutionPlotIrr[l]->SetLineColor(colors_irr[l]);
      resolutionPlotIrr[l]->SetMarkerStyle(marker_irr[l]);
      if(l==0){
	
	legAng->AddEntry(resolutionPlotIrr[l],"Optimal angle, 5.6 GeV",0);
	resolutionPlotIrr[0]->GetYaxis()->SetTitle("#sigma_{x} [#mum]");
	resolutionPlotIrr[0]->GetXaxis()->SetTitle("Threshold [%]");

        resolutionPlotIrr[0]->GetYaxis()->SetTitleFont(43);
        resolutionPlotIrr[0]->GetYaxis()->SetTitleSize(30); // labels will be 14 pixels

	resolutionPlotIrr[0]->GetYaxis()->SetLabelFont(43);
	resolutionPlotIrr[0]->GetYaxis()->SetLabelSize(30); // labels will be 14 pixels

	resolutionPlotIrr[0]->GetXaxis()->SetLabelFont(43);
	resolutionPlotIrr[0]->GetXaxis()->SetLabelSize(24); // labels will be 14 pixels

	resolutionPlotIrr[0]->GetXaxis()->SetTitleFont(43);
	resolutionPlotIrr[0]->GetXaxis()->SetTitleSize(24); // labels will be 14 pixels
	/*
	resolutionPlotIrr[0]->GetXaxis()->SetTitleSize(0.05);
	resolutionPlotIrr[0]->GetYaxis()->SetTitleSize(0.05);
	resolutionPlotIrr[0]->GetXaxis()->SetTitleOffset(0.9);
	resolutionPlotIrr[0]->GetYaxis()->SetTitleOffset(1.);
	resolutionPlotIrr[0]->GetXaxis()->SetLabelSize(0.05);
	resolutionPlotIrr[0]->GetYaxis()->SetLabelSize(0.05);
	*/
	resolutionPlotIrr[0]->GetXaxis()->SetLimits(0.,35.);
	resolutionPlotIrr[0]->GetYaxis()->SetRangeUser(0.,12.);
	resolutionPlotIrr[0]->Draw("AEP");
      }	  
      else
	resolutionPlotIrr[l]->Draw("EPsame");
      


    }
    legAng->Draw();
    legIrr->Draw();
    TString  outname = outputDir+"ResolutionSummaryPaper_angleScans25_"+name+"_allirr";
    c->SaveAs(outname+".eps");
    c->SaveAs(outname+".png");
    c->SaveAs(outname+".pdf");
    c->SaveAs(outname+".root");
    c->SaveAs(outname+".C");




  ///////////////////////////////
  //now compare irradiations and cluster size at best angle


  TCanvas *cc  = new TCanvas("cc", "FDB resolution", 600, 600);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  cc->SetLeftMargin(0.16);
  cc->SetTopMargin(0.1);
  cc->SetBottomMargin(0.03);
  cc->SetRightMargin(0.03);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetTextFont(43);
  gStyle->SetTextSize(10);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(24);

  cout << " initialized new canvas" << endl;
  TPad *pad2 = new TPad("pad2","",0,0.,1,0.45);
  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.25);
  pad2->SetRightMargin(0.03);
  pad2->SetLeftMargin(0.12);

  pad2->Draw();

  //pad2->SetFillStyle(4000); //will be transparent
  //#pad2->SetFrameFillStyle(0);
  TPad *pad1 = new TPad("pad1","",0,0.45,1,1);
  pad1->Draw();
  pad1->SetBottomMargin(0.03);
  pad1->SetTopMargin(0.04);
  pad1->SetRightMargin(0.03);
  pad1->SetLeftMargin(0.12);

  pad1->cd();
  gPad->SetTicks(1,1);
  
  resolutionPlotIrr[0]->GetXaxis()->SetTitleSize(0); // labels will be 14 pixels
  resolutionPlotIrr[0]->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels
  resolutionPlotIrr[0]->Draw("AEP");

  resolutionPlotIrr[1]->Draw("EPsame");
  resolutionPlotIrr[2]->Draw("EPsame");


  //// Write Txt files for Robert
  for( int k = 0; k< irradiations; k++){
    ofstream myres;
    myres.open ("/home/zoiirene/Output/TextFiles/Fig16_Resolution_"+SS_irr[k]+"_8p8deg.txt");
    myres << "Offline Threshold [%]  Resolution [um] \n";
    for(int i=0; i < dphcuts; i++){
      if( k == 0)    myres << i_dphcutPerc[i] << " " << resd_0[1][i] << "\n";
      if( k == 1)    myres << i_dphcutPerc[i] << " " << resd_1[1][i] << "\n";
      if( k == 2)    myres << i_dphcutPerc[i] << " " << res_3[i] << "\n";
  }
  myres.close();
  }

  legAng->Draw();
  legIrr->Draw();



		  pad2->cd();
  //  gROOT->SetStyle("Plain");
  TGraphErrors* clsizePlot[irradiations];
  gPad->SetTicks(1,1);
		  
  clsizePlot[0] = new TGraphErrors(dphcuts,i_dphcutPerc,clsize_0[11],err_0[1],clsizeerr_0[11]);
  clsizePlot[1] = new TGraphErrors(dphcuts,i_dphcutPerc,clsize_1[11],err_0[1],clsizeerr_1[11]);
  clsizePlot[2] = new TGraphErrors(dphcuts,i_dphcutPerc,clsize_3,err_0[1],clsizeerr_3);
  cout << " initialized new tgraps" << endl;



  clsizePlot[0]->GetYaxis()->SetTitle("Mean cluster size");
  clsizePlot[0]->GetXaxis()->SetTitle("Threshold [%]");

  clsizePlot[0]->GetYaxis()->SetNdivisions(5,0,5);
			       
  clsizePlot[0]->GetXaxis()->SetTitleFont(43);
  clsizePlot[0]->GetXaxis()->SetTitleSize(30); // labels will be 14 pixels
  clsizePlot[0]->GetXaxis()->SetTitleOffset(1.9); // labels will be 14 pixels
  clsizePlot[0]->GetXaxis()->SetLabelFont(43);
  clsizePlot[0]->GetXaxis()->SetLabelSize(30); // labels will be 14 pixels

  clsizePlot[0]->GetYaxis()->SetTitleFont(43);
  clsizePlot[0]->GetYaxis()->SetTitleSize(30); // labels will be 14 pixels
  clsizePlot[0]->GetYaxis()->SetLabelFont(43);
  clsizePlot[0]->GetYaxis()->SetLabelSize(30); // labels will be 14 pixels
  clsizePlot[0]->GetYaxis()->SetTitleOffset(0.8); // labels will be 14 pixels
  clsizePlot[0]->GetXaxis()->SetLimits(0.,35.);
  clsizePlot[0]->GetYaxis()->SetRangeUser(0.,5.);

  
  for (int l =0; l < irradiations; l++){


  //// Write Txt files for Robert
    ofstream myres;
    myres.open ("/home/zoiirene/Output/TextFiles/Fig16_ClusterSize_"+SS_irr[l]+"_8p8deg.txt");
    myres << "Offline Threshold [%]  Resolution [um] \n";
    for(int i=0; i < dphcuts; i++){
      if( l == 0)    myres << i_dphcutPerc[i] << " " << clsize_0[11][i] << "\n";
      if( l == 1)    myres << i_dphcutPerc[i] << " " << clsize_1[11][i] << "\n";
      if( l == 2)    myres << i_dphcutPerc[i] << " " << clsize_3[i] << "\n";
    }
    myres.close();


    clsizePlot[l]->SetTitle(" ");
    clsizePlot[l]->SetMarkerSize(1.);
    clsizePlot[l]->SetMarkerColor(colors_irr[l]);
    clsizePlot[l]->SetLineColor(colors_irr[l]);
    clsizePlot[l]->SetMarkerStyle(marker_irr[l]);
    cout << " initialized styles" << endl;
    if(l==0)  clsizePlot[0]->Draw("AEP");
    else
      clsizePlot[l]->Draw("EPsame");

    TLine *  lineb = new TLine( bestThr[l],0.,bestThr[l],5.);
    lineb->SetLineColor(colors_irr[l]);
    lineb->SetLineWidth(2);
    lineb->SetLineStyle(2);
    lineb->Draw("same");
    
  }
  TLine *  linea = new TLine( 0.,2.,35.,2.);
  linea->SetLineColor(kGray);
  linea->SetLineWidth(2);
  linea->SetLineStyle(2);
  linea->Draw("same");
  
  pad1->cd();
   cout << " back to pad 1" << endl;		  
   outname = outputDir+"ResolutionSummaryPaper_angleScans25_"+name+"_clsize_allirr";
   cc->SaveAs(outname+".eps");
   cc->SaveAs(outname+".png");
   cc->SaveAs(outname+".pdf");
   cc->SaveAs(outname+".root");
   cc->SaveAs(outname+".C");











  
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

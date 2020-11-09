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

#define BiasVoltages 4
#define Angles 8
bool print=true;
using namespace std;
#define dphcuts  1
#define comparisons 1
#define thresholds 3


void resolutionVSbias194n_TreeEntries(TString func = "RMSself")
{

  TString detectorA="148";
  TString detectorB="194i";
  TString detectorC="150";
  TString labelA="FDB150P_12_R4S100x25-P4_1, thr 12 ADC";
  TString labelB="FDB150P_12_R4S100x25-P1_1, thr 15 ADC";
  TString labelC="FDB150P_12_R4S100x25-P1_3, thr 12 ADC";
  TString info  ="beam momentum 5.2 GeV";  
  
  
  double BiasVoltage[BiasVoltages];
  BiasVoltage[0]=200;//V  
  BiasVoltage[1]=400;//V 3
  BiasVoltage[2]=600;//V 5
  BiasVoltage[3]=800;//V 7 

  double BiasVoltageError[BiasVoltages];

  double Angle[Angles]={8,10,12,14,16,18,20,22};
  TString ss_Angle[Angles];//={10};//,12};//,14,16,18,20};
  
  Double_t Resolution[Angles][BiasVoltages];
  Double_t Percentage[Angles][BiasVoltages];
  Double_t ResolutionError[Angles][BiasVoltages];
  TString ss_BiasVoltage[BiasVoltages];
  
  TString inputDir="/home/zoiirene/Output/TextFiles/";
  TString outputDir="/home/zoiirene/Output/Plots/";
  TString inputfile;
  TString testname="dycut_A12C13"; //"closest_A12C13";
  TString testnameUnf="RMSself_dycut_A13C14";
  
  TH1F * h_res;

  // Int_t run[BiasVoltages]={3763,3764,3765,3766,3767,3768,3769,3770};
  Int_t run[Angles][BiasVoltages];

  //8
  run[0][3] = 3823;
  run[0][2] = 3825;
  run[0][1] = 3827;
  run[0][0] = 3829;
  //10
  run[1][0] = 3770;
  run[1][1] = 3767;
  run[1][2] = 3765;
  run[1][3] = 3763;
  
  //12
  run[2][0] = 3783;
  run[2][1] = 3780;
  run[2][2] = 3777;
  run[2][3] = 3778;
  
  //14
  run[3][0] = 3791;
  run[3][1] = 3788;
  run[3][2] = 3786;
  run[3][3] = 3784;
  //16
  run[4][0] = 3799;
  run[4][1] = 3796;
  run[4][2] = 3794;
  run[4][3] = 3792;
  
  //18
  run[5][0] = 3813;
  run[5][1] = 3810;
  run[5][2] = 3808;
  run[5][3] = 3806;

  //20
  run[6][0] = 3821;
  run[6][1] = 3818;
  run[6][2] = 3816;
  run[6][3] = 3814;
  
//22
  run[7][3] = 3831;
  run[7][2] = 3833;
  run[7][1] = 3835;
  run[7][0] = 3838;
  

  TString Run[Angles][BiasVoltages];

  ostringstream strs[BiasVoltages];
  ostringstream strsa[Angles];
  TString pitch = "25";  
  if(print) cout << "Iniazializating bias voltages and runs "  << endl;

  MapTH1 res_map[Angles];
  MapTH1 res_map800;
  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  //  int i_dphcut[dphcuts] = {6, 10, 15, 19, 24, 29, 32, 40, 48}; //{15};
  //TString ss_dphcut[dphcuts] = {"6", "10", "15", "19", "24", "29", "32", "40", "48"}; //{"15"};
  int i_dphcut[dphcuts] = {15};
  TString ss_dphcut[dphcuts] = {"15"};

  TString Label[comparisons] = {"straightTracksY_isoAandCandB_straightTracksX"};
  TString Hist = "dx3_clphABC90evR";
  bool dphcut = true;
  
  TH1F * hdx3_clchargeABC90evR[Angles];

  for(int j = 0; j < Angles; j++)
    {
      strsa[j] << Angle[j];
      ss_Angle[j]=strsa[j].str();

      for(int i=0; i<BiasVoltages; i++)
	{ 

	  
	  if(print) cout << "Run: "<< i <<" " << run[j][i] << endl;
	  Run[j][i].Form("%d",run[j][i]);
	  if(print) cout << "Run: " << Run[j][i] << endl;
	  // ss_BiasVoltage[i].Form("%d",BiasVoltage[i]);	  
	  if(j==0)
	    {
	      BiasVoltageError[i]=0; //GeV
	      if(print) cout << "Beam energy " << i<< ": " << BiasVoltage[i] << " GeV" << endl;
	      strs[i] << BiasVoltage[i];
	      ss_BiasVoltage[i]=strs[i].str();
	      if(print) cout << "Beam energy " << i<< ": " << ss_BiasVoltage[i] << " GeV" << endl;
	    }
	}
      
      res_map[j] =  GetCheckHists(&(res_map[j]), BiasVoltages, dphcuts,dphcut, Run[j],ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR[j],true,testname);
      //      GetHists(&(res_map[j]),BiasVoltages, dphcuts, comparisons,dphcut, Run[j], i_dphcut, Label, Hist, h_res);
    }
  



  
  for(int k=0; k<dphcuts;k++){  
  

    if(print) cout << "#############         RMS:     ########################"  << endl;
  

    ofstream myfile[Angles];
    for(int j = 0; j < Angles; j++)
      {
      
	myfile[j].open ("/home/zoiirene/Output/TextFiles/Bscan_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+func+"_"+ss_Angle[j]+"deg_rmsTree_"+testname+".txt");
	myfile[j] << "A " << detectorA << "\n";
	myfile[j] << labelA << "\n";
	myfile[j] << "B " << detectorB<< "\n";
	myfile[j] << labelB << "\n";
	myfile[j] << "C " << detectorC<< "\n";
	myfile[j] << labelC << "\n";
	
	myfile[j] << "Momentum(GeV) RMS(um) Error\n";
      
	for(int i=0; i<BiasVoltages; i++)
	  {
	    auto it2 = res_map[j].find(std::make_pair(Run[j][i]+"_"+ss_dphcut[k]+"_"+pitch,Hist));
	    if(it2  != res_map[j].end())
	      {
		if(print)               cout << " found map " << endl;
		if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
		
		FitTH1(it2->second, &(Resolution[j][i]), &(ResolutionError[j][i]), ss_Angle[j]+"_"+ss_BiasVoltage[i]+"_dphcutB"+ss_dphcut[k], detectorA, detectorB, detectorC,func,&(Percentage[i][j]) );
		if(print) cout << "Beam energy " << i<< ": " << ss_BiasVoltage[i] << " GeV -> Resolution: " << Resolution[j][i] << " and res err: " << ResolutionError[j][i] << endl;
		myfile[j] << ss_BiasVoltage[i] << " " << Resolution[j][i] << " " << ResolutionError[j][i] << "\n";
	      }else{
	      cout << "********************               not found map " << Run[j][i]+"_"+Label[0] << " hist " << Hist << endl;
	    }
	  }
	
	myfile[j].close();
	DrawTGraphWithErrorDouble(BiasVoltages, BiasVoltage, Resolution[j], ResolutionError[j], Hist,Run[j][0], Label[0], "biasScan_dphcutB"+ss_dphcut[k]+"_"+ss_Angle[j]+"_rmsTree", "RMS 95% [#mum]",0.,7.,0.,8., "biasScan_"+ss_Angle[j], "Bias Voltage [V]"  );
      }
  




  ofstream myfileBias[BiasVoltages];
  ofstream myfileBiasPerc[BiasVoltages];
  for(int j = 0; j < BiasVoltages; j++)
    {
      
      myfileBias[j].open ("/home/zoiirene/Output/TextFiles/Ascan_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+func+"_"+ss_BiasVoltage[j]+"_rmsTree_"+testname+".txt");
      myfileBias[j] << "A " << detectorA << "\n";
      myfileBias[j] << labelA << "\n";
      myfileBias[j] << "B " << detectorB<< "\n";
      myfileBias[j] << labelB << "\n";
      myfileBias[j] << "C " << detectorC<< "\n";
      myfileBias[j] << labelC << "\n";

      myfileBias[j] << "Angle(deg) RMS(um) Error\n";


      myfileBiasPerc[j].open ("/home/zoiirene/Output/TextFiles/Ascan_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+func+"_"+ss_BiasVoltage[j]+"_Percentage_"+testname+".txt");
      myfileBiasPerc[j] << "A " << detectorA << "\n";
      myfileBiasPerc[j] << labelA << "\n";
      myfileBiasPerc[j] << "B " << detectorB<< "\n";
      myfileBiasPerc[j] << labelB << "\n";
      myfileBiasPerc[j] << "C " << detectorC<< "\n";
      myfileBiasPerc[j] << labelC << "\n";

      myfileBiasPerc[j] << "Angle(deg) Percentage\n";
      for(int i=0; i<Angles; i++)	{
	      myfileBias[j] << ss_Angle[i] << " " << Resolution[i][j] << " " << ResolutionError[i][j] << "\n";
	      myfileBiasPerc[j] << ss_Angle[i] << " " << Percentage[i][j] << "\n ";
	}
      myfileBias[j].close();
      myfileBiasPerc[j].close();
    }
      

  
  //////////////          RESOLUTION ///////////////////////
  double freshres = 2.70;
  double freshres_err = 0.01;

  for(int j = 0; j < Angles; j++)
    {
      
      for(int i=0; i<BiasVoltages; i++)
	{
	  
	  if(print) cout << "Energy " << i<< ": " << ss_BiasVoltage[i] << " GeV -> Resolution: " << Resolution[j][i] << " and res err: " << ResolutionError[j][i] << endl;
	  ExtractRes(&(Resolution[j][i]),&(ResolutionError[j][i]),true, freshres, freshres_err);
	  if(print) cout << "Energy " << i<< ": " << ss_BiasVoltage[i] << " GeV -> Resolution: " << Resolution[j][i] << " and res err: " << ResolutionError[j][i] << endl;
	  
	  
	}
      DrawTGraphWithErrorDouble(BiasVoltages, BiasVoltage, Resolution[j], ResolutionError[j], Hist,Run[j][0], Label[0], "biasScan_dphcutB"+ss_dphcut[k]+"_"+ss_Angle[j]+"_resTree", "Resolution [#mum]",0.,7.,0.,8., "biasScan_"+ss_Angle[j], "Bias Voltage [V]"  );
    }
  

  ofstream myfileBiasUnf[BiasVoltages];
  for(int j = 0; j < BiasVoltages; j++)
    {

      myfileBiasUnf[j].open ("/home/zoiirene/Output/TextFiles/Ascan_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+func+"_"+ss_BiasVoltage[j]+"_resTree_"+testname+".txt");
      myfileBiasUnf[j] << "A " << detectorA << "\n";
      myfileBiasUnf[j] << labelA << "\n";
      myfileBiasUnf[j] << "B " << detectorB<< "\n";
      myfileBiasUnf[j] << labelB << "\n";
      myfileBiasUnf[j] << "C " << detectorC<< "\n";
      myfileBiasUnf[j] << labelC << "\n";

      myfileBiasUnf[j] << "Angle(deg) RES(um) Error\n";

      for(int i=0; i<Angles; i++)
	{
	  myfileBiasUnf[j] << ss_Angle[i] << " " << Resolution[i][j] << " " << ResolutionError[i][j] << "\n";
	}
      myfileBiasUnf[j].close();
    }



  cout << "######## Angle dependent    resolution unfolding       ###########" << endl;
  TString  filenames = "AscanSelected_resTree_148_dphcutB12_"+testnameUnf+".txt";
  filenames      = inputDir+filenames;
  cout << filenames << endl;
  ifstream stream(filenames);
  double angle;
  double res,reserror;
  double freshresA[Angles];
  double freshresA_err[Angles];

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
  

  for(int j = 0; j < Angles; j++)
    {

      for(int i=0; i<BiasVoltages; i++)
	{

	  if(print) cout << "Energy " << i<< ": " << ss_BiasVoltage[i] << " GeV -> Resolution: " << Resolution[j][i] << " and res err: " << ResolutionError[j][i] << endl;
	  ExtractRes(&(Resolution[j][i]),&(ResolutionError[j][i]),true, freshres, freshres_err);
	  if(print) cout << "Energy " << i<< ": " << ss_BiasVoltage[i] << " GeV -> Resolution: " << Resolution[j][i] << " and res err: " << ResolutionError[j][i] << endl;


	}
      DrawTGraphWithErrorDouble(BiasVoltages, BiasVoltage, Resolution[j], ResolutionError[j], Hist,Run[j][0], Label[0], "biasScan_dphcutB"+ss_dphcut[k]+"_"+ss_Angle[j]+"_resTreeCorr", "Resolution [#mum]",0.,7.,0.,8., "biasScan_"+ss_Angle[j], "Bias Voltage [V]"  );
    }


  ofstream myfileBiasUnf2[BiasVoltages];
  for(int j = 0; j < BiasVoltages; j++)
    {

      myfileBiasUnf2[j].open ("/home/zoiirene/Output/TextFiles/Ascan_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+func+"_"+ss_BiasVoltage[j]+"_resTreeCorr_"+testname+".txt");
      myfileBiasUnf2[j] << "A " << detectorA << "\n";
      myfileBiasUnf2[j] << labelA << "\n";
      myfileBiasUnf2[j] << "B " << detectorB<< "\n";
      myfileBiasUnf2[j] << labelB << "\n";
      myfileBiasUnf2[j] << "C " << detectorC<< "\n";
      myfileBiasUnf2[j] << labelC << "\n";

      myfileBiasUnf2[j] << "Angle(deg) RES(um) Error\n";

      for(int i=0; i<Angles; i++)
	{
	  myfileBiasUnf2[j] << ss_Angle[i] << " " << Resolution[i][j] << " " << ResolutionError[i][j] << "\n";
	}
      myfileBiasUnf2[j].close();
    }




  cout << " doing other plots " << endl ; 



  
  Double_t resolution[BiasVoltages][Angles];
  Double_t resolutionError[BiasVoltages][Angles];
  


  double errx[Angles];
  for(int i =0; i< Angles;i++)
    {
      errx[i]=0;
      for(int j=0; j<BiasVoltages; j++)
	{
	  cout << Angle[i] << " " << BiasVoltage[j] << " " << Resolution[i][j] << endl;
	  resolution[j][i] = Resolution[i][j];
	  resolutionError[j][i] = ResolutionError[i][j];
	  cout << Angle[i] << " " << BiasVoltage[j] << " " << resolution[j][i] << endl;
	  
	}
    }

  TCanvas *c2 = new TCanvas("c2", "c2", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //  gStyle->SetOptStat();
  //gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  int color;

  TLegend* leg2 = new TLegend(0.3,0.75,0.85,0.85);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.04);
  

  

  TGraphErrors * gr[BiasVoltages];
  for(int j=0; j<BiasVoltages; j++)
    {
      gr[j]  = new TGraphErrors(Angles, Angle, resolution[j],errx,resolutionError[j]);//, errx, rms);
    }
  leg2->AddEntry(gr[3],"#phi_{eq} = 4#times10^{15} cm^{-2}, neutron, 5.2 GeV", "");
  for(int j=0; j<BiasVoltages; j++)
    {
      
      color = j+1;
      if(j+1 ==5) color = 95;
      if(j+1 ==10) color = 29;
      
      gr[j]->SetLineColor(color);
      gr[j]->SetLineWidth(1);
      gr[j]->SetLineStyle(2);
      gr[j]->SetMarkerColor(color);
      gr[j]->SetMarkerSize(1.5);
      gr[j]->SetMarkerStyle(21);
      gr[j]->GetXaxis()->SetRangeUser(8.,22.);
      gr[j]->SetMinimum(0.);
      gr[j]->SetMaximum(10.);
      
      //  gr[j]->SetTitle("Option ACP example");
      gr[j]->GetXaxis()->SetTitle("Angle [deg]");
      gr[j]->GetYaxis()->SetTitle("Resolution [#mum]");
      //gr[j]->SetFillStyle(3003);
      //gr[j]->SetFillColor(kRed-8);
            leg2->AddEntry(gr[j],ss_BiasVoltage[j], "lp");      
	    cout  << "try  " << BiasVoltage[j] << " " << resolution[j][0] <<" " << resolution[j][1] << endl;

	    if(j==0)   gr[j]->Draw("ACPE");
	    else gr[j]->Draw("CPEsame");      
    }

  leg2->SetNColumns(4);

  leg2->Draw("same");
  TString name;
  name = outputDir+"compare_angle_bias_Scan_dphcutB"+ss_dphcut[k]+"_"+detectorB+"_"+Label[0]+"_"+Hist;
  c2->SaveAs(name+".eps");
  c2->SaveAs(name+".pdf");
  c2->SaveAs(name+".png");
  c2->SaveAs(name+".root");


  ////// zoom in 
  
  TCanvas *c2z = new TCanvas("c2z", "c2z", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //  gStyle->SetOptStat();
  //gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  
  for(int j=0; j<BiasVoltages; j++)
    {
      gr[j]  = new TGraphErrors(Angles, Angle, resolution[j],errx,resolutionError[j]);//, errx, rms);
      color = j+1;
      if(j+1 ==5) color = 95;
      if(j+1 ==10) color = 29;

      gr[j]->SetLineColor(color);
      gr[j]->SetLineWidth(1);
      gr[j]->SetLineStyle(2);
      gr[j]->SetMarkerColor(color);
      gr[j]->SetMarkerSize(1.5);
      gr[j]->SetMarkerStyle(21);
      gr[j]->GetXaxis()->SetRangeUser(8.,22.);
      gr[j]->SetMinimum(4.);
      gr[j]->SetMaximum(7.);

      //  gr[j]->SetTitle("Option ACP example");
      gr[j]->GetXaxis()->SetTitle("Angle [deg]");
      gr[j]->GetYaxis()->SetTitle("Resolution [#mum]");
      //gr[j]->SetFillStyle(3003);
      //gr[j]->SetFillColor(kRed-8);
      //      leg2->AddEntry(gr[j],ss_BiasVoltage[j], "lp");
      cout  << "try  " << BiasVoltage[j] << " " << resolution[j][0] <<" " << resolution[j][1] << endl;

      if(j==0)   gr[j]->Draw("ACPE");
      else gr[j]->Draw("CPEsame");
    }

  
  leg2->SetTextSize(0.03);
  leg2->Draw("same");
  
  name = outputDir+"compare_angle_bias_Scan_zoom_dphcutB"+ss_dphcut[k]+"_"+detectorB+"_"+Label[0]+"_"+Hist;
  c2z->SaveAs(name+".eps");
  c2z->SaveAs(name+".pdf");
  c2z->SaveAs(name+".png");
  c2z->SaveAs(name+".root");

  }

  cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%% cluster size "<< endl;
  ////// cluster size
  Hist = "nrowB_clphABC90evR";
  for(int j = 0; j < Angles; j++)
    {
      res_map[j] =  GetCheckHists(&(res_map[j]), BiasVoltages, dphcuts,dphcut, Run[j],ss_dphcut,pitch, Hist,hdx3_clchargeABC90evR[j],true,testname);
    }

  TString myfile3name[dphcuts][BiasVoltages];
  for(int k =0; k<dphcuts;k++){
    ofstream myfile3[BiasVoltages];
    for(int i =0; i< BiasVoltages; i++){
      myfile3name[k][i]= "/home/zoiirene/Output/TextFiles/Ascan_clsizeB_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+ss_BiasVoltage[i]+"_"+testname+"_ABC90.txt";
      cout << myfile3name[k][i] << endl;
      myfile3[i].open (myfile3name[k][i]);
      myfile3[i] << "A " << detectorA << "\n";
      myfile3[i] << labelA << "\n";
      myfile3[i] << "B " << detectorB<< "\n";
      myfile3[i] << labelB << "\n";
      myfile3[i] << "C " << detectorC<< "\n";
      myfile3[i] << labelC << "\n";
      myfile3[i] << info << "\n";
      myfile3[i] << "Angle MeanNrowB Error\n";
    }

    double NrowB[Angles][BiasVoltages];
    double  NrowBerror[Angles][BiasVoltages];
    

    for(int j = 0; j < Angles; j++)
      {
	
	for(int i=0; i<BiasVoltages; i++)
	  {
	    //auto it2 = res_map[j].find(std::make_pair(Run[j][i]+"_"+ss_dphcut[k]+"_"+Label[0],Hist));
	    auto it2 = res_map[j].find(std::make_pair(Run[j][i]+"_"+ss_dphcut[k]+"_"+pitch,Hist));
	    if(it2  != res_map[j].end())
	      {
		if(print)  cout << " found map " << endl;
		if(print)  cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetMean() << endl;
		NrowB[j][i] = it2->second->GetMean();
		NrowBerror[j][i] = it2->second->GetMeanError();
		myfile3[i] << " " << Angle[j] << " " << NrowB[j][i]<< " " << NrowBerror[j][i] << "\n";
		cout << Angle[j] << " " << NrowB[j][i]<< " " << NrowBerror[j][i] << endl;
	      }else{
	      cout << "not found map " << Run[j][i]+"_"+Label[0] << " hist " << Hist << endl;
	    }

	  }
      }

    for(int i=0; i<BiasVoltages; i++)
      myfile3[i].close();

  }// k dphcuts
  ////// cluster size
  /*
  Hist = "hnrowB_ph95";
  TH1F * hnrowB_clphABC90evR[Angles]; // = new TH1I("dx3_clchargeABC90evR ", "triplet dx_clchargeABC90evR ; dx [mm];triplets", 500, -0.5, 0.5 ); //Cut at 90% events in Landau (only high tail)
  for(int j = 0; j < Angles; j++)
    hnrowB_clphABC90evR[j] =   GetCheckHists(&(res_map[j]), BiasVoltages, dphcuts,dphcut, Run[j],ss_dphcut,pitch, Hist,hnrowB_clphABC90evR[j],true,"testEntries");
  */
  Hist = "nrowB";
  for(int j = 0; j < Angles; j++)
    GetHists(&(res_map[j]),BiasVoltages, dphcuts, comparisons,dphcut, Run[j], i_dphcut, Label, Hist, h_res,true,testname);


  for(int k =0; k<dphcuts;k++){
  ofstream myfile3[BiasVoltages];
  for(int i =0; i< BiasVoltages; i++){
    myfile3[i].open ("/home/zoiirene/Output/TextFiles/Ascan_clsizeB_"+detectorB+"_dphcutB"+ss_dphcut[k]+"_"+ss_BiasVoltage[i]+"_"+testname+".txt");
    myfile3[i] << "A " << detectorA << "\n";
    myfile3[i] << labelA << "\n";
    myfile3[i] << "B " << detectorB<< "\n";
    myfile3[i] << labelB << "\n";
    myfile3[i] << "C " << detectorC<< "\n";
    myfile3[i] << labelC << "\n";
    myfile3[i] << info << "\n";
    myfile3[i] << "Angle MeanNrowB Error\n";
  }

  TCanvas *csza[Angles];/// = new TCanvas("csza", "csza", 1500, 900);

  double NrowB[Angles][BiasVoltages];
  double  NrowBerror[Angles][BiasVoltages];



  TLegend* leg2a = new TLegend(0.7,0.35,0.85,0.85);
  leg2a->SetLineColor(0);
  leg2a->SetTextSize(0.04);

  
  for(int j = 0; j < Angles; j++)
    {
      csza[j] = new TCanvas("csza", "csza", 1500, 900);  
      gPad->SetTicks(1,1);
      gROOT->SetStyle("Plain");
      gStyle->SetPadGridX(0);
      gStyle->SetPadGridY(0);
      gStyle->SetPalette(1);

      //gStyle->SetOptStat(1110);
      gStyle->SetOptTitle(0);
      
      for(int i=0; i<BiasVoltages; i++)
	{
	  gStyle->SetOptStat(0);
	  //auto it2 = res_map[j].find(std::make_pair(Run[j][i]+"_"+ss_dphcut[k]+"_"+pitch[0],Hist));
	  auto it2 = res_map[j].find(std::make_pair(Run[j][i]+"_"+ss_dphcut[k]+"_"+Label[0],Hist));
	  if(it2  != res_map[j].end())
	    {
	      if(print)  cout << " found map " << endl;
	      if(print)  cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetMean() << endl;
	      NrowB[j][i] = it2->second->GetMean();
	      NrowBerror[j][i] = it2->second->GetMeanError();
	      myfile3[i] << " " << Angle[j] << " " << NrowB[j][i]<< " " << NrowBerror[j][i] << "\n";
	    }else{
	    cout << "not found map " << Run[j][i]+"_"+Label[0] << " hist " << Hist << endl;
	  }
	  int color = i+1;
	  if(i+1 ==5) color = 95;
	  if(i+1 ==10) color = 29;
	  if(print)  cout << " color " << color << endl;
	  it2->second->SetLineColor(color);
	  it2->second->SetLineWidth(1);
	  it2->second->SetLineStyle(1);
	  //	  it2->second->SetMarkerColor(color);
	  //it2->second->SetMarkerSize(1.5);
	  //it2->second->SetMarkerStyle(21);
	  it2->second->GetXaxis()->SetRangeUser(0.,5.);
	  it2->second->SetMinimum(0.);
	  it2->second->SetMaximum(40000.);

	  //  it2->second->SetTitle("Option ACP example");
	  it2->second->GetXaxis()->SetTitle("cluster size B");
	  it2->second->GetYaxis()->SetTitle("Entries");
	  //it2->second->SetFillStyle(3003);
	  //it2->second->SetFillColor(kRed-8);
	  //      leg2->AddEntry(it2->second,ss_BiasVoltage[j], "lp");
	  if(j==0)          leg2a->AddEntry(it2->second,ss_BiasVoltage[i], "l");

	  it2->second->Draw("histsame");
	  //	  it2->second->Draw("CPEsame");

	  
	}
      leg2a->SetNColumns(2);
      leg2a->SetTextSize(0.03);
      leg2a->Draw("same");

     TString  name = outputDir+"compare_angle"+ss_Angle[j]+"_bias_Scan_clustersize_dphcutB"+ss_dphcut[k]+"_"+detectorB+"_"+Label[0]+"_"+Hist;
      csza[j]->SaveAs(name+".eps");
      csza[j]->SaveAs(name+".pdf");
      csza[j]->SaveAs(name+".png");
      csza[j]->SaveAs(name+".root");
    }


  for(int i =0; i< BiasVoltages; i++)  myfile3[i].close();  
  TCanvas *csz = new TCanvas("csz", "csz", 1500, 900);
  gPad->SetTicks(1,1);
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPalette(1);
  //  gStyle->SetOptStat();
  //gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);


  TGraph * grB[BiasVoltages];
  for(int j=0; j<BiasVoltages; j++)
    {
       
      grB[j]  = new TGraph(Angles, Angle, NrowB[j]);
      int color = j+1;
      if(j+1 ==5) color = 95;
      if(j+1 ==10) color = 29;

      grB[j]->SetLineColor(color);
      grB[j]->SetLineWidth(1);
      grB[j]->SetLineStyle(2);
      grB[j]->SetMarkerColor(color);
      grB[j]->SetMarkerSize(1.5);
      grB[j]->SetMarkerStyle(21);
      grB[j]->GetXaxis()->SetRangeUser(8.,22.);
      grB[j]->SetMinimum(0.);
      grB[j]->SetMaximum(4.);

      //  grB[j]->SetTitle("Option ACP example");
      grB[j]->GetXaxis()->SetTitle("Angle [deg]");
      grB[j]->GetYaxis()->SetTitle("Mean cluster size B");
      //grB[j]->SetFillStyle(3003);
      //grB[j]->SetFillColor(kRed-8);
      //      leg2->AddEntry(grB[j],ss_BiasVoltage[j], "lp");


      if(j==0)   grB[j]->Draw("ACPE");
      else grB[j]->Draw("CPEsame");
	
    }

  
  leg2a->SetTextSize(0.03);
  leg2a->Draw("same");
  
TString  name = outputDir+"compare_angle_bias_Scan_clustersize_dphcutB"+ss_dphcut[k]+"_"+detectorB+"_"+Label[0]+"_"+Hist;
  csz->SaveAs(name+".eps");
  csz->SaveAs(name+".pdf");
  csz->SaveAs(name+".png");
  csz->SaveAs(name+".root");
  }









  
}//resolution 


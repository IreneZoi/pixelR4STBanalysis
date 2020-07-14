#include <stdio.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <sstream>
#include "fileHandler.h"

bool print=true;
using namespace std;
#define measurements 7
#define irradiation 3
#define dphcuts  1
#define comparisons 1

void bestRes(TString func = "RMSself", bool unfolding = false)
{
  TString unf = "Unfolded";
  TString filename = "/home/zoiirene/Output/TextFiles/resolution_25gain1_"+func;
  if(unfolding)
    filename+=unf;
  filename+=".txt";
  if(print) cout << filename << endl;


  TString sensor[measurements] = {
    //  "108", 
    "120i",
    "150", 
    "130i",
    "148", 
    "146", 
    "163",
    "194i"
    // "152", 
    // "159", 
    // "102", 
    // "160", 
    // "133i"
  };
  Int_t irr[measurements]= {
    //0, 
    2, 
    0, 
    2, 
    0, 
    0, 
    0,
    4// , 
    // 0, 
    // 0, 
    // 0, 
    // 0, 
    // 2
  };
  int pitch[measurements]={
    //25,
    25,
    25,
    25,
    25,
    25,
    25,
    25
    // 50,
    // 50,
    // 50,
    // 50,
    // 50
  };
  TString bias[measurements]={
    //    "200", 
    "800", 
    "120", 
    "600", 
    "120", 
    "120", 
    "120",
    "800"//, 
    //"120"// , 
    // "?",   
    // "?",   
    // "120", 
    // "600"
  };
  TString angle[measurements]={
    //    "12.5", 
    "10 ", 
    "12.5",
    "11.25",
    " 8.75 ",
    "15 ",
    "11.25",
    "12"
    // "17,50",
    // "19.40",
    // "20.50",
    // "17.00",
    // "20"  
  };
 float beam[measurements]={
   //5.6,//6
   5.6,  // no 6 GeV data!!! 
   5.6,  //6 to be analysed.. is from different set of data 
   5.6, //6
   5.6, //6 to be analysed.. is from a different set of data       
   5.6,  //6 to be analysed.. is from a different set of data 
   5.6, //6
   5.6
   // 5.6,//6
   // 5.6,//6
   // 5.6, //5.6 available at not optimal angle (but close) but with only 40k events, 6 available at optimal angle     
   // 5.6
   //6
 };  
 TString runs[measurements] = {
   //   "1827", 
   "2801",      
   "2775",      
   "1820", 
   "2743",      
   "2773",      
   "2832",
   "3839"
   // "1027", 
   // "1894", 
   // "1870",      
   // "1037", 
   // "1846"
   
 };
 TString filelabel[measurements] = {
   //
   "thrScan_A12C15",
   "closest_A11C15",
   "closest_A32C42",
   "thrScan_A13C14",
   "closest_A11C14",
   "closest2_A11C12",
   "thrScan_A12C13"
   //
   //
   //
   //
   //
 };
 TString specific[measurements] = {
   //"FTH150P_6_R4S100x25-P1_2", 
   "FTH150P_6_R4S100x25-P1_2",
   "FDB150P_12_R4S100x25-P1_3",
   "FDB150P_12_R4S100x25-P1_4",
   "FDB150P_12_R4S100x25-P4_1",
   "FDB150Y_2_R4S100x25-Y2_2",
   "FDB150Y_2_R4S100x25-Y6_1",
   "FDB150P_12_R4S100x25-P1_1"
   // "FDD150P_22_R4S50x50-P1_1",
   // "FDB150P_12_R4S50x50-P3_1",
   // "FTH150P_6_R4S50x50-P3_2",
   // "FDB150P_12_R4S50x50-P4_1", 
   // "FDB150Y_2_R4S50x50-Y2_1" 
 };
 TString Short[measurements] = {
   //   "Pstop_default_FTH",
   "Pstop_default_FTH",   
   "Pstop_default_FDB",   
   "Pstop_default_FDB",   
   "Pstop_RD53Apads_FDB", 
   "Pspray_default_FDB",  
   "Pspray_RD53Apads_FDB",
   "Pstop_default_FDB"
   // "Pstop_default_FDD",   
   // "Pstop_bdotlarge_FDB", 
   // "Pstop_bdotlarge_FTH", 
   // "Pstop_bdotwiggle_FDB",
   // "Pspray_default_FDB" 
 };
 int thr[measurements] = {
   //15, 
   15, 
   12, 
   20, 
   12, 
   13, 
   14,
   15
   // 12, 
   // 20, 
   // 20, 
   // 12, 
   // 25
 };
 Double_t Resolution[measurements],ResolutionError[measurements],Percentage[measurements];
 Double_t RMS[measurements],RMSError[measurements];
 TString search;
 TString comment;
 MapTH1 res_map;

 TString ss_dphcut[dphcuts]; 
 
 TString Label[comparisons] = {"straightTracksY_isoAandCandB_straightTracksX"};
 TString Hist = "dx3_clchargeABC90evR";
 bool dphcut = true;
 TH1F * h_res;
 //GetHists(&res_map,measurements, 1, comparisons,dphcut, runs, thr, Label, Hist, h_res);
 TString ss_thr[measurements];

 TString ss_pitch = "25";
 TH1F * hdx3_clchargeABC90evR; 
 for(int i = 0; i<measurements;i++){
   cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     " << runs[i] << endl;
   TString Run[1] = {runs[i]};
   ss_thr[i].Form("%d",thr[i]);
   TString ss_dphcut[1] ={ss_thr[i]};
   TString label = filelabel[i];
   res_map =    GetCheckHists(&res_map, 1, dphcuts,dphcut, Run,ss_dphcut,ss_pitch, Hist,hdx3_clchargeABC90evR,true,label);
 }
 
 /*
 auto it2 = res_map.find(std::make_pair(runs[0]+"_"+ss_dphcut[0]+"_"+Label[0],Hist));
 if(print)               cout << " found map " << endl;
 if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;
 */
 cout << " ################################################################# " << endl;
 double freshres = 3.32;
 double freshres_err = 0.02;
 for(int i=0; i<measurements; i++)
   {
     cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     " << runs[i] << endl;     
     ss_thr[i].Form("%d",thr[i]);
     cout << ss_thr[i] << endl;
     auto it2 = res_map.find(std::make_pair(runs[i]+"_"+ss_thr[i]+"_"+ss_pitch,Hist)); //run dphcut pitch

     if(it2  != res_map.end())
       {

	 if(print)               cout << " found map " << endl;
	 if(print)                   cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;

	 FitTH1(it2->second, &(RMS[i]), &(RMSError[i]),  "dphcut"+ss_thr[i], runs[i], Label[0], Hist, func,&(Percentage[i]));

	 if(print) cout << "run " << i<< ": " << runs[i] << " RMS: " << RMS[i] << " and res err: " << RMSError[i] << endl;

	 Resolution[i] = RMS[i];
	 ResolutionError[i] = RMSError[i];


	 if(irr[i]!=0 )	 ExtractRes(&(Resolution[i]), &(ResolutionError[i]), true,freshres,freshres_err);
	 else if (i!=0)  ExtractRes(&(Resolution[i]), &(ResolutionError[i]));
	 if(print) cout << "run " << i<< ": " << runs[i] << " Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;

	 
       }
   }


 
  ofstream myfile;
  myfile.open (filename);
  myfile << "sensor  irr[e15 neq]    pitch[um]       bias[v] angle[deg]      beam[GeV]       res[um] reserr[um] runs    specific        short   thr[ADC]  \n";
  for(int i=0; i<measurements; i++)
    {
      myfile <<  sensor[i] << std::setw(3)
	     <<  irr[i]    << std::setw(3)
	     <<  pitch[i]  << "   " //std::setw(3)
	     <<  bias[i]   <<"   "
	     <<  angle[i]  <<"   "
	     <<  beam[i]   <<"   "
	     <<  Resolution[i] <<"   "
	     <<  ResolutionError[i] <<"   "
	     <<  runs[i]   <<"   "
	     <<  specific[i] <<"   "
	     <<  Short[i]  <<"   "
	     <<  thr[i]    <<" \n  ";
      
    }

  /*  
  myfile << "B " << detectorB<< "\n";
  myfile << labelB << "\n";
  myfile << "C " << detectorC<< "\n";
  myfile << labelC << "\n";
  myfile << info << "\n";
  myfile << "Momentum(GeV) Resolution(um) Error\n";

      FitTH1(h_res[i], &(Resolution[i]), &(ResolutionError[i]), "6", detectorA, detectorB, detectorC,func );
      if(print) cout << "Beam energy " << i<< ": 6 GeV -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      ResolutionSquare[i]=Resolution[i]*Resolution[i];
      ResolutionErrorSquare[i]=2*ResolutionError[i]*Resolution[i];
      myfile << "6  " << Resolution[i] << " " << ResolutionError[i] << "\n";

    }
  */
  myfile.close();

  




}  

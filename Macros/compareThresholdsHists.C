#include <stdio.h>
#include "TCanvas.h"
#include "TF1.h"
#include <TTreeFormula.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TVectorT.h>
#include <TMap.h>
#include <TAxis.h>
#include <sstream>
#include "fileHandler.h"

#define comparisons 1
#define runs 1
#define dphcuts 11
#define hists 6
bool print=true;
using namespace std;
TString outputDir = "/home/zoiirene/Output/Plots/";



void compareThresholdsHists()//TString hist1 = "clszBiii",TString hist2 = "clszBiii",TString run1 = "2743", TString label1 = "p-stop", TString run2 = "2832", TString label2 = "p-spray", bool norm = true)
{

  TString inputDir="/home/zoiirene/Output/TextFiles/";
  TString inputfile, Path;
 
  TFile * file[runs];
  TDirectory * histodir;
  TString Run[runs];
  Run[0] = "1827";
  TString Hist[hists];
  Hist[0] = "dx3";
  Hist[1] = "clchargeB";
  Hist[2] = "clphB";
  Hist[3] = "clsizeB";
  Hist[4] = "nrowvsxmB3";
  //  Hist[5] = "dx3_clchargeAC90evR";
  Hist[5] = "dx3_clchargeABC90evR";
  int i_dphcut[dphcuts];
  //  i_dphcut[0] =  10 ;
  // i_dphcut[1] =  11 ;
  //i_dphcut[2] =  12 ;
  i_dphcut[0] =  13 ;
  i_dphcut[1] =  14 ;
  i_dphcut[2] =  15 ;
  i_dphcut[3] =  18 ;
  i_dphcut[4] =  20 ;
  i_dphcut[5] =  22 ;
  i_dphcut[6] =  25 ;
  i_dphcut[7] =  27 ;
  i_dphcut[8] =  30 ;
  i_dphcut[9] =  35 ;
  i_dphcut[10] =  40 ;

  double sigma[dphcuts];
  double sigmaerr[dphcuts];
  double RMS[dphcuts];
  double RMSerr[dphcuts];

  TString ss_dphcut[dphcuts];
  for(int i =0; i< dphcuts;i++)
    ss_dphcut[i].Form("%d",i_dphcut[i]);

  TString Label[comparisons];
  /*
  Label[0] = "raw";
  Label[1] = "dxCAcut";
  */  
  //Label[0] = "beforeCorrections";
  //Label[1] = "nocuts";
  //Label[2] = "straightTracksY";
  //Label[3] = "straightTracksY_isoAandC";
  //Label[4] = "straightTracksY_isoAandCandB";
  //Label[5] = "straightTracksY_isoAandCandB_straightTracksX";
  Label[0] = "straightTracksY_isoAandCandB_straightTracksX";
  //Label[6] = "hitsOnTrack";
  /*
  Label[7] = "straightTracksY_isoAandCandB_chargerAandC";
  Label[8] = "straightTracksY_isoAandCandB_chargerAandCandB";
  Label[9] = "straightTracksY_isoAandCandB_chargeAandC";
  Label[10] = "straightTracksY_isoAandCandB_chargeAandCandB";
  */  
bool dphcut = true;

  MapTH1 m_hists;
  std::map<std::pair<TString, TString>, TH1F *>::iterator it;
  TH1F * h[hists];
  TString outputDir = "/home/zoiirene/Output/Plots/";

  TString name;

  for(int k=0; k<hists; k++)
    {
      GetHists(&m_hists,runs, dphcuts, comparisons,  dphcut, Run, i_dphcut, Label,Hist[k], h[k], true, "AC15");

      for ( it = m_hists.begin(); it != m_hists.end(); it++ )
	{
	  it->second->GetEntries();
	}
    }

  if(print) cout << " ***************** GOT ALL HISTS ********************** " << endl;
  ofstream myfile[hists];


  for(int k=0; k<hists; k++)
    {
      if(k == 5 || k ==3 )
	{

	  name = outputDir+"compare_thresholdsAC15_RMS_"+Run[0]+"_"+Label[0]+"_"+Hist[k];//+"_"+Hist;
	  myfile[k].open (name+".txt");
	  
	  myfile[k] << "\\begin{tabular}{|c|c|c|c|}\t\n";
	  myfile[k] << "\\hline \n";
	  myfile[k] << "dphcut &  \\multicolumn{3}{ |c| }{" << Hist[k] << "}\\\\ \n";
	  myfile[k] << "\\hline \n";
	  
	  myfile[k] << " & entries&  sigma [$\\mu m$]& sigma err [$\\mu m$]\\\\ \n";
	  myfile[k] << "\\hline \n";
	}
    }

  if(print) cout << " Initialized file " << name << endl;
  
  for(int i=0; i<runs; i++)
    {
      for(int l=0; l<dphcuts; l++)
	{
	  
	  for(int j=0; j<comparisons; j++)
	    {
	      
	      for(int k=0; k<hists; k++)
		{
		  if(k == 5 || k ==3 )
		    {
		      myfile[k] << i_dphcut[l] ;
		      
		      auto it2 = m_hists.find(std::make_pair(Run[i]+"_"+ss_dphcut[l]+"_"+Label[j],Hist[k]));
		      if(it2  != m_hists.end())
			{
			  if(print)               cout << " found map " << endl;
			  if(print)		      cout << "map key " << it2->first.first << " " << it2->first.second << " " << it2->second->GetEntries() << endl;  

			  FitTH1(it2->second, &(sigma[l]), &(sigmaerr[l]), "dphcut"+ss_dphcut[l], Run[i], Label[j], Hist[k], "RMS");
			  myfile[k] <<  "& " <<  setprecision(5) <<  it2->second->GetEntries() << "& " << setprecision(3) <<sigma[l] << "& " << setprecision(1) << sigmaerr[l];
                          if(k==5){
			    RMS[l] = sigma[l];
			    RMSerr[l] = sigmaerr[l];
			  }
			  if(print) cout << " ***************** done RMS ********************** " << endl;
			}
		      myfile[k] << "\\\\ \n";
		      myfile[k] << "\\hline \n";

		    }  
		      
		}//hists
	      
	    }//dphcuts
	}//comparisons
    }//runs
  
  for(int k=0; k<hists; k++)
    {
      if(k == 5 || k ==3)
	{
	  
	  myfile[k]<< "\\end{tabular} \n";
	  
	  
	  myfile[k].close();
	}
    }
  

  ///pasting

  TString detectorA="148";
  TString detectorB="108";
  TString detectorC="109";
  TString labelA="Pstop_RD53Apads_FDB, thr 30 ADC";
  TString labelB="Pstop_default_FTH, 200V";
  TString labelC="Pstop_default_FTH, thr 30 ADC";
  TString info = "beam energy 5.6 GeV";
  
  double Resolution[dphcuts];
  double ResolutionError[dphcuts];
  
  
  ofstream myfile2;
  myfile2.open ("/home/zoiirene/Output/TextFiles/ThrscanAC15_res_"+detectorB+"_"+Run[0]+".txt");
  myfile2 << "A " << detectorA << "\n";
  myfile2 << labelA << "\n";
  myfile2 << "B " << detectorB<< "\n";
  myfile2 << labelB << "\n";
  myfile2 << "C " << detectorC<< "\n";
  myfile2 << labelC << "\n";
  myfile2 << info << "\n";
  myfile2 << "dphcut(ADC) res(um) Error\n";

  for(int i=0; i<dphcuts; i++)
    {

      Resolution[i] = RMS[i];
      ResolutionError[i] = RMSerr[i];
      if(print) cout << "thr " << i<< ": " << ss_dphcut[i] << " ADC -> RMS: " << sigma[i] << " and res err: " << sigmaerr[i] << endl;

      ExtractRes(&(Resolution[i]), &(ResolutionError[i]), false); //,freshres,freshres_err);
      if(print) cout << "thr " << i<< ": " << ss_dphcut[i] << " ADC -> Resolution: " << Resolution[i] << " and res err: " << ResolutionError[i] << endl;
      myfile2 << " " << ss_dphcut[i] << " " << Resolution[i] << " " << ResolutionError[i] << "\n";
    }

  myfile2.close();



			    


  
  ////////


  
  

  if(print) cout << " ***************** wrote tab file ********************** " << endl;
  
  if(print) cout << "Plotting"  << endl;
  float xmin[hists]={-0.05,0,0,0,0,9};//,-0.05,-0.05};
  float xmax[hists]={0.05,50,600,10,25,41};//,0.05,0.05};
  float ymin[hists]={0,0,0,0,0,0.};//,0,0};
  float ymax[hists]={6.,4000,2500,50000,4.,4.};//,10000,10000};

  DrawThrScanHists(&m_hists,dphcuts, Run[0], i_dphcut, Label[0], hists,Hist,xmin,xmax,ymin,ymax);
  DrawThresholdLegend(&m_hists,dphcuts, Run[0], ss_dphcut, Label[0], Hist[0]);
  
    
  if(print) cout << " ***************** DONEEEE ********************** " << endl;
 
}//resolution 

  //  Run[0] = "2775";2773
/*
i_dphcut[0] =  0 ;
  i_dphcut[1] =  9 ;
  i_dphcut[2] =  10 ;
  i_dphcut[3] =  11 ;
  i_dphcut[4] =  12 ;
  i_dphcut[5] =  13 ;
  i_dphcut[6] =  14 ;
  i_dphcut[7] =  15 ;
  i_dphcut[8] =  18 ;
  i_dphcut[9] =  20 ;
  i_dphcut[10] =  22 ;
  i_dphcut[11] =  25 ;
  i_dphcut[12] =  27 ;
  i_dphcut[13] =  30 ;
  i_dphcut[14] =  32 ;
  i_dphcut[15] =  35 ;
  i_dphcut[16] =  37 ;
  i_dphcut[17] =  40 ;
*/

/*

  Run[0] = "2743";//2832
i_dphcut[0] =  0 ;
  i_dphcut[1] =  6 ;
  i_dphcut[2] =  9 ;
  i_dphcut[3] =  10 ;
  i_dphcut[4] =  11 ;
  i_dphcut[5] =  12 ;
  i_dphcut[6] =  13 ;
  i_dphcut[7] =  14 ;
  i_dphcut[8] =  15 ;
  i_dphcut[9] =  18 ;
  i_dphcut[10] =  20 ;
  i_dphcut[11] =  22 ;
  i_dphcut[12] =  25 ;
  i_dphcut[13] =  27 ;
  i_dphcut[14] =  30 ;
  i_dphcut[15] =  32 ;
  i_dphcut[16] =  35 ;
  i_dphcut[17] =  37 ;
  i_dphcut[18] =  40 ;
*/

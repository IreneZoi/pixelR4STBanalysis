//void FitTH1(TH1F* h1, Double_t *  sigma, Double_t *  sigmaerr, TString name, TString detectorA, TString detectorB, TString detectorC, TString func );
#include<bits/stdc++.h>
#include "TRandom3.h"

Double_t tp0Fit( Double_t *x, Double_t *par );
double G = 1E0;


double poissonsmearing(TH1* h, TString name, TString detectorA, TString detectorB, TString detectorC, int nexp = 10000){
  TRandom3 *myrnd = new TRandom3();
  TH1F *sigmas = new TH1F("sigmas","",200,0.,10.);

  //getting original sigma for reference
  double sigma = h->GetRMS() * 1000;
  double sigmaerr = h->GetRMSError() * 1000;
  cout << "original resolution " << sigma << " ± " << sigmaerr << endl;

  // Repeat this pseudo experiment 1k times:
  for(size_t i = 0; i < nexp; i++) {

    // new histogram, clone of the original:
    TH1F * smeared = (TH1F*)h->Clone("smeared");

    // loop over the bins:
    for(Int_t bin = 1; bin <= h->GetNbinsX(); bin++) {
      // Get random number from Poisson distribution with mean of original bin content:
      Double_t newcontent = myrnd->Poisson(h->GetBinContent(bin));
      //cout << " " << i << " bin " << bin << " orig " << original->GetBinContent(bin) << " rnd " << content << endl;
      smeared->SetBinContent(bin,newcontent);
    }
    Double_t res = smeared->GetRMS() * 1000;

    if(i%1000 == 0) cout << "i" << i << " sigma = " << res << endl;
    sigmas->Fill(res,1);

    delete smeared;
  }

    // plot width distribution
    gROOT->SetStyle("Plain");
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    //  gPad->SetTickx();
    //gPad->SetTicky();

    TCanvas * cs = new TCanvas("cs","cs",700,700);

    cs->cd();
    cs->SetFrameFillStyle(1000);
    cs->SetFrameFillColor(0);
    cs->SetLeftMargin(0.15);
    cs->SetRightMargin(0.15);
    cs->SetBottomMargin(0.2);
    gPad->SetTicks(1,1);


    sigmas->SetTitle("");
    sigmas->GetXaxis()->SetLabelFont(42);
    sigmas->GetXaxis()->SetLabelSize(0.04);
    sigmas->GetXaxis()->SetTitleSize(0.05);
    sigmas->GetXaxis()->SetTitleOffset(0.8);
    sigmas->GetXaxis()->SetTitleFont(42);

    sigmas->GetYaxis()->SetLabelFont(42);
    sigmas->GetYaxis()->SetLabelSize(0.04);
    sigmas->GetYaxis()->SetTitleSize(0.05);
    sigmas->GetYaxis()->SetTitleOffset(1.);
    sigmas->GetYaxis()->SetTitleFont(42);

    sigmas->GetYaxis()->SetTitle("Entries");
    sigmas->GetXaxis()->SetTitle("width [mm]");

    //TF1 * MyGaus = new TF1("MyGaus","gaus", -0.06,0.06);
    sigmas->SetMarkerStyle(20);
    sigmas->Draw("PZ");
    
    cs->Update();
    gStyle->SetOptFit(1111);
    TString outputDir="/home/zoiirene/Output/Plots/";

    //  TString outputDir = outputDir+"";
    TString outputFile = outputDir+"width_Fit_"+name+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps";
    cs->SaveAs(outputFile);

    cout << "statistical uncertainty: " << sigmas->GetRMS() << endl;
    return sigmas->GetRMS();
    
  }


Double_t conversion = TMath::Sqrt(2./3.);

void ExtractRes(Double_t *  sigma, Double_t *  sigmaerr, bool isIRR = false, Double_t sigma_fresh = 0., Double_t sigmaerr_fresh = 0.)
{

  cout << " before unfolding " << *sigma << " ± " << *sigmaerr << " um " << endl;
  if(!isIRR)
    {
      *sigma *= conversion;
      *sigmaerr *= conversion;
    }
  
  if(isIRR)
    {
      Double_t sigma_orig = *sigma;
      Double_t sigmaerr_orig = *sigmaerr;
      cout << " TMath::Sqrt((2*sigma_orig*sigma_orig-sigma_fresh*sigma_fresh)/2.) " << endl;
      *sigma = TMath::Sqrt((2*sigma_orig*sigma_orig-sigma_fresh*sigma_fresh)/2.);
      *sigmaerr = TMath::Sqrt( 1./(*sigma)*sigma_orig*sigma_orig*sigmaerr_orig*sigmaerr_orig  + 0.25/(*sigma)*sigma_fresh*sigma_fresh*sigmaerr_fresh*sigmaerr_fresh );
      
    }
    
  cout << "after unfolding " << *sigma << " ± " << *sigmaerr << " um " << endl;
  

}

void FitTH1(TH1* h1, Double_t *  sigma, Double_t *  sigmaerr, TString name, TString detectorA, TString detectorB, TString detectorC, TString func, Double_t * percentage, Double_t *min, Double_t * max, Double_t nsigma = 6)
{
  cout << " evaluating resolution with method " << func << endl;
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
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetBottomMargin(0.2);
  gPad->SetTicks(1,1);


  h1->SetTitle("");
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetLabelSize(0.04);
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetTitleOffset(0.8);
  h1->GetXaxis()->SetTitleFont(42);

  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelSize(0.04);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleOffset(1.);
  h1->GetYaxis()->SetTitleFont(42);

  h1->GetYaxis()->SetTitle("Entries");
  h1->GetXaxis()->SetTitle("residual [mm]");

  //TF1 * MyGaus = new TF1("MyGaus","gaus", -0.06,0.06);
  h1->SetMarkerStyle(20);
  h1->Draw("PZ");
  //h1->Fit("MyGaus","R");
  if(func == "gaus")
    {
      h1->GetXaxis()->SetRangeUser(-0.1,0.1);

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
      h1->GetXaxis()->SetRangeUser(-0.1,0.1);

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
  if(func == "FWHM")
    {
      cout <<func <<  " method" << endl;
      double maximum = h1->GetMaximum();
      double halfmaximum = h1->GetMaximum()/2.;
      int maxbin = h1->GetMaximumBin();
      cout <<  " maxbin " << maxbin << " at " << h1->GetBinCenter(h1->GetMaximumBin())<<endl;
      cout << "intital sigma = " << h1->GetRMS() * 1000 << " sigmaerr = " << h1->GetRMSError() * 1000 << endl;
      

      int bin1 = h1->FindFirstBinAbove(halfmaximum);
      int bin2 = h1->FindLastBinAbove(halfmaximum);
      cout << " max " << maximum << " max/2 " << halfmaximum << " bin1 " << h1->GetBinContent(bin1) << " bin2 " << h1->GetBinContent(bin2) << endl;
      double fwhm = h1->GetBinCenter(bin2) - h1->GetBinCenter(bin1);
      cout << "fwhm " << fwhm << " low " << bin1 << " high " << bin2 << endl;


      *sigma = fwhm * 1000;
      *sigmaerr = poissonsmearing( h1, name,  detectorA, detectorB, detectorC);
      cout << " resolution " << *sigma << " ± " << *sigmaerr << endl;
      TString mean;
      std::ostringstream sstream_mean;
      sstream_mean << setprecision(2) << h1->GetMean();// * 1000;
      mean = sstream_mean.str();
      TString meanerr;
      std::ostringstream sstream_meanerr;
      sstream_meanerr << setprecision(1) << h1->GetMeanError();// * 1000;
      meanerr = sstream_meanerr.str();
      
      TString rms;
      std::ostringstream sstream_rms;
      sstream_rms << setprecision(3) << h1->GetRMS() * 1000;
      rms = sstream_rms.str();
      
      TString rmserr;
      std::ostringstream sstream_rmserr;
      sstream_rmserr << setprecision(1) << h1->GetRMSError() * 1000;
      rmserr = sstream_rmserr.str();
      //      TGaxis::SetMaxDigits(2)      ;
      gStyle->SetOptStat(0);
      TLatex Tl2;
      Tl2.SetTextAlign(12);
      Tl2.SetTextSize(0.04);
      Tl2.DrawLatexNDC(0.23,0.7," #mu = "+mean+" mm");
      Tl2.DrawLatexNDC(0.23,0.62," #sigma = "+rms+" #pm "+rmserr+" #mum");
      //Tl2.DrawLatexNDC(0.6,0.54,"|#Delta #eta |< 1.3");
      //Tl2.DrawLatexNDC(0.6,0.51,"m_{jj} > 1050 GeV");
      //Tl2.DrawLatexNDC(0.6,0.48,"65 GeV < M_{SD} < 105 GeV");
      

    }
  if(func=="RMSself"){
    cout << func << endl;
    ///////////////////////////////
    // Self consistent RMS in N RMS
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>

    TH1* hSubR = (TH1*)h1->Clone();
    Double_t N = nsigma;
    Double_t NFull = h1->Integral(0,h1->GetNbinsX()+1);
    Double_t NSubR = NFull;
    Double_t prev_rms = 0;
    Double_t this_rms = 0;
    Double_t x1 = 0;
    Double_t x2 = 0;
    Int_t cnt = 0;

    do{

      cnt++;
      prev_rms = this_rms;
      this_rms = hSubR->GetRMS();
      x1 = hSubR->GetMean() - this_rms * N;
      x2 = hSubR->GetMean() + this_rms * N;


      hSubR->GetXaxis()->SetRangeUser(x1,x2);
      NSubR = hSubR->Integral();
      *percentage = NSubR / NFull;
          cout << "Iterations   " << cnt << endl
	             << "  RMS        " << this_rms << endl
               << "NSubR " << NSubR << " NFull " << NFull << endl 
	       << "  Percentage " << *percentage << endl
	       << "  New Range  " << x1 << " to " << x2 << endl << endl;


    } while( this_rms != prev_rms );



    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // Self consistent RMS in n RMS
    ///////////////////////////////
    h1->GetXaxis()->SetRangeUser(x1,x2);
    *min = x1*1000;
    *max = x2*1000;
    cout << " min " << *min << " max " << *max << endl;
    cout << "Getting RMS " << endl;
    *sigma = h1->GetRMS() * 1000;
    *sigmaerr = h1->GetRMSError() * 1000;
    cout << " resolution " << *sigma << " ± " << *sigmaerr << endl;
    *sigmaerr = poissonsmearing( h1, name,  detectorA, detectorB, detectorC);
    cout << " resolution updated " << *sigma << " ± " << *sigmaerr << " with precision " << *percentage<<  endl;

    TString mean;
    std::ostringstream sstream_mean;
    sstream_mean << setprecision(2) << h1->GetMean();// * 1000;
    mean = sstream_mean.str();
    TString meanerr;
    std::ostringstream sstream_meanerr;
    sstream_meanerr << setprecision(1) << h1->GetMeanError();// * 1000;
    meanerr = sstream_meanerr.str();

    TString rms;
    std::ostringstream sstream_rms;
    sstream_rms << setprecision(3) << h1->GetRMS() * 1000;
    rms = sstream_rms.str();

    TString rmserr;
    std::ostringstream sstream_rmserr;
    sstream_rmserr << setprecision(1) << *sigmaerr;
    cout << sstream_rmserr.str() << endl;
    rmserr = sstream_rmserr.str();
    cout << rmserr << endl;
    //      TGaxis::SetMaxDigits(2)      ;
    gStyle->SetOptStat(0);
    TLatex Tl2;
    Tl2.SetTextAlign(12);
    Tl2.SetTextSize(0.04);
    Tl2.DrawLatexNDC(0.23,0.7," #mu = "+mean+" mm");
    Tl2.DrawLatexNDC(0.23,0.62," #sigma = "+rms+" #pm "+rmserr+" #mum");
    //Tl2.DrawLatexNDC(0.6,0.54,"|#Delta #eta |< 1.3");
    //Tl2.DrawLatexNDC(0.6,0.51,"m_{jj} > 1050 GeV");
    //Tl2.DrawLatexNDC(0.6,0.48,"65 GeV < M_{SD} < 105 GeV");
    cout << "legend done " << endl;
    

    
  }
  if(func == "RMS") 
    {
      cout << " RMS method" << endl;
      double integral = h1->GetMaximum();
      int maxbin = h1->GetMaximumBin();
      cout <<  " maxbin " << maxbin << " at " << h1->GetBinCenter(h1->GetMaximumBin())<<endl;
      cout << "intital sigma = " << h1->GetRMS() * 1000 << " sigmaerr = " << h1->GetRMSError() * 1000 << endl;

      cout << " integral " << h1->Integral() << " entries " << h1->GetEntries() << endl;
      double integral95 = 0.95*h1->Integral();
      cout << " integral95 " << integral95 << endl;
      int i = 0;

      double low = 0;
      double high = 0;
      for(int i =0; i<h1->GetNbinsX()/2; i++)
	{
	  low = maxbin-i;
	  high = maxbin+i;

	  integral = h1->Integral(low,high);
	  cout << " integral " << integral << " low " << low << " high " << high << endl;
	  cout << " while "<< i << " fabs(integral-integral95)/integral95 " << fabs(integral-integral95)/integral95 << endl;

	  if(fabs(integral-integral95)/integral95 < 0.011 || integral>integral95)//integral>integral95)
	    break;

	 
	}
	  cout << "final integral " << integral << " low " << low << " high " << high << endl;


      cout << "final integral " << integral << " low " << low <<" high " << high << endl;
      h1->GetXaxis()->SetRange(low,high);

      *sigma = h1->GetRMS() * 1000;
      *sigmaerr = h1->GetRMSError() * 1000;
      cout << " resolution " << *sigma << " ± " << *sigmaerr << endl;
      TString mean;
      std::ostringstream sstream_mean;
      sstream_mean << setprecision(2) << h1->GetMean();// * 1000;
      mean = sstream_mean.str();
      TString meanerr;
      std::ostringstream sstream_meanerr;
      sstream_meanerr << setprecision(1) << h1->GetMeanError();// * 1000;
      meanerr = sstream_meanerr.str();
      
      TString rms;
      std::ostringstream sstream_rms;
      sstream_rms << setprecision(3) << h1->GetRMS() * 1000;
      rms = sstream_rms.str();
      
      TString rmserr;
      std::ostringstream sstream_rmserr;
      sstream_rmserr << setprecision(1) << h1->GetRMSError() * 1000;
      rmserr = sstream_rmserr.str();
      //      TGaxis::SetMaxDigits(2)      ;
      gStyle->SetOptStat(0);
      TLatex Tl2;
      Tl2.SetTextAlign(12);
      Tl2.SetTextSize(0.04);
      Tl2.DrawLatexNDC(0.23,0.7," #mu = "+mean+" mm");
      Tl2.DrawLatexNDC(0.23,0.62," #sigma = "+rms+" #pm "+rmserr+" #mum");
      //Tl2.DrawLatexNDC(0.6,0.54,"|#Delta #eta |< 1.3");
      //Tl2.DrawLatexNDC(0.6,0.51,"m_{jj} > 1050 GeV");
      //Tl2.DrawLatexNDC(0.6,0.48,"65 GeV < M_{SD} < 105 GeV");
      

    }
  if(func == "RMS95" || func == "RMS99")
    {
      cout << "Getting RMS " << endl;
      *sigma = h1->GetRMS() * 1000;
      *sigmaerr = h1->GetRMSError() * 1000;
      cout << " resolution " << *sigma << " ± " << *sigmaerr << endl;
      *sigmaerr = poissonsmearing( h1, name,  detectorA, detectorB, detectorC);
      cout << " resolution updated " << *sigma << " ± " << *sigmaerr << endl;
      
      TString mean;
      std::ostringstream sstream_mean;
      sstream_mean << setprecision(2) << h1->GetMean();// * 1000;
      mean = sstream_mean.str();
      TString meanerr;
      std::ostringstream sstream_meanerr;
      sstream_meanerr << setprecision(1) << h1->GetMeanError();// * 1000;
      meanerr = sstream_meanerr.str();
      
      TString rms;
      std::ostringstream sstream_rms;
      sstream_rms << setprecision(3) << h1->GetRMS() * 1000;
      rms = sstream_rms.str();
      
      TString rmserr;
      std::ostringstream sstream_rmserr;
      sstream_rmserr << setprecision(1) << *sigmaerr;
      cout << sstream_rmserr.str() << endl;
      rmserr = sstream_rmserr.str();
      cout << rmserr << endl;
      //      TGaxis::SetMaxDigits(2)      ;
      gStyle->SetOptStat(0);
      TLatex Tl2;
      Tl2.SetTextAlign(12);
      Tl2.SetTextSize(0.04);
      Tl2.DrawLatexNDC(0.23,0.7," #mu = "+mean+" mm");
      Tl2.DrawLatexNDC(0.23,0.62," #sigma = "+rms+" #pm "+rmserr+" #mum");
      //Tl2.DrawLatexNDC(0.6,0.54,"|#Delta #eta |< 1.3");
      //Tl2.DrawLatexNDC(0.6,0.51,"m_{jj} > 1050 GeV");
      //Tl2.DrawLatexNDC(0.6,0.48,"65 GeV < M_{SD} < 105 GeV");
      cout << "legend done " << endl;

    }
  //  h1->GetXaxis()->SetLimits(-0.02,0.02);
  c->Update();
  cout << " canvas updated " << endl;
  gStyle->SetOptFit(1111);
  TString outputDir="/home/zoiirene/Output/Plots/";

  //  TString outputDir = outputDir+"";
  TString outputFile = outputDir+"residual_"+func+"Fit_"+name+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps";
  c->SaveAs(outputFile);
  cout << " saved output file " << endl;
  
}//







Double_t tp0Fit( Double_t *x, Double_t *par )
{
  double G = 1E0;

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


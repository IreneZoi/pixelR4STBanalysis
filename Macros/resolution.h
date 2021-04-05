

    }
  if(func == "RMS")
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


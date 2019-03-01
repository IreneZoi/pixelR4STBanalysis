//void FitTH1(TH1F* h1, Double_t *  sigma, Double_t *  sigmaerr, TString name, TString detectorA, TString detectorB, TString detectorC, TString func );
Double_t tp0Fit( Double_t *x, Double_t *par );


double G = 1E0;

void FitTH1(TH1F* h1, Double_t *  sigma, Double_t *  sigmaerr, TString name, TString detectorA, TString detectorB, TString detectorC, TString func)
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
  gPad->SetTicks(1,1);


  h1->SetTitle("");
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetLabelSize(0.025);
  h1->GetXaxis()->SetTitleSize(0.035);
  h1->GetXaxis()->SetTitleOffset(0.8);
  h1->GetXaxis()->SetTitleFont(42);

  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelSize(0.025);
  h1->GetYaxis()->SetTitleSize(0.035);
  h1->GetYaxis()->SetTitleOffset(1.4);
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

	  if(fabs(integral-integral95)/integral95 < 0.01)//integral>integral95)
	    break;
	 
	}
	  cout << "final integral " << integral << " low " << low << " high " << high << endl;


      cout << "final integral " << integral << " low " << low <<" high " << high << endl;
      h1->GetXaxis()->SetRange(low,high);
      *sigma = h1->GetRMS() * 1000;
      *sigmaerr = h1->GetRMSError() * 1000;
      cout << " resolution " << *sigma << " ± " << *sigmaerr << endl;
      gStyle->SetOptStat(1111);
    }
  //  h1->GetXaxis()->SetLimits(-0.02,0.02);
  c->Update();
  gStyle->SetOptFit(1111);
  TString outputDir="/home/zoiirene/Output/Plots/";

  //  TString outputDir = outputDir+"";
  TString outputFile = outputDir+"residual_"+func+"Fit_"+name+"_A_"+detectorA+"_B_"+detectorB+"_C_"+detectorC+".eps";
  c->SaveAs(outputFile);
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


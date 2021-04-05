
// Daniel Pitzl, Dec 2011
// fit Gauss + p0
// .L fitgp0.C+
// .x fitgp0.C+("h012")
// .x fitgp0.C+("h012;cycle")
   
#include "TDirectory.h"
#include "TStyle.h"
#include "TH1.h"
#include "TMath.h" // Gamma
#include "TF1.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVector.h"
#include <iomanip> // setw

// global:

const int mdata = 256;
double xx[mdata];
double yy[mdata];
int Ndata = 0;

//----------------------------------------------------------------------
// Gauss function:

double gp0Fit( double *x, double *par )
{
  static int nn = 0;
  nn++;
  static double dx = 0.1;
  static double b1 = 0;
  if( nn == 1 ) b1 = x[0];
  if( nn == 2 ) dx = x[0] - b1;

  // Mean and width:

  double t = ( x[0] - par[0] ) / par[1];

  double pi = 3.14159265358979323846;
  double aa = dx / par[1] / sqrt(2*pi);

  double pk = par[2] * aa * exp( -0.5*t*t );

  return pk + par[3];
}

//------------------------------------------------------------------------------
// Minuit expects:

void FCN( Int_t &npar, double *grd, double &f, double *par, Int_t flag )
{
  bool ldb = 0;
  if( flag == 5 ) ldb = 1;
  if( ldb ) cout << endl << "FCN called with flag " << flag << endl;

  // select case using flag
  switch( flag ) {
  case 1: // Initialization
    cout << "FCN case 1: init" << endl;
    break;
  case 2: // Compute derivatives, store them in grd
    cout << "FCN case 2: derive" << endl;
    break;
  case 3: // after the fit is finished
    cout << "FCN case 3: done" << endl;
    break;
  default: // compute FCN

    double chisq = 0;

    if( ldb ) {
      cout << "par";
      for( int ip = 0; ip < npar; ++ip ) cout << "  " << par[ip];
    }

    for( int ii = 0; ii < Ndata; ++ii ) {

      double x = xx[ii];
      double y = yy[ii];
      double e = 1;
      if( y > 0.5 ) e = sqrt(y);

      double r = y - gp0Fit( &x, par ); // resid = data - fit
      double c = r/e;
      chisq += c*c;

    } // bins ii

    if( ldb ) cout << ", chisq " << chisq << endl;

    f = chisq;

  } // flag switch

} // FCN

//----------------------------------------------------------------------

void fitgp0( string hs, double x1 = 1, double x9 = 0 )
{
  TH1 *h = (TH1*)gDirectory->Get(hs.c_str());

  if( h == NULL ) {
    cout << hs << " does not exist\n";
    return;
  }
   
  h->SetMarkerStyle(21);
  h->SetMarkerSize(0.8);
  h->SetStats(1);
  gStyle->SetOptFit(101);

  //gROOT->ForceStyle();

  int nb = h->GetNbinsX();

  if( x9 < x1 ) { // no user input
    x1 = h->GetBinCenter(1);
    x9 = h->GetBinCenter(nb);
  }

  int i1 = h->FindBin(x1);
  int i9 = h->FindBin(x9);
  double n1 = h->GetBinContent(i1);
  double n9 = h->GetBinContent(i9);
  double bg = 0.5*(n1+n9);

  double dx = h->GetBinWidth(1);
  double nmax = h->GetBinContent( h->GetMaximumBin() );
  double xmax = h->GetBinCenter( h->GetMaximumBin() );

  double nn = 7*(nmax-bg);

  cout << "fit " << hs << " from " << x1 << " to " << x9 << endl;

  int jj = 0;
  for( int ii = i1; ii <= i9; ++ii ) {
    xx[jj] = h->GetBinCenter(ii);
    yy[jj] = h->GetBinContent(ii);
    jj++;
    if( jj == mdata ) {
      cout << "too many bins, increase mdata" << endl;
      break;
    }
  }
  Ndata = jj;

  // create a TF1 with the range from x1 to x9 and 4 parameters

  const int mpar = 4;

  TF1* gp0Fcn = new TF1( "gp0Fcn", gp0Fit, x1, x9, mpar );

  gp0Fcn->SetParName( 0, "mean" );
  gp0Fcn->SetParName( 1, "sigma" );
  gp0Fcn->SetParName( 2, "area" );
  gp0Fcn->SetParName( 3, "BG" );

  // set start values for some parameters:

  cout << "start dx " << dx
       << ", max " << nmax
       << " at " << xmax
       << ", bg " << bg
       << ", area " << nn
       << endl;

  gp0Fcn->SetParameter( 0, xmax ); // peak position
  gp0Fcn->SetParameter( 1, 4*dx ); // width
  gp0Fcn->SetParameter( 2, nn ); // N
  gp0Fcn->SetParameter( 3, bg );
  
  //gp0Fcn->FixParameter( 3, 0 );  

  gp0Fcn->SetNpx(500);
  gp0Fcn->SetLineWidth(4);
  gp0Fcn->SetLineColor(kMagenta);
  gp0Fcn->SetLineColor(kGreen);

  cout << endl << "Migrad: " << endl << endl;
  h->Fit( "gp0Fcn", "R", "ep" );
  // h->Fit( "gp0Fcn","V+","ep");

  cout << endl << "Improve:" << endl << endl;
  h->Fit( "gp0Fcn", "RM", "ep" );

  cout << endl << "Minos:" << endl << endl;
  h->Fit( "gp0Fcn", "RME", "ep" );
  h->GetXaxis()->SetTitle("residual [mm]");
  h->GetYaxis()->SetTitle("Tracks");
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.);
  h->GetYaxis()->SetTitleOffset(1.);
  h->GetXaxis()->SetTitleOffset(0.9);
  h->SetTitle("");
  h->Draw("histepsame");  // data again on top

  int npar = mpar;
  double grd[mpar];
  double chisq0;
  double par[mpar];
  for( int j = 0; j < mpar; ++j )
    par[j] = gp0Fcn->GetParameter( j );
  int flag = 5;
  FCN( npar, grd, chisq0, par, flag );
  cout << "FCN chisq " << chisq0 << endl;

  // step size:

  double stp[mpar];
  for( int j = 0; j < mpar; ++j )
    stp[j] = gp0Fcn->GetParError( j );

  // scan FCN around minimum:
  // chisq = chisq0 + ( (p-p0)/s )^2
  // => dchi = (dp/s)^2
  // first iter: dchi1 for dp1
  // we want dchi2 = 1
  // sqrt(dchi2/dchi1) = dp2/dp1
  // => dp2 = dp1 / sqrt(dchi1)

  double xar[mpar];
  for( int j = 0; j < mpar; ++j ) xar[j] = par[j];
  double dp1[mpar];

  flag = 4; // get chisq
  double chisq;

  for( int j = 0; j < mpar; ++j ) {

    double dp = stp[j];
    int iter = 0;
    bool again = 1;
    double dchi = 0;

    do {
      iter++;
      xar[j] = par[j] + dp;
      FCN( npar, grd, chisq, xar, flag );
      dchi = chisq - chisq0;
      cout << "par " << j << " iter " << iter << "  " << xar[j] << " dchi2 " << dchi << endl;
      if( dchi > 2 )
	dp = dp/sqrt(dchi);
      else if( dchi < 0.01 )
	dp = 2*dp;
      else if( dchi < 0.5 )
	dp = dp/sqrt(dchi);
      else
	again = 0;
      if( iter > 9 ) again = 0;
    }
    while( again );

    cout << "par " << j
	 << "  " << par[j]
	 << ", iter " << iter
	 << ", step " << dp
	 << ", dchi " << dchi
	 << endl;

    xar[j] = par[j]; // back to minimum
    dp1[j] = dp;

  } // j par

  // Hessian for 2nd order derivatives:
  // f''(x) = ( f(x + h) − 2f(x) + f(x − h) ) / h^2
  // d2f/dxdy =
  // ( f(a+h1, b+h2) - f(a+h1, b-h2) - f(a-h1, b+h2) + f(a-h1, b-h2) / (4 h1 h2)

  double H[mpar][mpar];

  for( int j = 0; j < mpar; ++j ) {
    double dpj = dp1[j];
    xar[j] = par[j] + dpj;
    FCN( npar, grd, chisq, xar, flag );
    double dchiup = chisq - chisq0;
    xar[j] = par[j] - dpj;
    FCN( npar, grd, chisq, xar, flag );
    double dchidn = chisq - chisq0;
    double f2nd = ( dchiup + dchidn ) / (dpj*dpj);
    H[j][j] = f2nd;
    for( int k = j+1; k < mpar; ++k ){
      double dpk = dp1[k];
      xar[j] = par[j] + dpj;
      xar[k] = par[k] + dpk;
      FCN( npar, grd, chisq, xar, flag );
      double lupup = chisq;

      xar[k] = par[k] - dpk;
      FCN( npar, grd, chisq, xar, flag );
      double lupdn = chisq;

      xar[j] = par[j] - dpj;
      FCN( npar, grd, chisq, xar, flag );
      double ldndn = chisq;

      xar[k] = par[k] + dpk;
      FCN( npar, grd, chisq, xar, flag );
      double ldnup = chisq;

      double df2didj = ( lupup - lupdn - ldnup + ldndn ) / ( 4*dpj*dpk);
      H[j][k] = df2didj;
      H[k][j] = df2didj;

      xar[k] = par[k]; // back to minimum

    } // k
    xar[j] = par[j]; // back to minimum
  } // j

  // invert H

  TMatrixD hesse( mpar, mpar );

  for( int j = 0; j < mpar; ++j )
    for( int k = 0; k < mpar; ++k )
      TMatrixDRow( hesse, j )(k) = 0.5*H[j][k]; // justify factor 1/2

  TDecompSVD svd( mpar, mpar );
  svd.SetMatrix( hesse );
  svd.Decompose();

  TVectorD eigen = svd.GetSig();
  cout << "Eigenvalues";
  for( int j = 0; j < mpar; ++j )
    cout << "  " << eigen(j);
  cout << endl;
  cout << "condition number  " << svd.Condition() << endl;

  double d1, d2;
  svd.Det( d1, d2 );
  cout << "det  " << d1 * pow( 2, d2 ) << endl;

  TMatrixD covar( mpar, mpar );
  covar = svd.Invert();
  //covar = hesse.Invert();

  double sigma[mpar];
  for( int j = 0; j < mpar; ++j ) {
    sigma[j] = sqrt( fabs( TMatrixDRow( covar, j )(j) ) );
    cout << "par " << j
	 << ":  " << par[j]
	 << " +- " << sigma[j]
	 << endl;
  }

  for( int k = 0; k < mpar; ++k )
    cout << setw(14) << k;
  cout << endl;

  for( int j = 0; j < mpar; ++j ) {
    cout << setw(9) << j << "  ";
    cout << setw(j*14) << "";
    for( int k = j; k < mpar; ++k ) {
      int iw = 14;
      if( k == j ) iw = 3;
      cout << setw(iw)
	   << TMatrixDRow( covar, j )(k) / sigma[j] / sigma[k];
    }
    cout << endl;
  }

}

#include <iostream>
#include <fstream>
#include <cmath>

#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/WrappedTF1.h"
#include "Math/WrappedParamFunction.h"

#include "LHAPDF/LHAPDF.h"

using namespace std;

LHAPDF::PDF * xf = LHAPDF::mkPDF("NNPDFpol11_100", 0);

double integrand0(const double * var, const double * par){
  double x = var[0];
  double Q = par[0];
  double result = xf->xfxQ(2, x, Q) - xf->xfxQ(-2, x, Q) - xf->xfxQ(1, x, Q) + xf->xfxQ(-1, x, Q);
  return result / x;
}

double integrand1(const double * var, const double * par){
  double x = var[0];
  double Q = par[0];
  double result = xf->xfxQ(2, x, Q) - xf->xfxQ(1, x, Q);
  return result / x;
}

double integrand2(const double * var, const double * par){
  double x = var[0];
  double Q = par[0];
  double result = xf->xfxQ(1, x, Q) + xf->xfxQ(-1, x, Q);
  return result / x;
}

double integrand3(const double * var, const double * par){
  double x = var[0];
  double Q = par[0];
  double result = xf->xfxQ(2, x, Q) + xf->xfxQ(-2, x, Q) - xf->xfxQ(1, x, Q) - xf->xfxQ(-1, x, Q);
  return result / x;
}
  
int main(const int argc, const char * argv[]){

  LHAPDF::setVerbosity(0);

  double Q = atof(argv[1]);

  TF1 f0("integrand0", &integrand0, 0.0, 1.0, 1);
  TF1 f1("integrand1", &integrand1, 0.0, 1.0, 1);
  TF1 f2("integrand2", &integrand2, 0.0, 1.0, 1);
  TF1 f3("integrand3", &integrand3, 0.0, 1.0, 1);

  ROOT::Math::WrappedTF1 wf0(f0);
  ROOT::Math::WrappedTF1 wf1(f1);
  ROOT::Math::WrappedTF1 wf2(f2);
  ROOT::Math::WrappedTF1 wf3(f3);

  wf0.SetParameters(&Q);
  wf1.SetParameters(&Q);
  wf2.SetParameters(&Q);
  wf3.SetParameters(&Q);

  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0, 1e-5, 10000);
  ig.SetFunction(wf3);
  double central = 0;
  double sum = 0;

  central = ig.Integral(1e-6, 1.0);

  cout << "central = " << central << endl;

  //return 0;
  
  for (int i = 1; i <= 100; i++){
    //cout  <<  i << endl;
    xf = LHAPDF::mkPDF("NNPDFpol11_100", i);
    //cout << ig.Integral(1e-5, 1.0) << endl;
    sum += pow(ig.Integral(1e-6, 1.0) - central, 2);
  }

  sum = sqrt(sum / 100);

  cout << "error = " << sum << endl;

      
  
  return 0;
}

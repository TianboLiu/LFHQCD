#include <iostream>
#include <fstream>
#include <cmath>

#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/WrappedTF1.h"
#include "Math/WrappedParamFunction.h"

#include "LHAPDF/LHAPDF.h"

using namespace std;

const LHAPDF::PDF * xf = LHAPDF::mkPDF("NNPDFpol11_100", 0);

double integrand0(const double * var, const double * par){
  double x = var[0];
  double Q = par[0];
  double result = xf->xfxQ(2, x, Q) + xf->xfxQ(-2, x, Q) + xf->xfxQ(1, x, Q) + xf->xfxQ(-1, x, Q) + xf->xfxQ(3, x, Q) + xf->xfxQ(-3, x, Q) + xf->xfxQ(4, x, Q) + xf->xfxQ(-4, x, Q) + xf->xfxQ(5, x, Q) + xf->xfxQ(-5, x, Q);
  return result / x;
}

double integrand1(const double * var, const double * par){
  double x = var[0];
  double Q = par[0];
  double result = xf->xfxQ(2, x, Q) + xf->xfxQ(-2, x, Q);
  return result / x;
}

double integrand2(const double * var, const double * par){
  double x = var[0];
  double Q = par[0];
  double result = xf->xfxQ(1, x, Q) + xf->xfxQ(-1, x, Q);
  return result / x;
}
  
  
int main(const int argc, const char * argv[]){

  double Q = atof(argv[1]);

  TF1 f0("integrand0", &integrand0, 0.0, 1.0, 1);
  TF1 f1("integrand1", &integrand1, 0.0, 1.0, 1);
  TF1 f2("integrand2", &integrand2, 0.0, 1.0, 1);

  ROOT::Math::WrappedTF1 wf0(f0);
  ROOT::Math::WrappedTF1 wf1(f1);
  ROOT::Math::WrappedTF1 wf2(f2);

  wf0.SetParameters(&Q);
  wf1.SetParameters(&Q);
  wf2.SetParameters(&Q);

  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0, 1e-6, 10000);
  ig.SetFunction(wf0);
  double r0 = ig.Integral(1e-7, 1);
  ig.SetFunction(wf1);
  double r1 = ig.Integral(1e-5, 1);
  ig.SetFunction(wf2);
  double r2 = ig.Integral(1e-5, 1);

  cout << "Delta Sigma = " << r0 << endl;
  cout << "Delta  u+   = " << r1 << endl;
  cout << "Delta  d+   = " << r2 << endl;
  cout << "Delta u - d = " << r1 - r2 << endl;
  
  return 0;
}

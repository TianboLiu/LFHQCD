#include "hoppet_v1.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>

#include "gsl/gsl_sf_gamma.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using namespace std;

static const int AllNum = 13;
static const int HalfNum = 6;
static const double AlphasMz=0.1181; //from world average value
static const double Mz = 91.1876;

double wx(const double x, const double a){
  return pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2));
}

double dwx(const double x, const double a){
  return 2.0 * a * (1.0 - x) * pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2)) + pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2)) * ((1.0 - x) / x - log(x));
}


double s = 0;
double Q0 = 1.0;
void model0(const double & x, const double & Q, double * pdf){//proton poly
  double a0 = 0.531;
  double q3 = 2.0 * (1.0 - wx(x, a0)) * dwx(x, a0);
  double q4 = 3.0 * pow(1.0 - wx(x, a0), 2) * dwx(x, a0);
  double gA = 1.27;

  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar 
  pdf[ 0+HalfNum] = 0; //gluon
  pdf[ 1+HalfNum] = -x * gA * s * q4; //d
  pdf[ 2+HalfNum] = x * gA * (1.0 - s) * q3; //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

const double ymax = 12.0;//max y=ln(1/x) value we want to access
const double dy = 0.05;//internal ln(1/x) grid spacing, usually 0.1~0.25
const double Qmin = 1.0;//lower limit of Q we want to access
const double Qmax = 2.0;//upper limit of Q we want to access
const double dlnlnQ = dy/4.0;//internal spacing in lnlnQ, usually dy/4.0
const int order = -6;//order of numerical interpolation
const int nloop = 2;//max number of loops we want to use, i.e. LO, NLO, NNLO
const int factscheme = 3;//the scheme, 1: factscheme_MSbar; 2: factscheme_DIS; 3: factscheme_PolMSbar
const double muR_Q = 1.0;//ratio between muR and Q
const double mc = 1.28;//charm quark mass in MSbar scheme
const double mb = 4.18;//bottom quark mass in MSbar scheme
const double mt = 173.1;//top quark mass

double FitFunction(const double * par){
  //cout << "Calling..." << endl;
  s = par[0];
  Q0 = par[1];
  hoppetEvolve(AlphasMz, Mz, nloop, muR_Q, model0, Q0);
  double pdf[AllNum];
  int points = 1000;
  double x = 0;
  double sum = 0;
  for (int i = 0; i < points; i++){
    x = 0.5 / points + 1.0 / points * i;
    hoppetEval(x, 1.087, pdf);
    sum += pdf[2+HalfNum] + pdf[-2+HalfNum] - pdf[1+HalfNum] - pdf[-1+HalfNum];
  }
  sum = sum / points;
  double chi2 = pow((sum - 0.247) / 0.0136, 2);
  return chi2;
}
  

int main(const int argc, const char * argv[]){

  //starting hoppet
  hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,factscheme);
  //Set thresholds with quark MSbar masses
  hoppetSetMSbarMassVFN(mc, mb, mt);
  
  double par[2] = {0.3, 0.59};
  
  cout << FitFunction(par) <<  endl;
  //return 0;
  
  const int Npar = 2;
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1.0e-4);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&FitFunction, Npar);
  min->SetFunction(f);
  min->SetLimitedVariable(0, "s", 0.3, 1e-6, 0.0, 1.0);
  //min->SetFixedVariable(0, "s", 0.0);
  min->SetLimitedVariable(1, "Q0", 1.0, 1e-6, 0.5, 1.1);
  //min->SetFixedVariable(1, "Q0", 1.087);

  min->Minimize();
  
  return 0; 
}


//-------------------------------------------------------------------------

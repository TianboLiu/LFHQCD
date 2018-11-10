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
static const double AlphasMz = 0.1181;
static const double Mz = 91.1876;

double U1(const double x, const double a){
  double b = 3.0 - 2.0 * a;
  return a * x + b * pow(x, 2) + (1.0 - a - b) * pow(x, 3);
}

double dU1(const double x, const double a){
  double b = 3.0 - 2.0 * a;
  return a + 2.0 * b * x + 3.0 * (1.0 - a - b) * pow(x, 2);
}

double qx1(const double x, const double tau, const double a){
  //double lambda = pow(0.5482, 2);
  double result = gsl_sf_gamma(tau - 0.5) / (sqrt(M_PI) * gsl_sf_gamma(tau - 1.0)) * pow(U1(x, a), -0.5) * dU1(x, a) * pow(1.0 - U1(x, a), tau - 2.0);
  return result;
}

double U2(const double x, const double a){
  return pow(x, a * (1.0 - x));
}

double dU2(const double x, const double a){
  return a * pow(x, a * (1.0 - x)) * ((1.0 - x) / x - log(x));
}

double qx2(const double x, const double tau, const double a){
  //double lambda = pow(0.5482, 2);
  double result = gsl_sf_gamma(tau - 0.5) / (sqrt(M_PI) * gsl_sf_gamma(tau - 1.0)) * pow(U2(x, a), -0.5) * dU2(x, a) * pow(1.0 - U2(x, a), tau - 2.0);
  return result;
}

double U3(const double x, const double a){
  return pow(x, a * (1.0 - x * x));
}

double dU3(const double x, const double a){
  return a * pow(x, a * (1.0 - x * x)) * ((1.0 - x * x) / x - 2.0 * x * log(x));
}

double qx3(const double x, const double tau, const double a){
  //double lambda = pow(0.5482, 2);
  double result = gsl_sf_gamma(tau - 0.5) / (sqrt(M_PI) * gsl_sf_gamma(tau - 1.0)) * pow(U3(x, a), -0.5) * dU3(x, a) * pow(1.0 - U3(x, a), tau - 2.0);
  return result;
}

double U4(const double x, const double a){
  return pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2));
}

double dU4(const double x, const double a){
  return 2.0 * a * (1.0 - x) * pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2)) + pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2)) * ((1.0 - x) / x - log(x));
}

double qx4(const double x, const double tau, const double a){
  //double lambda = pow(0.5482, 2);
  double result = gsl_sf_gamma(tau - 0.5) / (sqrt(M_PI) * gsl_sf_gamma(tau - 1.0)) * pow(U4(x, a), -0.5) * dU4(x, a) * pow(1.0 - U4(x, a), tau - 2.0);
  return result;
}

double (* qx)(const double x, const double tau, const double a);

double a0 = 1.0;
void model(const double & x, const double & Q, double * pdf){
  double r = 2.0;
  double s = 1.0 / 3.0;
  double q3 = qx(x, 3.0, a0);
  double q4 = qx(x, 4.0, a0);
  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar 
  pdf[ 0+HalfNum] = 0; //gluon
  pdf[ 1+HalfNum] = x * (((1.0 + s) - 2.0 / 3.0 * r) * q3 - (s - 2.0 / 3.0 * r) * q4); //d
  pdf[ 2+HalfNum] = x * ((2.0 * (1.0 + s) - 1.0 / 3.0 * r) * q3 - (2.0 * s -  1.0 / 3.0 * r) * q4); //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

const double ymax = 12.0;
const double dy = 0.05;
const double Qmin = 2.0;
const double Qmax = 4.0;
const double dlnlnQ = dy / 4.0;
const int order  = -6;
const int nloop = 3;
const int factscheme = 1;
const double muR_Q = 1.0;
const double mc = 1.28;
const double mb = 4.18;
const double mt = 173.1;

//const double Q0 = 0.866;
const double Q0 = 1.057;
const double Q1 = 3.162;

double FitFunction(const double * par){
  a0 = par[0];
  hoppetEvolve(AlphasMz, Mz, nloop, muR_Q, model, Q0);
  double pdf[AllNum];
  double xu = 0;
  double xd = 0;
  int points = 1000;
  double x = 0;
  for (int i = 0; i < points; i++){
    x = 0.5 / points + 1.0 / points * i;
    hoppetEval(x, Q1, pdf);
    xu += pdf[2+HalfNum] - pdf[-2+HalfNum];
    xd += pdf[1+HalfNum] - pdf[-1+HalfNum];
  }
  xu = xu / points;
  xd = xd / points;
  double chi2 = pow((xu - 0.2594) / 0.0038, 2) + pow((xd - 0.1083) / 0.0045, 2);
  return chi2;
}

int main(const int argc, const char * argv[]){

  const int opt = atoi(argv[1]);
  if (opt == 1)
    qx = & qx1;
  else if (opt == 2)
    qx = & qx2;
  else if (opt == 3)
    qx = & qx3;
  else if (opt == 4)
    qx = & qx4;
  else
    return 0;

  hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, factscheme);
  hoppetSetMSbarMassVFN(mc, mb, mt);
  
  const int Npar = 1;
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1.0e-4);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&FitFunction, Npar);
  min->SetFunction(f);
  min->SetVariable(0, "a", 1.0, 1.0e-4);
  min->Minimize();

  return 0;
}

#include "hoppet_v1.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>

#include "gsl/gsl_sf_gamma.h"
//#include "gsl/gsl_sf_beta.h"

using namespace std;

static const int AllNum = 13;
static const int HalfNum = 6;

double a0 = 1, r0 = 1, gamma0 = 1;
// definition of our initial condition function
void (* model)(const double & x, const double & Q, double * pdf);
void model0(const double & x, const double & Q, double * pdf);
void model1(const double & x, const double & Q, double * pdf);
void model2(const double & x, const double & Q, double * pdf);

double Beta(const double a, const double b){
  return gsl_sf_beta(a, b);
}

//----------------------------------------------------------------------
int main(const int argc, const char * argv[]){

  cout << Beta(3, -0.09) << endl;

  return 0;
  
  if (argc < 2){
    cout << "./evolution <input>" << endl;
    return 0;
  }

  const double ymax = 12.0;//max y=ln(1/x) value we want to access
  const double dy = 0.05;//internal ln(1/x) grid spacing, usually 0.1~0.25
  const double Qmin = 0.5;//lower limit of Q we want to access
  const double Qmax = 12.0;//upper limit of Q we want to access
  const double dlnlnQ = dy/4.0;//internal spacing in lnlnQ, usually dy/4.0
  const int order = -6;//order of numerical interpolation
  const double muR_Q = 1.0;//ratio between muR and Q
  const double mc = 1.28;//charm quark mass in MSbar scheme
  const double mb = 4.18;//bottom quark mass in MSbar scheme
  const double mt = 173.1;//top quark mass

  //const double Q0 = atof(argv[4]);//initial scale of the input xf(x)

  int opt, nloop, factscheme;
  double Q02, Q2, AlphaS, Q_AlphaS;

  char tmp[200];

  ifstream input(argv[1]);

  input >> tmp >> nloop;
  input >> tmp >> Q02;
  input >> tmp >> Q2;
  input >> tmp >> AlphaS;
  input >> tmp >> Q_AlphaS;
  input >> tmp >> factscheme;
  input >> tmp >> opt;
  input >> tmp >> a0;
  input >> tmp >> r0;
  input >> tmp >> gamma0;
  
  cout << AlphaS << endl;

  double Q0 = sqrt(Q02);
  double Q = sqrt(Q2);
  
  //starting hoppet
  hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,factscheme);

  //Set thresholds with quark MSbar masses
  hoppetSetMSbarMassVFN(mc, mb, mt);

  //Set input PDF
  if (opt == 0)//model choices
    model = & model0;
  else if (opt == 1)
    model = & model1;
  else if (opt == 2)
    model = & model2;
  else
    return 0;

  //a0 = atof(argv[6]);//parameter
  //Evolve PDF
  hoppetEvolve(AlphaS, Q_AlphaS, nloop, muR_Q, model, Q0);
  
  double X1[500], X2[500];
  for (int i = 0; i < 500; i++){
    X1[i] = pow(10.0, -4.0 + 0.004 * i);
    X2[i] = 0.01 + (1.0 - 0.01) / 499 * i;
  }
  
  FILE * fs = fopen("xf.dat", "w");
  fprintf(fs, "model:%d\tQ0^2=%.2f\tQ^2=%.2f GeV^2\n", opt, Q02, Q2);
  fprintf(fs, "x  xuv  xdv  xu  xd  xs  xc  xb  xg  xub  xdb  xsb  xcb  xbb\n");

  double x;
  double pdf[AllNum];

  for (int i = 0; i < 500; i++){
    x = X1[i];
    hoppetEval(x, Q, pdf);
    fprintf(fs, "%.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E\n",
	    x, pdf[2+HalfNum] - pdf[-2+HalfNum], pdf[1+HalfNum] - pdf[-1+HalfNum],
	    pdf[2+HalfNum], pdf[1+HalfNum], pdf[3+HalfNum], pdf[4+HalfNum], pdf[5+HalfNum], pdf[0+HalfNum],
	    pdf[-2+HalfNum], pdf[-1+HalfNum], pdf[-3+HalfNum], pdf[-4+HalfNum], pdf[-5+HalfNum]);
  }

  for (int i = 0; i < 500; i++){
    x = X2[i];
    hoppetEval(x, Q, pdf);
    fprintf(fs, "%.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E\n",
	    x, pdf[2+HalfNum] - pdf[-2+HalfNum], pdf[1+HalfNum] - pdf[-1+HalfNum],
	    pdf[2+HalfNum], pdf[1+HalfNum], pdf[3+HalfNum], pdf[4+HalfNum], pdf[5+HalfNum], pdf[0+HalfNum],
	    pdf[-2+HalfNum], pdf[-1+HalfNum], pdf[-3+HalfNum], pdf[-4+HalfNum], pdf[-5+HalfNum]);
  }
 
  fclose(fs);

  return 0; 
}

double wx(const double x){
  double a = a0;
  return pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2));
}

double dwx(const double x){
  double a = a0;
  return 2.0 * a * (1.0 - x) * pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2)) + pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2)) * ((1.0 - x) / x - log(x));
}

double qx(const double x, const double tau){
  double alpha0 = 0.48400716838522223;
  double a = a0;
  double result = gsl_sf_gamma(tau - 0.5) / (sqrt(M_PI) * gsl_sf_gamma(tau - 1.0)) * pow(wx(x), -0.5) * dwx(x) * pow(1.0 - wx(x), tau - 2.0);
  return result;
}

void model0(const double & x, const double & Q, double * pdf){//proton
  double r = r0;
  double q3 = qx(x, 3.0);
  double q4 = qx(x, 4.0);
  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar 
  pdf[ 0+HalfNum] = 0; //gluon
  pdf[ 1+HalfNum] = x * ((1.0 - 2.0 / 3.0 * r) * q3 + 2.0 / 3.0 * r * q4); //d
  pdf[ 2+HalfNum] = x * ((2.0 - 1.0 / 3.0 * r) * q3 + 1.0 / 3.0 * r * q4); //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void  model1(const double & x, const double & Q, double * pdf){//pion
  double gamma = gamma0;
  double q2 = qx(x, 2.0);
  double q4 = qx(x, 4.0);
  
  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = x * ((1.0 - gamma) * q2 + gamma * q4); //dbar
  
  pdf[ 0+HalfNum] = 0; //gluon

  pdf[ 1+HalfNum] = 0; //d
  pdf[ 2+HalfNum] = x * ((1.0 - gamma) * q2 + gamma * q4); //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}
//-------------------------------------------------------------------------


void model2(const double & x, const double & Q, double * pdf){//proton
  double s = 1.0 / 3.0;
  double r = 2.0;
  double q3 = qx(x, 3.0);
  double q4 = qx(x, 4.0);
  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar 
  pdf[ 0+HalfNum] = 0; //gluon
  pdf[ 1+HalfNum] = x * (((1.0 + s) - 2.0 / 3.0 * r) * q3 - (s - 2.0 / 3.0 * r) * q4); //d
  pdf[ 2+HalfNum] = x * ((2.0 * (1.0 + s) - 1.0 / 3.0 * r) * q3 - (2.0 * s - 1.0 / 3.0 * r) * q4); //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

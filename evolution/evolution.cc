#include "hoppet_v1.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>

#include "gsl/gsl_sf_gamma.h"

using namespace std;

static const int AllNum = 13;
static const int HalfNum = 6;
static const double AlphasMz=0.1181; //from world average value
static const double Mz = 91.1876;


// definition of our initial condition function
void (* model)(const double & x, const double & Q, double * pdf);
void model0(const double & x, const double & Q, double * pdf);
void model1(const double & x, const double & Q, double * pdf);
void model2(const double & x, const double & Q, double * pdf);
void model3(const double & x, const double & Q, double * pdf);
void model4(const double & x, const double & Q, double * pdf);
void model5(const double & x, const double & Q, double * pdf);
void model6(const double & x, const double & Q, double * pdf);

//----------------------------------------------------------------------
int main(const int argc, const char * argv[]){
  
  if (argc < 6){
    cout << "./evolution <scheme> <nloop> <model> <Q0> <Q>" << endl;
    return 0;
  }

  const double ymax = 12.0;//max y=ln(1/x) value we want to access
  const double dy = 0.05;//internal ln(1/x) grid spacing, usually 0.1~0.25
  const double Qmin = 0.5;//lower limit of Q we want to access
  const double Qmax = 12.0;//upper limit of Q we want to access
  const double dlnlnQ = dy/4.0;//internal spacing in lnlnQ, usually dy/4.0
  const int order = -6;//order of numerical interpolation
  const int nloop = atoi(argv[2]);//max number of loops we want to use, i.e. LO, NLO, NNLO
  const int factscheme = atoi(argv[1]);//the scheme, 1: factscheme_MSbar; 2: factscheme_DIS; 3: factscheme_PolMSbar
  const double muR_Q = 1.0;//ratio between muR and Q
  const double mc = 1.28;//charm quark mass in MSbar scheme
  const double mb = 4.18;//bottom quark mass in MSbar scheme
  const double mt = 173.1;//top quark mass

  const double Q0 = atof(argv[4]);//initial scale of the input xf(x)

  //starting hoppet
  hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,factscheme);

  //Set thresholds with quark MSbar masses
  hoppetSetMSbarMassVFN(mc, mb, mt);


  //Set input PDF
  const int opt = atoi(argv[3]);
  if (opt == 0)//model choices
    model = & model0;
  else if (opt == 1)
    model = & model1;
  else if (opt == 2)
    model = & model2;
  else if (opt == 3)
    model = & model3;
  else if (opt == 4)
    model = & model4;
  else if (opt == 5)
    model = & model5;
  else if (opt == 6)
    model = & model6;
  else
    return 0;

  //Evolve PDF
  hoppetEvolve(AlphasMz, Mz, nloop, muR_Q, model, Q0);
  
  double pdf[AllNum];
  double x, Q;

  Q = atof(argv[5]);//scale to evaluate the xf(x)

  FILE * fs = fopen("xf.dat", "w");
  
  //print out the headlines
  printf("Q0 = %.1f GeV\t Q = %.1f GeV\n", Q0, Q);
  fprintf(fs, "Q0 = %.1f GeV\t Q = %.1f GeV\n", Q0, Q);
  fprintf(fs, "x\t x*uv\t x*dv\t x*(uv-dv)\n");

  int Xspots = 1000;//X grid size
  for (int i = 0; i < Xspots; i++){
    x = 1.0 / Xspots * (i + 0.5);
    hoppetEval(x, Q, pdf);
    fprintf(fs, "%.3E\t%.3E\t%.3E\t%.3E\n",
	   x, pdf[2+HalfNum] - pdf[-2+HalfNum], pdf[1+HalfNum] - pdf[-1+HalfNum],
	    pdf[2+HalfNum] - pdf[-2+HalfNum] - pdf[1+HalfNum] + pdf[-1+HalfNum]);
  }

  fclose(fs);

  return 0; 
}




//-------------------------------------------------------------------------

double qx(const double x, const double tau, const double norm = 1.0, const double m1 = 0.0, const double m2 = 0.0){
  if (x >= 1.0) return 0;
  double lambda = pow(0.5682, 2);
  double result = 1.0 / norm * gsl_sf_gamma(tau - 0.5) / (sqrt(M_PI) * gsl_sf_gamma(tau - 1.0)) * pow(x, -0.5) * pow(1.0 - x, tau - 2.0) * exp(-1.0 / lambda * (m1 * m1 / x + m2 * m2 / (1.0 - x)) * x * log(1.0 / x) / (1.0 - x));
  return result;
}

void  model0(const double & x, const double & Q, double * pdf){
  double r = 1.5;
  double gamma_p = 0.0;
  double gamma_n = 0.0;
  double beta_n = 0.5;
  double norm = 1.0;
  double m1 = 0.0;
  double m2 = 0.0;
  double q3 = qx(x, 3.0, norm, m1, m2);
  double q4 = qx(x, 4.0, norm, m1, m2);
  double q5 = qx(x, 5.0, norm, m1, m2);
  double q6 = qx(x, 6.0, norm, m1, m2);
  
  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar
  
  pdf[ 0+HalfNum] = 0; //gluon

  pdf[ 1+HalfNum] = x * (((1.0 - gamma_p) - 2.0 * r / 3.0 * (1.0 - gamma_n)) * q3 + 2.0 * r / 3.0 * (1.0 - beta_n) * q4 + (gamma_p - 2.0 * r / 3.0 * gamma_n) * q5 + 2.0 * r / 3.0 * beta_n * q6); //d
  pdf[ 2+HalfNum] = x * ((2.0 * (1.0 - gamma_p) - r / 3.0 * (1.0 - gamma_n)) * q3 + r / 3.0 * (1.0 - beta_n) * q4 + (2.0 * gamma_p - r / 3.0 * gamma_n) * q5 + r / 3.0 * beta_n * q6); //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void  model1(const double & x, const double & Q, double * pdf){
  double r = 1.5;
  double gamma_p = 0.0;
  double gamma_n = 0.0;
  double beta_n = 0.5;
  double m1 = 0.05;
  double m2 = 0.1;
  double q3 = qx(x, 3.0, 0.957034, m1, m2);
  double q4 = qx(x, 4.0, 0.962172, m1, m2);
  double q5 = qx(x, 5.0, 0.963430, m1, m2);
  double q6 = qx(x, 6.0, 0.963694, m1, m2);
  
  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar
  
  pdf[ 0+HalfNum] = 0; //gluon

  pdf[ 1+HalfNum] = x * (((1.0 - gamma_p) - 2.0 * r / 3.0 * (1.0 - gamma_n)) * q3 + 2.0 * r / 3.0 * (1.0 - beta_n) * q4 + (gamma_p - 2.0 * r / 3.0 * gamma_n) * q5 + 2.0 * r / 3.0 * beta_n * q6); //d
  pdf[ 2+HalfNum] = x * ((2.0 * (1.0 - gamma_p) - r / 3.0 * (1.0 - gamma_n)) * q3 + r / 3.0 * (1.0 - beta_n) * q4 + (2.0 * gamma_p - r / 3.0 * gamma_n) * q5 + r / 3.0 * beta_n * q6); //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void  model2(const double & x, const double & Q, double * pdf){
  double r = 1.5;
  double gamma_p = 0.0;
  double gamma_n = 0.0;
  double beta_n = 0.5;
  double m1 = 0.1;
  double m2 = 0.2;
  double q3 = qx(x, 3.0, 0.843896, m1, m2);
  double q4 = qx(x, 4.0, 0.858396, m1, m2);
  double q5 = qx(x, 5.0, 0.862385, m1, m2);
  double q6 = qx(x, 6.0, 0.863247, m1, m2);
  
  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar
  
  pdf[ 0+HalfNum] = 0; //gluon

  pdf[ 1+HalfNum] = x * (((1.0 - gamma_p) - 2.0 * r / 3.0 * (1.0 - gamma_n)) * q3 + 2.0 * r / 3.0 * (1.0 - beta_n) * q4 + (gamma_p - 2.0 * r / 3.0 * gamma_n) * q5 + 2.0 * r / 3.0 * beta_n * q6); //d
  pdf[ 2+HalfNum] = x * ((2.0 * (1.0 - gamma_p) - r / 3.0 * (1.0 - gamma_n)) * q3 + r / 3.0 * (1.0 - beta_n) * q4 + (2.0 * gamma_p - r / 3.0 * gamma_n) * q5 + r / 3.0 * beta_n * q6); //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void  model3(const double & x, const double & Q, double * pdf){
  double r = 1.5;
  double gamma_p = 0.0;
  double gamma_n = 0.0;
  double beta_n = 0.5;
  double m1 = 0.2;
  double m2 = 0.4;
  double q3 = qx(x, 3.0, 0.529504, m1, m2);
  double q4 = qx(x, 4.0, 0.552687, m1, m2);
  double q5 = qx(x, 5.0, 0.560209, m1, m2);
  double q6 = qx(x, 6.0, 0.561944, m1, m2);
  
  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar
  
  pdf[ 0+HalfNum] = 0; //gluon

  pdf[ 1+HalfNum] = x * (((1.0 - gamma_p) - 2.0 * r / 3.0 * (1.0 - gamma_n)) * q3 + 2.0 * r / 3.0 * (1.0 - beta_n) * q4 + (gamma_p - 2.0 * r / 3.0 * gamma_n) * q5 + 2.0 * r / 3.0 * beta_n * q6); //d
  pdf[ 2+HalfNum] = x * ((2.0 * (1.0 - gamma_p) - r / 3.0 * (1.0 - gamma_n)) * q3 + r / 3.0 * (1.0 - beta_n) * q4 + (2.0 * gamma_p - r / 3.0 * gamma_n) * q5 + r / 3.0 * beta_n * q6); //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void model4(const double & x, const double & Q, double * pdf){//pion
  double gamma = 0.125;
  
  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = x * pow(x, -0.5) * ((1.0 - gamma) / 2.0 + 15.0 / 16.0 * gamma * pow(1.0 - x, 2)); //dbar
  
  pdf[ 0+HalfNum] = 0; //gluon

  pdf[ 1+HalfNum] = 0; //d
  pdf[ 2+HalfNum] = x * pow(x, -0.5) * ((1.0 - gamma) / 2.0 + 15.0 / 16.0 * gamma * pow(1.0 - x, 2)); //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void  model5(const double & x, const double & Q, double * pdf){
  double gamma = 0.125;
  double eta = 1.0;
  double norm = (pow(4.0, - 2.0 - eta) * sqrt(M_PI) * (30.0 + eta * (32.0 * (2.0 + eta) - gamma * (19.0 + 17.0 * eta))) * gsl_sf_gamma(1.0 + 2.0 * eta)) / gsl_sf_gamma(3.5 + 2.0 * eta);
  
  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = x * pow(x, -0.5) * ((1.0 - gamma) / 2.0 + 15.0 / 16.0 * gamma * pow(1.0 - x, 2)) * pow(x * (1.0 - x), eta) / norm; //dbar
  
  pdf[ 0+HalfNum] = 0; //gluon

  pdf[ 1+HalfNum] = 0; //d
  pdf[ 2+HalfNum] = x * pow(x, -0.5) * ((1.0 - gamma) / 2.0 + 15.0 / 16.0 * gamma * pow(1.0 - x, 2)) * pow(x * (1.0 - x), eta) / norm; //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void  model6(const double & x, const double & Q, double * pdf){
  double alpha = 0.70;
  double beta = 2.03;
  double gamma = 13.8;
  double delta = 2.0;
  double N = 1.0 / 1.193;
  double xv = N * pow(x, alpha) * pow(1.0 - x, beta) * (1.0 + gamma * pow(x, delta));
  
  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = xv; //dbar
  
  pdf[ 0+HalfNum] = 0; //gluon

  pdf[ 1+HalfNum] = 0; //d
  pdf[ 2+HalfNum] = xv; //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}


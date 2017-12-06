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
void model7(const double & x, const double & Q, double * pdf);
void model8(const double & x, const double & Q, double * pdf);
void model9(const double & x, const double & Q, double * pdf);
void model10(const double & x, const double & Q, double * pdf);

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
  else if (opt == 7)
    model = & model7;
  else if (opt == 8)
    model = & model8;
  else if (opt == 9)
    model = & model9;
  else if (opt == 10)
    model = & model10;
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
    x = pow(10.0, 4.0 / Xspots * i - 4.0);
    hoppetEval(x, Q, pdf);
    fprintf(fs, "%.3E\t%.3E\t%.3E\t%.3E\n",
	   x, pdf[2+HalfNum] - pdf[-2+HalfNum], pdf[1+HalfNum] - pdf[-1+HalfNum],
	    pdf[2+HalfNum] - pdf[-2+HalfNum] - pdf[1+HalfNum] + pdf[-1+HalfNum]);
  }

  fclose(fs);

  return 0; 
}




//-------------------------------------------------------------------------

double qx0(const double x, const double tau){//massless
  double result = gsl_sf_gamma(tau - 0.5) / (sqrt(M_PI) * gsl_sf_gamma(tau - 1.0)) * pow(x, -0.5) * pow(1.0 - x, tau - 2.0);
  return result;
}

double qxmass(const double x, const double tau, const double m){//exponential mass correction
  double lambda = pow(0.5482, 2);
  double result = qx0(x, tau) * exp(-m * m / (lambda * x * (1.0 - x)));
  return result;//not normalized
}

double qxlongi(const double x, const double tau, const double m){//longitudinal mass correction
  double lambda = pow(0.5482, 2);
  double eta = 4.0 * m * m / lambda;
  double result =  qx0(x, tau) * pow(x * (1.0 - x), eta);
  return result;//not normalized
}

double qxalpha(const double x, const double tau, const double alpha){//modified beta-function representation
  double lambda = pow(0.5482, 2);
  double result = gsl_sf_gamma(tau - 0.5) / (sqrt(M_PI) * gsl_sf_gamma(tau - 1.0)) * pow(1.0 - pow(1.0 - x, alpha), -0.5) * pow(1.0 - x, alpha * tau - alpha - 1.0);
  return result;
}

double ux(const double x){//u(x)=x e^[x(1-x)]
  return x * exp(x * (1.0 - x));
}

double dux(const double x){//u'(x)
  return (1.0 + x - 2.0 * x * x) * exp(x * (1.0 - x));
}

double qxu(const double x, const double tau){//modified by u(x)
  double lambda = pow(0.5482, 2);
  double result = gsl_sf_gamma(tau - 0.5) / (sqrt(M_PI) * gsl_sf_gamma(tau - 1.0)) * pow(ux(x), -0.5) * dux(x) * pow(1.0 - ux(x), tau - 2.0);
  return result;
}

double ux1(const double x, const double a = 1.0){//u(x)=ax+bx^2+(1-a-b)x^3
  //0<a<3
  double b = 3.0 - 2.0 * a;
  return a * x + b * x * x + (1.0 - a - b) * pow(x, 3);
}

double dux1(const double x, const double a = 1.0){
  double b = 3.0 - 2.0 * a;
  return a + 2.0 * b * x + 3.0 * (1.0 - a - b) * x * x;
}

double qxu1(const double x, const double tau, const double a = 1.0){//modified by u(x)
  double lambda = pow(0.5482, 2);
  double result = gsl_sf_gamma(tau - 0.5) / (sqrt(M_PI) * gsl_sf_gamma(tau - 1.0)) * pow(ux1(x, a), -0.5) * dux1(x, a) * pow(1.0 - ux1(x, a), tau - 2.0);
  return result;
}

double ux2(const double x, const double a = 1.5){//u(x)=x^[a(1-x)]
  return pow(x, a * (1.0 - x));
}

double dux2(const double x, const double a = 1.5){//
  return a * pow(x, a * (1.0 - x)) * ((1.0 - x) / x - log(x));
}

double qxu2(const double x, const double tau, const double a = 1.5){//modified by u(x)
  double lambda = pow(0.5482, 2);
  double result = gsl_sf_gamma(tau - 0.5) / (sqrt(M_PI) * gsl_sf_gamma(tau - 1.0)) * pow(ux2(x, a), -0.5) * dux2(x, a) * pow(1.0 - ux2(x, a), tau - 2.0);
  return result;
}

void  model0(const double & x, const double & Q, double * pdf){//proton massless
  double r = 1.5;
  double q3 = qx0(x, 3.0);
  double q4 = qx0(x, 4.0);
  
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

void  model1(const double & x, const double & Q, double * pdf){//proton exponential mass
  double r = 1.5;
  double m = 0.050;
  double q3 = qxmass(x, 3.0, m) / 0.770764;
  double q4 = qxmass(x, 4.0, m) / 0.727216;
  
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

void  model2(const double & x, const double & Q, double * pdf){//proton longitudinal mass
  double r = 1.5;
  double m = 0.050;
  double q3 = qxlongi(x, 3.0, m) / 0.908328;
  double q4 = qxlongi(x, 4.0, m) / 0.899495;
  
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

void  model3(const double & x, const double & Q, double * pdf){//proton massless with alpha 1.5
  double r = 1.5;
  double q3 = qxalpha(x, 3.0, 1.5);
  double q4 = qxalpha(x, 4.0, 1.5);
  
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

void  model4(const double & x, const double & Q, double * pdf){//proton massless with alpha 2.0
  double r = 1.5;
  double q3 = qxalpha(x, 3.0, 2.0);
  double q4 = qxalpha(x, 4.0, 2.0);
  
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

void  model5(const double & x, const double & Q, double * pdf){//proton massless with u(x) modification
  double r = 1.5;
  double q3 = qxu(x, 3.0);
  double q4 = qxu(x, 4.0);
  
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


void  model6(const double & x, const double & Q, double * pdf){//proton massless with u(x) modification
  double r = 1.5;
  double q3 = qxu1(x, 3.0, 0.0);
  double q4 = qxu1(x, 4.0, 0.0);
  
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

void  model7(const double & x, const double & Q, double * pdf){//proton massless with u(x) modification
  double r = 1.5;
  double q3 = qxu1(x, 3.0, 0.5);
  double q4 = qxu1(x, 4.0, 0.5);
  
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

void  model8(const double & x, const double & Q, double * pdf){//pion massless with u(x) modification
  double gamma = 0.125;
  double q2 = qxu1(x, 2.0, 0.5);
  double q4 = qxu1(x, 4.0, 0.5);
  
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

void  model9(const double & x, const double & Q, double * pdf){//proton massless with u(x) modification
  double r = 1.5;
  double q3 = qxu2(x, 3.0, 1.5);
  double q4 = qxu2(x, 4.0, 1.5);
  
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

void  model10(const double & x, const double & Q, double * pdf){//pion massless with u(x) modification
  double gamma = 0.125;
  double q2 = qxu2(x, 2.0, 1.5);
  double q4 = qxu2(x, 4.0, 1.5);
  
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

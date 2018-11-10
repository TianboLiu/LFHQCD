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

double a0 = 1.0;
// definition of our initial condition function
void (* model)(const double & x, const double & Q, double * pdf);
void model0(const double & x, const double & Q, double * pdf);
void model1(const double & x, const double & Q, double * pdf);
void model2(const double & x, const double & Q, double * pdf);
void model3(const double & x, const double & Q, double * pdf);
void model4(const double & x, const double & Q, double * pdf);
void model5(const double & x, const double & Q, double * pdf);
void model6(const double & x, const double & Q, double * pdf);
//void model7(const double & x, const double & Q, double * pdf);
//void model8(const double & x, const double & Q, double * pdf);
//void model9(const double & x, const double & Q, double * pdf);
//void model10(const double & x, const double & Q, double * pdf);

//----------------------------------------------------------------------
int main(const int argc, const char * argv[]){
  
  if (argc < 7){
    cout << "./evolution <scheme> <nloop> <model> <Q0> <Q> <a>" << endl;
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
  //else if (opt == 7)
  //  model = & model7;
  //else if (opt == 8)
  //  model = & model8;
  //else if (opt == 9)
  //  model = & model9;
  //else if (opt == 10)
  //  model = & model10;
  else
    return 0;

  a0 = atof(argv[6]);//parameter
  //Evolve PDF
  hoppetEvolve(AlphasMz, Mz, nloop, muR_Q, model, Q0);
  
  double pdf[AllNum];
  double x, Q;

  Q = atof(argv[5]);//scale to evaluate the xf(x)

  FILE * fs = fopen("xf.dat", "w");
  
  //print out the headlines
  printf("Q0 = %.3f GeV\t Q = %.3f GeV\n", Q0, Q);
  fprintf(fs, "Q0 = %.3f GeV\t Q = %.3f GeV\n", Q0, Q);
  fprintf(fs, "x\t x*(uv-dv)\t x*(u+-d+)\t xu+\t xd+\t xu\t xd\t xubar\t xdbar\n");

  double X[1000];
  for (int i = 0; i < 500; i++){
    X[i] = pow(10.0, -4.0 + 0.004 * i);
    X[i+500] = 0.01 + (1.0 - 0.01) / 499 * i;
  }
  for (int i = 0; i < 1000; i++){
    x = X[i];
    hoppetEval(x, Q, pdf);
    fprintf(fs, "%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\n",
	    x,
	    (pdf[2+HalfNum] - pdf[-2+HalfNum]) - (pdf[1+HalfNum] - pdf[-1+HalfNum]),
	    (pdf[2+HalfNum] + pdf[-2+HalfNum]) - (pdf[1+HalfNum] + pdf[-1+HalfNum]),
	    pdf[2+HalfNum] + pdf[-2+HalfNum],
	    pdf[1+HalfNum] + pdf[-1+HalfNum],
	    pdf[2+HalfNum],
	    pdf[1+HalfNum],
	    pdf[-2+HalfNum],
	    pdf[-1+HalfNum]);
  }


  fclose(fs);

  return 0; 
}

double wx(const double x, const double a){
  return pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2));
}

double dwx(const double x, const double a){
  return 2.0 * a * (1.0 - x) * pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2)) + pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2)) * ((1.0 - x) / x - log(x));
}

void model0(const double & x, const double & Q, double * pdf){//proton poly
  double qv = 2.0 * (1.0 - wx(x, a0)) * dwx(x, a0);

  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar 
  pdf[ 0+HalfNum] = 0; //gluon
  pdf[ 1+HalfNum] = 0; //d
  pdf[ 2+HalfNum] = x * qv; //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void model1(const double & x, const double & Q, double * pdf){//proton poly
  double q3 = 2.0 * (1.0 - wx(x, a0)) * dwx(x, a0);
  double q5 = 4.0 * pow(1.0 - wx(x, a0), 3) * dwx(x, a0);

  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar 
  pdf[ 0+HalfNum] = 0; //gluon
  pdf[ 1+HalfNum] = 0; //d
  pdf[ 2+HalfNum] = 0.7 * x * q3 + 0.3 * x * q5; //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void model2(const double & x, const double & Q, double * pdf){//proton poly
  double q3 = 2.0 * (1.0 - wx(x, a0)) * dwx(x, a0);
  double q5 = 4.0 * pow(1.0 - wx(x, a0), 3) * dwx(x, a0);

  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar 
  pdf[ 0+HalfNum] = 0; //gluon
  pdf[ 1+HalfNum] = 0; //d
  pdf[ 2+HalfNum] = 0.6 * x * q3 + 0.4 * x * q5; //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void model3(const double & x, const double & Q, double * pdf){//proton poly
  double q3 = 2.0 * (1.0 - wx(x, a0)) * dwx(x, a0);

  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar 
  pdf[ 0+HalfNum] = 0; //gluon
  pdf[ 1+HalfNum] = -1.0/3.0 * x * q3; //d
  pdf[ 2+HalfNum] = 2.0/3.0 * x * q3; //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void model4(const double & x, const double & Q, double * pdf){//proton poly
  double q3 = 2.0 * (1.0 - wx(x, 0.531)) * dwx(x, 0.531);
  double q4 = 3.0 * pow(1.0 - wx(x, 0.531), 2) * dwx(x, 0.531);
  double beta = a0;

  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar 
  pdf[ 0+HalfNum] = 0; //gluon
  pdf[ 1+HalfNum] = -beta * x * q4; //d
  pdf[ 2+HalfNum] = x * q3 - beta * x * q4; //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void model5(const double & x, const double & Q, double * pdf){//proton poly
  double q3 = 2.0 * (1.0 - wx(x, 0.531)) * dwx(x, 0.531);
  double q4 = 3.0 * pow(1.0 - wx(x, 0.531), 2) * dwx(x, 0.531);
  double s = a0;
  double p = 0;

  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = 0; //ubar
  pdf[-1+HalfNum] = 0; //dbar 
  pdf[ 0+HalfNum] = 0; //gluon
  pdf[ 1+HalfNum] = - x * (s + p) * q4; //d
  pdf[ 2+HalfNum] = x * ((1.0 - s) * q3 - p * q4); //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

void model6(const double & x, const double & Q, double * pdf){//proton poly
  double q3 = 2.0 * (1.0 - wx(x, a0)) * dwx(x, a0);
  double q4 = 3.0 * pow(1.0 - wx(x, a0), 2) * dwx(x, a0);
  double q5 = 4.0 * pow(1.0 - wx(x, a0), 3) * dwx(x, a0);
  double q6 = 5.0 * pow(1.0 - wx(x, a0), 4) * dwx(x, a0);
  double c4 = 0.25;
  double c5 = 0.03;
  double c6 = 0.06;
  double c3 = 1.0 - c4 - c5 - c6;

  pdf[-6+HalfNum] = 0; //tbar
  pdf[-5+HalfNum] = 0; //bbar
  pdf[-4+HalfNum] = 0; //cbar
  pdf[-3+HalfNum] = 0; //sbar
  pdf[-2+HalfNum] = x * c5 * q5; //ubar
  pdf[-1+HalfNum] = - x * c6 * q6; //dbar 
  pdf[ 0+HalfNum] = 0; //gluon
  pdf[ 1+HalfNum] = - x * c4 * q4; //d
  pdf[ 2+HalfNum] = x * c3 * q3; //u
  pdf[ 3+HalfNum] = 0; //s
  pdf[ 4+HalfNum] = 0; //c
  pdf[ 5+HalfNum] = 0; //b
  pdf[ 6+HalfNum] = 0; //t
}

//-------------------------------------------------------------------------

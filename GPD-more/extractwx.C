#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "Math/Functor.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/RootFinder.h"

using namespace std;

LHAPDF::PDF * xf;
LHAPDF::PDF * xpdf[101];

double Q0 = 1.057;
double q3(const double x){
  if (x >= 1.0) return 0;
  double xu = xf->xfxQ(2, x, Q0) - xf->xfxQ(-2, x, Q0);
  double xd = xf->xfxQ(1, x, Q0) - xf->xfxQ(-1, x, Q0);
  double result = (2.0 / 3.0 * xu - 1.0 / 3.0 * xd) / x;
  return result;
}

double xq3(const double x){
  if (x >= 1.0) return 0;
  double xu = xf->xfxQ(2, x, Q0) - xf->xfxQ(-2, x, Q0);
  double xd = xf->xfxQ(1, x, Q0) - xf->xfxQ(-1, x, Q0);
  double result = (2.0 / 3.0 * xu - 1.0 / 3.0 * xd);
  return result;
} 

double Qx(const double x){
  ROOT::Math::Functor1D wf(&q3);
  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-5, 1e-4, 2000);
  ig.SetFunction(wf);
  double result = ig.Integral(x, 1.0);
  return 2.0 * result - 2.0;
}

double DD;
double Equation(const double W){
  return pow(W, 3) - 3.0 * W - DD;
}

double solve_w(const double x){
  DD = Qx(x);
  ROOT::Math::Functor1D wf(&Equation);
  ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kBRENT);
  rf.SetFunction(wf, -1.0, 1.0);
  rf.Solve(1000, 1e-4, 1e-3);
  return pow(rf.Root(), 2);
}
     

int main(const int argc, const char * argv[]){

  if (argc < 4){
    cout << "./extractwx <pdfset> <NSET> <Q>" << endl;
    return 0;
  }

  const int NSET = atoi(argv[2]);
  Q0 = atof(argv[3]);

  for (int i = 0; i <= NSET; i++)
    xpdf[i] = LHAPDF::mkPDF(argv[1], i);

  double X[1000];
  for (int i = 0; i < 500; i++){
    X[i] = pow(10.0, -4.0 + 0.004 * i);
    X[i+500] = 0.01 + (1.0 - 0.01) / 499 * i;
  }
  double x;

  if (true){
    FILE * fs = fopen("wx.dat", "w");
    fprintf(fs, "%s\t%.3f GeV\n", argv[1], Q0);
    fprintf(fs, "x  w(x) delta_w\n");
    
    double w, dw;
    
    for (int i = 0; i < 1000; i++){
      x = X[i];
      xf = xpdf[0];
      w = solve_w(x);
      dw = 0.0;
      for (int j = 1; j <= NSET; j++){
	xf = xpdf[j];
	dw += pow(solve_w(x) - w, 2);
      }
      dw = sqrt(dw);
      fprintf(fs, "%.6E\t%.6E\t%.6E\n", x, w, dw);
    }
    
    fclose(fs);
  }

  if (true){
    FILE * fq = fopen("qx.dat", "w");
    fprintf(fq, "%s\t%.3f GeV\n", argv[1], Q0);
    fprintf(fq, "x\tq(x)\tdelta_q\txq(x)\tdelta_xq(x)\n");

    double q, dq, xq, dxq;
    for (int i = 0; i < 1000; i++){
      x = X[i];
      xf = xpdf[0];
      q = q3(x);
      xq = xq3(x);
      dq = 0;
      dxq = 0;
      for (int j = 1; j <= NSET; j++){
	xf = xpdf[j];
	dq += pow(q3(x) - q, 2);
	dxq += pow(xq3(x) - xq, 2);
      }
      dq = sqrt(dq);
      dxq = sqrt(dxq);
      fprintf(fq, "%.6E\t%.6E\t%.6E\t%.6E\t%.6E\n", x, q, dq, xq, dxq);
    }
    fclose(fq);
  }
    
  
  return 0;
}

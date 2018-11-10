#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "Math/GSLIntegrator.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/Interpolator.h"

using namespace std;

int points = 1000;

double xuv(const LHAPDF::PDF * xf);
double xdv(const LHAPDF::PDF * xf);

int main(const int argc, const char * argv[]){

  LHAPDF::PDF * xpdf;

  xpdf = LHAPDF::mkPDF("NNPDF30_nnlo_as_0118", 0);

  double Eu = xuv(xpdf);
  double Ed = xdv(xpdf);

  double Vu = 0;
  double Vd = 0;
  
  for (int i = 1; i <= 100; i++){
    xpdf = LHAPDF::mkPDF("NNPDF30_nnlo_as_0118", i);
    Vu += pow(xuv(xpdf) - Eu, 2);
    Vd += pow(xdv(xpdf) - Ed, 2);
  }
  Vu = Vu / 100.0;
  Vd = Vd / 100.0;

  cout << "xu:  " << Eu << "  " << sqrt(Vu) << endl;
  cout << "xd:  " << Ed << "  " << sqrt(Vd) << endl;

  return 0;
}

double Q2 = 1.057 * 1.057;
double xuv(const LHAPDF::PDF * xf){
  double step = 1.0 / points;
  double x0 = step / 2.0;
  double sum = 0;
  for (int i = 0; i < points; i++){
    sum += xf->xfxQ2(2, x0 + step * i, Q2) - xf->xfxQ2(-2, x0 + step * i, Q2);
  }
  sum *= step;
  return sum;
}

double xdv(const LHAPDF::PDF * xf){
  double step = 1.0 / points;
  double x0 = step / 2.0;
  double sum = 0;
  for (int i = 0; i < points; i++){
    sum += xf->xfxQ2(1, x0 + step * i, Q2) - xf->xfxQ2(-1, x0 + step * i, Q2);
  }
  sum *= step;
  return sum;
}
  

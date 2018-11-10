#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

using namespace std;

int main(const int argc, const char * argv[]){

  if (argc < 3){
    cout << "./genPDF <pdfset> <NEIG> <Q2>" << endl;
    return 0;
  }

  int NEIG = atoi(argv[2]);
  
  const LHAPDF::PDF * xf0 = LHAPDF::mkPDF(argv[1], 0);
  const LHAPDF::PDF * xf[100];
  for (int i = 0; i < NEIG; i++)
    xf[i] = LHAPDF::mkPDF(argv[1], i+1);

  double Q2 = atof(argv[3]);
  double X1[500], X2[500];
  for (int i = 0; i < 500; i++){
    X1[i] = pow(10.0, -4.0 + 0.004 * i);
    X2[i] = 0.01 + (1.0 - 0.01) / 499 * i;
  }

  FILE * fs = fopen("xfpara.dat", "w");
  fprintf(fs, "%s\t%.2f GeV^2\n", argv[1], Q2);
  fprintf(fs, "x  xuv  xdv  xu  xd  xs  xc  xb  xg  xub  xdb  xsb  xcb  xbb  Exuv  Exdv  Exu  Exd  Exs  Exc  Exb  Exg  Exub  Exdb  Exsb  Excb  Exbb\n");

  double x, xuv, xdv, xu, xd, xs, xc, xb, xg, xub, xdb, xsb, xcb, xbb;
  double Exuv, Exdv, Exu, Exd, Exs, Exc, Exb, Exg, Exub, Exdb, Exsb, Excb, Exbb;

  for (int i = 0; i < 500; i++){
    x = X1[i];
    xuv = xf0->xfxQ2(2, x, Q2) - xf0->xfxQ2(-2, x, Q2);
    xdv = xf0->xfxQ2(1, x, Q2) - xf0->xfxQ2(-1, x, Q2);
    xu = xf0->xfxQ2(2, x, Q2);
    xd = xf0->xfxQ2(1, x, Q2);
    xs = xf0->xfxQ2(3, x, Q2);
    xc = xf0->xfxQ2(4, x, Q2);
    xb = xf0->xfxQ2(5, x, Q2);
    xg = xf0->xfxQ2(21, x, Q2);
    xub = xf0->xfxQ2(-2, x, Q2);
    xdb = xf0->xfxQ2(-1, x, Q2);
    xsb = xf0->xfxQ2(-3, x, Q2);
    xcb = xf0->xfxQ2(-4, x, Q2);
    xbb = xf0->xfxQ2(-5, x, Q2);
    Exuv = 0; Exdv = 0; Exu = 0; Exd = 0; Exs = 0; Exc = 0; Exb = 0; Exg = 0; Exub = 0; Exdb = 0; Exsb = 0; Excb = 0; Exbb = 0;
    for (int j = 0; j < NEIG; j++){
      Exuv += pow(xf[j]->xfxQ2(2, x, Q2) - xf[j]->xfxQ2(-2, x, Q2) - xuv, 2) / 2.0;
      Exdv += pow(xf[j]->xfxQ2(1, x, Q2) - xf[j]->xfxQ2(-1, x, Q2) - xdv, 2) / 2.0;
      Exu += pow(xf[j]->xfxQ2(2, x, Q2) - xu, 2) / 2.0;
      Exd += pow(xf[j]->xfxQ2(1, x, Q2) - xd, 2) / 2.0;
      Exs += pow(xf[j]->xfxQ2(3, x, Q2) - xs, 2) / 2.0;
      Exc += pow(xf[j]->xfxQ2(4, x, Q2) - xc, 2) / 2.0;
      Exb += pow(xf[j]->xfxQ2(5, x, Q2) - xb, 2) / 2.0;
      Exg += pow(xf[j]->xfxQ2(21, x, Q2) - xg, 2) / 2.0;
      Exub += pow(xf[j]->xfxQ2(-2, x, Q2) - xub, 2) / 2.0;
      Exdb += pow(xf[j]->xfxQ2(-1, x, Q2) - xdb, 2) / 2.0;
      Exsb += pow(xf[j]->xfxQ2(-3, x, Q2) - xsb, 2) / 2.0;
      Excb += pow(xf[j]->xfxQ2(-4, x, Q2) - xcb, 2) / 2.0;
      Exbb += pow(xf[j]->xfxQ2(-5, x, Q2) - xbb, 2) / 2.0;
    }
    fprintf(fs, "%.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E\n",
	    x, xuv, xdv, xu, xd, xs, xc, xb, xg, xub, xdb, xsb, xcb, xbb,
	    sqrt(Exuv), sqrt(Exdv), sqrt(Exu), sqrt(Exd), sqrt(Exs), sqrt(Exc), sqrt(Exb), sqrt(Exg), sqrt(Exub), sqrt(Exdb), sqrt(Exsb), sqrt(Excb), sqrt(Exbb));
  }

  for (int i = 0; i < 500; i++){
    x = X2[i];
    xuv = xf0->xfxQ2(2, x, Q2) - xf0->xfxQ2(-2, x, Q2);
    xdv = xf0->xfxQ2(1, x, Q2) - xf0->xfxQ2(-1, x, Q2);
    xu = xf0->xfxQ2(2, x, Q2);
    xd = xf0->xfxQ2(1, x, Q2);
    xs = xf0->xfxQ2(3, x, Q2);
    xc = xf0->xfxQ2(4, x, Q2);
    xb = xf0->xfxQ2(5, x, Q2);
    xg = xf0->xfxQ2(21, x, Q2);
    xub = xf0->xfxQ2(-2, x, Q2);
    xdb = xf0->xfxQ2(-1, x, Q2);
    xsb = xf0->xfxQ2(-3, x, Q2);
    xcb = xf0->xfxQ2(-4, x, Q2);
    xbb = xf0->xfxQ2(-5, x, Q2);
    Exuv = 0; Exdv = 0; Exu = 0; Exd = 0; Exs = 0; Exc = 0; Exb = 0; Exg = 0; Exub = 0; Exdb = 0; Exsb = 0; Excb = 0; Exbb = 0;
    for (int j = 0; j < NEIG; j++){
      Exuv += pow(xf[j]->xfxQ2(2, x, Q2) - xf[j]->xfxQ2(-2, x, Q2) - xuv, 2) / 2.0;
      Exdv += pow(xf[j]->xfxQ2(1, x, Q2) - xf[j]->xfxQ2(-1, x, Q2) - xdv, 2) / 2.0;
      Exu += pow(xf[j]->xfxQ2(2, x, Q2) - xu, 2) / 2.0;
      Exd += pow(xf[j]->xfxQ2(1, x, Q2) - xd, 2) / 2.0;
      Exs += pow(xf[j]->xfxQ2(3, x, Q2) - xs, 2) / 2.0;
      Exc += pow(xf[j]->xfxQ2(4, x, Q2) - xc, 2) / 2.0;
      Exb += pow(xf[j]->xfxQ2(5, x, Q2) - xb, 2) / 2.0;
      Exg += pow(xf[j]->xfxQ2(21, x, Q2) - xg, 2) / 2.0;
      Exub += pow(xf[j]->xfxQ2(-2, x, Q2) - xub, 2) / 2.0;							     
      Exdb += pow(xf[j]->xfxQ2(-1, x, Q2) - xdb, 2) / 2.0;
      Exsb += pow(xf[j]->xfxQ2(-3, x, Q2) - xsb, 2) / 2.0;
      Excb += pow(xf[j]->xfxQ2(-4, x, Q2) - xcb, 2) / 2.0;
      Exbb += pow(xf[j]->xfxQ2(-5, x, Q2) - xbb, 2) / 2.0;
    }
    fprintf(fs, "%.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E\n",
	    x, xuv, xdv, xu, xd, xs, xc, xb, xg, xub, xdb, xsb, xcb, xbb,
	    sqrt(Exuv), sqrt(Exdv), sqrt(Exu), sqrt(Exd), sqrt(Exs), sqrt(Exc), sqrt(Exb), sqrt(Exg), sqrt(Exub), sqrt(Exdb), sqrt(Exsb), sqrt(Excb), sqrt(Exbb));
  }

  fclose(fs);
  
  return 0;
}

#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

using namespace std;

int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./NNPDFpol <Q>" << endl;
    return 0;
  }

  double Q = atof(argv[1]);

  LHAPDF::setVerbosity(0);
  
  LHAPDF::PDF * xf[101];
  for (int i = 0; i <= 100; i++)
    xf[i] = LHAPDF::mkPDF("NNPDFpol11_100", i);

  double X[1000];
  for (int i = 0; i < 500; i++){
    X[i] = pow(10.0, -4.0 + 0.004 * i);
    X[i+500] = 0.01 + (1.0 - 0.01) / 499 * i;
  }

  FILE * fs = fopen("fs.dat", "w");
  fprintf(fs, "NNPDFpol 1.1, Q = %.3f GeV\n", Q);
  fprintf(fs, "x\t x(u-d)\t error\t xu+\t error\t xd+\t error\t xs+\t error\t xu\t error\t xd\t error\t xs\t error\t xub\t error\t xdb\t error\t xsb\t error\n");

  double x;
  double central, diff;
  double uplus, dplus, splus, u, d, s, ubar, dbar, sbar;
  double Euplus, Edplus, Esplus, Eu, Ed, Es, Eubar, Edbar, Esbar;

  for (int j = 0; j < 1000; j++){
    x = X[j];
    u = xf[0]->xfxQ(2, x, Q);
    d = xf[0]->xfxQ(1, x, Q);
    s = xf[0]->xfxQ(3, x, Q);
    ubar = xf[0]->xfxQ(-2, x, Q);
    dbar = xf[0]->xfxQ(-1, x, Q);
    sbar = xf[0]->xfxQ(-3, x, Q);
    uplus = u + ubar;
    dplus = d + dbar;
    splus = s + sbar;
    central = uplus - dplus;
    Euplus = 0;
    Edplus = 0;
    Esplus = 0;
    Eu = 0;
    Ed = 0;
    Es = 0;
    Eubar = 0;
    Edbar = 0;
    Esbar = 0;
    diff = 0;
    for (int i = 1; i <= 100; i++){
      Euplus += pow( xf[i]->xfxQ(2, x, Q) + xf[i]->xfxQ(-2, x, Q) - uplus, 2);
      Edplus += pow( xf[i]->xfxQ(1, x, Q) + xf[i]->xfxQ(-1, x, Q) - dplus, 2);
      Esplus += pow( xf[i]->xfxQ(3, x, Q) + xf[i]->xfxQ(-3, x, Q) - splus, 2);
      Eu += pow(xf[i]->xfxQ(2, x, Q) - u, 2);
      Ed += pow(xf[i]->xfxQ(1, x, Q) - d, 2);
      Es += pow(xf[i]->xfxQ(3, x, Q) - s, 2);
      Eubar += pow(xf[i]->xfxQ(-2, x, Q) - ubar, 2);
      Edbar += pow(xf[i]->xfxQ(-1, x, Q) - dbar, 2);
      Esbar += pow(xf[i]->xfxQ(-3, x, Q) - sbar, 2);
      diff += pow(xf[i]->xfxQ(2, x, Q) + xf[i]->xfxQ(-2, x, Q) - xf[i]->xfxQ(1, x, Q) - xf[i]->xfxQ(-1, x, Q) - central, 2);
    }
    Euplus = sqrt(Euplus / 100);
    Edplus = sqrt(Edplus / 100);
    Esplus = sqrt(Esplus / 100);
    Eu = sqrt(Eu / 100);
    Ed = sqrt(Ed / 100);
    Es = sqrt(Es / 100);
    Eubar = sqrt(Eubar / 100);
    Edbar = sqrt(Edbar / 100);
    Esbar = sqrt(Esbar / 100);
    diff = sqrt(diff / 100);
    fprintf(fs, "%.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\t %.6E\n",
	    x, central, diff, uplus, Euplus, dplus, Edplus, splus, Esplus,
	    u, Eu, d, Ed, s, Es, ubar, Eubar, dbar, Edbar, s, Esbar);
  }

  fclose(fs);

  return 0;
}

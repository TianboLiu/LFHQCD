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

  LHAPDF::PDF * xf0 = LHAPDF::mkPDF("NNPDFpol11_100", 0);
  LHAPDF::PDF * xf[100];
  for (int i = 1; i <= 100; i++)
    xf[i] = LHAPDF::mkPDF("NNPDFpol11_100", i);

  double X[1000];
  for (int i = 0; i < 500; i++){
    X[i] = pow(10.0, -4.0 + 0.004 * i);
    X[i+500] = 0.01 + (1.0 - 0.01) / 499 * i;
  }

  FILE * fs = fopen("fs.dat", "w");
  fprintf(fs, "NNPDFpol 1.1, Q = %.3f GeV\n", Q);
  fprintf(fs, "x\t x(uv-dv)\t error\t x(u+-d+)\t error\n");

  double central1, diff1;
  double central2, diff2;

  for (int i = 0; i < 1000; i++){
    central1 = (xf0->xfxQ(2, X[i], Q) -  xf0->xfxQ(-2, X[i], Q)) - (xf0->xfxQ(1, X[i], Q) -  xf0->xfxQ(-1, X[i], Q));
    diff1 = 0;
    central2 = (xf0->xfxQ(2, X[i], Q) +  xf0->xfxQ(-2, X[i], Q)) - (xf0->xfxQ(1, X[i], Q) + xf0->xfxQ(-1, X[i], Q));
    diff2 = 0;
    for (int j = 1; j <= 100; j++){
      diff1 += pow( (xf[j]->xfxQ(2, X[i], Q) -  xf[j]->xfxQ(-2, X[i], Q)) - (xf[j]->xfxQ(1, X[i], Q) -  xf[j]->xfxQ(-1, X[i], Q)) - central1, 2);
      diff2 += pow( (xf[j]->xfxQ(2, X[i], Q) +  xf[j]->xfxQ(-2, X[i], Q)) - (xf[j]->xfxQ(1, X[i], Q) +  xf[j]->xfxQ(-1, X[i], Q)) - central1, 2);
    }
    diff1 = sqrt(diff1 / 100);
    diff2 = sqrt(diff2 / 100);
    fprintf(fs, "%.6E\t %.6E\t %.6E\t %.6E\t %.6E\n", X[i], central1, diff1, central2, diff2);
  }

  fclose(fs);

  return 0;
}

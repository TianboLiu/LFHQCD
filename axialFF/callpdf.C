#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

using namespace std;

int main(const int argc, const char * argv[]){

  if (argc < 4){
    cout << "./callpdf <unpol/pol> <flavor> <x> <Q2>" << endl;
    return 0;
  }
  
  LHAPDF::setVerbosity(0);
  LHAPDF::PDF * xf[101];
  
  if (strcmp(argv[1], "unpol") == 0){
    for (int i = 0; i <= 100; i++)
      xf[i] = LHAPDF::mkPDF("NNPDF30_nlo_as_0118", i);
  }
  else if (strcmp(argv[1], "pol") == 0){
    for (int i = 0; i <= 100; i++)
      xf[i] = LHAPDF::mkPDF("NNPDFpol11_100", i);
  }
  else
    return 0;

  double x = atof(argv[3]);
  double Q2 = atof(argv[4]);

  double val = 0;
  double err = 0;

  if (strcmp(argv[2], "u+") == 0){
    val = xf[0]->xfxQ2(2, x, Q2) + xf[0]->xfxQ2(-2, x, Q2);
    for (int i = 1; i <= 100; i++)
      err += pow( xf[i]->xfxQ2(2, x, Q2) + xf[i]->xfxQ2(-2, x, Q2) - val, 2);
    err = sqrt(err / 100);
  }
  else if (strcmp(argv[2], "d+") == 0){
    val = xf[0]->xfxQ2(1, x, Q2) + xf[0]->xfxQ2(-1, x, Q2);
    for (int i = 1; i <= 100; i++)
      err += pow( xf[i]->xfxQ2(1, x, Q2) + xf[i]->xfxQ2(-1, x, Q2) - val, 2);
    err = sqrt(err / 100);
  }
  else if (strcmp(argv[2], "uv") == 0){
    val = xf[0]->xfxQ2(2, x, Q2) - xf[0]->xfxQ2(-2, x, Q2);
    for (int i = 1; i <= 100; i++)
      err += pow( xf[i]->xfxQ2(2, x, Q2) - xf[i]->xfxQ2(-2, x, Q2) - val, 2);
    err = sqrt(err / 100);
  }
  else if (strcmp(argv[2], "dv") == 0){
    val = xf[0]->xfxQ2(1, x, Q2) - xf[0]->xfxQ2(-1, x, Q2);
    for (int i = 1; i <= 100; i++)
      err += pow( xf[i]->xfxQ2(1, x, Q2) - xf[i]->xfxQ2(-1, x, Q2) - val, 2);
    err = sqrt(err / 100);
  }
  else
    return 0;
  
  printf("x%s(%.3f,%.3f)\t%.6f\t%.6f\n",
	  argv[2], x, Q2,
	  val, err);
  
  return 0;
}

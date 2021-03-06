
#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/GaussIntegrator.h"

using namespace std;

const LHAPDF::PDF * CT14NNLO = LHAPDF::mkPDF("CT14nnlo", 0);
const LHAPDF::PDF * MMHT2014NNLO = LHAPDF::mkPDF("MMHT2014nnlo68cl", 0);
const LHAPDF::PDF * NNPDF30NNLO = LHAPDF::mkPDF("NNPDF30_nnlo_as_0118", 0);
const LHAPDF::PDF * NNPDF31NNLO = LHAPDF::mkPDF("NNPDF31_nnlo_as_0118", 0);

const LHAPDF::PDF * CT14NLO = LHAPDF::mkPDF("CT14nlo", 0);
const LHAPDF::PDF * MMHT2014NLO = LHAPDF::mkPDF("MMHT2014nlo68cl", 0);
const LHAPDF::PDF * NNPDF30NLO = LHAPDF::mkPDF("NNPDF30_nlo_as_0118", 0);
const LHAPDF::PDF * NNPDF31NLO = LHAPDF::mkPDF("NNPDF31_nlo_as_0118", 0);

const LHAPDF::PDF * pdf;

double xf(const double x, const double Q, const int flavor){
  if (flavor == 7)
    return pdf->xfxQ(2, x, Q) - pdf->xfxQ(-2, x, Q);
  if (flavor == 8)
    return pdf->xfxQ(1, x, Q) - pdf->xfxQ(-1, x, Q);
  return pdf->xfxQ(flavor, x, Q);
}

double Q0 = 1.0;
double xuv(double x){
  return xf(x, Q0, 7);
}

double xdv(double x){
  return xf(x, Q0, 8);
}

double M1_uv(){
  ROOT::Math::Functor1D wf(&xuv);
  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1e-6, 10000);
  ig.SetFunction(wf);
  return ig.Integral(0, 1);
}

double M1_dv(){
  ROOT::Math::Functor1D wf(&xdv);
  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1e-4, 10000);
  ig.SetFunction(wf);
  return ig.Integral(0, 1);
}

double _dbar_ubar(double x){
  return (xf(x, Q0, -1) - xf(x, Q0, -2)) / x;
}

double dbar_ubar(){
  ROOT::Math::Functor1D wf(&_dbar_ubar);
  ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1e-6, 10000);
  ig.SetFunction(wf);
  return ig.Integral(1e-3, 1);
}

int main(const int argc, const char * argv[]){

  const int opt = atoi(argv[1]);


  if (opt == 0){
    Q0 = 1.057;
    printf("Q0 = %.3f GeV\n", Q0); 
    printf("set             M1_uv\tM1_dv\n");
    pdf = CT14NNLO;
    printf("%s\t%.6f\t%.6f\n", "CT14NNLO", M1_uv(), M1_dv());
    pdf = MMHT2014NNLO;
    printf("%s\t%.6f\t%.6f\n", "MMHT2014NNLO", M1_uv(), M1_dv());
    pdf = NNPDF30NNLO;
    printf("%s\t%.6f\t%.6f\n", "NNPDF30NNLO", M1_uv(), M1_dv());
    pdf = NNPDF31NNLO;
    printf("%s\t%.6f\t%.6f\n", "NNPDF31NNLO", M1_uv(), M1_dv());
    
    Q0 = 1.087;
    printf("Q0 = %.3f GeV\n", Q0); 
    printf("set             M1_uv\tM1_dv\n");
    pdf = CT14NLO;
    printf("%s\t%.6f\t%.6f\n", "CT14NLO  ", M1_uv(), M1_dv());
    pdf = MMHT2014NLO;
    printf("%s\t%.6f\t%.6f\n", "MMHT2014NLO", M1_uv(), M1_dv());
    pdf = NNPDF30NLO;
    printf("%s\t%.6f\t%.6f\n", "NNPDF30NLO", M1_uv(), M1_dv());
    pdf = NNPDF31NLO;
    printf("%s\t%.6f\t%.6f\n", "NNPDF31NLO", M1_uv(), M1_dv());
  }

  if (opt == 1){
    Q0 = 1.057;
    pdf = NNPDF30NNLO;
    cout << dbar_ubar() << endl;
  }
  
  return 0;
}

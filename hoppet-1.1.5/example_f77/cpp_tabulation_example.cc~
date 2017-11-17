#include "../src/hoppet_v1.h"
#include<iostream>
#include<cmath>
#include<cstdio>

using namespace std;
static const int AllNum = 13;
static const int HalfNum = 6;
static const int Xspot=100;
static const double AlsMz=0.11784; //from world average value
static const double AlsMzError=0.0021; //from the same paper

static const double MassZ = 91.187;
static const double AlsMassZ = 0.125;

// definition of the initial condition function
// void  PoZXY_init(const double & x, const double & Q, double * pdf);
// void  chiralquark_init(const double & x, const double & Q, double * pdf) ;
// void  heralhc_init(const double & x, const double & Q, double * pdf) ;

// definition of our initial condition function
void  model1(const double & x, const double & Q, double * pdf);
void  model2(const double & x, const double & Q, double * pdf);

//------------------------------------------------
double GetAlphas(double Q, double AlphasMZ);
double GetAlphasGRV(double MU, double ALPSMZ);

//----------------------------------------------------------------------
int main()
{
  double dy    = 0.1;
  double ymax = 12.0;
  double Qmin = sqrt(0.3);
  double Qmax = 28000.0;
  double dlnlnQ = dy/4.0;
  int order = -6;
  int nloop = 2;//NNLO or NLO
  bool debug = true;//true for model 1, false for model 2
  bool evaluate = true;//true for AlsMassZ, false for GetAlphasGRV Q0
  // initialise with NNLO currently, VFN
  //hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,factscheme_PolMSbar);
  hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,factscheme_MSbar);

  // evolve the initial condition
  double Q0=sqrt(0.36);  //set initial Q2 
  double AsQ0, Q0Alphas;
  //AsQ0 = 0.35; Q0Alphas = Q0;   // set alphas initial value at Q0
  if(evaluate){
    AsQ0 = AlsMassZ; 
    Q0Alphas = MassZ;
  }// set alpha initial value at AsQ0=MassZ
  //double AsQ0 = GetAlphas(Q0, AlsMz); Q0Alphas = Q0;
  if(!evaluate){
    AsQ0 = GetAlphasGRV(Q0, AlsMz);
    Q0Alphas = Q0;
  }//Set Alphas initial value at Q0 with Alphas(Mz) value from World Average 
  cout<<AsQ0<<endl;
  //hoppetEvolve(AsQ0, Q0, nloop, 1.0, PoZXY_init, Q0);

  if (debug){
    hoppetEvolve(AsQ0, Q0Alphas, nloop, 1.0, model1, Q0);}
  
  if (!debug){
    hoppetEvolve(AsQ0, Q0Alphas, nloop, 1.0, model2, Q0);}
  // alternatively preprepare an evolution and then use its cached version.x
  //hoppetPreEvolve(AsQ0, Q0, nloop, 1.0, Q0);
  //hoppetCachedEvolve(heralhc_init);

  // output the results
  double pdf[AllNum];
  double xvals[Xspot];
  //xvals[Xspot]{1e-5,1e-4,1e-3,1e-2,0.1,0.3,0.5,0.7,0.9};
  for(int i=0;i<Xspot;i++){
    double series = i;
    xvals[i]=pow(10, 3 * (series/Xspot-1));
  }

  double Q = sqrt(4.0);     //set the point you would like to evolve to
 

  printf("           Evaluating PDFs at Q = %8.3f GeV\n",Q);

  
  FILE *pionpdf;
  if (debug){
    pionpdf = fopen("pionpdf.dat", "w");}

  if (!debug){
    pionpdf = fopen("pionpdf2.dat", "w");}
  
  for (int ix = 0; ix < Xspot; ix++){
    hoppetEval(xvals[ix], Q, pdf);
    fprintf(pionpdf, "%11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E\n",
	    xvals[ix], pdf[HalfNum+2], pdf[HalfNum-2], pdf[HalfNum+1], pdf[HalfNum-1], pdf[HalfNum+3], pdf[HalfNum-3], pdf[HalfNum+0]);
    //x u ubar d dbar s sbar gluon  
  }
  
  fclose(pionpdf);

 //  if(debug){
 // printf("    x         u         ubar        d           dbar         s          sbar        c          cbar         gluon         dbar/ubar  \n");
 //    for (int ix = 0; ix < Xspot; ix++) {
 //      hoppetEval(xvals[ix], Q, pdf);
 //      printf("%7.1E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E \n",xvals[ix],
 // 	     pdf[HalfNum+2],//u quark
 // 	     pdf[HalfNum-2],//ubar quark  
 // 	     pdf[HalfNum+1],//d quark
 // 	     pdf[HalfNum-1], //dbar quark
 // 	     pdf[HalfNum+3],//s quark
 // 	     pdf[HalfNum-3],//sbar quark
 // 	     pdf[HalfNum+4],//c quark
 // 	     pdf[HalfNum-4],//cbar quark
 // 	     pdf[HalfNum+0],//gluon  
 // 	     (pdf[HalfNum-1]/pdf[HalfNum-2]));
 //    }
 //  }

 
 //  if(evaluate){
 //    printf("    x      u-ubar       d-dbar   2(ubar+dbar)     s+sbar       c+cbar      gluon     s-sbar\n");
 //    for (int ix = 0; ix < Xspot; ix++) {
 //      hoppetEval(xvals[ix], Q, pdf);
 //      printf("%7.1E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E\n",xvals[ix],
 // 	     pdf[HalfNum+2]-pdf[HalfNum-2],//u-ubar quark  
 // 	     pdf[HalfNum+1]-pdf[HalfNum-1], //d-dbar quark
 // 	     2*(pdf[HalfNum-1]+pdf[HalfNum-2]),//2*(dbar+ubar)
 // 	     (pdf[HalfNum-3]+pdf[HalfNum+3]),//s+sbar
 // 	     (pdf[HalfNum-4]+pdf[HalfNum+4]),//c+cbar
 // 	     pdf[HalfNum+0],//gluon
 // 	     (pdf[HalfNum+3]-pdf[HalfNum-3]));//s-sbar
 //    }
 //  }
 //  printf("    x       s-sbar\n");
 //  for (int ix = 0; ix < Xspot; ix++) {
 //    hoppetEval(xvals[ix], Q, pdf);
 //    printf("%7.1E %11.4E\n",xvals[ix],
 // 	   (pdf[HalfNum+3]-pdf[HalfNum-3]));//s-sbar
 //  }
  
}

//----------------------------------------------------------------------
// the initial condition
double GetAlphas(double Q, double AlphasMZ){

  static const double TWOPI = 6.28318530717958647692528;
  static const double TWOPISQR = 39.47841760435743447533796;
  static const int NF = 5;
  static const double MZ = 91.1882;
  
  double BETA0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
  double BETA1 =  (51. - 19./3.*NF);
  
   // This is from NLOJET++, alpha.cc
  double res = AlphasMZ;
  double b0 = BETA0/TWOPI;
  double w = 1.0 + b0*AlphasMZ*log(Q/MZ);
  res /= w;
  double b1 = BETA1/TWOPISQR;
  res *= 1.0 - AlphasMZ*b1/b0*log(w)/w;
  return res;
}

double GetAlphasGRV(double MU, double ALPSMZ){
  // Implementation of Alpha_s evolution as function of Q.
 

  // c - initialize pi and beta functions
  const int NF	= 5; //it could be set to reset
  static const double TWOPI = 6.28318530717958647692528;
  const double PI4= TWOPI*2;
  const double B0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
  const double B1 =  (102. - 38./3.*NF);
  const double B10 = B1 / B0 / B0;
  const double ZMASS2 = pow ( 91.187 , 2 );

  // c - exact formula to extract Lambda from alpha_s(Mz)
  double Q2 = pow(MU,2);
  double LAM2 = ZMASS2 * exp( -1.*PI4/B0/ALPSMZ +  B10 * log( PI4/B0/ALPSMZ + B10) );

  // c - extract approx. alpha_s(mu) value 
  double LQ2 = log( Q2 / LAM2 ) ;
  double ASAPPROX = PI4/B0/LQ2 * (1. - B10*log(LQ2)/LQ2);
  double ALPHAS = ASAPPROX;
    
  // c - exact 2loop value by Newton procedure
  for ( int I = 1 ; I <=6 ; I++ ){
    double  F  = LQ2 - PI4/B0/ALPHAS + B10*log(PI4/B0/ALPHAS + B10);
    double FP = -1.*PI4/B0/(ALPHAS*1.01) + B10 * log(PI4/B0/(ALPHAS*1.01) + B10);
    double FM = -1.*PI4/B0/(ALPHAS*0.99) + B10 * log(PI4/B0/(ALPHAS*0.99) + B10);
    ALPHAS = ALPHAS - F/(FP-FM)*0.02*ALPHAS;
    // c      WRITE(*,*) ' LAMDA/a_s_approx/a_s = ',sqrt(lam2),ASAPPROX,ALPHAS
    // c        alpsmz = alpsmz + I*0.001
    // c       WRITE(*,*) ' alpsmz =', real(alpsmz)
  }
  
  printf("FNLO v1.4 as.  Q2=%7.4f  AlphasMZ_0=%7.4f, Alphas(Q)=%7.4f\n",Q2,ALPSMZ,ALPHAS);

  // c - that's it!
  return ALPHAS;

}

//-------------------------------------------------------------------------

void  model1(const double & x, const double & Q, double * pdf){
  double qv;
  double qsea;
  
  qv = 13.5604 * pow(x, 1.88002) * pow(1-x, 0.87976) * (1 - 2.37773 * pow (x * (1-x), 0.5) + 2.04689 * x * (1-x));
  qsea = 0./*558789 * pow(x, 1.0306) * pow(1-x, 1.2846) * (1 - 1.75194 * pow(x, 0.5) + 0.757313 * x)*/;

  pdf[ 0+HalfNum] = 0; //N_g * pow(x,0.8) * pow(1-x,3);//gluon 
  pdf[-3+HalfNum] = 0; //strange bar
  pdf[ 3+HalfNum] = 0; //strange
  pdf[ 2+HalfNum] = qv + qsea; //up
  pdf[-2+HalfNum] = qsea; //up bar
  pdf[ 1+HalfNum] = qsea; //down
  pdf[-1+HalfNum] = qv + qsea; //down bar

  pdf[ 4+HalfNum] = 0;// charm
  pdf[ 5+HalfNum] = 0;// bottom
  pdf[ 6+HalfNum] = 0;// top
  pdf[-4+HalfNum] = 0;//charm bar
  pdf[-5+HalfNum] = 0;//bottom bar
  pdf[-6+HalfNum] = 0;// top bar

}

void  model2(const double & x, const double & Q, double * pdf){
  double qv;
  double qsea;
  
  qv = 0.30983 * x * pow(1-x, 0.579534) * (2.80914 + 2.05078 * pow (x, 0.5) + 2.76892 * x);
  qsea = 0./*268644 * pow(x, 0.724103) * pow(1-x, 1.7175) * (1 - 1.8391 * pow(x, 0.5) + 0.835333 * x)*/;

  pdf[ 0+HalfNum] = 0; //N_g * pow(x,0.8) * pow(1-x,3);//gluon 
  pdf[-3+HalfNum] = 0; //strange bar
  pdf[ 3+HalfNum] = 0; //strange
  pdf[ 2+HalfNum] = qv + qsea; //up
  pdf[-2+HalfNum] = qsea; //up bar
  pdf[ 1+HalfNum] = qsea; //down
  pdf[-1+HalfNum] = qv + qsea; //down bar

  pdf[ 4+HalfNum] = 0;// charm
  pdf[ 5+HalfNum] = 0;// bottom
  pdf[ 6+HalfNum] = 0;// top
  pdf[-4+HalfNum] = 0;//charm bar
  pdf[-5+HalfNum] = 0;//bottom bar
  pdf[-6+HalfNum] = 0;// top bar

}




//--------------------------------------------------------------------------

// void  PoZXY_init(const double & x,
//                    const double & Q,
//                    double * pdf) {
//   double uv;
//   double dv;
//   double ubar; 
//   double dbar;
//   double N_ls=0.0; //dbar quark coefficience
//   double N_uv=13.6; //u quark coefficient
//   double N_dv=-3.37; //d quark coefficient
//   double N_db=N_ls/2; //dbar quark
//   double N_st=-2.4;//strange quark coefficience
//   //double N_g=0.0;   //set initial gluon coefficient
//   double N_g = 1.7;

//   uv = N_uv * pow(x,1.58) * pow((1-x),3.33); //set initial u quark
//   dv = N_dv * pow(x,1.53) * pow((1-x),2.58)* (1-1.99 * pow(x,0.5)+2.24 * x); // set initial d quark
//   dbar = N_db * pow(x,-0.1) * pow(1-x,6); // set initial dbar function on x
//   ubar = dbar * (1-x); //set initial ubar


//   pdf[ 0+HalfNum] = N_g * pow(x,0.8) * pow(1-x,3);//gluon 
//   pdf[-3+HalfNum] = 0; //sbar quark
//   pdf[ 3+HalfNum] =  N_st * pow(x,1.8) * pow(1-x,9.7) * (1 - 3.1 * pow(x, 0.5)+4.5 * x);//s quark
//   pdf[ 2+HalfNum] = uv + ubar; //u quark total
//   pdf[-2+HalfNum] = ubar; //u bar quark
//   pdf[ 1+HalfNum] = dv + dbar; //d quark total
//   pdf[-1+HalfNum] = dbar; //d bar quark total

//   pdf[ 4+HalfNum] = 0;// charm
//   pdf[ 5+HalfNum] = 0;// bottom
//   pdf[ 6+HalfNum] = 0;// top
//   pdf[-4+HalfNum] = 0;//charm bar
//   pdf[-5+HalfNum] = 0;//bottom bar
//   pdf[-6+HalfNum] = 0;// top bar
// }


// void  heralhc_init(const double & x,
//                    const double & Q,
//                    double * pdf) {
  
//   double uv;
//   double dv;
//   double ubar; 
//   double dbar;
  
//   double N_g = 1.7;
//   double N_ls = 0.387975;
//   double N_uv = 5.107200;
//   double N_dv = 3.064320;
//   double N_db = N_ls/2;
//   double N_st = 0.0;//strange quark coefficience

//   uv = N_uv * pow(x,0.8) * pow((1-x),3);
//   dv = N_dv * pow(x,0.8) * pow((1-x),4);
//   dbar = N_db * pow(x,-0.1) * pow(1-x,6);
//   ubar = dbar * (1-x);


//   pdf[ 0+HalfNum] = N_g * pow(x,0.8) * pow(1-x,3);//gluon 
//   pdf[-3+HalfNum] = 0.2 * (dbar+ ubar); //sbar quark
//   pdf[ 3+HalfNum] = N_st * pow(x,1.8) * pow(1-x,9.7) * (1 - 3.1 * pow(x, 0.5)+4.5 * x);//s quark
//   pdf[ 2+HalfNum] = uv + ubar; //u quark total
//   pdf[-2+HalfNum] = ubar; //u bar quark
//   pdf[ 1+HalfNum] = dv + dbar; //d quark total
//   pdf[-1+HalfNum] = dbar; //d bar quark total

//   pdf[ 4+HalfNum] = 0;// charm
//   pdf[ 5+HalfNum] = 0;// bottom
//   pdf[ 6+HalfNum] = 0;// top
//   pdf[-4+HalfNum] = 0;//charm bar
//   pdf[-5+HalfNum] = 0;//bottom bar
//   pdf[-6+HalfNum] = 0;// top bar
// }

// void  chiralquark_init(const double & x,
//                    const double & Q,
//                    double * pdf) {
  
//   double uu;
//   double dd;
//   double ubar; 
//   double dbar;
//   double st;
//   double stbar;
//   double N_g = 1.7;
//   double N_u = 1.66618;
//   double N_d = 0.634606;
//   double N_db = 0.878112;
//   double N_ub = 0.190835;
//   double N_stb = 0.567088;//strange quark coefficience
//   double N_st = 21.2255;

//   uu = N_u * pow(x,(1.0-0.413363)) * pow((1-x),1.74183) * (1 - 0.509391 * pow(x,0.5) + 3.92316 * x );
//   dd = N_d * pow(x,(1.0-0.459534)) * pow((1-x),6.32916) * (1 - 1.71424 * pow(x,0.5) + 1.29132 * x);
//   dbar = N_db * pow(x,(1.0-0.457605)) * pow((1-x),6.04274) * (1 - 1.11507 * pow(x,0.5) + 0.583017 * x);
//   ubar = N_ub *  pow(x,(1.0-0.613311)) * pow((1-x),1.95083) * (1 + 16.6968 * pow(x,0.5) - 7.19281 * x );
//   stbar =  N_stb *  pow(x,(1.0-0.462578)) * pow((1-x),6.35009) * (1 - 1.31012 * pow(x,0.5) + 0.816931 * x );
//   st = N_st * pow(x,(1.0+0.519606)) * pow((1-x),8.57117) * (1 - 2.93235 * pow(x,0.5) + 3.2472 * x ) ;
// //cout<<x<<"  "<<uu<<"  "<<ubar<<"  "<<dd<<"  "<<dbar<<"  "<<st<<"  "<<stbar<<endl;
//   pdf[ 0+HalfNum] = N_g * pow(x,0.8) * pow(1-x,3);//gluon 
//   pdf[-3+HalfNum] = stbar; //sbar quark
//   pdf[ 3+HalfNum] = st;//s quark
//   pdf[ 2+HalfNum] = uu; //u quark total
//   pdf[-2+HalfNum] = ubar; //u bar quark
//   pdf[ 1+HalfNum] = dd; //d quark total
//   pdf[-1+HalfNum] = dbar; //d bar quark total

//   pdf[ 4+HalfNum] = 0;// charm
//   pdf[ 5+HalfNum] = 0;// bottom
//   pdf[ 6+HalfNum] = 0;// top
//   pdf[-4+HalfNum] = 0;//charm bar
//   pdf[-5+HalfNum] = 0;//bottom bar
//   pdf[-6+HalfNum] = 0;// top bar
// }



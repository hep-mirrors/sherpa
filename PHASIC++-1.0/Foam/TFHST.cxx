///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//      Primitive Class of histograms                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include<cmath>
#include<iostream>
#include<fstream>
#include<iomanip>

#include"TFHST.h"

#define SW12 std::setprecision(7)<<std::setw(12)
#define SW8  std::setprecision(5)<<std::setw(8)



///////////////////////////////////////////////////////////////////////////////
TFHST::TFHST(){
// Empty histogram is legal!
  m_Entries=0.0;
  m_xmin=0.0;
  m_xmax=0.0;
  m_Bin1=NULL;
  m_Bin2=NULL;
  m_Nb  =0;
}
///////////////////////////////////////////////////////////////////////////////
TFHST::TFHST(double xmin, double xmax, int nb){
// PRINCIPAL Constructor for practical use
  m_Entries=0.0;
  int i;
  m_xmin=xmin;
  m_xmax=xmax;
  m_Nb=nb;
  m_Bin1 = NULL;
  m_Bin2 = NULL;
  if( m_xmin >= m_xmax || m_Nb<1 ){
    std::cout<<"TFHST constructor failed, xmax<xmin or Nb<1 "<<std::endl;
    exit(0);
  }
  if(nb>0){
    m_Bin1 = new double[m_Nb+2];
    m_Bin2 = new double[m_Nb+2];
    if(m_Bin1 == NULL || m_Bin2 == NULL){
      std::cout<<"TFHST constructor failed to allocate Bin1 or Bin2"<<std::endl;
      exit(0);
    }
    for (i=0; i<=m_Nb+1; i++) m_Bin1[i]=0.0;
    for (i=0; i<=m_Nb+1; i++) m_Bin2[i]=0.0;
  }
}
///////////////////////////////////////////////////////////////////////////////
TFHST::TFHST(TFHST &From){
// Explicit copy constructor (unused, so far?)
  long i;
  m_xmin = From.m_xmin;
  m_xmax = From.m_xmax;
  m_Nb=From.m_Nb;
  m_Bin1 = NULL;
  m_Bin2 = NULL;
  if(m_Nb>0){
    m_Bin1 = new double[m_Nb+2];
    m_Bin2 = new double[m_Nb+2];
    if(m_Bin1 == NULL || m_Bin2 == NULL){
      std::cout<<"TFHST constructor failed to allocate Bin1 or Bin2"<<std::endl;
      exit(0);
    }
    for (i=0; i<=m_Nb+1; i++) m_Bin1[i]=From.m_Bin1[i];
    for (i=0; i<=m_Nb+1; i++) m_Bin2[i]=From.m_Bin2[i];
  }
  std::cout << "TFHST::TFHST NEVER UESE Copy constructor !!!" << std::endl; exit(0);
}
///////////////////////////////////////////////////////////////////////////////
TFHST::~TFHST(){
//    Explicit Destructor
  delete [] m_Bin1;
  delete [] m_Bin2;
}


///////////////////////////////////////////////////////////////////////////////
TFHST& TFHST::operator =(TFHST &From){
// substitution TFHST = TFHST, used
  int i;
  if (&From == this) return *this;
  // allocate (and delete) ONLY if necessary
  if (m_Nb != From.m_Nb) {
    delete [] m_Bin1;
    m_Bin1 = new double[m_Nb+2];
  }
  if( m_Nb != From.m_Nb) {
    delete [] m_Bin2;
    m_Bin2 = new double[m_Nb+2];
  }
  m_Nb=From.m_Nb;
  m_xmin = From.m_xmin;
  m_xmax = From.m_xmax;
  for (i=0; i<=m_Nb+1; i++) m_Bin1[i] = From.m_Bin1[i];
  for (i=0; i<=m_Nb+1; i++) m_Bin2[i] = From.m_Bin2[i];
  if(m_chat)  std::cout<<" SUBSITUTE operator = "<<std::endl;
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
double& TFHST::operator[](int nb){
// This is non like normal [], because hist[5]=x not alowed
  if( (nb>=0) &&  (nb<=(m_Nb+1)) ){
    return m_Bin1[nb];
  }else{
    std::cout<<" STOP in TFHST[], wrong bin index "<<std::endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////////////
double& TFHST::operator()(int nb){
// getter for erreor entry (not used in Foam)
  if( (nb>=0) &&  (nb<=(m_Nb+1)) ){
    return m_Bin2[nb];
  }else{
    std::cout<<" STOP in TFHST(), wrong bin index "<<std::endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////////////
void TFHST::Fill(double xx, double wt){
//   fill histogram with single weighted event
  int ib;
  ib = (int)( (xx-m_xmin)/(m_xmax-m_xmin) *m_Nb) +1;
  if(ib<1)    ib=0;       //unterflow
  if(ib>m_Nb) ib=m_Nb+1;  //overrflow
  m_Bin1[ib] += wt;
  m_Bin2[ib] += wt*wt;
  m_Entries  +=1.0;
}


///////////////////////////////////////////////////////////////////////////////
double TFHST::GetBinContent(int nb){
// Get content of bin nb.    nb=1 for 1-st bin, undeflow nb=0, overflow nb=Nb+1
  if( (nb>=0) &&  (nb<m_Nb+2) ){
    return m_Bin1[nb];
  }else{
    std::cout<<" STOP in TFHST[], wrong bin index "<<std::endl;
    exit(0);
  }
}

///////////////////////////////////////////////////////////////////////////////
double TFHST::GetBinError(int nb){
// Get content of bin nb.    nb=1 for 1-st bin, undeflow nb=0, overflow nb=Nb+1
  if( (nb>=0) &&  (nb<m_Nb+2) ){
    return sqrt(m_Bin2[nb]);
  }else{
    std::cout<<" STOP in TFHST[], wrong bin index "<<std::endl;
    exit(0);
  }
}
///////////////////////////////////////////////////////////////////////////////
void TFHST::Print(char *ch){
  // do nothing
}

///////////////////////////////////////////////////////////////////////////////
void TFHST::Print(void){
  int i,k;
  double bi,bi2,err,xx;
  double Lmax =80.0;
  std::cout << "-----------------------------------Histogram-------------------------------------------"<<std::endl;
  double bmin= +1e150;
  double bmax= -1e150;
  for(i=1; i<=m_Nb;i++){
    if(m_Bin1[i]>bmax) bmax=m_Bin1[i];
    if(m_Bin1[i]<bmin) bmin=m_Bin1[i];
  }
  if(bmax==bmin){
    std::cout<<" ++++ EMPTY Histogram "<<std::endl;
    return;
  }
  std::cout <<" min_bin ="<<SW12<<bmin<<"      ";
  std::cout <<" max_bin ="<<SW12<<bmax<<std::endl;
  for(i=0; i<=m_Nb+1; i++){
    bi  = m_Bin1[i];
    bi2 = m_Bin2[i];
    err =sqrt(bi2);
    int line =0;
    if(i>0 && i<m_Nb+1) {
      xx=m_xmin+(i-1)*((m_xmax-m_xmin)/m_Nb);
      line= (int)(Lmax*bi/(bmax-bmin));
      std::cout <<SW8<<xx <<"  ";
      std::cout <<SW12<<bi <<"  +- " <<SW12 << err << " ";
      std::cout<<"I";
      for(k=0;k<line;k++) std::cout<<"X";
    }else if(i==0){
      std::cout <<" Underflow";
      std::cout <<SW12<<bi <<"  +- " <<SW12 << err << " ";
    }else if(i==m_Nb+1){    
      std::cout <<" Overlow  ";
      std::cout <<SW12<<bi <<"  +- " <<SW12 << err << " ";
    }
    std::cout <<std::endl;
  }
  std::cout << "---------------------------------------------------------------------------------------"<<std::endl;
  std::cout <<std::endl;
}
///////////////////////////////////////////////////////////////////////////////
void TFHST::Reset(){
// Reset bin content
  int i;
  for (i=0; i<=m_Nb+1; i++) m_Bin1[i]=0.0;
  for (i=0; i<=m_Nb+1; i++) m_Bin2[i]=0.0;
}
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         End of Class of histograms                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

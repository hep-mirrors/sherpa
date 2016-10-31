#include "AMEGIC++/Amplitude/Zfunc.H"
#include "AMEGIC++/Amplitude/Pfunc.H"
#include "AMEGIC++/Amplitude/Zfunctions/Zfunc_Calc.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

#define new_calc

Zfunc::Zfunc() 
{
  p_equal  = this;
  m_narg   = 0;
  m_ncoupl = 0;
  m_nprop  = 0;
}

Zfunc::Zfunc(const Zfunc& z) 
{
  m_narg   = 0;
  m_ncoupl = 0;
  m_nprop  = 0;
  *this = z;
}
 
Zfunc::~Zfunc() {
  if (m_narg!=0)   delete[] p_arguments;
  if (m_ncoupl!=0) delete[] p_couplings;
  if (m_nprop!=0)  delete[] p_propagators;
  for (CL_Iterator clit=m_calclist.begin();clit!=m_calclist.end();++clit) {delete[] (*clit).p_args;}
}
    	
Zfunc& Zfunc::operator=(const Zfunc& z) {
  if (this!=&z) {
    if (m_narg!=0)   delete[] p_arguments;
    if (m_ncoupl!=0) delete[] p_couplings;
    if (m_nprop!=0)  delete[] p_propagators;
    
    m_type   = z.m_type;
    m_narg   = z.m_narg;
    m_ncoupl = z.m_ncoupl;
    m_nprop  = z.m_nprop;
    
    if (m_narg>0) {
      p_arguments = new int[m_narg];
      for (int i=0;i<m_narg;i++) p_arguments[i] = z.p_arguments[i];
    }
    if (m_ncoupl>0) {
      p_couplings = new Complex[m_ncoupl];
      for (int i=0;i<m_ncoupl;i++) p_couplings[i] = z.p_couplings[i];
    }
    if (m_nprop>0) {
      p_propagators = new Argument[m_nprop];
      for (int i=0;i<m_nprop;i++) p_propagators[i] = z.p_propagators[i];
    }
    
    m_sign   = z.m_sign;
    m_str    = z.m_str;
    
    for (CL_Iterator clit=m_calclist.begin();clit!=m_calclist.end();++clit) delete[] (*clit).p_args;
    m_calclist.clear();

    for (Calc_List::const_iterator clit=z.m_calclist.begin();clit!=z.m_calclist.end();++clit) {
      int* newargs = new int[2*m_narg];
      for (int j=0;j<2*m_narg;j++) newargs[j] = (*clit).p_args[j];
      m_calclist.push_back(CValue(newargs,(*clit).m_value));
    }
    
    p_equal      = z.p_equal;
    p_calculator = z.p_calculator;
  }
  return *this;
}

void Zfunc::ReplaceProp(vector<Pair>* pairlist)
{
  for(int i=0;i<m_narg;i++){
    for (size_t k=0;k<pairlist->size();k++) {
      if ((*pairlist)[k].pold==p_arguments[i]) {
	p_arguments[i] = (*pairlist)[k].pnew;
	break;
      }
    }
  }
  for (int i=0;i<m_nprop;i++) {
    for (size_t k=0;k<pairlist->size();k++) {
      if ((*pairlist)[k].pold==p_propagators[i].numb) {
	p_propagators[i].numb = (*pairlist)[k].pnew;
	break;
      }
    }
  }
}

void Zfunc::ClearCalcList()
{
  if (!m_calclist.empty()){
    for (CL_Iterator clit=m_calclist.begin();clit!=m_calclist.end();++clit) delete[] (*clit).p_args;
    m_calclist.clear();
  }
}

void Zfunc::Print() 
{
  if (!msg_LevelIsTracking()) return;

  msg_Out()<<"Z(["<<m_type<<"],";
  msg_Out()<<"[";
  for (int i=0;i<m_narg-1;i++) msg_Out()<<p_arguments[i]<<";";
  
  if (m_narg>0) msg_Out()<<p_arguments[m_narg-1];
  msg_Out()<<"][";
  msg->Out().precision(2);
  for (int i=0;i<m_ncoupl-1;i++) {
    if ( real(p_couplings[i])!=0.0 &&
	 imag(p_couplings[i])==0.0 )
      msg_Out()<<real(p_couplings[i])<<";";
    if (  real(p_couplings[i])==0.0 &&
	  imag(p_couplings[i])!=0.0 )
      msg_Out()<<imag(p_couplings[i])<<" I;";
    if ( real(p_couplings[i])!=0.0 &&
	 imag(p_couplings[i])!=0.0 )
      msg_Out()<<real(p_couplings[i])<<"+"<<imag(p_couplings[i])<<" I;";
    if (  real(p_couplings[i])==0.0 &&
	  imag(p_couplings[i])==0.0 )
      msg_Out()<<"0;";
  }
  if ( real(p_couplings[m_ncoupl-1])!=0.0 &&
       imag(p_couplings[m_ncoupl-1])==0.0 )
      msg_Out()<<real(p_couplings[m_ncoupl-1])<<"])";
  if (  real(p_couplings[m_ncoupl-1])==0.0 &&
	imag(p_couplings[m_ncoupl-1])!=0.0 )
    msg_Out()<<imag(p_couplings[m_ncoupl-1])<<" I])";
  if ( real(p_couplings[m_ncoupl-1])!=0.0 &&
       imag(p_couplings[m_ncoupl-1])!=0.0 )
    msg_Out()<<real(p_couplings[m_ncoupl-1])<<"+"<<imag(p_couplings[m_ncoupl-1])<<" I])";
  if (  real(p_couplings[m_ncoupl-1])==0.0 &&
	imag(p_couplings[m_ncoupl-1])==0.0 )
	msg_Out()<<"0])";
  msg_Out()<<endl;
  msg->Out().precision(6);
}




Zfunc_Group::Zfunc_Group(const Zfunc& z) : Zfunc(z) 
{
  m_op = '+';
  m_sign = 1;
  if (m_nprop!=0) delete[] p_propagators;
  m_nprop = 0;
  p_equal = this;
}

Zfunc_Group::Zfunc_Group(Zfunc& z1,Zfunc& z2,int si,Pfunc_List* pl)
{
  int n99 = 0;                    //Marker 99
  int i1  = 0;
  for(int i=0;i<z1.m_narg;i++){
    if (z1.p_arguments[i]==si) i1++;
    if (z1.p_arguments[i]==99) n99++;   //polarisation dummys not needed in Zfunc_Group
  }
  int i2  = 0;
  for(int i=0;i<z2.m_narg;i++){
    if (z2.p_arguments[i]==si) i2++;
    if (z2.p_arguments[i]==99) n99++;
  }
  if ( (i1==0&&i2>0) || (i2==0&&i1>0) ) {
    msg_Error()<<"Error in Zfunc_Group(Z*Z-Constructor): sum index, will abort."<<endl;
    Abort();
  }

  m_type   = "";
  m_narg   = z1.m_narg   + z2.m_narg - i1 - i2 - n99;
  m_ncoupl = z1.m_ncoupl + z2.m_ncoupl;
  m_nprop  = 0;

  int ii = 0;
  if (m_narg>0) {
    p_arguments = new int[m_narg];
    while (ii<m_narg) {
      for (int i=0;i<z1.m_narg;i++) 
	if ( z1.p_arguments[i]!=si && z1.p_arguments[i]!=99 ){
	  p_arguments[ii] = z1.p_arguments[i];
	  ii++;
	}
      for (int i=0;i<z2.m_narg;i++) 
	if( z2.p_arguments[i]!=si && z2.p_arguments[i]!=99 ){
	  p_arguments[ii] = z2.p_arguments[i];
	  ii++;
	}
    }
  }

  if (z1.GetOp()=='*') m_nprop += z1.m_nprop;
  if (z1.GetOp()==0) {
    for (int i=0;i<z1.m_nprop;i++){
      for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
	Pfunc*  p = *pit;
	if (p->arg[0]==z1.p_propagators[i].numb && p->on==0){
	  m_nprop++;
	  break;
	}
      } 
    }
  }
  int k = m_nprop;

  if (z2.GetOp()=='*') m_nprop += z2.m_nprop;
  if (z2.GetOp()==0) {
    for (int i=0;i<z2.m_nprop;i++){
      for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
	Pfunc*  p = *pit;
	if (p->arg[0]==z2.p_propagators[i].numb && p->on==0){
	  m_nprop++;
	  break;
	}
      } 
    }
  }

  m_nprop++;

  if (m_nprop>0) {
    p_propagators = new Argument[m_nprop];
    if ( (z1.GetOp()=='*') && (z1.m_nprop>0) ) { 
      for (int i=0;i<k;i++) p_propagators[i] = z1.p_propagators[i];
      k = z1.m_nprop;
      delete[] z1.p_propagators;
      z1.m_nprop = 0;
    }
    if (z1.GetOp()==0) {
      int j = 0;
      for (int i=0;i<z1.m_nprop;i++) {
	for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
	  Pfunc*  p = *pit;
	  if (p->arg[0]==z1.p_propagators[i].numb && p->on==0){
	    p_propagators[j] = z1.p_propagators[i];
	    j++;
	    break;
	  }
	} 
      }
    }
    if ( (z2.GetOp()=='*') && (z2.m_nprop>0) ) {
      for (int i=0;i<z2.m_nprop;i++) p_propagators[i+k] = z2.p_propagators[i];
      delete[] z2.p_propagators;
      z2.m_nprop = 0;      
    }
    if (z2.GetOp()==0) {
      int j = k;
      for (int i=0;i<z2.m_nprop;i++){
	for (Pfunc_Iterator pit=pl->begin();pit!=pl->end();++pit) {
	  Pfunc*  p = *pit;
	  if (p->arg[0]==z2.p_propagators[i].numb && p->on==0){
	    p_propagators[j] = z2.p_propagators[i];
	    j++;
	    break;
	  }
	} 
      }
    }
    p_propagators[m_nprop-1].numb = si;
  }
	
  if (m_ncoupl>0) {
    p_couplings = new Complex[m_ncoupl];
    for (int i=0;i<z1.m_ncoupl;i++) p_couplings[i] = z1.p_couplings[i];
    for (int i=z1.m_ncoupl;i<z1.m_ncoupl+z2.m_ncoupl;i++) p_couplings[i] = z2.p_couplings[i-z1.m_ncoupl];
  }

  p_equal = this;
  m_sign     = 1;
  m_sumindex = si;
  m_op       = '*';
  m_zlist.push_back(&z1);
  m_zsigns.push_back(1);
  m_zlist.push_back(&z2);
  m_zsigns.push_back(1);
}

void Zfunc_Group::ReplaceProp(vector<Pair>* pairlist)
{
  for (size_t k=0;k<pairlist->size();k++) {
    if ((*pairlist)[k].pold==m_sumindex) {
      m_sumindex = (*pairlist)[k].pnew;
      break;
    }
  }
  Zfunc::ReplaceProp(pairlist);
  for (size_t i=0;i<m_zlist.size();i++) m_zlist[i]->ReplaceProp(pairlist);
}

void Zfunc_Group::ClearCalcList()
{
  Zfunc::ClearCalcList();
  for (size_t i=0;i<m_zlist.size();i++) m_zlist[i]->ClearCalcList();
}

void Zfunc_Group::KillZList()
{
  for (Zfunc_Iterator zit=m_zlist.begin();zit!=m_zlist.end();++zit) {
    (*zit)->KillZList();
    delete (*zit);
  }  
}

void Zfunc_Group::Print() 
{
  if (!msg_LevelIsTracking()) return;
  msg_Out()<<"SZ(["<<m_type<<"],";
  msg_Out()<<"[";
  for (int i=0;i<m_narg-1;i++) msg_Out()<<p_arguments[i]<<";";
  
  if (m_narg>0) msg_Out()<<p_arguments[m_narg-1];
  msg_Out()<<"])";
  msg_Out()<<endl;

  if (m_op=='+'){  
    for (size_t i=0;i<m_zlist.size();i++) {
      if (m_zsigns[i]==-1) {
	msg_Out()<<"   - "; //<<m_zlist[i]->p_propagators[0].numb<<" * ";
	m_zlist[i]->Print();
      }
      else {
	if (m_zlist[i]->p_propagators!=NULL) {
	  msg_Out()<<"   + ";
	  m_zlist[i]->Print(); 
	}
	else {
	  msg_Out()<<" ??? "<<" * ";m_zlist[i]->Print();
	}
      }
    }
  }
  if (m_op=='*'){
    for (size_t i=0;i<m_zlist.size();i++) {
      if (i>0) msg_Out()<<"  *";else msg_Out()<<" ->";
      m_zlist[i]->Print();
    }
    msg_Out()<<"Sum over "<<m_sumindex<<endl;
  }  
}

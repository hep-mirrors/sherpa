#include "Lorentz_Function.H"

using namespace AMEGIC;

void Lorentz_Function::AddPermutation(int sign,int a,int b=-1,int c=-1,int d=-1)
{
  int* newperm = new int[NofIndex()];
  newperm[0] = partarg[a];
  if (NofIndex()>1) newperm[1] = partarg[b];
  if (NofIndex()>2) newperm[2] = partarg[c];
  if (NofIndex()>3) newperm[3] = partarg[d];

  m_permlist.push_back(newperm);
  m_signlist.push_back(sign);
}

void Lorentz_Function::InitPermutation()
{
  if (!m_permlist.empty()) {
    for (short int i=0;i<m_permlist.size();i++) delete[] m_permlist[i]; 
    m_permlist.clear();
    m_signlist.clear();
  }

  switch (type) {
  case lf::Gab   : 
    AddPermutation(1,0,1);
    AddPermutation(1,1,0);  
    break;
  case lf::VVSS   : 
    AddPermutation(1,0,1);
    AddPermutation(1,1,0);  
    break;
  case lf::Gauge3:
    AddPermutation( 1,0,1,2);
    AddPermutation(-1,0,2,1);  
    AddPermutation(-1,1,0,2);
    AddPermutation(-1,2,1,0);  
    AddPermutation( 1,1,2,0);
    AddPermutation( 1,2,0,1);  
    break;
  case lf::Gauge4: 
    AddPermutation( 1,0,1,2,3);
    AddPermutation( 1,1,0,2,3);
    AddPermutation( 1,0,1,3,2);
    AddPermutation( 1,1,0,3,2);
    AddPermutation( 1,2,3,1,0);    
    AddPermutation( 1,3,2,1,0);
    AddPermutation( 1,2,3,0,1);
    AddPermutation( 1,3,2,0,1);
    break;
  case lf::Gluon4: 
    AddPermutation( 1,0,1,2,3);
    AddPermutation(-1,2,1,0,3);
    AddPermutation(-1,0,3,2,1);
    AddPermutation( 1,2,3,0,1);
    AddPermutation( 1,1,0,3,2);
    AddPermutation(-1,3,0,1,2);
    AddPermutation(-1,1,2,3,0);
    AddPermutation( 1,3,2,1,0);
    break;
  case lf::SSV   : 
    AddPermutation( 1,0,1,2);
    AddPermutation(-1,1,0,2);
    break;
  case lf::VVT:
    AddPermutation( 1,0,1,2);
    AddPermutation( 1,1,0,2);
    break;
  case lf::SST:
    AddPermutation( 1,0,1,2);
    AddPermutation( 1,1,0,2);
    break;    
  case lf::VVGS:
    AddPermutation( 1,0,1,2);
    AddPermutation( 1,1,0,2);
    break;
  case lf::SSGS:
    AddPermutation( 1,0,1);
    AddPermutation( 1,1,0);
    break;
  case lf::VVVT:
    AddPermutation( 1,0,1,2,3);
    AddPermutation(-1,0,2,1,3);  
    AddPermutation(-1,1,0,2,3);
    AddPermutation(-1,2,1,0,3);  
    AddPermutation( 1,1,2,0,3);
    AddPermutation( 1,2,0,1,3);  
    break;    
  }
  m_permcount = 0;
}

int Lorentz_Function::ResetPermutation() 
{
  m_permcount=0;
  for (short int i=0;i<NofIndex();i++) partarg[i]  = m_permlist[m_permcount][i];
}

int Lorentz_Function::NextPermutation()
{
  if (NofIndex()<2) return 0;
  m_permcount++;
  if (m_permcount==m_permlist.size()) return 0;
  
  for (short int i=0;i<NofIndex();i++) partarg[i]  = m_permlist[m_permcount][i];
					 
}

int Lorentz_Function::GetSign() 
{
  if (m_signlist.empty()) return 1;
  return m_signlist[m_permcount];
}

void AMEGIC::Lorentz_Function2MPI(const Lorentz_Function * lf , MPI_Lorentz_Function & mpi_lf) {
  
  mpi_lf.m_type =  lf->type;
  for (int i=0; i<4; ++i)  
    mpi_lf.m_partarg[i] = lf->partarg[i];
}

Lorentz_Function * AMEGIC::MPI2Lorentz_Function(const MPI_Lorentz_Function & mpi_lf ) {

  Lorentz_Function * lf ;
  
  lf = new Lorentz_Function((AMEGIC::lf::code)(mpi_lf.m_type));
  for (int i=0; i<4; ++i)
    lf->partarg[i] = mpi_lf.m_partarg[i];
    
  return lf;
}

std::ostream & AMEGIC::operator<<(std::ostream & s, const MPI_Lorentz_Function & lf) {
  s<<lf.m_type<<",";
  s<<lf.m_partarg[0]<<","<<lf.m_partarg[1]<<","<<lf.m_partarg[2]<<","<<lf.m_partarg[3];
  return s;
}

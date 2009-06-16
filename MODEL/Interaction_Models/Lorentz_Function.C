#include "MODEL/Interaction_Models/Lorentz_Function.H"
#include "ATOOLS/Org/Message.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE MODEL::LF_Key
#define OBJECT_TYPE MODEL::Lorentz_Function
#include "ATOOLS/Org/Getter_Function.C"

using namespace MODEL;
#include "ATOOLS/Org/Exception.H"
Lorentz_Function::Lorentz_Function(const std::string &type): 
  m_type(type), p_next(NULL) 
{
  for (short int i=0;i<4;i++) m_partarg[i] = -1;
}

Lorentz_Function::~Lorentz_Function() 
{
  for (size_t i=0; i<m_permlist.size();++i) delete [] m_permlist[i];
  if (p_next) delete p_next;
}

void Lorentz_Function::SetParticleArg(int a,int b,int c,int d) 
{
  m_partarg[0] = a;
  m_partarg[1] = b;
  m_partarg[2] = c;
  m_partarg[3] = d;
}

void Lorentz_Function::AddPermutation(int sign,int a,int b,int c,int d)
{
  int* newperm = new int[NofIndex()];
  newperm[0] = m_partarg[a];
  if (NofIndex()>1) newperm[1] = m_partarg[b];
  if (NofIndex()>2) newperm[2] = m_partarg[c];
  if (NofIndex()>3) newperm[3] = m_partarg[d];

  m_permlist.push_back(newperm);
  m_signlist.push_back(sign);
}

void Lorentz_Function::InitPermutation()
{
  if (!m_permlist.empty()) {
    for (size_t i=0;i<m_permlist.size();i++) delete[] m_permlist[i]; 
    m_permlist.clear();
    m_signlist.clear();
  }
  m_permcount = 0;
}

int Lorentz_Function::ResetPermutation() 
{
  m_permcount=0;
  for (short int i=0;i<NofIndex();i++) m_partarg[i]  = m_permlist[m_permcount][i];
  return 1;
}

int Lorentz_Function::NextPermutation()
{
  if (NofIndex()<2) return 0;
  m_permcount++;
  if (m_permcount==(int)m_permlist.size()) return 0;
  
  for (short int i=0;i<NofIndex();i++) m_partarg[i]  = m_permlist[m_permcount][i];
  return 1;
}

int Lorentz_Function::GetSign() 
{
  if (m_signlist.empty()) return 1;
  return m_signlist[m_permcount];
}

void Lorentz_Function2MPI(const Lorentz_Function * lf , MPI_Lorentz_Function & mpi_lf) 
{  
  mpi_lf.m_type =  lf->Type();
  for (int i=0; i<4; ++i)  
    mpi_lf.m_partarg[i] = lf->ParticleArg(i);
}

Lorentz_Function * MPI2Lorentz_Function(const MPI_Lorentz_Function & mpi_lf ) 
{
  Lorentz_Function * lf = LF_Getter::GetObject(mpi_lf.m_type,LF_Key());
  lf->SetParticleArg(mpi_lf.m_partarg[0],mpi_lf.m_partarg[1],mpi_lf.m_partarg[2],mpi_lf.m_partarg[3]);
  return lf;
}

std::string Lorentz_Function::Str(int a) const
{
  MyStrStream sstr;
  sstr<<m_partarg[a];
  std::string help;
  sstr>>help;
  return help;
} 


std::ostream & MODEL::operator<<(std::ostream & s, const MPI_Lorentz_Function & lf) {
  s<<lf.m_type<<",";
  s<<lf.m_partarg[0]<<","<<lf.m_partarg[1]<<","<<lf.m_partarg[2]<<","<<lf.m_partarg[3];
  return s;
}


Lorentz_Function & Lorentz_Function::operator=(const Lorentz_Function & l)
{
  if (this!=&l) {
    if (m_type!=l.m_type) THROW(fatal_error,"Internal error");
    m_permcount=l.m_permcount;
    int noi=l.NofIndex();

    for (size_t i=0; i<m_permlist.size();++i) delete [] m_permlist[i];
    m_permlist.clear();
    m_signlist.clear();
    if (p_next) delete p_next;

    for (size_t i=0; i<l.m_permlist.size();++i) {
      m_signlist[i]=l.m_signlist[i];
      m_permlist.push_back(new int[noi]);
      for (int j=0; j<noi;++j) {
	m_permlist[i][j]=l.m_permlist[i][j];
      }
    }
    for (int i=0; i<4;++i) 
      m_partarg[i]=l.m_partarg[i];
    if (l.p_next!=0) p_next = l.p_next->GetCopy();
    else p_next = 0;
  }
  return *this;
}

bool Lorentz_Function::CutVectors() 
{ 
  return false; 
}

#include "Lorentz_Function.H"
#include "Message.H"

using namespace AMEGIC;

void Lorentz_Function::AddPermutation(int sign,int a,int b=-1,int c=-1,int d=-1)
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

  switch (m_type) {
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
  case lf::Triangle   : 
    AddPermutation(1,0,1);
    AddPermutation(1,1,0);  
    break;
  case lf::Box:
    AddPermutation( 1,0,1,2);
    AddPermutation(-1,0,2,1);  
    AddPermutation(-1,1,0,2);
    AddPermutation(-1,2,1,0);  
    AddPermutation( 1,1,2,0);
    AddPermutation( 1,2,0,1);  
    break;
  case lf::C4GS   : 
     AddPermutation(1,0,1);
     AddPermutation(1,1,0);  
     break;
  default : break;
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

void AMEGIC::Lorentz_Function2MPI(const Lorentz_Function * lf , MPI_Lorentz_Function & mpi_lf) 
{  
  mpi_lf.m_type =  lf->Type();
  for (int i=0; i<4; ++i)  
    mpi_lf.m_partarg[i] = lf->ParticleArg(i);
}

Lorentz_Function * AMEGIC::MPI2Lorentz_Function(const MPI_Lorentz_Function & mpi_lf ) 
{
  Lorentz_Function * lf = new Lorentz_Function((AMEGIC::lf::code)(mpi_lf.m_type));
  lf->SetParticleArg(mpi_lf.m_partarg[0],mpi_lf.m_partarg[1],mpi_lf.m_partarg[2],mpi_lf.m_partarg[3]);
  return lf;
}

std::string Lorentz_Function::String(int shortversion) const 
{
  if (m_type==lf::SSS)  return std::string("1");
  if (m_type==lf::FFS)  return std::string("1");
  if (m_type==lf::SSSS) return std::string("1");
  std::string help;
  switch (m_type) {
  case lf::Gamma:  
    // Gam[0]
    help = std::string("Gam[") + Str(0) + std::string("]");break;
  case lf::Pol:  
    // Eps[0]
    help = std::string("Eps[") + Str(0) + std::string("]");break;
  case lf::Gab:   
    // G[0,1]
    help = std::string("G[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
  case lf::VVSS:   
    // G(2V2S)[0,1]
    help = std::string("G(2V2S)[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
  case lf::Gauge3: 
    // (P[0,2]-P[1,2])*G(0,1)+(P[1,0]-P[2,0])*G(1,2)+(P[2,1]-P[0,1])*G(2,0)
    if (shortversion) {
      help += std::string("V3[") + Str(0) + std::string(",") + 
	Str(1) + std::string(",") + 
	Str(2) + std::string("]");
    }
    else {
      help  = std::string("(P[") + Str(0) + std::string(",") + Str(2) + std::string("]-");
      help += std::string("P[")  + Str(1) + std::string(",") + Str(2) + std::string("])*");
      help += std::string("G[")  + Str(0) + std::string(",") + Str(1) + std::string("]");
	  
      help += std::string("+");
	  
      help += std::string("(P[") + Str(1) + std::string(",") + Str(0) + std::string("]-");
      help += std::string("P[")  + Str(2) + std::string(",") + Str(0) + std::string("])*");
      help += std::string("G[")  + Str(1) + std::string(",") + Str(2) + std::string("]");
	  
      help += std::string("+");
	  
      help += std::string("(P[") + Str(2) + std::string(",") + Str(1) + std::string("]-");
      help += std::string("P[")  + Str(0) + std::string(",") + Str(1) + std::string("])*");
      help += std::string("G[")  + Str(2) + std::string(",") + Str(0) + std::string("]");
    }
    break;
  case lf::SSV:
    //P[0,2]-P[1,2]
    help = std::string("P[") + Str(0) + std::string(",") + Str(2) +std::string("]-"); 
    help += std::string("P[") + Str(1) + std::string(",") + Str(2) +std::string("]");
    break;
  case lf::Gauge4: 
    //(2G(0,1)*G(2,3)-G(0,2)*G(1,3)-G(0,3)*G(1,2))
    help  = std::string("(2*G[")  + Str(0) + std::string(",") + Str(1) + std::string("]*");
    help += std::string("G[")  + Str(2) + std::string(",") + Str(3) + std::string("]-");
    help += std::string("G[")  + Str(0) + std::string(",") + Str(2) + std::string("]*");
    help += std::string("G[")  + Str(1) + std::string(",") + Str(3) + std::string("]-");
    help += std::string("G[")  + Str(0) + std::string(",") + Str(3) + std::string("]*");
    help += std::string("G[")  + Str(1) + std::string(",") + Str(2) + std::string("])");
    break;
  case lf::Gluon4:
    //G(0,1)*G(2,3)-G(0,3)*G(2,1)
    if (shortversion) {
      help += std::string("G4[") + Str(0) + std::string(",") + 
	Str(1) + std::string(",") + 
	Str(2) + std::string(",") + 
	Str(3) + std::string("]");
    }
    else {
      help  = std::string("(G[")  + Str(0) + std::string(",") + Str(1) + std::string("]*");
      help += std::string("G[")  + Str(2) + std::string(",") + Str(3) + std::string("]-");
      help += std::string("G[")  + Str(0) + std::string(",") + Str(3) + std::string("]*");
      help += std::string("G[")  + Str(2) + std::string(",") + Str(1) + std::string("])");
    }
    break;
  case lf::FFT:
    help = std::string("FFT[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
  case lf::FFVT:
    help = std::string("FFVT[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
  case lf::FFVGS:
    help = std::string("FFVGS[") + Str(0) + std::string("]");break;
  case lf::VVT:
    help = std::string("VVT[") + Str(0) + std::string(",") + Str(1) 
      + std::string(",") + Str(2) + std::string("]");break;    
  case lf::VVGS:
    help = std::string("VVGS[") + Str(0) + std::string(",") + Str(1) 
      + std::string(",") + Str(2) + std::string("]");break;
  case lf::Triangle:   
    // G[0,1]
    help = std::string("T[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
  case lf::Box:   
    // G[0,1]
    help = std::string("B[") + Str(0) + std::string(",") + Str(1) + std::string(",") + 
      Str(2) + std::string("]");break;
  case lf::C4GS:   
    // G[0,1]
    help = std::string("AddOn5Vertex[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
  default :
    return std::string("1");
  }     
  return help;
}

std::string Lorentz_Function::Str(int a) const
{
  MyStrStream sstr;
  sstr<<m_partarg[a];
  std::string help;
  sstr>>help;
  return help;
} 


std::ostream & AMEGIC::operator<<(std::ostream & s, const MPI_Lorentz_Function & lf) {
  s<<lf.m_type<<",";
  s<<lf.m_partarg[0]<<","<<lf.m_partarg[1]<<","<<lf.m_partarg[2]<<","<<lf.m_partarg[3];
  return s;
}


Lorentz_Function & Lorentz_Function::operator=(const Lorentz_Function & l)
{
  if (this!=&l) {
    m_type     =l.m_type;
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
    if (l.p_next!=0) 
      p_next = new Lorentz_Function(*(l.p_next));
    else
      p_next = 0;
  }
  return *this;
}

Lorentz_Function::~Lorentz_Function() {
  for (size_t i=0; i<m_permlist.size();++i) delete [] m_permlist[i];
  if (p_next) delete p_next;
}

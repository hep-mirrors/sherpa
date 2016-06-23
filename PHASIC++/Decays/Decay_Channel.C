#include "PHASIC++/Decays/Decay_Channel.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Color.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/My_MPI.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"
#include "METOOLS/SpinCorrelations/Spin_Density.H"
#include <algorithm>

using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

Decay_Channel::Decay_Channel(const Flavour & _flin,
                             const ATOOLS::Mass_Selector* ms) :
  m_width(0.), m_deltawidth(-1.), m_minmass(0.), m_max(0.), m_symfac(-1.0),
  m_iwidth(0.), m_ideltawidth(-1.), m_active(1),
  p_channels(NULL), p_amps(NULL), p_ms(ms)
{
  m_flavours.push_back(_flin);
}

Decay_Channel::~Decay_Channel()
{
  for (size_t i(0); i<m_diagrams.size(); ++i) {
    delete m_diagrams[i];
  }
  if (p_channels) delete p_channels;
  if (p_amps) delete p_amps;
}

bool Decay_Channel::FlavourSort(const Flavour &fl1,const Flavour &fl2)
{
  // TODO: Get rid of this custom sorting, but then the hadron decay channel
  // files have to be changed as well (order mapping in MEs)
  kf_code kf1(fl1.Kfcode()), kf2(fl2.Kfcode());
  if (kf1>kf2) return true;
  if (kf1<kf2) return false;
  /*
      anti anti -> true
      anti part -> false
      part anti -> true
      anti anti -> true
      */
  return !(fl1.IsAnti()&&!fl2.IsAnti());
}


void Decay_Channel::AddDecayProduct(const ATOOLS::Flavour& flout,
				    const bool & sort)
{
  m_flavours.push_back(flout);
  if (!sort) return;
  // sort
  Flavour flin=m_flavours[0];
  Flavour_Vector flouts(m_flavours.size()-1);
  for (size_t i=1; i<m_flavours.size(); ++i) {
    flouts[i-1]=m_flavours[i];
  }
  std::sort(flouts.begin(), flouts.end(),Decay_Channel::FlavourSort);
  m_flavours.clear();
  m_flavours.resize(flouts.size()+1);
  m_flavours[0]=flin;
  for (size_t i=0; i<flouts.size(); ++i) {
    m_flavours[i+1]=flouts[i];
  }

  m_minmass += p_ms->Mass(flout);
}

void Decay_Channel::AddDiagram(METOOLS::Spin_Amplitudes* amp) {
  m_diagrams.push_back(amp);
}

void Decay_Channel::AddChannel(PHASIC::Single_Channel* chan)
{
  p_channels->Add(chan);
}

void Decay_Channel::ResetChannels()
{
  p_channels->DropAllChannels(false);
  p_channels->Reset();
}

void Decay_Channel::Output() const
{
  msg_Out()<<(*this);
}

namespace PHASIC {
  std::ostream &operator<<(std::ostream &os,const Decay_Channel &dc)
  {
    os<<left<<setw(18)<<dc.IDCode();
    os<<setw(25)<<dc.Name();
    os<<setw(10)<<dc.m_width;
    if (dc.m_deltawidth>0.) os<<"("<<setw(10)<<dc.m_deltawidth<<")";
    os<<" GeV";
    if (dc.Active()<1) {
      os<<" [disabled]";
    }
    return os;
  }
}

string Decay_Channel::Name() const
{
  string name=m_flavours[0].IDName()+string(" --> ");
  for (size_t i=1; i<m_flavours.size(); ++i) {
    name+=m_flavours[i].IDName()+string(" ");
  }
  return name;
}

string Decay_Channel::IDCode() const
{
  string code=ToString((long int)m_flavours[0]);
  for (size_t i=1; i<m_flavours.size(); ++i) {
    code+=","+ToString((long int)m_flavours[i]);
  }
  return code;
}

double Decay_Channel::Lambda(const double& a,
                             const double& b,
                             const double& c) const
{
  double L = (sqr(a-b-c)-4.*b*c);
  if (L>0.0) return sqrt(L)/2/sqrt(a);
  if (L>-Accu()) return 0.0;
  msg_Error()<<"passed impossible mass combination:"<<std::endl;
  msg_Error()<<"m_a="<<sqrt(a)<<" m_b="<<sqrt(b)<<" m_c="<<sqrt(c)<<endl;
  msg_Error()<<"L="<<L<<endl;
  return 0.;
}

double Decay_Channel::MassWeight(const double& s,
                                 const double& sp,
                                 const double& b,
                                 const double& c) const
{
  return Lambda(sp,b,c)/Lambda(s,b,c)*s/sp;
}

double Decay_Channel::GenerateMass(const double& max,
                                   const double& width) const
{
  double mass=-1.0;
  double decaymin = MinimalMass();
  DEBUG_FUNC(decaymin<<" < m["<<GetDecaying()<<"] < "<<max);
  if(decaymin>max) mass=-1.0;
  else if (decaymin==0.0) {
    mass=m_flavours[0].RelBWMass(decaymin, max,
                                 p_ms->Mass(m_flavours[0]), width);
  }
  else {
    double s=sqr(p_ms->Mass(GetDecaying()));
    double mb(0.0), mc(0.0);
    for (int i=0; i<NOut(); ++i) {
      mc+=p_ms->Mass(GetDecayProduct(i));
      if(p_ms->Mass(GetDecayProduct(i))>mb)
        mb=p_ms->Mass(GetDecayProduct(i));
    }
    mc-=mb;
    double b=sqr(mb);
    double c=sqr(mc);
    double spmax=2.0*b+2.0*c+sqrt(sqr(b)+14.0*b*c+sqr(c));
    double wmax=MassWeight(s,spmax,b,c);
    double w=0.0;
    int trials(0);
    do {
      mass=m_flavours[0].RelBWMass(decaymin, max,
                                   p_ms->Mass(m_flavours[0]), width);
      double sp=sqr(mass);
      w=MassWeight(s,sp,b,c);
      ++trials;
      if (w>wmax+Accu())
        msg_Error()<<METHOD<<" w="<<w<<" > wmax="<<wmax<<std::endl;
    } while (w<ran->Get()*wmax && trials<1000);
  }
  DEBUG_VAR(mass);
  return mass;
}

double Decay_Channel::SymmetryFactor()
{
  if (m_symfac<0.0) {
    std::map<Flavour,size_t> fc;
    for (size_t i=1; i<m_flavours.size(); ++i) {
      std::map<Flavour,size_t>::iterator fit(fc.find(m_flavours[i]));
      if (fit==fc.end()) fit=fc.insert(make_pair(m_flavours[i],0)).first;
      ++fit->second;
    }
    m_symfac=1.0;
    for (std::map<Flavour,size_t>::const_iterator fit(fc.begin());
         fit!=fc.end();++fit) {
      m_symfac*=Factorial(fit->second);
    }
  }
  return m_symfac;
}

void Decay_Channel::CalculateWidth(double acc, double ref, int iter)
{
  p_channels->Reset();
  int maxopt    = p_channels->Number()*int(pow(2.,2*(int(NOut())-2)));

  int      opt=0;
  double   value, sum=0., sum2=0., result=1., disc;
  double   n=0., mv[3]={0.,0.,0.};

  double flux(1./(2.*p_ms->Mass(GetDecaying())));
  std::vector<Vec4D> momenta(1+NOut());
  momenta[0] = Vec4D(p_ms->Mass(GetDecaying()),0.,0.,0.);
  ref/=flux;
  double crit = (ref>0.0?ref:result);
  m_ideltawidth=crit;

  while(opt<maxopt && m_ideltawidth>acc*crit) {
    for (int ln=1;ln<iter+1;ln++) {
      value = Differential(momenta, false, NULL);
      mv[0] += 1.0;
      mv[1] += value;
      mv[2] += ATOOLS::sqr(value);
      p_channels->AddPoint(value);
      if (value>m_max) {
        m_max = value;
      }
    }
    opt++;
#ifdef USING__MPI
    if (MPI::COMM_WORLD.Get_size()) {
      mpi->MPIComm()->Allreduce(MPI_IN_PLACE,mv,3,MPI::DOUBLE,MPI::SUM);
      mpi->MPIComm()->Allreduce(MPI_IN_PLACE,&m_max,1,MPI::DOUBLE,MPI::MAX);
    }
#endif
    n+=mv[0];
    sum+=mv[1];
    sum2+=mv[2];
    mv[0]=mv[1]=mv[2]=0.0;

    p_channels->MPISync();
    p_channels->Optimize(0.01);

    result = sum/n;
    disc   = sqr(sum/n)/((sum2/n - sqr(sum/n))/(n-1));
    if (disc!=0.0) m_ideltawidth  = result/sqrt(abs(disc));
    crit = (ref>0.0?ref:result);
  }

  m_iwidth  = flux*sum/n;
  m_ideltawidth *= flux;
  disc   = sqr(m_iwidth)/((sum2*sqr(flux)/n - sqr(m_iwidth))/(n-1));
  if (disc!=0.0) m_ideltawidth  = m_iwidth/sqrt(abs(disc));
  if(abs(m_ideltawidth)/m_iwidth<1e-6) m_ideltawidth=0.0;
}

double Decay_Channel::Differential(ATOOLS::Vec4D_Vector& momenta, bool anti,
                                   METOOLS::Spin_Density* sigma,
                                   const std::vector<ATOOLS::Particle*>& p)
{
  Poincare labboost(momenta[0]);
  labboost.Boost(momenta[0]);
  Channels()->GeneratePoint(&momenta.front(),NULL);
  Channels()->GenerateWeight(&momenta.front(),NULL);

  labboost.Invert();
  for (size_t i(0); i<momenta.size(); ++i) labboost.Boost(momenta[i]);
  double dsigma_lab=ME2(momenta, anti, sigma, p);
  return dsigma_lab*Channels()->Weight();
}

double Decay_Channel::ME2(const ATOOLS::Vec4D_Vector& momenta, bool anti,
                          METOOLS::Spin_Density* sigma,
                          const std::vector<ATOOLS::Particle*>& p)
{
  if (GetDiagrams().size()<1) return 0.0;

  for(size_t i(0); i<GetDiagrams().size(); ++i) {
    GetDiagrams()[i]->Calculate(momenta, anti);
  }

  Complex sumijlambda_AiAj(0.0,0.0);

  if (sigma) {
    for (size_t i(0); i<m_diagrams.size(); ++i) DEBUG_VAR(*m_diagrams[i]);
    if (p_amps) delete p_amps;
    vector<int> spin_i(p.size(), -1), spin_j(p.size(), -1);
    p_amps=new Amplitude2_Tensor(p,0,m_diagrams,spin_i, spin_j);
    DEBUG_VAR(*p_amps);
    sumijlambda_AiAj=(*sigma)*p_amps->ReduceToMatrix(sigma->Particle());
  }
  else {
    for (size_t i(0); i<GetDiagrams().size(); ++i) {
      Spin_Amplitudes* Ai=GetDiagrams()[i];
      for (size_t j(0); j<GetDiagrams().size(); ++j) { // 0?
        Spin_Amplitudes* Aj=GetDiagrams()[j];
  
        // for debugging:
        if (Ai->size()!=Aj->size())
          THROW(fatal_error,"Trying to multiply two amplitudes with different "+
                string("number of helicity combinations."));
  
        for (size_t lambda=0; lambda<Ai->size(); ++lambda) {
          sumijlambda_AiAj+=(*Ai)[lambda]*conj((*Aj)[lambda]);
        }
      }
    }
  }
  if (!IsZero(sumijlambda_AiAj.imag(),1.0e-6)) {
    PRINT_INFO("Sum-Squaring matrix element yielded imaginary part.");
    PRINT_VAR(sumijlambda_AiAj);
  }

  double value=sumijlambda_AiAj.real();
  value /= double(GetDecaying().IntSpin()+1);
  if (GetDecaying().StrongCharge())
    value/=double(abs(GetDecaying().StrongCharge()));
  value /= SymmetryFactor();
  return value;
}

void Decay_Channel::
GenerateKinematics(ATOOLS::Vec4D_Vector& momenta, bool anti,
		   METOOLS::Spin_Density* sigma,
		   const std::vector<ATOOLS::Particle*>& parts)
{
  static std::string mname(METHOD);
  Return_Value::IncCall(mname);
  if(momenta.size()==2) {
    momenta[1]=momenta[0];
    if (sigma) {
      if (p_amps) delete p_amps;
      p_amps=new Amplitude2_Tensor(parts, 0);
    }
    return;
  }
  double value(0.);
  int trials(0);
  do {
    if(trials>10000) {
      msg_Error()<<METHOD<<"("<<Name()<<"): "
                 <<"Rejected decay kinematics 10000 times. "
                 <<"This indicates a wrong maximum. "
                 <<"Will accept kinematics."
                 <<endl;
      Return_Value::IncRetryMethod(mname);
      break;
    }
    value = Differential(momenta,anti,sigma, parts);
    if(value/m_max>1.05 && m_max>1e-30) {
      if(value/m_max>1.3) {
        msg_Info()<<METHOD<<"("<<Name()<<") warning:"<<endl
                  <<"  d\\Gamma(x)="<<value<<" > max(d\\Gamma)="<<m_max
                  <<std::endl;
      }
      m_max=value;
      Return_Value::IncRetryMethod(mname);
      break;
    }
    trials++;
  } while( ran->Get() > value/m_max );
}

namespace ATOOLS {
  template <> Blob_Data<PHASIC::Decay_Channel*>::~Blob_Data() {}
  template class Blob_Data<PHASIC::Decay_Channel*>;
  template PHASIC::Decay_Channel*&Blob_Data_Base::Get<PHASIC::Decay_Channel*>();
}

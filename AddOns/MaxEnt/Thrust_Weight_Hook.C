#include "SHERPA/Tools/Userhook_Base.H"
#include "ATOOLS/Org/Message.H"
#include "SHERPA/Single_Events/Event_Handler.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include <algorithm>
#include "MODEL/Main/Running_AlphaS.H"
#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "CSSHOWER++/Main/CS_Shower.H"
#include "PDF/Main/Shower_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"

/*
 * counter-terms applied to the X+bb events in a fused sample.
 */

using namespace ATOOLS;
using namespace SHERPA;

class Thrust_Weight_Hook : public Userhook_Base, public Tag_Replacer {

private:
  MODEL::Running_AlphaS *p_as;
  Sherpa* p_sherpa;
  double m_thrust;
  std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;
public:

  Thrust_Weight_Hook(const Userhook_Arguments args) :
    Userhook_Base("Thrust_Weight"), p_as(MODEL::as),
    p_sherpa(args.p_sherpa)
  {
    msg_Debugging()<<"Thrust_Weight user hook active."<<std::endl;
    Settings& s = Settings::GetMainSettings();
    std::vector<std::string> params = s["THRUST_WEIGHTS"].
      SetDefault(std::vector<std::string>(1,"1")).GetVector<std::string>();
    for (size_t i(0);i<params.size();++i) {
      m_calcs.push_back(new Algebra_Interpreter());
      m_calcs.back()->SetTagReplacer(this);
      m_calcs.back()->AddTag("Tau","1.0");
      DEBUG_VAR(params[i]);
      m_calcs.back()->Interprete(params[i]);
      if (msg_LevelIsDebugging()) m_calcs.back()->PrintEquation();
    }
  }

  ~Thrust_Weight_Hook()
  {
    for (size_t i(0);i<m_calcs.size();++i) delete m_calcs[i];
  }

  inline void RotateMoms(vector<Vec3D> &p,const Vec3D &ref)
  {
    for(std::vector<Vec3D>::iterator
	  i=p.begin();i!=p.end();++i) *i=*i-ref*(ref**i);
  }

  static bool bigger(const ATOOLS::Vec3D &lhs,const ATOOLS::Vec3D &rhs)
  {
    return lhs.Sqr()>rhs.Sqr(); 
  }

  Vec3D NewAxis(const vector<Vec3D> &p,const Vec3D &ref)
  {
    Vec3D nextref = Vec3D(0.,0.,0.);
    int addsign;
    for (unsigned int i=0;i<p.size();++i) {
      addsign = 1;
      if (ref*p[i]<0.) addsign = -1;
      nextref = nextref+addsign*p[i];
    }
    return nextref/nextref.Abs();  
  }

  double SumP(const vector<Vec3D> &p)
  { 
    double sum_p = 0.;
    for (unsigned int i=0;i<p.size();i++) sum_p+=p[i].Abs();
    return sum_p;
  }

  double SumNP(const vector<Vec3D> &p,const Vec3D &n)
  { 
    double sum_np = 0.;
    for (unsigned int i=0;i<p.size();i++) sum_np+=dabs(p[i]*n);
    return sum_np;
  }

  void Calculate(const Blob_List *blobs)
  {
    std::vector<Vec3D> m_vectors, initialaxes;
    Particle_List fs(blobs->ExtractLooseParticles(1));
    for (Particle_List::const_iterator
	   pit(fs.begin());pit!=fs.end();++pit) {
      m_vectors.push_back((*pit)->Momentum());
    }
    m_thrust = 0.;
    Vec3D maxthrustaxis, lastaxis, curraxis, m_thrustaxis;
    double maxthrust=0., lastthrust , currthrust;
    unsigned int min_generators(std::min(4,(int)m_vectors.size()));
    int addsign;
    for (int pass=0; pass<2; pass++) {
      initialaxes.clear();
      if (pass==1) RotateMoms(m_vectors,m_thrustaxis);
      sort(m_vectors.begin(),m_vectors.end(),&bigger);
      for(unsigned int i=1;i<=intpow(2,min_generators-1);++i) {
	Vec3D axis;
	for(unsigned int j=1;j<=min_generators;++j) {
	  addsign = -1;
	  if (intpow(2,j)*((i+intpow(2,j-1)-1)/intpow(2,j)) >= i) addsign = 1;
	  axis = axis+addsign*m_vectors[j-1];
	}
	initialaxes.push_back(axis);
      }
      sort(initialaxes.begin(),initialaxes.end(), &bigger);
      for(unsigned int j=0;j<initialaxes.size();j++) 
	initialaxes[j] = initialaxes[j]/initialaxes[j].Abs();
      unsigned int ident = 0;
      double sump        = SumP(m_vectors);
      maxthrust          = 0.;
      for(unsigned int j=0; (j<initialaxes.size()) && (ident<2); j++) {
	curraxis         = initialaxes[j];
	currthrust       = SumNP(m_vectors,curraxis)/sump;
	lastthrust       = 0.;
	while (currthrust > lastthrust+1.e-4) {
	  lastthrust     = currthrust;
	  lastaxis       = curraxis;
	  curraxis       = NewAxis(m_vectors,curraxis);
	  currthrust     = SumNP(m_vectors,curraxis)/sump;
	}
	if (lastthrust < maxthrust-1.e-4) break;
	if (lastthrust > maxthrust+1.e-4) {
	  ident          = 0;
	  maxthrustaxis  = lastaxis;
	  maxthrust      = lastthrust;
	}
	ident++;
      }
      if (pass==0) { 
	m_thrustaxis = maxthrustaxis; 
	m_thrust     = maxthrust; 
      }
    }
  }

  std::string ReplaceTags(std::string &expr) const
  {
    return m_calcs.front()->ReplaceTags(expr);
  }

  Term *ReplaceTags(Term *term) const
  {
    switch (term->Id()) {
    case 1:
      term->Set(1.-m_thrust);
      return term;
    }
    return term;
  }

  void AssignId(Term *term)
  {
    if (term->Tag()=="Tau") term->SetId(1);
  }

  ATOOLS::Return_Value::code Run(ATOOLS::Blob_List* blobs)
  {
    DEBUG_FUNC(p_sherpa->GetInitHandler()->
	       GetMatrixElementHandler()->
	       Process()->Parent()->Name());
    Calculate(blobs);
    msg_Debugging()<<"T = "<<m_thrust<<", \\tau = "<<1-m_thrust<<"\n";
    if (m_thrust>1-0.0045) return Return_Value::Nothing;
    auto me_w_info = (*blobs->FindFirst(btp::Signal_Process))
      ["MEWeightInfo"]->Get<ME_Weight_Info*>();
    Weights_Map &wmap = (*blobs->FindFirst(btp::Signal_Process))
      ["WeightsMap"]->Get<Weights_Map>();
    for (size_t i(0);i<m_calcs.size();++i) {
      double w=m_calcs[i]->Calculate()->Get<double>();
      msg_Debugging()<<i<<": w = "<<w<<"\n";
      wmap *= w;
      *me_w_info *= w;
    }
    return Return_Value::Nothing;
  }

  void Finish() {}

};

DECLARE_GETTER(Thrust_Weight_Hook,"Thrust_Weight",
               Userhook_Base,Userhook_Arguments);

Userhook_Base *ATOOLS::Getter<Userhook_Base,Userhook_Arguments,Thrust_Weight_Hook>::
operator()(const Userhook_Arguments &args) const
{
  return new Thrust_Weight_Hook(args);
}

void ATOOLS::Getter<Userhook_Base,Userhook_Arguments,Thrust_Weight_Hook>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Thrust_Weight userhook";
}

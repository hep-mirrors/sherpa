#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"

#include "AddOns/OpenLoops/OpenLoops_Interface.H"
#include "AddOns/OpenLoops/GGH_KFactor_Setter.H"

#include "AMEGIC++/Main/Single_Process_External.H"
#include "AMEGIC++/Main/Single_Process.H"

#include <algorithm>

using namespace ATOOLS;

namespace PHASIC {

  typedef std::vector<size_t> Size_t_Vec;
  typedef std::vector<Size_t_Vec> SizeTVec_Vec;
  typedef std::map<Flavour, Size_t_Vec > Decay_Map;
  
  class GGH_Decay_KFactor_Setter : public GGH_KFactor_Setter {

  private:

    SizeTVec_Vec       m_index_vec;
    Decay_Map          m_decay_map;
    size_t             m_higgs_index;
    size_t             GetHiggsIndex() const;
    Decay_Map          GetDecayMap() const;
    SizeTVec_Vec       GetIndexVec() const;
    Vec4D_Vector       MappedMomenta(const std::vector<Vec4D>& p) const;
    void               Check() const;
    static inline bool CompareByValue(const Decay_Info* lhs, const Decay_Info* rhs){return (*lhs)<(*rhs);}

  public:

    GGH_Decay_KFactor_Setter(const KFactor_Setter_Arguments &args) : GGH_KFactor_Setter(args) { 
      m_decay_map   = GetDecayMap();  
      p_ampl        = GetAmpl();
      m_index_vec   = GetIndexVec();  
      p_ampl        = GetAmpl();
      m_higgs_index = GetHiggsIndex();
      Check();
    }
    ~GGH_Decay_KFactor_Setter();

  protected:

    double             MassCorrectionFactor(const Vec4D_Vector& p);
    Cluster_Amplitude* GetAmpl() const;
    Vec4D_Vector       GetMomenta() const;
    void               SetNextAmplitude();
  };
  
}

using namespace PHASIC;

GGH_Decay_KFactor_Setter::~GGH_Decay_KFactor_Setter() {
  if(p_next_ampl)
    p_next_ampl->Delete();
}

double GGH_Decay_KFactor_Setter::MassCorrectionFactor(const Vec4D_Vector& p)
{
  OpenLoops::s_interface->SetParameter("mass(25)",p[m_higgs_index].Mass());
  OpenLoops::s_interface->SetParameter("width(25)",0);
  return GGH_KFactor_Setter::MassCorrectionFactor(p);
}

void GGH_Decay_KFactor_Setter::Check() const
{
  if(p_proc->NIn()!=2)
    THROW(not_implemented, "Number of incoming flavours different from two");
  if(!ContainsDecays(*p_proc))
    THROW(fatal_error, "This KFactor is intended for hard ME decays, use the 'GGH' KFactor instead");
}

Vec4D_Vector GGH_Decay_KFactor_Setter::GetMomenta() const 
{
  return MappedMomenta(p_proc->Integrator()->Momenta());
}

void GGH_Decay_KFactor_Setter::SetNextAmplitude() {
  if(p_next_ampl)
    p_next_ampl->Delete();
  p_next_ampl = p_proc->ScaleSetter(1)->Amplitudes().front()->Next()->Copy();
  if(!p_next_ampl)
    THROW(fatal_error, "Scale setter returns invalid amplitude");
  DecayInfo_Vector& decays = p_next_ampl->Decays();
  // sorting might not be neccessary, seem to come sorted already
  std::sort(decays.begin(), decays.end(), GGH_Decay_KFactor_Setter::CompareByValue);
  for(DecayInfo_Vector::iterator it=decays.begin(); it!=decays.end(); ++it){
    ClusterLeg_Vector legs = p_next_ampl->Legs();
    size_t id = (**it).m_id;
    for(ClusterLeg_Vector::const_iterator lit=legs.begin(); lit!=legs.end(); ++lit){
      size_t ida = (*lit)->Id();
      if ( (ida & id) ){
	ClusterLeg_Vector::const_iterator ljt = lit; ++ljt;
	for(; ljt!=legs.end(); ++ljt){
	  size_t idb  = (*ljt)->Id();
	  if ( (ida | idb) == (id)){
	    p_next_ampl->CombineLegs((*lit),(*ljt),(*it)->m_fl);
	    lit=ljt=legs.end()-1; //break both loops
	  }
	}
      }
    }
  }
}

Cluster_Amplitude* GGH_Decay_KFactor_Setter::GetAmpl() const {
  Cluster_Amplitude* ret(Cluster_Amplitude::New());
  ret->SetNIn(2);
  ret->CreateLeg(Vec4D(), p_proc->Flavours()[0].Bar());
  ret->CreateLeg(Vec4D(), p_proc->Flavours()[1].Bar());
  for(std::vector<Subprocess_Info>::const_iterator it = p_proc->Info().m_fi.m_ps.begin(); it != p_proc->Info().m_fi.m_ps.end(); ++it)
    ret->CreateLeg(Vec4D(), it->m_fl);
  Process_Base::SortFlavours(ret);
  return ret;
}

// build a map: keys   = decaying flavours in final state
//            : values = indices of decay products as they appear in p_proc->Flavours
std::map<Flavour, Size_t_Vec > GGH_Decay_KFactor_Setter::GetDecayMap() const{
  std::map<Flavour, Size_t_Vec > map;
  Flavour_Vector fs_flavs = p_proc->Flavours();
  if(fs_flavs.size()<3)
    THROW(fatal_error, "Internal error, too few flavours in process");
  fs_flavs.erase(fs_flavs.begin()); fs_flavs.erase(fs_flavs.begin());
  for(std::vector<Subprocess_Info>::const_iterator it = p_proc->Info().m_fi.m_ps.begin(); it != p_proc->Info().m_fi.m_ps.end(); ++it)
    {
      Flavour_Vector decay_flavs = it->GetExternal();
      if(decay_flavs.size()<2) continue;
      if(!map.insert(std::make_pair(it->m_fl, Size_t_Vec())).second)
	THROW(not_implemented, "Cannot handle multiple decays of the same flavour");
      Size_t_Vec& indvec = map[it->m_fl];
      for(Flavour_Vector::const_iterator kt(decay_flavs.begin()); kt!=decay_flavs.end(); ++kt) {
	if( std::count(decay_flavs.begin(), decay_flavs.end(), *kt) != std::count(fs_flavs.begin(), fs_flavs.end(), *kt) )
	  THROW(not_implemented, "Cannot unambiguously map final state flavours to decay products");
	for(size_t i(0); i<fs_flavs.size(); i++) {
	  if( (fs_flavs[i] == *kt) && std::find(indvec.begin(), indvec.end(), i+2)==indvec.end() ) {
	    indvec.push_back(i+2);
	    break;
	  }
	}
      }
      if(indvec.size()!=decay_flavs.size())
	THROW(fatal_error, "Internal error");
    }
  return map;
}

size_t GGH_Decay_KFactor_Setter::GetHiggsIndex() const {
  for(size_t i(2); i<p_ampl->Legs().size(); i++){
    if (p_ampl->Leg(i)->Flav().Kfcode() == kf_h0)
      return i;
  }
  THROW(fatal_error, "Internal error");
}

std::vector<Size_t_Vec > GGH_Decay_KFactor_Setter::GetIndexVec() const {
  std::vector<std::pair<Flavour,size_t > > fs_flavs;
  for(size_t i(2); i<p_proc->Flavours().size(); i++)
    fs_flavs.push_back(std::make_pair(Flavour(p_proc->Flavours()[i]), i));
  SizeTVec_Vec ret(p_ampl->Legs().size());
  ret[0].push_back(0);
  ret[1].push_back(1);
  for(size_t i(2); i<p_ampl->Legs().size(); i++){
    Flavour flav = p_ampl->Leg(i)->Flav();
    std::map<Flavour, Size_t_Vec >::const_iterator it = (m_decay_map.find(flav));
    if(it!=m_decay_map.end())
      ret[i] = it->second;
    else{
      ret[i] = Size_t_Vec(1);
      for(std::vector<std::pair<Flavour,size_t> >::iterator jt = fs_flavs.begin();jt!=fs_flavs.end();++jt){
	if(jt->first == flav){
	  ret[i][0] = (jt->second);
	  fs_flavs.erase(jt);
	  break;
	}
      }
    }
  }
  for (SizeTVec_Vec::const_iterator it=ret.begin(); it!=ret.end(); ++it){
    for(Size_t_Vec::const_iterator jt=it->begin(); jt!=it->end(); ++jt){
      for(SizeTVec_Vec::const_iterator kt=++it; kt!=ret.end(); ++kt){
   	if(std::find(kt->begin(),kt->end(), *jt)!=kt->end())
   	  THROW(fatal_error, "Internal error");
      }
      --it;
    }
  }
  return ret;
}

Vec4D_Vector GGH_Decay_KFactor_Setter::MappedMomenta(const Vec4D_Vector& p) const {
  Vec4D_Vector ret(m_index_vec.size());
  ret[0] = p[0];  ret[1] = p[1];
  for(size_t i(2); i<m_index_vec.size(); i++){
    for(size_t j(0); j<m_index_vec[i].size(); j++){
      ret[i] += p[m_index_vec[i][j]];
    }
  }
  return ret;
}
 

DECLARE_GETTER(GGH_Decay_KFactor_Setter,
	       "GGH_Decay",
	       KFactor_Setter_Base,
	       KFactor_Setter_Arguments);

KFactor_Setter_Base *Getter<KFactor_Setter_Base,
					    KFactor_Setter_Arguments,
					    GGH_Decay_KFactor_Setter>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new GGH_Decay_KFactor_Setter(args);
}

void Getter<KFactor_Setter_Base,
		    KFactor_Setter_Arguments,
		    GGH_Decay_KFactor_Setter>::
PrintInfo(std::ostream &str,const size_t width) const 
{ str<<"GGH_Decay K-Factor\n"; }

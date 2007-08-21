#include "XS_QCD_CDBG.H"

#include "N_Parton_CDBG.H"
#include "Phase_Space_Handler.H"
#include "Color_Integrator.H"
#include "Helicity_Integrator.H"
#include "XS_Group.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {

  template <> 
  Single_XS *Single_XS::GetProcess<CDXS_pp_np>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (nqcd!=nout || nqed!=0) return NULL;
    std::vector<int> nquark(7,0);
    for (size_t i(0);i<nin+nout;++i) {
      if (!flavours[i].IsGluon()) {
	if (!flavours[i].IsQuark()) return NULL;
	int n(flavours[i].IsAnti()?1:-1);
	nquark[flavours[i].Kfcode()]+=i<2?-n:n;
      }
    }
    for (size_t i(1);i<=6;++i) if (nquark[i]!=0) return NULL;
    return new CDXS_pp_np(nin,nout,flavours,model); 
  }
  
}

CDXS_pp_np::CDXS_pp_np(const size_t &nin,const size_t &nout,
		       const ATOOLS::Flavour *const fl,
		       XS_Model_Base *const model): 
  Single_XS(nin,nout,fl,model), p_map(NULL)
{
  SetColorScheme(cls::sample);
  SetHelicityScheme(hls::sample);
  for (short unsigned int i(0);i<4;++i) 
    p_colours[i][0]=p_colours[i][1]=0;
  m_nstrong=m_nin+m_nout;
  std::vector<Flavour> flavs(nin+nout);
  for (size_t i(0);i<nin+nout;++i) flavs[i]=p_flavours[i];
  p_bg = new N_Parton_CDBG(m_nin,m_nout,flavs,model);
  m_p.resize(m_nvector);
}

CDXS_pp_np::~CDXS_pp_np()
{
  if (p_bg!=NULL) delete p_bg;
}

bool CDXS_pp_np::MapProcess()
{
  if (p_parent==NULL || p_parent==this) return false;
  std::vector<Flavour> mapfl(m_nin+m_nout);
  std::map<Flavour,Flavour> flmap;
  int ckf(1);
  for (size_t i(0);i<m_nin+m_nout;++i) {
    if (flmap.find(p_flavours[i])==flmap.end()) {
      if (p_flavours[i].IsQuark()) {
	Flavour qr(p_flavours[i].Kfcode());
	flmap[qr]=Flavour((kf::code)ckf++);
	flmap[qr.Bar()]=flmap[qr].Bar();
      }
      else {
	flmap[p_flavours[i]]=p_flavours[i];
      }
    }
    mapfl[i]=flmap[p_flavours[i]];
  }
  std::string mapname(GenerateName(m_nin,m_nout,&mapfl.front()));
  XS_Base *mp(((XS_Group*)p_parent)->Matching(mapname));
  if (mp==NULL || mp==this) return false;
  p_map=dynamic_cast<CDXS_pp_np*>(mp);
  if (p_map==NULL) return false;
  msg_Debugging()<<METHOD<<"(): Map '"<<m_name<<"' -> '"
		 <<p_map->Name()<<"'."<<std::endl;
  return true;
}

bool CDXS_pp_np::Tests()
{
  if (MapProcess()) {
    delete p_bg;
    p_bg=NULL;
    return true;
  }
  if (m_gpath.length()>0) p_bg->PrintGraphs(m_gpath);
  msg_Info()<<METHOD<<"(): Test '"<<m_name<<"' <- '"
	    <<p_activepshandler->Process()->Name()<<"'."<<std::endl;
  p_activepshandler->InitIntegrators();
  Helicity_Integrator *helint(p_activepshandler->HelicityIntegrator());
  if (helint!=NULL) {
    p_bg->SetHelicityIntegrator(helint);
    Flavour_Vector fl(m_nin+m_nout);
    for (size_t i(0);i<fl.size();++i) fl[i]=p_flavours[i];
    if (!helint->Construct(fl)) return false;
  }
  Color_Integrator *colint(p_activepshandler->ColorIntegrator());
  p_bg->SetColorIntegrator(colint);
  Idx_Vector ids(m_nin+m_nout,0);
  Int_Vector types(m_nin+m_nout,0), acts(ids.size(),1);
  for (size_t i(0);i<ids.size();++i) {
    ids[i]=i;
    if (p_flavours[i].IsGluon()) types[i]=0;
    else if (p_flavours[i].IsAnti()) types[i]=i<2?1:-1;
    else types[i]=i<2?-1:1;
  }
  if (!colint->ConstructRepresentations(ids,types,acts)) return false;
  p_activepshandler->TestPoint(p_momenta);
  for (size_t i(0);i<m_nvector;++i) m_p[i]=p_momenta[i];
  return p_bg->GaugeTest(m_p);
}

double CDXS_pp_np::operator()(const double s,const double t,const double u)
{
  if (p_map!=NULL) return p_map->LastXS();
  for (size_t i(0);i<m_nvector;++i) m_p[i]=p_momenta[i];
  m_p[0]=Vec4D(m_p[0][0],0.0,0.0,m_p[0][3]<0.0?-m_p[0][0]:m_p[0][0]);
  m_p[1]=Vec4D(m_p[1][0],0.0,0.0,m_p[1][3]<0.0?-m_p[1][0]:m_p[1][0]);
  return p_bg->Differential(m_p);
}

bool CDXS_pp_np::SetColours(const double s,const double t,const double u)
{ 
  return false;
}

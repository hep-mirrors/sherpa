#include "EXTRA_XS/Cluster/Cluster_Algorithm.H"

#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "EXTRA_XS/Main/Single_Process.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Message.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

Cluster_Algorithm::Cluster_Algorithm():
  p_ampl(NULL), p_clus(NULL) {}

Cluster_Algorithm::~Cluster_Algorithm()
{
}

bool Cluster_Algorithm::Cluster(Single_Process *const xs)
{
  Selector_Base *jf=xs->Selector()
    ->GetSelector("Jetfinder");
  ME2_Base *me(xs->GetME());
  bool swap(xs->Integrator()->InSwaped());
  if (me==NULL) THROW(not_implemented,"Non-ME-specified process");
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetJF(jf);
  const Vec4D_Vector &moms(xs->Integrator()->Momenta());
  me->SetColours(moms);
  PHASIC::Process_Base *pb(xs->IsMapped()?xs->MapProc():xs);
  double muf2(pb->ScaleSetter()->Scale(stp::fac));
  double mur2(pb->ScaleSetter()->Scale(stp::ren));
  for (size_t i(0);i<xs->NIn()+xs->NOut();++i) {
    size_t id(1<<p_ampl->Legs().size());
    size_t idx(i<2?(swap?1-i:i):i);
    ColorID col(me->Colours()[idx][0],me->Colours()[idx][1]);
    if (i<2) col=col.Conj();
    Flavour flav(i<2?xs->Flavours()[i].Bar():
		 xs->Flavours()[i]);
    Vec4D mom(i<2?-moms[i]:moms[i]);
    p_ampl->CreateLeg(mom,flav,col,id);
    p_ampl->Legs().back()->SetStat(1);
    p_ampl->Legs().back()->SetNMax
      (xs->Info().m_fi.NMaxExternal());
  }
  // set colour partners
  p_ampl->SetNIn(xs->NIn());
  p_ampl->SetMuR2(mur2);
  p_ampl->SetMuF2(muf2);
  p_ampl->SetProcs(xs);
  double kt2(xs->Shower()->GetClusterDefinitions()->CoreScale(p_ampl));
  p_ampl->SetKT2(kt2);
  p_ampl->SetMu2(kt2);
  p_ampl->SetOrderEW(xs->OrderEW());
  p_ampl->SetOrderQCD(xs->OrderQCD());
  msg_Debugging()<<*p_ampl<<"\n";
  return true;
}


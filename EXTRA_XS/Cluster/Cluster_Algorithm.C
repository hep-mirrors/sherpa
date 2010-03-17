#include "EXTRA_XS/Cluster/Cluster_Algorithm.H"

#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "EXTRA_XS/Main/Single_Process.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
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
  bool swap(xs->Integrator()->InSwaped()), trig(true);
  if (jf) {
    Vec4D_Vector moms(xs->Integrator()->Momenta());
    if (swap) {
      std::swap<Vec4D>(moms[0],moms[1]);
      for (size_t i(0);i<moms.size();++i)
	moms[i]=Vec4D(moms[i][0],-moms[i]);
    }
    trig=jf->Trigger(moms);
  }
  msg_Debugging()<<METHOD<<"(): trig = "<<trig<<"\n";
  if (me==NULL) THROW(not_implemented,"Non-ME-specified process");
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetJF(jf);
  const Vec4D_Vector &moms(xs->Integrator()->Momenta());
  me->SetColours(moms);
  double muf2(xs->ScaleSetter()->Scale(stp::fac));
  double mur2(xs->ScaleSetter()->Scale(stp::ren));
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
      (trig?xs->Info().m_fi.NMaxExternal():
       xs->Info().m_fi.NExternal());
  }
  // set colour partners
  p_ampl->SetNIn(xs->NIn());
  p_ampl->SetMuR2(mur2);
  p_ampl->SetMuF2(muf2);
  double kt2(0.0);
  if (p_ampl->Leg(0)->Flav().Resummed() && 
      p_ampl->Leg(1)->Flav().Resummed() &&
      p_ampl->Leg(2)->Flav().Resummed() &&
      p_ampl->Leg(3)->Flav().Resummed()) {
    kt2=Max(moms[2].MPerp2(),moms[3].MPerp2());
  }
  else {
    int sintt(xs->GetME()->SIntType());
    kt2=std::numeric_limits<double>::max();
    if (sintt&1) kt2=Min(kt2,(moms[0]+moms[1]).Abs2());
    if (sintt&2) kt2=Min(kt2,dabs((moms[0]-moms[2]).Abs2()));
    if (sintt&4) kt2=Min(kt2,dabs((moms[0]-moms[3]).Abs2()));
  }
  p_ampl->SetKT2QCD(kt2);
  p_ampl->SetX1(xs->Integrator()->ISR()->X1());
  p_ampl->SetX2(xs->Integrator()->ISR()->X2());
  p_ampl->SetOrderEW(xs->OrderEW());
  p_ampl->SetOrderQCD(xs->OrderQCD());
  if (msg_LevelIsDebugging()) p_ampl->Print();
  return true;
}


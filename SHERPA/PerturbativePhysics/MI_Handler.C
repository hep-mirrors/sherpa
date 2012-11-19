#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "EXTRA_XS/Main/Single_Process.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"

#ifdef USING__Amisic
#include "AMISIC++/Main/Amisic.H"
#endif

using namespace SHERPA;
using namespace ATOOLS;

MI_Handler::MI_Handler(std::string path,std::string file,MODEL::Model_Base *model,
		       BEAM::Beam_Spectra_Handler *beam,PDF::ISR_Handler *isr) :
  p_beam(beam),
  p_isr(isr),
#ifdef USING__Amisic
  p_amisic(NULL),
#endif
  p_ampl(NULL),
  m_type(None),
  m_ycut(1.0e-7)
{
#ifdef USING__Amisic
  std::string mihandler="None";
  ATOOLS::Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.AddWordSeparator("\t");
  read.SetInputPath(path);
  read.SetInputFile(file);
  mihandler=read.GetValue<std::string>("MI_HANDLER","Amisic");
  path+=read.GetValue<std::string>("INPUT_PATH","");
  file=read.GetValue<std::string>("INPUT_FILE",file);
  if (!ATOOLS::rpa->gen.Beam1().IsHadron() ||
      !ATOOLS::rpa->gen.Beam2().IsHadron()) mihandler="None";
  if (mihandler==std::string("Amisic")) {
    p_amisic = new AMISIC::Amisic(model,beam,isr);
    p_amisic->SetInputPath(path);
    p_amisic->SetOutputPath(ATOOLS::rpa->gen.Variable("SHERPA_RUN_PATH")+"/");
    p_amisic->SetInputFile(file);
    if (!p_amisic->Initialize()) {
      msg_Error()<<METHOD<<"(): Cannot initialize MPI generator. "
		 <<"Continue without."<<std::endl;
      delete p_amisic;
      p_amisic=NULL;
      return;
    }
    m_ycut=p_amisic->HardBase()->Stop(0);
    m_ycut=ATOOLS::sqr(m_ycut/ATOOLS::rpa->gen.Ecms());
    m_type=Amisic;
  }
#endif
}

MI_Handler::~MI_Handler() 
{
#ifdef USING__Amisic
  if (p_amisic!=NULL) delete p_amisic;
#endif
}

bool MI_Handler::GenerateHardProcess(ATOOLS::Blob *blob)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return p_amisic->GenerateHardProcess(blob);
#endif
  default    : break;
  }
  return false;
}

bool MI_Handler::GenerateSoftProcess(ATOOLS::Blob *blob)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return p_amisic->GenerateSoftProcess(blob);
#endif
  default    : break;
  }
  return false;
}

bool MI_Handler::GenerateEvent(ATOOLS::Blob_List *bloblist)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: break; // p_amisic->GenerateEvent(bloblist);
#endif
  default    : break;
  }
  return false;
}

bool MI_Handler::VetoHardProcess(ATOOLS::Blob *blob)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return p_amisic->VetoHardProcess(blob);
#endif
  default    : break;
  }
  return false;
}

void MI_Handler::SetScaleMin(double scalemin,unsigned int i)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic:
    p_amisic->HardBase()->SetStop(scalemin,i);
    p_amisic->SoftBase()->SetStart(scalemin,i);
    break;
#endif
  default:
    break;
  }
}

void MI_Handler::SetScaleMax(double scalemax,unsigned int i)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic:
    p_amisic->HardBase()->SetStart(scalemax,i);
    break;
#endif
  default:
    break;
  }
}

double MI_Handler::ScaleMin(unsigned int i)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return p_amisic->HardBase()->Stop(i);
#endif
  default    : break;
  }
  return 0.;
}

double MI_Handler::ScaleMax(unsigned int i)
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return p_amisic->HardBase()->Start(i);
#endif
  default    : break;
  }
  return 0.;
}

void MI_Handler::Reset()
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: p_amisic->Reset();
#endif
  default    : break;
  }
}

void MI_Handler::CleanUp()
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: p_amisic->CleanUp();
#endif
  default    : break;
  }
}

std::string MI_Handler::MIGenerator() 
{
  return Name();
}

MI_Handler::TypeID MI_Handler::Type() 
{
  return m_type;
}

std::string MI_Handler::Name() 
{
  switch (m_type) {
#ifdef USING__Amisic
  case Amisic: return std::string("Amisic");
#endif
  case None  : return std::string("None");
  default    : break;
  }
  return std::string("Unknown");
}

PDF::ISR_Handler *MI_Handler::ISRHandler()
{
  return p_isr;
}

ATOOLS::Cluster_Amplitude *MI_Handler::ClusterConfiguration()
{
#ifdef USING__Amisic
  PHASIC::Process_Base *xs(p_amisic->HardBase()->XS());
  if (xs->Get<EXTRAXS::Single_Process>()==NULL) return NULL;
  EXTRAXS::ME2_Base *me(xs->Get<EXTRAXS::Single_Process>()->GetME());
  if (me==NULL) THROW(fatal_error,"Cannot handle non-generic ME's.");
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  p_ampl = Cluster_Amplitude::New();
  const Vec4D_Vector &moms(xs->Integrator()->Momenta());
  me->SetColours(moms);
  double muf2(xs->ScaleSetter()->Scale(stp::fac));
  double mur2(xs->ScaleSetter()->Scale(stp::ren));
  for (size_t i(0);i<xs->NIn()+xs->NOut();++i) {
    size_t id(1<<p_ampl->Legs().size());
    size_t idx(i);
    ColorID col(me->Colours()[idx][0],me->Colours()[idx][1]);
    if (i<2) col=col.Conj();
    Flavour flav(i<2?xs->Flavours()[i].Bar():
		 xs->Flavours()[i]);
    Vec4D mom(i<2?-moms[i]:moms[i]);
    p_ampl->CreateLeg(mom,flav,col,id);
    p_ampl->Legs().back()->SetStat(1);
    p_ampl->Legs().back()->SetNMax(xs->Info().m_fi.m_nmax);
  }
  // set colour partners
  p_ampl->SetNIn(xs->NIn());
  p_ampl->SetMuR2(mur2);
  p_ampl->SetMuF2(muf2);
  p_ampl->SetKT2(muf2);
  p_ampl->SetMu2(mur2);
  p_ampl->SetOrderEW(xs->OrderEW());
  p_ampl->SetOrderQCD(xs->OrderQCD());
  p_ampl->SetMS(this);
  msg_Debugging()<<*p_ampl<<"\n";
#endif
  return p_ampl;
}

double MI_Handler::Mass(const ATOOLS::Flavour &fl) const
{
  return fl.Mass();
}

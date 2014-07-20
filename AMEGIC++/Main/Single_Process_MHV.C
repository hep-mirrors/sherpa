#include "AMEGIC++/Main/Single_Process_MHV.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "AMEGIC++/Phasespace/Phase_Space_Generator.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace AMEGIC;
using namespace MODEL;
using namespace PHASIC;
using namespace PDF;
using namespace BEAM;
using namespace ATOOLS;
using namespace std;

/*-------------------------------------------------------------------------------

  Constructors

  ------------------------------------------------------------------------------- */

AMEGIC::Single_Process_MHV::Single_Process_MHV():
  m_gen_str(2), m_ownamps(false), p_hel(0), p_BS(0), p_ampl(0), p_shand(0), p_psgen(0),  
  p_MHVamp(0), p_momlist(0), p_partner(this)
{
}

AMEGIC::Single_Process_MHV::~Single_Process_MHV()
{
  if (p_hel)      {delete p_hel; p_hel=0;}
  if (p_BS)       {delete p_BS;   p_BS=0;}
  if (p_shand)    {delete p_shand;p_shand=0;}
  if (p_ampl)     {delete p_ampl; p_ampl=0;}
  if (p_psgen)    {delete p_psgen; p_psgen=0;}
  if (p_MHVamp&&m_ownamps)     {delete p_MHVamp; p_MHVamp=0;}
#ifndef Basic_Sfuncs_In_MHV 
  if (p_momlist)    {delete p_momlist; p_momlist=0;}
#endif
}

/*------------------------------------------------------------------------------
  
  Generic stuff for initialization of Single_Process_MHVes
      
  ------------------------------------------------------------------------------*/

void AMEGIC::Single_Process_MHV::PolarizationNorm() {
  m_Norm = SymmetryFactors() * m_pol.Spin_Average(m_nin,&m_flavs.front());
}


int AMEGIC::Single_Process_MHV::InitAmplitude(Model_Base * model,Topology* top,
					 vector<Process_Base *> & links,
					 vector<Process_Base *> & errs)
{
  Init();
  model->GetCouplings(m_cpls);
  if (!model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;
  m_newlib   = false;
  m_libnumb  = 0;
  m_pslibname = m_libname = ToString(m_nin)+"_"+ToString(m_nout);
  if (m_gen_str>1) m_ptypename = "P"+m_libname;
  else m_ptypename = "N"+m_libname;
  PolarizationNorm();
  if (m_gen_str>1) {
    ATOOLS::MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename);
  }
  string newpath=rpa->gen.Variable("SHERPA_CPP_PATH");
  ATOOLS::MakeDir(newpath);

  p_hel    = new Helicity(m_nin,m_nout,&m_flavs.front(),p_pl);
  p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nin+m_nout,&m_flavs.front(),p_b);  
  p_BS->Setk0(s_gauge);

  //////////////////////////////////////////////// 
#ifdef Basic_Sfuncs_In_MHV
  p_momlist = p_BS; 
#else
  p_momlist = new MomentumList(m_nin,m_nout); 
#endif 

  Flavour* fl = &m_flavs.front();
  int *plist = new int[m_nin+m_nout];
  for (size_t i=0;i<m_nin;i++) plist[i] = fl[i];
  for (size_t i=m_nin;i<m_nin+m_nout;i++) plist[i]=-fl[i];
  p_MHVamp = FullAmplitude_MHV_Handler(model,&m_cpls,m_nin+m_nout,plist,p_momlist,m_ownamps); 

  delete [] plist;
  //////////////////////////////////////////////

  p_shand  = new String_Handler(m_gen_str,p_BS,model->GetVertex()->GetCouplings());
  int oew(m_oew), oqcd(m_oqcd), ntchanmin(m_ntchanmin);
  p_ampl   = new Amplitude_Handler(m_nin+m_nout,&m_flavs.front(),p_b,p_pinfo,model,top,oqcd,oew,ntchanmin,
				   &m_cpls,p_BS,p_shand,m_print_graphs,0);
  m_oew=oew;
  m_oqcd=oqcd;
  if (p_ampl->GetGraphNumber()==0) {
    msg_Tracking()<<"AMEGIC::Single_Process_MHV::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return 0;
  }
  map<string,Complex> cplmap;
  for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
    cplmap.clear();
    if (FlavCompare(links[j]) && p_ampl->CompareAmplitudes(links[j]->GetAmplitudeHandler(),m_sfactor,cplmap)) {
      if (p_hel->Compare(links[j]->GetHelicity(),m_nin+m_nout)) {
	m_sfactor = sqr(m_sfactor);
	msg_Tracking()<<"AMEGIC::Single_Process_MHV::InitAmplitude : Found compatible process for "<<Name()<<" : "<<links[j]->Name()<<endl;
	  
	p_mapproc = p_partner = (Single_Process_MHV*)links[j];
	m_iresult = p_partner->Result()*m_sfactor;

	InitFlavmap(p_partner);
	
	Minimize();
	return 1;
      }
    }
  }

  p_ampl->FillPointlist();
  p_BS->Initialize();


  int result(Tests());
  switch (result) {
  case 1 :
    if (p_partner==this) links.push_back(this);
    msg_Info()<<".";
    
    if (m_gen_str<2) return 1;
    if (p_partner==this && Result()>0.) SetUpIntegrator();
    return 1;
  case -3: return 0;
  default :
    msg_Error()<<"ERROR in AMEGIC::Single_Process_MHV::InitAmplitude : "<<std::endl
	       <<"   Failed for "<<m_name<<"."<<endl;
    errs.push_back(this);
    return 0;
  }

  return 1;
}



int AMEGIC::Single_Process_MHV::Tests() 
{
  int number      = 1;
  int gauge_test  = 1;

  /* ---------------------------------------------------
     
     The reference result for momenta moms

     --------------------------------------------------- */

  string testname = string("");

  double M2 = 0.;
  double helvalue;
  if (gauge_test) {
#ifdef Basic_Sfuncs_In_MHV
    p_BS->Setk0(0);
    p_BS->CalcEtaMu(p_testmoms); 
#else
    p_momlist->PutMomenta(p_testmoms);
#endif    
 
    msg_Tracking()<<"AMEGIC::Single_Process_MHV::Tests for "<<m_name<<std::endl
		  <<"   Prepare gauge test and init helicity amplitudes. This may take some time."
	      <<std::endl;
    for (size_t i=0;i<p_hel->MaxHel();i++) { 
      if (p_hel->On(i)) {
	helvalue = p_MHVamp->MSquare((*p_hel)[i],p_BS)*p_hel->PolarizationFactor(i); 
	M2      +=  helvalue;
      } 
    }
    M2     *= p_MHVamp->ParticlesNorm();
    m_iresult  = M2;
  }
  /* ---------------------------------------------------
     
  First test : gauge test
  
  --------------------------------------------------- */
#ifdef Basic_Sfuncs_In_MHV
  p_BS->Setk0(s_gauge);
  p_BS->CalcEtaMu(p_testmoms); 
#else
  p_momlist->PutMomenta(p_testmoms);
#endif  
  number++;

  double M2g = 0.;
  double * M_doub = new double[p_hel->MaxHel()];
  for (size_t i=0; i<p_hel->MaxHel(); ++i) { 
    if (p_hel->On(i)) {
      M_doub[i]  = p_MHVamp->MSquare((*p_hel)[i],p_BS)*p_hel->PolarizationFactor(i);
      M2g       += M_doub[i];
    }
  }

  //shorten helicities
  int switchhit = 0;
  for (size_t i=0;i<p_hel->MaxHel();i++) {
    if (M_doub[i]==0. || dabs(M_doub[i]/M2g)<(ATOOLS::Accu()*1.e-2)) {
      p_hel->SwitchOff(i);
      switchhit++;
    }
  }
  msg_Tracking()<<"AMEGIC::Single_Process_MHV::Tests for "<<m_name<<std::endl
		<<"   Switched off or mapped "<<switchhit<<" helicities."<<std::endl;

  M2g    *= p_MHVamp->ParticlesNorm();
  m_iresult  = M2g;
  p_BS->StartPrecalc();

  if (gauge_test) {
    if (!ATOOLS::IsEqual(M2,M2g)) {
      msg_Out()<<"WARNING:  Gauge test not satisfied: "
	       <<M2<<" vs. "<<M2g<<" : "<<dabs(M2/M2g-1.)*100.<<"%"<<endl
	       <<"Gauge(1): "<<abs(M2)<<endl
	       <<"Gauge(2): "<<abs(M2g)<<endl;
    }
  }

  for (size_t i=0;i<p_hel->MaxHel();i++) {
    if (p_hel->On(i)) {
      for (size_t j=i+1;j<p_hel->MaxHel();j++) {
	if (p_hel->On(j)) {
#ifdef FuckUp_Helicity_Mapping
	  if (ATOOLS::IsEqual(M_doub[i],M_doub[j])) {
	    p_hel->SwitchOff(j);
	    p_hel->SetPartner(i,j);
	    p_hel->IncMultiplicity(i);
	  }
#endif
	}
      }
    }
  }
  
  delete[] M_doub;
  return 1;
}

int AMEGIC::Single_Process_MHV::PerformTests()
{
  return 1;
}

/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/

bool AMEGIC::Single_Process_MHV::FillIntegrator
(PHASIC::Phase_Space_Handler *const psh)
{
  if (!SetUpIntegrator()) THROW(fatal_error,"No integrator");
  return Process_Base::FillIntegrator(psh);
}

bool AMEGIC::Single_Process_MHV::SetUpIntegrator() 
{  
  if (m_nin==2) {
    if ( (m_flavs[0].Mass() != p_int->ISR()->Flav(0).Mass()) ||
	 (m_flavs[1].Mass() != p_int->ISR()->Flav(1).Mass()) ) p_int->ISR()->SetPartonMasses(m_flavs);
    if (CreateChannelLibrary()) return 1;
  }
  if (m_nin==1) if (CreateChannelLibrary()) return 1;
  m_newlib=true;
  return 0;
}

bool AMEGIC::Single_Process_MHV::CreateChannelLibrary()
{
  p_psgen     = new Phase_Space_Generator(m_nin,m_nout);
  bool newch  = 0;
  if (m_nin>=1)  newch = p_psgen->Construct(p_channellibnames,m_ptypename,m_pslibname,&m_flavs.front(),this); 
  if (newch>0) return 0;
  return 1;
}

/*------------------------------------------------------------------------------
  
  Process management
  
  ------------------------------------------------------------------------------*/
void AMEGIC::Single_Process_MHV::Minimize()
{
  if (p_partner==this) return;
  if (p_hel)      {delete p_hel; p_hel=0;}
  if (p_BS)       {delete p_BS;   p_BS=0;}
  if (p_shand)    {delete p_shand;p_shand=0;}
  if (p_ampl)     {delete p_ampl; p_ampl=0;}
  if (p_psgen)    {delete p_psgen; p_psgen=0;}

  m_oqcd      = p_partner->OrderQCD();
  m_oew       = p_partner->OrderEW();
}

double AMEGIC::Single_Process_MHV::Partonic(const Vec4D_Vector &moms,const int mode) 
{ 
  if (mode==1) return m_lastxs;
  if (!Selector()->Result()) return m_lastxs = 0.0;
  if (!(IsMapped() && LookUp())) {
    p_partner->ScaleSetter()->CalculateScale(moms,m_cmode);
  }
  return DSigma(moms,m_lookup); 
}

double AMEGIC::Single_Process_MHV::DSigma(const ATOOLS::Vec4D_Vector &_moms,bool lookup)
{
  m_lastxs = 0.;
  if (p_partner == this) {
    m_lastxs = m_Norm * operator()((ATOOLS::Vec4D*)&_moms.front());
  }
  else {
    if (lookup && p_partner->m_lookup) m_lastxs = p_partner->LastXS()*m_sfactor;
    else m_lastxs = m_Norm * p_partner->operator()((ATOOLS::Vec4D*)&_moms.front())*m_sfactor;
  }
  return m_lastxs;
}

double AMEGIC::Single_Process_MHV::operator()(const ATOOLS::Vec4D* mom)
{
  double M2(0.);

#ifdef Basic_Sfuncs_In_MHV
  p_BS->CalcEtaMu((ATOOLS::Vec4D*)mom); 
#else
  p_momlist->PutMomenta(mom);
#endif  

  double helvalue;
  for (size_t i=0; i<p_hel->MaxHel(); ++i) { 
      if (p_hel->On(i)) { 
	  helvalue = p_MHVamp->MSquare((*p_hel)[i],p_BS) * p_hel->Multiplicity(i) * p_hel->PolarizationFactor(i);  
	  M2       += helvalue;
      }
  }
  return M2*p_MHVamp->ParticlesNorm() * KFactor();
}


void AMEGIC::Single_Process_MHV::FillAmplitudes(vector<METOOLS::Spin_Amplitudes>& amps,
                                                std::vector<std::vector<Complex> >& cols)
{
  if (p_partner==this) p_ampl->FillAmplitudes(amps, cols, p_hel, 1.0);
  else p_partner->FillAmplitudes(amps, cols, sqrt(m_sfactor));
}

void AMEGIC::Single_Process_MHV::FillAmplitudes(vector<METOOLS::Spin_Amplitudes>& amps,
                                                std::vector<std::vector<Complex> >& cols,
                                                double sfactor)
{
  if (p_partner==this) p_ampl->FillAmplitudes(amps, cols, p_hel, sfactor);
  else p_partner->FillAmplitudes(amps, cols, sfactor*sqrt(m_sfactor));
}

int AMEGIC::Single_Process_MHV::NumberOfDiagrams() { 
  if (p_partner==this) return p_ampl->GetGraphNumber(); 
  return p_partner->NumberOfDiagrams();
}

Point * AMEGIC::Single_Process_MHV::Diagram(int i) { 
  if (p_partner==this) return p_ampl->GetPointlist(i); 
  return p_partner->Diagram(i);
} 


void AMEGIC::Single_Process_MHV::AddChannels(std::list<std::string>* tlist) 
{
  if (p_partner==this) {    
    list<string>* clist = p_channellibnames;
    for (list<string>::iterator it=clist->begin();it!=clist->end();++it) {
      bool hit = 0;
      for (list<string>::iterator jt=tlist->begin();jt!=tlist->end();++jt) {
	if ((*it)==(*jt)) {
	  hit = 1;
	  break;
	}
      }
      if (!hit) tlist->push_back((*it));
    }
  }
}

void AMEGIC::Single_Process_MHV::FillCombinations
(Point *const p,size_t &id)
{
  if (p->middle) return;
  if (p->left==NULL || p->right==NULL) {
    id=1<<p->number;
    return;
  }
  size_t ida(id), idb(id);
  FillCombinations(p->left,ida);
  FillCombinations(p->right,idb);
  id=ida+idb;
  size_t idc((1<<(m_nin+m_nout))-1-id);
#ifdef DEBUG__Fill_Combinations
  msg_Debugging()<<"  comb "<<ID(ida)
		 <<" "<<ID(idb)<<" "<<ID(idc)<<"\n";
#endif
  m_ccombs.insert(std::pair<size_t,size_t>(ida,idb));
  m_ccombs.insert(std::pair<size_t,size_t>(idb,ida));
  m_ccombs.insert(std::pair<size_t,size_t>(idb,idc));
  m_ccombs.insert(std::pair<size_t,size_t>(idc,idb));
  m_ccombs.insert(std::pair<size_t,size_t>(idc,ida));
  m_ccombs.insert(std::pair<size_t,size_t>(ida,idc));
  if (idc!=1) {
    bool in(false);
    Flavour fl(ReMap(p->fl,p->GetPropID()));
    Flavour_Vector cf(m_cflavs[id]);
    for (size_t i(0);i<cf.size();++i)
      if (cf[i]==fl) {
	in=true;
	break;
      }
    if (!in) {
      m_cflavs[idc].push_back(fl.Bar());
      m_cflavs[id].push_back(fl);
#ifdef DEBUG__Fill_Combinations
      msg_Debugging()<<"  flav "<<ID(idc)<<" / "
		     <<ID(id)<<" -> "<<fl<<"\n";
#endif
    }
  }
}

void AMEGIC::Single_Process_MHV::FillCombinations()
{
#ifdef DEBUG__Fill_Combinations
  msg_Debugging()<<METHOD<<"(): '"<<m_name<<"' {\n";
#endif
  size_t nd(p_partner->NumberOfDiagrams());
  for (size_t i(0);i<nd;++i) {
    Point *p(p_partner->Diagram(i));
    size_t id(1<<p->number);
    FillCombinations(p,id);
  }
#ifdef DEBUG__Fill_Combinations
  msg_Debugging()<<"  } -> "<<m_cflavs.size()
		 <<" flavours, "<<m_ccombs.size()
		 <<" combinations\n";
  msg_Debugging()<<"}\n";
#endif
}

bool AMEGIC::Single_Process_MHV::Combinable
(const size_t &idi,const size_t &idj)
{
  if (m_ccombs.empty()) FillCombinations();
  Combination_Set::const_iterator 
    cit(m_ccombs.find(std::pair<size_t,size_t>(idi,idj)));
  return cit!=m_ccombs.end();
}

const Flavour_Vector &AMEGIC::Single_Process_MHV::
CombinedFlavour(const size_t &idij)
{
  if (m_cflavs.empty()) FillCombinations();
  CFlavVector_Map::const_iterator fit(m_cflavs.find(idij));
  if (fit==m_cflavs.end()) THROW(fatal_error,"Invalid request");
  return fit->second;
}

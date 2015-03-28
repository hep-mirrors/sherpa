#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_MHV.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_MHV.H"
#include "AMEGIC++/Amplitude/FullAmplitude_MHV_Base.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "AMEGIC++/Phasespace/Phase_Space_Generator.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"

#include <unistd.h>

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


Single_LOProcess_MHV::Single_LOProcess_MHV(const Process_Info &pi,
                                           BEAM::Beam_Spectra_Handler *const beam,
                                           PDF::ISR_Handler *const isr) :  
  Single_LOProcess(pi, beam, isr)
{
  m_ownamps = false;
  m_emitgluon = false;
}


Single_LOProcess_MHV::~Single_LOProcess_MHV()
{
  if (p_hel)      {delete p_hel; p_hel=0;}
  if (p_BS)       {delete p_BS;   p_BS=0;}
  if (p_shand)    {delete p_shand;p_shand=0;}
  if (p_ampl)     {delete p_ampl; p_ampl=0;}
#ifndef Basic_Sfuncs_In_MHV 
  if (p_momlist)    {delete p_momlist; p_momlist=0;}
#endif
  if (p_MHVamp && m_ownamps) {delete p_MHVamp; p_MHVamp=0;}
}


/*------------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ------------------------------------------------------------------------------*/



int Single_LOProcess_MHV::InitAmplitude(Amegic_Model * model,Topology* top,
					vector<Process_Base *> & links,
					vector<Process_Base *> & errs,int checkloopmap)
{
  m_type = 21;
  if (!model->p_model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;
  model->p_model->GetCouplings(m_cpls);
  
  m_partonlist.clear();
  for (size_t i=0;i<m_nin;i++) if (m_flavs[i].Strong()) m_partonlist.push_back(i);
  for (size_t i=m_nin;i<m_nin+m_nout;i++) if (m_flavs[i].Strong()) m_partonlist.push_back(i);

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
  p_MHVamp = FullAmplitude_MHV_Handler(model->p_model,&m_cpls,m_nin+m_nout,plist,p_momlist,m_ownamps,127,127); 

  delete [] plist;
  //////////////////////////////////////////////

  p_shand  = new String_Handler(m_gen_str,p_BS,model->p_model->GetCouplings());

  int ntchanmin(m_ntchanmin);
  p_ampl   = new Amplitude_Handler(m_nin+m_nout,fl,p_b,p_pinfo,model,top,m_maxcpl,m_mincpl,ntchanmin,
                                   &m_cpls,p_BS,p_shand,m_print_graphs,0,true);
  m_ntchanmin=ntchanmin;
  if (p_ampl->GetGraphNumber()==0) {
    msg_Tracking()<<"Single_LOProcess_MHV::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return 0;
  }
  map<string,Complex> cplmap;
  for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
    cplmap.clear();
    if (FlavCompare(links[j]) && p_ampl->CompareAmplitudes(links[j]->GetAmplitudeHandler(),m_sfactor,cplmap)) {
      if (p_hel->Compare(links[j]->GetHelicity(),m_nin+m_nout)) {
	m_sfactor = sqr(m_sfactor);
	msg_Tracking()<<"AMEGIC::Single_Process_MHV::InitAmplitude : Found compatible process for "<<Name()<<" : "<<links[j]->Name()<<endl;
	  
	p_partner = (Single_LOProcess_MHV*)links[j];
	m_iresult = p_partner->Result()*m_sfactor;
	InitFlavmap(p_partner);

	Minimize();
	return 1;
      }
    }
  }

  p_ampl->FillPointlist();
  p_BS->Initialize();

  switch (Tests()) {
  case 1 :
    if (p_partner==this) links.push_back(this);
    msg_Info()<<".";
    
    return 1;
  case -3: return 0;
  default :
    msg_Error()<<"ERROR in Single_Fin_Process_MHV::InitAmplitude : "<<std::endl
	       <<"   Failed for "<<m_name<<"."<<endl;
    errs.push_back(this);
    return 0;
  }
  return 1;
}


int Single_LOProcess_MHV::InitAmplitude(Amegic_Model * model,Topology* top,
					vector<Process_Base *> & links,
					vector<Process_Base *> & errs,
					std::vector<ATOOLS::Vec4D>* epol,std::vector<double> * pfactors)
{
  m_type = 11;
  if (!model->p_model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;
  model->p_model->GetCouplings(m_cpls);
  
  if (m_gen_str>1) {
    ATOOLS::MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename);
  }
  string newpath=rpa->gen.Variable("SHERPA_CPP_PATH");
  ATOOLS::MakeDir(newpath);
  int cnt=0;

  vector<int> fi_tags;
  m_pinfo.m_fi.GetTags(fi_tags);
  if (fi_tags.size()!=m_nout) THROW(fatal_error, "Internal error.");
  for (size_t i(0);i<m_pinfo.m_ii.m_ps.size();++i) {
    if (m_pinfo.m_ii.m_ps[i].m_tag==-1) {
      m_emit=i;
      cnt++;
    }
    if (m_pinfo.m_ii.m_ps[i].m_tag==-2) {
      m_spect=i;
      cnt+=10;
    }
  }
  for (size_t i(0);i<m_nout;++i) {
    if (fi_tags[i]==-1) {
      m_emit=i+NIn();
      cnt++;
    }
    if (fi_tags[i]==-2) {
      m_spect=i+NIn();
      cnt+=10;
    }
  }
  if (cnt!=11) THROW(critical_error,"misstaged process "+m_name);
  m_emitgluon = m_flavs[m_emit].IsGluon();
  m_name+= "_S"+ToString((int)m_emit)+"_"+ToString((int)m_spect);

  if (m_flavs[m_emit].IsGluon()) {
    p_pl[m_emit]=Pol_Info();
    p_pl[m_emit].Init(2);
    p_pl[m_emit].pol_type = 'e'; 
    p_pl[m_emit].type[0] = 90;
    p_pl[m_emit].type[1] = 91;
    p_pl[m_emit].factor[0] = 1.;
    p_pl[m_emit].factor[1] = 1.;
  }

  Flavour* fl = &m_flavs.front();
  p_hel    = new Helicity(m_nin,m_nout,fl,p_pl);
  p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nin+m_nout,fl,p_b);  
  p_BS->Setk0(s_gauge);
  p_epol   = epol;
  //////////////////////////////////////////////// 

#ifdef Basic_Sfuncs_In_MHV
  p_momlist = p_BS; 
#else
  p_momlist = new MomentumList(m_nin,m_nout); 
#endif 

  int *plist = new int[m_nin+m_nout];
  for (size_t i=0;i<m_nin;i++) plist[i]=fl[i];
  for (size_t i=m_nin;i<m_nin+m_nout;i++) plist[i]=-fl[i];
  p_MHVamp = FullAmplitude_MHV_Handler(model->p_model,&m_cpls,m_nin+m_nout,plist,p_momlist,m_ownamps,m_emit,m_spect); 

  delete [] plist;
  //////////////////////////////////////////////

  p_shand  = new String_Handler(m_gen_str,p_BS,model->p_model->GetCouplings());

  int ntchanmin(m_ntchanmin);
  p_ampl   = new Amplitude_Handler(m_nin+m_nout,&m_flavs.front(),p_b,p_pinfo,model,top,m_maxcpl,m_mincpl,ntchanmin,
                                   &m_cpls,p_BS,p_shand,m_print_graphs,0,true);
  m_ntchanmin=ntchanmin;
  if (p_ampl->GetGraphNumber()==0) {
    msg_Tracking()<<"Single_LOProcess_MHV::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return 0;
  }

  map<string,Complex> cplmap;
  for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
    cplmap.clear();
    if (FlavCompare(links[j]) && p_ampl->CompareAmplitudes(links[j]->GetAmplitudeHandler(),m_sfactor,cplmap)) {
      if (p_hel->Compare(links[j]->GetHelicity(),m_nin+m_nout)) {
	Single_LOProcess_MHV *pp=dynamic_cast<Single_LOProcess_MHV*>(links[j]);
	if (m_emit!=pp->m_emit || m_spect!=pp->m_spect ||
	    p_sub->m_ijt!=pp->p_sub->m_ijt || p_sub->m_kt!=pp->p_sub->m_kt ||
	    p_sub->m_i!=pp->p_sub->m_i || p_sub->m_j!=pp->p_sub->m_j || p_sub->m_k!=pp->p_sub->m_k) continue;
	m_sfactor = sqr(m_sfactor);
	msg_Tracking()<<"AMEGIC::Single_Process_MHV::InitAmplitude : Found compatible process for "<<Name()<<" : "<<links[j]->Name()<<endl;
	  
	p_partner = (Single_LOProcess_MHV*)links[j];
	m_iresult = p_partner->Result()*m_sfactor;
	InitFlavmap(p_partner);

	Minimize();
	return 1;
      }
    }
  }

  p_ampl->FillPointlist();
  p_BS->Initialize();

  int tr=Tests(pfactors);
//   PRINT_INFO("Tests Result: "<<tr);
  switch (tr) {
  case 1 :
    if (p_partner==this) links.push_back(this);
    Minimize();
   
    return 1;
  case -3: return 0;
  default :
    msg_Error()<<"ERROR in Single_LOProcess_MHV::InitAmplitude : "<<std::endl
	       <<"   Failed for "<<m_name<<"."<<endl;
//     errs.push_back(this);
    return 0;
  }
  return 1;
}



int Single_LOProcess_MHV::Tests(std::vector<double> * pfactors) {

  int number      = 1;
  int gauge_test  = 1;
  int string_test = 1;

  /* ---------------------------------------------------
     
     The reference result for momenta moms

     --------------------------------------------------- */

  string testname = string("");
  if (FoundMappingFile(testname,m_pslibname)) {
    if (testname != string("")) {
      gauge_test = string_test = 0;
    }
  }
  
  double M2 = 0.;
  double helvalue;

  if (gauge_test) {
#ifdef Basic_Sfuncs_In_MHV
     p_BS->Setk0(0);
     p_BS->CalcEtaMu(p_testmoms); 
#else
     p_momlist->PutMomenta(p_testmoms);
#endif    

    msg_Tracking()<<"Single_LOProcess_MHV::Tests for "<<m_name<<std::endl
		  <<"   Prepare gauge test and init helicity amplitudes. This may take some time."
		  <<std::endl;
    if (m_emitgluon) p_MHVamp->SetSqMatrix((*pfactors)[1],CalculateHelicityPhase(p_testmoms));
    for (size_t i=0;i<p_hel->MaxHel();i++) { 
      if (p_hel->On(i) && p_hel->GetEPol(i)==90) {
	helvalue = p_MHVamp->MSquare((*p_hel)[i],p_BS)*p_hel->PolarizationFactor(i); 
	M2      +=  helvalue;
      } 
    }
     
    M2  *= p_MHVamp->ParticlesNorm();
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
 for (size_t i=0; i<p_hel->MaxHel(); ++i) M_doub[i]=0.;
 if (m_emitgluon) p_MHVamp->SetSqMatrix((*pfactors)[1],CalculateHelicityPhase(p_testmoms));
 for (size_t i=0; i<p_hel->MaxHel(); ++i) { 
     if (p_hel->On(i) && p_hel->GetEPol(i)==90) {
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
  msg_Tracking()<<"Single_LOProcess_MHV::Tests for "<<m_name<<std::endl
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
  
  m_libname    = testname;

  /* ---------------------------------------------------
     
     Second test : string test

     --------------------------------------------------- */

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



Complex Single_LOProcess_MHV::CalculateHelicityPhase(const ATOOLS::Vec4D * mom)
{
#ifndef Basic_Sfuncs_In_MHV
  msg_Error()<<"Must use Basic_Sfuncs to generate dipole subtraction terms!"<<endl;
  abort();
#endif
  ATOOLS::Vec4D p=mom[(int)GetEmit()];
    ATOOLS::Vec4D q(p[0],-1.*Vec3D(p));
    ATOOLS::Vec4D e1=(*p_epol)[0];
    double mu12=-1./(2.*(e1*p_BS->Getk0()));
    Complex eta13=2.*csqrt(p*p_BS->Getk0())*csqrt(q*p_BS->Getk0());
    
    Complex eip1 = Complex(0.,1.)*(-p_BS->CalcS(q,e1)*conj(p_BS->CalcS(e1,p))+mu12*eta13)/p_BS->CalcS(q,p);

    return eip1;
}


double Single_LOProcess_MHV::operator()(const ATOOLS::Vec4D_Vector &labmom,const ATOOLS::Vec4D *mom,
					std::vector<double> * pfactors,std::vector<ATOOLS::Vec4D>* epol,const int mode)
{
  if (p_partner!=this) {
    return m_lastxs = p_partner->operator()(labmom,mom,pfactors,epol,mode)*m_sfactor;
  }
  p_int->SetMomenta(labmom);
  p_scale->CalculateScale(labmom,m_cmode);

#ifdef Basic_Sfuncs_In_MHV
     p_BS->CalcEtaMu((ATOOLS::Vec4D*)mom); 
#else
     p_momlist->PutMomenta(mom);
#endif  
  double M2(0.);

  if (m_emitgluon) p_MHVamp->SetSqMatrix((*pfactors)[1],CalculateHelicityPhase(mom));

  for (size_t i=0;i<p_hel->MaxHel();i++) {
    if (p_hel->On(i) && p_hel->GetEPol(i)==90) {
      double mh=p_MHVamp->MSquare((*p_hel)[i],p_BS) * p_hel->Multiplicity(i) * p_hel->PolarizationFactor(i);
      
      M2 += mh;
    }
  }

  m_lastxs = M2*p_MHVamp->ParticlesNorm();
  return m_lastxs;
}



void Single_LOProcess_MHV::Calc_AllXS(const ATOOLS::Vec4D_Vector &labmom,
				      const ATOOLS::Vec4D *mom,std::vector<std::vector<double> > &dsij,const int mode) 
{
  p_int->SetMomenta(labmom);
  p_scale->CalculateScale(labmom,m_cmode);
#ifdef Basic_Sfuncs_In_MHV
  p_BS->CalcEtaMu((ATOOLS::Vec4D*)mom); 
#else
  p_momlist->PutMomenta(mom);
#endif  

  dsij[0][0] =0.;
  for (size_t i=0;i<m_partonlist.size();i++) {
    for (size_t k=i+1;k<m_partonlist.size();k++) {
      dsij[k][i] =0.;
    }
  }

  for (size_t i=0;i<p_hel->MaxHel();i++) {
    if (p_hel->On(i)) {
      double fac = p_hel->Multiplicity(i) * p_hel->PolarizationFactor(i) * p_MHVamp->ParticlesNorm();
      p_MHVamp->CalculateAmps((*p_hel)[i],p_BS);
      dsij[0][0] += p_MHVamp->MSquare(0,0)*fac;
      for (size_t i=0;i<m_partonlist.size();i++) {
	for (size_t k=i+1;k<m_partonlist.size();k++) {
	  dsij[i][k] = dsij[k][i] += p_MHVamp->MSquare(m_partonlist[i],m_partonlist[k])*fac;
	}
      }
    }
  }
//   PRINT_INFO(Name());
//   cout<<"Result 00: "<<p_dsij[0][0]<<endl;
//   for (size_t i=0;i<m_partonlist.size();i++) {
//     for (size_t k=i+1;k<m_partonlist.size();k++) {
//       cout<<"Result "<<i<<k<<": "<<p_dsij[i][k]<<endl;
//     }
//   }
//   cout<<"------------------------------------------------"<<endl;
}

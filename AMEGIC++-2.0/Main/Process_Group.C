#include "Process_Group.H"
#include "Single_Process.H"
#include "Single_Process_MHV.H"

#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"

#include "Running_AlphaS.H"
#include "Combined_Selector.H"
#include "MathTools.H"
#include "Shell_Tools.H"


#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace AMEGIC;
using namespace MODEL;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

/*----------------------------------------------------------------------------------
  
  Constructors

  ----------------------------------------------------------------------------------*/

Process_Group::Process_Group() :
  Process_Base(NULL,0,0,NULL,NULL,NULL,0,0,0,0,0,-1.),
  m_resetted(false), m_weventmode(0), m_enable_mhv(false)
{ 
  m_name  = "Empty_Group"; 
  p_pl    = 0;
  m_procs.clear();
  m_efunc="1";
}

Process_Group::Process_Group(Process_Info* pinfo,int _nin,int _nout,Flavour *& _fl,
			     PDF::ISR_Handler * _isr,BEAM::Beam_Spectra_Handler * _beam,Selector_Data * _seldata,
			     int _gen_str,int _orderQCD, int _orderEW,
			     int _kfactorscheme,int _scalescheme,double _scale,
			     Pol_Info * _pl,int _nex,Flavour * _ex_fl,int usepi,double ycut, double error,
			     std::string e_func, bool enable_mhv) :
  Process_Base(pinfo,_nin,_nout,_fl,_isr,_beam,_gen_str,_orderQCD,_orderEW,
	       _scalescheme,_kfactorscheme,_scale,_pl,_nex,_ex_fl,ycut,error),
  m_resetted(false), m_weventmode(0), m_enable_mhv(enable_mhv)
{
  p_selected  = NULL;

  string stan,oli;
  m_efunc=e_func;
  GenerateNames(m_nin,p_flin,p_plin,m_name,stan,oli);

  p_flavours   = new Flavour[m_nvector];
  p_pl   = new Pol_Info[m_nvector];
  p_b    = new int[m_nvector];
  for (size_t i=0;i<m_nin;i++) { 
    p_flavours[i] = p_flin[i]; 
    p_pl[i] = p_plin[i]; 
    p_b[i]  = -1; 
  }
  for (size_t i=m_nin;i<m_nin+m_nout;i++)  { 
    p_flavours[i] = p_flout[i-m_nin]; 
    p_pl[i] = p_plout[i-m_nin]; 
    p_b[i]  = 1; 
  } 
  for (size_t i=m_nin+m_nout;i<m_nvector;i++) { 
    p_flavours[i] = Flavour(kf::pol); 
    p_b[i]  = 1; 
  }

  m_usepi = usepi;

  ConstructProcesses(_seldata);
  GroupProcesses();

  if (_seldata) p_selector = new Combined_Selector(m_nin,m_nout,p_flavours,_seldata,ycut);
  else {
    if (m_nout>2) 
      msg.Out()<<"WARNING in Process_Group "<<m_name<<endl
	       <<"   No selection cuts specified. Init No_Selector !"<<endl;
    p_selector = new No_Selector();
  }
  SetSelector(p_selector);

  if (m_scalescheme==65 && m_updatescales) {
    double dr   = rpa.gen.DeltaR();
    double ycut = rpa.gen.Ycut();
    m_scale[stp::fac] = ycut*sqr(rpa.gen.Ecms());
    ycut=Min(ycut,ycut*sqr(dr));
    m_scale[stp::as] = ycut*sqr(rpa.gen.Ecms());
    SetScales(m_scale[stp::fac],m_scale[stp::as]);
  }
}



Process_Group::~Process_Group()
{
  for(int i=m_procs.size();i>0;i--) {
    if (m_procs[i-1]) delete m_procs[i-1];
  }
}




/*----------------------------------------------------------------------------------
  
  Management of the processes included in the Process_Group

  ----------------------------------------------------------------------------------*/

void Process_Group::ConstructProcesses(ATOOLS::Selector_Data * _seldata) {
  int bsh_pol=p_beamhandler->Polarisation();
  int beam_is_poled[2]={bsh_pol&1,bsh_pol&2};
  // ====

  p_pinfo->Reshuffle();
  p_pinfo->Expand();
  int nsproc=p_pinfo->NProcs();
  int nout = p_pinfo->Nout();
  int  * flindex = new int[m_nin+nout];
  for (size_t i=0;i<m_nin+nout;i++) flindex[i] = 0;
  char * plindex = new char[m_nin+nout];
  for (size_t i=0;i<m_nin+nout;i++) plindex[i] = ' ';
  Flavour  * _fl  = new Flavour[m_nin+nout];
  Flavour  * ofl  = new Flavour[m_nin+nout];
  Pol_Info  * _pl = new Pol_Info[m_nin+nout];
  Pol_Info  * opl = new Pol_Info[m_nin+nout];
  for (size_t i=0;i<m_nin;i++) {
    ofl[i]=p_flin[i];
    opl[i]=p_plin[i];
  }
  p_pinfo->GetFlavList(ofl+m_nin);
  p_pinfo->GetPolList(opl+m_nin);
  Process_Info *pi = NULL;

  string _name,_stan,_oli;
  bool flag = 1;
  bool take = 0;
  int overflow = 0;
  for (;;) {
    if (!flag) break;
    for (size_t i=0;i<m_nin+nout;++i) {
      if (ofl[i].Size() != 1) {
	_fl[i] = ofl[i][flindex[i]]; 
	_pl[i] = Pol_Info(_fl[i]);
	if (_pl[i].pol_type==opl[i].pol_type) {
	  _pl[i] = opl[i];
	}
	else {
	  if (opl[i].GetPol()!=' ' && opl[i].GetPol()!='s') {
	    msg.Out()<<" WARNING in Process_Group::ConstructProcesses : "<<std::endl
		     <<"   wrong polarisation state in Particle.dat."<<endl
		     <<" Polarisation ignored, run will continue."<<endl;
	  }
	}
      }
      else {
	_fl[i] = ofl[i];
	_pl[i] = opl[i];
      }
    }
    if (CF.ValidProcess(m_nin,_fl,nout,_fl+m_nin)) {
      overflow = SetPolarisations(plindex,_pl,beam_is_poled);
      for (size_t i=0;i<m_nin+nout;++i) {
	if (plindex[i]!=' ') _pl[i].SetPol(plindex[i]);
      }
      for (int j=0;j<nsproc;j++) {
	if (nsproc==1) pi = p_pinfo->GetSubProcess(-1);
	else pi = p_pinfo->GetSubProcess(j);
 	pi->ResetSubList(nout,_fl+m_nin,_pl+m_nin);
	pi->Reshuffle();
	GenerateNames(m_nin,_fl,_pl,_name,_stan,_oli,pi);
	take = 1;
	for (size_t k=0;k<m_procs.size();k++) {
	  if (_name == m_procs[k]->Name()) { 
	    take = 0; 
	    delete pi;
	    break; 
	  }
	}
	if (take) {
	  if (m_enable_mhv && CF.PureGluonic(m_nin,_fl,m_nout,&_fl[m_nin]))
	    Add(new Single_Process_MHV(pi,m_nin,m_nout,_fl,p_isrhandler,p_beamhandler,_seldata,m_gen_str,m_orderQCD,m_orderEW,
				       m_kfactorscheme,m_scalescheme,m_scale[stp::as],_pl,m_nex,p_ex_fl,m_usepi,m_ycut,m_maxerror,m_efunc));
	  else
	    Add(new Single_Process(pi,m_nin,m_nout,_fl,p_isrhandler,p_beamhandler,_seldata,m_gen_str,m_orderQCD,m_orderEW,
				   m_kfactorscheme,m_scalescheme,m_scale[stp::as],_pl,m_nex,p_ex_fl,m_usepi,m_ycut,m_maxerror,m_efunc));
	}
      }
    }
    else take=0;

    if (overflow || take==0) {
      for (size_t i=0; i<m_nin+nout; ++i) plindex[i]=' ';
      for (size_t i=m_nin+nout-1;i>=0;--i) {
	if (ofl[i].Size()-1>flindex[i]) {
	  ++flindex[i];
	  break;
	}
	else {
	  if (i==0) flag = 0;
	  flindex[i] = 0;
	}
	if (i==0) break;
      }
    }
  }
  delete [] _fl;
  delete [] _pl;
  delete [] ofl;
  delete [] opl;
  delete [] flindex;
  delete [] plindex;
}

int Process_Group::SetPolarisations(char * plindex, Pol_Info * pl, int * beam_is_poled) 
{
  for (size_t i=m_nin;i<m_nin+m_nout;++i) {
    if (pl[i].DoFNumber()==1) {
      plindex[i]=pl[i].GetPol();
    }
  }

  for (size_t i=0;i<m_nin;++i) {
    if (!beam_is_poled[i])   plindex[i]=pl[i].GetPol();
  }

  for (size_t i=0;i<m_nin;++i) {
    int over=0;
    if ( beam_is_poled[i]) {
      switch (plindex[i]) {
      case ' ': plindex[i]='+'; break;
      case '+': plindex[i]='-'; break;
      case '-': plindex[i]='0'; break;
      case '0': plindex[i]='+';
	over=1;
      } 
      if (pl[i].DoFNumber()==2 && plindex[i]=='0') {
	plindex[i]='+';
	over=1;
      }
    }
    else 
      over=1;

    if (!over) {
      if (beam_is_poled[1] && plindex[1]==' ') plindex[1]='+';
      return 0;
    }
  }
  return 1;
}

void Process_Group::GroupProcesses() {
  // First : Check for identical masses.
  double * massin  = new double[m_nin];
  double * massout = new double[m_nout];
  double   sum_massin  = 0.;
  double   sum_massout = 0.;
  Flavour* flout = new Flavour[m_nout];
  p_pinfo->GetTotalFlavList(flout);

  for (size_t i=0;i<m_nin;i++)  {
    massin[i]   = p_flin[i].Mass();
    sum_massin += massin[i];
  }
  for (size_t i=0;i<m_nout;i++) {
    massout[i]   = p_flout[i].Mass();
    sum_massout += massout[i];
  }

  bool massok = 1;
  for (size_t i=0;i<m_procs.size();i++) {
    for (size_t j=0;j<m_procs[i]->NIn();j++) {
      if (!(ATOOLS::IsEqual(massin[j],(m_procs[i]->Flavours()[j]).Mass()))) {
	msg.Error()<<"Error in Process_Group::GroupProcesses : "<<std::endl
		   <<"   Incoming masses "<<massin[j]<<" vs. "<<(m_procs[i]->Flavours()[j]).Mass()
		   <<" for "<<p_flin[j]<<" "<<m_procs[i]->Flavours()[j]<<endl
		   <<"   Continue run and hope for the best."<<std::endl;
	massok = 0; break;
      }
    }
    if (!massok) break;
    for (size_t j=0;j<m_procs[i]->NOut();j++) {
      if (!(ATOOLS::IsEqual(massout[j],(m_procs[i]->Flavours()[j+m_procs[i]->NIn()]).Mass()))) {
	msg.Error()<<"Error in Process_Group::GroupProcesses : "<<std::endl
		   <<"   Outgoing masses "<<massout[j]<<" vs. "
		   <<(m_procs[i]->Flavours()[j+m_procs[i]->NIn()]).Mass()
		   <<" for "<<p_flout[j]<<" "<<m_procs[i]->Flavours()[j+m_procs[i]->NIn()]<<endl
		   <<"   Continue run and hope for the best."<<std::endl;
	massok = 0; break;
      }
    }
    if (!massok) break;
  }
  if (massok) {
    SetISRThreshold(ATOOLS::Max(sum_massin,sum_massout));
  }
  else {
    msg.Error()<<"ERROR in Process_Group::GroupProcesses : "<<m_name<<endl
	       <<"   Processes do not have equal masses. Abort the run."<<endl;
    abort();
  }
  delete [] massin;
  delete [] massout;

  std::vector<Process_Base *> singleprocs = m_procs;
  while (m_procs.size()>0) m_procs.pop_back();
  
  Process_Base  * sproc;
  Process_Group * group;
  string         help;
  Flavour        flav1,flav2;
  for (size_t i=0;i<singleprocs.size();i++) {
    sproc      = singleprocs[i];
    help       = string("SG_");
    if (sproc->NIn()==2) {
      flav1    = sproc->Flavours()[0]; 
      flav2    = sproc->Flavours()[1]; 
      if ( (flav1.IsVector()) && (flav2.IsVector()) ) {
	if ( (flav1.Bar() == flav1) && (flav2.Bar() == flav2) ) {
	  if (flav1 == flav2) help += string("V_V_->_");
	                 else help += string("V_V'_->_");
	}
	else if (flav1.Bar() == flav1) { 
	  if (flav2.Charge() > 0) help += string("V0_V+_->_"); 
	                     else help += string("V0_V-_->_"); 
	}
	else if (flav2.Bar() == flav2) { 
	  if (flav1.Charge() > 0) help += string("V+_V0_->_"); 
	                     else help += string("V-_V0_->_"); 
	}
	else {
	  if (flav1.Charge() > 0) help += string("V+_"); 
	                     else help += string("V-_"); 
	  if (flav2.Charge() > 0) help += string("V+_->_"); 
	                     else help += string("V-_->_"); 
	}
      }
      else if ( (flav1.IsFermion()) && (flav2.IsVector()) ) {
	if (flav1.IsAnti())          help += string("fb_");
	                        else help += string("f_");
	if (flav2==flav2.Bar())      help += string("V_->_");
	else if (flav2.Charge() > 0) help += string("V+_->_"); 
	                        else help += string("V-_->_"); 
      }
      else if ( (flav1.IsVector()) && (flav2.IsFermion()) ) {
	if (flav1==flav1.Bar())      help += string("V_");
	else if (flav1.Charge() > 0) help += string("V+_"); 
	                        else help += string("V-_"); 
	if (flav2.IsAnti())          help += string("fb_->_");
	                        else help += string("f_->_");
      }
      else if ( (flav1.IsFermion()) && (flav2.IsFermion()) ) {
	if ( (flav1.IsAnti()) && (flav2==flav1) )             help += string ("fb_fb_->_");
	else if ( !(flav1.IsAnti()) && (flav2==flav1) )       help += string ("f_f_->_");
	else if ( !(flav1.IsAnti()) && (flav2==flav1.Bar()) ) help += string ("f_fb_->_");
	else if ( (flav1.IsAnti()) && (flav2.IsAnti()) )      help += string ("fb_fb'_->_");
	else if ( !(flav1.IsAnti()) && (flav2.IsAnti()) )     help += string ("fb_f'_->_");
	else if ( (flav1.IsAnti()) && !(flav2.IsAnti()) )     help += string ("f_fb'_->_");
	else if ( !(flav1.IsAnti()) && !(flav2.IsAnti()) )    help += string ("f_f'_->_");
      }
    }
    else {
      flav1    = sproc->Flavours()[0]; 
      if (flav1.IsVector())        help += string("V_->_");
      else if (flav1.IsFermion())  help += string("f_->_");
      else                         help += string("S_->_");
    }
    int scalars = 0,fermions = 0,vectors = 0;
    for (size_t j=0;j<sproc->NOut();j++) {
      if ((sproc->Flavours()[sproc->NIn()+j]).IsScalar())  scalars++;
      if ((sproc->Flavours()[sproc->NIn()+j]).IsFermion()) fermions++;
      if ((sproc->Flavours()[sproc->NIn()+j]).IsVector())  vectors++;
    }
    help += ToString(scalars)+"S_"+ToString(fermions)+"F_"+ToString(vectors)+"V";
    
    bool found = 0;
    for (size_t j=0;j<m_procs.size();j++) {
      if (m_procs[j]->Name() == help) {
	m_procs[j]->Add(sproc);
	found = 1;
	break;
      }
    }
    if (!found) {
      group = new Process_Group();
      group->SetName(help);
      group->SetAtoms(0);
      group->SetBeam(p_beamhandler);
      group->SetISR(p_isrhandler);
      group->Add(sproc);
      group->SetParent(this);
      m_procs.push_back(group);
    }
  }

  for (size_t i=0;i<m_procs.size();i++) {
    msg_Tracking()<<"Process_Group::GroupProcesses "<<m_procs[i]->Name()<<" : "<<m_procs[i]->Size()<<endl;
    for (size_t j=0;j<m_procs[i]->Size();j++) 
      msg_Tracking()<<"    "<<((*m_procs[i])[j])->Name()<<endl;
    msg_Tracking()<<"--------------------------------------------------"<<endl;
  }
}

void Process_Group::Add(Process_Base * _proc) 
{
  if (m_procs.size()==0) {
    m_nin     = _proc->NIn();
    m_nout    = _proc->NOut();
    m_nvector    = _proc->NVector();
    m_nstrong = _proc->NStrong();
    m_neweak  = _proc->NEWeak();
    if (p_flavours==NULL) {
      p_flavours = new Flavour[m_nin+m_nout];
      for (size_t i=0;i<m_nin+m_nout;i++) p_flavours[i] = (_proc->Flavours())[i];
    }
  }
  else {
    if (_proc->NVector() > m_nvector) m_nvector = _proc->NVector();
  }
  if ( (m_nin != _proc->NIn()) || (m_nout != _proc->NOut())) {
    msg.Error()<<"Error : Cannot add process "<<_proc->Name()
	       <<" to group "<<m_name<<" ! "<<endl
	       <<"   Inconsistent number of external legs."<<endl
	       <<"  Before : ("<<m_nin<<" -> "<<m_nout<<" )"<<endl
	       <<"  Now    : ("<<_proc->NIn()<<" -> "<<_proc->NOut()<<" )"<<endl;
    return;
  }
  _proc->SetParent(this);
  m_procs.push_back(_proc);
}

bool Process_Group::Find(string _name,Process_Base *& _proc) 
{
  if (m_name==_name) {
    _proc = this;
    return 1;
  }
  for (size_t i=0;i<m_procs.size();i++) {
    if (m_procs[i]->Find(_name,_proc)) return 1;
  }
  return 0;
}

void Process_Group::WriteOutXSecs(std::ofstream & _to)
{
//   _to<<m_name<<"  "<<m_totalxs<<"  "<<m_max<<"  "<<m_totalerr<<" "
//      <<m_totalsum<<" "<<m_totalsumsqr<<" "<<m_n<<endl;
  _to.precision(12);
  _to<<m_name<<"  "<<m_totalxs<<"  "<<m_max<<"  "<<m_totalerr<<" "
     <<m_totalsum<<" "<<m_totalsumsqr<<" "<<m_n<<" "
     <<m_ssum<<" "<<m_ssumsqr<<" "<<m_ssigma2<<" "<<m_sn<<" "<<m_wmin<<" "<<m_son<<endl; 
  for (size_t i=0;i<m_procs.size();i++) m_procs[i]->WriteOutXSecs(_to);
}

void Process_Group::WriteOutHistogram(std::string filename)
{
  Integrable_Base::WriteOutHistogram(filename);
  for (size_t i=0;i<m_procs.size();i++) m_procs[i]->WriteOutHistogram(filename);
}

void Process_Group::SetEvents(const double number) 
{
  m_expevents=m_dicedevents=0;
  m_anasum=m_validanasum=0.0;
  for (size_t i=0;i<m_procs.size();++i) {
    m_procs[i]->SetEvents(number/TriggerEfficiency());
    m_expevents+=m_procs[i]->ExpectedEvents();
//     PRINT_INFO(Parent()->Name()<<" "<<Name()<<" "<<m_procs[i]->Name()<<" "
//  	       <<m_procs[i]->ExpectedEvents()<<" "<<m_expevents);
  }
  if (Parent()==this)
    rpa.gen.SetNumberOfEstimatedEvents(m_expevents);
}

bool Process_Group::SelectOne()
{
  if (m_weventmode<0) return SelectOneFromList();
  DeSelect();
  if (m_totalxs==0) p_selected = m_procs[int(ran.Get()*m_procs.size())];
  else {
    double disc;
    if (m_atoms) {
      // select according to total xsecs.
      disc = m_totalxs * ran.Get();

      /*
      double test=0.;
      for (size_t i=0;i<m_procs.size();i++) {
	test+=m_procs[i]->TotalXS();
      }
      cout<<"m_total(group) = "<<m_totalxs<<" vs. "<<test<<" "<<m_totalxs/test-1<<endl;
      */

      for (size_t i=0;i<m_procs.size();i++) {
	disc -= m_procs[i]->TotalXS();
	if (disc<0.) {
	  p_selected = m_procs[i];
	  p_selected->SelectOne();
	  return true;
	}
      }
      if (disc>0.) { 
	msg.Error()<<"ERROR in Process_Group::SelectOne() : "
		   <<"   Total xsec, max = "<<m_totalxs<<", "<<m_max<<endl;
	return false;
      }
    }
    else {
      double m=0.;
      for (size_t i=0;i<m_procs.size();i++) {
	m+= m_procs[i]->Max();
      }
      if (!ATOOLS::IsEqual(m,m_max)) {
	SetMax(0.);
      }
      disc = m_max * ran.Get();

      /*
      double test=0.;
      for (size_t i=0;i<m_procs.size();i++) {
	test+=m_procs[i]->TotalXS();
      }
      cout<<"m_total(group (atom)) = "<<m_totalxs<<" vs. "<<test<<" "<<m_totalxs/test-1<<endl;
      */

      for (size_t i=0;i<m_procs.size();i++) {
	disc -= m_procs[i]->Max();
	if (disc<0.) {
	  p_selected = m_procs[i];
	  p_selected->SelectOne();
	  return true;
	}
      }
      if (disc>0.) { 
	msg.Error()<<"ERROR in Process_Group::SelectOne() : "
		   <<"   Total xsec, max = "<<m_totalxs<<", "<<m_max<<endl;
	return false;
      }
    }
  }
  return true;
}

bool Process_Group::SelectOneFromList()
{
  DeSelect();
  for (size_t i=0;i<m_procs.size();i++) {
    if (m_procs[i]->DicedEvents()<m_procs[i]->ExpectedEvents()) {
      p_selected=m_procs[i];
      if (p_selected->SelectOneFromList()) {
	if (Parent()==this) {
	  AddEvent(0.0,0.0,1);
	}
	break;
      }
    }
  }
  if (p_selected==NULL) return false;
  return true;
}

void Process_Group::DeSelect() {
  p_selected = 0;
  for (size_t i=0;i<m_procs.size();i++) m_procs[i]->DeSelect();
}

bool Process_Group::ReSelect(int i) {
  if (m_weventmode>=0) return SelectOne();
  if (i==1) return true;
  return Parent()->SelectOne();
}

void Process_Group::Empty() {
  for (size_t i=0;i<m_procs.size();i++) m_procs[i]->Empty();
}



void Process_Group::SetResDir(std::string _resdir) {
  m_resdir = _resdir;
  for (size_t i=0;i<m_procs.size();i++) m_procs[i]->SetResDir(m_resdir);
}


void Process_Group::SetScale(const double _scale)
{
  Process_Base::SetScale(_scale);
  for (size_t i=0;i<m_procs.size();i++) m_procs[i]->SetScale(_scale); 
} 

void Process_Group::SetScales(double q2_fac, double q2_ren)
{ 
  Process_Base::SetScales(q2_fac,q2_ren);
  for (size_t i=0;i<m_procs.size();i++) m_procs[i]->SetScales(q2_fac,q2_ren); 
} 


void Process_Group::SetISRThreshold(double _isrth)
{
  m_threshold = _isrth;
  for (size_t i=0;i<m_procs.size();i++) m_procs[i]->SetISRThreshold(m_threshold); 
} 


void Process_Group::SetTables(bool _tables)
{
  m_tables = _tables;
  for (size_t i=0;i<m_procs.size();i++) m_procs[i]->SetTables(m_tables);
} 

void Process_Group::SetTotal(int flag, int depth)  { 
  if (flag!=2) {
    m_totalxs  = TotalResult(); 
    m_totalerr = TotalVar();
    //     m_totalerr = sqrt( (m_totalsumsqr/m_n - 
    // 			(ATOOLS::sqr(m_totalsum)-m_totalsumsqr)/(m_n*(m_n-1.)) )  / m_n); 
    if ((m_nin==1 && m_nout==2) || m_n==1) m_totalerr = 0.;
    if (p_selector) p_selector->Output();
    m_max = 0.;
    for (size_t i=0;i<m_procs.size();i++) {
      m_procs[i]->SetTotal(flag, depth+1);
      m_max += m_procs[i]->Max(); // naive sum, probably unneccessary large
    }
  }
  else {
    //   flag==2  means  check xs with sum of subprocesses
    //               update maximum to sum of maximum
    SetMax(0., depth+1);
  }

  //  RescaleXSec(1.);
  
  if (m_nin==2 && flag==0) {
    if ( (depth<=0 && msg.LevelIsInfo()) || msg.LevelIsTracking()) {
      for (int i=0;i<depth;++i) msg.Out()<<"  ";
      msg_Info()<<om::bold<<m_name<<om::reset<<" : "
		<<om::blue<<om::bold<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"<<om::reset
		<<" +/- "<<om::reset<<om::blue<<m_totalerr/m_totalxs*100.<<" %,"<<om::reset
		<<om::bold<<" exp. eff: "<<om::red<<(100.*m_totalxs/m_max)<<" %."<<om::reset<<endl;    
    }
  }
  if (m_nin==1) {
    msg_Info()<<"Total width for "<<m_name<<" : "
	      <<m_totalxs<<" GeV"
	      <<" +/- "<<m_totalerr/m_totalxs*100.<<"%, max : "<<m_max<<endl;
  }
}

void Process_Group::SetMax(const double max, int depth) {
  if (max>0.) {
    m_max=max;
    return;
  }
  // parameter is dummy!
  double sum = 0.;
  m_max = 0.;
  for (size_t i=0;i<m_procs.size();i++) {
    m_procs[i]->SetTotal(2,depth);
    sum   += m_procs[i]->TotalXS();
    m_max += m_procs[i]->Max(); // naive sum, probably unneccessary large
  }
  if (m_totalxs!=0.) {
    if (!ATOOLS::IsEqual(sum,m_totalxs)) {
      int io = msg.Out().precision(12);
      /*
      msg.Out()<<"WARNING in Process_Group::SetMax :"<<std::endl
	       <<"   In group "<<Name()<<": xs and sum of daughters does not agree ! "<<endl
	       <<" sum="<<sum<<"  total:"<<m_totalxs
	       <<"  ("<<((sum-m_totalxs)/m_totalxs)<<")"<<endl;
      */
      msg.Out().precision(io);
    }
    if (m_atoms) m_totalxs=sum;
  }
}

void Process_Group::ResetMax(int flag) {
  m_max = 0.;
  for (size_t i=0;i<m_procs.size();++i) {
    m_procs[i]->ResetMax(flag);
  }
}
void Process_Group::OptimizeResult() 
{
  double ssigma2 = ATOOLS::sqr(m_ssum/m_sn)/((m_ssumsqr/m_sn - ATOOLS::sqr(m_ssum/m_sn))/(m_sn-1));
  if (ssigma2>m_wmin) {
    m_ssigma2  += ssigma2; 
    m_totalsum += m_ssum*ssigma2/m_sn;
    m_totalsumsqr+= m_ssumsqr*ssigma2/m_sn;
    m_ssum     = 0.;
    m_ssumsqr  = 0.;
    m_sn       = 0;
    if (ssigma2/m_son>m_wmin) m_wmin = ssigma2/m_son;
    m_son      = 0;
  }
//   cout<<"Weights (actual/min) "<<ssigma2<<" "<<m_wmin<<endl;
  m_son++;
  for (size_t i=0;i<m_procs.size();++i) m_procs[i]->OptimizeResult();
}

void Process_Group::InitWeightHistogram() {
  Integrable_Base::InitWeightHistogram();
  for (size_t i=0;i<m_procs.size();++i) m_procs[i]->InitWeightHistogram();
}

void Process_Group::SetMaxJetNumber(int max) {
  for (size_t i=0;i<m_procs.size();i++) {
    m_procs[i]->SetMaxJetNumber(max);
  }  
  m_maxjetnumber = max;
}

void Process_Group::SetAtoms(bool _atoms) { m_atoms = _atoms; }


/*----------------------------------------------------------------------------------
  
  Initialization of the processes

  ----------------------------------------------------------------------------------*/

int Process_Group::InitAmplitude(Interaction_Model_Base * model,Topology * top,Vec4D *& testmoms,
				 vector<Process_Base *> & links,vector<Process_Base *> & errs,
				 int & totalsize, int & procs, int & current_atom)
{
  int okay = 1;
  vector <string> deletethem;

  for (size_t i=0;i<m_procs.size();i++) {
    if (m_atoms) { 
      delete [] testmoms; testmoms = 0; 
      current_atom = links.size();
    }
    switch (m_procs[i]->InitAmplitude(model,top,testmoms,links,errs,totalsize,procs,current_atom)) {
    case -3 :
      msg_Tracking()<<"Amplitude is zero for "<<m_procs[i]->Name()<<endl
		     <<"   delete it."<<endl;
      deletethem.push_back(m_procs[i]->Name());
      break;
    case -2 : 
      msg.Error()<<"Error in creation of amplitude "<<m_procs[i]->Name()<<endl;
      return -2;
    case -1 : 
      msg_Tracking()<<"No diagrams or amplitudes for "<<m_procs[i]->Name()<<endl
		    <<"   delete it."<<endl;
      deletethem.push_back(m_procs[i]->Name());
      break;
    case 0 :
      okay = 0;      
      break;
    default :
      break;
    }
  }

  bool flag = 1;
  if (errs.size()>0){
    double sum = 0;
    for (size_t i=0;i<m_procs.size();i++) sum+=m_procs[i]->Result();

    for (size_t i=0;i<errs.size();i++) 
      if (!ATOOLS::IsZero(errs[i]->Result()/sum)) flag = 0;
    if (!flag) return -2;
    
    //delete
    for (size_t i=0;i<errs.size();i++) 
      for (size_t j=0;j<m_procs.size();j++) 
	if (errs[i]->Name() == m_procs[j]->Name()) {
	  deletethem.push_back(m_procs[j]->Name());
	  msg.Out()<<"Faulty process "<<m_procs[j]->Name()<<" is negligible"<<endl
		   <<"   delete it."<<endl;
	}
    errs.clear();
  }

  for (size_t i=0;i<m_procs.size();i++) {
    flag = 0;
    if (m_procs[i]->Size() == 0) flag = 1;
    if (!flag) {
      for (size_t j=0;j<deletethem.size();j++) {
	if (m_procs[i]->Name() == deletethem[j]) flag = 1;
      }
    }
    if (flag) {
      delete m_procs[i];
      for (size_t j=i;j<m_procs.size()-1;j++) m_procs[j] = m_procs[j+1];
      //just to continue at position i in the next try
      i--;
      m_procs.pop_back();
    }
  }

  if (okay==0) {
    links.clear();
    for (size_t i=0;i<m_procs.size();i++) if (m_procs[i]) delete (m_procs[i]); 
    m_procs.clear();
  }

  return okay;
}



bool Process_Group::SetUpIntegrator()
{
  bool okay = 1;
  if (m_atoms) {
    PRINT_INFO(m_procs.size());
    for (size_t i=0;i<m_procs.size();i++) {
      if (m_procs[i]->Partner()==NULL) {
	if (!(m_procs[i]->SetUpIntegrator())) okay = 0;
      }
    }
    return okay;
  }
  
  if (m_nin==2) {
    if ( (p_flavours[0].Mass() != p_isrhandler->Flav(0).Mass()) ||
	 (p_flavours[1].Mass() != p_isrhandler->Flav(1).Mass()) ) p_isrhandler->SetPartonMasses(p_flavours);
  }
  p_pshandler  = new Phase_Space_Handler(this,p_isrhandler,p_beamhandler,m_maxerror);
  p_pshandler->SetUsePI(m_usepi);

  //  if (m_nin==2 ) 
  AddChannels(this);
//  if (m_nin==2) 
  { for (size_t i=0;i<m_procs.size();i++) m_procs[i]->Empty(); }
  return 1;
}




/*----------------------------------------------------------------------------------
  
  Evaluation of the processes included in All_Processes

  ----------------------------------------------------------------------------------*/

bool Process_Group::CalculateTotalXSec(std::string _resdir)
{
  msg_Info()<<"Process_Group::CalculateTotalXSec("<<_resdir<<")"<<endl;
  if (m_atoms) {
    bool okay = 1;
    for (size_t i=0;i<m_procs.size();i++) {
      if (!(m_procs[i]->CalculateTotalXSec(_resdir))) okay = 0;
    }
    return okay;
  }
  else {
    std::string filename =_resdir+"/"+m_name+".xs_tot";
    std::string histofile =_resdir+string("/WD_")+m_name;
    string _name;
    double _totalxs,_totalerr,_max,sum,sqrsum,ssum,ssqrsum,ss2,wmin;
    long int n,sn,son;
    if (_resdir!=string("")) {
      if (IsFile(filename)) {
	ifstream from;
	bool okay=1;
	from.open(filename.c_str());
	while (from) {
	  from>>_name>>_totalxs>>_max>>_totalerr>>sum>>sqrsum>>n>>ssum>>ssqrsum>>ss2>>sn
	      >>wmin>>son;
	  Process_Base * proc = NULL;
	  if (Find(_name,proc)) {
	    proc->SetTotalXS(_totalxs);
	    proc->SetTotalError(_totalerr);
	    proc->SetMax(_max);
	    proc->SetSum(sum);
	    proc->SetSumSqr(sqrsum);
	    proc->SetPoints(n);
	    proc->SetSSum(ssum);
	    proc->SetSSumSqr(ssqrsum);
	    proc->SetSigmaSum(ss2);
	    proc->SetSPoints(sn);
	    proc->SetWMin(wmin);
	    proc->SetOptCounter(son);
	    proc->ReadInHistogram(histofile);
	  }
	  else {
	    okay = 0;
	  }
	}
	from.close();
	p_pshandler->ReadIn(_resdir+string("/MC_")+m_name);
	if (p_pshandler->BeamIntegrator() != 0) p_pshandler->BeamIntegrator()->Print();
	if (p_pshandler->ISRIntegrator() != 0)  p_pshandler->ISRIntegrator()->Print();
	if (p_pshandler->FSRIntegrator() != 0)  p_pshandler->FSRIntegrator()->Print();
	p_pshandler->InitIncoming();
	ReadInHistogram(histofile);
	if (m_totalxs<=0.) {
	  msg.Error()<<"ERROR in Process_Group::CalculateTotalXSec :"
		     <<"   In "<<m_name<<"::CalculateTotalXSec("<<_resdir<<")"<<endl
		     <<"   Something went wrong : Negative xsec : "<<m_totalxs<<endl;
	  return 0;
	}
	else {
	  if (okay) {
	    msg_Tracking()<<"In "<<m_name<<"::CalculateTotalXSec("<<_resdir<<")"<<endl
			  <<"   Found all xsecs. Continue"<<endl;
	    SetTotal(2);
	  }
	}
      }
    }
    if (m_nin==2) {
      if (p_flavours[0].Mass()!=p_isrhandler->Flav(0).Mass() ||
	  p_flavours[1].Mass() != p_isrhandler->Flav(1).Mass()) 
	p_isrhandler->SetPartonMasses(p_flavours);
    }
    m_tables = 0;
    m_resultpath=_resdir;
    m_resultfile=filename;
    m_histofile=histofile;
    ATOOLS::Exception_Handler::AddTerminatorObject(this);
    double var=TotalVar();
    m_totalxs = p_pshandler->Integrate();
    
    if (m_nin==2) m_totalxs /= ATOOLS::rpa.Picobarn(); 
    
    SetTotal(0);
        
    if (!(ATOOLS::IsZero((m_totalxs-TotalResult())/(m_totalxs+TotalResult())))) {
      msg.Error()<<"ERROR in Process_Group::CalculateTotalXSec :"
		 <<"Result of PS-Integrator and internal summation do not coincide!"<<endl
		 <<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<TotalResult()<<endl;
    }

    if (m_totalxs>0.) {
      if (ATOOLS::IsEqual(var,TotalVar())) {
	ATOOLS::Exception_Handler::RemoveTerminatorObject(this);
	return 1;
      }
      if (_resdir!=string("")) {
	std::ofstream to;
	to.open(filename.c_str(),ios::out);
	to.precision(12);
	msg_Info()<<"Store result : xs for "<<m_name<<" : ";
	if (m_nin==2) msg_Info()<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb";
	if (m_nin==1) msg_Info()<<m_totalxs<<" GeV";
	msg_Info()<<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
		  <<"       max : "<<m_max<<endl;
	WriteOutXSecs(to);
	int  mode_dir = 493;
	ATOOLS::MakeDir(histofile.c_str(),mode_dir,0); 
	WriteOutHistogram(histofile);
	p_pshandler->WriteOut(_resdir+string("/MC_")+m_name);
	to.close();
      }
      ATOOLS::Exception_Handler::RemoveTerminatorObject(this);
      return 1;
    }
    ATOOLS::Exception_Handler::RemoveTerminatorObject(this);
  }
  return 0;
}

void Process_Group::PrepareTerminate()
{
  if (m_resultpath.length()==0 && m_resultfile.length()==0) return;
  SetTotal(0);
  if (m_totalxs<=0.) return;
  std::ofstream to;
  to.open(m_resultfile.c_str(),ios::out);
  to.precision(12);
  msg_Info()<<"Store result : xs for "<<m_name<<" : ";
  if (m_nin==2) msg_Info()<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb";
  if (m_nin==1) msg_Info()<<m_totalxs<<" GeV";
  msg_Info()<<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
	    <<"       max : "<<m_max<<endl;
  WriteOutXSecs(to);
  int  mode_dir = 493;
  ATOOLS::MakeDir(m_histofile.c_str(),mode_dir,0); 
  WriteOutHistogram(m_histofile);
  p_pshandler->WriteOut(m_resultpath+string("/MC_")+m_name);
  to.close();
}

void  Process_Group::RescaleXSec(double fac) {
  double sumxs=0.;
  for (size_t i=0;i<m_procs.size();++i) sumxs +=m_procs[i]->TotalXS();
  if (sumxs!=0.) m_totalxs = sumxs;
  
  for (size_t i=0;i<m_procs.size();i++) {
    m_procs[i]->RescaleXSec(fac);
  }
  Process_Base::RescaleXSec(fac);
}

void Process_Group::SetupEnhance() {
  if (m_enhancefac==1. && m_maxfac==1. && m_maxeps==0.) return;

  double xs=TotalXS();
  for (size_t i=0;i<m_procs.size();++i) {
    m_procs[i]->SetEnhance(m_enhancefac,m_maxfac,m_maxeps);
    m_procs[i]->SetupEnhance();
  }
  if (m_enhancefac!=1.) {
    SetTotalXS(xs*m_enhancefac);
  }
}


bool Process_Group::LookUpXSec(double ycut,bool calc,string obs) {
  bool okay = 1;
  if (m_atoms) {
    for (size_t i=0;i<m_procs.size();i++) {
      if (!(m_procs[i]->LookUpXSec(ycut,1,obs))) okay = 0;
    }
    if (msg.LevelIsTracking() && okay) {
      msg.Out()<<"Process_Group::LookUpXSec() : "<<std::endl
	       <<"   Read in cross section for "<<m_name<<" from file in directory "<<m_resdir<<endl;
    }
    return okay;
  }
  else {
    for (size_t i=0;i<m_procs.size();i++) {
      if (!(m_procs[i]->LookUpXSec(ycut,0,obs))) okay = 0;
    }
    if (okay) {                    
      m_totalxs = 0; m_max = 0;
      for (size_t i=0;i<m_procs.size();i++) {
	m_totalxs += m_procs[i]->TotalXS();
	m_max     += m_procs[i]->Max();
      }
      p_pshandler->ReadIn(m_resdir+string("/MC_")+m_name);
      if (msg.LevelIsTracking()) {
	msg.Out()<<"Process_Group::LookUpXSec() : "<<std::endl
		 <<"   Read in cross sections for "<<m_name<<" from file in directory "<<m_resdir<<endl;
	if (p_pshandler->BeamIntegrator() != 0) p_pshandler->BeamIntegrator()->Print();
	if (p_pshandler->ISRIntegrator() != 0)  p_pshandler->ISRIntegrator()->Print();
	if (p_pshandler->FSRIntegrator() != 0)  p_pshandler->FSRIntegrator()->Print();
      }
      return 1;
    }
    if (calc) {
      if (!(PrepareXSecTables())) return 0;
      okay = 1;
      for (size_t i=0;i<m_procs.size();i++) {
	if (!(m_procs[i]->LookUpXSec(ycut,0,obs))) okay = 0;
      }
      if (okay) {                    
	m_totalxs = 0; m_max = 0;
	for (size_t i=0;i<m_procs.size();i++) {
	  m_totalxs += m_procs[i]->TotalXS();
	  m_max     += m_procs[i]->Max();
	}
	msg_Tracking()<<"Process_Group::LookUpXSec() : "<<std::endl
		      <<"   Read in cross sections for "<<m_name<<" from file in directory "<<m_resdir<<endl;
	return 1;
      }
    }
    return 0;
  }
}

bool Process_Group::PrepareXSecTables()
{
  msg_Info()<<"Process_Group::PrepareXSecTables()"<<std::endl;
  if (m_atoms) {
    bool okay = 1;
    for (size_t i=0;i<m_procs.size();i++) {
      if (!(m_procs[i]->PrepareXSecTables())) okay = 0;
    }
    return okay;
  }
  else {
    if (m_nin==2) {
      if ( (p_flavours[0].Mass() != p_isrhandler->Flav(0).Mass()) ||
	   (p_flavours[1].Mass() != p_isrhandler->Flav(1).Mass()) ) p_isrhandler->SetPartonMasses(p_flavours);
    }
    m_totalxs = p_pshandler->Integrate()/ATOOLS::rpa.Picobarn(); 
    if (!(ATOOLS::IsZero((m_totalxs-TotalResult())/(m_totalxs+TotalResult())))) {
      msg.Error()<<"ERROR in Process_Group::PrepareXSecTables :"<<std::endl
		 <<"   Result of PS-Integrator and internal summation do not coincide for"<<endl
		 <<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<TotalResult()<<endl;
    }
    SetTotal(1);
    p_pshandler->WriteOut(m_resdir+string("/MC_")+m_name);
    if (m_totalxs>0.) return 1;
  }
  return 0;
}






void Process_Group::AddPoint(const double value) 
{
  Integrable_Base::AddPoint(value);

  for (size_t i=0;i<m_procs.size();i++) {
    if (dabs(m_last)>0.) {
      m_procs[i]->AddPoint(value*m_procs[i]->Last()/m_last);
    }
    else {
      m_procs[i]->AddPoint(0.);
    }  
  }
}






double Process_Group::Differential(const Vec4D * p)
{
  m_last = 0;
  for (size_t i=0;i<m_procs.size();i++) {
    m_last += m_procs[i]->DSigma(p,1);
  }
  //ControlOutput(p);

  if ((!(m_last<=0)) && (!(m_last>0))) {
    msg.Error()<<"ERROR in Process_Group::Differential :"<<endl;
    for (size_t i=0;i<m_nin+m_nout;i++) 
      msg.Error()<<"   "<<i<<" th Momentum "<<p[i]<<endl;
    PrintDifferential();
  }
  return m_last;
}



double Process_Group::Differential2()
{
  if (p_isrhandler->On()==0) return 0.;
  double tmp = 0.;
  for (size_t i=0;i<m_procs.size();i++) tmp += m_procs[i]->DSigma2();

  if ((!(tmp<=0)) && (!(tmp>0))) {
    msg.Error()<<"ERROR in Process_Group::Differential2 :"<<endl;
    PrintDifferential();
  }
  m_last += tmp;
  return tmp;
}



double Process_Group::DSigma(const Vec4D * p,bool lookup)
{
  m_last = 0;
  for (size_t i=0;i<m_procs.size();i++) m_last += m_procs[i]->DSigma(p,lookup);
  //ControlOutput(p);
  return m_last;
}



double Process_Group::DSigma2()
{
  double tmp = 0.;
  for (size_t i=0;i<m_procs.size();i++) tmp += m_procs[i]->DSigma2();
  m_last += tmp;
  return tmp;
}


bool Process_Group::OneEvent(double _mass) {
  if (m_atoms) {
    SelectOne();
    return dynamic_cast<Process_Base*>(p_selected)->OneEvent(_mass);
  }
  return p_pshandler->OneEvent(_mass);
}

bool Process_Group::SameEvent() {
  if (m_atoms) {
    if (p_selected)
      return p_selected->SameEvent();
    msg.Error()<<" ERROR in bool Process_Group::SameEvent() "<<endl;
    return 0;
  }
  return p_pshandler->SameEvent();
}

ATOOLS::Blob_Data_Base *Process_Group::WeightedEventNS(const int mode) 
{
  if (m_atoms) {
    bool oldresetted=false;
    while (!oldresetted) {
      while (SelectOneFromList()) {
	ATOOLS::Blob_Data_Base *info = p_selected->WeightedEventNS(mode);
	if (rpa.gen.NumberOfDicedEvents()==rpa.gen.NumberOfEvents()-1)
	  rpa.gen.SetNumberOfEvents(rpa.gen.NumberOfEvents()+1);
	if (info!=NULL) return info;
      }
      oldresetted=m_resetted;
      if (!m_resetted) {
	double gmin=10.;
	ResetEvents(gmin);
	m_resetted=true;
      }
    }
    rpa.gen.SetNumberOfEvents(rpa.gen.NumberOfDicedEvents());
    return NULL;
  }
  return p_pshandler->WeightedEvent(mode);
}

ATOOLS::Blob_Data_Base *  Process_Group::WeightedEvent(const int mode) {
  if (m_weventmode<0) return WeightedEventNS(mode);
  if (m_atoms) {
    SelectOne();
    return p_selected->WeightedEvent(mode);
  }
  return p_pshandler->WeightedEvent(mode);
}

ATOOLS::Blob_Data_Base *  Process_Group::SameWeightedEvent() {
  if (m_atoms) {
    return p_selected->SameWeightedEvent();
  }
  return p_pshandler->SameWeightedEvent();
}




/*----------------------------------------------------------------------------------
  
  Helpers

  ----------------------------------------------------------------------------------*/

void Process_Group::PrintDifferential()
{
  if (!(msg.LevelIsDebugging())) return;
  m_last = 0;
  msg.Out()<<"--------------------------------------------------------"<<endl
	   <<"--------------------------------------------------------"<<endl
	   <<"--------------------------------------------------------"<<endl;
  for (size_t i=0;i<m_procs.size();i++) {
    msg.Out()<<"====================================================="<<endl;
    m_procs[i]->PrintDifferential();
    m_last += m_procs[i]->Last();
    msg.Out()<<"====================================================="<<endl;
  }
  msg.Out()<<"--------------------------------------------------------"<<endl
	   <<" Total : "<<m_last<<endl
	   <<"--------------------------------------------------------"<<endl
	   <<"--------------------------------------------------------"<<endl
	   <<"--------------------------------------------------------"<<endl;
}

void Process_Group::ControlOutput(Vec4D * p)
{ 
  msg.Out()<<"***************************************************************************"<<endl
	   <<"***************************************************************************"<<endl;
  double s   = (p[0]+p[1]).Abs2();
  double t   = (p[0]-p[2]).Abs2();
  double u   = (p[0]-p[3]).Abs2();
  m_scale[stp::as]    = ATOOLS::sqr(ATOOLS::rpa.gen.Ecms());
  double a_s = as->AlphaS(m_scale[stp::as]);
  msg.Out()<<"-------- Process_Group : "<<m_name<<" : DSigma -------------------"<<endl
	   <<"         scale = "<<m_scale[stp::as]<<" = "<<2.*s*t*u/(s*s+u*u+t*t)<<endl
	   <<"         s,t,u = "<<s<<", "<<t<<", "<<u<<" : "<<sqrt(4.*M_PI*a_s)<<endl
	   <<"-----------------------------------------------------------------------"<<endl;
  double g4  = sqr(4.*M_PI*a_s);
  if (m_name == string("gg -> gg"))
    msg.Out()<<"gg   -> gg   : "
	     <<g4 * ((3.-(t*u)/(s*s)-(s*u)/(t*t)-(t*s)/(u*u)))<<endl;
  if (m_name == string("qg -> qg"))
    msg.Out()<<"qg   -> qg   : "
	     <<g4 * (-4./9.*(s*s+u*u)/(s*u)+(s*s+u*u)/(t*t))<<endl;
  if (m_name == string("qbg -> qbg"))
    msg.Out()<<"qg   -> qg   : "
	     <<g4 * (-4./9.*(s*s+u*u)/(s*u)+(s*s+u*u)/(t*t))<<endl;
  if (m_name == string("qbq'b -> qbq'b"))
    msg.Out()<<"qq'  -> qq'   : "
	     <<g4 * (4.*(s*s+u*u)/(9.*t*t))<<endl;

  // if (m_name == string("SG_f f' -> 0S 2F 0V"))
    msg.Out()<<"qq'  -> qq'   : "
	     <<g4 * (4.*(s*s+u*u)/(9.*t*t))<<endl;
  if (m_name == string("qq'b -> qq'b"))
    msg.Out()<<"qq'b -> qq'b  : "
	     <<g4 * (4.*(s*s+u*u)/(9.*t*t))<<endl;
  // if (m_name == string("qq -> qq"))
    msg.Out()<<"qq   -> qq    : "
	     <<g4 * ((4./9.*((s*s+u*u)/(t*t)+(s*s+t*t)/(u*u))-
		       8./27.*(s*s)/(u*t))/2.)<<endl;
  if (m_name == string("SG_f f -> 0S 2F 0V"))
    msg.Out()<<"qq   -> qq    : "
	     <<g4 * ((4./9.*((s*s+u*u)/(t*t)+(s*s+t*t)/(u*u))-
		       8./27.*(s*s)/(u*t))/2.)<<endl;
  if (m_name == string("qqb -> q'q'b"))
    msg.Out()<<"qqb  -> q'q'b : "
	     <<g4 * (4.*(t*t+u*u)/(9.*s*s))<<endl;
  if (m_name == string("qqb -> qqb"))
    msg.Out()<<"qqb  -> qqb   : "
	     <<g4 * (4./9.*((s*s+u*u)/(t*t)+(u*u+t*t)/(s*s))-
		      8./27.*(u*u)/(s*t))<<endl;
  if (m_name == string("qqb -> gg"))
    msg.Out()<<"qqb  -> gg    : "
	     <<g4 * ((32./27.*(t*t+u*u)/(t*u)-
		       8./3.*(t*t+u*u)/(s*s))/2.)<<endl;
  if (m_name == string("gg -> qqb"))
    msg.Out()<<"gg   -> qqb   : "
	     <<g4 * (1./6.*(t*t+u*u)/(t*u)-3./8.*(t*t+u*u)/(s*s))<<endl;
  msg.Out()<<"-----------------------------------------------------------------------"<<endl;
}

void Process_Group::SetPrintGraphs(bool print_graphs) 
{
 m_print_graphs=print_graphs; 
 std::cout<<" "<<Name()<<" setprintgraphs : "<<print_graphs<<std::endl;
 for (size_t i=0;i<m_procs.size();i++) m_procs[i]->SetPrintGraphs(print_graphs);
}

void Process_Group::AddEvent(const double xs,const double validxs,const int ncounts)
{ 
  m_dicedevents+=ncounts;
  if (xs!=0.) m_accevents++;
  p_selected->AddEvent(xs,validxs,ncounts);
}

void Process_Group::ResetEvents(double gmin)
{
  m_gmin = gmin;
//   PRINT_INFO(Name());
  if (Parent()==this){
    double number = 0.;
    GetGMin(gmin,number);
    number/=(double)m_dicedevents*rpa.gen.NumberOfEvents()/(double)m_accevents;
//     PRINT_INFO(gmin<<" "<<number);
    gmin=ATOOLS::Min(gmin,number);
    rpa.gen.SetWAnaScale(gmin);
  }
  m_expevents=0;
  for (size_t i=0;i<m_procs.size();++i) {
    m_procs[i]->ResetEvents(gmin);
    m_expevents+=m_procs[i]->ExpectedEvents();
//     PRINT_INFO(m_dicedevents<<" "<<m_expevents<<" "<<gmin);    
  }
  if (Parent()==this)
    rpa.gen.SetNumberOfEstimatedEvents(m_expevents);
}

void Process_Group::GetGMin(double &g, double &meff)
{
  for (size_t i=0;i<m_procs.size();++i) m_procs[i]->GetGMin(g,meff);
}

void Process_Group::SetWEventMode(int mode) 
{ 
  for (size_t i=0;i<m_procs.size();++i) m_procs[i]->SetWEventMode(mode);
  m_weventmode=mode; 
}

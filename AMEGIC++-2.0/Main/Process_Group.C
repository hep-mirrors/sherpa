#include "Process_Group.H"
#include "Single_Process.H"

#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"

#include "Running_AlphaS.H"
#include "Combined_Selector.H"
#include "Primitive_Observable_Base.H"
#include "MathTools.H"


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
  Process_Base(0,0,NULL,NULL,NULL,0,0,0,0,0,1.,-1.)
{ 
  m_name  = "Empty_Group"; 
  p_fl    = 0;
  p_pl    = 0;

  p_sel   = 0; 
  p_cuts  = 0;
  p_ps    = 0;
  p_moms  = 0;
  m_procs.clear();
}



Process_Group::Process_Group(int _nin,int _nout,Flavour *& _fl,
			     PDF::ISR_Handler * _isr,BEAM::Beam_Spectra_Handler * _beam,Selector_Data * _seldata,
			     int _gen_str,int _orderQCD, int _orderEW,
			     int _kfactorscheme,int _scalescheme,double _scalefactor,double _scale,
			     Pol_Info * _pl,int _nex,Flavour * _ex_fl) :
  Process_Base(_nin,_nout,_fl,_isr,_beam,_gen_str,_orderQCD,_orderEW,
	       _scalescheme,_kfactorscheme,_scalefactor,_scale,_pl,_nex,_ex_fl)
{
  p_selected  = NULL;

  string stan,oli;
  GenerateNames(m_nin,p_flin,p_plin,m_nout,p_flout,p_plout,m_name,stan,oli);

  p_fl   = new Flavour[m_nvec];
  p_pl   = new Pol_Info[m_nvec];
  p_b    = new int[m_nvec];
  for (short int i=0;i<m_nin;i++) { 
    p_fl[i] = p_flin[i]; 
    p_pl[i] = p_plin[i]; 
    p_b[i]  = -1; 
  }
  for (short int i=m_nin;i<m_nin+m_nout;i++)  { 
    p_fl[i] = p_flout[i-m_nin]; 
    p_pl[i] = p_plout[i-m_nin]; 
    p_b[i]  = 1; 
  } 
  for (short int i=m_nin+m_nout;i<m_nvec;i++) { 
    p_fl[i] = Flavour(kf::pol); 
    p_b[i]  = 1; 
  }

  ConstructProcesses(_seldata);
  GroupProcesses();

  InitCuts();
  if (_seldata) p_sel = new Combined_Selector(m_nin,m_nout,p_fl,_seldata);
  else {
    if (m_nout>2) 
      msg.Error()<<"Potential Error in Process_Group "<<m_name<<endl
		 <<"   No selection cuts specified. Init No_Selector !"<<endl;
    p_sel = new No_Selector();
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
  int bsh_pol=p_beam->Polarisation();
  int beam_is_poled[2]={bsh_pol&1,bsh_pol&2};
  // ====

  int  * flindex = new int[m_nin+m_nout];
  for (int i=0;i<m_nin+m_nout;i++) flindex[i] = 0;
  char * plindex = new char[m_nin+m_nout];
  for (int i=0;i<m_nin+m_nout;i++) plindex[i] = ' ';
  Flavour  * _fl  = new Flavour[m_nin+m_nout];
  Pol_Info  * _pl = new Pol_Info[m_nin+m_nout];

  string _name,_stan,_oli;
  bool flag = 1;
  bool take;
  int overflow;
  for (;;) {
    if (!flag) break;
    for (int i=0;i<m_nin+m_nout;++i) {
      if (p_fl[i].Size() != 1) {
	_fl[i] = p_fl[i][flindex[i]]; 
	_pl[i] = Pol_Info(_fl[i]);
	if (_pl[i].pol_type==p_pl[i].pol_type) {
	  _pl[i] = p_pl[i];
	}
	else {
	  if (p_pl[i].GetPol()!=' ' && p_pl[i].GetPol()!='s') {
	    msg.Out()<<" WARNING: wrong polarisation state in Particle.dat !!!!"<<endl;
	    msg.Out()<<"          Polarisation ignored. "<<endl;
	  }
	}
      }
      else {
	_fl[i] = p_fl[i];
	_pl[i] = p_pl[i];
      }
    }
    overflow = SetPolarisations(plindex,_pl,beam_is_poled);
    for (int i=0;i<m_nin+m_nout;++i) {
      if (plindex[i]!=' ') _pl[i].SetPol(plindex[i]);
    }
    GenerateNames(m_nin,_fl,_pl,m_nout,_fl+m_nin,_pl+m_nin,_name,_stan,_oli);
    take = 1;
    for (int k=0;k<m_procs.size();k++) {
      if (_name == m_procs[k]->Name()) { take = 0; break; }
    }
    if (take) {
      if (CheckExternalFlavours(m_nin,_fl,m_nout,_fl+m_nin)) {
	Add(new Single_Process(m_nin,m_nout,_fl,p_isr,p_beam,_seldata,m_gen_str,m_orderQCD,m_orderEW,
			       m_kfactorscheme,m_scalescheme,m_scalefactor,m_asscale,_pl,m_nex,p_ex_fl));
      }
      else {
	take=0;
      }
    }
    if (overflow || take==0) {
      for (int i=0; i<m_nin+m_nout; ++i) plindex[i]=' ';
      for (int i=m_nin+m_nout-1;i>=0;--i) {
	if (p_fl[i].Size()-1>flindex[i]) {
	  ++flindex[i];
	  break;
	}
	else {
	  if (i==0) flag = 0;
	  flindex[i] = 0;
	}
      }
    }
  }
  delete [] _fl;
  delete [] _pl;
}

int Process_Group::SetPolarisations(char * plindex, Pol_Info * pl, int * beam_is_poled) 
{
  for (int i=m_nin;i<m_nin+m_nout;++i) {
    if (pl[i].DoFNumber()==1) {
      plindex[i]=pl[i].GetPol();
    }
  }

  for (int i=0;i<m_nin;++i) {
    if (!beam_is_poled[i])   plindex[i]=pl[i].GetPol();
  }

  for (int i=0;i<m_nin;++i) {
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

  for (int i=0;i<m_nin;i++)  {
    massin[i]   = p_flin[i].Mass();
    sum_massin += massin[i];
  }
  for (int i=0;i<m_nout;i++) {
    massout[i]   = p_flout[i].Mass();
    sum_massout += massout[i];
  }

  bool massok = 1;
  for (int i=0;i<m_procs.size();i++) {
    for (int j=0;j<m_procs[i]->Nin();j++) {
      if (!(ATOOLS::IsEqual(massin[j],(m_procs[i]->Flavs()[j]).Mass()))) {
	msg.Error()<<"Error in Incoming masses ; "<<massin[j]<<" vs. "<<(m_procs[i]->Flavs()[j]).Mass()
		   <<" for "<<p_flin[j]<<" "<<m_procs[i]->Flavs()[j]<<endl;
	massok = 0; break;
	}
    }
    if (!massok) break;
    for (int j=0;j<m_procs[i]->Nout();j++) {
      if (!(ATOOLS::IsEqual(massout[j],(m_procs[i]->Flavs()[j+m_procs[i]->Nin()]).Mass()))) {
	msg.Error()<<"Error in outgoing masses ; "<<massout[j]<<" vs. "<<(m_procs[i]->Flavs()[j+m_procs[i]->Nin()]).Mass()
		   <<"for "<<p_flout[j]<<" "<<m_procs[i]->Flavs()[j+m_procs[i]->Nin()]<<endl;
	massok = 0; break;
      }
    }
    if (!massok) break;
  }
  if (massok) {
    SetISRThreshold(ATOOLS::Max(sum_massin,sum_massout));
  }
  else {
    msg.Error()<<"Error in Process_Group : "<<m_name<<endl
	       <<"   Processes do not have equal masses. Abort."<<endl;
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
  for (int i=0;i<singleprocs.size();i++) {
    sproc      = singleprocs[i];
    help       = string("SG_");
    if (sproc->Nin()==2) {
      flav1    = sproc->Flavs()[0]; 
      flav2    = sproc->Flavs()[1]; 
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
      flav1    = sproc->Flavs()[0]; 
      if (flav1.IsVector())        help += string("V_->_");
      else if (flav1.IsFermion())  help += string("f_->_");
      else                         help += string("S_->_");
    }
    int scalars = 0,fermions = 0,vectors = 0;
    for (int j=0;j<sproc->Nout();j++) {
      if ((sproc->Flavs()[sproc->Nin()+j]).IsScalar())  scalars++;
      if ((sproc->Flavs()[sproc->Nin()+j]).IsFermion()) fermions++;
      if ((sproc->Flavs()[sproc->Nin()+j]).IsVector())  vectors++;
    }
    char numb[20];
    sprintf(numb,"%i",scalars);
    help += string(numb) + string("S_");
    sprintf(numb,"%i",fermions);
    help += string(numb) + string("F_");
    sprintf(numb,"%i",vectors);
    help += string(numb) + string("V");


    bool found = 0;
    for (int j=0;j<m_procs.size();j++) {
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
      group->SetBeam(p_beam);
      group->SetISR(p_isr);
      group->Add(sproc);
      m_procs.push_back(group);
    }
  }

  for (int i=0;i<m_procs.size();i++) {
    msg.Tracking()<<"Process_Group "<<m_procs[i]->Name()<<" : "<<m_procs[i]->Size()<<endl;
    for (int j=0;j<m_procs[i]->Size();j++) 
      msg.Tracking()<<"    "<<((*m_procs[i])[j])->Name()<<endl;
    msg.Tracking()<<"--------------------------------------------------"<<endl;
  }
}

void Process_Group::Add(Process_Base * _proc) 
{
  if (m_procs.size()==0) {
    m_nin     = _proc->Nin();
    m_nout    = _proc->Nout();
    m_nvec    = _proc->Nvec();
    m_nstrong = _proc->NStrong();
    m_neweak  = _proc->NEWeak();
    if (p_fl==NULL) {
      p_fl = new Flavour[m_nin+m_nout];
      for (int i=0;i<m_nin+m_nout;i++) p_fl[i] = (_proc->Flavs())[i];
    }
  }
  else {
    if (_proc->Nvec() > m_nvec) m_nvec = _proc->Nvec();
  }
  if ( (m_nin != _proc->Nin()) || (m_nout != _proc->Nout())) {
    msg.Error()<<"Error : Cannot add process "<<_proc->Name()
	       <<" to group "<<m_name<<" ! "<<endl
	       <<"   Inconsistent number of external legs."<<endl
	       <<"  Before : ("<<m_nin<<" -> "<<m_nout<<" )"<<endl
	       <<"  Now    : ("<<_proc->Nin()<<" -> "<<_proc->Nout()<<" )"<<endl;
    return;
  }
  m_procs.push_back(_proc);
}

bool Process_Group::Find(string _name,Process_Base *& _proc) 
{
  if (m_name==_name) {
    _proc = this;
    return 1;
  }
  for (int i=0;i<m_procs.size();i++) {
    if (m_procs[i]->Find(_name,_proc)) return 1;
  }
  return 0;
}

void Process_Group::WriteOutXSecs(std::ofstream & _to)
{
  _to<<m_name<<"  "<<m_totalxs<<"  "<<m_max<<"  "<<m_totalerr<<endl;
  for (int i=0;i<m_procs.size();i++) m_procs[i]->WriteOutXSecs(_to);
}


void Process_Group::SelectOne()
{
  DeSelect();
  if (m_totalxs==0) p_selected = m_procs[int(ran.Get()*m_procs.size())];
  else {
    double disc;
    if (m_atoms) {
      // select according to total xsecs.
      disc = m_totalxs * ran.Get();
      for (int i=0;i<m_procs.size();i++) {
	disc -= m_procs[i]->Total();
	if (disc<0.) {
	  p_selected = m_procs[i];
	  p_selected->SelectOne();
	  return;
	}
      }
      if (disc>0.) { 
	msg.Error()<<"Error in Process_Group::SelectOne() : "
		   <<"Total xsec, max = "<<m_totalxs<<", "<<m_max<<endl;
	return;
      }
    }
    else {
      double m=0.;
      for (int i=0;i<m_procs.size();i++) {
	m+= m_procs[i]->Max();
      }
      if (!ATOOLS::IsEqual(m,m_max)) {
	SetMax(0.);
      }
      disc = m_max * ran.Get();
      for (int i=0;i<m_procs.size();i++) {
	disc -= m_procs[i]->Max();
	if (disc<0.) {
	  p_selected = m_procs[i];
	  p_selected->SelectOne();
	  return;
	}
      }
      if (disc>0.) { 
	msg.Error()<<"Error in Process_Group::SelectOne() : "
		   <<"Total xsec, max = "<<m_totalxs<<", "<<m_max<<endl;
	return;
      }
    }
  }
}



void Process_Group::DeSelect() {
  p_selected = 0;
  for (int i=0;i<m_procs.size();i++) m_procs[i]->DeSelect();
}



Process_Base * Process_Group::Selected() { 
  if (p_selected==this) return this;
  return p_selected->Selected(); 
}    



void Process_Group::Empty() {
  for (int i=0;i<m_procs.size();i++) m_procs[i]->Empty();
}



void Process_Group::SetResDir(std::string _resdir) {
  m_resdir = _resdir;
  for (int i=0;i<m_procs.size();i++) m_procs[i]->SetResDir(m_resdir);
}


void Process_Group::SetScale(double _scale)
{
  Process_Base::SetScale(_scale);
  for (int i=0;i<m_procs.size();i++) m_procs[i]->SetScale(_scale); 
} 


void Process_Group::SetISRThreshold(double _isrth)
{
  m_isrthreshold = _isrth;
  for (int i=0;i<m_procs.size();i++) m_procs[i]->SetISRThreshold(m_isrthreshold); 
} 


void Process_Group::SetTables(bool _tables)
{
  m_tables = _tables;
  for (int i=0;i<m_procs.size();i++) m_procs[i]->SetTables(m_tables);
} 

void Process_Group::SetTotalXS(int tables)  { 
  if (tables!=2) {
    m_totalxs  = m_totalsum/m_n; 
    m_totalerr = sqrt( (m_totalsumsqr/m_n - 
			(ATOOLS::sqr(m_totalsum)-m_totalsumsqr)/(m_n*(m_n-1.)) )  / m_n); 
    if ((m_nin==1 && m_nout==2) || m_n==1) m_totalerr = 0.;
    if (p_sel) p_sel->Output();
    m_max = 0.;
    for (int i=0;i<m_procs.size();i++) {
      m_procs[i]->SetTotalXS(tables);
      m_max += m_procs[i]->Max(); // naive sum, probably unneccessary large
    }
  }
  else {
    //   _tables==2  means  check xs with sum of subprocesses
    //               update maximum to sum of maximum
    SetMax(0.);
  }
  msg.Events()<<"-----------------------------------------------------------------------"<<endl;
  if (m_nin==2) 
    msg.Events()<<"Total XS for "<<m_name<<"("<<m_procs.size()<<") : "
		<<m_totalxs*rpa.Picobarn()<<" pb"<<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl;
  
  if (m_nin==1) 
    msg.Events()<<"Total Width for "<<m_name<<"("<<m_procs.size()<<") : "
		<<m_totalxs<<" GeV"<<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl;
}

void Process_Group::SetMax(double max) {
  if (max>0.) {
    m_max=max;
    return;
  }
  // paramter is dummy!
  double sum = 0.;
  m_max = 0.;
  for (int i=0;i<m_procs.size();i++) {
    m_procs[i]->SetTotalXS(2);
    sum   += m_procs[i]->Total();
    m_max += m_procs[i]->Max(); // naive sum, probably unneccessary large
  }
  if (m_totalxs!=0.) {
    if (!ATOOLS::IsEqual(sum,m_totalxs)) {
      msg.Events().precision(12);
      msg.Events()<<" WARNING: group "<<Name()<<": xs and sum of daughters does not agree ! "<<endl
		  <<" sum="<<sum<<"  total:"<<m_totalxs
		  <<"  ("<<((sum-m_totalxs)/m_totalxs)<<")"<<endl;
    }
    m_totalxs=sum;
  }
}

void Process_Group::SetMaxJetNumber(int max) {
  for (int i=0;i<m_procs.size();i++) {
    m_procs[i]->SetMaxJetNumber(max);
  }  
  m_maxjetnumber = max;
}

void Process_Group::SetAtoms(bool _atoms) { m_atoms = _atoms; }


/*----------------------------------------------------------------------------------
  
  Initialization of the processes

  ----------------------------------------------------------------------------------*/

int Process_Group::InitAmplitude(Interaction_Model_Base * model,Topology * top,Vec4D *& testmoms,
				 vector<Single_Process *> & links,vector<Single_Process *> & errs,
				 int & totalsize, int & procs)
{
  int okay = 1;
  vector <string> deletethem;

  for (int i=0;i<m_procs.size();i++) {
    msg.Debugging()<<"========================================================="<<endl
		   <<"========================================================="<<endl
		   <<"Process_Group::InitAmplitude for "<<m_procs[i]->Name()<<endl;
    if (m_atoms) { delete [] testmoms; testmoms = 0; }

    switch (m_procs[i]->InitAmplitude(model,top,testmoms,links,errs,totalsize,procs)) {
    case -3 :
      msg.Debugging()<<"Amplitude is zero: "<<m_procs[i]->Name()<<endl
		     <<"   delete it."<<endl;
      deletethem.push_back(m_procs[i]->Name());
      break;
    case -2 : 
      msg.Error()<<"Error in creation of amplitude "<<m_procs[i]->Name()<<endl;
      return -2;
    case -1 : 
      msg.Debugging()<<"No diagrams or amplitudes for "<<m_procs[i]->Name()<<endl
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
    for (int i=0;i<m_procs.size();i++) sum+=m_procs[i]->Result();

    for (int i=0;i<errs.size();i++) 
      if (!ATOOLS::IsZero(errs[i]->Result()/sum)) flag = 0;
    if (!flag) return -2;
    
    //delete
    for (int i=0;i<errs.size();i++) 
      for (int j=0;j<m_procs.size();j++) 
	if (errs[i]->Name() == m_procs[j]->Name()) {
	  deletethem.push_back(m_procs[j]->Name());
	  msg.Out()<<"Faulty process "<<m_procs[j]->Name()<<" is negligible"<<endl
		   <<"   delete it."<<endl;
	}
    errs.clear();
  }

  for (int i=0;i<m_procs.size();i++) {
    flag = 0;
    if (m_procs[i]->Size() == 0) flag = 1;
    if (!flag) {
      for (int j=0;j<deletethem.size();j++) {
	if (m_procs[i]->Name() == deletethem[j]) flag = 1;
      }
    }
    if (flag) {
      delete m_procs[i];
      for (int j=i;j<m_procs.size()-1;j++) m_procs[j] = m_procs[j+1];
      //just to continue at position i in the next try
      i--;
      m_procs.pop_back();
    }
  }

  if (okay==0) {
    links.clear();
    for (int i=0;i<m_procs.size();i++) if (m_procs[i]) delete (m_procs[i]); 
    m_procs.clear();
  }

  msg.Debugging()<<"Process_Group::Initialize Amplitude for "<<m_name;
  if (okay) msg.Debugging()<<" successful."<<endl;
       else msg.Debugging()<<" failed."<<endl;
  return okay;
}



bool Process_Group::SetUpIntegrator()
{
  bool okay = 1;
  if (m_atoms) {
    for (int i=0;i<m_procs.size();i++) {
      if (m_procs[i]->Partner()==NULL) {
	if (!(m_procs[i]->SetUpIntegrator())) okay = 0;
      }
    }
    return okay;
  }
  
  p_sel->BuildCuts(p_cuts);
  if (m_nin==2) {
    if ( (p_fl[0].Mass() != p_isr->Flav(0).Mass()) ||
	 (p_fl[1].Mass() != p_isr->Flav(1).Mass()) ) p_isr->SetPartonMasses(p_fl);
  }
  p_ps  = new Phase_Space_Handler(this,p_isr,p_beam);
  AddChannels(this,p_ps->FSRIntegrator(),p_ps->BeamParameters(),p_ps->ISRParameters());
  if (!p_ps->CreateIntegrators()) return 0;
  if (m_nin==2) { for (int i=0;i<m_procs.size();i++) m_procs[i]->Empty(); }
  return 1;
}




/*----------------------------------------------------------------------------------
  
  Evaluation of the processes included in All_Processes

  ----------------------------------------------------------------------------------*/

bool Process_Group::CalculateTotalXSec(std::string _resdir)
{
  msg.Tracking()<<"Process_Group::CalculateTotalXSec("<<_resdir<<")"<<endl;
  if (m_atoms) {
    bool okay = 1;
    for (int i=0;i<m_procs.size();i++) {
      if (!(m_procs[i]->CalculateTotalXSec(_resdir))) okay = 0;
    }
    return okay;
  }
  else {
    char filename[100];
    sprintf(filename,"%s.xstotal",(_resdir+string("/")+m_name).c_str());
    string _name;
    double _totalxs,_totalerr,_max;
    if (_resdir!=string("")) {
      if (IsFile(filename)) {
	ifstream from;
	bool okay=1;
	from.open(filename);
	while (from) {
	  from>>_name>>_totalxs>>_max>>_totalerr;
	  if (_name==m_name) m_totalxs += _totalxs;
	  msg.Events()<<"Found result : xs for "<<_name<<" : "
		      <<_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
		      <<" +/- "<<_totalerr/_totalxs*100.<<"%,"<<endl
		      <<"       max : "<<_max<<endl;
	  Process_Base * _proc = NULL;
	  if (Find(_name,_proc)) {
	    _proc->SetTotal(_totalxs);
	    _proc->SetMax(_max);
	  }
	  else {
	    okay = 0;
	  }
	}
	from.close();
	p_ps->ReadIn(_resdir+string("/MC_")+m_name);

	if (p_ps->BeamIntegrator() != 0) p_ps->BeamIntegrator()->Print();
	if (p_ps->ISRIntegrator() != 0)  p_ps->ISRIntegrator()->Print();
	if (p_ps->FSRIntegrator() != 0)  p_ps->FSRIntegrator()->Print();
	if (m_totalxs>0.) {
	  if (okay) {
	    msg.Debugging()<<"In "<<m_name<<"::CalculateTotalXSec("<<_resdir<<")"<<endl
			   <<"   Found all xsecs. Continue"<<endl;
	    SetTotalXS(2);
	    return 1;
	  }
	}
	else {
	  msg.Error()<<"In "<<m_name<<"::CalculateTotalXSec("<<_resdir<<")"<<endl
		     <<"   Something went wrong : Negative xsec : "<<m_totalxs<<endl;
	  return 0;
	}
      }
    }

    if (m_nin==2) {
      if ( (p_fl[0].Mass() != p_isr->Flav(0).Mass()) ||
	   (p_fl[1].Mass() != p_isr->Flav(1).Mass()) ) p_isr->SetPartonMasses(p_fl);
    }
    p_sel->BuildCuts(p_cuts);
    m_tables  = 0;
    
    m_totalxs = p_ps->Integrate();
    if (m_nin==2) m_totalxs /= ATOOLS::rpa.Picobarn(); 
    if (!(ATOOLS::IsZero((m_n*m_totalxs-m_totalsum)/(m_n*m_totalxs+m_totalsum)))) {
      msg.Error()<<"Result of PS-Integrator and internal summation do not coincide!"<<endl
		 <<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<m_totalsum/m_n<<endl;
    }
    SetTotalXS(0);
    if (m_totalxs>0.) {
      if (_resdir!=string("")) {
	std::ofstream to;
	to.open(filename,ios::out);
	to.precision(12);
	msg.Events()<<"Store result : xs for "<<m_name<<" : ";
	if (m_nin==2) msg.Events()<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb";
	if (m_nin==1) msg.Events()<<m_totalxs<<" GeV";
	msg.Events()<<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
		    <<"       max : "<<m_max<<endl;
	WriteOutXSecs(to);
	p_ps->WriteOut(_resdir+string("/MC_")+m_name);
	to.close();
      }
      return 1;
    }
  }
  return 0;
}



void  Process_Group::RescaleXSec(double fac) {
  Process_Base::RescaleXSec(fac);
  for (int i=0;i<m_procs.size();i++) {
    m_procs[i]->RescaleXSec(fac);
  }
}

bool Process_Group::LookUpXSec(double ycut,bool calc,string obs) {
  msg.Tracking()<<"Process_Group::LookUpXSec() for "<<m_name<<" in "<<m_resdir<<endl;
  bool okay = 1;
  if (m_atoms) {
    for (int i=0;i<m_procs.size();i++) {
      if (!(m_procs[i]->LookUpXSec(ycut,1,obs))) okay = 0;
    }
    return okay;
  }
  else {
    for (int i=0;i<m_procs.size();i++) {
      if (!(m_procs[i]->LookUpXSec(ycut,0,obs))) okay = 0;
    }
    if (okay) {                    
      m_totalxs = 0; m_max = 0;
      for (int i=0;i<m_procs.size();i++) {
	m_totalxs += m_procs[i]->Total();
	m_max     += m_procs[i]->Max();
      }
      p_ps->ReadIn(m_resdir+string("/MC_")+m_name);
      if (p_ps->BeamIntegrator() != 0) p_ps->BeamIntegrator()->Print();
      if (p_ps->ISRIntegrator() != 0)  p_ps->ISRIntegrator()->Print();
      if (p_ps->FSRIntegrator() != 0)  p_ps->FSRIntegrator()->Print();
      return 1;
    }
    if (calc) {
      if (!(PrepareXSecTables())) return 0;
      okay = 1;
      for (int i=0;i<m_procs.size();i++) {
	if (!(m_procs[i]->LookUpXSec(ycut,0,obs))) okay = 0;
      }
      if (okay) {                    
	m_totalxs = 0; m_max = 0;
	for (int i=0;i<m_procs.size();i++) {
	  m_totalxs += m_procs[i]->Total();
	  m_max     += m_procs[i]->Max();
	}
	return 1;
      }
    }
    return 0;
  }
}

bool Process_Group::PrepareXSecTables()
{
  if (m_atoms) {
    bool okay = 1;
    for (int i=0;i<m_procs.size();i++) {
      if (!(m_procs[i]->PrepareXSecTables())) okay = 0;
    }
    return okay;
  }
  else {
    if (m_nin==2) {
      if ( (p_fl[0].Mass() != p_isr->Flav(0).Mass()) ||
	   (p_fl[1].Mass() != p_isr->Flav(1).Mass()) ) p_isr->SetPartonMasses(p_fl);
    }
    p_sel->BuildCuts(p_cuts);
    m_totalxs = p_ps->Integrate()/ATOOLS::rpa.Picobarn(); 
    if (!(ATOOLS::IsZero((m_n*m_totalxs-m_totalsum)/(m_n*m_totalxs+m_totalsum)))) {
      msg.Error()<<"Result of PS-Integrator and internal summation do not coincide!"<<endl;
      msg.Error()<<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<m_totalsum/m_n<<endl;
    }
    SetTotalXS(1);
    p_ps->WriteOut(m_resdir+string("/MC_")+m_name);
    if (m_totalxs>0.) return 1;
  }
  return 0;
}






void Process_Group::AddPoint(const double value) 
{
  m_n++;
  m_totalsum    += value;
  m_totalsumsqr += value*value;
  if (value>m_max) m_max = value;

  for (int i=0;i<m_procs.size();i++) {
    if (dabs(m_last)>0.) {
      m_procs[i]->AddPoint(value*m_procs[i]->Last()/m_last);
    }
    else {
      m_procs[i]->AddPoint(0.);
    }  
  }
}






double Process_Group::Differential(Vec4D * p)
{
  m_last = 0;
  for (int i=0;i<m_procs.size();i++) {
    m_last += m_procs[i]->DSigma(p,1);
  }
  //ControlOutput(p);

  if ((!(m_last<=0)) && (!(m_last>0))) {
    msg.Error()<<"---- Process_Group::Differential -------------------"<<endl;
    for (int i=0;i<m_nin+m_nout;i++) 
      msg.Error()<<i<<" th Momentum "<<p[i]<<endl;
    PrintDifferential();
  }
  return m_last;
}



double Process_Group::Differential2()
{
  if (p_isr->On()==0) return 0.;
  double tmp = 0.;
  for (int i=0;i<m_procs.size();i++) tmp += m_procs[i]->DSigma2();

  if ((!(tmp<=0)) && (!(tmp>0))) {
    msg.Error()<<"---- Process_Group::Differential -------------------"<<endl;
    PrintDifferential();
  }
  m_last += tmp;
  return tmp;
}



double Process_Group::DSigma(Vec4D * p,bool lookup)
{
  m_last = 0;
  for (int i=0;i<m_procs.size();i++) m_last += m_procs[i]->DSigma(p,lookup);
  //ControlOutput(p);
  return m_last;
}



double Process_Group::DSigma2()
{
  double tmp = 0.;
  for (int i=0;i<m_procs.size();i++) tmp += m_procs[i]->DSigma2();
  m_last += tmp;
  return tmp;
}


bool Process_Group::OneEvent(double _mass) {
  if (m_atoms) {
    SelectOne();
    return p_selected->OneEvent(_mass);
  }
  return p_ps->OneEvent(_mass);
}

bool Process_Group::SameEvent() {
  if (m_atoms) {
    if (p_selected)
      return p_selected->SameEvent();
    msg.Error()<<" ERROR in bool Process_Group::SameEvent() "<<endl;
    return 0;
  }
  return p_ps->SameEvent();
}

ATOOLS::Blob_Data_Base *  Process_Group::WeightedEvent() {
  if (m_atoms) {
    SelectOne();
    return p_selected->WeightedEvent();
  }
  return p_ps->WeightedEvent();
}

ATOOLS::Blob_Data_Base *  Process_Group::SameWeightedEvent() {
  if (m_atoms) {
    return p_selected->SameWeightedEvent();
  }
  return p_ps->SameWeightedEvent();
}




/*----------------------------------------------------------------------------------
  
  Helpers

  ----------------------------------------------------------------------------------*/

void Process_Group::PrintDifferential()
{
  if (!(rpa.gen.Debugging())) return;
  m_last = 0;
  msg.Out()<<"--------------------------------------------------------"<<endl
	   <<"--------------------------------------------------------"<<endl
	   <<"--------------------------------------------------------"<<endl;
  for (int i=0;i<m_procs.size();i++) {
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
  m_asscale    = ATOOLS::sqr(ATOOLS::rpa.gen.Ecms());
  double a_s = as->AlphaS(m_asscale);
  msg.Out()<<"-------- Process_Group : "<<m_name<<" : DSigma -------------------"<<endl
	   <<"         scale = "<<m_asscale<<" = "<<2.*s*t*u/(s*s+u*u+t*t)<<endl
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


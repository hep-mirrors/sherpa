#include "Combined_Selector.H"
#include "Standard_Selector.H"
#include "Jet_Finder.H"
#include "Dipole_Jet_Finder.H"
#include "Cone_Finder.H"
#include "Run_Parameter.H"
#include "Selector_Bias.H"
#include "Variable_Selector.H"
#include "Message.H"
#include "MyStrStream.H"

using namespace ATOOLS;
using namespace std;

Combined_Selector::Combined_Selector(int _nin,int _nout, Flavour * _fl,
				     Selector_Data * _seldata, std::string ycut) {

  m_name  = std::string("Combined_Selector"); 
  m_nin   = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl    = _fl;
  m_count = 0;
  m_smin = 0.;
  m_update = 0;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();

  if (_seldata==NULL) return; 
  int                               type;
  std::vector<int>                  activetypes;
  std::vector<Flavour>              critflavs;
  double                            rmin,rmax;
  int                               help;
  std::string                       helps;
  bool                              init;
  Selector_Base                   * sel;
  std::vector<std::pair<double,double> > bounds;
  for (int i=0;i<_seldata->Size();i++) {
    _seldata->FillData(i,type,critflavs,bounds,help);
    helps=_seldata->GetData(i).helps;
    rmin=bounds.front().first;
    rmax=bounds.front().second;
    init = 1;
    for (size_t j=0;j<activetypes.size();j++) {
      if (type==activetypes[j]) {
	  if (type!=14) m_sels[j]->SetRange(critflavs,rmin,rmax);
	  else          m_sels[j]->SetRange(critflavs,help,rmin,rmax);
	init = 0;
      }
    }
    if (init) {
      int jettype = 0;
      switch (type) {
      case 1 : 
	if (_nin==2) {
	  int instrong(0);
	  for (int j=0;j<_nin;j++) { if (_fl[j].Strong()) instrong++; }
	  if (instrong==0) jettype = 1;
	  if (instrong==1) jettype = 2;
	  if (instrong==2) jettype = 4;
	}
	{
	  std::string ccut(ToString(rmin));
	  if (rpa.gen.Variable("Y_CUT")!="")
	    ccut=rpa.gen.Variable("Y_CUT");
	  else rpa.gen.SetVariable("Y_CUT",ccut);
	  if (ycut.length()>0 && ycut!="-1.") ccut=ycut;
	  Jet_Finder * jf = new Jet_Finder(_nin+_nout,_fl,ccut,jettype);
	  jf->SetDeltaR(ToType<double>(rpa.gen.Variable("DELTA_R")));
	  sel = jf;
	}
	activetypes.push_back(type);
	break;
      case 2 :
	sel = new Cone_Finder(_nin+_nout,_fl,rmin);
	activetypes.push_back(type);
	break;
      case 3 : 
      {
	if (_nin==1) jettype = 0;
	if (_nin==2) {
	  int instrong(0);
	  for (int j=0;j<_nin;j++) { if (_fl[j].Strong()) instrong++; }
	  if (instrong==0) jettype = 1;
	  if (instrong==1) jettype = 2;
	  if (instrong==2) jettype = 4;
	}
	rmin=Max(rmin,ToType<double>(rpa.gen.Variable("Y_CUT")));
	rpa.gen.SetVariable("Y_CUT",ToString(rmin));
	rpa.gen.SetVariable("DELTA_R",ToString(rmax));
	double ccut(ToType<double>(ycut));
	if (ccut>0.) rmin=ccut;
	Dipole_Jet_Finder *djf
	  (new Dipole_Jet_Finder(_nin+_nout,_fl,rmin,dipjet_type::ariadne,
				 dipjet_mode::code(jettype),false));
	sel = djf;
	activetypes.push_back(type);
	break;
      }
      case 11 :
	sel = new Energy_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	activetypes.push_back(type);	
	break;
      case 12 :
	sel = new PT_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	activetypes.push_back(type);
	break;
      case 13 :
	sel = new Rapidity_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	activetypes.push_back(type);
	break;
      case 14 : 
	sel = new Angle_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,help,rmin,rmax);
	activetypes.push_back(type);
	break;
      case 15 : 
	sel = new ET_Selector(_nin,_nout,_fl);
 	sel->SetRange(critflavs,rmin,rmax);
	activetypes.push_back(type);
	break;
      case 16 :
	sel = new PseudoRapidity_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	activetypes.push_back(type);
	break;
      case 21 : 
	sel = new Mass_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	activetypes.push_back(type);
	break;
      case 22 : 
	sel = new Angle_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	activetypes.push_back(type);	
	break;
      case 24 :
	sel = new X_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	activetypes.push_back(type);
	break;
      case 25 :
	sel = new BFKL_PT_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	activetypes.push_back(type);
	break;
      case 26 : 
	sel = new Delta_Eta_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	activetypes.push_back(type);
	break;
      case 27 : 
	sel = new Delta_Phi_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	activetypes.push_back(type);
	break;
      case 28 : 
	sel = new Delta_R_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	activetypes.push_back(type);
	break;
      case 101 :
	sel = new ET_Bias(_nin,_nout,_fl,helps);
	sel->SetRange(critflavs,bounds);
	activetypes.push_back(type);
	break;
      case 102 :
	sel = new PT_Bias(_nin,_nout,_fl,helps);
	sel->SetRange(critflavs,bounds);
	activetypes.push_back(type);
	break;
      case 103 :
	sel = new Eta_Bias(_nin,_nout,_fl,helps);
	sel->SetRange(critflavs,bounds);
	activetypes.push_back(type);
	break;
      case 113 :
	sel = new Delta_Eta_Bias(_nin,_nout,_fl,helps);
	sel->SetRange(critflavs,bounds);
	activetypes.push_back(type);
	break;
      case 114 :
	sel = new Delta_Phi_Bias(_nin,_nout,_fl,helps);
	sel->SetRange(critflavs,bounds);
	activetypes.push_back(type);
	break;
      case 115 :
	sel = new Delta_R_Bias(_nin,_nout,_fl,helps);
	sel->SetRange(critflavs,bounds);
	activetypes.push_back(type);
	break;
      case 116 :
	sel = new Mass_Bias(_nin,_nout,_fl,helps);
	sel->SetRange(critflavs,bounds);
	activetypes.push_back(type);
	break;
      case 126 :
	sel = new Variable_Selector(_nin,_nout,_fl,helps);
	sel->SetRange(critflavs,bounds);
	break;
      default :
	sel = NULL;
	msg_Error()<<"Error in Combined_Selector::Combined_Selector."<<endl
			      <<"  Unknown type of selector-data : "<<type<<endl; 
      }
      m_update += sel->NeedUpdate(); 
      Add(sel);
    }
  }

  /*for (int i=0;i<m_sels.size();i++) {
    if (m_sels[i]->Smin()>m_smin) m_smin = m_sels[i]->Smin();
    }*/
}

Combined_Selector::~Combined_Selector()
{
  while (m_sels.size()>0) {
    delete *m_sels.begin();
    m_sels.erase(m_sels.begin());
  }
}

int  Combined_Selector::NeedUpdate() 
{ 
  return m_update>0; 
}

void Combined_Selector::Add(Selector_Base * sel) { 
  if (sel != NULL) m_sels.push_back(sel); 
}

bool Combined_Selector::Trigger(const Vec4D* p) 
{
  m_count++;
  if (!(m_count%1000000)) Output(); 
  for (size_t i=0; i<m_sels.size(); ++i) {
    if (!(m_sels[i]->Trigger(p))) return 0;
  }
  return 1;
}

void Combined_Selector::BuildCuts(Cut_Data * cuts)
{
  for (size_t i=0; i<m_sels.size(); ++i) 
    if (!m_sels[i]->IsConditional()) m_sels[i]->BuildCuts(cuts);

  //smin update!!!

  for (size_t i=0; i<m_sels.size(); ++i) 
    if (m_sels[i]->IsConditional()) m_sels[i]->BuildCuts(cuts);
  for (size_t i=0; i<m_sels.size(); ++i) 
    if (m_sels[i]->IsConditional()) m_sels[i]->BuildCuts(cuts);

  for (size_t i=0; i<m_osc.size(); ++i) cuts->Setscut(m_osc[i].first,m_osc[i].second);
  cuts->Complete();
  for (size_t i=0; i<m_osc.size(); ++i) cuts->Setscut(m_osc[i].first,m_osc[i].second);
}

void Combined_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts)
{
  cuts->Reset(NeedUpdate());
  if (NeedUpdate()) 
    for (size_t i=0; i<m_sels.size(); ++i) m_sels[i]->UpdateCuts(sprime,y,cuts);
  cuts->Update(sprime,y);      
  for (size_t i=0; i<m_osc.size(); ++i) cuts->Setscut(m_osc[i].first,m_osc[i].second);
}

void Combined_Selector::AddOnshellCondition(std::string s,double d)
{
  m_osc.push_back(std::pair<std::string,double>(s,d));
}

void Combined_Selector::Output()
{
  msg_Debugging()<<"========================================="<<std::endl
			    <<"Efficiency of the Selector : "<<m_name<<std::endl;
  for (size_t i=0; i<m_sels.size(); ++i) m_sels[i]->Output();
  msg_Debugging()<<"========================================="<<std::endl;
}

Selector_Base * Combined_Selector::GetSelector(std::string name)
{
  for (size_t i=0; i<m_sels.size(); ++i) 
    if (m_sels[i]->Name()==name) return m_sels[i];
  return 0;
  
}

void Combined_Selector::SetProcessName(const std::string &name)   
{ 
  m_procname=name; 
  for (size_t i=0; i<m_sels.size(); ++i)
    m_sels[i]->SetProcessName(m_procname);
}

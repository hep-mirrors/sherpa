#include "Combined_Selector.H"
#include "Standard_Selector.H"
#include "Jet_Finder.H"
#include "Cone_Finder.H"
#include "Message.H"

using namespace APHYTOOLS;
using namespace AMATOOLS;

Combined_Selector::Combined_Selector(int _nin,int _nout, Flavour * _fl,
				     Selector_Data * _seldata) {
  m_name  = std::string("Combined_Selector"); 
  m_nin   = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl    = _fl;
  m_count = 0;

  if (_seldata==NULL) return; 
  int                               type;
  std::vector<int>                  activetypes;
  std::vector<APHYTOOLS::Flavour>   critflavs;
  double                            rmin,rmax;
  int                               help;
  bool                              init;
  Selector_Base                   * sel;

  for (int i=0;i<_seldata->Size();i++) {
    _seldata->Data(i,type,critflavs,help,rmin,rmax);
    init = 1;
    for (int j=0;j<activetypes.size();j++) {
      if (type==activetypes[j]) {
	  if (type!=14) m_sels[j]->SetRange(critflavs,rmin,rmax);
	  else          m_sels[j]->SetRange(critflavs,help,rmin,rmax);
	init = 0;
      }
    }
    if (init) {
      int jettype = 0;
      switch (type) {
      case 1  : 
	if (_nin==2) {
	  int instrong = 0;
	  for (int j=0;j<_nin;j++) { if (_fl[j].Strong()) instrong++; }
	  if (instrong==0) jettype = 1;
	  if (instrong==1) jettype = 2;
	  if (instrong==2) jettype = 4;
	}
	sel = new Jet_Finder(_nin+_nout,_fl,rmin,1,jettype);
	break;
      case 2 :
	sel = new Cone_Finder(_nin+_nout,rmin);
	break;
      case 11 :
	sel = new Energy_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	break;
      case 12 :
	sel = new PT_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	break;
      case 13 :
	sel = new Rapidity_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	break;
      case 14 : 
	sel = new Angle_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,help,rmin,rmax);
	break;
      case 21 : 
	sel = new Mass_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	break;
      case 22 : 
	sel = new Angle_Selector(_nin,_nout,_fl);
	sel->SetRange(critflavs,rmin,rmax);
	break;
      default :
	sel = NULL;
	AORGTOOLS::msg.Error()<<"Error in Combined_Selector::Combined_Selector."<<endl
			      <<"  Unknown type of selector-data : "<<type<<endl; 
      } 
      Add(sel);
    }
  }
}

Combined_Selector::~Combined_Selector()
{
  int count=0;
  for(int i=0; i<m_sels.size(); ++i) {
    if (m_sels[i-count]) delete m_sels[i-count];
    count++;
  }
}

void Combined_Selector::Add(Selector_Base * sel) { 
  if (sel != NULL) m_sels.push_back(sel); 
}

bool Combined_Selector::Trigger(const Vec4D* p) 
{
  m_count++;
  if (!(m_count%1000000)) Output(); 
  for (short int i=0; i<m_sels.size(); ++i) {
    if (!(m_sels[i]->Trigger(p))) return 0;
  }
  return 1;
}

void Combined_Selector::BuildCuts(Cut_Data * cuts)
{
  for (short int i=0; i<m_sels.size(); ++i) m_sels[i]->BuildCuts(cuts);
}

void Combined_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts)
{
  for (short int i=0; i<m_sels.size(); ++i) m_sels[i]->UpdateCuts(sprime,y,cuts);
}

void Combined_Selector::Output()
{
  AORGTOOLS::msg.Tracking()<<"========================================="<<std::endl
			   <<"Efficiency of the Selector : "<<m_name<<std::endl;
  for (short int i=0; i<m_sels.size(); ++i) m_sels[i]->Output();
  AORGTOOLS::msg.Tracking()<<"========================================="<<std::endl;
}

#include "Full_Decay_Table.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Full_Decay_Channel::Full_Decay_Channel(Decay_Channel & _dec) :
  p_proc(NULL)
{
  m_dec = _dec;
}


Full_Decay_Channel::Full_Decay_Channel(const Flavour _fl) :
  p_proc(NULL)
{
  m_dec =  Decay_Channel(_fl);
}

void Full_Decay_Channel::AddDecayProduct(const Flavour & _fl) 
{
  m_dec.AddDecayProduct(_fl);
}

void Full_Decay_Channel::Output()
{
  FlavourSet decs = m_dec.GetDecayProducts();
  msg.Out()<<"  "<<(m_dec.GetDecaying())<<" -> ";
  for (FlSetIter fl=decs.begin();fl!=decs.end();++fl) msg.Out()<<(*fl)<<" ";
  msg.Out()<<" : "<<m_dec.Width()<<" GeV.";
  if (p_proc) msg.Out()<<"  : "<<p_proc->Name();
  msg.Out()<<endl;
}
 

bool Full_Decay_Channel::CreateDecay()
{
  FlavourSet decs = m_dec.GetDecayProducts();
  Flavour * flavs = new Flavour[1+decs.size()];
  flavs[0] = (m_dec.GetDecaying());
  int i    = 0;
  for (FlSetIter fl=decs.begin();fl!=decs.end();++fl) {
    i++;
    flavs[i] = (*fl);
  }
  if (decs.size()==2) p_proc = new Single_Process(1,decs.size(),flavs);
                 else p_proc = new Single_Process(1,decs.size(),flavs,NULL,NULL,NULL,2);
  m_dec.SetProcessName(p_proc->Name());
}

bool Full_Decay_Channel::CalculateWidth() 
{
}


void Full_Decay_Channel::SetWidth(double _w ) 
{
  if (_w<0.) m_dec.SetWidth(p_proc->Total()); 
  else m_dec.SetWidth(_w);
  if (rpa.gen.Events()) Output();
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

Full_Decay_Table::Full_Decay_Table(const Flavour _fl) :
  m_flin(_fl), m_width(0.), m_isevaluated(0) { }


void Full_Decay_Table::AddDecayChannel(Decay_Channel * _dc)
{
  if (m_isevaluated) return;
  m_channels.push_back(new Full_Decay_Channel(*_dc));
  m_width += _dc->Width();
}
void Full_Decay_Table::AddDecayChannel(Full_Decay_Channel * _dc)
{
  if (m_isevaluated) return;
  m_channels.push_back(_dc);
  m_width += _dc->Width();
}




bool Full_Decay_Table::InitAllDecays(Interaction_Model_Base * _model,Topology * _top)
{
  if (m_isevaluated) return 1;
  Vec4D * moms = NULL;
  vector<Single_Process *> links,errs;
  int totalsize = 0;
  int procs     = 0;

  bool okay = 1;
  for (int i=0;i<m_decaymodes.size();i++) {
    msg.Tracking()<<"============================================================"<<endl;
    if (moms) { delete [] moms; moms = NULL; }
    links.clear();
    okay = okay && m_decaymodes[i]->InitAmplitude(_model,_top,moms,links,errs,totalsize,procs);
    for (int j=0;j<links.size();j++) {
      msg.Tracking()<<"Set up integrator of "<<j<<" : "<<links[j]->Name()<<endl;
      if (!(links[j]->SetUpIntegrator())) okay = 0;
    }
    if (okay) {
      msg.Tracking()<<"Set Selector of "<<i<<" : "<<(*m_decaymodes[i])[0]->Name()<<endl;
      m_decaymodes[i]->SetSelector((*m_decaymodes[i])[0]->Selector());
      okay = okay && m_decaymodes[i]->SetUpIntegrator();
    }
  }
  return okay;
}

void Full_Decay_Table::ArrangeDecays()
{
  if (m_isevaluated) return;
  Process_Group * group;
  Process_Base  * proc, * test;
  int             nout;
  Flavour       * flavs;
  bool newgroup = 1;
  bool alreadyin,match,initgroup;
  while (newgroup) {
    newgroup = 0;
    for (int i=0;i<m_channels.size();i++) {
      proc     = m_channels[i]->GetProcessBase();
      nout     = proc->Nout();
      initgroup = 1;
      if (m_decaymodes.size()>0) {
	for (int j=0;j<m_decaymodes.size();j++) {
	  alreadyin = 0;
	  for (int k=0;k<m_decaymodes[j]->Size();k++) {
	    if (proc->Name()==(*m_decaymodes[j])[k]->Name()) {
	      alreadyin = 1;
	      break;
	    } 
	  }
	  if (alreadyin) { initgroup = 0; break; }
	  if (m_decaymodes[j]->Nout()==nout) {
	    match = 1;
	    for (int k=0;k<m_decaymodes[j]->Nout();k++) {
	      if ((((*m_decaymodes[j])[0]->Flavs())[k+1]).Mass()!=((proc->Flavs())[k+1]).Mass()) {
		match = 0;
		break;
	      }
	    }
	    if (match==1) {
	      m_decaymodes[j]->Add(proc);
	      initgroup = 0;
	      break;
	    }
	  }
	}
      }
      if (initgroup==1) {
	newgroup = 1;
	flavs = proc->Flavs();
	group = new Process_Group();
	group->SetName("Decay for : "+flavs[0].TexName());
	group->Add(proc);
	group->SetAtoms(0);
	m_decaymodes.push_back(group);
	break;
      }
    }
  }
}



void Full_Decay_Table::CalculateWidths()
{
  if (m_isevaluated) return;
  for (int i=0;i<m_decaymodes.size();i++) {
    msg.Tracking()<<"Full_Decay_Table::CalculateWidths for "<<m_decaymodes[i]->Size()<<" decay(s)."<<endl;
    m_decaymodes[i]->CalculateTotalXSec(string(""));
  }
  for (int i=0;i<m_channels.size();i++) {
    m_channels[i]->SetWidth();
    m_width += m_channels[i]->Width();
  }
  m_isevaluated = 1;
  Output();
  m_flin.SetWidth(m_width);
}

Decay_Channel * Full_Decay_Table::GetChannel(int _ch) 
{
  if (_ch<0 || _ch>=m_channels.size()) {
    msg.Error()<<"Error in Full_Decay_Table::Channel("<<_ch<<")."<<endl
	       <<"    Out of bounds : 0 ... "<<m_channels.size()-1<<"."<<endl
	       <<"    Return NULL."<<endl;
    return new Decay_Channel();
  }
  return m_channels[_ch]->GetDecayChannel();
}

Full_Decay_Channel * Full_Decay_Table::GetFullChannel(int _ch) 
{
  if (_ch<0 || _ch>=m_channels.size()) {
    msg.Error()<<"Error in Full_Decay_Table::GetFullChannel("<<_ch<<")."<<endl
	       <<"    Out of bounds : 0 ... "<<m_channels.size()-1<<"."<<endl
	       <<"    Return NULL."<<endl;
    return NULL;
  }
  return m_channels[_ch];
}

void Full_Decay_Table::Output() {
  msg.Out()<<"Decay table for : "<<m_flin<<", total width : "<<m_width<<" GeV ("
	   <<m_flin.Width()<<" GeV)."<<endl
	   <<"----------------------------------------------------------------"<<endl;
  for (int i=0;i<m_channels.size();i++) m_channels[i]->Output();
  msg.Out()<<"----------------------------------------------------------------"<<endl;
}






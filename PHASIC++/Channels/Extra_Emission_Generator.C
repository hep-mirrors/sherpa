#include "PHASIC++/Channels/Extra_Emission_Generator.H"

#include "PHASIC++/Channels/CS_Dipoles.H"
#include "PHASIC++/Process/POWHEG_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Integration_Info.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Data_Writer.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace ATOOLS;
using namespace PHASIC;

const size_t s_noptmin(10);

Extra_Emission_Generator::Extra_Emission_Generator():
  m_opt(5), m_weight(0.0), m_asum(0.0), m_numtrig(false)
{
  Data_Reader read(" ",";","!","=");
  read.SetInputPath(rpa->GetPath());
  read.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
  if (!read.ReadFromFile(m_omode,"EEG_OMODE")) m_omode=3;
  else msg_Info()<<METHOD<<"(): Set mode "<<m_omode<<".\n";
  if (!read.ReadFromFile(m_opt,"EEG_OSTEP")) m_opt=5;
  else msg_Info()<<METHOD<<"(): Set steps "<<m_opt<<".\n";
  if (!read.ReadFromFile(m_amin,"EEG_AMIN")) m_amin=1.0e-6;
  else msg_Info()<<METHOD<<"(): Set \\alpha_{min} = "<<m_amin<<".\n";
  if (!read.ReadFromFile(m_Q2min,"EEG_Q2MIN")) m_Q2min=1.0e-6;
  else msg_Info()<<METHOD<<"(): Set Q^2_{min} = "<<m_Q2min<<".\n";
}

Extra_Emission_Generator::~Extra_Emission_Generator() 
{
  for (size_t i(0);i<m_dipoles.size();++i) delete m_dipoles[i];
}

bool Extra_Emission_Generator::AddDipole
(POWHEG_Process *const proc,CS_Dipole *const dip)
{
  for (size_t i(0);i<m_dipoles.size();++i)
    if (dip->IsMapped(m_dipoles[i])) {
      delete dip;
      return false;
    }
  dip->InitVegas("");
  dip->SetAMin(m_amin);
  dip->SetQ2Min(m_Q2min);
  m_dipoles.push_back(dip);
  dip->SetIdx(proc->DipoleMap()->find(dip->GetSubEvt()->p_id)->second);
  m_dmap[dip->Idx()]=dip;
  return true;
}

bool Extra_Emission_Generator::InitDipoles
(POWHEG_Process *const proc,Process_Base *const sproc,
 Phase_Space_Handler *const psh)
{
  p_proc=proc;
  p_info=psh->GetInfo();
  m_nin=sproc->NIn();
  for (size_t i(0);i<sproc->Size();++i) {
    NLO_subevtlist *subs((*sproc)[i]->GetSubevtList());
    for (size_t j(0);j<subs->size()-1;++j) {
      NLO_subevt *sub((*subs)[j]);
      if (sub->m_i<m_nin) {
        if (sub->m_k<m_nin) AddDipole(proc,new II_Dipole(sub,this,psh));
        else AddDipole(proc,new IF_Dipole(sub,this,psh));
      }
      else {
        if (sub->m_k<m_nin) AddDipole(proc,new FI_Dipole(sub,this,psh));
        else AddDipole(proc,new FF_Dipole(sub,this,psh));
      }
    }
  }
  for (size_t i(0);i<m_dipoles.size();++i)
    m_dipoles[i]->SetAlpha(1.0/m_dipoles.size());
  // output dipoles
  if (msg_LevelIsDebugging()) {
    DEBUG_FUNC("");
    for (size_t i=0;i<m_dipoles.size();++i)
	msg_Debugging()<<*m_dipoles[i]<<std::endl;
  }
  return true;
}

Vec4D_Vector Extra_Emission_Generator::GeneratePoint
(const Vec4D_Vector &p,Cut_Data *const cuts)
{
  DEBUG_FUNC("");
  double rns[4];
  for (size_t i(0);i<4;++i) rns[i]=ran->Get();
  msg_Debugging()<<"in EEG: ";
  msg_Debugging()<<"#1 = "<<rns[1]<<", #2 = "<<rns[2]
                 <<", phi = "<<rns[3]<<"\n";
  msg_Debugging()<<"Born point {\n";
  for (size_t j(0);j<p.size();++j)
    msg_Debugging()<<"  "<<p[j]<<"\n";
  msg_Debugging()<<"}\n";
  m_asum=0.0;
  CSDipole_Vector cdips;
  for (size_t i(0);i<m_dipoles.size();++i) {
    if (m_dipoles[i]->ValidPoint(p)) {
      cdips.push_back(m_dipoles[i]);
      m_asum+=cdips.back()->Alpha(1);
      m_dipoles[i]->SetOn(true);
    }
    else {
      m_dipoles[i]->SetOn(false);
      msg_Debugging()<<"invalid for "<<m_dipoles[i]->Id()<<"\n";
    }
  }
  if (!(m_numtrig=cdips.size())) return Vec4D_Vector();
  p_active=NULL;
  double disc(rns[0]*m_asum), psum(0.0);
  for (size_t i(0);i<cdips.size();++i)
    if ((psum+=cdips[i]->Alpha(1))>=disc) {
      p_active=cdips[i];
      break;
    }
  if (p_active==NULL) THROW(fatal_error,"Internal error");
  msg_Debugging()<<"selected "<<p_active->Id()<<"\n";
  return p_active->GeneratePoint(p,cuts,&rns[1]);
}

bool Extra_Emission_Generator::GenerateWeight
(const Vec4D_Vector &p,Cut_Data *const cuts,bool activeonly)
{
  DEBUG_FUNC("");
  m_weight=0.0;
  p_info->ResetAll();
  if (p.empty()) return false;
  double asum(0.0);
  for (size_t i(0);i<m_dipoles.size();++i) {
    CS_Dipole *cdip(m_dipoles[i]);
    if (activeonly || cdip==p_active) continue;
    msg_Debugging()<<"add "<<cdip->Id()<<" {\n";
    double wgt(cdip->GenerateWeight(p,cuts));
    double alpha(cdip->Alpha(1));
    p_info->ResetAll();
    msg_Debugging()<<"} -> w = "<<wgt<<", a = "<<alpha<<"\n";
    if (wgt!=0.0) m_weight+=alpha/wgt;
    asum+=alpha;
  }
  msg_Debugging()<<"add "<<p_active->Id()<<" {\n";
  double wgt(p_active->GenerateWeight(p,cuts));
  double alpha(p_active->Alpha(1));
  msg_Debugging()<<"} -> w = "<<wgt<<", a = "<<alpha<<"\n";
  if (wgt!=0.0) m_weight+=alpha/wgt;
  else {
    msg_Error()<<METHOD<<"(): No weight from active dipole !"<<std::endl; 
    return false;
  }
  asum+=alpha;
  if (IsBad(asum/m_weight))
    msg_Error()<<METHOD<<"(): Bad weight "<<asum
	       <<" / "<<m_weight<<"."<<std::endl;
  m_weight=asum/m_weight;
  msg_Debugging()<<"m_weight = "<<m_weight<<"\n";
  if (!(m_weight>0.)) msg_Out()<<m_weight<<std::endl;
  return true;
}

double Extra_Emission_Generator::SelectionWeight(const size_t &idx) const
{
  return m_dmap.find(idx)->second->Alpha();
}

void Extra_Emission_Generator::AddPoint(const double &value)
{ 
  double bme=p_proc->LastB()+p_proc->LastVI();
  double rme=p_proc->LastRS();
  p_active->AddPoint(value,m_weight,bme,rme,1);
  for (size_t i(0);i<m_dipoles.size();++i) {
    if (m_dipoles[i]==p_active) continue;
    m_dipoles[i]->AddPoint(value,m_weight,bme,rme,0);
  }
}

void Extra_Emission_Generator::Optimize()  
{
  msg_Tracking()<<"Optimize EEG ("<<m_opt<<") {\n";
  size_t off(0);
  {
    msg_Indent();
    bool aopt(true);
    double csum(0.0), wmean(0.0), nc(0.0);
    for (size_t i(0);i<m_dipoles.size();++i)
      if (m_dipoles[i]->N()<s_noptmin) {
	aopt=false;
	msg_Tracking()<<"Too few points in channel "<<i
		      <<" ( "<<m_dipoles[i]->N()<<" vs. "<<s_noptmin
		      <<" ).\nSkip \\alpha optimization in this step.\n";
	break;
      }
    for (size_t i(0);i<m_dipoles.size();++i) {
      CS_Dipole *v(m_dipoles[i]);
      msg_Debugging()<<v->Id()<<" : alpha = "<<v->Alpha()<<std::endl;
      if (v->Alpha()<=0.0) ++off;
      else {
	if (m_opt==1 && (m_omode&2)) v->Optimize();
	if ((m_omode&1) && aopt) {
	  v->SetOldAlpha(v->Alpha());
  	  v->SetAlpha(v->Alpha()*sqrt(dabs(v->Mean())));
	  csum+=v->Alpha();
	  wmean+=dabs(v->Mean());
	  ++nc;
	}
      }
    }
    wmean/=nc;
    if (aopt) {
    msg_Tracking()<<std::string(116,'-')<<"\n";
    for (size_t i(0);i<m_dipoles.size();++i) {
      CS_Dipole *v(m_dipoles[i]);
      if (v->Alpha()>0.0) {
 	if (nc>0.0) v->SetAlpha(v->Alpha()/csum);
	double dev(int((v->Alpha()/v->OldAlpha()-1.0)*10000)/100.0);
	double re(int(v->Sigma()/v->Mean()*10000)/100.0);
	if (v->N()<2) re=100.0;
	if (v->Alpha()) {
	msg_Tracking()<<std::left<<std::setw(14)<<v->Id()
		      <<": n = "<<std::right<<std::setw(5)
		      <<v->N()<<"  w' = "<<std::setw(15)
		      <<(nc>0.0?v->Mean()/wmean:v->Mean())
		      <<" +- "<<std::setw(6)<<re
		      <<" %  =>  a = "<<std::setw(15)<<v->OldAlpha()
		      <<" -> "<<std::setw(15)<<v->Alpha()
		      <<" ( "<<std::setw(6)<<std::right
		      <<dev<<std::left<<" % )\n";
	}
 	v->Reset();
      }
    }
    msg_Tracking()<<std::string(116,'-')<<"\n";
    }
  }
  msg_Tracking()<<"}";
  if (off) msg_Tracking()<<" "<<off<<" channels off";
  msg_Tracking()<<"\n";
  if (m_opt>1) --m_opt;
} 

void Extra_Emission_Generator::EndOptimize()  
{
  for (size_t i(0);i<m_dipoles.size();++i)
    m_dipoles[i]->EndOptimize();
  m_opt=0;
} 

void Extra_Emission_Generator::MPISync()
{
  for (size_t i(0);i<m_dipoles.size();++i)
    m_dipoles[i]->MPISync();
} 

void Extra_Emission_Generator::WriteOut(std::string pid)
{ 
  MakeDir(pid,false);
  pid+="_CS";
  std::vector<std::vector<std::string> > pvds(m_dipoles.size());
  for (size_t i(0);i<m_dipoles.size();++i)
    m_dipoles[i]->WriteOut(pid,pvds[i]);
  pvds.push_back(std::vector<std::string>(1,ToString(m_opt)));
  Data_Writer writer;
  writer.SetOutputPath(pid);
  writer.SetOutputFile("_EEG_PV");
  writer.MatrixToFile(pvds);
}

void Extra_Emission_Generator::ReadIn(std::string pid)
{
  pid+="_CS";
  Data_Reader reader;
  reader.SetAddCommandLine(false);
  reader.SetInputPath(pid);
  reader.SetInputFile("_EEG_PV");
  std::vector<std::vector<std::string> > pvds;
  reader.MatrixFromFile(pvds);
  if (m_dipoles.size()>pvds.size()-1)
    THROW(fatal_error,"Corrupted input file");
  for (size_t i(0);i<m_dipoles.size();++i)
    m_dipoles[i]->ReadIn(pid,pvds[i]);
  if (pvds.back().size()!=1)
    THROW(fatal_error,"Corrupted input file");
  m_opt=ToType<int>(pvds.back().front());
}

void Extra_Emission_Generator::Print()  
{
  msg_Tracking()<<"EEG with "<<m_dipoles.size()<<" dipoles\n";
  for (size_t i(0);i<m_dipoles.size();++i) {
    msg_Tracking()<<"  "<<m_dipoles[i]->Id()<<" : "
		  <<m_dipoles[i]->Alpha()<<"\n";
  }
  msg_Tracking()<<"----------------------------------------------\n";
}

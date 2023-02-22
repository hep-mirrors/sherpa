#include "PHASIC++/Channels/RS_Emission_Generator.H"

#include "PHASIC++/Channels/CS_Dipoles.H"
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
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Data_Writer.H"
#include "PHASIC++/Channels/Vegas.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace ATOOLS;
using namespace PHASIC;

const size_t s_noptmin(10);

RS_Emission_Generator::RS_Emission_Generator():
  m_opt(5)
{
  Default_Reader reader;
  reader.SetInputPath(rpa->GetPath());
  reader.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
  if (reader.Read<size_t>(m_omode,"EEG_OMODE", 2)) {
    msg_Info()<<METHOD<<"(): Set mode "<<m_omode<<"."<<std::endl;
  }
  if (reader.Read<size_t>(m_opt,"EEG_OSTEP", 5)) {
    msg_Info()<<METHOD<<"(): Set steps "<<m_opt<<"."<<std::endl;
  }
  if (reader.Read(m_Q2min,"EEG_Q2MIN", 1.0e-6)) {
    msg_Info()<<METHOD<<"(): Set Q^2_{min} = "<<m_Q2min<<"."<<std::endl;
  }
  reader.CloseInFile(0);
  reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  if (reader.Read(m_amin,"DIPOLE_AMIN",Max(ATOOLS::Accu(),1.0e-8))) {
    msg_Info()<<METHOD<<"(): Set \\alpha_{min} = "<<m_amin<<".\n";
  }
}

RS_Emission_Generator::~RS_Emission_Generator() 
{
  for (size_t i(0);i<m_dipoles.size();++i) delete m_dipoles[i];
}

bool RS_Emission_Generator::AddDipole(CS_Dipole *const dip)
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
  return true;
}

bool RS_Emission_Generator::InitDipoles
(Process_Base *const sproc,Phase_Space_Handler *const psh)
{
  m_nin=sproc->NIn();
  for (size_t i(0);i<sproc->Size();++i) {
    NLO_subevtlist *subs((*sproc)[i]->GetSubevtList());
    for (size_t j(0);j<subs->size()-1;++j) {
      NLO_subevt *sub((*subs)[j]);
      if (sub->m_i<m_nin) {
	if (sub->m_k<m_nin) AddDipole(new II_Dipole(sub,psh,1));
        else AddDipole(new IF_Dipole(sub,psh,1));
      }
      else {
	if (sub->m_k<m_nin) AddDipole(new FI_Dipole(sub,psh,1));
	else AddDipole(new FF_Dipole(sub,psh,1));
      }
    }
  }
  for (size_t i(0);i<m_dipoles.size();++i)
    m_dipoles[i]->SetAlpha(1.0/m_dipoles.size());
  if (msg_LevelIsDebugging()) {
    for (size_t i(0);i<m_dipoles.size();++i) {
      msg_Debugging()<<*m_dipoles[i]<<std::endl;
    }
  }
  return true;
}

bool RS_Emission_Generator::GeneratePoint
(const Vec4D_Vector &p,Cut_Data *const cuts)
{
  DEBUG_FUNC("");
  m_p.clear();
  p_active=NULL;
  double rns[4];
  for (size_t i(0);i<4;++i) rns[i]=ran->Get();
  msg_Debugging()<<"in EEG: ";
  msg_Debugging()<<"#1 = "<<rns[1]<<", #2 = "<<rns[2]
                 <<", #3 = "<<rns[3]<<"\n";
  msg_Debugging()<<"Born point {\n";
  for (size_t j(0);j<p.size();++j)
    msg_Debugging()<<"  "<<p[j]<<"\n";
  msg_Debugging()<<"}\n";
  double asum(0.0);
  CSDipole_Vector cdips;
  for (size_t i(0);i<m_dipoles.size();++i) {
    if (!m_dipoles[i]->ValidPoint(p)) {
      msg_Debugging()<<"invalid for "<<m_dipoles[i]->Id()<<"\n";
      continue;
    }
    cdips.push_back(m_dipoles[i]);
    asum+=cdips.back()->Alpha(1);
  }
  if (cdips.empty()) return false;
  double disc(rns[0]*asum), psum(0.0);
  for (size_t i(0);i<cdips.size();++i)
    if ((psum+=cdips[i]->Alpha(1))>=disc) {
      p_active=cdips[i];
      break;
    }
  if (p_active==NULL) THROW(fatal_error,"Internal error");
  msg_Debugging()<<"selected "<<p_active->Id()<<"\n";
  m_p=p_active->GeneratePoint(p,cuts,&rns[1]);
  return true;
}

bool RS_Emission_Generator::GenerateWeight
(Cut_Data *const cuts,bool activeonly)
{
  DEBUG_FUNC("");
  if (p_active==NULL) {
    msg_Debugging()<<"Invalid Born\n";
    return false;
  }
  double wsum(0), asum(0);
  for (size_t i(0);i<m_dipoles.size();++i) {
    if (!m_dipoles[i]->On() || m_dipoles[i]==p_active) continue;
    msg_Debugging()<<"Dipole "<<m_dipoles[i]->Id()<<" {\n";
    double wgt(m_dipoles[i]->GenerateWeight(m_p,cuts));
    msg_Debugging()<<"} -> w = "<<wgt
  		   <<" ( a = "<<m_dipoles[i]->Alpha(1)<<" )\n";
    wsum+=m_dipoles[i]->Alpha(1)/wgt;
    asum+=m_dipoles[i]->Alpha(1);
  }
  msg_Debugging()<<"Dipole "<<p_active->Id()<<" {\n";
  double wgt(p_active->GenerateWeight(m_p,cuts));
  msg_Debugging()<<"} -> w = "<<wgt
		 <<" ( a = "<<p_active->Alpha(1)<<" )\n";
  wsum+=p_active->Alpha(1)/wgt;
  asum+=p_active->Alpha(1);
  m_weight=asum/wsum;
  return true;
}

void RS_Emission_Generator::AddPoint(const double &value)
{ 
  p_active->AddPoint(value,m_weight,1);
  for (size_t i(0);i<m_dipoles.size();++i)
    if (m_dipoles[i]!=p_active)
      m_dipoles[i]->AddPoint(value,m_weight,0);
}

void RS_Emission_Generator::Optimize()  
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

void RS_Emission_Generator::EndOptimize()  
{
  for (size_t i(0);i<m_dipoles.size();++i)
    m_dipoles[i]->EndOptimize();
  m_opt=0;
} 

void RS_Emission_Generator::MPISync()
{
#ifdef USING__MPI
  size_t i(0), jc(0), jr(0);
  std::vector<double> sv;
  for (size_t i(0);i<m_dipoles.size();++i)
    m_dipoles[i]->MPICollect(sv,jc);
  if (mpi->Size())
    mpi->Allreduce(&sv[0],sv.size(),MPI_DOUBLE,MPI_SUM);
  for (size_t i(0);i<m_dipoles.size();++i)
    m_dipoles[i]->MPIReturn(sv,jr);
#endif
  for (size_t i(0);i<m_dipoles.size();++i)
    m_dipoles[i]->MPISync();
} 

void RS_Emission_Generator::WriteOut(std::string pid)
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

bool RS_Emission_Generator::ReadIn(std::string pid)
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
  return true;
}

void RS_Emission_Generator::Print()  
{
  msg_Tracking()<<"EEG with "<<m_dipoles.size()<<" dipoles\n";
  for (size_t i(0);i<m_dipoles.size();++i) {
    msg_Tracking()<<"  "<<m_dipoles[i]->Id()<<" : "
		  <<m_dipoles[i]->Alpha()<<"\n";
  }
  msg_Tracking()<<"----------------------------------------------\n";
}

// namespace PHASIC
// {
//   std::ostream &operator<<(std::ostream &ostr,const Dipole_Params &dp)
//   {
//     ostr<<*dp.p_dip<<"\n";
//     for (size_t i(0);i<dp.m_procs.size();++i)
//       ostr<<"  "<<dp.m_procs[i]->Name()<<"\n";
//     for (size_t i(0);i<dp.m_p.size();++i)
//       ostr<<"  "<<dp.m_p[i]<<"\n";
//     ostr<<"-> "<<dp.m_weight<<"\n";
//     return ostr;
//   }
// }


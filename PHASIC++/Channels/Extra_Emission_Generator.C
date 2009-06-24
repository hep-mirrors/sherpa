#include "PHASIC++/Channels/Extra_Emission_Generator.H"

#include "PHASIC++/Channels/CS_Dipoles.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Data_Writer.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Channels/Vegas.H"

using namespace ATOOLS;
using namespace PHASIC;

const double s_alphamin(1.0e-10);

Extra_Emission_Generator::Extra_Emission_Generator
(Process_Integrator *const pint, bool isqcd):
  p_int(pint), m_opt(5), m_weight(0.0), m_isqcd(isqcd)
{
  double amin(0.0);
  Data_Reader read(" ",";","!","=");
  read.SetInputFile("Integration.dat");
  if (!read.ReadFromFile(m_omode,"EEG_OMODE")) m_omode=1;
  else msg_Info()<<METHOD<<"(): Set mode "<<m_omode<<".\n";
  if (!read.ReadFromFile(m_opt,"EEG_OSTEP")) m_opt=5;
  else msg_Info()<<METHOD<<"(): Set steps "<<m_opt<<".\n";
  if (!read.ReadFromFile(amin,"EEG_AMIN")) amin=1.0e-6;
  else msg_Info()<<METHOD<<"(): Set \\alpha_{min} = "<<amin<<".\n";
  m_nin=p_int->Process()->NIn();
  if (isqcd){
    InitDipoles(p_int->Process());
  }else{
    InitEWDipoles(p_int->Process());
  }
  SetAMin(amin);
}

Extra_Emission_Generator::~Extra_Emission_Generator() 
{
}

bool Extra_Emission_Generator::AddDipole(CS_Dipole *dip)
{
  for (size_t i=0;i<m_dipoles.size();++i)
    if (dip->LIJ().Id()==m_dipoles[i]->LIJ().Id() &&
      	dip->LK().Id()==m_dipoles[i]->LK().Id() &&
	dip->LI().Id()==m_dipoles[i]->LI().Id() &&
      	dip->LJ().Id()==m_dipoles[i]->LJ().Id()) {
      delete dip;
      return false;
    }
  dip->InitVegas(p_int->Process());
  m_dipoles.push_back(dip);
  return true;
}

bool Extra_Emission_Generator::InitDipoles(Process_Base *const proc)
{
  std::vector<Cluster_Leg> legs;
  for (size_t i(0);i<proc->Flavours().size();++i) {
    Cluster_Leg lij(NULL,Vec4D(),proc->Flavours()[i]);
    if (i<m_nin) {
      lij.SetFlav(lij.Flav().Bar());
      lij.SetStat(1);
    }
    lij.SetId(1<<i);
    for (size_t j(0);j<legs.size();++j) {
      const Cluster_Leg &lk(legs[j]);
      if (lij.Flav().Strong() && lk.Flav().Strong()) {
	for (size_t lid=p_int->Process()->NIn();
	     lid<p_int->Process()->NIn()+p_int->Process()->NOut()+1;++lid) {
	  Cluster_Leg li(NULL,Vec4D(),kf_jet);
	  li.SetId(1<<lid);
	  for (size_t ljd=lid+1;
	       ljd<p_int->Process()->NIn()+p_int->Process()->NOut()+1;++ljd) {
	    Cluster_Leg lj(NULL,Vec4D(),kf_jet);
	    lj.SetId(1<<ljd);
	if (lij.Stat()) {
	  if (lk.Stat()) {//II
	  }
	  else {//IF
	  }
	}
	else {
	  if (lk.Stat()) {//FI
	  }
	  else {//FF
	    AddDipole(new FF_Dipole(lk,lij,li,lj));
	    AddDipole(new FF_Dipole(lij,lk,li,lj));
	  }
	}
	  }
	}
      }
    }
    legs.push_back(lij);
  }
  // output dipoles
  if (msg_LevelIsDebugging()) {
    DEBUG_FUNC("");
    msg_Debugging()<<"Legs:\n";
    for (size_t i=0;i<legs.size();++i)
      msg_Debugging()<<legs[i]<<std::endl;
    msg_Debugging()<<"Dipoles:\n";
    for (size_t i=0;i<m_dipoles.size();++i)
	msg_Debugging()<<*m_dipoles[i]<<std::endl;
  }
  return true;
}


bool Extra_Emission_Generator::InitEWDipoles(Process_Base *const proc)
{
  std::vector<Cluster_Leg> legs;
  for (size_t i(0);i<proc->Flavours().size();++i) {
    Cluster_Leg lij(NULL,Vec4D(),proc->Flavours()[i]);
    if (i<m_nin) {
      lij.SetFlav(lij.Flav().Bar());
      lij.SetStat(1);
    }
    lij.SetId(1<<i);
    for (size_t j(0);j<legs.size();++j) {
      const Cluster_Leg &lk(legs[j]);
      if (lij.Flav().Electromagnetic() && lij.Flav()==lk.Flav().Bar()) {
	if (lij.Stat()) {
	  if (lk.Stat()) {//II
	  }
	  else {//IF
	  }
	}
	else {
	  if (lk.Stat()) {//FI
	  }
	  else {//FF
// 	    AddDipole(new FF_Dipole(lk,lij));
// 	    AddDipole(new FF_Dipole(lij,lk));
	  }
	}
      }
    }
    legs.push_back(lij);
  }
  // output dipoles
  if (msg_LevelIsDebugging()) {
    DEBUG_FUNC("");
    msg_Debugging()<<"Legs:\n";
    for (size_t i=0;i<legs.size();++i)
      msg_Debugging()<<legs[i]<<std::endl;
    msg_Debugging()<<"Dipoles:\n";
    for (size_t i=0;i<m_dipoles.size();++i)
	msg_Debugging()<<*m_dipoles[i]<<std::endl;
  }
  return true;
}

void Extra_Emission_Generator::SetAMin(const double &amin)
{
  for (size_t i=0;i<m_dipoles.size();++i)
      m_dipoles[i]->SetAMin(amin);
}

bool Extra_Emission_Generator::GeneratePoint(Cluster_Amplitude *const ampl)
{
  DEBUG_FUNC("");
  double rns[4];
  for (size_t i(0);i<4;++i) rns[i]=ran.Get();
  msg_Debugging()<<"in EEG: ";
  msg_Debugging()<<"#y = "<<rns[1]<<", #z = "<<rns[2]
                 <<", #phi = "<<rns[3]<<"\n";
  // choose dipole
  double asum(0.0);
  for (size_t i(0);i<m_dipoles.size();++i)
    asum+=m_dipoles[i]->Alpha();
  CS_Dipole *active(NULL);
  double disc(rns[0]*asum), psum(0.0);
  for (size_t i(0);i<m_dipoles.size();++i)
    if ((psum+=m_dipoles[i]->Alpha())>=disc) {
      active=m_dipoles[i];
      break;
    }
  if (active==NULL) THROW(fatal_error,"Internal error");
  msg_Debugging()<<"selected "<<active->Id()<<"\n";
  // copy generate extra emission
  Cluster_Amplitude *next(ampl->InitNext());
  next->CopyFrom(ampl);
  if (m_isqcd){
    next->CreateLeg(Vec4D(),Flavour(kf_jet),
		  ColorID(),1<<next->Legs().size());
  }else{
    next->CreateLeg(Vec4D(),Flavour(kf_photon),
		  ColorID(),1<<next->Legs().size());
  }
  next->Legs().back()->SetStat(1);
  next->SetIdNew(next->Legs().back()->Id());
  msg_Debugging()<<*ampl<<"\n";
  active->GeneratePoint(next,&rns[1]);
  msg_Debugging()<<*next<<"\n";
  return 1;
}

bool Extra_Emission_Generator::GenerateWeight(const Cluster_Amplitude *ampl)
{
  DEBUG_FUNC("");
  msg_Debugging()<<*ampl<<"\n";
  m_weight=0.0;
  double asum(0.0);
  for (size_t i(0);i<m_dipoles.size();++i) {
    CS_Dipole *cdip(m_dipoles[i]);
    msg_Debugging()<<"add "<<cdip->Id()<<"{\n";
    double wgt(cdip->GenerateWeight(ampl));
    msg_Debugging()<<"} -> w = "<<wgt<<", a = "<<cdip->Alpha()<<"\n";
    if (wgt==0.0) {
      m_weight=0.0;
      for (size_t j(0);j<m_dipoles.size();++j)
	m_dipoles[j]->GetVegas()->SetCheck(0);
      return false;
    }
    m_weight+=cdip->Alpha()/wgt;
    asum+=cdip->Alpha();
  }
  if (IsBad(asum/m_weight))
    msg_Error()<<METHOD<<"(): Bad weight "<<asum
	       <<" / "<<m_weight<<"."<<std::endl;
  m_weight=asum/m_weight;
  msg_Debugging()<<"m_weight = "<<m_weight<<"\n";
  if (!(m_weight>0.)) msg_Out()<<m_weight<<std::endl;
  return true;
}

void Extra_Emission_Generator::AddPoint(const double &value)
{ 
  for (size_t i(0);i<m_dipoles.size();++i) {
    m_dipoles[i]->AddPoint(value);
    msg_Debugging()<<m_dipoles[i]->Id()<<": <w> = "
		   <<m_dipoles[i]->Mean()<<" +- "
		   <<m_dipoles[i]->Sigma()<<", max = "
		   <<m_dipoles[i]->Max()<<"\n";
  }
}

void Extra_Emission_Generator::Optimize()  
{
  msg_Tracking()<<"Optimize EEG for '"
		<<p_int->Process()->Name()<<"' ("<<m_opt<<") {\n";
  size_t off(0);
  {
    msg_Indent();
    double csum(0.0), wmean(0.0), nc(0.0);
    for (size_t i(0);i<m_dipoles.size();++i) {
      CS_Dipole *v(m_dipoles[i]);
      if (v->Alpha()<=0.0) ++off;
      else {
	if (m_opt==1 && m_omode&2) v->GetVegas()->Optimize();
	if (m_omode&1) {
	  v->SetOldAlpha(v->Alpha());
	  v->SetAlpha(v->Alpha()*sqrt(dabs(v->Mean())));
	  csum+=v->Alpha();
	  wmean+=dabs(v->Mean());
	  ++nc;
	}
      }
    }
    wmean/=nc;
    csum/=nc;
    msg_Tracking()<<std::string(116,'-')<<"\n";
    for (size_t i(0);i<m_dipoles.size();++i) {
      CS_Dipole *v(m_dipoles[i]);
      if (v->Alpha()>0.0) {
	if (nc>0.0) v->SetAlpha(v->Alpha()/csum);
	if (v->Alpha()<s_alphamin) v->SetAlpha(0.0);
	double dev(int((v->Alpha()/v->OldAlpha()-1.0)*10000)/100.0);
	double re(int(v->Sigma()/v->Mean()*10000)/100.0);
	if (v->N()<2) re=100.0;
	msg_Tracking()<<std::left<<std::setw(14)<<v->Id()
		      <<": n = "<<std::right<<std::setw(5)
		      <<v->N()<<"  w' = "<<std::setw(15)
		      <<(nc>0.0?v->Mean()/wmean:v->Mean())
		      <<" +- "<<std::setw(6)<<re
		      <<" %  =>  a = "<<std::setw(15)<<v->OldAlpha()
		      <<" -> "<<std::setw(15)<<v->Alpha();
	if (v->Alpha()<s_alphamin) {
	  msg_Tracking()<<std::left<<" (      off )\n";
	}
	else {
	  msg_Tracking()<<" ( "<<std::setw(6)
			<<std::right<<dev<<std::left<<" % )\n";
	}
	v->Reset();
      }
    }
    msg_Tracking()<<std::string(116,'-')<<"\n";
  }
  msg_Tracking()<<"}";
  if (off) msg_Tracking()<<" "<<off<<" channels off";
  msg_Tracking()<<"\n";
  if (m_opt>1) --m_opt;
} 

void Extra_Emission_Generator::EndOptimize()  
{
  for (size_t i(0);i<m_dipoles.size();++i)
    m_dipoles[i]->GetVegas()->EndOptimize();
  m_opt=0;
} 

void Extra_Emission_Generator::WriteOut(std::string pid)
{ 
  MakeDir(pid,false);
  pid+="/CS";
  std::vector<std::vector<std::string> > 
    pvds(m_dipoles.size(),std::vector<std::string>(7));
  for (size_t i(0);i<m_dipoles.size();++i) {
    CS_Dipole *v(m_dipoles[i]);
    v->GetVegas()->WriteOut(pid);
    pvds[i][0]=v->Id();
    pvds[i][1]=ToString(v->Alpha(),12);
    pvds[i][2]=ToString(v->OldAlpha(),12);
    pvds[i][3]=ToString(v->N(),12);
    pvds[i][4]=ToString(v->Sum(),12);
    pvds[i][5]=ToString(v->Sum2(),12);
    pvds[i][6]=ToString(v->Max(),12);
  }
  pvds.push_back(std::vector<std::string>(1,ToString(m_opt)));
  Data_Writer writer;
  writer.SetOutputPath(pid);
  writer.SetOutputFile("_"+p_int->Process()->Name()+"_EEG_PV");
  writer.MatrixToFile(pvds);
}

void Extra_Emission_Generator::ReadIn(std::string pid)
{
  pid+="/CS";
  Data_Reader reader;
  reader.SetAddCommandLine(false);
  reader.SetInputPath(pid);
  reader.SetInputFile("_"+p_int->Process()->Name()+"_EEG_PV");
  std::vector<std::vector<std::string> > pvds;
  reader.MatrixFromFile(pvds);
  if (m_dipoles.size()>pvds.size()-1)
    THROW(fatal_error,"Corrupted input file");
  for (size_t i(0);i<m_dipoles.size();++i) {
    CS_Dipole *v(m_dipoles[i]);
    if (v->Id()!=pvds[i][0])
      THROW(fatal_error,"Corrupted input file");
    v->GetVegas()->ReadIn(pid);
    v->SetAlpha(ToType<double>(pvds[i][1],12));
    v->SetOldAlpha(ToType<double>(pvds[i][2],12));
    v->SetN(ToType<double>(pvds[i][3],12));
    v->SetSum(ToType<double>(pvds[i][4],12));
    v->SetSum2(ToType<double>(pvds[i][5],12));
    v->SetMax(ToType<double>(pvds[i][6],12));
  }
  if (pvds.back().size()!=1)
    THROW(fatal_error,"Corrupted input file");
  m_opt=ToType<int>(pvds.back().front());
}

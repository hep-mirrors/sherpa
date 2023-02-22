#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/My_MPI.H"
#ifdef USING__HDF5
#ifdef USING__MPI

#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PDF/Main/PDF_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Phys/Variations.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;
using namespace HighFive;

namespace SHERPA {

  class Output_HDF5 : public Output_Base {

    File *p_file;
    std::map<std::string,DataSet> m_dss;
    std::vector<std::vector<double> > m_ecache, m_pcache;
    double m_xs, m_xserr, m_max, m_trials;
    std::string m_basename;
    int m_setcol, m_ncache, m_events, m_unweight, m_offset, m_nmax;
    size_t m_nweights, m_neprops, m_npprops, m_clabel;
    Matrix_Element_Handler *p_me;
    Variations *p_vars;

  public:

    Output_HDF5(const Output_Arguments &args):
      Output_Base("HDF5"), m_xs(1.0), m_xserr(1.0), m_max(1.0),
      m_trials(0.0), m_events(0), m_offset(0), m_nmax(0)
    {
      m_basename=args.m_outpath+"/"+args.m_outfile;
      m_setcol=args.p_reader->GetValue<int>("HDF5_SET_COLORS",1);
      m_ncache=args.p_reader->GetValue<int>("HDF5_CACHE_SIZE",10000);
      m_unweight=args.p_reader->GetValue<int>("HDF5_UNWEIGHT",0);
      m_clabel=args.p_reader->GetValue<int>("HDF5_START_COLOR_LABEL",500);
      m_ncache=std::min(m_ncache,(int)rpa->gen.NumberOfEvents());
      p_me=args.p_init->GetMatrixElementHandler();
      p_vars=args.p_init->GetVariations();
      std::vector<std::string> params;
      args.p_reader->VectorFromFile(params,"HDF5_MPIIO_PARAMS");
      MPI_Info info;
      MPI_Info_create(&info);
      for (size_t i(0);i+1<params.size();i+=2) {
	msg_Info()<<METHOD<<"(): Add MPIIO parameters '"
		  <<params[i]<<"' -> '"<<params[i+1]<<"'\n";
	MPI_Info_set(info,params[i].c_str(),params[i+1].c_str());
      }
      p_file = new File(m_basename+".hdf5",
			File::ReadWrite|File::Create|File::Truncate,
			MPIOFileDriver(MPI_COMM_WORLD,MPI_INFO_NULL));
    }

    ~Output_HDF5()
    {
      m_dss.clear();
      delete p_file;
    }

    void SetXS(const double& xs, const double& xserr)
    {
      m_xs = xs;
      m_xserr = xserr;
      m_max =1.;
    }

    bool SetSumSqrColors(Cluster_Amplitude *const ampl)
    {
      Process_Base *proc(ampl->Proc<Process_Base>());
      SP(Color_Integrator) colint(proc->Integrator()->ColorIntegrator());
      colint->GenerateOrders();
      const Idx_Matrix &orders(colint->Orders());
      std::vector<double> psum(orders.size());
      double csum(0.0);
      for (int i(0);i<orders.size();++i) {
	int last(0), first(-1), counter(m_clabel);
	for (size_t j(0);j<orders[i].size();++j) {
	  Cluster_Leg *cl(ampl->Leg((size_t)orders[i][j]));
	  if (cl->Flav().StrongCharge()==3) {
	    last=++counter;
	    cl->SetCol(ColorID(last,0));
	  }
	  else if (cl->Flav().StrongCharge()==-3) {
	    cl->SetCol(ColorID(0,last));
	    last=0;
	  }
	  else if (cl->Flav().StrongCharge()==8) {
	    int nlast(last);
	    last=++counter;
	    cl->SetCol(ColorID(last,nlast));
	    if (nlast==0) first=j;
	  }
	  else {
	    cl->SetCol(ColorID(0,0));
	  }
	}
	if (first>=0) {
	  Cluster_Leg *fl(ampl->Leg((size_t)orders[i][first]));
	  fl->SetCol(ColorID(fl->Col().m_i,last));
	}
	msg_Debugging()<<"odering "<<orders[i]<<", first = "<<first<<"\n";
	msg_Debugging()<<*ampl<<"\n";
	if ((m_setcol&2) && first) continue;
	int valid(true);
	for (size_t j(0);j<ampl->Legs().size();++j) {
	  Cluster_Leg *fl(ampl->Leg(j));
	  if (fl->Flav().Strong() &&
	      fl->Col().m_i==fl->Col().m_j) valid=false;
	}
	if (!valid) continue;
	csum+=psum[i]=dabs(proc->Differential(*ampl,1|2|4));
	msg_Debugging()<<"sc: csum = "<<psum[i]<<"\n";
      }
      if (csum==0.0) return false;
      double disc(csum*ran->Get()), sum(0.0);
      for (size_t i(0);i<orders.size();++i)
	if ((sum+=psum[i])>=disc) {
	  msg_Debugging()<<"selected ordering "<<i<<" -> "<<orders[i]<<"\n";
	  int last(0), first(-1), counter(0);
	  for (size_t j(0);j<orders[i].size();++j) {
	    Cluster_Leg *cl(ampl->Leg((size_t)orders[i][j]));
	    if (cl->Flav().StrongCharge()==3) {
	      last=++counter;
	      cl->SetCol(ColorID(last,0));
	    }
	    else if (cl->Flav().StrongCharge()==-3) {
	      cl->SetCol(ColorID(0,last));
	      last=0;
	    }
	    else if (cl->Flav().StrongCharge()==8) {
	      int nlast(last);
	      last=++counter;
	      cl->SetCol(ColorID(last,nlast));
	      if (nlast==0) first=j;
	    }
	    else {
	      cl->SetCol(ColorID(0,0));
	    }
	  }
	  if (first>=0) {
	    Cluster_Leg *fl(ampl->Leg((size_t)orders[i][first]));
	    fl->SetCol(ColorID(fl->Col().m_i,last));
	  }
	  return true;
	}
      THROW(fatal_error,"Internal error");
    }

    void SetColors(Blob_List* blobs)
    {
      Blob *sp(blobs->FindFirst(btp::Signal_Process));
      Process_Base *proc((*sp)["Process"]->Get<Process_Base*>());
      SP(Color_Integrator) ci(proc->Integrator()->ColorIntegrator());
      Cluster_Amplitude *ampl(Cluster_Amplitude::New());
      ampl->SetNIn(sp->NInP());
      ampl->SetProc(proc);
      ampl->SetMuF2((*sp)["Factorisation_Scale"]->Get<double>());
      ampl->SetMuR2((*sp)["Renormalization_Scale"]->Get<double>());
      for (size_t i(0);i<sp->NInP();++i) {
	Particle *p(sp->InParticle(i));
	ColorID col(0,0);
	if (ci!=NULL) col=ColorID(ci->I()[i],ci->J()[i]);
	ampl->CreateLeg(-p->Momentum(),p->Flav().Bar(),col);
      }
      for (size_t i(0);i<sp->NOutP();++i) {
	Particle *p(sp->OutParticle(i));
	ColorID col(0,0);
	if (ci!=NULL) col=ColorID(ci->I()[sp->NInP()+i],ci->J()[sp->NInP()+i]);
	ampl->CreateLeg(p->Momentum(),p->Flav(),col);
      }
      msg_Debugging()<<"before color setting "<<*ampl<<"\n";
      while (!SetSumSqrColors(ampl)) {
	msg_Debugging()<<"color setting failed. generate new point\n";
	while (!ci->GeneratePoint());
	const PHASIC::Int_Vector &ni(ci->I()), &nj(ci->J());
	for (size_t i(0);i<ampl->Legs().size();++i)
	  ampl->Leg(i)->SetCol(ColorID(ni[i],nj[i]));
      }
      msg_Debugging()<<"after color setting "<<*ampl<<"\n";
      for (size_t i(0);i<sp->NInP();++i) {
	Particle *p(sp->InParticle(i));
	p->SetFlow(1,ampl->Leg(i)->Col().m_j);
	p->SetFlow(2,ampl->Leg(i)->Col().m_i);
      }
      for (size_t i(0);i<sp->NOutP();++i) {
	Particle *p(sp->OutParticle(i));
	p->SetFlow(1,ampl->Leg(sp->NInP()+i)->Col().m_i);
	p->SetFlow(2,ampl->Leg(sp->NInP()+i)->Col().m_j);
      }      
      ampl->Delete();
    }
    
    void Header()
    {
      std::vector<int> versionno = {2, 0, 0};
      m_dss["version"]=p_file->createDataSet<int>
	("version",DataSpace::From(versionno));
      m_dss["version"].write(versionno);
      std::vector<double> idata;
      std::vector<std::string> inames;
      inames.push_back("beamA");
      idata.push_back((long int)rpa->gen.Beam1());
      inames.push_back("beamB");
      idata.push_back((long int)rpa->gen.Beam2());
      inames.push_back("energyA");
      idata.push_back(rpa->gen.PBeam(0)[0]);
      inames.push_back("energyB");
      idata.push_back(rpa->gen.PBeam(1)[0]);
      inames.push_back("PDFgroupA");
      idata.push_back(0);
      inames.push_back("PDFgroupB");
      idata.push_back(0);
      inames.push_back("PDFsetA");
      idata.push_back(rpa->gen.PDF(0)?rpa->gen.PDF(0)->LHEFNumber():-1);
      inames.push_back("PDFsetB");
      idata.push_back(rpa->gen.PDF(1)?rpa->gen.PDF(1)->LHEFNumber():-1);
      inames.push_back("weightingStrategy");
      idata.push_back(ToType<int>(rpa->gen.Variable("EVENT_GENERATION_MODE"))==0?1:3);
      const PHASIC::Process_Vector procs(p_me->AllProcesses());
      inames.push_back("numProcesses");
      idata.push_back(procs.size());
      DataSetCreateProps props;
      m_dss["init"]=p_file->createDataSet<double>
	("init",DataSpace::From(idata));
      m_dss["init"].write(idata);
      m_dss["init"].createAttribute<std::string>
	("properties",DataSpace::From(inames)).write(inames);
      std::vector<std::vector<double> > pdata
	(procs.size(),std::vector<double>(6,0));
      std::vector<std::string> pnames(6);
      pnames[0]="procId";
      pnames[1]="npLO";
      pnames[2]="npNLO";
      pnames[3]="xSection";
      pnames[4]="error";
      pnames[5]="unitWeight";
      for (size_t i(0);i<pdata.size();++i) {
	pdata[i][0]=i+1;
	pdata[i][1]=pdata[i][2]=-1;
	if (procs[i]->Info().Has(nlo_type::real))
	  pdata[i][2]=procs[i]->NIn()+procs[i]->NOut()-1;
	else pdata[i][1]=procs[i]->NIn()+procs[i]->NOut();
	if (procs[i]->Get<MCatNLO_Process>()) {
	  pdata[i][3]=(*procs[i])[0]->Integrator()->TotalXS()+
	    (*procs[i])[1]->Integrator()->TotalXS();
	  pdata[i][4]=sqrt(sqr((*procs[i])[0]->Integrator()->TotalError())+
		      sqr((*procs[i])[1]->Integrator()->TotalError()));
	}
	else {
	  pdata[i][3]=procs[i]->Integrator()->TotalXS();
	  pdata[i][4]=procs[i]->Integrator()->TotalError();
	}
	pdata[i][3]*=rpa->Picobarn();
	pdata[i][4]*=rpa->Picobarn();
	pdata[i][5]=procs[i]->Integrator()->Max();
	m_nmax=std::max(m_nmax,int(procs[i]->NIn()+procs[i]->NOut()));
      }
      m_dss["procInfo"]=p_file->createDataSet<double>
	("procInfo",DataSpace::From(pdata));
      m_dss["procInfo"].write(pdata);
      m_dss["procInfo"].createAttribute<std::string>
	("properties",DataSpace::From(pnames)).write(pnames);
      Initialize(m_nmax);
    }

    void ChangeFile()
    {
    }

    void Footer()
    {
      Write(1);
    }

    void Initialize(const size_t &nup)
    {
      size_t size(MPI::COMM_WORLD.Get_size());
      std::vector<size_t> min(1,size);
      min.front()*=m_unweight?1:rpa->gen.NumberOfEvents();
      std::vector<size_t> max(min);
      DataSetCreateProps props;
      if (m_unweight) max.front()=DataSpace::UNLIMITED;
      std::vector<std::string> wnames(1,"NOMINAL");
      if (p_vars) {
	const Variations::Parameters_Vector *params
	  (p_vars->GetParametersVector());
	for (size_t i(0);i<params->size();++i)
	  wnames.push_back((*params)[i]->m_name);
      }
      m_nweights=wnames.size();
      // LHEF event information
      std::vector<std::string> enames((m_neprops=9)+m_nweights);
      enames[0]="pid";
      enames[1]="nparticles";
      enames[2]="start";
      enames[3]="trials";
      enames[4]="scale";
      enames[5]="fscale";
      enames[6]="rscale";
      enames[7]="aqed";
      enames[8]="aqcd";
      for (size_t i(0);i<wnames.size();++i)
	enames[m_neprops+i]=wnames[i];
      min.push_back(enames.size());
      max.push_back(enames.size());
      if (m_unweight) {
	props=DataSetCreateProps();
	props.add(Chunking({m_unweight*size,m_neprops+m_nweights}));
      }
      m_dss["events"]=p_file->createDataSet<double>
	("events",DataSpace(min,max),props);
      m_ecache.reserve(m_ncache);
      for (size_t i(0);i<m_ecache.size();++i)
	m_ecache[i].reserve(m_neprops+m_nweights);
      m_dss["events"].createAttribute<std::string>
	("events",DataSpace::From(enames)).write(enames);
      // LHEF particle information
      std::vector<std::string> pnames(m_npprops=13);
      pnames[0]="id";
      pnames[1]="status";
      pnames[2]="mother1";
      pnames[3]="mother2";
      pnames[4]="color1";
      pnames[5]="color2";
      pnames[6]="px";
      pnames[7]="py";
      pnames[8]="pz";
      pnames[9]="e";
      pnames[10]="m";
      pnames[11]="lifetime";
      pnames[12]="spin";
      min.front()*=nup;
      if (!m_unweight) max.front()*=nup;
      min.back()=pnames.size();
      max.back()=pnames.size();
      if (m_unweight) {
	props=DataSetCreateProps();
	props.add(Chunking({m_unweight*size*nup,m_npprops}));
      }
      m_dss["particles"]=p_file->createDataSet<double>
	("particles",DataSpace(min,max),props);
      m_pcache.reserve(m_ncache);
      for (size_t i(0);i<m_pcache.size();++i)
	m_pcache[i].reserve(m_npprops);
      m_dss["particles"].createAttribute<std::string>
	("properties",DataSpace::From(pnames)).write(pnames);
    }

    void Write(const int &mode=0)
    {
      if (mode!=1 && m_events<m_ncache) return;
      if (m_trials>0) {
	std::vector<double> eprops(m_neprops+m_nweights,-1.);
	eprops[0]=1;
	eprops[1]=m_nmax;
	eprops[2]=m_offset*m_nmax;
	for (size_t i(0);i<m_nweights;++i) eprops[m_neprops+i]=0.;
	eprops[3]=m_trials;
	m_ecache.push_back(eprops);
	m_trials=0.0;
	for (size_t i(0);i<m_nmax;++i)
	  m_pcache.push_back(std::vector<double>(m_npprops,0.));
      }
      m_events=0;
      size_t ncache(m_ecache.size());
      std::size_t size(MPI::COMM_WORLD.Get_size());
      std::size_t crank(MPI::COMM_WORLD.Get_rank());
      std::vector<int> ncaches(size,ncache);
      mpi->Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&ncaches[0],1,MPI::INT);
      size_t sumcache(0), rank(0);
      for (size_t i(0);i<size;++i) sumcache+=ncaches[i];
      for (size_t i(0);i<crank;++i) rank+=ncaches[i];
      if (sumcache==0) return;
      for (size_t i(0);i<m_ecache.size();++i) m_ecache[i][2]+=rank*m_nmax;
      if (m_unweight) {
	m_dss["events"].resize({m_offset+sumcache,m_neprops+m_nweights});
	m_dss["particles"].resize({(m_offset+sumcache)*m_nmax,m_npprops});
      }
      rank+=m_offset;
      m_offset+=sumcache;
      m_dss["events"].select({rank,0},{ncache,m_ecache.front().size()}).write(m_ecache);
      m_ecache.clear();
      m_dss["particles"].select({rank*m_nmax,0},{ncache*m_nmax,m_pcache.front().size()}).write(m_pcache);
      m_pcache.clear();
    }

    void Output(Blob_List* blobs)
    {
      DEBUG_FUNC(m_basename);
      ++m_events;
      double wratio(1.0);
      const auto weight(blobs->Weight());
      Blob *sp(blobs->FindFirst(btp::Signal_Process));
      const auto trials((*sp)["Trials"]->Get<double>());
      Process_Base *proc((*sp)["Process"]->Get<Process_Base*>());
      if (m_unweight) {
	double weight((*sp)["Weight"]->Get<double>());
	double max((*sp)["Max"]->Get<double>());
	const auto disc = max * ran->Get();
	if (std::abs(weight) < disc) {
	  m_trials+=trials;
	  Write();
	  return;
	}
	wratio=max/weight;
      }
      if (weight!=0.0 && m_setcol) SetColors(blobs);
      size_t nup(sp->NInP()+sp->NOutP());
      m_ecache.push_back(std::vector<double>(m_neprops,-1));
      m_ecache.back().push_back(weight*wratio);
      if (p_vars) {
	Variation_Weights wgts
	  ((*sp)["Variation_Weights"]->Get<Variation_Weights>());
	for (size_t i(0);i<wgts.NumberOfParameters();++i)
	  m_ecache.back().push_back(wgts.GetVariationWeightAt(i)*wratio);
      }
      m_ecache.back()[0]=0;
      if (proc) {
	while (proc->Parent()!=proc) proc=proc->Parent();
	const PHASIC::Process_Vector procs(p_me->AllProcesses());
	for (size_t i(0);i<procs.size();++i)
	  if (procs[i]==proc) m_ecache.back()[0]=i+1;
      }
      m_ecache.back()[1]=nup;
      m_ecache.back()[3]=m_trials+trials;
      m_trials=0.0;
      double mur2=(*sp)["Renormalization_Scale"]->Get<double>();
      double muf2=(*sp)["Factorisation_Scale"]->Get<double>();
      double muq2=muf2?muf2:mur2;
      if ((*sp)["Resummation_Scale"])
	muq2=(*sp)["Resummation_Scale"]->Get<double>();
      m_ecache.back()[4]=sqrt(muq2);
      m_ecache.back()[5]=sqrt(muf2);
      m_ecache.back()[6]=sqrt(mur2);
      m_ecache.back()[8]=mur2?(*MODEL::as)(mur2):-1.0;
      for (int i=0;i<sp->NInP();++i) {
	Vec4D p(sp->InParticle(i)->Momentum());
	m_pcache.push_back(std::vector<double>(13,0));
	m_pcache.back()[9]=p[0];
	m_pcache.back()[6]=p[1];
	m_pcache.back()[7]=p[2];
	m_pcache.back()[8]=p[3];
	m_pcache.back()[10]=sp->InParticle(i)->FinalMass();
	m_pcache.back()[0]=(long int)sp->InParticle(i)->Flav();
	m_pcache.back()[12]=sp->InParticle(i)->Flav().IntSpin();
	m_pcache.back()[4]=sp->InParticle(i)->GetFlow(1);
	m_pcache.back()[5]=sp->InParticle(i)->GetFlow(2);
	m_pcache.back()[1]=-1;
	m_pcache.back()[2]=0;
	m_pcache.back()[3]=0;
	m_pcache.back()[11]=0.;
      }
      for (int k=0;k<sp->NOutP();++k) {
	int i(k+sp->NInP());
	Vec4D p(sp->OutParticle(k)->Momentum());
	m_pcache.push_back(std::vector<double>(13,0));
	m_pcache.back()[9]=p[0];
	m_pcache.back()[6]=p[1];
	m_pcache.back()[7]=p[2];
	m_pcache.back()[8]=p[3];
	m_pcache.back()[10]=sp->OutParticle(k)->FinalMass();
	m_pcache.back()[0]=(long int)sp->OutParticle(k)->Flav();
	m_pcache.back()[12]=sp->OutParticle(k)->Flav().IntSpin();
	m_pcache.back()[4]=sp->OutParticle(k)->GetFlow(1);
	m_pcache.back()[5]=sp->OutParticle(k)->GetFlow(2);
	m_pcache.back()[1]=1;
	m_pcache.back()[2]=1;
	m_pcache.back()[3]=2;
	m_pcache.back()[11]=0.;
      }
      size_t rank(m_ecache.size()-1);
      m_ecache.back()[2]=(m_offset+rank)*m_nmax;
      for (size_t i(nup);i<m_nmax;++i)
	m_pcache.push_back(std::vector<double>(13,0));
      Write();
    }

  };// end of class Output_HDF5

}// end of namespace SHERPA

DECLARE_GETTER(Output_HDF5,"HDF5",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_HDF5>::
operator()(const Output_Arguments &args) const
{
  return new Output_HDF5(args);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_HDF5>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"HDF5 output (version 2)";
}

#endif
#endif

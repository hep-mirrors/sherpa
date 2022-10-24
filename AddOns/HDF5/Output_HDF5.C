#ifdef XXXXXXXXXXXXXX
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
    std::map<std::string,std::vector<int> > m_icache;
    std::map<std::string,std::vector<double> > m_dcache;
    std::map<std::string,std::vector<size_t> > m_scache;
    std::map<std::string,int> m_nup;
    std::vector<std::vector<double> > m_wcache;
    double m_xs, m_xserr, m_max, m_trials;
    std::string m_basename;
    int m_setcol, m_ncache, m_events, m_unweight, m_offset, m_nmax;
    size_t m_nweights, m_clabel;
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
      Group init = p_file->createGroup("init");
      std::vector<int> beamid = { int((long int)rpa->gen.Beam1()) } ;
      m_dss["beamA"]=init.createDataSet<int>("beamA",DataSpace::From(beamid));
      m_dss["beamA"].write(beamid);
      beamid[0]=(long int)rpa->gen.Beam2();
      m_dss["beamB"]=init.createDataSet<int>("beamB",DataSpace::From(beamid));
      m_dss["beamB"].write(beamid);
      std::vector<double> ebeam = { rpa->gen.PBeam(0)[0] };
      m_dss["energyA"]=init.createDataSet<double>("energyA",DataSpace::From(ebeam));
      m_dss["energyA"].write(ebeam);
      ebeam[0]=rpa->gen.PBeam(1)[0];
      m_dss["energyB"]=init.createDataSet<double>("energyB",DataSpace::From(ebeam));
      m_dss["energyB"].write(ebeam);
      std::vector<int> pdfgroup = { 0 };
      m_dss["PDFgroupA"]=init.createDataSet<int>("PDFgroupA",DataSpace::From(pdfgroup));
      m_dss["PDFgroupA"].write(pdfgroup);
      m_dss["PDFgroupB"]=init.createDataSet<int>("PDFgroupB",DataSpace::From(pdfgroup));
      m_dss["PDFgroupB"].write(pdfgroup);
      std::vector<int> pdfset = { rpa->gen.PDF(0)?rpa->gen.PDF(0)->LHEFNumber():-1 };
      m_dss["PDFsetA"]=init.createDataSet<int>("PDFsetA",DataSpace::From(pdfset));
      m_dss["PDFsetA"].write(pdfset);
      pdfset[0]=rpa->gen.PDF(1)?rpa->gen.PDF(1)->LHEFNumber():-1;
      m_dss["PDFsetB"]=init.createDataSet<int>("PDFsetB",DataSpace::From(pdfset));
      m_dss["PDFsetB"].write(pdfset);
      std::vector<int> weighting = { ToType<int>(rpa->gen.Variable("EVENT_GENERATION_MODE"))==0?1:3 };
      m_dss["weightingStrategy"]=init.createDataSet<int>("weightingStrategy",DataSpace::From(weighting));
      m_dss["weightingStrategy"].write(weighting);
      const PHASIC::Process_Vector procs(p_me->AllProcesses());
      std::vector<int> nprocs(1,procs.size());
      m_dss["numProcesses"]=init.createDataSet<int>("numProcesses",DataSpace::From(nprocs));
      m_dss["numProcesses"].write(nprocs);
      std::vector<int> versionno = {1, 0, 0};
      m_dss["version"]=init.createDataSet<int>("version",DataSpace::From(versionno));
      m_dss["version"].write(versionno);
      Group process = p_file->createGroup("procInfo");
      std::vector<int> procid(procs.size()), nplo(procs.size()), npnlo(procs.size());
      std::vector<double> xsec(procs.size()), err(procs.size()), max(procs.size());
      for (size_t i(0);i<procid.size();++i) {
	procid[i]=i+1;
	npnlo[i]=nplo[i]=-1;
	if (procs[i]->Info().Has(nlo_type::real))
	  npnlo[i]=procs[i]->NIn()+procs[i]->NOut()-1;
	else nplo[i]=procs[i]->NIn()+procs[i]->NOut();
	if (procs[i]->Get<MCatNLO_Process>()) {
	  xsec[i]=(*procs[i])[0]->Integrator()->TotalXS()+
	    (*procs[i])[1]->Integrator()->TotalXS();
	  err[i]=sqrt(sqr((*procs[i])[0]->Integrator()->TotalError())+
		      sqr((*procs[i])[1]->Integrator()->TotalError()));
	}
	else {
	  xsec[i]=procs[i]->Integrator()->TotalXS();
	  err[i]=procs[i]->Integrator()->TotalError();
	}
	xsec[i]*=rpa->Picobarn();
	err[i]*=rpa->Picobarn();
	max[i]=procs[i]->Integrator()->Max();
	m_nmax=std::max(m_nmax,int(procs[i]->NIn()+procs[i]->NOut()));
      }
      m_dss["procId"]=process.createDataSet<int>("procId",DataSpace::From(procid));
      m_dss["procId"].write(procid);
      m_dss["npLO"]=process.createDataSet<int>("npLO",DataSpace::From(nplo));
      m_dss["npLO"].write(nplo);
      m_dss["npNLO"]=process.createDataSet<int>("npNLO",DataSpace::From(npnlo));
      m_dss["npNLO"].write(npnlo);
      m_dss["xSection"]=process.createDataSet<double>("xSection",DataSpace::From(xsec));
      m_dss["xSection"].write(xsec);
      m_dss["error"]=process.createDataSet<double>("error",DataSpace::From(err));
      m_dss["error"].write(err);
      m_dss["unitWeight"]=process.createDataSet<double>("unitWeight",DataSpace::From(max));
      m_dss["unitWeight"].write(max);
      Initialize(m_nmax);
    }

    void ChangeFile()
    {
    }

    void Footer()
    {
      Write(1);
    }

    void AddSizeTDataSet(Group &group,const std::string &name,
			 const std::vector<size_t> &min,
			 const std::vector<size_t> &max,
			 const DataSetCreateProps &props)
    {
      m_dss[name]=group.createDataSet<size_t>(name,DataSpace(min,max),props);
      m_scache[name].reserve(m_ncache);
    }

    void AddIntDataSet(Group &group,const std::string &name,
		       std::vector<size_t> &min,std::vector<size_t> &max,
		       const DataSetCreateProps &props,const int nup=1)
    {
      min.front()*=nup;
      if (max.front()!=DataSpace::UNLIMITED) max.front()*=nup;
      m_dss[name]=group.createDataSet<int>(name,DataSpace(min,max),props);
      min.front()/=nup;
      if (max.front()!=DataSpace::UNLIMITED) max.front()/=nup;
      m_icache[name].reserve(m_ncache);
      m_nup[name]=nup;
    }

    void AddDoubleDataSet(Group &group,const std::string &name,
			  std::vector<size_t> &min,std::vector<size_t> &max,
			  const DataSetCreateProps &props,const int nup=1)
    {
      min.front()*=nup;
      if (max.front()!=DataSpace::UNLIMITED) max.front()*=nup;
      m_dss[name]=group.createDataSet<double>(name,DataSpace(min,max),props);
      min.front()/=nup;
      if (max.front()!=DataSpace::UNLIMITED) max.front()/=nup;
      m_dcache[name].reserve(m_ncache);
      m_nup[name]=nup;
    }

    void Initialize(const size_t &nup)
    {
      size_t size(MPI::COMM_WORLD.Get_size());
      std::vector<size_t> min(1,size);
      min.front()*=m_unweight?1:rpa->gen.NumberOfEvents();
      std::vector<size_t> max(min);
      DataSetCreateProps props;
      if (m_unweight) {
	max.front()=DataSpace::UNLIMITED;
	props.add(Chunking(std::vector<hsize_t>(1,m_unweight*size)));
      }
      // LHEF event information
      Group event = p_file->createGroup("event");
      AddIntDataSet(event,"pid",min,max,props);
      AddIntDataSet(event,"nparticles",min,max,props);
      AddSizeTDataSet(event,"start",min,max,props);
      AddDoubleDataSet(event,"trials",min,max,props);
      AddDoubleDataSet(event,"scale",min,max,props);
      AddDoubleDataSet(event,"fscale",min,max,props);
      AddDoubleDataSet(event,"rscale",min,max,props);
      AddDoubleDataSet(event,"aqed",min,max,props);
      AddDoubleDataSet(event,"aqcd",min,max,props);
      // LHEF particle information
      Group particle = p_file->createGroup("particle");
      AddIntDataSet(particle,"id",min,max,props,nup);
      AddIntDataSet(particle,"status",min,max,props,nup);
      AddIntDataSet(particle,"mother1",min,max,props,nup);
      AddIntDataSet(particle,"mother2",min,max,props,nup);
      AddIntDataSet(particle,"color1",min,max,props,nup);
      AddIntDataSet(particle,"color2",min,max,props,nup);
      AddDoubleDataSet(particle,"px",min,max,props,nup);
      AddDoubleDataSet(particle,"py",min,max,props,nup);
      AddDoubleDataSet(particle,"pz",min,max,props,nup);
      AddDoubleDataSet(particle,"e",min,max,props,nup);
      AddDoubleDataSet(particle,"m",min,max,props,nup);
      AddDoubleDataSet(particle,"lifetime",min,max,props,nup);
      AddDoubleDataSet(particle,"spin",min,max,props,nup);
      std::vector<std::string> wnames(1,"NOMINAL");
      if (p_vars) {
	const Variations::Parameters_Vector *params
	  (p_vars->GetParametersVector());
	for (size_t i(0);i<params->size();++i)
	  wnames.push_back((*params)[i]->m_name);
      }
      m_nweights=wnames.size();
      min.push_back(m_nweights);
      max.push_back(m_nweights);
      if (m_unweight) {
	props=DataSetCreateProps();
	props.add(Chunking({m_unweight*size,m_nweights}));
      }
      m_dss["weight"]=event.createDataSet<double>
	("weight",DataSpace(min,max),props);
      m_wcache.reserve(m_ncache);
      for (size_t i(0);i<m_wcache.size();++i)
	m_wcache[i].reserve(m_nweights);
      m_dss["weight"].createAttribute<std::string>
	("weight_names",DataSpace::From(wnames)).write(wnames);
    }

    void Write(const int &mode=0)
    {
      if (mode!=1 && m_events<m_ncache) return;
      if (m_events==0) {
	if (m_trials)
	  std::cout<<METHOD<<"(rank="
		   <<MPI::COMM_WORLD.Get_rank()
		   <<": Error! # trials = "<<m_trials<<".\n";
	return;
      }
      if (m_scache["start"].empty()) {
	m_icache["pid"].push_back(1);
	m_icache["nparticles"].push_back(m_nmax);
	m_dcache["scale"].push_back(-1.);
	m_dcache["fscale"].push_back(-1.);
	m_dcache["rscale"].push_back(-1.);
	m_dcache["aqed"].push_back(-1.);
	m_dcache["aqcd"].push_back(-1.);
	m_scache["start"].push_back(m_offset*m_nmax);
	for (size_t i(0);i<m_nmax;++i) {
	  m_dcache["e"].push_back(0.);
	  m_dcache["px"].push_back(0.);
	  m_dcache["py"].push_back(0.);
	  m_dcache["pz"].push_back(0.);
	  m_dcache["m"].push_back(0.);
	  m_icache["id"].push_back(0);
	  m_dcache["spin"].push_back(0);
	  m_icache["color1"].push_back(0);
	  m_icache["color2"].push_back(0);
	  m_icache["status"].push_back(0);
	  m_icache["mother1"].push_back(0);
	  m_icache["mother2"].push_back(0);
	  m_dcache["lifetime"].push_back(0.);
	}
	m_wcache.push_back(std::vector<double>(m_nweights,0.));
	m_dcache["trials"].push_back(m_trials);
	m_trials=0.0;
      }
      m_events=0;
      size_t ncache(m_scache["start"].size());
      std::size_t size(MPI::COMM_WORLD.Get_size());
      std::size_t crank(MPI::COMM_WORLD.Get_rank());
      std::vector<int> ncaches(size,ncache);
      mpi->MPIComm()->Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,&ncaches[0],1,MPI::INT);
      size_t sumcache(0), rank(0);
      for (size_t i(0);i<size;++i) sumcache+=ncaches[i];
      for (size_t i(0);i<crank;++i) rank+=ncaches[i];
      if (sumcache==0) return;
      std::vector<size_t> &start(m_scache["start"]);
      for (size_t i(0);i<start.size();++i) start[i]+=rank*m_nup["id"];
      if (m_unweight) {
	for (std::map<std::string,std::vector<size_t> >::iterator
	       it=m_scache.begin();it!=m_scache.end();++it)
	  m_dss[it->first].resize(std::vector<size_t>(1,m_offset+sumcache));
	for (std::map<std::string,std::vector<int> >::iterator
	       it=m_icache.begin();it!=m_icache.end();++it)
	  m_dss[it->first].resize(std::vector<size_t>(1,(m_offset+sumcache)*m_nup[it->first]));
	for (std::map<std::string,std::vector<double> >::iterator
	       it=m_dcache.begin();it!=m_dcache.end();++it)
	  m_dss[it->first].resize(std::vector<size_t>(1,(m_offset+sumcache)*m_nup[it->first]));
	m_dss["weight"].resize({m_offset+sumcache,m_nweights});
      }
      rank+=m_offset;
      m_offset+=sumcache;
      for (std::map<std::string,std::vector<size_t> >::iterator
	     it=m_scache.begin();it!=m_scache.end();++it) {
	m_dss[it->first].select({rank},{ncache}).write(it->second);
	it->second.clear();
      }
      for (std::map<std::string,std::vector<int> >::iterator
	     it=m_icache.begin();it!=m_icache.end();++it) {
	int nup(it->second.size()/ncache);
	m_dss[it->first].select({rank*nup},{ncache*nup}).write(it->second);
	it->second.clear();
      }
      for (std::map<std::string,std::vector<double> >::iterator
	     it=m_dcache.begin();it!=m_dcache.end();++it) {
	int nup(it->second.size()/ncache);
	m_dss[it->first].select({rank*nup},{ncache*nup}).write(it->second);
	it->second.clear();
      }
      m_dss["weight"].select({rank,0},{ncache,m_wcache.front().size()}).write(m_wcache);
      m_wcache.clear();
    }

    void Output(Blob_List* blobs)
    {
      DEBUG_FUNC(m_basename);
      ++m_events;
      double wratio(1.0);
      const auto weight(blobs->Weight());
      Blob *sp(blobs->FindFirst(btp::Signal_Process));
      const auto trials((*sp)["Trials"]->Get<double>());
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
      m_icache["pid"].push_back(1);
      m_icache["nparticles"].push_back(nup);
      m_dcache["trials"].push_back(m_trials+trials);
      m_trials=0.0;
      double mur2=(*sp)["Renormalization_Scale"]->Get<double>();
      double muf2=(*sp)["Factorisation_Scale"]->Get<double>();
      double muq2=muf2?muf2:mur2;
      if ((*sp)["Resummation_Scale"])
	muq2=(*sp)["Resummation_Scale"]->Get<double>();
      m_dcache["scale"].push_back(sqrt(muq2));
      m_dcache["fscale"].push_back(sqrt(muf2));
      m_dcache["rscale"].push_back(sqrt(mur2));
      Poincare cms(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
      std::vector<double> weights(1,weight*wratio);
      if (p_vars) {
	Variation_Weights wgts
	  ((*sp)["Variation_Weights"]->Get<Variation_Weights>());
	for (size_t i(0);i<wgts.NumberOfParameters();++i)
	  weights.push_back(wgts.GetVariationWeightAt(i)*wratio);
      }
      m_wcache.push_back(weights);
      m_dcache["aqed"].push_back(-1.);
      m_dcache["aqcd"].push_back(mur2?(*MODEL::as)(mur2):-1.0);
      for (int i=0;i<sp->NInP();++i) {
	Vec4D p(sp->InParticle(i)->Momentum());
	m_dcache["e"].push_back(p[0]);
	m_dcache["px"].push_back(p[1]);
	m_dcache["py"].push_back(p[2]);
	m_dcache["pz"].push_back(p[3]);
	m_dcache["m"].push_back(sp->InParticle(i)->FinalMass());
	m_icache["id"].push_back((long int)sp->InParticle(i)->Flav());
	m_dcache["spin"].push_back(sp->InParticle(i)->Flav().IntSpin());
	m_icache["color1"].push_back(sp->InParticle(i)->GetFlow(1));
	m_icache["color2"].push_back(sp->InParticle(i)->GetFlow(2));
	m_icache["status"].push_back(-1);
	m_icache["mother1"].push_back(0);
	m_icache["mother2"].push_back(0);
	m_dcache["lifetime"].push_back(0.);
      }
      for (int k=0;k<sp->NOutP();++k) {
	int i(k+sp->NInP());
	Vec4D p(sp->OutParticle(k)->Momentum());
	m_dcache["e"].push_back(p[0]);
	m_dcache["px"].push_back(p[1]);
	m_dcache["py"].push_back(p[2]);
	m_dcache["pz"].push_back(p[3]);
	m_dcache["m"].push_back(sp->OutParticle(k)->FinalMass());
	m_icache["id"].push_back((long int)sp->OutParticle(k)->Flav());
	m_dcache["spin"].push_back(sp->OutParticle(k)->Flav().IntSpin());
	m_icache["color1"].push_back(sp->OutParticle(k)->GetFlow(1));
	m_icache["color2"].push_back(sp->OutParticle(k)->GetFlow(2));
	m_icache["status"].push_back(1);
	m_icache["mother1"].push_back(1);
	m_icache["mother2"].push_back(2);
	m_dcache["lifetime"].push_back(0.);
      }
      size_t rank(m_scache["start"].size());
      m_scache["start"].push_back((m_offset+rank)*m_nmax);
      for (size_t i(nup);i<m_nmax;++i) {
	m_dcache["e"].push_back(0.);
	m_dcache["px"].push_back(0.);
	m_dcache["py"].push_back(0.);
	m_dcache["pz"].push_back(0.);
	m_dcache["m"].push_back(0.);
	m_icache["id"].push_back(0);
	m_dcache["spin"].push_back(0);
	m_icache["color1"].push_back(0);
	m_icache["color2"].push_back(0);
	m_icache["status"].push_back(0);
	m_icache["mother1"].push_back(0);
	m_icache["mother2"].push_back(0);
	m_dcache["lifetime"].push_back(0.);
      }
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
  str<<"HDF5 output";
}

#endif
#endif
#endif

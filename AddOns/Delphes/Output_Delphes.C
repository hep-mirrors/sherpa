#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__DELPHES
#ifndef USING__ROOT
#error Delphes cannot be used without Root. \
  Please reconfigure with --enable-root
#endif
#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Tools/HepMC2_Interface.H"
#include "Utilities/ExRootAnalysis/interface/ExRootTreeWriter.h"
#include "Utilities/ExRootAnalysis/interface/BlockClasses.h"

namespace HepMC { class IO_GenEvent; class GenCrossSection; class GenEvent; }

namespace DELPHES {

  class Output_Delphes: public SHERPA::Output_Base {
  private:

    SHERPA::HepMC2_Interface m_hepmc2;

    HepMC::IO_GenEvent *p_iogenevent;
    HepMC::GenCrossSection *p_xs;
    HepMC::GenEvent *p_event;

    ExRootTreeWriter *p_treeWriter;
    ExRootTreeBranch *p_branchGenEvent;
    ExRootTreeBranch *p_branchGenParticle;

    std::string m_basename, m_ext;

    int m_mode;

    void AnalyseParticles(ExRootTreeBranch *branch,const HepMC::GenEvent& evt);
    void AnalyseEvent(ExRootTreeBranch *branch,const HepMC::GenEvent& evt,
		      const Long64_t eventNumber,const double weight);
    void ReadStats();
    void getStatsFromTuple(int &mo1, int &mo2, int &da1, int &da2,
			   int &status, int &pid, int j) const;

    std::vector<HepMC::GenParticle*> index_to_particle;
    std::map<HepMC::GenParticle*,int> particle_to_index;

    inline int find_in_map(HepMC::GenParticle *p) const
    {
      std::map<HepMC::GenParticle*,int>::const_iterator iter = particle_to_index.find(p);
      return (iter == particle_to_index.end()) ? 0 : iter->second;
    }

  public:

    Output_Delphes(const SHERPA::Output_Arguments &args,int mode);

    ~Output_Delphes();

    void SetXS(const double& xs, const double& xserr);
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void ChangeFile();

  };// end of class Output_Delphes

  class Output_Delphes_GenEvent: public Output_Delphes {};

}// end of namespace DELPHES

#include "HepMC/GenEvent.h"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

#ifdef USING__HEPMC2__DEFS
#include "HepMC/HepMCDefs.h"
#ifdef HEPMC_HAS_CROSS_SECTION
#include "HepMC/GenCrossSection.h"
#endif
#endif

#include "Utilities/ExRootAnalysis/interface/ExRootTreeWriter.h"
#include "Utilities/ExRootAnalysis/interface/BlockClasses.h"
#include "Utilities/ExRootAnalysis/interface/ExRootTreeBranch.h"

using namespace DELPHES;
using namespace SHERPA;
using namespace ATOOLS;

Output_Delphes::Output_Delphes(const Output_Arguments &args,int mode):
  Output_Base(mode?"Delphes":"DelphesS"), p_event(NULL), m_mode(mode)
{
  m_basename=args.m_outpath+"/"+args.m_outfile;
  m_ext=m_mode?"_ge.root":".root";
  int precision       = args.p_reader->GetValue<int>("OUTPUT_PRECISION",12);
#ifdef HEPMC_HAS_CROSS_SECTION
  p_xs = new HepMC::GenCrossSection();
#endif
  if (m_mode==0) p_event = new HepMC::GenEvent();
  p_treeWriter = new ExRootTreeWriter(m_basename+m_ext,"GEN");
  p_branchGenEvent = p_treeWriter->NewBranch
    ("Event", TRootLHEFEvent::Class());
  p_branchGenParticle = p_treeWriter->NewBranch
    ("Particle", TRootC::GenParticle::Class());
}

Output_Delphes::~Output_Delphes()
{
  p_treeWriter->Write();
  delete p_treeWriter;
  p_treeWriter=NULL;
  if (m_mode==0) delete p_event;
}

void Output_Delphes::SetXS(const double& xs, const double& xserr)
{
#ifdef HEPMC_HAS_CROSS_SECTION
  p_xs->set_cross_section(xs, xserr);
#endif
}

void Output_Delphes::Output(Blob_List* blobs, const double weight)
{
#ifdef USING__HEPMC2__IOGENEVENT
  if (m_mode==0) {
    p_event->clear();
    m_hepmc2.Sherpa2ShortHepMC(blobs,*p_event,weight);
  }
  else {
    m_hepmc2.Sherpa2HepMC(blobs,weight);
    p_event = m_hepmc2.GenEvent();
  }
#ifdef HEPMC_HAS_CROSS_SECTION
  p_event->set_cross_section(*p_xs);
#endif
  p_treeWriter->Clear();
  AnalyseEvent(p_branchGenEvent,*p_event,
	       p_event->event_number()+1,weight);
  AnalyseParticles(p_branchGenParticle,*p_event);      
  p_treeWriter->Fill();
#endif
}

void Output_Delphes::ChangeFile()
{
}

void Output_Delphes::AnalyseEvent
(ExRootTreeBranch *branch,const HepMC::GenEvent& evt,
 const Long64_t eventNumber,const double weight)
{
  TRootLHEFEvent *element;
  element = static_cast<TRootLHEFEvent*>(branch->NewEntry());
  element->Number = eventNumber;
  element->Nparticles = evt.particles_size();
  element->ProcessID = evt.signal_process_id();
  element->Weight = weight;
  element->ScalePDF = evt.event_scale();
  element->CouplingQED = evt.alphaQED();
  element->CouplingQCD = evt.alphaQCD();
}

void  Output_Delphes::ReadStats()
{
 unsigned int particle_counter=0;
 index_to_particle.clear(); 
 particle_to_index.clear();
 index_to_particle.reserve(p_event->particles_size()); 
 index_to_particle[0] = 0; 
 HepMC::GenEvent::vertex_const_iterator v;
 for (v = p_event->vertices_begin(); v != p_event->vertices_end(); ++v ) {
   // making a list of incoming particles of the vertices
   // so that the mother indices in HEPEVT can be filled properly
   HepMC::GenVertex::particles_out_const_iterator p1;
   for (p1 = (*v)->particles_in_const_begin();p1 != (*v)->particles_in_const_end(); ++p1 ) {
     ++particle_counter;
     //particle_counter can be very large for heavy ions
     if(particle_counter >= index_to_particle.size() ) {
       //make it large enough to hold up to this index
       index_to_particle.resize(particle_counter+1);
     }                 
     index_to_particle[particle_counter] = *p1;
     particle_to_index[*p1] = particle_counter;
   }   
   // daughters are entered only if they aren't a mother of
   // another vertex
   HepMC::GenVertex::particles_out_const_iterator p2;
   for (p2 = (*v)->particles_out_const_begin();p2 != (*v)->particles_out_const_end(); ++p2) {
     if (!(*p2)->end_vertex()) {
       ++particle_counter;
       //particle_counter can be very large for heavy ions
       if(particle_counter  >= index_to_particle.size() ) {        
	 //make it large enough to hold up to this index
	 index_to_particle.resize(particle_counter+1);
       }                   
       index_to_particle[particle_counter] = *p2;
       particle_to_index[*p2] = particle_counter;
     }
   }
 } 
}

void Output_Delphes::getStatsFromTuple(int &mo1, int &mo2, int &da1, int &da2, int &status, int &pid, int j) const
{
  status =  index_to_particle[j]->status();
  pid = index_to_particle[j]->pdg_id();
  if ( index_to_particle[j]->production_vertex() ) {
    int num_mothers = index_to_particle[j]->production_vertex()->particles_in_size();
    if (num_mothers ==0) { 
      mo1 = 0; 
      mo2 = 0; 
    }
    else {
      int first_mother = find_in_map(*(index_to_particle[j]->production_vertex()->particles_in_const_begin()));
      int last_mother = first_mother + num_mothers - 1;
      if ( first_mother == 0 ) last_mother = 0;
      mo1=first_mother;
      mo2=last_mother;
    } // if num_mothers !=0
  } 
  else // no data on production_vertex
    {
      mo1 =0;
      mo2 =0;
    }   
  if (index_to_particle[j]->end_vertex()) {
    // make sure first and last daughter are indeed the first and last
    int first_daughter = find_in_map(*(index_to_particle[j]->end_vertex()->particles_begin(HepMC::children)));
    int last_daughter = 0;
    HepMC::GenVertex::particle_iterator ic;
    for (ic = index_to_particle[j]->end_vertex()->particles_begin(HepMC::children);ic != index_to_particle[j]->end_vertex()->particles_end(HepMC::children); ++ic) {
      int current_daughter = find_in_map(*ic) ; 
      if (current_daughter < first_daughter)
	first_daughter = current_daughter;
      if (current_daughter > last_daughter)
	last_daughter = current_daughter;
    }
    if (first_daughter== 0) last_daughter = 0;
    da1=first_daughter;
    da2=last_daughter;
  } 
  else {
    da1=0;
    da2=0;
  }
} 

void Output_Delphes::AnalyseParticles(ExRootTreeBranch *branch, const HepMC::GenEvent& evt)
{
  TRootC::GenParticle *element;
  TLorentzVector momentum;
  Double_t signPz;
  ReadStats();
  for(int n=1; n<=evt.particles_size(); n++) {
    int mo1, mo2, da1, da2, status, pid;
    getStatsFromTuple(mo1,mo2,da1,da2,status,pid,n);
    element = static_cast<TRootC::GenParticle*>(branch->NewEntry());
    element->PID = pid;
    element->Status = status;
    element->M1 = mo1 - 1; // added -1 as the numbering in the tree starts from 0
    element->M2 = mo2 - 1;
    element->D1 = da1 - 1;
    element->D2 = da2 - 1;
    element->E = index_to_particle[n]->momentum().e();
    element->Px = index_to_particle[n]->momentum().px();
    element->Py = index_to_particle[n]->momentum().py();
    element->Pz = index_to_particle[n]->momentum().pz();
    element->PT = sqrt(pow(element->Px,2)+pow(element->Py,2));
    momentum.SetPxPyPzE(element->Px, element->Py, element->Pz, element->E);
    signPz = (element->Pz >= 0.0) ? 1.0 : -1.0;
    element->Eta = element->PT < 1e-6 ? signPz*999.9 : momentum.Eta();
    element->Phi = index_to_particle[n]->momentum().phi();
    HepMC::GenVertex* vrtI = (index_to_particle[n])->production_vertex();
    HepMC::GenVertex::particles_in_const_iterator partI;
    if(vrtI) {
      element->T = vrtI->position().t();
      element->X = vrtI->position().x();
      element->Y = vrtI->position().y();
      element->Z = vrtI->position().z();
    }
    else {
      element->T = 0.;
      element->X = 0.;
      element->Y = 0.;
      element->Z = 0.;
    }  
  }
}

DECLARE_GETTER(Output_Delphes,"Delphes_Short",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_Delphes>::
operator()(const Output_Arguments &args) const
{
  return new Output_Delphes(args,0);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_Delphes>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Delphes short output";
}

DECLARE_GETTER(Output_Delphes_GenEvent,"Delphes_GenEvent",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_Delphes_GenEvent>::
operator()(const Output_Arguments &args) const
{
  return new Output_Delphes(args,1);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_Delphes_GenEvent>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Delphes long output";
}

#endif

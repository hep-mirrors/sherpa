#include <cassert>
#include "AHADIC++/Main/Ahadic.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "AHADIC++/Tools/Dipole.H"
#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "AHADIC++/Tools/Hadron_Init.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "AHADIC++/Formation/Cluster_Formation_Handler.H"
#include "AHADIC++/Decays/Cluster_Decay_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Ahadic::Ahadic(string path, string file)  :
  m_fullinfo(false), m_maxtrials(3), m_clulist()
{
  
  Data_Reader dr(" ",";","!","=");
  dr.AddComment("#");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(path);
  dr.SetInputFile(file);
  m_anadir = dr.GetValue<string>("FRAGMANADIR", "Fragmentation_Analysis");
  hadpars =  new Hadronisation_Parameters();
  hadpars->Init(path,file);

  p_cformhandler = new Cluster_Formation_Handler(&m_clulist,hadpars->AnaOn());
  p_cdechandler  = new Cluster_Decay_Handler(&m_clulist,hadpars->AnaOn());



  msg_Tracking()<<"Initialisation of Ahadic complete."<<endl;
}


Ahadic::~Ahadic() 
{
  writeHistos();
  CleanUp();
  if (p_cdechandler)  { delete p_cdechandler;  p_cdechandler=NULL;  }
  if (p_cformhandler) { delete p_cformhandler; p_cformhandler=NULL; }
  delete hadpars;
}

// One write function for all histos in AHADIC analysis
void Ahadic::writeHistos() {
    std::map < std::string, std::map< std::string, ATOOLS::Histogram * > > megaMap;
    //if (p_cdechandler)
        //megaMap[string("ClusterDecay")]     = p_cdechandler->getHistos();
    if (p_cformhandler) {
      megaMap[string("ClusterFormation")] = p_cformhandler->GetHistos();
      megaMap[string("GluonDecayer")]     = p_cformhandler->GetGDHistos();
      megaMap[string("SoftClusters")]     = p_cformhandler->GetSCHistos();
    }
   
    /// Write out to disk, create subdirectories if necessary
    string outdir = m_anadir;
    if (!megaMap.empty()) ATOOLS::MakeDir(outdir, 493);
    // Outer Loop over analyses
    for (map<std::string, std::map< std::string, ATOOLS::Histogram * > >::iterator oit=megaMap.begin();
         oit!=megaMap.end();oit++) {
        msg_Tracking() << "Processing " << oit->first << endl;
        if (!oit->second.empty()) ATOOLS::MakeDir(string(m_anadir + '/'+ oit->first), 493);
        // Inner loop over histograms
        for (map<string,Histogram *>::iterator hit=oit->second.begin();
             hit!=oit->second.end();hit++) {
            #ifdef USING__MPI
            hit->second->MPISync();
            #endif
            if (hit->second->Integral() > 0 ) hit->second->Scale(1./hit->second->Integral());
            hit->second->Output(m_anadir + '/'+ oit->first + "/" + hit->first+std::string(".dat"));
        }
    }
    msg_Info() << "Output written to " << m_anadir << endl;
    //megaMap.clear();
}

Return_Value::code Ahadic::Hadronize(Blob_List * blobs)
{
  static std::string mname(METHOD);
  Return_Value::IncCall(mname);

  // Filter BlobList for a reduced list of Blobs that need hadronisation
  std::vector<Blob*> v_hadronize;
  for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();) {
    if ((*blit)->Has(blob_status::needs_hadronization) &&
	(*blit)->Type()==btp::Fragmentation) {
        v_hadronize.push_back(*blit);
    }
    blit++;
  }

  // Loop over blobs that need hadronisation
  Blob * blob(NULL);
  Return_Value::code result;
  bool moveon(false);
  for (unsigned int ib=0;ib<v_hadronize.size();++ib) {

      blob   = v_hadronize[ib];
      moveon = false;
      Reset(); // Ahadic Reset --- calls Cluster::ResetClusterNumber

      /// m_maxtrials is set in the ctor --- the retry argument, i, is
      /// ONLY used for debug output
      for (short int i=0;i<m_maxtrials;i++) {
	try {
            // Calls Formation and Decay, second arg only for debug output
            result = Hadronize(blob ,i); // result contains return value 
	} 
        catch (Return_Value::code ret) { // Seems to catch a throw from Cluster_Former::ConstructClusters
	  msg_Error()<<"ERROR in "<<METHOD<<" : \n"
		     <<"   Hadronization for blob "<<(*blob)<<"\n"
		     <<"("<<blob<<"; "<<blob->NInP()<<" -> "
		     <<blob->NOutP()<<") "
		     <<"did not work out,\n"
		     <<"   will trigger retrying the event.\n";
	  CleanUp(blob);
	  return Return_Value::Retry_Event;
	}

        /// Decision on how to proceed based on return value of Hadronize        
	switch (result) {
            case Return_Value::Success : 
              moveon = true;
              break;
            case Return_Value::Retry_Event : 
              {
              blobs->ColorConservation();
              msg_Tracking()<<"ERROR in "<<METHOD<<" :\n"
                            <<"   Hadronization for blob "
                            <<"("<<blob<<"; "<<blob->NInP()<<" -> "
                            <<blob->NOutP()<<") "
                            <<"did not work out,"<<std::endl
                            <<"   will trigger retrying the event.\n";
              CleanUp(blob);
              if (rpa->gen.Beam1().IsLepton() ||
                  rpa->gen.Beam2().IsLepton()) {
                msg_Tracking()<<METHOD<<": "
                              <<"Non-hh collision. Request new event instead.\n";
                return Return_Value::New_Event;
              }
              return result;
              }
            case Return_Value::Retry_Method :
              {
              msg_Tracking()<<"Warning in "<<METHOD<<" : "<<std::endl
                            <<"   Hadronization for blob "
                            <<"("<<blob<<"; "<<blob->NInP()<<" -> "
                            <<blob->NOutP()<<") "
                            <<"did not work out properly in the "
                            <<(i+1)<<"th attempt,"<<std::endl
                            <<"   retry it "<<m_maxtrials<<" times."<<std::endl;
              Return_Value::IncRetryMethod(mname);
              CleanUp(blob);
              }
              break;
            case Return_Value::Nothing :
            default:
              msg_Tracking()<<"Warning in "<<METHOD<<":"<<std::endl
                            <<"   Calling Hadronization for Blob("<<blob<<") "
                            <<"yields "<<int(result)<<"."<<std::endl
                            <<"   Continue and hope for the best."<<std::endl;
              moveon = true;
              break;
	} // End switch
	
        if (moveon) break;
      } // End trials

      CleanUp(); // NOT Ahadic::Cleanup(Blob)!!!
     
      if (moveon && SanityCheck(blob, sqr(rpa->gen.Ecms()))) {
	blob->SetStatus(blob_status::needs_hadrondecays);
	blob->SetType(btp::Fragmentation);
      }
      else {
	CleanUp(blob);
	return Return_Value::Retry_Event;
      }
    }
  return Return_Value::Success;
}  



Return_Value::code Ahadic::Hadronize(Blob * blob, int retry) {
  assert(m_clulist.empty());
  blob->SetType(btp::Cluster_Formation);
  blob->SetTypeSpec("AHADIC-1.0");


  switch (p_cformhandler->FormClusters(blob)) {
  case -1 : 
    msg_Tracking()<<"ERROR in "<<METHOD<<" (formation):"<<std::endl
		  <<"   Will retry event."<<std::endl;
    p_cformhandler->Reset();
    return Return_Value::Retry_Event;
  case  0 :
    msg_Tracking()<<"ERROR in "<<METHOD<<" (formation):"<<std::endl
		  <<"   Will retry method."<<std::endl;
    p_cformhandler->Reset();
    return Return_Value::Retry_Method;
  case 1 :
    if (retry>0) 
      msg_Tracking()<<"   Passed cluster form now ("<<retry<<"th trial).\n";
    break;
  }
  
  switch (p_cdechandler->DecayClusters(blob)) {
  case -1 : 
    msg_Tracking()<<"ERROR in "<<METHOD<<" (decay):"<<std::endl
		  <<"   Will retry event."<<std::endl;
    return Return_Value::Retry_Event;
  case  0 :
    msg_Tracking()<<"ERROR in "<<METHOD<<" (decay):"<<std::endl
		  <<"   Will retry method."<<std::endl;
    return Return_Value::Retry_Method;
  case  1 :
    if (retry>0) 
      msg_Tracking()<<"   Passed cluster decays now ("<<retry<<"th trial).\n";
    break;
  }


  
  if (!m_clulist.empty()) {
    msg_Error()<<"Error in "<<METHOD<<" (deletion):"<<std::endl
	       <<"   "<<m_clulist.size()<<" clusters undeleted:"<<std::endl
	       <<m_clulist<<std::endl;
  }
  assert(m_clulist.empty());

  if (blob->CheckMomentumConservation().Abs2()>1.e-6) {
    msg_Error()<<METHOD<<"(): Momentum imbalance: "
	       <<blob->CheckMomentumConservation()<<"\n";
    msg_Debugging()<<(*blob)<<"\n";
    return Return_Value::Retry_Event;
  }
  return Return_Value::Success;
}

/// Reset clusters ???
void Ahadic::Reset() {
  if (Cluster::RemainingClusters()) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   "<<Cluster::RemainingClusters()<<" clusters undeleted.\n";
  }
  assert(!Cluster::RemainingClusters());
  Cluster::ResetClusterNumber();
  control::s_AHAparticles=0;
}

bool Ahadic::SanityCheck(Blob * blob,double norm2) {
  Vec4D checkmom(blob->CheckMomentumConservation());
  if (dabs(checkmom.Abs2())/norm2>1.e-12 ||
      (norm2<0. && norm2>0.) ||
      Cluster::RemainingClusters()!=0 ||
      control::s_AHAparticles!=blob->NOutP()
      /* || control::s_AHAprotoparticles!=0*/) {
    msg_Error()<<"ERROR in "<<METHOD<<" : "<<endl
	       <<"   Momentum/particle-blob number violation for "
	       <<(Cluster::RemainingClusters())
	       <<" remaining clusters (norm2 = "<<norm2<<")."<<endl
	       <<"   Protoparticles = "<<control::s_AHAprotoparticles
	       <<"/ parts = "<<control::s_AHAparticles<<" vs. "<<blob->NOutP()
	       <<"   : "<<checkmom<<" ("<<sqrt(Max(0.,checkmom.Abs2()))<<")\n"
	       <<(*blob)<<endl;
    return false;
  }
  msg_Tracking()<<"Passed "<<METHOD<<" with "
		<<"protoparticles = "<<control::s_AHAprotoparticles
		<<"/ parts = "<<control::s_AHAparticles<<" vs. "<<blob->NOutP()
		<<"   : "<<blob->CheckMomentumConservation()<<endl;
  msg_Debugging()<<(*blob)<<endl;
  return true;
}

void Ahadic::CleanUp(Blob * blob) {
  if (Cluster::RemainingActives()>0) {
    msg_Tracking()<<METHOD<<": "<<Cluster::RemainingActives()
		  <<" remaining Clusters found:"<<std::endl;
    //Cluster::PrintActives(msg_Tracking());
    Cluster::DeleteActives();
  }
  if (!m_clulist.empty()) m_clulist.clear();
  
  if (Dipole::RemainingActives()>0) {
    msg_Tracking()<<METHOD<<": "<<Dipole::RemainingActives()
		  <<" remaining Dipoles found:"<<std::endl;
    //Dipole::PrintActives(msg_Tracking());
    Dipole::DeleteActives();
  }

  if (Proto_Particle::RemainingActives()>0) {
    msg_Tracking()<<METHOD<<": "<<Proto_Particle::RemainingActives()
		  <<" remaining Proto_Particles found:"<<std::endl;
    //Proto_Particle::PrintActives(msg_Tracking());
    Proto_Particle::DeleteActives();
  }

  if (Proto_Particle_List::RemainingActives()>0) {
    msg_Tracking()<<METHOD<<": "<<Proto_Particle_List::RemainingActives()
		  <<" remaining Proto_Particle_Lists found:"<<std::endl;
    //Proto_Particle_List::PrintActives(msg_Tracking());
    Proto_Particle_List::DeleteActives();
  }


  if(blob) blob->DeleteOutParticles(0);
}

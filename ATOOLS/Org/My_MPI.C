#include "ATOOLS/Org/My_MPI.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

#include <stddef.h>
#include <cstring>
#include <unistd.h>

using namespace ATOOLS;

My_MPI *ATOOLS::mpi(NULL);

My_MPI::My_MPI():
  m_hassend(false), m_hasrecv(false)
{
}

My_MPI::~My_MPI()
{
#ifdef USING__MPI
  if (m_hassend) m_send.Free();
  if (m_hasrecv) m_recv.Free();
#endif  
}

void My_MPI::SetUpSendRecv(Data_Reader *const read)
{
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    std::string rname=read->GetValue<std::string>
      ("MPI_NODE_NAME",std::string(".*"));
    double starttime=rpa->gen.Timer().RealTime();
    msg_Info()<<METHOD<<"(): Analyzing MPI environment {\n";
    int rank=MPI::COMM_WORLD.Get_rank(), pid(getpid()), hlen;
    char host[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(host,&hlen);
    std::vector<std::string> hres(RegExMatch(host,rname));
    if (hres.size()) strcpy(host,hres.front().c_str());
    if (rank==0) {
      msg_Info()<<"  Rank "<<rank<<", pid "<<pid
		<<" running on "<<host<<".\n";
      char mhost[MPI_MAX_PROCESSOR_NAME];
      MPI_Get_processor_name(mhost,&hlen);
      std::vector<std::string> hres(RegExMatch(mhost,rname));
      if (hres.size()) strcpy(mhost,hres.front().c_str());
      std::vector<std::string> hosts(size,host);
      std::map<std::string,size_t> procs;
      std::vector<std::vector<int> > recv(size);
      std::map<std::string,std::vector<int> > recvmap;
      std::vector<int> mrecv(1,-1);
      procs[mhost]=1;
      for (int tag=1;tag<size;++tag) {
	MPI::COMM_WORLD.Recv(&pid,1,MPI::INT,MPI::ANY_SOURCE,tag);
	MPI::COMM_WORLD.Recv(host,MPI_MAX_PROCESSOR_NAME,
			     MPI::CHAR,MPI::ANY_SOURCE,tag);
	msg_Info()<<"  Rank "<<tag<<", pid "<<pid
		  <<" running on "<<host<<"."<<std::endl;
	hosts[tag]=host;
	recvmap[hosts[tag]].push_back(tag);
	if (procs.find(hosts[tag])!=procs.end()) {
	  if (strcmp(host,mhost)==0 && procs[mhost]==1)
	    mrecv.push_back(tag);
	  ++procs[hosts[tag]];
	}
	else {
	  mrecv.push_back(tag);
	  procs[hosts[tag]]=1;
	}
      }
      if (procs.size()==1) {
	mrecv.insert(mrecv.end(),
		     recvmap[mhost].begin(),
		     recvmap[mhost].end());
	recvmap[mhost].resize(1);
      }
      std::set<std::string> send;
      for (int tag=1;tag<size;++tag) {
	std::vector<int> recv(1,recvmap[hosts[tag]][0]);
	if (send.find(hosts[tag])==send.end()) {
	  send.insert(hosts[tag]);
	  recv=recvmap[hosts[tag]];
	  recv[0]=0;
	}
	if (recv.size()>1)
	  msg_Info()<<"  Rank "<<tag<<" send/recv "<<recv<<".\n";
	int nrecv(recv.size());
	MPI::COMM_WORLD.Send(&nrecv,1,MPI::INT,tag,size+tag);
	MPI::COMM_WORLD.Send(&recv.front(),nrecv,MPI::INT,tag,size+tag);
      }
      SetMPIRecv(mrecv);
      double diff=rpa->gen.Timer().RealTime()-starttime;
      msg_Info()<<"} -> "<<FormatTime(size_t(diff))<<" elapsed"<<std::endl;
    }
    else {
      MPI::COMM_WORLD.Send(&pid,1,MPI::INT,0,rank);
      MPI::COMM_WORLD.Send(host,MPI_MAX_PROCESSOR_NAME,MPI::CHAR,0,rank);
      int nrecv;
      MPI::COMM_WORLD.Recv(&nrecv,1,MPI::INT,0,size+rank);
      std::vector<int> recv(nrecv);
      MPI::COMM_WORLD.Recv(&recv.front(),nrecv,MPI::INT,0,size+rank);
      SetMPIRecv(recv);
    }
  }
#endif
}

void My_MPI::SetMPIRecv(std::vector<int> r)
{
#ifdef USING__MPI
  int rank=MPI::COMM_WORLD.Get_rank();
  if (rank==0) {
    m_hasrecv=true;
    m_recv=MPI::COMM_WORLD.Split(rank,rank);
    m_send=MPI::COMM_WORLD.Split(MPI_UNDEFINED,rank);
  }
  else {
    if (r[0]==0) {
      m_hassend=m_hasrecv=true;
      m_send=MPI::COMM_WORLD.Split(r[0],rank);
      m_recv=MPI::COMM_WORLD.Split(rank,rank);
    }
    else {
      m_hassend=true;
      m_recv=MPI::COMM_WORLD.Split(MPI_UNDEFINED,rank);
      m_send=MPI::COMM_WORLD.Split(r[0],rank);
    }
  }
#endif
}

bool My_MPI::HasMPISend() const
{
  return m_hassend;
}

bool My_MPI::HasMPIRecv() const
{
#ifdef USING__MPI
  if (m_hasrecv) return m_recv.Get_size()>1;
#endif
  return false;
}


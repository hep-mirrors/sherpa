#include "HD_PS_Base.H"
#include "Hadron_Decay_Channel.H"
#include "Two_Body_PSs.H"
#include "Three_Body_PSs.H"
#include "Rambo.H"
#include "Data_Reader.H"
#include "Message.H"

using namespace HADRONS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Single_Channel * HD_Channel_Selector::GetChannel(int nin,int nout,
						 const Flavour * flavs,string name)
{
  if (nin>1 || nout<2) {
    msg.Error()<<"Error in HD_Channel_Selector::GetChannel : "<<endl
	       <<"   No PS for channel ("<<nin<<" -> "<<nout<<" )"<<endl
	       <<"   Return nothing and hope for the best."<<endl;
    return NULL;
  }
  if (name==string("Isotropic")) return new Rambo(1,nout,flavs);
  if (nout==2) {
    if (name==string("Iso2")) return new Iso2Channel(flavs);
  }
  if (nout==3) {
    if (name==string("Dalitz_photon_23")) return new Dalitz(flavs,Flavour(kf::photon),2,3);
    if (name==string("Dalitz_rho_770_23"))    return new Dalitz(flavs,Flavour(kf::rho_770),2,3);
  }

  msg.Error()<<"Error in HD_Channel_Selector::GetChannel : "<<endl
	     <<"   No channel for ("<<nin<<" -> "<<nout<<" )"<<endl
	     <<"   Return nothing and hope for the best."<<endl;
  return NULL;
}


HD_PS_Base::HD_PS_Base(Hadron_Decay_Channel * hdc,vector<string> & _pst,
		       bool & mustinit) :
  Multi_Channel(_pst[2]), p_hdc(hdc),
  p_channelselector(new HD_Channel_Selector), m_file(string("")),
  m_res(-1.), m_error(0.), m_max(-1.), m_flux(1./(2.*hdc->Flavours()[0].Mass()))
{
  if (_pst.size()>2) m_file = _pst[3];
  mustinit = Construct();
  delete p_channelselector;
}


HD_PS_Base::~HD_PS_Base() {}

void HD_PS_Base::Initialise() 
{
  CalculateNormalisedWidth();
  WriteOut();
  //abort();
}

bool HD_PS_Base::Construct()
{
  if (m_file!=string("")) {
    vector<vector<string> > helpsvv;
    Data_Reader reader = Data_Reader(string("|"),string(";"),string("!"));
    reader.AddComment("#");
    reader.AddComment("//");
    reader.SetInputPath("./");
    reader.SetInputFile(m_file);
    reader.SetMatrixType(reader.MTransposed);
    reader.MatrixFromFile(helpsvv,"");
    string name;
    double weight;
    for (int i=0;i<helpsvv.size();i++) {
      if (helpsvv[i][0]==string("Channels")) {
	i++;
	while (helpsvv[i][0]!=string("}")) {
	  weight=1.;
	  if (helpsvv[i].size()>1) {
	    weight=atof(helpsvv[i][1].c_str());
	  }
	  AddChannel(helpsvv[i][0],weight);
	  i++;
	}
      }
      if (helpsvv[i][0]==string("Dalitz-Parameters")) {
	i++;
	while (helpsvv[i][0]!=string("}")) {
	  if (helpsvv[i].size()==5) {
	    vector<double> dals;
	    for (int j=0;j<5;j++) dals.push_back(atof(helpsvv[i][j].c_str()));
	    //p_hdc->GetME()->SetDalitzParameters(dals);
	  }
	  i++;
	}
      }
      if (helpsvv[i][0]==string("Result")) {
	i++;
	while (helpsvv[i][0]!=string("}")) {
	  m_res   = atof(helpsvv[i][0].c_str());
	  m_error = atof(helpsvv[i][1].c_str());
	  m_max   = atof(helpsvv[i][2].c_str());
	  i++;
	}
	return false;
      }
    }
  }
  return true;
}

void HD_PS_Base::AddChannel(string name,double weight) {
  Single_Channel * sc = p_channelselector->GetChannel(1,p_hdc->NOut(),
						      p_hdc->Flavours(),name);
  sc->SetAlpha(weight);
  Add(sc);
}

void HD_PS_Base::CalculateNormalisedWidth() {
  cout<<"HD_PS_Base::CalculateNormalisedWidth()"<<endl;
  Reset();
  long int iter = Number()*5000*int(pow(2.,int(p_hdc->NOut())-2));
  int maxopt    = Number()*int(pow(2.,2*(int(p_hdc->NOut())-2)));
  
  long int n;
  int      opt=0;
  double   value, oldvalue=0., sum=0., sum2=0., result=-1., disc;
  bool     maxincrease, simple=false;
  while(opt<maxopt || m_error>0.01) {
    maxincrease = false;
    for (n=1;n<iter+1;n++) {
      value = p_hdc->Differential();
      sum  += value;
      sum2 += ATOOLS::sqr(value);
      AddPoint(value);
      if (value>m_max) { m_max = value; maxincrease = true; }
      if (value!=0. && value==oldvalue) { simple = true; break; }
      oldvalue = value;
    }
    opt++;
    Optimize(0.01);

    if (simple) break;
    n      = opt*iter;
    result = sum/n;
    disc   = sqr(sum/n)/((sum2/n - sqr(sum/n))/(n-1));
    if (disc>0) m_error  = result/disc;
  } 
  m_res  = m_flux*sum/n;
  m_max *= m_flux;
  disc   = sqr(m_res)/((sum2*sqr(m_flux)/n - sqr(m_res))/(n-1));
  if (disc>0) m_error  = result/disc;
}


bool HD_PS_Base::WriteOut() {
  system((string("mv ")+m_file+string(" ")+m_file+string(".old")).c_str());

  ofstream to;
  to.open(m_file.c_str(),ios::out);
  to<<"Channels {"<<endl;
  for (int i=0;i<channels.size();i++) {
    if (channels[i]->Name()==string("Rambo"))
      to<<"    Isotropic"<<" "<<channels[i]->Alpha()<<";"<<endl;
    else
      to<<"    "<<channels[i]->ChID()<<" "<<channels[i]->Alpha()<<";"<<endl;
  }
  to<<"}"<<endl;
  
  char buffer[100];
  ifstream from;
  from.open((m_file+string(".old")).c_str());
  while (from.getline(buffer,100)) {
    if (buffer==string("Dalitz-Parameters {")) {
      to<<"Dalitz-Parameters {"<<endl;
      while (buffer!=string("}")) {
	from.getline(buffer,100);
	to<<buffer<<endl;
      }
      break;
    }
  }
  from.close();
  to<<"Result {"<<endl;
  to<<"   "<<m_res<<" "<<m_error<<" "<<m_max<<";"<<endl;
  to<<"}"<<endl;
  to.close();
}

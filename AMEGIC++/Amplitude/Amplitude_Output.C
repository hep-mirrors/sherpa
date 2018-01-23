#include "AMEGIC++/Amplitude/Amplitude_Output.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Amplitude_Output::Amplitude_Output(std::string pid, Topology * _top,
                                   std::string gpath, int gmode)
{
  m_graphmode=gmode;
  std::string script("/plot_graphs");
  if (!FileExists(rpa->gen.Variable("SHERPA_CPP_PATH")+script))
    Copy(rpa->gen.Variable("SHERPA_SHARE_PATH")+script,
         rpa->gen.Variable("SHERPA_CPP_PATH")+script);
  gpath+=std::string("/Amegic/");
  MakeDir(gpath);
  pid=pid.substr(pid.rfind('/')+1);
  std::string fname=gpath+pid+std::string(".tex");
  pios.open(fname.c_str());
  top = _top;
  ampl=0;
  counter=0;
  maincounter=1;
  subcounter=0;
  super_amplitude=false;
  for (int i=0; i<3;++i) captions.push_back("");
  WriteHeader(pid);
}

void Amplitude_Output::WriteHeader(const std::string &name)
{
  pios<<"\\documentclass{article} "<<endl;
  pios<<"\\usepackage{feynmp} "<<endl;
  pios<<"\\unitlength=1mm "<<endl;
  pios<<"\\newcommand{\\m}{-}"<<endl;
  pios<<"\\newcommand{\\p}{+}"<<endl;
  pios<<"\\newcommand{\\ti}{*}"<<endl;

  pios<<"\\setlength{\\textwidth}{25cm}"<<endl;
  pios<<"\\setlength{\\textheight}{25cm}"<<endl;
  pios<<"\\setlength{\\topmargin}{0cm}"<<endl;
  pios<<"\\setlength{\\headsep}{0pt}"<<endl;
  pios<<"\\setlength{\\headheight}{0pt}"<<endl;
  pios<<"\\setlength{\\oddsidemargin}{0pt}"<<endl;
  pios<<"\\setlength{\\evensidemargin}{0pt} "<<endl;

  pios<<"\\setlength{\\tabcolsep}{5mm}  "<<endl;

  pios<<"\\begin{document} "<<endl;
  pios<<"\\pagestyle{empty}"<<endl;
  pios<<"\\begin{fmffile}{"<<name<<"_fg} "<<endl;

}

string Amplitude_Output::Int2String(const int i) {
  MyStrStream str;
  str<<i;
  string o;
  str>>o;
  return o;
}

void Amplitude_Output::LegCount(Point * mo) {
  if (!mo) {
    msg_Error()<<METHOD<<"(): ERROR: no point found, continue run."<<endl;
    return;
  }

  if (mo->left==0) {
    if (mo->b==1) ++nout;
    else ++nin;
    return;
  }
  ++nmed;
  LegCount(mo->left);
  LegCount(mo->right);
  if (mo->middle)
    LegCount(mo->middle);
}

int Amplitude_Output::InclInComming(Point * mo) {
  if (mo==0) return 0;
  if (mo->b==-1 && mo->left==0)
    return 1;
  if (mo->left==0)
    return 0;

  int test = 4*InclInComming(mo->right);
  test += 2*InclInComming(mo->middle);
  test += InclInComming(mo->left);
  
  if (test==0) return 0;

  if (test==4) {
    Point * help=mo->left;
    mo->left = mo->right;
    mo->right = help;
  } 
  if (test==2) {
    Point * help=mo->left;
    mo->left = mo->middle;
    mo->middle = help;
  } 
  return 1;
}

void Amplitude_Output::WriteOut(Point * start) {
  // count orders
  std::vector<int> cpls;
  start->FindOrder(cpls);

  // make working copy
  if (ampl==0) ampl=new Point[12*2 + 1];  //  up to 10 outgoing particles
  int count_all=0;
  top->Copy(start,ampl,count_all);

  InclInComming(ampl);  

  ostream & s= pios;
  nin=1;
  nout=0;
  nmed=0;
  LegCount(ampl);

  // number, fl
  for (int i=0; i<nin; ++i) {
    ins.push_back(string("i")+Int2String(i));
  }
  for (int i=0; i<nout; ++i) {
    outs.push_back(string("o")+Int2String(i));
  }
  for (int i=0; i<nmed; ++i) {
    meds.push_back(string("v")+Int2String(i));
  }

  // the writing out:
  // start graph environment
  s<<endl;
  if (counter%3==0) {
    s<<"\\begin{tabular}{ccc}"<<endl;
  }

  // write caption with graph number
  MyStrStream str;
  if (m_graphmode==1) {
    if (super_amplitude)
      str<<maincounter<<"("<<subcounter++<<")";
    else
      str<<maincounter++;
  }
  else {
    str<<maincounter++;
  }
  str>>captions[counter%3];
  captions[counter%3]=std::string(" Graph ")+captions[counter%3];

  // add orders to caption
  captions[counter%3]+=std::string(" $\\mathcal{O}(g_s^")+ToString(cpls[0])
                       +std::string("\\,g^")+ToString(cpls[1]);
  for (size_t i(2);i<cpls.size();++i)
    captions[counter%3]+=std::string("\\,g_\text{BSM")+ToString(i-1)
                         +std::string("}^")+ToString(cpls[i]);
  captions[counter%3]+=std::string(")$");

  s<<" % Graph "<<++counter<<endl;
  s<<"\\begin{fmfgraph*}(40,40) "<<endl;

  if (m_graphmode==1) {
    // define incoming points at bottom
    s<<"  \\fmfbottom{"<<ins[0];
    for (int i=1; i<nin; ++i) s<<","<<ins[i];
    s<<"} "<<endl;
    // define outgoing points at top
    s<<"  \\fmftop{"<<outs[nout-1];
    for (int i=nout-2; i>=0; --i) s<<","<<outs[i];
    s<<"} "<<endl;
  }
  else {
    // define incoming and outgoints points in a circle
    s<<"  \\fmfsurround{"<<outs[0];
    for (int i=1; i<nout; ++i) s<<","<<outs[i];
    for (int i=0; i<nin; ++i) s<<","<<ins[i];
    s<<"} "<<endl;
  }

  // draw start line (left incoming)
  s<<"  \\fmf{";
  if (ampl->fl.IsPhoton()) s<<"photon";
  else if (ampl->fl.IsGluon()) s<<"gluon";
  else if (ampl->fl.IsVector()) s<<"boson";
  else if (ampl->fl.IsFermion()) s<<"fermion";
  else if (ampl->fl.IsScalar()) s<<"dashes";
  else s<<"dots";

  bool flip(ampl->fl.IsAnti());
  std::string begin(ins[0]);
  std::string end(meds[0]);
  if (flip) std::swap(begin,end);
  s<<",label=$"<<ampl->fl.TexName()<<"$}{"<<begin<<","<<end<<"} "<<endl;

  s<<"  \\fmfv{label="<<ampl->number<<"}{"<<ins[0]<<"} "<<endl;

  oc=0;
  ic=1;
  mc=1;

  DrawLine(meds[0],ampl->left,flip);
  DrawLine(meds[0],ampl->middle,flip);
  DrawLine(meds[0],ampl->right,flip);
  // draw dots

  // draw numbers
  

  // close graph environment
  s<<"\\end{fmfgraph*} "<<endl<<endl;
  if (counter%3==0) {
    s<<"\\\\[15pt]"<<endl;
    for (int i=0;;++i) {
      s<<captions[i];
      if (i==2) break;
      s<<" & "<<endl;
    }
    s<<"\\\\[15mm]"<<endl;

    s<<"\\end{tabular}"<<endl;
  }
  else {
    s<<"&"<<endl;
  }

}


void Amplitude_Output::DrawLine(string from, Point * d, bool flip) {
  if (d==0) return;

  ostream & s= pios;

  string to;
  if (d->left==0 && d->b==1) {
    // if (d->left==0) {
    to=outs[oc++];
    s<<"  \\fmfv{label="<<d->number<<"}{"<<to<<"} "<<endl;
  }
  else if (d->left==0) {
    to=ins[ic++];
    s<<"  \\fmfv{label="<<d->number<<"}{"<<to<<"} "<<endl;
  }
  else 
    to=meds[mc++];

  // draw line
  s<<"  \\fmf{";
  if (d->fl.IsPhoton()) s<<"photon";
  else if (d->fl.IsGluon()) s<<"gluon";
  else if (d->fl.IsVector()) s<<"boson";
  else if (d->fl.IsFermion()) s<<"fermion";
  else if (d->fl.IsScalar()) s<<"dashes";
  else s<<"dots";
  s<<",label=$"<<d->fl.TexName()<<"$";

  if (!flip)
    s<<"}{"<<from<<","<<to<<"} "<<endl;
  else
    s<<"}{"<<to<<","<<from<<"} "<<endl;

  DrawLine(to,d->left,flip);
  DrawLine(to,d->middle,flip);
  DrawLine(to,d->right,flip);
}


Amplitude_Output::~Amplitude_Output()
{
  if (counter%3!=0) {
    pios<<"\\\\[12pt]"<<endl;
    for (int i=0;;++i) {
      pios<<captions[i];
      if (i==counter%3-1) break;
      pios<<" & "<<endl;
    }
    pios<<endl;
    pios<<"\\end{tabular}"<<endl;
  }
  pios<<"\\end{fmffile} "<<endl;
  pios<<"\\end{document} "<<endl;
  pios.close();

  if (ampl) delete [] ampl;  ampl=0;
}

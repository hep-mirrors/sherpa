#include "AMEGIC++/Amplitude/Amplitude_Output.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Amplitude_Output::Amplitude_Output(std::string pid, Topology * _top)
{
  std::string fname=rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+pid+std::string(".tex");
  pios.open(fname.c_str());
  top = top;
  ampl=0;
  counter=0;
  maincounter=1;
  subcounter=0;
  super_amplitude=false;
  for (int i=0; i<3;++i) captions.push_back("");
  WriteHeader();
}

void Amplitude_Output::WriteHeader()
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
  pios<<"\\begin{fmffile}{graphs} "<<endl;

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
    msg_Error()<<"ERROR in Amplitude_Output::LegCount : no point found, continue run."<<endl;
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
  MyStrStream str;
  if (super_amplitude) 
    str<<maincounter<<"("<<subcounter++<<")";
  else 
    str<<maincounter++;
  str>>captions[counter%3];
  captions[counter%3]=std::string(" Graph ")+captions[counter%3];

  s<<" % Graph "<<++counter<<endl;
  s<<"\\begin{fmfgraph*}(40,40) "<<endl;

  // define incoming points
  s<<"  \\fmfbottom{"<<ins[0];
  for (int i=1; i<nin; ++i) s<<","<<ins[i];
  s<<"} "<<endl;
  // define outgoing points
  s<<"  \\fmftop{"<<outs[nout-1];
  for (int i=nout-2; i>=0; --i) s<<","<<outs[i];
  s<<"} "<<endl;
  
  /*
  // define incoming and outgoints points in a circle
  s<<"  \\fmfsurround{"<<outs[0];
  for (int i=1; i<nout; ++i) s<<","<<outs[i];
  for (int i=0; i<nin; ++i) s<<","<<ins[i];
  s<<"} "<<endl;
  */

  // draw start line
  s<<"  \\fmf{";
  if (ampl->fl.IsPhoton()) s<<"photon";
  else if (ampl->fl.IsGluon()) s<<"gluon";
  else if (ampl->fl.IsBoson()) s<<"dashes";
  else s<<"fermion";
  if (!(ampl->fl.IsGluon())&& !(ampl->fl.IsPhoton()))
    s<<",label=$"<<ampl->fl.TexName()<<"$}{"<<ins[0]<<","<<meds[0]<<"} "<<endl;
  else
    s<<"}{"<<ins[0]<<","<<meds[0]<<"} "<<endl;

  s<<"  \\fmfv{label="<<ampl->number<<"}{"<<ins[0]<<"} "<<endl;

  oc=0;
  ic=1;
  mc=1;
  
  DrawLine(meds[0],ampl->left);
  DrawLine(meds[0],ampl->middle);
  DrawLine(meds[0],ampl->right);
  // draw dots

  // draw numbers
  

  // close graph environment
  s<<"\\end{fmfgraph*} "<<endl<<endl;
  if (counter%3==0) {
    s<<"\\\\"<<endl;
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


void Amplitude_Output::DrawLine(string from, Point * d) {
  ostream & s= pios;
  if (d==0) return;

  int in_or_out=0;

  string to;
  if (d->left==0 && d->b==1) {
    // if (d->left==0) {
    to=outs[oc++];
    s<<"  \\fmfv{label="<<d->number<<"}{"<<to<<"} "<<endl;
    in_or_out=1;
  }
  else if (d->left==0) {
    to=ins[ic++];
    s<<"  \\fmfv{label="<<d->number<<"}{"<<to<<"} "<<endl;
    in_or_out=1;
  }
  else 
    to=meds[mc++];

  // draw line
  s<<"  \\fmf{";
  if (d->fl.IsPhoton()) s<<"photon";
  else if (d->fl.IsGluon()) s<<"gluon";
  else if (d->fl.IsBoson()) s<<"dashes";
  else if (in_or_out) s<<"fermion";
  else s<<"plain";
  if (!(d->fl.IsGluon())&& !(d->fl.IsPhoton()))
    s<<",label=$"<<d->fl.TexName()<<"$";
  if (d->b==-1)
    s<<"}{"<<to<<","<<from<<"} "<<endl;
  else 
    s<<"}{"<<from<<","<<to<<"} "<<endl;
  
  DrawLine(to,d->left);
  DrawLine(to,d->middle);
  DrawLine(to,d->right);
}


Amplitude_Output::~Amplitude_Output()
{
  if (counter%3!=0) {
    pios<<"\\\\"<<endl;
    for (int i=0;;++i) {
      pios<<captions[i];
      if (i==counter%3-1) break;
      pios<<" & "<<endl;
    }
    pios<<"\\end{tabular}"<<endl;
  }
  pios<<"\\end{fmffile} "<<endl;
  pios<<"\\end{document} "<<endl;
  pios.close();

  if (ampl) delete [] ampl;  ampl=0;
}

#include "ISR_Channel.H"
#include "Channel_Elements.H"
#include "Channel_Basics.H"
#include "Run_Parameter.H"
#include "Message.H"

#include <stdio.h>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

SimplePoleUniform::SimplePoleUniform(const double exponent,const std::string cinfo,
				     Integration_Info *info):
  m_exponent(exponent)
{
  char help[3];
  sprintf(help,"%i",int(100.*exponent));
  name=std::string("Simple_Pole_Uniform_")+std::string(help);
  m_spkey.SetInfo(std::string("Simple_Pole_")+std::string(help));
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}
  
void SimplePoleUniform::GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode) 
{
  CalculateLimits(spkey,ykey);
  m_spkey[3]=CE.MasslessPropMomenta(m_exponent,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]+=CE.DiceYUniform(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}
  
void SimplePoleUniform::GenerateWeight(const int mode) 
  {
  if (m_spkey.Weight()==UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_exponent,m_spkey[0],m_spkey[1],m_spkey[3]);
      }
    }
  if (m_ykey.Weight()==UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYUniform(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),mode);
      }
    }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
  }
  
bool SimplePoleUniform::CalculateLimits(Info_Key &spkey,Info_Key &ykey) 
  {
  m_spkey[2]=spkey[2];
  if (!m_zchannel) {
    m_spkey[0]=spkey[0];
    m_spkey[1]=spkey[3];
    }
    else {
    m_spkey[0]=spkey[3];
    m_spkey[1]=spkey[1];
    m_ykey[0]=ykey[0];
    m_ykey[1]=ykey[1];
    double logtau=.5*log(spkey[3]/spkey[2]);
    m_xkey[0]=logtau+ykey[2];
    m_xkey[2]=logtau-ykey[2];
    m_xkey[1]=0.;
    m_xkey[3]=0.;
  }
  return true;
}

SimplePoleUniform::SimplePoleUniform(double _sprimeexp, double _deltay1, double _deltay2) :
  sprimeexp(_sprimeexp)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  char help[3];
  sprintf(help,"%i",int(100.*sprimeexp));
  name     = std::string("SimplePoleUniform"+std::string(help));
  ms = rans = 0;
}


void SimplePoleUniform::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  sprime = CE.MasslessPropMomenta(sprimeexp,sprimerange[0],sprimerange[1],rans[0]);
  y      = CE.DiceYUniform(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void SimplePoleUniform::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  weight  = 1./CE.MasslessPropWeight(sprimeexp,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYUniform(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


SimplePoleCentral::SimplePoleCentral(double _sprimeexp, double _deltay1, double _deltay2) :
  sprimeexp(_sprimeexp)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  char help[3];
  sprintf(help,"%i",int(100.*sprimeexp));
  name     = std::string("SimplePoleCentral"+std::string(help));
  ms = rans = 0;
}


void SimplePoleCentral::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  sprime = CE.MasslessPropMomenta(sprimeexp,sprimerange[0],sprimerange[1],rans[0]);
  y      = CE.DiceYCentral(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void SimplePoleCentral::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  weight  = 1./CE.MasslessPropWeight(sprimeexp,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYCentral(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 



SimplePoleForward::SimplePoleForward(double _sprimeexp,double _yexp,double _deltay1, double _deltay2) : 
  sprimeexp(_sprimeexp), yexp(_yexp)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  char help[3];
  sprintf(help,"%i",int(100.*sprimeexp));
  name     = std::string("SimplePoleForward"+std::string(help));
  ms = rans = 0;
}


void SimplePoleForward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  sprime = CE.MasslessPropMomenta(sprimeexp,sprimerange[0],sprimerange[1],rans[0]);
  y = CE.DiceYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void SimplePoleForward::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  weight  = 1./CE.MasslessPropWeight(sprimeexp,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 

SimplePoleBackward::SimplePoleBackward(double _sprimeexp,double _yexp,double _deltay1, double _deltay2) :
  sprimeexp(_sprimeexp), yexp(_yexp)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  char help[3];
  sprintf(help,"%i",int(100.*sprimeexp)); 
  name     = std::string("SimplePoleBackward"+std::string(help));
  ms = rans = 0;
}

void SimplePoleBackward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  sprime = CE.MasslessPropMomenta(sprimeexp,sprimerange[0],sprimerange[1],rans[0]);
  y = CE.DiceYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void SimplePoleBackward::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  weight  = 1./CE.MasslessPropWeight(sprimeexp,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
}

ResonanceUniform::ResonanceUniform(double _mass,double _width,double _deltay1, double _deltay2) :
  mass(_mass), width(_width)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  name     = std::string("ResonanceUniform");
  ms = rans = 0;
};


void ResonanceUniform::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  sprime = CE.MassivePropMomenta(mass,width,1,sprimerange[0],sprimerange[1],rans[0]);
  y = CE.DiceYUniform(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void ResonanceUniform::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  weight  = 1./CE.MassivePropWeight(mass,width,1,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYUniform(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


ResonanceCentral::ResonanceCentral(double _mass,double _width,double _deltay1, double _deltay2) :
  mass(_mass), width(_width)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  name     = std::string("ResonanceCentral");
  ms = rans = 0;
}


void ResonanceCentral::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  sprime = CE.MassivePropMomenta(mass,width,1,sprimerange[0],sprimerange[1],rans[0]);
  y      = CE.DiceYCentral(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void ResonanceCentral::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  weight  = 1./CE.MassivePropWeight(mass,width,1,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYCentral(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 



ResonanceForward::ResonanceForward(double _mass,double _width,double _yexp,double _deltay1, double _deltay2) :
  mass(_mass), width(_width), yexp(_yexp)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  name     = std::string("ResonanceForward");
  ms = rans = 0;
}


void ResonanceForward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  sprime = CE.MassivePropMomenta(mass,width,1,sprimerange[0],sprimerange[1],rans[0]);
  y = CE.DiceYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void ResonanceForward::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  weight  = 1./CE.MassivePropWeight(mass,width,1,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 

ResonanceBackward::ResonanceBackward(double _mass,double _width,double _yexp,double _deltay1, double _deltay2) :
  mass(_mass), width(_width), yexp(_yexp)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  name     = std::string("ResonanceBackward");
  ms = rans = 0;
};


void ResonanceBackward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  sprime = CE.MassivePropMomenta(mass,width,1,sprimerange[0],sprimerange[1],rans[0]);
  y = CE.DiceYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void ResonanceBackward::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  weight  = 1./CE.MassivePropWeight(mass,width,1,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


ThresholdUniform::ThresholdUniform(double _mass,double _deltay1, double _deltay2) :
  mass(_mass)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  name     = std::string("ThresholdUniform");
  ms = rans = 0;
}


void ThresholdUniform::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  sprime = CE.ThresholdMomenta(mass,sprimerange[0],sprimerange[1],rans[0]);
  y = CE.DiceYUniform(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void ThresholdUniform::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  weight  = 1./CE.ThresholdWeight(mass,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYUniform(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


ThresholdCentral::ThresholdCentral(double _mass,double _deltay1, double _deltay2) :
  mass(_mass)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  name     = std::string("ThresholdCentral");
  ms = rans = 0;
}


void ThresholdCentral::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  sprime = CE.ThresholdMomenta(mass,sprimerange[0],sprimerange[1],rans[0]);
  y      = CE.DiceYCentral(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void ThresholdCentral::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  weight  = 1./CE.ThresholdWeight(mass,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYCentral(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 



ThresholdForward::ThresholdForward(double _mass,double _yexp,double _deltay1, double _deltay2) :
  mass(_mass), yexp(_yexp)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  name     = std::string("ThresholdForward");
  ms = rans = 0;
}


void ThresholdForward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  sprime = CE.ThresholdMomenta(mass,sprimerange[0],sprimerange[1],rans[0]);
  y = CE.DiceYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void ThresholdForward::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  weight  = 1./CE.ThresholdWeight(mass,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 

ThresholdBackward::ThresholdBackward(double _mass,double _yexp,double _deltay1, double _deltay2) :
  mass(_mass), yexp(_yexp)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  name     = std::string("ThresholdBackward");
  ms = rans = 0;
}


void ThresholdBackward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  sprime = CE.ThresholdMomenta(mass,sprimerange[0],sprimerange[1],rans[0]);
  y = CE.DiceYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void ThresholdBackward::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  weight  = 1./CE.ThresholdWeight(mass,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


void LLUniform::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double pole = sprimerange[2];
  if (ATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  sprime = CE.LLPropMomenta(1.-beta,pole,sprimerange[0],sprimerange[1],rans[0]);
  y = CE.DiceYUniform(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void LLUniform::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  double pole = sprimerange[2];
  if (ATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  weight  = 1./CE.LLPropWeight(1.-beta,pole,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYUniform(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


void LLCentral::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double pole = sprimerange[2];
  if (ATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  sprime = CE.LLPropMomenta(1.-beta,pole,sprimerange[0],sprimerange[1],rans[0]);
  y      = CE.DiceYCentral(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void LLCentral::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  double pole = sprimerange[2];
  if (ATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  weight  = 1./CE.LLPropWeight(1.-beta,pole,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYCentral(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


void LLForward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double pole = sprimerange[2];
  if (ATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  sprime = CE.LLPropMomenta(1.-beta,pole,sprimerange[0],sprimerange[1],rans[0]);
  y      = CE.DiceYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void LLForward::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  double pole = sprimerange[2];
  if (ATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  weight  = 1./CE.LLPropWeight(1.-beta,pole,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


void LLBackward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double pole = sprimerange[2];
  if (ATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  sprime = CE.LLPropMomenta(1.-beta,pole,sprimerange[0],sprimerange[1],rans[0]);
  y      = CE.DiceYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void LLBackward::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  double pole = sprimerange[2];
  if (ATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  weight  = 1./CE.LLPropWeight(1.-beta,pole,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 




void LBSComptonPeakUniform::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double help   = CE.LLPropMomenta(sexp, sprimerange[2], sprimerange[0], sprimerange[1], rans[0]);
  if (sprimerange[0]<sprimerange[2]*pole && sprimerange[2]*pole<sprimerange[1]) {
    sprime = help - sprimerange[1] + sprimerange[2]*pole;
    if (sprime<sprimerange[0]) 
      sprime = help + (sprimerange[2]*pole - sprimerange[0]);
  }
  else {
    sprime = help;
  }
  y = CE.DiceYUniform(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void LBSComptonPeakUniform::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  double help=sprime;
  if (sprimerange[0]<sprimerange[2]*pole && sprimerange[2]*pole<sprimerange[1]) {
    if (sprime > pole*sprimerange[2]) 
      help = sprime - (sprimerange[2]*pole - sprimerange[0]);
    else
      help = sprime + sprimerange[1] - sprimerange[2]*pole;
  }
  weight  = 1./CE.LLPropWeight(sexp,sprimerange[2], sprimerange[0], sprimerange[1] ,help);
  weight *= 1./sprimerange[2];  
  weight *= CE.WeightYUniform(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 

void LBSComptonPeakCentral::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double help   = CE.LLPropMomenta(sexp, sprimerange[2], sprimerange[0], sprimerange[1], rans[0]);
  if (sprimerange[0]<sprimerange[2]*pole && sprimerange[2]*pole<sprimerange[1]) {
    sprime = help - sprimerange[1] + sprimerange[2]*pole;
    if (sprime<sprimerange[0]) 
      sprime = help + (sprimerange[2]*pole - sprimerange[0]);
  }
  else {
    sprime = help;
  }
  y      = CE.DiceYCentral(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void LBSComptonPeakCentral::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  double help=sprime;
  if (sprimerange[0]<sprimerange[2]*pole && sprimerange[2]*pole<sprimerange[1]) {
    if (sprime > pole*sprimerange[2]) 
      help = sprime - (sprimerange[2]*pole - sprimerange[0]);
    else
      help = sprime + sprimerange[1] - sprimerange[2]*pole;
  }
  weight  = 1./CE.LLPropWeight(sexp,sprimerange[2], sprimerange[0], sprimerange[1] ,help);
  weight *= 1./sprimerange[2];  
  weight *= CE.WeightYCentral(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


void LBSComptonPeakForward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double help   = CE.LLPropMomenta(sexp, sprimerange[2], sprimerange[0], sprimerange[1], rans[0]);
  if (sprimerange[0]<sprimerange[2]*pole && sprimerange[2]*pole<sprimerange[1]) {
    sprime = help - sprimerange[1] + sprimerange[2]*pole;
    if (sprime<sprimerange[0]) 
      sprime = help + (sprimerange[2]*pole - sprimerange[0]);
  }
  else {
    sprime = help;
  }
  y = CE.DiceYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void LBSComptonPeakForward::GenerateWeight(double sprime,double y,int mode)
{
  weight      = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) {
    return;
  }
  double help=sprime;
  if (sprimerange[0]<sprimerange[2]*pole && sprimerange[2]*pole<sprimerange[1]) {
    if (sprime > pole*sprimerange[2]) 
      help = sprime - (sprimerange[2]*pole - sprimerange[0]);
    else
      help = sprime + sprimerange[1] - sprimerange[2]*pole;
  }
  weight  = 1./CE.LLPropWeight(sexp,sprimerange[2], sprimerange[0], sprimerange[1] ,help);
  weight *= 1./sprimerange[2];  
  weight *= CE.WeightYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


void LBSComptonPeakBackward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double help   = CE.LLPropMomenta(sexp, sprimerange[2], sprimerange[0], sprimerange[1], rans[0]);
  if (sprimerange[0]<sprimerange[2]*pole && sprimerange[2]*pole<sprimerange[1]) {
    sprime = help - sprimerange[1] + sprimerange[2]*pole;
    if (sprime<sprimerange[0]) 
      sprime = help + (sprimerange[2]*pole - sprimerange[0]);
  }
  else {
    sprime = help;
  }
  y = CE.DiceYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void LBSComptonPeakBackward::GenerateWeight(double sprime,double y,int mode)
{
  weight      = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) {
    return;
  }
  double help=sprime;
  if (sprimerange[0]<sprimerange[2]*pole && sprimerange[2]*pole<sprimerange[1]) {
    if (sprime > pole*sprimerange[2]) 
      help = sprime - (sprimerange[2]*pole - sprimerange[0]);
    else
      help = sprime + sprimerange[1] - sprimerange[2]*pole;
  }
  weight  = 1./CE.LLPropWeight(sexp,sprimerange[2], sprimerange[0], sprimerange[1] ,help);
  weight *= 1./sprimerange[2];  
  weight *= CE.WeightYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 

#include "ISR_Channel.H"
#include "Channel_Elements.H"
#include "Channel_Basics.H"
#include "Run_Parameter.H"
#include "Message.H"

#include <stdio.h>

using namespace PHASIC;
using namespace AORGTOOLS;
using namespace std;

SimplePoleUniform::SimplePoleUniform(double _sprimeexp, double _deltay1, double _deltay2) :
  sprimeexp(_sprimeexp)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  char help[3];
  sprintf(help,"%i",int(100.*sprimeexp));
  name     = std::string("SimplePoleUniform"+std::string(help));
  ms = rans = 0;
  msg.Out()<<"Init Simple_Pole : "<<sprimeexp<<" / "<<deltay[0]<<" / "<<deltay[1]<<endl;
};


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
  msg.Out()<<"Init Simple_Pole : "<<sprimeexp<<" / "<<deltay[0]<<" / "<<deltay[1]<<endl;
};


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
  /*
    if (weight<0.) {
    msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight.";
    msg.Error()<<"with  sprime , y = "<<sprime<<" , "<<y<<endl;
    }
  */
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
};


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
  /*
    if (weight<0.) {
    msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight.";
    msg.Error()<<"with  sprime , y = "<<sprime<<" , "<<y<<endl;
    }
  */
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
};

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
  /*
    if (weight<0.) {
    msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight.";
    msg.Error()<<"with  sprime , y = "<<sprime<<" , "<<y<<endl;
    }
  */
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
};


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
  /*
    if (weight<0.) {
    msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight.";
    msg.Error()<<"with  sprime , y = "<<sprime<<" , "<<y<<endl;
    }
  */
} 



ResonanceForward::ResonanceForward(double _mass,double _width,double _yexp,double _deltay1, double _deltay2) :
  mass(_mass), width(_width), yexp(_yexp)
{
  deltay[0] = _deltay1;
  deltay[1] = _deltay2;
  name     = std::string("ResonanceForward");
  ms = rans = 0;
};


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
  /*
    if (weight<0.) {
    msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight.";
    msg.Error()<<"with  sprime , y = "<<sprime<<" , "<<y<<endl;
    }
  */
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
  /*
    if (weight<0.) {
    msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight." ;
    }
    msg.Error()<<"****  sprime , y = "<<sprime<<" , "<<y<<endl;
  */
} 


void LLUniform::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double pole = sprimerange[2];
  if (AMATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  sprime = CE.LLPropMomenta(1.-beta,pole,sprimerange[0],sprimerange[1],rans[0]);
  y = CE.DiceYUniform(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void LLUniform::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  double pole = sprimerange[2];
  if (AMATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  weight  = 1./CE.LLPropWeight(1.-beta,pole,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYUniform(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


void LLCentral::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double pole = sprimerange[2];
  if (AMATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  sprime = CE.LLPropMomenta(1.-beta,pole,sprimerange[0],sprimerange[1],rans[0]);
  y      = CE.DiceYCentral(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void LLCentral::GenerateWeight(double sprime,double y,int mode)
{
  //AORGTOOLS::msg.Out()<<"Cen : "<<endl;GetRange(); 
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  double pole = sprimerange[2];
  if (AMATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  weight  = 1./CE.LLPropWeight(1.-beta,pole,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYCentral(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


void LLForward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double pole = sprimerange[2];
  if (AMATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  sprime = CE.LLPropMomenta(1.-beta,pole,sprimerange[0],sprimerange[1],rans[0]);
  y      = CE.DiceYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void LLForward::GenerateWeight(double sprime,double y,int mode)
{
  //AORGTOOLS::msg.Out()<<"FW  : "<<endl;GetRange(); 
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  double pole = sprimerange[2];
  if (AMATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  weight  = 1./CE.LLPropWeight(1.-beta,pole,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


void LLBackward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double pole = sprimerange[2];
  if (AMATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  sprime = CE.LLPropMomenta(1.-beta,pole,sprimerange[0],sprimerange[1],rans[0]);
  y      = CE.DiceYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void LLBackward::GenerateWeight(double sprime,double y,int mode)
{
  //AORGTOOLS::msg.Out()<<"BW  : "<<endl;GetRange(); 
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  double pole = sprimerange[2];
  if (AMATOOLS::IsEqual(sprimerange[2],sprimerange[1])) pole *= factor;
  weight  = 1./CE.LLPropWeight(1.-beta,pole,sprimerange[0],sprimerange[1],sprime);
  weight *= 1./sprimerange[2];
  weight *= CE.WeightYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 




void LBSComptonPeakUniform::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  /*! former routine
    ! doesn't work since the pole is far off the region of integration
    offset      = (-1.+pole)*sprimerange[2];
    double help = CE.LLPropMomenta(1.,sprimerange[2],
    sprimerange[0]+offset,sprimerange[1]+offset,ran[0]);
    if (help<sprimerange[0]) help += (sprimerange[1]-sprimerange[0]);
    sprime      = help;
  */
  /*! new version
    // tested; ok, but doesn't comply with beam handler
    // will be all right, if 'sprimerange[2]' is chosen to be the end of the spectrum 
    // and 'pole' denotes the coordinate of the compton peak relative to current 'sprimerange[2]'
    double help   = CE.LLPropMomenta(1., sprimerange[2]*1.000001, sprimerange[0], sprimerange[2], rans[0]);
    double shift  = sprimerange[2] * (1. - pole);
    if (help<shift) sprime = (sprimerange[2] - shift) + help
    else sprime = help - shift;
  */
  double help   = sprimerange[2] * pole; 
  sprime        = CE.LLPropMomenta(sexp, sprimerange[2], sprimerange[0], help, rans[0]);
  y = CE.DiceYUniform(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void LBSComptonPeakUniform::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) return;
  /*! former routine
    offset      = (-1.+pole)*sprimerange[2];
    double help = sprime;
    if (help>sprimerange[1]+offset) help -= (sprimerange[1]-sprimerange[0]);
    weight  = 1./CE.LLPropWeight(1.,sprimerange[2],
    sprimerange[0]+offset,sprimerange[1]+offset,help);
  */
  /*! new version
    // tested; ok, but doesn't comply with beam handler
    // will be all right, if 'sprimerange[2]' is chosen to be the end of the spectrum 
    // and 'pole' denotes the coordinate of the compton peak relative to current 'sprimerange[2]'
    double help   = sprime;
    double shift  = sprimerange[2] * pole;
    if (sprime>shift) help -= shift
    else help += sprimerange[2]-shift;
    weight = 1./CE.LLPropWeight(1., sprimerange[2]*1.000001, sprimerange[0], sprimerange[2], help);
  */
  double help = sprimerange[2] * pole;
  weight  = 1./CE.LLPropWeight(sexp,sprimerange[2], sprimerange[0], help ,sprime);
  weight *= 1./help;
  weight *= CE.WeightYUniform(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 

void LBSComptonPeakCentral::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double help   = sprimerange[2] * pole; 
  sprime        = CE.LLPropMomenta(sexp, sprimerange[2], sprimerange[0], help, rans[0]);
  y      = CE.DiceYCentral(sprime/sprimerange[2], yrange, deltay, mode, rans[1]);
}

void LBSComptonPeakCentral::GenerateWeight(double sprime,double y,int mode)
{
  weight  = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) {
    //AORGTOOLS::msg.Out()<<" Cen : sprime, y0, y : "<<sprime/sprimerange[2]<<" , "<<sprimerange[1]/sprimerange[2]<<","
    //			<<0.5*log(sprime/sprimerange[2])<<" , "<<y<<" -> "<<weight<<endl;
    return;
  }
  double help = sprimerange[2] * pole;
  weight  = 1./CE.LLPropWeight(sexp,sprimerange[2], sprimerange[0], help ,sprime);
  weight *= 1./help;
  weight *= CE.WeightYCentral(sprime/sprimerange[2], yrange, deltay, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


void LBSComptonPeakForward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double help   = sprimerange[2] * pole; 
  sprime        = CE.LLPropMomenta(sexp, sprimerange[2], sprimerange[0], help, rans[0]);
  y = CE.DiceYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void LBSComptonPeakForward::GenerateWeight(double sprime,double y,int mode)
{
  weight      = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) {
    return;
  }
  double help = sprimerange[2] * pole;
  weight  = 1./CE.LLPropWeight(sexp,sprimerange[2], sprimerange[0], help ,sprime);
  weight *= 1./help;
  weight *= CE.WeightYForward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 


void LBSComptonPeakBackward::GeneratePoint(double & sprime,double & y,int mode,double * rans)
{
  double help   = sprimerange[2] * pole; 
  sprime        = CE.LLPropMomenta(sexp, sprimerange[2], sprimerange[0], help, rans[0]);
  y = CE.DiceYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, rans[1]);
}

void LBSComptonPeakBackward::GenerateWeight(double sprime,double y,int mode)
{
  weight      = 0.;
  if ((sprime<sprimerange[0]) || (sprime>sprimerange[1])) {
    return;
  }
  double help = sprimerange[2] * pole;
  weight  = 1./CE.LLPropWeight(sexp,sprimerange[2], sprimerange[0], help ,sprime);
  weight *= 1./help;
  weight *= CE.WeightYBackward(sprime/sprimerange[2], yrange, deltay, yexp, mode, y);
  if (weight<0.) msg.Error()<<"Negative weight in "<<name<<"::GenerateWeight."<<endl;
} 

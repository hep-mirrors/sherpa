#include "PHASIC++/Scales/Tag_Setter.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Particle_Dresser.H"
#include "PHASIC++/Process/Process_Base.H"

using namespace PHASIC;
using namespace ATOOLS;

std::string Tag_Setter::ReplaceTags(std::string &expr) const
{
  return p_calc->ReplaceTags(expr);
}

Term *Tag_Setter::ReplaceTags(Term *term) const
{
  if (term->Id()>=10) {
    if (term->Id()>=100) {
      term->Set(p_setter->Momenta()[term->Id()-100]);
      return term;
    }
    term->Set(p_setter->Scales()[term->Id()-10]);
    return term;
  }
  switch (term->Id()) {
  case 1:
    term->Set(p_setter->Scale(stp::fac));
    return term;
  case 2:
    term->Set(p_setter->Scale(stp::ren));
    return term;
  case 3:
    term->Set(p_setter->Scale(stp::res));
    return term;
  case 4:
    term->Set(sqr(p_setter->HTM()));
    return term;
  case 5:
    term->Set(sqr(p_setter->HT()));
    return term;
  case 6:
    term->Set(sqr(p_setter->HTMprime()));
    return term;
  case 7:
    term->Set(sqr(p_setter->HTprime()));
    return term;
  case 8:
    term->Set(p_setter->PSum());
    return term;
  case 9:
    term->Set(p_setter->PSum());
    return term;
  }
  return term;
}

void Tag_Setter::AssignId(Term *term)
{
  if (term->Tag()=="MU_F2") term->SetId(1);
  else if (term->Tag()=="MU_R2") term->SetId(2);
  else if (term->Tag()=="MU_Q2") term->SetId(3);
  else if (term->Tag()=="H_TM2") term->SetId(4);
  else if (term->Tag()=="H_T2")  term->SetId(5);
  else if (term->Tag()=="H_TMp2") term->SetId(6);
  else if (term->Tag()=="H_Tp2") term->SetId(7);
  else if (term->Tag()=="P_SUM") term->SetId(8);
  else if (term->Tag()=="TAUB") term->SetId(9);
  else {
    term->SetId(100+ToType<int>
		(term->Tag().substr
		 (2,term->Tag().length()-3)));
  }
}

namespace PHASIC {

  class H_TY2: public Function {
  private:

    Scale_Setter_Base *p_setter;

  public:

    inline H_TY2(Scale_Setter_Base *const setter):
      Function("H_TY2"), p_setter(setter) {}

    Term *Evaluate(const std::vector<Term*> &args) const
    {
      double htyfac(args[0]->Get<double>()), htyexp(args[1]->Get<double>());
      Vec4D psum(0.,0.,0.,0.);
      const Vec4D_Vector &p(p_setter->Momenta());
      for (size_t i(p_setter->NIn());i<p.size();++i) psum+=p[i];
      double yboost((psum/(double)(p.size()-p_setter->NIn())).Y());
      double hty(0.0);
      for (size_t i(p_setter->NIn());i<p.size();++i) 
        hty+=p[i].PPerp()*exp(htyfac*pow(std::abs(p[i].Y()-yboost),htyexp));
      Term *res(Term::New(sqr(hty)));
      p_interpreter->AddTerm(res);
      return res;
    }

  };// end of class H_TY2

  class Dressed_H_Tp2: public Function {
  private:

    size_t m_l1,m_l2;
    std::vector<size_t> m_photons,m_charges;
    Particle_Dresser *p_dresser;
    Scale_Setter_Base *p_setter;

  public:

    Dressed_H_Tp2(Scale_Setter_Base *const setter):
      Function("DH_Tp2"), m_l1(0), m_l2(0), p_dresser(NULL), p_setter(setter)
    {
      DEBUG_FUNC(p_setter->Process()->Name());
      size_t nl(0);
      for (size_t i(p_setter->Process()->NIn());
           i<p_setter->Process()->Flavours().size();++i) {
        if (p_setter->Process()->Flavours()[i].IsLepton()) {
          nl++;
          if      (nl==1) m_l1=i;
          else if (nl==2) m_l2=i;
          else           {m_l1=m_l2=0;}
        }
        if      (p_setter->Process()->Flavours()[i].Charge())
          m_charges.push_back(i);
        else if (p_setter->Process()->Flavours()[i].IsPhoton())
          m_photons.push_back(i);
      }
      msg_Debugging()<<"Found "<<nl<<" leptons: "<<m_l1<<" "<<m_l2<<std::endl;
    }

    Vec4D_Vector ConeDress(const Vec4D_Vector& p,
                           std::vector<double>& dR2s) const
    {
      DEBUG_FUNC("photons: "<<m_photons<<", charges: "<<m_charges);
      if (!m_photons.size() || !m_charges.size()) return p;
      Vec4D_Vector pp(p);
      std::vector<bool> valid(m_photons.size(),true);
      std::vector<std::vector<double> > dij;
      dij.resize(m_charges.size());
      double maxd(std::numeric_limits<double>::max()),dmin(maxd);
      size_t ii(0),jj(0),max(std::numeric_limits<size_t>::max());
      // calculate initial dijs=dR(i,j)^2
      for (size_t i(0);i<m_charges.size();++i) {
        dij[i].resize(m_photons.size());
        for (size_t j(0);j<m_photons.size();++j) {
          dij[i][j]=pp[m_charges[i]].DR2(pp[m_photons[j]])/dR2s[i];
          if (dij[i][j]<dmin) { dmin=dij[i][j]; ii=i; jj=j; }
        }
      }
      msg_Debugging()<<"dmin = "<<dmin<<std::endl;
      while (dmin<1.) {
        if (msg_LevelIsDebugging()) {
          msg_Out()<<"dij: ";
          for (size_t i(0);i<dij.size();++i) {
            msg_Out()<<dij[i]<<"\n     ";
          }
          msg_Out()<<"-> i: "<<ii<<" , j: "<<jj<<" , dmin="<<dmin<<std::endl;
        }
        // mark photon that is recombined
        valid[jj]=false;
        // recombine, do not recompute always with respect to bare axis
        pp[m_charges[ii]]+=pp[m_photons[jj]];
        pp[m_photons[jj]]=Vec4D(0.,0.,0.,0.);
        for (size_t i(0);i<m_charges.size();++i) dij[i][jj]=maxd;
        // find new dmin
        dmin=maxd;
        for (size_t i(0);i<m_charges.size();++i) {
          for (size_t j(0);j<m_photons.size();++j) if (valid[j]) {
            if (dij[i][j]<dmin) { dmin=dij[i][j]; ii=i; jj=j; }
          }
        }
      }
      return pp;
    }

    Vec4D_Vector RecombinationDress(const Vec4D_Vector& p,
                                    const double& exp,
                                    std::vector<double>& dR2s) const
    {
      THROW(not_implemented,"Not implemented.");
      return p;
    }

    Term *Evaluate(const std::vector<Term*> &args) const
    {
      DEBUG_FUNC(p_setter->Process()->Name()<<" "<<m_l1<<" "<<m_l2);
      if (m_l1==0 || m_l2==0) THROW(fatal_error,"Lepton indices not set.");
      Vec4D_Vector p(p_setter->Momenta());
      if (m_photons.size() && m_charges.size() && args.size()) {
        msg_Debugging()<<"Reading arguments."<<std::endl;
        std::string method(args[0]->Get<std::string>());
        double dRglobal(args[1]->Get<double>());
        std::vector<double> dR2s(m_charges.size(),sqr(dRglobal));
        for (size_t i(2);i<args.size();i+=2) {
          kf_code kf(args[i]->Get<double>());
          double dR(args[i+1]->Get<double>());
          msg_Debugging()<<"Setting dR="<<dR<<" for kf="<<kf<<std::endl;
          for (size_t j(0);j<m_charges.size();++j) {
            if (kf==p_setter->Process()->Flavours()[m_charges[j]].Kfcode()) {
              dR2s[j]=sqr(dR);
            }
          }
        }
        if (msg_LevelIsDebugging()) {
          msg_Out()<<"Cone radii:"<<std::endl;
          for (size_t i(0);i<m_charges.size();++i)
            msg_Out()<<i<<": "<<p_setter->Process()->Flavours()[m_charges[i]]
                     <<" -> dR="<<sqrt(dR2s[i])<<std::endl;
          msg_Out()<<"Original momenta:"<<std::endl;
          for (size_t i(0);i<p_setter->Momenta().size();++i)
            msg_Out()<<p_setter->Process()->Flavours()[i]<<" "
                     <<p[i]<<std::endl;
        }
        if      (method=="Cone")   p=ConeDress(p,dR2s);
        else if (method=="kt")     p=RecombinationDress(p, 1.,dR2s);
        else if (method=="CA")     p=RecombinationDress(p, 0.,dR2s);
        else if (method=="antikt") p=RecombinationDress(p,-1.,dR2s);
        else THROW(fatal_error,"Dressing method not implemented.");
        if (msg_LevelIsDebugging())
          for (size_t i(0);i<p_setter->Momenta().size();++i)
            msg_Out()<<p_setter->Process()->Flavours()[i]<<" "
                     <<p[i]<<std::endl;
      }
      double htp((p[m_l1]+p[m_l2]).MPerp());
      for (size_t i(p_setter->Process()->NIn());i<p.size();++i)
        if (i!=m_l1 && i!=m_l2) htp+=p[i].PPerp();
      Term *res(Term::New(sqr(htp)));
      p_interpreter->AddTerm(res);
      return res;
    }

  };// end of class Dressed_H_Tp2

  class Dressed_H_Tln2: public Function {
  private:

    std::vector<size_t> m_leptons,m_neutrinos,m_photons,m_charges;
    Particle_Dresser *p_dresser;
    Scale_Setter_Base *p_setter;

  public:

    Dressed_H_Tln2(Scale_Setter_Base *const setter):
      Function("DH_Tln2"), p_dresser(NULL), p_setter(setter)
    {
      DEBUG_FUNC(p_setter->Process()->Name());
      for (size_t i(p_setter->Process()->NIn());
           i<p_setter->Process()->Flavours().size();++i) {
        if      (p_setter->Process()->Flavours()[i].IsChargedLepton())
          m_leptons.push_back(i);
        else if (p_setter->Process()->Flavours()[i].IsNeutrino())
          m_neutrinos.push_back(i);
        if      (p_setter->Process()->Flavours()[i].Charge())
          m_charges.push_back(i);
        else if (p_setter->Process()->Flavours()[i].IsPhoton())
          m_photons.push_back(i);
      }
      msg_Debugging()<<"Found "<<m_leptons.size()<<" leptons and "
                               <<m_neutrinos.size()<<" neutrinos."<<std::endl;
    }

    Vec4D_Vector ConeDress(const Vec4D_Vector& p,
                           std::vector<double>& dR2s) const
    {
      DEBUG_FUNC("photons: "<<m_photons<<", charges: "<<m_charges);
      if (!m_photons.size() || !m_charges.size()) return p;
      Vec4D_Vector pp(p);
      std::vector<bool> valid(m_photons.size(),true);
      std::vector<std::vector<double> > dij;
      dij.resize(m_charges.size());
      double maxd(std::numeric_limits<double>::max()),dmin(maxd);
      size_t ii(0),jj(0),max(std::numeric_limits<size_t>::max());
      // calculate initial dijs=dR(i,j)^2
      for (size_t i(0);i<m_charges.size();++i) {
        dij[i].resize(m_photons.size());
        for (size_t j(0);j<m_photons.size();++j) {
          dij[i][j]=pp[m_charges[i]].DR2(pp[m_photons[j]])/dR2s[i];
          if (dij[i][j]<dmin) { dmin=dij[i][j]; ii=i; jj=j; }
        }
      }
      msg_Debugging()<<"dmin = "<<dmin<<std::endl;
      while (dmin<1.) {
        if (msg_LevelIsDebugging()) {
          msg_Out()<<"dij: ";
          for (size_t i(0);i<dij.size();++i) {
            msg_Out()<<dij[i]<<"\n     ";
          }
          msg_Out()<<"-> i: "<<ii<<" , j: "<<jj<<" , dmin="<<dmin<<std::endl;
        }
        // mark photon that is recombined
        valid[jj]=false;
        // recombine, do not recompute always with respect to bare axis
        pp[m_charges[ii]]+=pp[m_photons[jj]];
        pp[m_photons[jj]]=Vec4D(0.,0.,0.,0.);
        for (size_t i(0);i<m_charges.size();++i) dij[i][jj]=maxd;
        // find new dmin
        dmin=maxd;
        for (size_t i(0);i<m_charges.size();++i) {
          for (size_t j(0);j<m_photons.size();++j) if (valid[j]) {
            if (dij[i][j]<dmin) { dmin=dij[i][j]; ii=i; jj=j; }
          }
        }
      }
      return pp;
    }

    Vec4D_Vector RecombinationDress(const Vec4D_Vector& p,
                                    const double& exp,
                                    std::vector<double>& dR2s) const
    {
      THROW(not_implemented,"Not implemented.");
      return p;
    }

    Term *Evaluate(const std::vector<Term*> &args) const
    {
      if (m_leptons.size()==0) THROW(fatal_error,"No leptons available.");
      DEBUG_FUNC(p_setter->Process()->Name());
      Vec4D_Vector p(p_setter->Momenta());
      if (m_photons.size() && m_charges.size() && args.size()) {
        msg_Debugging()<<"Reading arguments."<<std::endl;
        std::string method(args[0]->Get<std::string>());
        double dRglobal(args[1]->Get<double>());
        std::vector<double> dR2s(m_charges.size(),sqr(dRglobal));
        for (size_t i(2);i<args.size();i+=2) {
          kf_code kf(args[i]->Get<double>());
          double dR(args[i+1]->Get<double>());
          msg_Debugging()<<"Setting dR="<<dR<<" for kf="<<kf<<std::endl;
          for (size_t j(0);j<m_charges.size();++j) {
            if (kf==p_setter->Process()->Flavours()[m_charges[j]].Kfcode()) {
              dR2s[j]=sqr(dR);
            }
          }
        }
        if (msg_LevelIsDebugging()) {
          msg_Out()<<"Cone radii:"<<std::endl;
          for (size_t i(0);i<m_charges.size();++i)
            msg_Out()<<i<<": "<<p_setter->Process()->Flavours()[m_charges[i]]
                     <<" -> dR="<<sqrt(dR2s[i])<<std::endl;
          msg_Out()<<"Original momenta:"<<std::endl;
          for (size_t i(0);i<p_setter->Momenta().size();++i)
            msg_Out()<<p_setter->Process()->Flavours()[i]<<" "
                     <<p[i]<<std::endl;
        }
        if      (method=="Cone")   p=ConeDress(p,dR2s);
        else if (method=="kt")     p=RecombinationDress(p, 1.,dR2s);
        else if (method=="CA")     p=RecombinationDress(p, 0.,dR2s);
        else if (method=="antikt") p=RecombinationDress(p,-1.,dR2s);
        else THROW(fatal_error,"Dressing method not implemented.");
        if (msg_LevelIsDebugging())
          for (size_t i(0);i<p_setter->Momenta().size();++i)
            msg_Out()<<p_setter->Process()->Flavours()[i]<<" "
                     <<p[i]<<std::endl;
      }
      double htp(0);
      for (size_t i(0);i<m_leptons.size();++i) htp += p[m_leptons[i]].PPerp();
      Vec4D ptmis(0.,0.,0.,0.);
      for (size_t i(0);i<m_neutrinos.size();++i) ptmis += p[m_neutrinos[i]];
      msg_Debugging()<<"H_{T,ln} = "<<htp+ptmis.PPerp()<<std::endl;
      Term *res(Term::New(sqr(htp+ptmis.PPerp())));
      p_interpreter->AddTerm(res);
      return res;
    }

  };// end of class Dressed_H_Tln2

  class Dressed_MPerp2: public Function {
  private:

    size_t m_l1,m_l2;
    std::vector<size_t> m_photons,m_charges;
    Scale_Setter_Base *p_setter;

  public:

    Dressed_MPerp2(Scale_Setter_Base *const setter):
      Function("DMPerp2"), m_l1(0), m_l2(0), p_setter(setter)
    {
      DEBUG_FUNC(p_setter->Process()->Name());
      size_t nl(0);
      for (size_t i(p_setter->Process()->NIn());
           i<p_setter->Process()->Flavours().size();++i) {
        if (p_setter->Process()->Flavours()[i].IsLepton()) {
          nl++;
          if      (nl==1) m_l1=i;
          else if (nl==2) m_l2=i;
          else           {m_l1=m_l2=0;}
        }
        if      (p_setter->Process()->Flavours()[i].Charge())
          m_charges.push_back(i);
        else if (p_setter->Process()->Flavours()[i].IsPhoton())
          m_photons.push_back(i);
      }
    }

    Vec4D_Vector ConeDress(const Vec4D_Vector& p,
                           std::vector<double>& dR2s) const
    {
      DEBUG_FUNC("photons: "<<m_photons<<", charges: "<<m_charges);
      if (!m_photons.size() || !m_charges.size()) return p;
      Vec4D_Vector pp(p);
      std::vector<bool> valid(m_photons.size(),true);
      std::vector<std::vector<double> > dij;
      dij.resize(m_charges.size());
      double maxd(std::numeric_limits<double>::max()),dmin(maxd);
      size_t ii(0),jj(0),max(std::numeric_limits<size_t>::max());
      // calculate initial dijs=dR(i,j)^2
      for (size_t i(0);i<m_charges.size();++i) {
        dij[i].resize(m_photons.size());
        for (size_t j(0);j<m_photons.size();++j) {
          dij[i][j]=pp[m_charges[i]].DR2(pp[m_photons[j]])/dR2s[i];
          if (dij[i][j]<dmin) { dmin=dij[i][j]; ii=i; jj=j; }
        }
      }
      msg_Debugging()<<"dmin = "<<dmin<<std::endl;
      while (dmin<1.) {
        if (msg_LevelIsDebugging()) {
          msg_Out()<<"dij: ";
          for (size_t i(0);i<dij.size();++i) {
            msg_Out()<<dij[i]<<"\n     ";
          }
          msg_Out()<<"-> i: "<<ii<<" , j: "<<jj<<" , dmin="<<dmin<<std::endl;
        }
        // mark photon that is recombined
        valid[jj]=false;
        // recombine, do not recompute always with respect to bare axis
        pp[m_charges[ii]]+=pp[m_photons[jj]];
        pp[m_photons[jj]]=Vec4D(0.,0.,0.,0.);
        for (size_t i(0);i<m_charges.size();++i) dij[i][jj]=maxd;
        // find new dmin
        dmin=maxd;
        for (size_t i(0);i<m_charges.size();++i) {
          for (size_t j(0);j<m_photons.size();++j) if (valid[j]) {
            if (dij[i][j]<dmin) { dmin=dij[i][j]; ii=i; jj=j; }
          }
        }
      }
      return pp;
    }

    Vec4D_Vector RecombinationDress(const Vec4D_Vector& p,
                                    const double& exp,
                                    std::vector<double>& dR2s) const
    {
      THROW(not_implemented,"Not implemented.");
      return p;
    }

    Term *Evaluate(const std::vector<Term*> &args) const
    {
      if (m_l1==0 || m_l2==0) THROW(fatal_error,"Lepton indices not set.");
      DEBUG_FUNC(p_setter->Process()->Name()<<" "<<m_l1<<" "<<m_l2);
      Vec4D_Vector p(p_setter->Momenta());
      if (m_photons.size() && m_charges.size() && args.size()) {
        msg_Debugging()<<"Reading arguments."<<std::endl;
        std::string method(args[0]->Get<std::string>());
        double dRglobal(args[1]->Get<double>());
        std::vector<double> dR2s(m_charges.size(),sqr(dRglobal));
        for (size_t i(2);i<args.size();i+=2) {
          kf_code kf(args[i]->Get<double>());
          double dR(args[i+1]->Get<double>());
          msg_Debugging()<<"Setting dR="<<dR<<" for kf="<<kf<<std::endl;
          for (size_t j(0);j<m_charges.size();++j) {
            if (kf==p_setter->Process()->Flavours()[m_charges[j]].Kfcode()) {
              dR2s[j]=sqr(dR);
            }
          }
        }
        if (msg_LevelIsDebugging()) {
          msg_Out()<<"Cone radii:"<<std::endl;
          for (size_t i(0);i<m_charges.size();++i)
            msg_Out()<<i<<": "<<p_setter->Process()->Flavours()[m_charges[i]]
                     <<" -> dR="<<sqrt(dR2s[i])<<std::endl;
          msg_Out()<<"Original momenta:"<<std::endl;
          for (size_t i(0);i<p_setter->Momenta().size();++i)
            msg_Out()<<p_setter->Process()->Flavours()[i]<<" "
                     <<p[i]<<std::endl;
        }
        if      (method=="Cone")   p=ConeDress(p,dR2s);
        else if (method=="kt")     p=RecombinationDress(p, 1.,dR2s);
        else if (method=="CA")     p=RecombinationDress(p, 0.,dR2s);
        else if (method=="antikt") p=RecombinationDress(p,-1.,dR2s);
        else THROW(fatal_error,"Dressing method not implemented.");
        if (msg_LevelIsDebugging())
          for (size_t i(0);i<p_setter->Momenta().size();++i)
            msg_Out()<<p_setter->Process()->Flavours()[i]<<" "
                     <<p[i]<<std::endl;
      }
      Term *res(Term::New((p[m_l1]+p[m_l2]).MPerp2()));
      p_interpreter->AddTerm(res);
      return res;
    }

  };// end of class Dressed_MPerp2
}

void Tag_Setter::SetTags(Algebra_Interpreter *const calc)
{
  calc->AddTag("MU_F2","1.0");
  calc->AddTag("MU_R2","1.0");
  calc->AddTag("MU_Q2","1.0");
  calc->AddTag("H_T2","1.0");
  calc->AddTag("H_TM2","1.0");
  calc->AddTag("H_Tp2","1.0");
  calc->AddTag("H_TMp2","1.0");
  calc->AddTag("P_SUM","(1.0,0.0,0.0,0.0)");
  calc->AddTag("TAUB","1.0");
  calc->AddFunction(new H_TY2(p_setter));
  calc->AddFunction(new Dressed_H_Tp2(p_setter));
  calc->AddFunction(new Dressed_H_Tln2(p_setter));
  calc->AddFunction(new Dressed_MPerp2(p_setter));
  for (size_t i=0;i<p_setter->Scales().size();++i)
    calc->AddTag("MU_"+ToString(i)+"2","1.0");
  for (size_t i=0;i<p_setter->NIn()+p_setter->NOut();++i) 
    calc->AddTag("p["+ToString(i)+"]",ToString(Vec4D()));
}

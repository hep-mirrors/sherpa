#include "Calculator.H"
#include "String_Generator.H"

using namespace AMEGIC;

VVVV_Calc::VVVV_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type     = zl::VVVV;
  ncoupl=16;narg=8;pn=5;
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gauge3));
  lorentzlist.push_back(Lorentz_Function(lf::Gauge3));
  for (short int i=0;i<4;i++) lorentzlist[i].SetParticleArg(i);
  lorentzlist[4].SetParticleArg(1,0,4);           
  lorentzlist[5].SetParticleArg(-4,2,3);     
}

Kabbala VVVV_Calc::GGGG() 
{   
  return 
    +Z(1,0)*( (X(3,1)-X(3,0))*(X(2,3)+X(2,4))
	      +(X(3,2)+X(3,4))*(X(2,0)-X(2,1))
	      )
    +Z(3,2)*( (X(1,3)-X(1,2))*(X(0,1)-X(0,4))
	      +(X(1,0)-X(1,4))*(X(0,2)-X(0,3))
	      )
    +(X(3,2)+X(3,4))*(Z(2,0)*(X(1,4)-X(1,0))
		      -Z(2,1)*(X(0,4)-X(0,1))
		      )
    +(X(2,3)+X(2,4))*(Z(3,0)*(X(1,0)-X(1,4))
		      -Z(3,1)*(X(0,1)-X(0,4))
		      )   
    +Z(3,2)*Z(1,0)*(V(0,3)-V(0,2)+V(1,2)-V(1,3));	    
}

Kabbala VVVV_Calc::Do() 
{
  Kabbala factor = sgen->GetEnumber(coupl[8])*sgen->GetEnumber(coupl[9]);

  if (IsZero(M(0)) &&
      IsZero(M(1)) &&
      IsZero(M(2)) &&
      IsZero(M(3)) &&
      IsZero(M(4))) {
    return factor*GGGG();
  }

  return 
    factor*(
	    M(0)*M(1)*M(2)*M(3)*X(0,0)*X(1,1)*X(2,2)*X(3,3)*(+V(0,4)*(V(1,3)*V(2,4)-V(1,2)*V(3,4))
							     -V(1,4)*(V(0,3)*V(2,4)-V(0,2)*V(3,4))
							     ) //% 4
	    +Triple_M1(0,1,2,3,4)-Triple_M1(0,1,3,2,4)+Triple_M1(0,2,3,1,4)-Triple_M1(1,2,3,0,4)    
	    +Triple_M2(0,1,2,3,4)-Triple_M2(1,0,2,3,4)+Triple_M2(1,0,3,2,4)-Triple_M2(0,1,3,2,4)
	    +Double_M1(0,1,2,3,4)-Double_M1(2,3,0,1,4)+Double_M2(0,1,2,3,4)-Double_M2(1,0,2,3,4)
	    +Double_M2(1,0,3,2,4)-Double_M2(0,1,3,2,4)+Double_M3(0,1,2,3,4)-Double_M3(1,0,2,3,4)    
	    -Double_M3(3,2,0,1,4)+Double_M3(2,3,0,1,4)+Single_M(0,1,2,3,4)-Single_M(1,0,2,3,4)
	    +Single_M(3,2,1,0,-4)-Single_M(2,3,1,0,-4)
	    
	    +M(4)*( (V(2,4)-V(3,4))*Z(3,2)*(X(0,1)*X(1,4)-X(1,0)*X(0,4))
		   +(V(0,4)-V(1,4))*Z(1,0)*(X(2,3)*X(3,4)-X(3,2)*X(2,4))
		   +Z(1,0)*Z(3,2)*( V(0,4)*(V(2,4)-V(3,4))
				   +V(1,4)*(V(3,4)-V(2,4))
				   )
		   +(X(0,1)*X(1,4)-X(0,4)*X(1,0))*(X(2,3)*X(3,4)-X(2,4)*X(3,2))
		   )
	    
	    +Z(1,0)*( (X(3,1)-X(3,0))*(X(2,3)+X(2,4))
		     +(X(3,2)+X(3,4))*(X(2,0)-X(2,1))
		     )
	    +Z(3,2)*( (X(1,3)-X(1,2))*(X(0,1)-X(0,4))
		     +(X(1,0)-X(1,4))*(X(0,2)-X(0,3))
		     )
	    +(X(3,2)+X(3,4))*(Z(2,0)*(X(1,4)-X(1,0))
			      -Z(2,1)*(X(0,4)-X(0,1))
			      )
	    
	    +(X(2,3)+X(2,4))*(Z(3,0)*(X(1,0)-X(1,4))
			      -Z(3,1)*(X(0,1)-X(0,4))
			      )
	    
	    +Z(3,2)*Z(1,0)*(V(0,3)-V(0,2)+V(1,2)-V(1,3)));	    
}
  
Kabbala VVVV_Calc::Triple_M1(const int &a,const int &b,
			     const int &c,const int &d,
			     const int &e)
{
  return 
    M(a)*M(b)*M(c)*X(a,a)*X(b,b)*X(c,c)*(X(d,e)*(V(a,e)*V(b,c)-V(a,c)*V(b,e))
					 +X(d,c)*(V(a,e)*V(b,d)-V(a,d)*V(b,e))
					 +(V(c,d)+V(c,e))*(X(d,a)*V(b,e)-X(d,b)*V(a,e))
					 );
  /*
    +M(0)*M(1)*M(2)*X(0,0)*X(1,1)*X(2,2)*(+X(3,4)*(V(0,4)*V(1,2)-V(0,2)*V(1,4))
					  +X(3,2)*(V(0,4)*V(1,3)-V(0,3)*V(1,4))
					  +(V(2,3)+V(2,4))*(X(3,0)*V(1,4)-X(3,1)*V(0,4))
					  ) //%8
    */      
}

Kabbala VVVV_Calc::Triple_M2(const int &a,const int &b,
			     const int &c,const int &d,
			     const int &e)
{
  return 
    M(a)*M(c)*M(e)*X(a,a)*X(c,c)*(X(d,c)*V(d,e)-X(d,e)*V(c,d))*(V(b,e)*X(b,a)-V(a,b)*X(b,e));
  //+M(0)*M(2)*M(4)*X(0,0)*X(2,2)*(X(3,2)*V(3,4)-X(3,4)*V(2,3))*(V(1,4)*X(1,0)-V(0,1)*X(1,4))
}

Kabbala VVVV_Calc::Double_M1(const int &a,const int &b,
			     const int &c,const int &d,
			     const int &e)
{
  return 
    M(a)*M(b)*X(a,a)*X(b,b)*( (X(d,b)*V(a,e)-X(d,a)*V(b,e))*(X(c,d)+X(c,e))
			     +(X(c,a)*V(b,e)-X(c,b)*V(a,e))*(X(d,c)+X(d,e))
			     +Z(d,c)*(+V(a,e)*(V(b,c)-V(b,d))
				      -V(b,e)*(V(a,c)-V(a,d))
				      )
			     );
  /*
    +M(0)*M(1)*X(0,0)*X(1,1)*(+(X(3,1)*V(0,4)-X(3,0)*V(1,4))*(X(2,3)+X(2,4))
                              -(X(2,1)*V(0,4)-X(2,0)*V(1,4))*(X(3,2)+X(3,4))
                              +Z(3,2)*(+V(0,4)*(V(1,2)-V(1,3))
                                       -V(1,4)*(V(0,2)-V(0,3))
    )
    )
  */
}

Kabbala VVVV_Calc::Double_M2(const int &a,const int &b,
			     const int &c,const int &d,
			     const int &e)
{
  return 
    M(b)*M(c)*X(b,b)*X(c,c)*((X(d,e)*X(a,c)+X(d,c)*X(a,d))*(V(b,e)-V(a,b))
			     +(X(a,e)*X(d,b)-X(d,a)*X(a,b))*(V(c,d)+V(c,e))   
			     +X(d,e)*(V(a,c)*X(a,b)-V(b,c)*X(a,e))
			     +X(d,c)*(V(a,d)*X(a,b)-V(b,d)*X(a,e))
			     +Z(d,a)*(V(c,e)+V(c,d))*(V(b,a)-V(b,e))
			     );
    /*
+M(1)*M(2)*X(1,1)*X(2,2)*(+(X(3,4)*X(0,2)+X(3,2)*X(0,3))*(V(1,4)-V(0,1))
			  +(X(0,4)*X(3,1)-X(3,0)*X(0,1))*(V(2,3)+V(2,4))   
			  +X(3,4)*(V(0,2)*X(0,1)-V(1,2)*X(0,4))
			  +X(3,2)*(V(0,3)*X(0,1)-V(1,3)*X(0,4))
			  +Z(3,0)*(V(2,4)+V(2,3))*(V(1,0)-V(1,4))
			  )
    */
    }

Kabbala VVVV_Calc::Double_M3(const int &a,const int &b,
			     const int &c,const int &d,
			     const int &e)
{
  return 
    +M(b)*M(e)*X(b,b)*(Z(d,c)*(V(a,e)*X(a,b)-V(a,b)*X(a,e))*(V(d,e)-V(c,e))
		       +(X(c,d)*X(d,e)-X(c,e)*X(d,c))*(V(a,b)*X(a,e)-V(a,e)*X(a,b))
		       );
    /*
      +M(1)*M(4)*X(1,1)*(+Z(3,2)*(V(0,4)*X(0,1)-V(0,1)*X(0,4))*(V(3,4)-V(2,4))
                         +(X(2,3)*X(3,4)-X(3,2)*X(2,4))*(V(0,1)*X(0,4)-V(0,4)*X(0,1))
      )
    */
}
  
Kabbala VVVV_Calc::Single_M(const int &a,const int &b,
			    const int &c,const int &d,
			    const int &e)
{
  return 
    M(a)*X(a,a)*( (X(d,a)*X(b,e)-X(d,b)*X(b,a))*(X(c,e)+X(c,d))
		 -(X(c,a)*X(b,e)-X(c,b)*X(b,a))*(X(d,e)+X(d,c))
		 +(V(a,b)-V(a,e))*(Z(d,b)*(X(c,d)+X(c,e))-Z(c,b)*(X(d,e)+X(d,c)))		   
		 +Z(d,c)*(X(b,a)*(V(b,d)-V(b,c))
			  -X(b,e)*(V(a,d)-V(a,c))
			  +(X(b,c)-X(b,d))*(V(a,b)-V(a,e))
			  )
		 );
  /*
+M(0)*X(0,0)*(+(X(3,0)*X(1,4)-X(3,1)*X(1,0))*(X(2,4)+X(2,3))
	      -(X(2,0)*X(1,4)-X(2,1)*X(1,0))*(X(3,4)+X(3,2))
	      +(V(0,1)-V(0,4))*(Z(3,1)*(X(2,3)+X(2,4))
	                       -Z(2,1)*(X(3,2)+X(3,4)))		   
	      +Z(3,2)*(+X(1,0)*(V(1,3)-V(1,2))
		       -X(1,4)*(V(0,3)-V(0,2))
		       +(X(1,2)-X(1,3))*(V(0,1)-V(0,4))
		       )
	      )
  */
}  
  











#include "AddOns/Analysis/Triggers/ConeMCFM.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

#define _VERBOSE 1

using namespace ATOOLS;
using namespace std;

ConeMCFM::ConeMCFM(double Rmin, double Rsep) : m_Rmin(Rmin), m_Rsep(Rsep) {}


double ConeMCFM::getet(Vec4D p1)
{
	return(sqr(p1[1]*p1[1]+p1[2]*p1[2])*p1[0]/sqr(p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3])); //MCFM definition
}


double ConeMCFM::deltarq(Vec4D p1, Vec4D p2)
{
	double pt1, pt2, dphi;
	double eta1, eta2;

	pt1=sqr(p1[1]*p1[1]+p1[2]*p1[2]);
	pt2=sqr(p2[1]*p2[1]+p2[2]*p2[2]);
	dphi=acos(Min(1.0,Max(-1.0,((p1[1]*p2[1]+p1[2]*p2[2])/pt1/pt2))));

	eta1=etarap(p1);
	eta2=etarap(p2);

	return(sqr(pow(eta1-eta2,2)+pow(dphi,2)));
}



double ConeMCFM::etarap(Vec4D p1)
{
	double energy;
	
	energy=sqr(Vec3D(p1).Abs());
	return(0.5*log(energy+p1[3])/(energy-p1[3]));
}


void ConeMCFM::ConstructJets(std::vector<Vec4D> & qjet) 
{
  	std::vector<Vec4D> qfinal;
  	std::vector<Vec4D> protoq;
  	
//	int mxpart(m_nout); 

  	int jets(0);
  	int maxjet;
  	int maxproto(0);

	//CHECK vector labels from here!
	std::vector<std::vector<int> > protoc;
  
	maxjet=qjet.size();

	switch(maxjet) {
		case 0:
			{
			jets=0;		
			};break;

// skip clustering if we only have one parton
		case 1:
			{
			qfinal.push_back(qjet[0]);
			jets=1;
			};break;
// set up proto-jets
		default:
			{
     			maxproto=-1;
      			for(int i=1; i<=maxjet; i++)
			{
        			maxproto++;
        			protoc[maxproto].push_back(1);
        			protoc[maxproto].push_back(i);
				protoq.push_back(qjet[i-1]);
      			};
			
			for(int i1=1; i1<=maxjet; i1++)
				{
				for (int i2=i1+1; i2<=maxjet; i2++)
					{
					maxproto++;
        				protoc[maxproto].push_back(2);
        				protoc[maxproto].push_back(i1);
        				protoc[maxproto].push_back(i2);
					protoq.push_back(qjet[i1-1]+qjet[i2-1]);
					
					if(  (deltarq(protoq[maxproto],protoq[i1-1])>m_Rmin)
					    ||(deltarq(protoq[maxproto],protoq[i2-1])>m_Rmin)
     				 	    ||(deltarq(protoq[i1-1],protoq[i2-1])>m_Rmin*m_Rsep))
						{
						protoq.erase(protoq.begin()+maxproto);
						maxproto=maxproto-1;
						};
      				
					}
				}
			if(maxjet > 2)
				{	
				for(int i1=1; i1<=maxjet; i1++)
					{
					for (int i2=i1+1; i2<=maxjet; i2++)
						{
						for (int i3=i2+1; i3<=maxjet; i3++)
							{
							maxproto++;
        						protoc[maxproto].push_back(3);
        						protoc[maxproto].push_back(i1);
        						protoc[maxproto].push_back(i2);
        						protoc[maxproto].push_back(i3);
							protoq.push_back(qjet[i1-1]+qjet[i2-1]+qjet[i3-1]);
						
							if(   (deltarq(protoq[maxproto],protoq[i1-1])>m_Rmin)
						    		||(deltarq(protoq[maxproto],protoq[i2-1])>m_Rmin)
						    		||(deltarq(protoq[maxproto],protoq[i3-1])>m_Rmin)
     					 	    		||(deltarq(protoq[i1-1],protoq[i2-1])>m_Rmin*m_Rsep)
     					 	    		||(deltarq(protoq[i1-1],protoq[i3-1])>m_Rmin*m_Rsep)
     					 	    		||(deltarq(protoq[i2-1],protoq[i3-1])>m_Rmin*m_Rsep))
								{
								protoq.erase(protoq.begin()+maxproto);
								maxproto=maxproto-1;
								};
      							}
						}
					}
				}
			if(maxjet > 3)
				{	
				for(int i1=1; i1<=maxjet; i1++)
					{
					for (int i2=i1+1; i2<=maxjet; i2++)
						{
						for (int i3=i2+1; i3<=maxjet; i3++)
							{
							for (int i4=i3+1; i4<=maxjet; i4++)
								{
								maxproto++;
        							protoc[maxproto].push_back(4);
        							protoc[maxproto].push_back(i1);
        							protoc[maxproto].push_back(i2);
        							protoc[maxproto].push_back(i3);
        							protoc[maxproto].push_back(i4);
								protoq.push_back(qjet[i1-1]+qjet[i2-1]+qjet[i3-1]+qjet[i4-1]);
						
								if(   (deltarq(protoq[maxproto],protoq[i1-1])>m_Rmin)
							    		||(deltarq(protoq[maxproto],protoq[i2-1])>m_Rmin)
							    		||(deltarq(protoq[maxproto],protoq[i3-1])>m_Rmin)
							    		||(deltarq(protoq[maxproto],protoq[i4-1])>m_Rmin)
     						 	    		||(deltarq(protoq[i1-1],protoq[i2-1])>m_Rmin*m_Rsep)
     						 	    		||(deltarq(protoq[i1-1],protoq[i3-1])>m_Rmin*m_Rsep)
     						 	    		||(deltarq(protoq[i1-1],protoq[i4-1])>m_Rmin*m_Rsep)
     						 	    		||(deltarq(protoq[i2-1],protoq[i3-1])>m_Rmin*m_Rsep)
     						 	    		||(deltarq(protoq[i2-1],protoq[i4-1])>m_Rmin*m_Rsep)
     						 	    		||(deltarq(protoq[i3-1],protoq[i4-1])>m_Rmin*m_Rsep))
									{
									protoq.erase(protoq.begin()+maxproto);
									maxproto=maxproto-1;
									}
      								}
							}
						}
					}
				}
      			if (maxjet>4) 
				{
				 msg_Error()<<"cannot do more than 4 jets in ConeMCFM jet algorithm"<<endl
					                           <<"Abort."<<endl;
				     abort();
				}
			
// loops through all iterations of the algorithm.
			if(maxproto>0)
				{
				MergeSplit(qjet, protoq, protoc, qfinal);
				jets=qfinal.size();
				}
			};break;
	}			
#if _VERBOSE
	                cout<< "Finished finding jets: got "<<jets <<endl;
#endif
			for(int i=0; i<jets; i++)	
				{
				m_pjets.push_back(qfinal[i]);
				}
			return;
}

void ConeMCFM::MergeSplit(std::vector<Vec4D> & qjet ,std::vector<Vec4D> & protoq, std::vector<std::vector<int> > & protoc, std::vector<Vec4D> & qfinal)
{
		int maxproto;
		maxproto=protoq.size()-1;

		int eti(-1);
		int ni(-1);
	while(maxproto>-1)
		{
// find highest Et proto-jet
		double etmax(-0.0);
		double et;
		for(int i=0; i<=maxproto; i++)
			{
			et=getet(protoq[i]);
			if(et>etmax) {eti=i; etmax=et;};
			}

// check to see if any partons are shared by this proto-jet
		int shared(0);
		double sharedet(0.0);
		vector<int> sharedc;
		for(int i=0;i<=maxproto; i++)
			{
			sharedc[i]=0;
			if(i!=eti)
				{
				for(int j=1;j<=protoq[i][0]; j++)
					{
					for(int k=1;k<=protoq[eti][0]; k++)
						{
						if(protoc[i][j]==protoc[eti][k]) 
							{
							shared++;
							sharedc[i]=1;
							}
						}
					}
				}
			}
		switch(shared)
			{
			case 0:	
			{
			qfinal.push_back(protoq[eti]);
		//shuffle down proto-jets
			maxproto=maxproto-1;
			protoc.erase(protoc.begin()+eti);
			protoq.erase(protoq.begin()+eti);		
			};break;
			default:
			{
// a parton is shared: perform split/merge procedure
		double net(0.0);

		for(int i=0; i<=maxproto;i++)
			{
			et=getet(protoq[i]);
			if(et>net && sharedc[i]==1)
				{
				ni=i;
				net=et;
				};
			}
// calculate the shared Et
		Vec4D qshared;
		qshared=Vec4D(0.,0.,0.,0.);
		
		for(int i=1; i<=protoc[eti][0]; i++)
			{
			for(int j=1; j<=protoc[ni][0];j++)
				{
				if(protoc[eti][i]==protoc[ni][j])
					{
					qshared=qshared+qjet[protoc[eti][i]];
					}
				}
			}
		sharedet=getet(qshared);

#if _VERBOSE
		cout<< "Proto-jet is " << eti <<endl;
		cout<< "Highest Et neighbor is " << ni <<endl;
		cout<< "Shared Et is " << sharedet <<endl;
		cout<< "Neighbor Et is " << net <<endl;
#endif

// condition for merge of jets
		if(sharedet/net > 0.5)
			{
			for(int i=1;i<=protoc[ni][0];i++)
				{
				shared=0;
				for(int j=1; j<=protoc[eti][0];j++)
					{
					if(protoc[ni][i]==protoc[eti][j]) {shared=1;};
					}
//add cells that are not shared
				if(shared==0)
					{
					protoc[eti][0]=protoc[eti][0]+1;
					protoc[eti][protoc[eti][0]]=protoc[ni][i];
					protoq[eti]=protoq[eti]+qjet[protoc[ni][i]];
					}
				}
//shuffle down the proto-jet data
			protoq.erase(protoq.begin()+ni);
			protoc.erase(protoc.begin()+ni);	
#if _VERBOSE
			cout<< "merged proto-jets "<< eti <<" and "<< ni <<endl;
#endif	

			}
		else
			{
//we should split the proto-jets
			for(int i=1;i<=protoc[ni][0];i++)
                                {
	                        shared=0;
	                        for(int j=1; j<=protoc[eti][0];j++)
		                     {
			             if(protoc[ni][i]==protoc[eti][j]) {shared=j;};
			             }
//if cell is shared decide where to put it bases on distance in Delta_R
//update the contents list and momentum of the protojet it's removed from
				if(shared >0)
					{
					if (deltarq(qjet[protoc[ni][i]-1],protoq[ni])<deltarq(qjet[protoc[ni][i]-1],protoq[eti]))
						{
// shared cell closer to neighbor
						if(protoc[eti][0]==1)
							{
							protoc.erase(protoc.begin()+eti);
							protoq.erase(protoq.begin()+eti);
							}
						else
							{
						protoc[eti].erase(protoc[eti].begin()+shared);
						protoc[eti][0]=protoc[eti][0]-1;
						protoq[eti]=protoq[eti]-qjet[protoc[ni][i]-1];	
							}
						
						}	
					else	
						{
// shared cell closer to original proto-jet
						if(protoc[ni][0]==1)
							{	
							protoc.erase(protoc.begin()+ni);
							protoq.erase(protoq.begin()+ni);
							}
						else	
							{
							protoc[ni][0]=protoc[ni][0]-1;
							protoq[ni]=protoq[ni]-qjet[protoc[ni][i]-1];
							protoc[ni].erase(protoc[ni].begin()+i);
							}
						}
					}
				}
#if _VERBOSE
			cout<< "split proto-jets "<< eti <<" and "<< ni <<endl;
#endif	
			}
			};break;
			}
// iterate merging and splitting until all protojets are gone
		maxproto=protoq.size()-1;
		}
	return;
}

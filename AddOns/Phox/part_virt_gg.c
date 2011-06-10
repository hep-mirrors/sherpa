// q(p1) + qb(p2) --> gamma(p3) + gamma(p4) 
// with p1+p2=p3+p4 
// sc = 2*p1.p2, tc=-2*p1.p3, uc=-2*p2.p3 
// eqi : electric charge of quark in unit of e
// q2 : arbitrary scale > 0 
// N : nb of colour of quarks 
// 
// Normalisation
// 1/(4 N^2) \alpha^2 \alpha_s/sc^2 \Gamma(1+\epsilon) \Gamma(1-\epsilon)^2/\Gamma(1-2 \epsilon) (4 \pi \mu^2/q2)^{\epsilon} 
// 
// tc2 :coefficient of 1/epsilon^2 
// tc1 : coefficient of 1/epsilon 
// tc0 : coefficient O(epsilon^0) 
      tc2 = -8.0*(N+1.0)*(N-1.0)*eqi*eqi*eqi*eqi*(uc*uc+tc*tc)/tc/uc;
      tc1 = 4.0*eqi*eqi*eqi*eqi*(N-1.0)*(N+1.0)*(2.0*sc*sc-uc*uc-tc*tc+2.0*log(
sc/q2)*uc*uc+2.0*log(sc/q2)*tc*tc)/tc/uc;
      tc0 = 4.0*eqi*eqi*eqi*eqi*(-2.0*log(sc/q2)*sc*sc-2.0*log(sc/q2)*tc*tc-2.0
*log(sc/q2)*uc*uc+2.0*sc*sc*pow(log(sc/q2),2.0)+pow(log(-tc/q2),2.0)*tc*tc+pow(
log(-tc/q2),2.0)*sc*sc+3.0*uc*uc*log(-tc/q2)+0.3141592653589793E1*
0.3141592653589793E1*tc*tc+0.3141592653589793E1*0.3141592653589793E1*uc*uc+pow(
log(-uc/q2),2.0)*uc*uc+pow(log(-uc/q2),2.0)*sc*sc+3.0*tc*tc*log(-uc/q2)-2.0*log
(sc/q2)*log(-tc/q2)*tc*tc-2.0*log(sc/q2)*log(-tc/q2)*sc*sc-4.0*tc*uc*log(sc/q2)
+2.0*uc*log(-tc/q2)*tc-2.0*log(sc/q2)*log(-uc/q2)*uc*uc-2.0*log(sc/q2)*log(-uc/
q2)*sc*sc+2.0*tc*log(-uc/q2)*uc-4.0*tc*tc-4.0*uc*uc+sc*sc)*(N+1.0)*(N-1.0)/tc/
uc;

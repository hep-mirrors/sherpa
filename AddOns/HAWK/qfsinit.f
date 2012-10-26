************************************************************************
        subroutine qfsinit(inqu,inqd,inql,inqn,inmu,inmd,xsw,xcw)
************************************************************************
        implicit real*8 (a-z)
        complex*16 xsw,xcw,guu,gdd            
        common/qf/qu,qd,ql,qn,qf(4),guu(-1:1),gdd(-1:1),mu,mu2,md,md2

        qu  = inqu
        qd  = inqd
        qf(1) = qu
        qf(2) = qd
        qf(3) = qd
        qf(4) = qu
        ql  = inql
        qn  = inqn
        mu  = inmu
        mu2 = mu*mu
        md  = inmd
        md2 = md*md
        guu(-1) = -qu*xsw/xcw; 
        guu(1)  = -qu*xsw/xcw+0.5d0/(xsw*xcw); 
        gdd(-1) = -qd*xsw/xcw; 
        gdd(1)  = -qd*xsw/xcw-0.5d0/(xsw*xcw); 
        guu(0)  = (0d0,0d0)
        gdd(0)  = (0d0,0d0)

        write (*,*) '@@@ Test qf  = (',qf(1),qf(2),qf(3),qf(4),')'
        write (*,*) '@@@ Test guu = (',guu(-1),guu(0),guu(1),')'
        write (*,*) '@@@ Test gdd = (',gdd(-1),gdd(0),gdd(1),')'


        end

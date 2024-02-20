!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                             int2e.f90                          !!!!
!!!!        This file calculats the general two electron intgral.   !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!1. ERINT
!description: subroutine
!  Calculat the general overlap intgral, Kinetic intgral and Nuclear 
!   attraction  intgral  of shell.    
!   input:   (1) katom        integer(natom)
!               Z of atoms
!            (2)coor        double precision(natom,3)
!               (coordinates of atom )
!            (3) ntype     integer(nshell)
!             (ntype  is  type of gauss function for all shells)
!               (S=1,P=2,D=3,F=4,G=5,H=6,I=7)
!            (4) ngf     integer(nshell)
!             (The gauss function number of all shells)
!            (5) gfc    double precision(ngfsum)
!             (gfc is coefficient  of primitive gauss functions)
!               ngfsum=SUM(ngf(:))
!            (6) gfe    double precision(ngfsum)
!             (gfc is orbital exponent  of primitive gauss functions)
!   output:  (1) S     double precision
!            (2) T     double precision  
!            (3) V     double precision  
!2. recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,
!     ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,Px,Qx,Ay,By,Cy,Dy,Py,
!              Qy,Az,Bz,Cz,Dz,Pz,Qz,T)  result(tem)
!description: recursive function
!  Calculat the ERI <ac|1/r12|bd>(m).
!  use DPK method
!  ref:F W.Chen,Computational Methods in Quantum Chemistry
!   input:   (1) m   integer
!               (m=0,<ac|1/r12|bd>)
!            (2) nax,nay,naz      integer
!            (3) nbx,nby,nbz      integer
!            (4) ncx,ncy,ncz      integer
!            (5) ndx,ndy,ndz      integer
!             (2),(3),(4),(5) gauss function parameter
!            (6) ep             double precision
!            (7) eq             double precision
!            (8) ABK              double precision
!                 ABK=K_{AB}  
!            (9) CDK              double precision
!                 CDK=K_{CD}  
!            (10) Ax,Ay,Az         double precision
!            (11) Bx,By,Bz         double precision
!            (12) Cx,Cy,Cz         double precision
!            (13) Dx,Dy,Dz         double precision
!          (10),(11),(12),(13) Coordinates of center a,b,c,d
!            (14) Px,Py,Pz         double precision
!                P is the new center after combining A,B.
!            (15) Qx,Qy,Qz         double precision
!                Q is the new center after combining C,D.
!            (16) T               double precision
!                 T=(ea*eb)/(ea+eb)*((Px-Qx)^2+(Py-Qy)^2+(Pz-Qz)^2)
!   output:  (1) result           double precision 
!3. recERI2(m,nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,
!     ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,Px,Qx,Ay,By,Cy,Dy,Py,
!              Qy,Az,Bz,Cz,Dz,Pz,Qz,T)  result(tem)
!description: recursive function
!  Calculat the ERI <ac|1/r12|bd>(m).
!  use S.Obara;A.saika method
!  ref:F W.Chen,Computational Methods in Quantum Chemistry
!   input:   (1) m   integer
!               (m=0,<ac|1/r12|bd>)
!            (2) nax,nay,naz      integer
!            (3) nbx,nby,nbz      integer
!            (4) ncx,ncy,ncz      integer
!            (5) ndx,ndy,ndz      integer
!             (2),(3),(4),(5) gauss function parameter
!            (6) ep             double precision
!            (7) eq             double precision
!            (8) ABK              double precision
!                 ABK=K_{AB}  
!            (9) CDK              double precision
!                 CDK=K_{CD}  
!            (10) Ax,Ay,Az         double precision
!            (11) Bx,By,Bz         double precision
!            (12) Cx,Cy,Cz         double precision
!            (13) Dx,Dy,Dz         double precision
!          (10),(11),(12),(13) Coordinates of center a,b,c,d
!            (14) Px,Py,Pz         double precision
!                P is the new center after combining A,B.
!            (15) Qx,Qy,Qz         double precision
!                Q is the new center after combining C,D.
!            (16) T               double precision
!                 T=(ea*eb)/(ea+eb)*((Px-Qx)^2+(Py-Qy)^2+(Pz-Qz)^2)
!   output:  (1) result           double precision 
!4. anaERI1(nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,
!     ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,Px,Qx,Ay,By,Cy,Dy,Py,
!              Qy,Az,Bz,Cz,Dz,Pz,Qz,T)  result(tem)
!description: function
!  Calculat the ERI <ac|1/r12|bd>(m).
!  use (Analytical method) O.Ohata method 
!  ref:F W.Chen,Computational Methods in Quantum Chemistry
!   input:  
!            (1) nax,nay,naz      integer
!            (2) nbx,nby,nbz      integer
!            (3) ncx,ncy,ncz      integer
!            (4) ndx,ndy,ndz      integer
!             (1),(2),(3),(4) gauss function parameter
!            (5) ep             double precision
!            (6) eq             double precision
!            (7) ABK              double precision
!                 ABK=K_{AB}  
!            (8) CDK              double precision
!                 CDK=K_{CD}  
!            (9) Ax,Ay,Az         double precision
!            (10) Bx,By,Bz         double precision
!            (11) Cx,Cy,Cz         double precision
!            (12) Dx,Dy,Dz         double precision
!          (9),(10),(11),(12) Coordinates of center a,b,c,d
!            (13) Px,Py,Pz         double precision
!                P is the new center after combining A,B.
!            (14) Qx,Qy,Qz         double precision
!                Q is the new center after combining C,D.
!            (15) T               double precision
!                 T=(ea*eb)/(ea+eb)*((Px-Qx)^2+(Py-Qy)^2+(Pz-Qz)^2)
!   output:  (1) result           double precision 
!5. BFanaERI1
!description:  function
!  BFanaERI1 is a part of anaERI1.
!6. recERIfast
!description:  function
!   recERIfast is the efficient version of recERI1.
!7. recERI1xyz
!description:  function
!  recERI1xyz is a part of recERIfast.
!======================================================================!
  Module Int2e
    contains
    subroutine ERINT(nshell,nbasis,natom,coor,neri,&
                    kcenter,gfc,gfe,ntype,ngf,ngfall,&
                    kbasst,kbased,knbas,kpgfst,&
                       kpgfed,naxyz,ERIall,methodERI)
        use Global
        implicit double precision (A-H,O-Z)
        Integer*16 neri 
        dimension::kbasst(nshell),kbased(nshell),knbas(nshell),&
                 kpgfst(nshell),kpgfed(nshell),naxyz(3*nbasis),&
                 ERIall(neri),coor(natom,3)   
        dimension::kcenter(Natom),gfc(ngfall),gfe(ngfall),ntype(nshell),ngf(nshell)     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!=================================================================!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  SHELL  !!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        write(*,*)'*** Computing two electron intgral ***'
        write(*,*)'neri=',neri   !! number of stored ERI
        nsize=neri*8
        write(*,*)'The storage space of storing ERI  requires approximately ' 
        write(*,*) 'size',nsize,'B',nsize/1024**1d0,'Kb',nsize/1024**2d0,'Mb',nsize/1024**3d0,'Gb'          
        keri=1
        do Ish=1,nshell   ! I shell
            Iatom=kcenter(Ish)
            cix=coor(Iatom,1)
            ciy=coor(Iatom,2)
            ciz=coor(Iatom,3)
            do Jsh=1,Ish  ! J shell  
            ! do Jsh=1,nshell  ! J shell 
                Jatom=kcenter(Jsh)
                cJx=coor(Jatom,1)
                cJy=coor(Jatom,2)
                cJz=coor(Jatom,3)
                Rij=(cix-cjx)**2+(ciy-cjy)**2+(ciz-cjz)**2
                do ksh=1,Ish   ! K shell
                    kkatom=kcenter(Ksh)
                    ckx=coor(Kkatom,1)
                    cky=coor(Kkatom,2)
                    ckz=coor(Kkatom,3)
                    do Lsh=1,ksh  ! L shell  
                        latom=kcenter(lsh)
                        clx=coor(latom,1)
                        cly=coor(latom,2)
                        clz=coor(latom,3)
                        Rkl=(ckx-clx)**2+(cky-cly)**2+(ckz-clz)**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! BASIS FUNCTION  !!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            do ibas=kbasst(ish),kbased(ish)
                nax=naxyz(ibas)
                nay=naxyz(ibas+nbasis)
                naz=naxyz(ibas+2*nbasis)
                jbased=kbased(jsh) 
                ipgfst=kpgfst(ish)
                ipgfed=kpgfed(ish)
!!!!! only for sp shell                    
                if(ntype(ish)==12)then
                    if(nax+nay+naz==0) then
                        ipgfed=kpgfed(ish)-ngf(ish)/2
                    else  
                        ipgfst=kpgfst(ish)+ngf(ish)/2
                    end if
                end if                                        
                do jbas=kbasst(jsh),jbased
                    if(ibas<jbas) cycle
                    nbx=naxyz(jbas)
                    nby=naxyz(jbas+nbasis)
                    nbz=naxyz(jbas+2*nbasis)
                    jpgfst=kpgfst(jsh)
                    jpgfed=kpgfed(jsh)                        
!!!!! only for sp shell                        
                    if(ntype(jsh)==12)then
                        if(nbx+nby+nbz==0)then
                            jpgfed=kpgfed(jsh)-ngf(jsh)/2
                        else  
                            jpgfst=kpgfst(jsh)+ngf(jsh)/2
                        end if
                    end if                      
                    do kbas=kbasst(ksh),kbased(ksh)
                        ncx=naxyz(kbas)
                        ncy=naxyz(kbas+nbasis)
                        ncz=naxyz(kbas+2*nbasis)
                        lbased=kbased(lsh)      
                        kkpgfst=kpgfst(ksh)
                        kkpgfed=kpgfed(ksh)
!!!!! only for sp shell                    
                        if(ntype(ksh)==12)then
                            if(ncx+ncy+ncz==0) then
                                kkpgfed=kpgfed(ksh)-ngf(ksh)/2
                            else  
                                kkpgfst=kpgfst(ksh)+ngf(ksh)/2
                            end if
                        end if                                           
                        do lbas=kbasst(lsh),lbased
                            ndx=naxyz(lbas)
                            ndy=naxyz(lbas+nbasis)
                            ndz=naxyz(lbas+2*nbasis) 
                            lpgfst=kpgfst(lsh)
                            lpgfed=kpgfed(lsh)
!!!!! only for sp shell                    
                            if(ntype(lsh)==12)then
                                if(ndx+ndy+ndz==0) then
                                    lpgfed=kpgfed(lsh)-ngf(lsh)/2
                                else  
                                    lpgfst=kpgfst(lsh)+ngf(lsh)/2
                                end if
                            end if                              
                            ij=ibas*(ibas+1)/2+Jbas
                            kl=kbas*(kbas+1)/2+lbas  
                            if(ibas<jbas .or.kbas<lbas .or.ij<kl) cycle
                        eri=0d0 
                        eri1=0d0
                        eri2=0d0
                        eri3=0d0
                        eri4=0d0
!!!!!!!!!!!!!!!  <ac|bd>=(ab|cd)
!!! primitive gauss functions  (a|                   
                        do ipgf=ipgfst,ipgfed     
                                ea=gfe(ipgf)
                                ca=gfc(ipgf)
                            dNa=(2*ea/pi)**0.75*((4*ea)**(nax+nay+naz)&
                                /(Dfact(2*nax-1)*Dfact(2*nay-1)*Dfact(2*naz-1)))**0.5
!!! primitive gauss functions (b|
                                do jpgf=jpgfst,jpgfed 
                                    eb=gfe(jpgf)
                                    cb=gfc(jpgf)
                                    ep=ea+eb
                                    ep1=1d0/ep
                                    ABK=exp(-ea*eb*ep1*Rij)
                                    PX=ep1*(ea*cix+eb*cjx)
                                    PY=ep1*(ea*ciy+eb*cjy)
                                    PZ=ep1*(ea*ciz+eb*cjz)                                    
                                    if(abs(Rij)<1D-15) then
                                        ABK=1d0
                                        PX=cix
                                        PY=ciy
                                        PZ=ciz                                    
                                    end if   
                    dNb=(2*eb/pi)**0.75*((4*eb)**(nbx+nby+nbz)&
                            /(Dfact(2*nbx-1)*Dfact(2*nby-1)*Dfact(2*nbz-1)))**0.5                                       
!!! primitive gauss functions  c|)                   
                        do kpgf=kkpgfst,kkpgfed    
                                ec=gfe(kpgf)
                                cc=gfc(kpgf)
                        dNc=(2*ec/pi)**0.75*((4*ec)**(ncx+ncy+ncz)&
                            /(Dfact(2*ncx-1)*Dfact(2*ncy-1)*Dfact(2*ncz-1)))**0.5 
!!! primitive gauss functions |d)
                                do lpgf=lpgfst,lpgfed
                                    ed=gfe(lpgf)
                                    cd=gfc(lpgf)
                                    eq=ec+ed
                                    eq1=1d0/eq
                                    CDK=exp(-ec*ed*eq1*Rkl)
                                    QX=eq1*(ec*ckx+ed*clx)
                                    QY=eq1*(ec*cky+ed*cly)
                                    QZ=eq1*(ec*ckz+ed*clz)                                    
                                    if(abs(Rkl)<1D-15) then
                                        CDK=1d0
                                        QX=ckx
                                        QY=cky
                                        QZ=ckz                                    
                                    end if   
                                    rho=ep*eq/(ep+eq)
                                    PQ=(Px-Qx)**2+(Py-Qy)**2+(Pz-Qz)**2
                                    TT=rho*PQ
                                dNd=(2*ed/pi)**0.75*((4*ed)**(ndx+ndy+ndz)&
                                    /(Dfact(2*ndx-1)*Dfact(2*ndy-1)*Dfact(2*ndz-1)))**0.5 
                                naxx=nax
                                nbxx=nbx
                                ncxx=ncx
                                ndxx=ndx
                                Axx=cix
                                Bxx=cjx  
                                Cxx=ckx
                                Dxx=clx 
                                pxx=Px
                                qxx=Qx                 
                                nayy=nay
                                nbyy=nby
                                ncyy=ncy
                                ndyy=ndy
                                Ayy=ciy
                                Byy=cjy  
                                Cyy=cky
                                Dyy=cly 
                                pyy=Py
                                qyy=Qy                                  
                                nazz=naz
                                nbzz=nbz
                                nczz=ncz
                                ndzz=ndz
                                Azz=ciz
                                Bzz=cjz  
                                Czz=ckz
                                Dzz=clz 
                                pzz=Pz
                                qzz=Qz
                                epp=ep
                                eqq=eq                                 
!!!!!!!
                        if(methodERI==1)then
!!!!!!!!! recursive method 1 
!!!!!!! Dupuis-Rys-King method
                        
                        eri=eri+ca*cb*cc*cd*dNa*dNb*dNc*dNd*recERI1(0,naxx,nayy,nazz,&
                                nbxx,nbyy,nbzz,ncxx,ncyy,nczz,ndxx,ndyy,ndzz,&
                                epp,eqq,ABK,CDK,Axx,Bxx,Cxx,Dxx,Pxx,Qxx,Ayy,Byy,Cyy,Dyy,Pyy,&
                                Qyy,Azz,Bzz,Czz,Dzz,Pzz,Qzz,TT) 

                        else if(methodERI==2)then
!!!!!!!!!!!!! Recursive method 2   
!!!!!!! S.Obara;A.saika  The Journal of Chemical Physics 84, 3963 (1986);
!!!!! doi: 10.1063/1.450106                        
                        methodERI=2                        
                        eri=eri+ca*cb*cc*cd*dNa*dNb*dNc*dNd*recERI2(0,naxx,nayy,nazz,&
                                nbxx,nbyy,nbzz,ncxx,ncyy,nczz,ndxx,ndyy,ndzz,&
                                epp,eqq,ABK,CDK,Axx,Bxx,Cxx,Dxx,Pxx,Qxx,Ayy,Byy,Cyy,Dyy,Pyy,&
                                Qyy,Azz,Bzz,Czz,Dzz,Pzz,Qzz,TT)     

                        else if(methodERI==3)then                        
!!!!!!! 3.(Analytical method) O.Ohata method (1966)
!!!! https://doi.org/10.1143/JPSJ.21.2313                          
                        methodERI=3
                   eri=eri+ca*cb*cc*cd*dNa*dNb*dNc*dNd&
                         *anaERI1(naxx,nayy,nazz,nbxx,nbyy,nbzz,ncxx,ncyy,nczz,ndxx,ndyy,ndzz,&
                         epp,eqq,ABK,CDK,Axx,Bxx,Cxx,Dxx,Pxx,Qxx,Ayy,Byy,Cyy,Dyy,Pyy,&
                         Qyy,Azz,Bzz,Czz,Dzz,Pzz,Qzz,TT)  
                        else if(methodERI==4)then    
!!!!! fast method                                                                 
                        methodERI=4
                   eri=eri+ca*cb*cc*cd*dNa*dNb*dNc*dNd&
                        *recERIfast(naxx,nayy,nazz,nbxx,nbyy,nbzz,ncxx,ncyy,nczz,ndxx,ndyy,ndzz,&
                         epp,eqq,ABK,CDK,Axx,Bxx,Cxx,Dxx,Pxx,Qxx,Ayy,Byy,Cyy,Dyy,Pyy,&
                         Qyy,Azz,Bzz,Czz,Dzz,Pzz,Qzz,TT)                        
!!!!                          
                        else if(methodERI==-1)then
!!!! for test                                           
                        eri1=eri1+ca*cb*cc*cd*dNa*dNb*dNc*dNd*recERI1(0,naxx,nayy,nazz,&
                                nbxx,nbyy,nbzz,ncxx,ncyy,nczz,ndxx,ndyy,ndzz,&
                                epp,eqq,ABK,CDK,Axx,Bxx,Cxx,Dxx,Pxx,Qxx,Ayy,Byy,Cyy,Dyy,Pyy,&
                                Qyy,Azz,Bzz,Czz,Dzz,Pzz,Qzz,TT) 
!!!! for test                                           
                        eri2=eri2+ca*cb*cc*cd*dNa*dNb*dNc*dNd*recERI2(0,naxx,nayy,nazz,&
                                nbxx,nbyy,nbzz,ncxx,ncyy,nczz,ndxx,ndyy,ndzz,&
                                epp,eqq,ABK,CDK,Axx,Bxx,Cxx,Dxx,Pxx,Qxx,Ayy,Byy,Cyy,Dyy,Pyy,&
                                Qyy,Azz,Bzz,Czz,Dzz,Pzz,Qzz,TT)      
!!!! for test                                           
                        eri3=eri3+ca*cb*cc*cd*dNa*dNb*dNc*dNd&
                         *anaERI1(naxx,nayy,nazz,nbxx,nbyy,nbzz,ncxx,ncyy,nczz,ndxx,ndyy,ndzz,&
                         epp,eqq,ABK,CDK,Axx,Bxx,Cxx,Dxx,Pxx,Qxx,Ayy,Byy,Cyy,Dyy,Pyy,&
                         Qyy,Azz,Bzz,Czz,Dzz,Pzz,Qzz,TT) 
!!!! for test                            
                        eri4=eri4+ca*cb*cc*cd*dNa*dNb*dNc*dNd&
                         *recERIfast(naxx,nayy,nazz,nbxx,nbyy,nbzz,ncxx,ncyy,nczz,ndxx,ndyy,ndzz,&
                         epp,eqq,ABK,CDK,Axx,Bxx,Cxx,Dxx,Pxx,Qxx,Ayy,Byy,Cyy,Dyy,Pyy,&
                         Qyy,Azz,Bzz,Czz,Dzz,Pzz,Qzz,TT)                                                                        
                        else 
                            stop "bad ERI method"
                        end if  
                        end do  !!!  d gf
                    end do  !!!  c gf                   
                end do  !!! b gf                 
            end do  !!! a gf 
            ! if(abs(eri)>=1d-9)write(99,*)'I=',ibas,'J=',Jbas,'K=',Kbas,'L=',Lbas,'ERI=',eri
            ! if(kprinteri==2)then
            !     write(99,*)'I=',ibas,'J=',Jbas,'K=',Kbas,'L=',Lbas,'ERI=',eri
            ! end if
                    ijkl=indexeri(ibas,jbas,kbas,lbas)
                    if(abs(eri)<1d-16)eri=0d0
                    ERIall(ijkl)=eri
!!! test 
                    if  (methodERI==-1)then
                        msum=nax+nay+naz+nbx+nby+nbz+ncx+ncy+ncz+ndx+ndy+ndz
                        if (abs(eri1-eri2)>1d-16 .or.abs(eri1-eri3)>1d-16 .or.abs(eri2-eri3)>1d-16&
                                    .OR. abs(eri1-eri4)>1d-16)then
                            write(99,*)'I=',ibas,'J=',Jbas,'K=',Kbas,'L=',Lbas,'delt ERI12=',eri1-eri2,&
                            'delt ERI13=',eri1-eri3,'delt ERI23=',eri2-eri3,'delt ERI14=',eri1-eri4,&
                            'TT=',TT,'m=',msum,'ERI1=',eri1,'ERI2=',eri2,'ERI3=',eri3,'ERI4=',eri4,&
                            'nx1=',nax+nbx,'nx2=',ncx+ndx,'ny1=',nay+nby,'ny2=',ncy+ndy,'nz1=',naz+nbz,&
                            'nz2=',ncz+ndz,nax,nbx,ncx,ndx,nay,nby,ncy,ndy,naz,nbz,ncz,ndz
                        end if
                    end if
!!!!!!                    
                    keri=keri+1 
                    if(mod(keri,10000)==0)write(*,*)keri,'ERI done,',real(keri)/neri*100d0,'% done'
                        end do  !!!  d basis
                    end do  !!!  c basis                   
                end do  !!! b basis                 
            end do  !!! a basis                    
!!!!!!!!!!!!!!!                                                
                    end do  !!!  L shell 
                end do  !!!  K shell                   
            end do  !!!  J shell                  
        end do  !!!  I shell      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!! now in IO.f90                                   
        ! if(methodERI==1)then 
        !     write(99,'(23x,A40)')'ERI: Dupuis-Rys-King method '
        ! else if(methodERI==2) then
        !     write(99,'(23x,A40)')'ERI: S.Obara Recursive Method '
        ! else if(methodERI==3) then
        !     write(99,'(23x,A40)')'ERI: (Analytical) O.Ohata method '        
        ! else
        !     write(99,'(23x,A40)')'ERI: Something wrong !!! '
        ! end if
!!!!!!!      
    write(*,*)'*** Two electron intgral Done ***'
!!!!
    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive function recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)  result(tem)
!!!!!                                
    use Global,only:PI,F,F_pade,F_ser
    implicit double precision (A-H,O-Z)    
    tem=0d0
    if(nax<0.or.nbx<0.or.nay<0.or.nby<0.or.naz<0.or.nbz<0.or.&
        ncx<0.or.ndx<0.or.ncy<0.or.ndy<0.or.ncz<0.or.ndz<0)then
        tem=0d0
        return
    end if
    epq1=1d0/(ep+eq)
!!!!!!!!! recursive method 1 
!!!!!!! Dupuis-Rys-King method
!!!! F W.Chen,Computational Methods in Quantum Chemistry
        if (nax==0.and.nbx==0.and.ncx==0.and.ndx==0)then
            if (nay==0.and.nby==0.and.ncy==0.and.ndy==0)then
                if (naz==0.and.nbz==0.and.ncz==0.and.ndz==0)then
                    tem=2d0*ABK*CDK*pi**2.5d0/(ep*eq)/sqrt(ep+eq)*F(m,TT)
                    ! tem=2*ABK*CDK/(ep*eq)*pi**2.5/sqrt(ep+eq)*F_pade(m,TT)
                    ! tem=2d0*ABK*CDK/(ep*eq)*pi**2.5/sqrt(ep+eq)*F_ser(m,TT)
                    return
!!!! z                    
                else if(ndz==0) then            !!!! (ab|c0)z
                    if(nbz==0) then             !!!! (a0|c0)z
                        if(ncz==0) then         !!!! (a0|00)z
                            tem=0.5d0*(naz-1)/ep*(&
                            recERI1(m,nax,nay,naz-2,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                            -eq*epq1*&
                            recERI1(m+1,nax,nay,naz-2,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) )&
                            +(Ppz-Az)*recERI1(m,nax,nay,naz-1,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                            +eq*(QQz-PPz)*epq1*recERI1(m+1,nax,nay,naz-1,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                            +0.5d0*ncz*epq1* recERI1(m+1,nax,nay,naz-1,nbx,nby,nbz,ncx,ncy,ncz-1,ndx,ndy,ndz,&
                                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)                                                              
                            return
                        else                !!!! (ab|00)z  nbz/=0
                        tem=0.5d0*(ncz-1)/eq*(&
                        recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz-2,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                        -ep*epq1*&
                        recERI1(m+1,nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz-2,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) )&
                        +(QQz-Cz)*recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz-1,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                        +ep*(PPz-QQz)*epq1*recERI1(m+1,nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz-1,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                        +0.5d0*naz*epq1* recERI1(m+1,nax,nay,naz-1,nbx,nby,nbz,ncx,ncy,ncz-1,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)                           
                        return
                        end if
                    else        !!! ncz/=0  !!!! (ab|c0)z
                    tem=recERI1(m,nax,nay,naz+1,nbx,nby,nbz-1,ncx,ncy,ncz,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                        +(Az-Bz)*recERI1(m,nax,nay,naz,nbx,nby,nbz-1,ncx,ncy,ncz,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)                                                   
                        return
                    end if
                else           !!!! ndz /= 0       (ab|cd)z
                    tem=recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz+1,ndx,ndy,ndz-1,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                        +(Cz-Dz)*recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz-1,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) 
                    return
                end if    
!!!! y                    
            else if(ndy==0) then                !!!! (ab|c0)y
                    if(nby==0) then         !!!! (ab|00)y
                        if(ncy==0) then     !!!! (a0|00)y
                         tem=0.5d0*(nay-1)/ep*(&
                            recERI1(m,nax,nay-2,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                            -eq*epq1*&
                            recERI1(m+1,nax,nay-2,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) )&
                            +(Ppy-Ay)*recERI1(m,nax,nay-1,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                            +eq*(QQy-PPy)*epq1*recERI1(m+1,nax,nay-1,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                            +0.5d0*ncy*epq1* recERI1(m+1,nax,nay-1,naz,nbx,nby,nbz,ncx,ncy-1,ncz,ndx,ndy,ndz,&
                                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)                                                              
                            return
                        else            !!!! (ab|00)y  nby/=0
                        tem=0.5d0*(ncy-1)/eq*(&
                        recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx,ncy-2,ncz,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                        -ep*epq1*&
                        recERI1(m+1,nax,nay,naz,nbx,nby,nbz,ncx,ncy-2,ncz,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) )&
                        +(QQy-Cy)*recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx,ncy-1,ncz,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                        +ep*(PPy-QQy)*epq1*recERI1(m+1,nax,nay,naz,nbx,nby,nbz,ncx,ncy-1,ncz,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                        +0.5d0*nay*epq1* recERI1(m+1,nax,nay-1,naz,nbx,nby,nbz,ncx,ncy-1,ncz,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)    
                    return
                        end if
                    else            !!! ncy/=0  !!!! (ab|c0)y

                        tem=recERI1(m,nax,nay+1,naz,nbx,nby-1,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                        +(Ay-By)*recERI1(m,nax,nay,naz,nbx,nby-1,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)                                                                                   
                        return
                    end if
                else           !!!! ndy /= 0       (ab|cd)y
                    tem=recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx,ncy+1,ncz,ndx,ndy-1,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                        +(Cy-Dy)*recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy-1,ndz,&
                                        ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                        QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) 
                    return
                end if                          
!!!! x                    
        else if(ndx==0) then       !!!! (ab|c0)x
            if(nbx==0) then        !!!! (ab|00)x
                if(ncx==0) then    !!!! (a0|00)x
                 tem=0.5d0*(nax-1)/ep*(&
                    recERI1(m,nax-2,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                    -eq*epq1*&
                    recERI1(m+1,nax-2,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) )&
                    +(Ppx-Ax)*recERI1(m,nax-1,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                    +eq*(QQx-PPx)*epq1*recERI1(m+1,nax-1,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                    +0.5d0*ncx*epq1* recERI1(m+1,nax-1,nay,naz,nbx,nby,nbz,ncx-1,ncy,ncz,ndx,ndy,ndz,&
                                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)                                                              
                 return
                else    !!!! (ab|00)x  nbx/=0
                tem=0.5d0*(ncx-1)/eq*(&
                recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx-2,ncy,ncz,ndx,ndy,ndz,&
                                ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                -ep*epq1*&
                recERI1(m+1,nax,nay,naz,nbx,nby,nbz,ncx-2,ncy,ncz,ndx,ndy,ndz,&
                                ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) )&
                +(QQx-Cx)*recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx-1,ncy,ncz,ndx,ndy,ndz,&
                                ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                +ep*(PPx-QQx)*epq1*recERI1(m+1,nax,nay,naz,nbx,nby,nbz,ncx-1,ncy,ncz,ndx,ndy,ndz,&
                                ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                +0.5d0*nax*epq1* recERI1(m+1,nax-1,nay,naz,nbx,nby,nbz,ncx-1,ncy,ncz,ndx,ndy,ndz,&
                                ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)   
                 return
                end if
            else !!! nbx/=0  !!!! (ab|c0)x
                    tem=recERI1(m,nax+1,nay,naz,nbx-1,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                    +(Ax-Bx)*recERI1(m,nax,nay,naz,nbx-1,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)                                                                                        
                return
            end if
        else    !!!! ndx /= 0       (ab|cd)x
            tem=recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx+1,ncy,ncz,ndx-1,ndy,ndz,&
                                ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                +(Cx-Dx)*recERI1(m,nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx-1,ndy,ndz,&
                                ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                                QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) 
            return
        end if           
    end function recERI1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                      
    recursive function recERI2(m,naxx,nayy,nazz,nbxx,nbyy,nbzz,ncxx,ncyy,nczz,ndxx,ndyy,ndzz,&
                                epp,eqq,ABK,CDK,Axx,Bxx,Cxx,Dxx,Px,Qx,Ayy,Byy,Cyy,Dyy,Py,&
                                Qy,Azz,Bzz,Czz,Dzz,Pz,Qz,TT)  result(tem)
!!!!!                                
        use Global,only:PI,F,F_pade,F_ser,dexchange,iexchange
        implicit double precision (A-H,O-Z)  
        logical amax  
!!!! copy values for exchange
        nax=naxx
        nay=nayy
        naz=nazz
        nbx=nbxx
        nby=nbyy
        nbz=nbzz
        ncx=ncxx
        ncy=ncyy
        ncz=nczz
        ndx=ndxx
        ndy=ndyy
        ndz=ndzz
!!!!
        ep=epp
        eq=eqq
!!!!
        ax=Axx
        bx=Bxx
        cx=Cxx
        dx=Dxx
        ppx=Px
        QQx=qx
        ay=Ayy
        by=Byy
        cy=Cyy
        dy=Dyy
        ppy=Py
        QQy=qy
        az=Azz
        bz=Bzz
        cz=Czz
        dz=Dzz
        ppz=Pz
        QQz=qz
!!!!!!!!!!!!!!!!!    
        if(nax<0.or.nbx<0.or.nay<0.or.nby<0.or.naz<0.or.nbz<0.or.&
         ncx<0.or.ndx<0.or.ncy<0.or.ndy<0.or.ncz<0.or.ndz<0)then
            tem=0d0
            return
        end if
!!!!      
        tem=0d0       
        amax=.true.
        epq1=1d0/(ep+eq)
        Wx=(ep*PPx+eq*QQx)*epq1
        Wy=(ep*Ppy+eq*Qqy)*epq1
        Wz=(ep*Ppz+eq*Qqz)*epq1
!!!!!!!!!!!!! Recursive method 2   
!!!!!!!       S.Obara;A.saika  The Journal of Chemical Physics 84, 3963 (1986); doi: 10.1063/1.450106
    if (nax==0.and.nbx==0.and.ncx==0.and.ndx==0)then
        if (nay==0.and.nby==0.and.ncy==0.and.ndy==0)then
            if (naz==0.and.nbz==0.and.ncz==0.and.ndz==0)then
                tem=2d0*ABK*CDK*pi**2.5d0/(ep*eq)/sqrt(ep+eq)*F(m,TT)
                ! tem=2*ABK*CDK/(ep*eq)*pi**2.5/sqrt(ep+eq)*F_pade(m,TT)
                ! tem=2d0*ABK*CDK/(ep*eq)*pi**2.5/sqrt(ep+eq)*F_ser(m,TT)
                return
            else    !!!!z
!!!! make sure naz is max  
                if(naz>=nbz) then 
                    amax=naz<ncz .or. naz<ndz
                else 
                    amax=nbz<ncz .or. nbz<ndz
                end if    
                if(amax)then !!! exchange ac,bd  a-><-c  b-><-d
!!!!!!! Direct exchange
                    ! chtem=ep
                    ! ep=eq
                    ! eq=chtem
!!!!!!! Use function Iexchange/Dexchange                    
                    call dexchange(ep,eq)                              
!!!! exchange z
!!!!!!!! Direct exchange
                    ! nchtem=naz
                    ! naz=ncz
                    ! ncz=nchtem
                    ! chtem=Az
                    ! Az=Cz
                    ! Cz=chtem    
                    ! nchtem=nbz
                    ! nbz=ndz
                    ! ndz=nchtem
                    ! chtem=Bz
                    ! Bz=Dz
                    ! Dz=chtem 
                    ! chtem=Ppz
                    ! Ppz=QQz
                    ! Qqz=chtem  
!!!!!!!! Use function Iexchange/Dexchange
                    call Iexchange(naz,ncz)                              
                    call dexchange(Az,Cz)                              
                    call Iexchange(nbz,ndz)                              
                    call dexchange(Bz,Dz)  
                    call dexchange(PPz,QQz)                                
                end if 
!!!! make sure nax>=nbx 
                if(naz<nbz)then !!!! exchange a,b
!!!!!!!!! Direct exchange                
                    ! nchtem=naz
                    ! naz=nbz
                    ! nbz=nchtem
                    ! chtem=Az
                    ! Az=Bz
                    ! Bz=chtem
!!!!!!!! Use function Iexchange/Dexchange
                    call Iexchange(naz,nBz)                              
                    call dexchange(Az,Bz)   
                end if                
                tem=(Ppz-Az)*recERI2(m,nax,nay,naz-1,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                    +(Wz-Ppz)*recERI2(m+1,nax,nay,naz-1,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                    +0.5d0/ep*(naz-1)*(recERI2(m,nax,nay,naz-2,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                            -eq*epq1*recERI2(m+1,nax,nay,naz-2,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)   )&
                    + 0.5d0/ep*nbz*(recERI2(m,nax,nay,naz-1,nbx,nby,nbz-1,ncx,ncy,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                            -eq*epq1*recERI2(m+1,nax,nay,naz-1,nbx,nby,nbz-1,ncx,ncy,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)   )& 
                    +0.5d0*ncz*epq1 *recERI2(m+1,nax,nay,naz-1,nbx,nby,nbz,ncx,ncy,ncz-1,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) & 
                    +0.5d0*ndz*epq1 *recERI2(m+1,nax,nay,naz-1,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz-1,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)      
                return                         
            end if
        else    !!!!   y
!!!! make sure nay is max  
            if(nay>=nby) then 
                amax=nay<ncy .or. nay<ndy
            else 
                amax=nby<ncy .or. nby<ndy
            end if 
            if(amax)then !!! exchange ac,bd  a-><-c  b-><-d 
!!!!!!!! Direct exchange            
                ! chtem=ep
                ! ep=eq
                ! eq=chtem   
!!!!!!!! Use function Iexchange/Dexchange
                call dexchange(ep,eq)                                             
!!!! change y
!!!!!!!! Direct exchange
                ! nchtem=nay
                ! nay=ncy
                ! ncy=nchtem
                ! chtem=Ay
                ! Ay=Cy
                ! Cy=chtem      
                ! nchtem=nby
                ! nby=ndy
                ! ndy=nchtem
                ! chtem=By
                ! By=Dy
                ! Dy=chtem  
                ! chtem=PPy
                ! PPy=Qqy
                ! Qqy=chtem  
!!!!!!!! Use function Iexchange/Dexchange
                call Iexchange(nay,ncy)                              
                call dexchange(Ay,Cy)                              
                call Iexchange(nby,ndy)                              
                call dexchange(By,Dy)  
                call dexchange(PPy,QQy)                                
!!!! change z
                ! nchtem=naz
                ! naz=ncz
                ! ncz=nchtem
                ! chtem=Az
                ! Az=Cz
                ! Cz=chtem      
                ! nchtem=nbz
                ! nbz=ndz
                ! ndz=nchtem
                ! chtem=Bz
                ! Bz=Dz
                ! Dz=chtem 
                ! chtem=Ppz
                ! Ppz=QQz
                ! QQz=chtem   
!!!!!!!! Use function Iexchange/Dexchange
                call Iexchange(naz,ncz)                              
                call dexchange(Az,Cz)                              
                call Iexchange(nbz,ndz)                              
                call dexchange(Bz,Dz)  
                call dexchange(PPz,QQz)                      
            end if
!!!! make sure nay>=nby
                if(nay<nby)then
!!!!! exchange ab y
!!!!!!!! Direct exchange
                    ! nchtem=nay
                    ! nay=nby
                    ! nby=nchtem
                    ! chtem=Ay
                    ! Ay=By
                    ! By=chtem
!!!!!!!! Use function Iexchange/Dexchange
                    call Iexchange(nay,nBy)                              
                    call dexchange(Ay,By)   
!!!!! exchange ab z
!!!!!!!! Direct exchange
                    ! nchtem=naz
                    ! naz=nbz
                    ! nbz=nchtem
                    ! chtem=Az
                    ! Az=Bz
                    ! Bz=chtem
!!!!!!!! Use function Iexchange/Dexchange
                    call Iexchange(naz,nBz)                              
                    call dexchange(Az,Bz)   
                end if                 
                tem=(Ppy-Ay)*recERI2(m,nax,nay-1,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                    +(Wy-Ppy)*recERI2(m+1,nax,nay-1,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                    ! +eq*epq1*(QQy-Ppy)*recERI2(m+1,nax,nay,naz-1,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                    +0.5d0/ep*(nay-1)*(recERI2(m,nax,nay-2,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                            -eq*epq1*recERI2(m+1,nax,nay-2,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)   )&
                    + 0.5d0/ep*nby*(recERI2(m,nax,nay-1,naz,nbx,nby-1,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                            -eq*epq1*recERI2(m+1,nax,nay-1,naz,nbx,nby-1,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)   )& 
                    +0.5d0*ncy*epq1 *recERI2(m+1,nax,nay-1,naz,nbx,nby,nbz,ncx,ncy-1,ncz,ndx,ndy,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) & 
                    +0.5d0*ndy*epq1 *recERI2(m+1,nax,nay-1,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy-1,ndz,&
                            ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                            QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)     
                return                                
            end if 
    else    !!!!x
!!!! make sure nax is max  
        if(nax>=nbx) then 
            amax=nax<ncx .or. nax<ndx
        else 
            amax=nbx<ncx .or. nbx<ndx
        end if 
        if(amax)then !!! exchange ab,cd  a-><-c  b-><-d
!!!!!!!! Direct exchange        
            ! chtem=ep
            ! ep=eq
            ! eq=chtem  
!!!!!!!! Use function Iexchange/Dexchange                     
            call dexchange(ep,eq)                                 
!!!! change x  
!!!!!!!! Direct exchange                          
            ! nchtem=nax
            ! nax=ncx
            ! ncx=nchtem
            ! chtem=Ax
            ! Ax=Cx
            ! Cx=chtem      
            ! nchtem=nbx
            ! nbx=ndx
            ! ndx=nchtem
            ! chtem=Bx
            ! Bx=Dx
            ! Dx=chtem   
            ! chtem=Ppx
            ! Ppx=QQx
            ! QQx=chtem      
!!!!!!!! Use function Iexchange/Dexchange
            call Iexchange(nax,ncx)                              
            call dexchange(Ax,Cx)                              
            call Iexchange(nbx,ndx)                              
            call dexchange(Bx,Dx)  
            call dexchange(ppx,qqx)                       
!!!! change y
!!!!!!!! Direct exchange
            ! nchtem=nay
            ! nay=ncy
            ! ncy=nchtem
            ! chtem=Ay
            ! Ay=Cy
            ! Cy=chtem      
            ! nchtem=nby
            ! nby=ndy
            ! ndy=nchtem
            ! chtem=By
            ! By=Dy
            ! Dy=chtem   
            ! chtem=Ppy
            ! Ppy=QQy
            ! QQy=chtem  
!!!!!!!! Use function Iexchange/Dexchange
            call Iexchange(nay,ncy)                              
            call dexchange(Ay,Cy)                              
            call Iexchange(nby,ndy)                              
            call dexchange(By,Dy)              
!!!! change z
!!!!!!!! Direct exchange
            ! nchtem=naz
            ! naz=ncz
            ! ncz=nchtem
            ! chtem=Az
            ! Az=Cz
            ! Cz=chtem      
            ! nchtem=nbz
            ! nbz=ndz
            ! ndz=nchtem
            ! chtem=Bz
            ! Bz=Dz
            ! Dz=chtem 
            ! chtem=Ppz
            ! Ppz=QQz
            ! QQz=chtem      
!!!!!!!! Use function Iexchange/Dexchange
            call Iexchange(naz,ncz)                              
            call dexchange(Az,Cz)                              
            call Iexchange(nbz,ndz)                              
            call dexchange(Bz,Dz)  
            call dexchange(PPz,QQz)                         
        end if   
!!!! make sure nax>=nbx 
        if(nax<nbx)then
!!!! exchange ab
!!!!!!!! Direct exchange
            ! nchtem=nax
            ! nax=nbx
            ! nbx=nchtem
            ! chtem=Ax
            ! Ax=Bx
            ! Bx=chtem
!!!!!!!! Use function Iexchange/Dexchange
            call Iexchange(nax,nBx)                              
            call dexchange(Ax,Bx)   
!!!! change y
!!!!!!!! Direct exchange
            ! nchtem=nay
            ! nay=nby
            ! nby=nchtem
            ! chtem=Ay
            ! Ay=By
            ! By=chtem
!!!!!!!! Use function Iexchange/Dexchange
            call Iexchange(nay,nBy)                              
            call dexchange(Ay,By)   
!!!! change z
!!!!!!!! Direct exchange
            ! nchtem=naz
            ! naz=nbz
            ! nbz=nchtem
            ! chtem=Az
            ! Az=Bz
            ! Bz=chtem
!!!!!!!! Use function Iexchange/Dexchange
            call Iexchange(naz,nBz)                              
            call dexchange(Az,Bz)   
        end if               
        tem=(Ppx-Ax)*recERI2(m,nax-1,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
            +(Wx-Ppx)*recERI2(m+1,nax-1,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
            +0.5d0/ep*(nax-1)*(recERI2(m,nax-2,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                    -eq*epq1*recERI2(m+1,nax-2,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)   )&
            + 0.5d0/ep*nbx*(recERI2(m,nax-1,nay,naz,nbx-1,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) &
                    -eq*epq1*recERI2(m+1,nax-1,nay,naz,nbx-1,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)   )& 
            +0.5d0*ncx*epq1 *recERI2(m+1,nax-1,nay,naz,nbx,nby,nbz,ncx-1,ncy,ncz,ndx,ndy,ndz,&
                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT) & 
            +0.5d0*ndx*epq1 *recERI2(m+1,nax-1,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx-1,ndy,ndz,&
                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,PPx,QQx,Ay,By,Cy,Dy,PPy,&
                    QQy,Az,Bz,Cz,Dz,PPz,QQz,TT)  
            return
    end if       
    end function recERI2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function anaERI1(nax,nay,naz,nbx,nby,nbz,ncx,ncy,ncz,ndx,ndy,ndz,&
                    ep,eq,ABK,CDK,Ax,Bx,Cx,Dx,Px,Qx,Ay,By,Cy,Dy,Py,&
                    Qy,Az,Bz,Cz,Dz,Pz,Qz,TT)                             
!!!!!                                
        use Global
        implicit double precision (A-H,O-Z)  
!!!!!!! 3.(Analytical method) O.Ohata method (1966): https://doi.org/10.1143/JPSJ.21.2313
        anaERI1=0d0
        delt=1d0/(4*ep)+1d0/(4*eq)
        erixyz=0d0
        erix=0d0
        eriy=0d0
        eriz=0d0
!!!  x                          
        do ii1=0,nax+nbx 
        do ii2=0,ncx+ndx
            do ir1=0,floor(ii1*0.5d0)
            do ir2=0,floor(ii2*0.5d0)                
                do iu=0,floor((ii1+ii2-2*ir1-2*ir2)*0.5d0)
!!!
            ! erix=(-1d0)**(ii2+iu)*fi(ii1,nax,nbx,Px-Ax,Px-Bx)&
            !         *fi(ii2,ncx,ndx,Qx-Cx,Qx-Dx)*ifact(ii1)*ifact(ii2)&
            !         *(4d0*ep)**(ir1-ii1)*(4d0*eq)**(ir2-ii2)&
            !         *delt**(2*(ir1+ir2)+iu-ii1-ii2)&
            !         /(ifact(ir1)*ifact(ir2)*ifact(ii1-2*ir1)*ifact(ii2-2*ir2))&
            !         *ifact(ii1+ii2-2*(ir1+ir2))*(Qx-Px)**(ii1+ii2-2*(iu+ir1+ir2))&
            !         /(ifact(iu)*ifact(ii1+ii2-2*(ir1+ir2+iu)))
!!!! use function BFanaERI1
            erix=BFanaERI1(nax,nbx,ncx,ndx,Ax,Bx,Cx,Dx,Px,Qx,ep,eq,delt,&
                    ii1,ii2,ir1,ir2,iu)
!!!  y   
            do jj1=0,nay+nby 
            do jj2=0,ncy+ndy 
                do is1=0,floor(jj1*0.5d0)
                do is2=0,floor(jj2*0.5d0)
                    do iv=0,floor((jj1+jj2-2*is1-2*is2)*0.5d0)
!!!!                    
                        ! eriy=(-1d0)**(jj2+iv)*fi(jj1,nay,nby,Py-Ay,Py-By)&
                        !         *fi(jj2,ncy,ndy,Qy-Cy,Qy-Dy)*ifact(jj1)*ifact(jj2)&
                        !         *(4d0*ep)**(is1-jj1)*(4d0*eq)**(is2-jj2)&
                        !         *delt**(2*(is1+is2)+iv-jj1-jj2)&
                        !         /(ifact(is1)*ifact(is2)*ifact(jj1-2*is1)*ifact(jj2-2*is2))&
                        !         *ifact(jj1+jj2-2*(is1+is2))*(Qy-Py)**(jj1+jj2-2*(iv+is1+is2))&
                        !         /(ifact(iv)*ifact(jj1+jj2-2*(is1+is2+iv))) 
!!!! use function BFanaERI1
                    eriy=BFanaERI1(nay,nby,ncy,ndy,Ay,By,Cy,Dy,Py,qy,ep,eq,delt,&
                    jj1,jj2,is1,is2,iv)    
!!!  z              
                        do kk1=0,naz+nbz 
                        do kk2=0,ncz+ndz
                            do it1=0,floor(kk1*0.5d0)
                            do it2=0,floor(kk2*0.5d0)
                                do iw=0,floor((kk1+kk2-2*it1-2*it2)*0.5d0)
!!!!                                
                                ! eriz=(-1d0)**(kk2+iw)*fi(kk1,naz,nbz,Pz-Az,Pz-Bz)&
                                !         *fi(kk2,ncz,ndz,Qz-Cz,Qz-Dz)*ifact(kk1)*ifact(kk2)&
                                !         *(4d0*ep)**(it1-kk1)*(4d0*eq)**(it2-kk2)&
                                !         *delt**(2*(it1+it2)+iw-kk1-kk2)&
                                !         /(ifact(it1)*ifact(it2)*ifact(kk1-2*it1)*ifact(kk2-2*it2))&
                                !         *ifact(kk1+kk2-2*(it1+it2))*(Qz-Pz)**(kk1+kk2-2*(iw+it1+it2))&
                                !         /(ifact(iw)*ifact(kk1+kk2-2*(it1+it2+iw)))        
!!!! use function BFanaERI1
                    eriz=BFanaERI1(naz,nbz,ncz,ndz,Az,Bz,Cz,Dz,Pz,qz,ep,eq,delt,&
                    kk1,kk2,it1,it2,iw)                               
                        !!!!
                                    mu=ii1+jj1+kk1+ii2+jj2+kk2-2*(ir1+is1+it1+ir2+is2+it2)-(iu+iv+iw)
                                    erixyz=erixyz+erix*eriy*eriz*F(mu,TT)
                                    ! erixyz=erixyz+erix*eriy*eriz*F_pade(mu,TT)                                    
                                    ! erixyz=erix*eriy*eriz*f(mu,TT)                                        
                                        end do  !iw
                                    end do  !it2                                    
                                    end do  !it1
                                end do  !kk2                                 
                                end do  !kk1             
                            end do  !iv
                        end do  !is2                         
                        end do  !is1 
                    end do  !jj2                        
                    end do  !jj1                                      
                end do !iu
            end do  !ir2            
            end do  !ir1
        end do  !ii2
        end do  !ii1
        anaERI1=erixyz*2*pi**2*ABK*CDK/(ep*eq)*sqrt(pi/(eq+ep))
        return
    end function anaERI1
!!!!!!!
    function BFanaERI1(nax,nbx,ncx,ndx,Ax,Bx,Cx,Dx,Px,Qx,ep,eq,delt,&
                        i1,i2,ir1,ir2,iu)
            use Global ,only:ifact,fi
            implicit double precision (A-H,O-Z)  
            tem=(-1)**(i1)*fi(i1,nax,nbx,Px-Ax,Px-Bx)*fi(i2,ncx,ndx,Qx-Cx,Qx-Dx)&
                *ifact(i1)*ifact(i2)/(4*delt)**(i1+i2)&
                *(ep)**(ir1-i1)*(eq)**(ir2-i2)*(2*delt)**(2*ir1+2*ir2)&
                /ifact(ir1)/ifact(ir2)/ifact(i1-2*ir1)/ifact(i2-2*ir2)&
                *ifact(i1+i2-2*(ir1+ir2))*(-1)**iu*(px-qx)**(i1+i2-2*(ir1+ir2)-2*iu)&
                *delt**iu/ifact(iu)/ifact(i1+i2-2*(ir1+ir2)-2*iu)
            BFanaERI1=tem
            return
    end function BFanaERI1
!======================================================================!
 function recERIfast(naxx,nayy,nazz,nbxx,nbyy,nbzz,ncxx,ncyy,nczz,ndxx,ndyy,ndzz,&
                epp,eqq,ABK,CDK,Axx,Bxx,Cxx,Dxx,Px,Qx,Ayy,Byy,Cyy,Dyy,Py,&
                 Qy,Azz,Bzz,Czz,Dzz,Pz,Qz,TT)
!!!!!                                
        use Global,only:PI,F,iexchange,dexchange,F_pade
        implicit double precision (A-H,O-Z)    
        allocatable::temp1(:)
!!!! copy values for exchange
        nax=naxx
        nay=nayy
        naz=nazz
        nbx=nbxx
        nby=nbyy
        nbz=nbzz
        ncx=ncxx
        ncy=ncyy
        ncz=nczz
        ndx=ndxx
        ndy=ndyy
        ndz=ndzz
!!!!
        ep=epp
        eq=eqq
!!!!
        ax=Axx
        bx=Bxx
        cx=Cxx
        dx=Dxx
        ppx=Px
        QQx=qx
        ay=Ayy
        by=Byy
        cy=Cyy
        dy=Dyy
        ppy=Py
        QQy=qy
        az=Azz
        bz=Bzz
        cz=Czz
        dz=Dzz
        ppz=Pz
        QQz=qz
!!!!!!!!!!!!!!!!!    
!!!!        
        if(nax<0.or.nbx<0.or.nay<0.or.nby<0.or.naz<0.or.nbz<0.or.&
         ncx<0.or.ndx<0.or.ncy<0.or.ndy<0.or.ncz<0.or.ndz<0)then
           stop 'bad input for recERIfast'
        end if
        Lsum=nax+nay+naz+nbx+nby+nbz+ncx+ncy+ncz+ndx+ndy+ndz
        allocate(temp1(0:Lsum))  
        do i=0,Lsum
            temp1(i)=F(i,TT)
            ! temp1(i)=F_pade(i,TT)
        end do        
!!!!           
        if(nax<nbx)  then
            call iexchange(nax,nbx)
            call dexchange(ax,bx)
        end if 
        if(nay<nby)  then
            call  iexchange(nay,nby)
            call  dexchange(ay,by)
        end if   
        if(naz<nbz)  then
            call   iexchange(naz,nbz)
            call   dexchange(az,bz)
        end if    
!!!
        if(ncx<ndx)  then
            call iexchange(ncx,ndx)
            call dexchange(cx,dx)
        end if 
        if(ncy<ndy)  then
            call  iexchange(ncy,ndy)
            call  dexchange(cy,dy)
        end if   
        if(ncz<ndz)  then
            call   iexchange(ncz,ndz)
            call   dexchange(cz,dz)
        end if       
!!!
        if(nax+nbx<ncx+ndx)then
            call dexchange(ep,eq)            
            call iexchange(nax,ncx)
            call dexchange(ax,cx)   
            call iexchange(nbx,ndx)
            call dexchange(bx,dx)  
            call dexchange(ppx,qqx)   
            call iexchange(nay,ncy)
            call dexchange(ay,cy)   
            call iexchange(nby,ndy)
            call dexchange(by,dy)     
            call dexchange(ppy,qqy)   
            call iexchange(naz,ncz)
            call dexchange(az,cz)   
            call iexchange(nbz,ndz)
            call dexchange(bz,dz)     
            call dexchange(ppz,qqz)                         
        end if  
        call recERI1xyz(nax,nbx,ncx,ndx,Lsum,ep,eq,Ax,Bx,Cx,Dx,Ppx,Qqx,temp1,1)
        if(nay+nby<ncy+ndy)then
            call dexchange(ep,eq)            
            call iexchange(nay,ncy)
            call dexchange(ay,cy)   
            call iexchange(nby,ndy)
            call dexchange(by,dy)     
            call dexchange(ppy,qqy)   
            call iexchange(naz,ncz)
            call dexchange(az,cz)   
            call iexchange(nbz,ndz)
            call dexchange(bz,dz)     
            call dexchange(ppz,qqz)                
        end if              
        call recERI1xyz(nay,nby,ncy,ndy,Lsum,ep,eq,Ay,By,Cy,Dy,Ppy,Qqy,temp1,1)
        if(naz+nbz<ncz+ndz)then
            call dexchange(ep,eq)            
            call iexchange(naz,ncz)
            call dexchange(az,cz)   
            call iexchange(nbz,ndz)
            call dexchange(bz,dz)     
            call dexchange(ppz,qqz)  
        end if                           
        call recERI1xyz(naz,nbz,ncz,ndz,Lsum,ep,eq,Az,Bz,Cz,Dz,Ppz,Qqz,temp1,0)
!!!!         
        recERIfast=temp1(0)*ABK*CDK*2d0*pi**2.5d0/(eq*ep)/sqrt(eq+ep)
        deallocate(temp1)
    end function recERIfast
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine recERI1xyz(nax,nbx,ncx,ndx,Lsum,ep,eq,Ax,Bx,Cx,Dx,Px,Qx,Finreout,Notz)
        implicit double precision (A-H,O-Z)    
        dimension::Finreout(0:Lsum),temp1(0:nax+nbx+1,0:ncx+ndx+1,0:Lsum),&
                temp2(0:nax+nbx+1,0:nbx),temp3(0:nax+nbx+1,0:nbx,0:ncx+ndx+1,0:ndx),&
                 temp4(0:ncx+ndx+1,0:ndx)
!!!! init
        if(Lsum==0) return
        nx1=nax+nbx
        nx2=ncx+ndx
        nx=nx1+nx2
        if(nx==0) return        
        epq1=1d0/(ep+eq)
        Wx=(ep*Px+eq*Qx)*epq1
!!!!
        temp1(0,0,:)=Finreout(:)
        temp1(1,0,Lsum)=(Px-Ax)*Finreout(Lsum)
        temp1(0,1,Lsum)=(Qx-Cx)*Finreout(Lsum)
        temp1(1,1,Lsum)=(Qx-Cx)*temp1(1,0,Lsum)
!!!! Lsum=1
        if(Lsum==1) then
            if(nx1==1)then
                temp1(1,0,0)=(Px-Ax)*Finreout(0)+(Wx-Px)*Finreout(1)
                Finreout(0:1)=temp1(1,0,0:1)
                return
            else
                temp1(0,1,0)=(Qx-Cx)*Finreout(0)+(Wx-Qx)*Finreout(1)
                Finreout(0:1)=temp1(0,1,0:1)
                return                
            end if
        end if
!!!! nx1=1
        if(nx1==1)then
            if(nx2==0)then
                do m=Lsum-1,0,-1
                    temp1(0,0,m)=Finreout(m)
                    temp1(1,0,m)=(Px-Ax)*Finreout(m)+(Wx-Px)*Finreout(m+1)
                end do
                if(Notz==0)then
                    Finreout(0)=temp1(1,0,0)
                    return
                end if
                do m=0,Lsum
                    Finreout(m)=temp1(1,0,m)
                end do
                return      
            else        !!!! nx2=1, because nx1>=nx2 
                do m=Lsum-1,0,-1
                    temp1(0,0,m)=Finreout(m)
                    temp1(1,0,m)=(Px-Ax)*temp1(0,0,m)+(Wx-Px)*temp1(0,0,m+1)
                    temp1(0,1,m)=(Qx-Cx)*temp1(0,0,m)+(Wx-Qx)*temp1(0,0,m+1)
                    temp1(1,1,m)=0.5d0*epq1*temp1(0,0,m+1)+(Qx-Cx)*temp1(1,0,m)&
                        +(Wx-Qx)*temp1(1,0,m+1)
                end do    
                if(Notz==0)then
                    Finreout(0)=temp1(1,1,0)
                    return
                end if
                do m=0,Lsum
                    Finreout(m)=temp1(1,1,m)
                end do
                return   
            end if
        end if
!!!! nx1 >= 2
!!! nx2=0
        if(nx2==0)then
            do i=1,nx1
                temp1(i+1,0,Lsum)=i*0.5d0/ep*temp1(i-1,0,Lsum)&
                        +(Px-Ax)*temp1(i,0,Lsum)
            end do
            do m=Lsum-1,0,-1
                temp1(0,0,m)=Finreout(m)
                temp1(1,0,m)=(Px-Ax)*Finreout(m)+(Wx-Px)*Finreout(m+1)
                do i=1,nx1
                    temp1(i+1,0,m)=0.5d0*i/ep*(temp1(i-1,0,m)-eq*epq1*temp1(i-1,0,m+1))&
                    +(Px-Ax)*temp1(i,0,m)+(Wx-Px)*temp1(i,0,m+1)
                end do
            end do
            if(nbx .ne. 0)then 
                if (Notz==1)then
                    do m=0,Lsum
                        do i=0,nx1+1
                            temp2(i,0)=temp1(i,nx2,m)
                        end do            
                        do  j=1,nbx
                            do  i=nx1,nax,-1
                                temp2(i,j)=temp2(i+1,j-1)+(Ax-Bx)*temp2(i,j-1)
                            end do
                        end do 
                        Finreout(m)=temp2(nax,nbx)
                    end do
                    return
                else
                    do i=0,nx1+1
                        temp2(i,0)=temp1(i,nx2,0)
                    end do  
                    do  j=1,nbx
                        do  i=nx1,nax,-1
                            temp2(i,j)=temp2(i+1,j-1)+(Ax-Bx)*temp2(i,j-1)
                        end do
                    end do           
                    Finreout(0)=temp2(nax,nbx) 
                    return
                end if 
            else 
                if(Notz==0)then
                    Finreout(0)=temp1(nx1,0,0)
                    return
                end if
                do m=0,Lsum
                    Finreout(m)=temp1(nx1,0,m)
                end do
                return                       
            end if                
!!! nx2 = 1            
        else if(nx2==1)then
            do i=1,nx1
                temp1(i+1,0,Lsum)=i*0.5d0/ep*temp1(i-1,0,Lsum)&
                        +(Px-Ax)*temp1(i,0,Lsum)
                temp1(i+1,1,Lsum)=(qx-cx)*temp1(i+1,0,Lsum)                    
            end do   
            do m=Lsum-1,0,-1
                temp1(0,0,m)=Finreout(m)
                temp1(1,0,m)=(Px-Ax)*Finreout(m)+(Wx-Px)*Finreout(m+1)
                temp1(0,1,m)=(Qx-Cx)*Finreout(m)+(Wx-Qx)*Finreout(m+1)
                temp1(1,1,m)=0.5d0*epq1*temp1(0,0,m+1)+(Qx-Cx)*temp1(1,0,m)&
                        +(Wx-Qx)*temp1(1,0,m+1)
                do i=1,nx1
                    temp1(i+1,0,m)=0.5d0*i/ep*(temp1(i-1,0,m)-eq*epq1*temp1(i-1,0,m+1))&
                    +(Px-Ax)*temp1(i,0,m)+(Wx-Px)*temp1(i,0,m+1)
                    temp1(i+1,1,m)=(i+1)*0.5d0*epq1*temp1(i,0,m+1)+(Qx-Cx)*temp1(i+1,0,m)&
                        +(Wx-Qx)*temp1(i+1,0,m+1)
                end do   
            end do      
            if(nbx .ne. 0)then 
                if (Notz==1)then
                    do m=0,Lsum
                        do i=0,nx1+1
                            temp2(i,0)=temp1(i,nx2,m)
                        end do            
                        do  j=1,nbx
                            do  i=nx1,nax,-1
                                temp2(i,j)=temp2(i+1,j-1)+(Ax-Bx)*temp2(i,j-1)
                            end do
                        end do 
                        Finreout(m)=temp2(nax,nbx)
                    end do
                    return
                else
                    do i=0,nx1+1
                        temp2(i,0)=temp1(i,nx2,0)
                    end do  
                    do  j=1,nbx
                        do  i=nx1,nax,-1
                            temp2(i,j)=temp2(i+1,j-1)+(Ax-Bx)*temp2(i,j-1)
                        end do
                    end do           
                    Finreout(0)=temp2(nax,nbx) 
                    return
                end if 
            else 
                if(Notz==0)then
                    Finreout(0)=temp1(nx1,1,0)
                    return
                end if
                do m=0,Lsum
                    Finreout(m)=temp1(nx1,1,m)
                end do
                return                       
            end if                                                
    else !!!!  nx2>=2
!!!! general
        do i=1,nx1
            temp1(i+1,0,Lsum)=i*0.5d0/ep*temp1(i-1,0,Lsum)&
                    +(Px-Ax)*temp1(i,0,Lsum)
            temp1(i+1,1,Lsum)=(qx-cx)*temp1(i+1,0,Lsum)             
            do j=1,nx2 
                temp1(0,j+1,Lsum)=j*0.5d0/eq*temp1(0,j-1,Lsum)&
                        +(Qx-Cx)*temp1(0,j,Lsum) 
                temp1(1,j+1,Lsum)=j*0.5d0/eq*temp1(1,j-1,Lsum)&
                        +(Qx-Cx)*temp1(1,j,Lsum)
                temp1(i+1,j+1,Lsum)=j*0.5d0/eq*temp1(i+1,j-1,Lsum)&
                    +(Qx-Cx)*temp1(i+1,j,Lsum)                       
            end do      
        end do                 
        do m=Lsum-1,0,-1
            temp1(1,0,m)=(Px-Ax)*Finreout(m)+(Wx-Px)*Finreout(m+1)
            temp1(0,1,m)=(Qx-Cx)*Finreout(m)+(Wx-Qx)*Finreout(m+1)
            temp1(1,1,m)=0.5d0*epq1*temp1(0,0,m+1)+(Qx-Cx)*temp1(1,0,m)&
                    +(Wx-Qx)*temp1(1,0,m+1)  
            do i=1,nx1
                temp1(i+1,0,m)=0.5d0*i/ep*(temp1(i-1,0,m)-eq*epq1*temp1(i-1,0,m+1))&
                    +(Px-Ax)*temp1(i,0,m)+(Wx-Px)*temp1(i,0,m+1)
                temp1(i+1,1,m)=(i+1)*0.5d0*epq1*temp1(i,0,m+1)+(Qx-Cx)*temp1(i+1,0,m)&
                        +(Wx-Qx)*temp1(i+1,0,m+1)                                    
                do j=1,nx2
                    temp1(0,j+1,m)=0.5d0*j/eq*(temp1(0,j-1,m)-ep*epq1*temp1(0,j-1,m+1))&
                       +(Qx-Cx)*temp1(0,j,m)+(Wx-Qx)*temp1(0,j,m+1)
                    temp1(1,j+1,m)= 0.5d0*j/eq*(temp1(1,j-1,m)-ep*epq1*temp1(1,j-1,m+1))&
                        +0.5d0*epq1*temp1(0,j,m+1)+(Qx-Cx)*temp1(1,j,m)&
                        +(Wx-Qx)*temp1(1,j,m+1)
                    temp1(i+1,j+1,m)=0.5d0*j/eq*(temp1(i+1,j-1,m)-ep*epq1*temp1(i+1,j-1,m+1))&
                        +0.5d0*(i+1)*epq1*temp1(i,j,m+1)+(Qx-Cx)*temp1(i+1,j,m)&
                        +(Wx-Qx)*temp1(i+1,j,m+1)                       
                end do              
            end do                 
        end do          
    end if    
!!!
    if(nbx==0 .and. ndx==0)then
        if(Notz==1)then
            Finreout(0:Lsum)=temp1(nx1,nx2,0:Lsum)
            return
        else
            Finreout(0)=temp1(nx1,nx2,0)
            return                
        end if    
    else if(ndx==0) then
        if (Notz==1)then
            do m=0,Lsum
                do i=0,nx1+1
                    temp2(i,0)=temp1(i,nx2,m)
                end do            
                do  j=1,nbx
                    do  i=nx1,nax,-1
                        temp2(i,j)=temp2(i+1,j-1)+(Ax-Bx)*temp2(i,j-1)
                    end do
                end do 
                Finreout(m)=temp2(nax,nbx)
            end do
            return
        else
            do i=0,nx1+1
                temp2(i,0)=temp1(i,nx2,0)
            end do  
            do  j=1,nbx
                do  i=nx1,nax,-1
                    temp2(i,j)=temp2(i+1,j-1)+(Ax-Bx)*temp2(i,j-1)
                end do
            end do           
            Finreout(0)=temp2(nax,nbx) 
            return
        end if   
    else  if(nbx==0) then      
        if (Notz==1)then
            do m=0,Lsum
                do i=0,nx2+1
                    temp4(i,0)=temp1(nx1,i,m)
                end do            
                do  j=1,ndx
                    do  i=nx2,ncx,-1
                        temp4(i,j)=temp4(i+1,j-1)+(Cx-Dx)*temp4(i,j-1)
                    end do
                end do 
                Finreout(m)=temp4(ncx,ndx)
            end do
            return
        else
            do i=0,nx2+1
                temp4(i,0)=temp1(nx1,i,0)
            end do  
            do  j=1,ndx
                do  i=nx2,ncx,-1
                    temp4(i,j)=temp4(i+1,j-1)+(Cx-Dx)*temp4(i,j-1)
                end do
            end do           
            Finreout(0)=temp4(ncx,ndx) 
            return
        end if                          
    else     !!!! nbx>0,ndx>0
        if (Notz==1)then
            do m=0,Lsum
                do i=0,nx1+1
                    do j=0,nx2+1
                        temp3(i,0,j,0)=temp1(i,j,m)
                    end do      
                end do      
                do  j=1,nbx
                    do  i=nx1,nax,-1
                        temp3(i,j,:,0)=temp3(i+1,j-1,:,0)+(Ax-Bx)*temp3(i,j-1,:,0)
                    end do
                end do 
                do  j=1,ndx
                    do  i=nx2,ncx,-1
                        temp3(nax,nbx,i,j)=temp3(nax,nbx,i+1,j-1)+(Cx-Dx)*temp3(nax,nbx,i,j-1)
                    end do
                end do                 
                Finreout(m)=temp3(nax,nbx,ncx,ndx)
            end do
            return
!!!!!!            
        else 
            do i=0,nx1+1
                do j=0,nx2+1
                    temp3(i,0,j,0)=temp1(i,j,0)
                end do      
            end do           
            do  j=1,nbx
                do  i=nx1,nax,-1
                    temp3(i,j,:,0)=temp3(i+1,j-1,:,0)+(Ax-Bx)*temp3(i,j-1,:,0)
                end do
            end do 
            do  j=1,ndx
                do  i=nx2,ncx,-1
                    temp3(nax,nbx,i,j)=temp3(nax,nbx,i+1,j-1)+(Cx-Dx)*temp3(nax,nbx,i,j-1)
                end do
            end do            
            Finreout(0)=temp3(nax,nbx,ncx,ndx)
            return
        end if  
    end if         
!!!!!!!!!!!!!!!!!!!!!!!!        
    end subroutine recERI1xyz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!
subroutine   CINT_ERI
        use Global
        use init
        implicit double precision (A-H,O-Z)
!!!!
        ! integer,parameter :: CHARGE_OF  = 1
        ! integer,parameter :: PTR_COORD  = 2
        ! integer,parameter :: NUC_MOD_OF = 3
        ! integer,parameter :: PTR_ZETA   = 4
        ! integer,parameter :: ATM_SLOTS  = 6
        
        ! integer,parameter :: ATOM_OF    = 1
        ! integer,parameter :: ANG_OF     = 2
        ! integer,parameter :: NPRIM_OF   = 3
        ! integer,parameter :: NCTR_OF    = 4
        ! integer,parameter :: KAPPA_OF   = 5
        ! integer,parameter :: PTR_EXP    = 6
        ! integer,parameter :: PTR_COEFF  = 7
        ! integer,parameter :: BAS_SLOTS  = 8
        
        integer,parameter :: PTR_ENV_START = 21
        
        integer,allocatable :: atm(:,:)
        integer,allocatable :: bas(:,:)
        double precision,allocatable :: env(:)
        integer,external :: CINTcgto_cart,cint2e_ip1_cart
        double precision,external :: CINTgto_norm
        
        integer ::  off
        integer :: di, dj, dk, dl
        integer :: shls(4)
        double precision,allocatable :: buf2e(:,:,:,:),buf1e(:,:)
        ! On 32-bit machine, integer(4) :: opt
        integer(8) :: opt
!!!!!!! only for Pople's basis set (sp)
        allocatable::kcentersp(:),ntypesp(:),ngfsp(:),&
            kbasstsp(:),knbassp(:),kpgfstsp(:),kpgfedsp(:),kbasedsp(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nsp=count(ntype==12)
        if(nsp>0)then
            nshellsp=nshell+nsp
            allocate(kcentersp(nshellsp),ntypesp(nshellsp),ngfsp(nshellsp),&
                    kbasstsp(nshellsp),kbasedsp(nshellsp),knbassp(nshellsp),&
                    kpgfstsp(nshellsp),kpgfedsp(nshellsp))
            ncountsp=0
            do i=1,nshell
                kcentersp(i+ncountsp)=kcenter(i)
                ntypesp(i+ncountsp)=ntype(i)
                ngfsp(i+ncountsp)=ngf(i)
                kbasstsp(i+ncountsp)=kbasst(i)
                kbasedsp(i+ncountsp)=kbased(i)
                knbassp(i+ncountsp)=knbas(i)
                kpgfstsp(i+ncountsp)=kpgfst(i)
                kpgfedsp(i+ncountsp)=kpgfed(i)
                ngfsp(i+ncountsp)=ngf(i)
                if(ntype(i)==12)then !!! SP
!!!! S           
                    ntypesp(i+ncountsp)=1
                    kbasstsp(i+ncountsp)=kbasst(i)
                    kbasedsp(i+ncountsp)=kbasst(i)
                    knbassp(i+ncountsp)=1
                    ngfsp(i+ncountsp)=ngf(i)/2 
                    kpgfstsp(i+ncountsp)=kpgfst(i)
                    kpgfedsp(i+ncountsp)=kpgfstsp(i+ncountsp)+ngfsp(i+ncountsp)-1    
!!!! P                
                    kcentersp(i+ncountsp+1)=kcentersp(i+ncountsp)
                    ntypesp(i+ncountsp+1)=2
                    knbassp(i+ncountsp+1)=3   
                    kbasstsp(i+ncountsp+1)=kbasstsp(i+ncountsp)+1
                    kbasedsp(i+ncountsp+1)=kbasstsp(i+ncountsp)+3
!!!
                    ngfsp(i+ncountsp+1)=ngf(i)/2   
                    kpgfstsp(i+ncountsp+1)=kpgfedsp(i+ncountsp)+1
                    kpgfedsp(i+ncountsp+1)=kpgfstsp(i+ncountsp+1)&
                            +ngfsp(i+ncountsp+1)-1    
                    ncountsp=ncountsp+1          
                end if
            end do   
            write(*,*) nshellsp
            write(*,*) kcentersp
            write(*,*) ntypesp
            write(*,*) knbassp
            write(*,*) kbasstsp
            write(*,*) kbasedsp
            write(*,*) ngfsp 
            write(*,*) kpgfstsp
            write(*,*) kpgfedsp
!!!!!!!!!!            
            ERIall=0d0
            allocate (atm(6,natom))
            allocate (bas(8,nshellsp))
            allocate (env(100000))
            off = PTR_ENV_START
            do i=1,natom       
                atm(1,i) = KATOM(i)
                atm(2,i) = off ! note the 0-based index
                env(off +1 ) = coor(i,1) ! x (Bohr)
                env(off + 2) = coor(i,2) ! y (Bohr)
                env(off + 3) = coor(i,3) ! z (Bohr)
                off = off + 3 
                atm(3,i) = 1
                atm(5,i) = 0
            end do     
            do n=1,nshellsp         
! basis #1
                bas(1 ,n)  = kcentersp(n)-1 ! note that it's the first atom, the index is 0-based
                bas(2 ,n)  = Ntypesp(n)-1
                bas(3 ,n)  = ngfsp(n)
                bas(4 ,n)  = 1
                bas(6  ,n)  = off ! note the 0-based index
                do i=1,ngfsp(n)
                    mpgf=kpgfstsp(n)+i-1
                    env(off + i) = gfe(mpgf)       
                end do 
                off = off + ngfsp(n)
                bas(7,n) = off
                do i=1,ngfsp(n)
                    env(off + i) = gfc(kpgfstsp(n)+i-1)&
                                   *CINTgto_norm(bas(2,n), env(bas(6,n)+i))
                end do
                off = off + ngfsp(n)
            end do        
!!!!!!!!!!!!!!  1e
            do i=1,nshellsp
                shls(1) = i-1
                kdi = CINTcgto_cart(i-1, bas)
                do j=1,nshellsp
                    shls(2) = j-1
                    kdj = CINTcgto_cart(j-1, bas)
                    allocate (buf1e(kdi,kdj))
                    do ibas=kbasstsp(i),kbasedsp(i)                          
                        do jbas=kbasstsp(j),kbasedsp(j)
                            ibu2ei=ibas-kbasstsp(i)+1
                            ibu2ej=jbas-kbasstsp(j)+1
                            call cint1e_ovlp_cart(buf1e, shls, atm, natom, bas, nshellsp, env, opt) 
                            S(ibas,jbas)=buf1e(ibu2ei,ibu2ej)
                            S(jbas,ibas)=S(ibas,jbas)
                            call cint1e_kin_cart(buf1e, shls, atm, natom, bas, nshellsp, env, opt) 
                            T(ibas,jbas)=buf1e(ibu2ei,ibu2ej)
                            T(jbas,ibas)=T(ibas,jbas)
                            call cint1e_nuc_cart(buf1e, shls, atm, natom, bas, nshellsp, env, opt) 
                            V(ibas,jbas)=buf1e(ibu2ei,ibu2ej)
                            V(jbas,ibas)=V(ibas,jbas)
                        end do
                    end do            
                    deallocate (buf1e)
                end do
            end do        
            Hcore=T+V
!!!!!!!!!!!!!!! 2e     
            nshelldone=0   
            call cint2e_cart_optimizer(opt,atm,natom,bas,nshellsp,env)
            do i=1,nshellsp
                nshelldone=nshelldone+1        
                if(mod(nshelldone,5)==0)then
                    write(*,*)real(nshelldone)/nshellsp*100d0,'% ERI done'    
                end if 
                shls(1) = i-1
                kdi = CINTcgto_cart(i-1, bas)
                do j=1,nshellsp
                    shls(2) = j-1
                    kdj = CINTcgto_cart(j-1, bas)
                    do k=1,nshellsp
                        shls(3) = k-1
                        kdk = CINTcgto_cart(k-1, bas)
                        do l=1,nshellsp
                            shls(4) = l-1
                            kdl = CINTcgto_cart(l-1, bas)     
                            allocate (buf2e(kdi,kdj,kdk,kdl))
                            call cint2e_cart(buf2e, shls, atm, natom, bas, nshellsp, env, opt)
!!!                            
            do ibas=kbasstsp(i),kbasedsp(i)                          
                do jbas=kbasstsp(j),kbasedsp(j)
                    do kbas=kbasstsp(k),kbasedsp(k)
                        do lbas=kbasstsp(l),kbasedsp(l)
                            ijkl=indexeri(ibas,jbas,kbas,lbas)
                            ibu2ei=ibas-kbasstsp(i)+1
                            ibu2ej=jbas-kbasstsp(j)+1
                            ibu2ek=kbas-kbasstsp(k)+1
                            ibu2el=lbas-kbasstsp(l)+1
                            ERIall(ijkl)=buf2e(ibu2ei,ibu2ej,ibu2ek,ibu2el)
                        end do
                    end do
                end do
            end do
                        deallocate (buf2e)
                    end do 
                end do 
            end do 
        end do      
            deallocate (atm, bas, env)    
            deallocate(kcentersp,ntypesp,ngfsp,kbasstsp,&
                        kbasedsp,knbassp,kpgfstsp,kpgfedsp)
            return 
        end if !!!! only for Pople basis set.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!  general        
        ERIall=0d0
        allocate (atm(6,natom))
        allocate (bas(8,nshell))
        allocate (env(100000))
        off = PTR_ENV_START
        do i=1,natom       
            atm(1,i) = KATOM(i)
            atm(2,i) = off ! note the 0-based index
            env(off +1 ) = coor(i,1) ! x (Bohr)
            env(off + 2) = coor(i,2) ! y (Bohr)
            env(off + 3) = coor(i,3) ! z (Bohr)
            off = off + 3 
            atm(3,i) = 1
            atm(5,i) = 0
        end do     
        do n=1,nshell         
!!! basis #1
!!! note that it's the first atom, the index is 0-based
            bas(1 ,n)  = kcenter(n)-1 
            bas(2 ,n)  = Ntype(n)-1
            bas(3 ,n)  = ngf(n)
            bas(4 ,n)  = 1
! note the 0-based index            
            bas(6  ,n)  = off 
            do i=1,ngf(n)
                mpgf=kpgfst(n)+i-1
                env(off + i) = gfe(mpgf)       
            end do 
            off = off + ngf(n)
            bas(7,n) = off
            do i=1,ngf(n)
             env(off + i) = gfc(kpgfst(n)+i-1)*CINTgto_norm(bas(2,n), env(bas(6,n)+i))
            end do
            off = off + ngf(n)
        end do        
!!!!!!!!!!!!!!  1e
        do i=1,nshell
                shls(1) = i-1
                kdi = CINTcgto_cart(i-1, bas)
            do j=1,nshell
                shls(2) = j-1
                kdj = CINTcgto_cart(j-1, bas)
                allocate (buf1e(kdi,kdj))
                do ibas=kbasst(i),kbased(i)                          
                    do jbas=kbasst(j),kbased(j)
                        ibu2ei=ibas-kbasst(i)+1
                        ibu2ej=jbas-kbasst(j)+1
                        call cint1e_ovlp_cart(buf1e, shls, atm, natom, bas, nshell, env, opt) 
                        S(ibas,jbas)=buf1e(ibu2ei,ibu2ej)
                        S(jbas,ibas)=S(ibas,jbas)
                        call cint1e_kin_cart(buf1e, shls, atm, natom, bas, nshell, env, opt) 
                        T(ibas,jbas)=buf1e(ibu2ei,ibu2ej)
                        T(jbas,ibas)=T(ibas,jbas)
                        call cint1e_nuc_cart(buf1e, shls, atm, natom, bas, nshell, env, opt) 
                        V(ibas,jbas)=buf1e(ibu2ei,ibu2ej)
                        V(jbas,ibas)=V(ibas,jbas)
                    end do
                end do            
                deallocate (buf1e)
            end do
        end do        
        Hcore=T+V
!!!!!!!!!!!!!!! 2e        
        nshelldone=0
        call cint2e_cart_optimizer(opt,atm,natom,bas,nshell,env)
        do i=1,nshell
            nshelldone=nshelldone+1        
            if(mod(nshelldone,5)==0)write(*,*)real(nshelldone)/nshell*100d0,'% ERI done'     
                shls(1) = i-1
                kdi = CINTcgto_cart(i-1, bas)
            do j=1,nshell
                shls(2) = j-1
                kdj = CINTcgto_cart(j-1, bas)
                do k=1,nshell
                    shls(3) = k-1
                    kdk = CINTcgto_cart(k-1, bas)
                    do l=1,nshell
                        shls(4) = l-1
                        kdl = CINTcgto_cart(l-1, bas)     
                        allocate (buf2e(kdi,kdj,kdk,kdl))
                        call cint2e_cart(buf2e, shls, atm, natom, bas, nshell, env, opt)
                        do ibas=kbasst(i),kbased(i)             
                            do jbas=kbasst(j),kbased(j)
                                do kbas=kbasst(k),kbased(k)
                                    do lbas=kbasst(l),kbased(l)
                                ijkl=indexeri(ibas,jbas,kbas,lbas)
                                ibu2ei=ibas-kbasst(i)+1
                                ibu2ej=jbas-kbasst(j)+1
                                ibu2ek=kbas-kbasst(k)+1
                                ibu2el=lbas-kbasst(l)+1
                                ERIall(ijkl)=buf2e(ibu2ei,ibu2ej,ibu2ek,ibu2el)
                                    end do
                                end do
                            end do
                        end do
                        deallocate (buf2e)
                    end do 
                end do 
            end do 
        end do      
        deallocate (atm, bas, env)     
!!!!!        
    end subroutine   CINT_ERI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    End  Module Int2e

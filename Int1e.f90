!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                          int1e.f90                             !!!!
!!!!      This file calculats the general one electron intgral.     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!! 1. STVINT(nshell,nbasis,natom,kbasst,kbased,knbas,kpgfst,
!!            kpgfed,naxyz, S,T,V,MethodS,MethodT,MethodV)
!description: subroutine
!  Calculate the general overlap intgral, Kinetic intgral and Nuclear 
!   attraction  intgral  of shell.    
!   input:   (1) nshell        integer
!               (Number of shells)
!            (2)nbasis        integer
!               (Number of basis)
!            (3) natom    integer
!               (Number of atoms)
!            (4) kbasst     integer(nshell)
!             (The  first number of basis function in one shell)
!            (5) kbased    integer(nshell)
!             (The  last number of basis function in one shell)
!            (6) knbas    integer(nshell)
!             (Number of basis function in one shell)
!            (7) kpgfst    integer(nshell)
!             (The  first number of primitive gaussian function in one shell)
!            (8) kpgfed    integer(nshell)
!             (The  last number of primitive gaussian function in one shell)
!            (9) naxyz    integer(3*nbasis)
!             (Type of basis function,x in [1,nbasis],y in[1+nbasis,2nbasis]
!                   ,z in[1+2nbasis,3nbasis]    )
!   output:  (1) S     double precision
!            (2) T     double precision  
!            (3) V     double precision  
!            (4) MethodS     integer (Describe the type of used method )
!            (5) MethodT     integer (Describe the type of used method )
!            (6) MethodV     integer (Describe the type of used method )
!2. GHINTtest(nax,nbx,zetap,cp,cA,cb,result) (Abandoned now)
!description: subroutine
! Test for the Gauss–Hermite numerical integration.
!3. GHINTS(nax,nbx,zetap,cp,cA,cb,result)
!description: subroutine
!  Calculat the Gauss–Hermite numerical integration.
!   input:   (1) nax      integer
!            (2) nbx      integer
!            (3) zetap    integer
!            (4) cp    double precision
!                P_i ,i={x,y,z}
!            (5) cA    double precision
!                A_i，i={x,y,z}
!               ngfsum=SUM(ngf(:))
!            (6) cB    double precision
!                B_i，i={x,y,z}
!   output:  (1) result    double precision 
!4. recNAI1(m,nax,nay,naz,nbx,nby,nbz,&
!        ep1,ABK,Ax,Bx,Px,Cx,Ay,By,&
!        Py,Cy,Az,Bz,Pz,Cz,T) 
!description: recursive function
!  Calculat the NAI <a|1/rC|b>.
!  use DPK method
!  ref:F W.Chen,Computational Methods in Quantum Chemistry
!   input:   (1) m   integer
!               (m=0)
!            (2) nax,nay,naz      integer
!            (3) nbx,nby,nbz      integer
!               (2),(3) gauss function parameter
!            (4) ep1              double precision
!                (ep1=1/(ea+eb)) 
!            (5) ABK              double precision
!                 ABK=K_{AB}  
!            (6) Ax,Ay,Az         double precision
!            (7) Bx,By,Bz         double precision
!               (6),(7) Coordinates of center a,b
!            (8) Cx,Cy,Cz         double precision
!                Coordinates of nuclear C
!            (9) Px,Py,Pz         double precision
!                P is the new center after combining A,B.
!            (10) T               double precision
!                   T=(ea+eb)*((Px-Cx)^2+(Py-Cy)^2+(Pz-Cz)^2)
!   output:  (1) result           double precision 
!5. recNAI2(m,nax,nay,naz,nbx,nby,nbz,&
!        ep1,ABK,Ax,Bx,Px,Cx,Ay,By,&
!        Py,Cy,Az,Bz,Pz,Cz,T) 
!description: recursive function
!  Calculat the NAI <a|1/rC|b>.
!  use  S.Obara;A.saika method 
!  ref:The Journal of Chemical Physics 84, 3963 (1986); doi: 10.1063/1.450106
!   input:   (1) m   integer
!               (m=0)
!            (2) nax,nay,naz      integer
!            (3) nbx,nby,nbz      integer
!               (2),(3) gauss function parameter
!            (4) ep1              double precision
!                (ep1=1/(ea+eb)) 
!            (5) ABK              double precision
!                 ABK=K_{AB}  
!            (6) Ax,Ay,Az         double precision
!            (7) Bx,By,Bz         double precision
!               (6),(7) Coordinates of center a,b
!            (8) Cx,Cy,Cz         double precision
!                Coordinates of nuclear C
!            (9) Px,Py,Pz         double precision
!                P is the new center after combining A,B.
!            (10) T               double precision
!                   T=(ea+eb)*((Px-Cx)^2+(Py-Cy)^2+(Pz-Cz)^2)
!   output:  (1) result           double precision 
!6. anaNAI1(nax,nay,naz,nbx,nby,nbz,&
!        ep1,ABK,Ax,Bx,Px,Cx,Ay,By,&
!        Py,Cy,Az,Bz,Pz,Cz,T) 
!description: recursive function
!  Calculat the NAI <a|1/rC|b>.
!  use (Analytical method) O.Ohata method 
!  ref: O.Ohata method (1966): https://doi.org/10.1143/JPSJ.21.2313
!   input:   
!            (1) nax,nay,naz      integer
!            (2) nbx,nby,nbz      integer
!               (1),(2) gauss function parameter
!            (3) ep1              double precision
!                (ep1=1/(ea+eb)) 
!            (4) ABK              double precision
!                 ABK=K_{AB}  
!            (5) Ax,Ay,Az         double precision
!            (6) Bx,By,Bz         double precision
!               (5),(6) Coordinates of center a,b
!            (7) Cx,Cy,Cz         double precision
!                Coordinates of nuclear C
!            (8) Px,Py,Pz         double precision
!                P is the new center after combining A,B.
!            (9) T               double precision
!                   T=(ea+eb)*((Px-Cx)^2+(Py-Cy)^2+(Pz-Cz)^2)
!   output:  (1) result           double precision 
!7. recNAI1fast,recNAI2fast
!description:  function
!   recNAI1fast is the efficient version of recNAI1.
!   recNAI2fast is the efficient version of recNAI2.
!8. recNAI1xyz,recNAI2xyz
!description:  function
!   recNAI1xyz is a part of recNAI1fast.
!   recNAI2xyz is a part of recNAI2fast.
!======================================================================!
  Module Int1e
  contains
    subroutine STVINT(nshell,nbasis,natom,coor,katom,&
            kcenter,gfc,gfe,ntype,ngf,ngfall,&     
            kbasst,kbased,knbas,kpgfst,&
            kpgfed,naxyz,S,T,V,hcore,dipx,dipy,dipz,&
            MethodS,MethodT,MethodV)
        use Global
        implicit double precision (A-H,O-Z)
        dimension::kbasst(nshell),kbased(nshell),knbas(nshell),&
                 kpgfst(nshell),kpgfed(nshell),naxyz(3*nbasis),&
                 S(nbasis,nbasis),T(nbasis,nbasis),V(nbasis,nbasis),&
                 Hcore(nbasis,nbasis),coor(natom,3),katom(natom),&
                 dipx(nbasis,nbasis),dipy(nbasis,nbasis),dipz(nbasis,nbasis)
        dimension::kcenter(Natom),gfc(ngfall),gfe(ngfall),ntype(nshell),ngf(nshell)     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!=====================================================================!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  SHELL  !!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        write(*,*)'*** Computing one electron intgral ***'
        num=0
 !!!!!!       
        do Ish=1,nshell   ! I shell
            Iatom=kcenter(Ish)
            cix=coor(Iatom,1)
            ciy=coor(Iatom,2)
            ciz=coor(Iatom,3)
            do Jsh=1,Ish  ! J shell  
                Jatom=kcenter(Jsh)
                cJx=coor(Jatom,1)
                cJy=coor(Jatom,2)
                cJz=coor(Jatom,3)
                Rij=(cix-cjx)**2+(ciy-cjy)**2+(ciz-cjz)**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! BASIS FUNCTION  !!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                do ibas=kbasst(ish),kbased(ish)
                    nax=naxyz(ibas)
                    nay=naxyz(ibas+nbasis)
                    naz=naxyz(ibas+2*nbasis)
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
!!!!!!!!!                    
                    jbased=kbased(jsh)                       
                    ! if(ish==jsh)   jbased=ibas  
                    do jbas=kbasst(jsh),jbased
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
!!! primitive gauss functions  <a|
                        sij=0d0
                        sijn=0d0
                        tijn=0d0
                        sij2=0d0
                        sij2n=0d0
                        vij=0d0
                        vijn=0d0  
                        dipijx=0d0  
                        dipijy=0d0  
                        dipijz=0d0               
                        do ipgf=ipgfst,ipgfed     
                                ea=gfe(ipgf)
                                ca=gfc(ipgf)
!!! primitive gauss functions |b>
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
!!!!!!! Overlap
!!!!!!!!!!!!!!!!!!!!!!!!!! <a|b> 
!!!!!!! normalization
                    ! dNa=(2*ea/pi)**0.75*((4*ea)**(naxyz(ibas)+naxyz(ibas+nbasis)&
                    !             +naxyz(ibas+2*nbasis))/(Dfact(2*naxyz(ibas)-1)&
                    !             *Dfact(2*naxyz(ibas+nbasis)-1)&
                    !             *Dfact(2*naxyz(ibas+2*nbasis)-1)))**0.5
                    ! dNb=(2*eb/pi)**0.75*((4*eb)**(naxyz(jbas)+naxyz(jbas+nbasis)&
                    !             +naxyz(jbas+2*nbasis))/(Dfact(2*naxyz(jbas)-1)&
                    !             *Dfact(2*naxyz(jbas+nbasis)-1)&
                    !             *Dfact(2*naxyz(jbas+2*nbasis)-1)))**0.5
                    ! dNa=(2*ea/pi)**0.75*((4*ea)**(nax+nay&
                    !             +naz)/(Dfact(2*nax-1)&
                    !             *Dfact(2*nay-1)&
                    !             *Dfact(2*naz-1)))**0.5
                    ! dNb=(2*eb/pi)**0.75*((4*eb)**(nbx+nby&
                    !             +nbz)/(Dfact(2*nbx-1)&
                    !             *Dfact(2*nby-1)&
                    !             *Dfact(2*nbz-1)))**0.5  
                   dna=dnormgf(ea,nax,nay,naz)
                   dnb=dnormgf(eb,nbx,nby,nbz)
!!!! 
                    call GHINTS(nax,nbx,ep,PX,cix,cjx,rstx)
                    call GHINTS(nay,nby,ep,PY,ciy,cjy,rsty)
                    call GHINTS(naz,nbz,ep,PZ,ciz,cjz,rstz)
                    if(MethodS==1)then                                                    
!!!!!!!! 1. Gauss–Hermite     
!!!!! T needs rstx,rsty,rstz.                                             
                        sab=rstx*rsty*rstz*ABK*ep1**1.5
                    else  if(MethodS==2)then 
!!!!!!!! 2. Analytical method
                        sx=0d0
                        sy=0d0
                        sz=0d0
                        do ii=0,floor((nax+nbx)/2d0)
                            sx=sx+fi(2*ii,nax,nbx,Px-cix,Px-cjx)*ep**(-ii-0.5d0)*gamma(ii+0.5d0)                     
                        end do
                        do jj=0,floor((nay+nby)*0.5d0)
                            sy=sy+fi(2*jj,nay,nby,Py-ciy,Py-cjy)*ep**(-jj-0.5d0)*gamma(jj+0.5d0)                                             
                        end do
                        do kk=0,floor((naz+nbz)*0.5d0)
                            sz=sz+fi(2*kk,naz,nbz,Pz-ciz,Pz-cjz)*ep**(-kk-0.5d0)*gamma(kk+0.5d0)
                        end do
                        sab=sx*sy*sz*ABK  !!!  use gamma function
                    else 
                            stop "bad Overlap method"    
                    end if 
!============
!!!!!!!   S(i,j) (normalization)
                    sijn=sijn+ca*cb*sab*dNa*dNb   
!!!!!!!!!!! Kinetic energy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    if(methodT==1)then
!=============== 1. direct ===================================!
!===============  +2,0,-2 
!!!!!!!!!!!!!!!!!!!!!!!!  <a|b+2>   
                        call GHINTS(nax,nbx+2,ep,PX,cix,cjx,rstx2)
                        call GHINTS(nay,nby+2,ep,PY,ciy,cjy,rsty2)
                        call GHINTS(naz,nbz+2,ep,PZ,ciz,cjz,rstz2)                        
                        sab2x=rstx2*rsty*rstz*ABK*ep1**1.5
                        sab2y=rstx*rsty2*rstz*ABK*ep1**1.5
                        sab2z=rstx*rsty*rstz2*ABK*ep1**1.5                   
! !!!!!!!!!!!!!!    
! !!!!!!!!!!!!!!!!!!!!!!!!  <a|b-2> 
                        if(nbx-2<0) then
                            sab_2x=0
                        else
                            call GHINTS(nax,nbx-2,ep,PX,cix,cjx,rstx_2)
                            sab_2x=rstx_2*rsty*rstz*ABK*ep1**1.5
                        end if
                        if(nby-2<0) then
                            sab_2y=0
                        else
                            call GHINTS(nay,nby-2,ep,Py,ciy,cjy,rsty_2)
                            sab_2y=rstx*rsty_2*rstz*ABK*ep1**1.5
                        end if                   
                        if(nbz-2<0) then 
                            sab_2z=0
                        else
                            call GHINTS(naz,nbz-2,ep,Pz,ciz,cjz,rstz_2)
                            sab_2z=rstx*rsty*rstz_2*ABK*ep1**1.5
                        end if
! !!!!!!           
! !!!!!!!!  sum
                        tix=-0.5d0*nbx*(nbx-1)*sab_2x+eb*(2*nbx+1)*sab-2*eb**2*sab2x
                        tiy=-0.5d0*nby*(nby-1)*sab_2y+eb*(2*nby+1)*sab-2*eb**2*sab2y
                        tiz=-0.5d0*nbz*(nbz-1)*sab_2z+eb*(2*nbz+1)*sab-2*eb**2*sab2z
! !!!!!!!
! !!!!!!!!!!!!! T(i,j)
                        tijn=tijn+(tix+tiy+tiz)*ca*cb*dNa*dNb   
!!!!
                    else if(methodT==2)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                 
!=====================2.integral by parts ========================!
!===============  -1-1,-1+1,+1-1,+1+1
!!!!!!!!!!!!!!!!!!!!!!!!  <a-1|b-1>  
                        if(nax-1<0.or.nbx-1<0) then
                            sab_1_1x=0
                        else
                            call GHINTS(nax-1,nbx-1,ep,PX,cix,cjx,rstx_1_1)
                            sab_1_1x=rstx_1_1*rsty*rstz*ABK*ep1**1.5
                        end if
                        if(nay-1<0.or.nby-1<0) then
                            sab_1_1y=0
                        else
                            call GHINTS(nay-1,nby-1,ep,Py,ciy,cjy,rsty_1_1)
                            sab_1_1y=rsty_1_1*rstx*rstz*ABK*ep1**1.5
                        end if
                        if(naz-1<0.or.nbz-1<0) then
                            sab_1_1z=0
                        else
                            call GHINTS(naz-1,nbz-1,ep,Pz,ciz,cjz,rstz_1_1)
                            sab_1_1z=rstz_1_1*rsty*rstx*ABK*ep1**1.5
                        end if                
!!!!!!!!!!!!!!!!!!!!!!!!  <a-1|b+1>   
                        if(nax-1<0) then
                            sab_11x=0
                        else
                            call GHINTS(nax-1,nbx+1,ep,PX,cix,cjx,rstx_11)
                            sab_11x=rstx_11*rsty*rstz*ABK*ep1**1.5
                        end if
                        if(nay-1<0.) then
                            sab_11y=0
                        else
                            call GHINTS(nay-1,nby+1,ep,Py,ciy,cjy,rsty_11)
                            sab_11y=rsty_11*rstx*rstz*ABK*ep1**1.5
                        end if
                        if(naz-1<0) then
                            sab_11z=0
                        else
                            call GHINTS(naz-1,nbz+1,ep,Pz,ciz,cjz,rstz_11)
                            sab_11z=rstz_11*rsty*rstx*ABK*ep1**1.5
                        end if                
!!!!!!!!!!!!!!!!!!!!!!!!  <a+1|b-1>   
                        if(nbx-1<0) then
                            sab1_1x=0
                        else
                            call GHINTS(nax+1,nbx-1,ep,PX,cix,cjx,rstx1_1)
                            sab1_1x=rstx1_1*rsty*rstz*ABK*ep1**1.5
                        end if
                        if(nby-1<0) then
                            sab1_1y=0
                        else
                            call GHINTS(nay+1,nby-1,ep,Py,ciy,cjy,rsty1_1)
                            sab1_1y=rsty1_1*rstx*rstz*ABK*ep1**1.5
                        end if
                        if(nbz-1<0) then
                            sab1_1z=0
                        else
                            call GHINTS(naz+1,nbz-1,ep,Pz,ciz,cjz,rstz1_1)
                            sab1_1z=rstz1_1*rsty*rstx*ABK*ep1**1.5
                        end if                
!!!!!!!!!!!!!!!!!!!!!!!!  <a+1|b+1>   
                        call GHINTS(nax+1,nbx+1,ep,PX,cix,cjx,rstx11)
                        call GHINTS(nay+1,nby+1,ep,PY,ciy,cjy,rsty11)
                        call GHINTS(naz+1,nbz+1,ep,PZ,ciz,cjz,rstz11)                        
                        sab11x=rstx11*rsty*rstz*ABK*ep1**1.5
                        sab11y=rstx*rsty11*rstz*ABK*ep1**1.5
                        sab11z=rstx*rsty*rstz11*ABK*ep1**1.5                   
!!!  sum                   
                        tix=nax*nbx*sab_1_1x-2*eb*nax*sab_11x-2*ea*nbx*sab1_1x+4*ea*eb*sab11x
                        tiy=nay*nby*sab_1_1y-2*eb*nay*sab_11y-2*ea*nby*sab1_1y+4*ea*eb*sab11y
                        tiz=naz*nbz*sab_1_1z-2*eb*naz*sab_11z-2*ea*nbz*sab1_1z+4*ea*eb*sab11z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    
!!!!!!!!!!! T(i,j) for methodT=2 ,which is diff from 1 for factor 1/2
                        tijn=tijn+(tix+tiy+tiz)*ca*cb*dNa*dNb*0.5
                    else 
                        stop "bad Kinetic method"                    
                    end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    
!!!!!! NUCLEAR ATTRACTION
!============== <a|1/rc|b> ==============! vab
!!!!!!! 1.(Analytical method) O.Ohata method (1966): https://doi.org/10.1143/JPSJ.21.2313
!!!!!!!!!!! Now this method is transferred to the function anaNAI1.!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Direct computing don't be used anymore,but the code remains for reference.!!!!
!!!!!  Basic information
!                     MethodV=1
!                     vab=0d0
!                     vabn=0d0
!                     do  kc=1,natom                     
!                         CX=coor(kc,1)
!                         CY=coor(kc,2)
!                         CZ=coor(kc,3)
!                         NZ=katom(kc)
!                         TT=ep*((PX-CX)**2+(PY-CY)**2+(PZ-CZ)**2)
!                         if(TT<1d-012) TT=0d0
! !!!!!!!!!!!!!!!!!!!!!!!!         
!                         yip=0.25d0*ep1
!                         vxyz=0d0
!     !!!  x                          
!             do ii=0,nax+nbx 
!                 do ir=0,floor(ii*0.5d0)
!                     do iu=0,floor((ii-2*ir)*0.5d0)

!                         vx=(-1)**ii*fi(ii,nax,nbx,Px-cix,Px-cjx)*(-1)**iu&
!                             *ifact(ii)*(PX-CX)**(ii-2*ir-2*iu)*yip**(ir+iu)/&
!                             (ifact(ir)*ifact(iu)*ifact(ii-2*ir-2*iu)  )              
!     !!!  y   
!                         do jj=0,nay+nby 
!                             do is=0,floor(jj*0.5d0)
!                                 do iv=0,floor((jj-2*is)*0.5d0)

!                                     vy=(-1)**jj*fi(jj,nay,nby,Py-ciy,Py-cjy)*(-1)**iv&
!                                       *ifact(jj)*(Py-Cy)**(jj-2*is-2*iv)*yip**(is+iv)/&
!                                       (ifact(is)*ifact(iv)*ifact(jj-2*is-2*iv))      
!     !!!  z              
!                                     do kk=0,naz+nbz 
!                                         do it=0,floor(kk*0.5d0)
!                                             do iw=0,floor((kk-2*it)*0.5d0)

           
!                                                 vz=(-1)**kk*fi(kk,naz,nbz,Pz-ciz,Pz-cjz)*(-1)**iw&
!                                                 *ifact(kk)*(Pz-Cz)**(kk-2*it-2*iw)*yip**(it+iw)/&
!                                                 (ifact(it)*ifact(iw)*ifact(kk-2*it-2*iw))                                                     
!                                     !!!!
!                                                 mu=ii+jj+kk-2*(ir+is+it)-(iu+iv+iw)
!                                                 ! vxyz=vxyz+vx*vy*vz*f(mu,TT)
!                                                 vxyz=vxyz+vx*vy*vz*F_pade(mu,TT)

!                                             end do  !iw
!                                         end do  !it
!                                     end do  !kk             
!                                 end do  !iv
!                             end do  !is 
!                         end do  !jj                                      
!                     end do !iu
!                 end do  !ir
!             end do  !ii
!                                 vabn=vabn+vxyz*nz*2*pi*ep1*ABK*dNa*dNb
!                     end do  ! natom         
!!!=====================================================================!!!!!
!!!!!!! Compute NAI by using functions (Mainly used now)
!!!!!!! 1.(Analytical method) O.Ohata method (1966): https://doi.org/10.1143/JPSJ.21.2313
!!!!!!! 2.Recursion (Recursive function method) Maybe low efficiency.
!!!!!!!!!! 
                    naxx=nax
                    nbxx=nbx
                    Axx=cix
                    Bxx=cjx                    
                    nayy=nay
                    nbyy=nby
                    Ayy=ciy
                    Byy=cjy
                    nazz=naz
                    nbzz=nbz
                    Azz=ciz
                    Bzz=cjz  
                    vab=0d0 
                    do  kc=1,natom
                        CX=coor(kc,1)
                        CY=coor(kc,2)
                        CZ=coor(kc,3)
                        NZ=katom(kc)
                        TT=ep*((PX-CX)**2+(PY-CY)**2+(PZ-CZ)**2)
                        ! if(abs(TT)<1d-12) TT=0d0
!!!!
                        if (methodV==1)then
!!!! Method 1 
!!! Ohata method          
                        vab=vab+anaNAI1(naxx,nayy,nazz,nbxx,nbyy,nbzz,&
                                    ep1,ABK,Axx,Bxx,Px,Cx,Ayy,Byy,Py,Cy,&
                                    Azz,Bzz,Pz,Cz,TT)*nz*2*pi*ep1*ABK
                        else if (methodV==2)then
!!!! Method 2     
!!! DRK  method            
                            ! vab=vab+nz*recNAI1(0,naxx,nayy,nazz,nbxx,nbyy,nbzz,&
                            !         ep1,ABK,Axx,Bxx,Px,Cx,Ayy,Byy,Py,Cy,&
                            !         Azz,Bzz,Pz,Cz,TT)  
!!!! FAST CODE
                            vab=vab+nz*recNAI1fast(naxx,nayy,nazz,nbxx,nbyy,nbzz,&
                                    ep1,ABK,Axx,Bxx,Px,Cx,Ayy,Byy,Py,Cy,&
                                    Azz,Bzz,Pz,Cz,TT)                                                                  
                        else if (methodV==3)then                                    
!!!!! Method 3     
!!!  Obara,Saika method                        
                            ! vab=vab+nz*recNAI2(0,naxx,nayy,nazz,nbxx,nbyy,nbzz,&
                            !         ep1,ABK,Axx,Bxx,Px,Cx,Ayy,Byy,Py,Cy,&
                            !         Azz,Bzz,Pz,Cz,TT)  
!!!!!!!!!!!!! FAST CODE
                        vab=vab+nz*recNAI2fast(naxx,nayy,nazz,nbxx,nbyy,nbzz,&
                                    ep1,ABK,Axx,Bxx,Px,Cx,Ayy,Byy,Py,Cy,&
                                    Azz,Bzz,Pz,Cz,TT)                                                                    
                        else 
                            stop "bad NAI method"                                    
                        end if                                           
                    end do
                    vabn=vab*dNa*dNb
!================= <a|1/rc|b> -> sum(ca*cb*<a|1/rc|b>)==============!
!!  V(i,j)                   
                    vijn=vijn+vabn*ca*cb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! Multipole matrices
!!! dipole moment
                    call GHINTS(nax+1,nbx,ep,PX,cix,cjx,drstx)
                    call GHINTS(nay+1,nby,ep,PY,ciy,cjy,drsty)
                    call GHINTS(naz+1,nbz,ep,PZ,ciz,cjz,drstz) 
!!!! x                                                                                                                 
                    sa1bx=drstx*rsty*rstz*ABK*ep1**1.5
                    dipxab=sab*cix+sa1bx
                    dipijx=dipijx+dipxab*dNa*dNb*ca*cb
!!!! y                    
                    sa1by=rstx*drsty*rstz*ABK*ep1**1.5
                    dipyab=sab*ciy+sa1by
                    dipijy=dipijy+dipyab*dNa*dNb*ca*cb  
!!!! z                    
                    sa1bz=rstx*rsty*drstz*ABK*ep1**1.5
                    dipzab=sab*ciz+sa1bz
                    dipijz=dipijz+dipzab*dNa*dNb*ca*cb                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            end do                          ! b pgf
                        end do                              ! a pgf
!!!!!!!!!!!!  Matrix assignment                    
                        if(abs(sijn)<1d-16) sijn=0d0
                        S(ibas,jbas)=sijn     
                        S(jbas,ibas)=sijn  
                        if(abs(tijn)<1d-16) tijn=0d0
                        T(ibas,jbas)=tijn
                        T(jbas,ibas)=tijn
                        if(abs(vijn)<1d-16) vijn=0d0                       
                        V(ibas,jbas)=vijn
                        V(jbas,ibas)=vijn
                        dipx(ibas,jbas)=dipijx
                        dipx(jbas,ibas)=dipijx
                        dipy(ibas,jbas)=dipijy
                        dipy(jbas,ibas)=dipijy
                        dipz(ibas,jbas)=dipijz
                        dipz(jbas,ibas)=dipijz                                                
!!!!!!!!!!!
                    end do  ! J basis             
                end do    ! I basis  
!!!!!!!!!!!
            end do      ! J shell
        end do         ! I shell
        Hcore=T-V
        write(*,*)'*** One electron intgral Done***'        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end  subroutine
!======================================================================!
    subroutine GHINTtest(npoint,resul)
!!! $\int_{-\infity}^{+\infity}e^(-x^2)*x^2 dx$
!!! Only test.
        use Global ,only :pi
        use gausshermite
        implicit double precision (A-H,O-Z)
        temp=0d0
        if (npoint.eq.5) then
            write(*,*)'test2'
            do ii=1,npoint
                temp=temp+gaussHerw5(ii)*gaussHerx5(ii)**2
            end do
        else if (npoint<=1) then
            result=sqrt(pi)
            return        
        else   
            write(*,*)'test1'
            do ii=1,npoint
                temp=temp+gaussHerw(ii)*gaussHerx(ii)**2
            end do
        end if
        resul=temp
    end  subroutine
!======================================================================!
    subroutine GHINTS(nax,nbx,zetap,cp,cA,cb,result)
! $\int_{-\infity}^{+\infity}e^(-x^2)*(x/sprt(p)+px-ax)^nax*
! (x/sprt(p)+px-bx)^nbx dx$
        use Global ,only :pi
        use gausshermite
        implicit double precision (A-H,O-Z)
        if (nax+nbx==0)then
            result=sqrt(pi)
            return
        end if
        p1=1d0/sqrt(zetap)
        ! npoint=floor((nax+nbx)/2d0)+1
        npoint=ceiling((nax+nbx+1)/2d0)
        npointst=npoint*(npoint-1)/2+1
        npointed=npointst+npoint-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        temp=0d0
        do ii=npointst,npointed
            temp=temp+gauHerw(ii)*(gauHerx(ii)*p1+cp-ca)**nax &
                    *(gauHerx(ii)*p1+cp-cb)**nbx
        end do
        result=temp
    end  subroutine
!======================================================================!
!======================================================================!
    recursive function recNAI1(m,nax,nay,naz,nbx,nby,nbz,&
                                ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,&
                                Az,Bz,PPz,Cz,TT)  result(tem)
!!!!!                                
        use Global,only:PI,F,F_pade
        implicit double precision (A-H,O-Z)    
        tem=0d0
        if(nax<0.or.nbx<0.or.nay<0.or.nby<0.or.naz<0.or.nbz<0)then
            tem=0d0
            return
        end if
!!!!!!!!! recursive method 1 
!!!! DPK method
!!!! F W.Chen,Computational Methods in Quantum Chemistry
    MethodV=2
        if (nax==0.and.nbx==0)then
            if (nay==0.and.nby==0)then
                if (naz==0.and.nbz==0)then
                    tem=2*pi*ep1*ABK*F(m,TT)
                    ! tem=2*pi*ep1*ABK*F_pade(m,TT)                    
                    return
                else  if(nbz==0) then
                tem=(naz-1)*ep1/2d0*recNAI1(m,nax,nay,naz-2,nbx,nby,0,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                    +(ppz-az)*recNAI1(m,nax,nay,naz-1,nbx,nby,0,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                    -(ppz-az)*recNAI1(m+1,nax,nay,naz-1,nbx,nby,0,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                    -(naz-1)*ep1/2d0*recNAI1(m+1,nax,nay,naz-2,nbx,nby,0,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                    +(Cz-Az)*recNAI1(m+1,nax,nay,naz-1,nbx,nby,0,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) 
                    return
                 else           
                    tem=recNAI1(m,nax,nay,naz+1,nbx,nby,nbz-1,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                        +(Az-Bz)*recNAI1(m,nax,nay,naz,nbx,nby,nbz-1,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) 
                    return
                end if
            else  if(nby==0) then
                tem=(nay-1)*ep1/2d0*recNAI1(m,nax,nay-2,naz,nbx,0,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                    +(ppy-ay)*recNAI1(m,nax,nay-1,naz,nbx,0,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                    -(ppy-ay)*recNAI1(m+1,nax,nay-1,naz,nbx,0,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                    -(nay-1)*ep1/2d0*recNAI1(m+1,nax,nay-2,naz,nbx,0,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                    +(Cy-Ay)*recNAI1(m+1,nax,nay-1,naz,nbx,0,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT)                     
                    return
            else 
                tem=recNAI1(m,nax,nay+1,naz,nbx,nby-1,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                    +(Ay-By)*recNAI1(m,nax,nay,naz,nbx,nby-1,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) 
                return
            end if
        else  if(nbx==0) then
            tem=(nax-1)*ep1/2d0*recNAI1(m,nax-2,nay,naz,0,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                +(ppx-ax)*recNAI1(m,nax-1,nay,naz,0,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                -(ppx-ax)*recNAI1(m+1,nax-1,nay,naz,0,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                -(nax-1)*ep1/2d0*recNAI1(m+1,nax-2,nay,naz,0,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                +(Cx-Ax)*recNAI1(m+1,nax-1,nay,naz,0,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) 
            return
        else 
            tem=recNAI1(m,nax+1,nay,naz,nbx-1,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
                +(Ax-Bx)*recNAI1(m,nax,nay,naz,nbx-1,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) 
            return
        end if
    end function recNAI1
!======================================================================!
!======================================================================!
    recursive function recNAI2(m,naxx,nayy,nazz,nbxx,nbyy,nbzz,&
                                ep11,ABK,Axx,Bxx,Px,Cxx,Ayy,Byy,Py,Cyy,&
                                Azz,Bzz,Pz,Czz,TT)  result(tem)
!!!!!
        use Global,only:PI,F,F_pade,F_ser,dexchange,iexchange
        implicit double precision (A-H,O-Z)
!!!! copy values for exchange
        nax=naxx
        nay=nayy
        naz=nazz
        nbx=nbxx
        nby=nbyy
        nbz=nbzz
        ep1=ep11
        ax=Axx
        bx=Bxx
        cx=Cxx        
        ppx=Px
        ay=Ayy
        by=Byy
        Cy=Cyy
        ppy=Py
        az=Azz
        bz=Bzz
        Cz=Czz
        ppz=Pz        
!!!!!!            
        tem=0d0
        if(nax<0.or.nbx<0.or.nay<0.or.nby<0.or.naz<0.or.nbz<0)then
            tem=0d0
            return
        end if
!!!!!!!!!!!!! Recursive method 2   
!!!!!!!       S.Obara;A.saika  The Journal of Chemical Physics 84, 3963 (1986); doi: 10.1063/1.450106
        MethodV=3
        if (nax==0.and.nbx==0)then
            if (nay==0.and.nby==0)then
                if (naz==0.and.nbz==0)then
                    tem=2*ABK*pi*ep1*F(m,TT)
                    ! tem=2*ABK*pi*ep1*F_ser(m,TT)
                    ! tem=2*ABK*pi*ep1*F_pade(m,TT)
                    return
                else 
!!!! make sure naz>=nbz
                if(naz<nbz)then
                    ! nchtem=naz
                    ! naz=nbz
                    ! nbz=nchtem
                    ! chtem=Az
                    ! Az=Bz
                    ! Bz=chtem
!!!!!!!! Use function Iexchange/Dexchange
                    call Iexchange(naz,nbz)                              
                    call dexchange(Az,Bz)                     
                end if                 
!!!! z
        tem=(Ppz-Az)*recNAI2(m,nax,nay,naz-1,nbx,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
            -(Ppz-Cz)*recNAI2(m+1,nax,nay,naz-1,nbx,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
         +(naz-1)*ep1*0.5d0*(recNAI2(m,nax,nay,naz-2,nbx,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
            -recNAI2(m+1,nax,nay,naz-2,nbx,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) )&
         +(nbz)*ep1*0.5d0*(recNAI2(m,nax,nay,naz-1,nbx,nby,nbz-1,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
            -recNAI2(m+1,nax,nay,naz-1,nbx,nby,nbz-1,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) )
                    return
                end if
        else 
!!!! make sure nay>=nby
        if(nay<nby)then
            ! nchtem=nay
            ! nay=nby
            ! nby=nchtem
            ! chtem=Ay
            ! Ay=By
            ! By=chtem
!!!!!!!! Use function Iexchange/Dexchange
        call Iexchange(nay,nby)                              
        call dexchange(Ay,By)             
!!!!   change z
            ! nchtem=naz
            ! naz=nbz
            ! nbz=nchtem
            ! chtem=Az
            ! Az=Bz
            ! Bz=chtem
!!!!!!!! Use function Iexchange/Dexchange
            call Iexchange(naz,nbz)                              
            call dexchange(Az,Bz)             
        end if
!!!! y
        tem=(Ppy-Ay)*recNAI2(m,nax,nay-1,naz,nbx,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
            -(Ppy-Cy)*recNAI2(m+1,nax,nay-1,naz,nbx,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
            +(nay-1)*ep1*0.5d0*(recNAI2(m,nax,nay-2,naz,nbx,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
            -recNAI2(m+1,nax,nay-2,naz,nbx,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) )&
            +(nby)*ep1*0.5d0*(recNAI2(m,nax,nay-1,naz,nbx,nby-1,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
            -recNAI2(m+1,nax,nay-1,naz,nbx,nby-1,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) )            
            return
        end if
    else 
!!!! make sure nax>=nbx 
    if(nax<nbx)then
        ! nchtem=nax
        ! nax=nbx
        ! nbx=nchtem
        ! chtem=Ax
        ! Ax=Bx
        ! Bx=chtem
!!!!!!!! Use function Iexchange/Dexchange
        call Iexchange(nax,nbx)                              
        call dexchange(Ax,Bx)         
!!!! change y
        ! nchtem=nay
        ! nay=nby
        ! nby=nchtem
        ! chtem=Ay
        ! Ay=By
        ! By=chtem
!!!!!!!! Use function Iexchange/Dexchange
        call Iexchange(nay,nby)                              
        call dexchange(Ay,By)         
!!!! change z
        ! nchtem=naz
        ! naz=nbz
        ! nbz=nchtem
        ! chtem=Az
        ! Az=Bz
        ! Bz=chtem
!!!!!!!! Use function Iexchange/Dexchange
        call Iexchange(naz,nbz)                              
        call dexchange(Az,Bz)         
    end if
!!!! x 
      tem=(Ppx-Ax)*recNAI2(m,nax-1,nay,naz,nbx,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
        -(Ppx-Cx)*recNAI2(m+1,nax-1,nay,naz,nbx,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
        +(nax-1)*ep1*0.5d0*(recNAI2(m,nax-2,nay,naz,nbx,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
        -recNAI2(m+1,nax-2,nay,naz,nbx,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) )&
        +(nbx)*ep1*0.5d0*(recNAI2(m,nax-1,nay,naz,nbx-1,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) &
        -recNAI2(m+1,nax-1,nay,naz,nbx-1,nby,nbz,ep1,ABK,Ax,Bx,PPx,Cx,Ay,By,PPy,Cy,Az,Bz,PPz,Cz,TT) )
      return
     end if
    end function recNAI2
!======================================================================!
    function anaNAI1(nax,nay,naz,nbx,nby,nbz,ep1,ABK,Ax,Bx,Px,Cx,Ay,By,Py,Cy,&
                    Az,Bz,Pz,Cz,TT)  
!!!!!                                
        use Global
        implicit double precision (A-H,O-Z)    
!!!!!!!!! recursive method 1 
!!!!!!! 1.(Analytical method) O.Ohata method (1966): https://doi.org/10.1143/JPSJ.21.2313
        yip=0.25d0*ep1
        vxyz=0d0
!!!!!  x                          
        do ii=0,nax+nbx 
            do ir=0,floor(ii*0.5d0)
                do iu=0,floor((ii-2*ir)*0.5d0)
!!!!
                    vx=(-1)**ii*fi(ii,nax,nbx,Px-Ax,Px-Bx)*(-1)**iu&
                        *ifact(ii)*(PX-CX)**(ii-2*ir-2*iu)*yip**(ir+iu)/&
                        (ifact(ir)*ifact(iu)*ifact(ii-2*ir-2*iu)  )              
!!!!!  y   
                    do jj=0,nay+nby 
                        do is=0,floor(jj*0.5d0)
                            do iv=0,floor((jj-2*is)*0.5d0)
!!!!
                                vy=(-1)**jj*fi(jj,nay,nby,Py-Ay,Py-By)*(-1)**iv&
                                    *ifact(jj)*(Py-Cy)**(jj-2*is-2*iv)*yip**(is+iv)/&
                                    (ifact(is)*ifact(iv)*ifact(jj-2*is-2*iv))      
!!!!!  z              
                                do kk=0,naz+nbz 
                                    do it=0,floor(kk*0.5d0)
                                        do iw=0,floor((kk-2*it)*0.5d0)
!!!!                                        
                                            vz=(-1)**kk*fi(kk,naz,nbz,Pz-Az,Pz-Bz)*(-1)**iw&
                                            *ifact(kk)*(Pz-Cz)**(kk-2*it-2*iw)*yip**(it+iw)/&
                                            (ifact(it)*ifact(iw)*ifact(kk-2*it-2*iw))                                                     
!!!
                                            mu=ii+jj+kk-2*(ir+is+it)-(iu+iv+iw)
                                            vxyz=vxyz+vx*vy*vz*F(mu,TT)
                                            ! vxyz=vxyz+vx*vy*vz*F_pade(mu,TT)

                                        end do  !iw
                                    end do  !it
                                end do  !kk             
                            end do  !iv
                        end do  !is 
                    end do  !jj                                      
                end do !iu
            end do  !ir
        end do  !ii
        anaNAI1=vxyz
    end function anaNAI1
!======================================================================!
!======================================================================!
    function recNAI1fast (naxx,nayy,nazz,nbxx,nbyy,nbzz,&
             ep11,ABK,Axx,Bxx,Px,Cxx,Ayy,Byy,Py,Cyy,Azz,Bzz,Pz,Czz,TT) 
!!!!!                                
        use Global,only:PI,F,iexchange,dexchange
        implicit double precision (A-H,O-Z)    
        allocatable::temp1(:)
!!!! copy values for exchange
        nax=naxx
        nay=nayy
        naz=nazz
        nbx=nbxx
        nby=nbyy
        nbz=nbzz
        ep1=ep11
        ax=Axx
        bx=Bxx
        cx=Cxx        
        ppx=Px
        ay=Ayy
        by=Byy
        Cy=Cyy
        ppy=Py
        az=Azz
        bz=Bzz
        Cz=Czz
        ppz=Pz         
!!!!        
        if(naxx<0.or.nbxx<0.or.nayy<0.or.nbyy<0.or.nazz<0.or.nbzz<0)then
           stop 'bad input for recNAI1fast'
        end if
        Lsum=nax+nay+naz+nbx+nby+nbz
        allocate(temp1(0:Lsum))  
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
!!!!         
        do i=0,Lsum
            temp1(i)=F(i,TT)
        end do
        call recNAI1xyz(nax,nbx,Lsum,ep1,Ax,Bx,Ppx,Cx,temp1,1)
        call recNAI1xyz(nay,nby,Lsum,ep1,Ay,By,Ppy,Cy,temp1,1)
        call recNAI1xyz(naz,nbz,Lsum,ep1,Az,Bz,Ppz,Cz,temp1,0)
        recNAI1fast=temp1(0)*ABK*2d0*pi*ep1
        deallocate(temp1)
    end function recNAI1fast
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine recNAI1xyz(nax,nbx,Lsum,ep1,Ax,Bx,Px,Cx,Finreout,Notz)
        implicit double precision (A-H,O-Z)    
        dimension::Finreout(0:Lsum),temp1(0:nax+nbx+1,0:Lsum),&
                    temp2(0:nax+nbx+1,0:nbx)
!!!! init
        nx=nax+nbx
        if(Lsum==0 .or. nx==0) return
!!!!
        temp1(0,Lsum)=Finreout(Lsum)
        temp1(1,Lsum)=(Px-Ax)*Finreout(Lsum)
!!!! Lsum=1
        if(Lsum==1) then
            temp1(1,0)=(Px-Ax)*Finreout(0)-(Px-Ax)*Finreout(1)+(Cx-Ax)*Finreout(1)
            Finreout(0:1)=temp1(1,0:1)
            return
        end if
!!!! nax+nbx=1
        if(nx==1)then
            do m=Lsum-1,0,-1
                temp1(0,m)=Finreout(m)
                temp1(1,m)=(Px-Ax)*Finreout(m)-(Px-Ax)*Finreout(1+m)+(Cx-Ax)*Finreout(m+1)
            end do
            if(nbx==0)then
                if(Notz==0)then
                    Finreout(0)=temp1(1,0)
                    return
                end if
                do m=0,Lsum
                    Finreout(m)=temp1(1,m)
                end do
                return
            end if
        end if
!!!! nax+nbx >= 2
        do i=1,nx
            temp1(i+1,Lsum)=i*0.5d0*ep1*temp1(i-1,Lsum)+(Px-Ax)*temp1(i,Lsum)
        end do
        do m=Lsum-1,0,-1
            temp1(0,m)=Finreout(m)
            temp1(1,m)=(Px-Ax)*Finreout(m)-(Px-Ax)*Finreout(1+m)+(Cx-Ax)*Finreout(m+1)
            do i=1,nx
                temp1(i+1,m)=i*0.5d0*ep1*(temp1(i-1,m)-temp1(i-1,m+1))&
                        +(Px-Ax)*temp1(i,m)-(Px-Ax)*temp1(i,m+1)&
                        +(Cx-Ax)*temp1(i,m+1)
            end do
        end do
!!!!    
        if(nbx==0)then
            if(Notz==1)then
                Finreout(0:Lsum)=temp1(nax,0:Lsum)
                return
            else
                Finreout(0)=temp1(nax,0)
                return                
            end if
        end if
!!!!
        if(nbx==1)then
            if(Notz==1)then
                do m=0,Lsum
                    do i=0,nax+2
                        temp2(i,0)=temp1(i,m)
                    end do
                    temp2(nax+1,1)=temp2(nax+2,0)+(Ax-Bx)*temp2(nax+1,0)
                    Finreout(m)=temp2(nax+1,0)+(Ax-Bx)*temp2(nax,0)
                end do        
                return
            else
                do i=0,nax+2
                     temp2(i,0)=temp1(i,0)
                end do 
                temp2(nax+1,1)=temp2(nax+2,0)+(Ax-Bx)*temp2(nax+1,0)
                Finreout(0)=temp2(nax+1,0)+(Ax-Bx)*temp2(nax,0)                
                return                
            end if
        end if
!!!!
        if (Notz==1)then
            do m=0,Lsum
                do i=0,nx+1
                    temp2(i,0)=temp1(i,m)
                end do            
                do  j=1,nbx
                    do  i=nx,nax,-1
                        temp2(i,j)=temp2(i+1,j-1)+(Ax-Bx)*temp2(i,j-1)
                    end do
                end do 
                Finreout(m)=temp2(nax,nbx)
            end do
        else
            do i=0,nx+1
                temp2(i,0)=temp1(i,0)
            end do  
            do  j=1,nbx
                do  i=nx,nax,-1
                    temp2(i,j)=temp2(i+1,j-1)+(Ax-Bx)*temp2(i,j-1)
                end do
            end do           
            Finreout(0)=temp2(nax,nbx) 
        end if       
!!!!!!!!!!!!!!!!!!!!!!!!        
    end subroutine
!======================================================================!
    function recNAI2fast (naxx,nayy,nazz,nbxx,nbyy,nbzz,&
             ep11,ABK,Axx,Bxx,Px,Cxx,Ayy,Byy,Py,Cyy,Azz,Bzz,Pz,Czz,TT) 
!!!!!                                
        use Global,only:PI,F,iexchange,dexchange
        implicit double precision (A-H,O-Z)    
        allocatable::temp1(:)
!!!! copy values for exchange
        nax=naxx
        nay=nayy
        naz=nazz
        nbx=nbxx
        nby=nbyy
        nbz=nbzz
        ep1=ep11
        ax=Axx
        bx=Bxx
        cx=Cxx        
        ppx=Px
        ay=Ayy
        by=Byy
        Cy=Cyy
        ppy=Py
        az=Azz
        bz=Bzz
        Cz=Czz
        ppz=Pz         
!!!!        
        if(naxx<0.or.nbxx<0.or.nayy<0.or.nbyy<0.or.nazz<0.or.nbzz<0)then
           stop 'bad input for recNAI1fast'
        end if
        Lsum=nax+nay+naz+nbx+nby+nbz
        allocate(temp1(0:Lsum))  
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
!!!!         
        do i=0,Lsum
            temp1(i)=F(i,TT)
        end do
        call recNAI2xyz(nax,nbx,Lsum,ep1,Ax,Bx,Ppx,Cx,temp1,1)
        call recNAI2xyz(nay,nby,Lsum,ep1,Ay,By,Ppy,Cy,temp1,1)
        call recNAI2xyz(naz,nbz,Lsum,ep1,Az,Bz,Ppz,Cz,temp1,0)
        recNAI2fast=temp1(0)*ABK*2d0*pi*ep1
        deallocate(temp1)
    end function recNAI2fast
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine recNAI2xyz(nax,nbx,Lsum,ep1,Ax,Bx,Px,Cx,Finreout,Notz)
        implicit double precision (A-H,O-Z)    
        dimension::Finreout(0:Lsum),temp1(0:nax+nbx+1,0:Lsum),&
                    temp2(0:nax+nbx+1,0:nbx)
!!!! init
        nx=nax+nbx
        if(Lsum==0 .or. nx==0) return
!!!!
        temp1(0,Lsum)=Finreout(Lsum)
        temp1(1,Lsum)=(Px-Ax)*Finreout(Lsum)
!!!! Lsum=1
        if(Lsum==1) then
            temp1(1,0)=(Px-Ax)*Finreout(0)-(Px-Cx)*Finreout(1)
            Finreout(0:1)=temp1(1,0:1)
            return
        end if
!!!! nax+nbx=1
        if(nx==1)then
            do m=Lsum-1,0,-1
                temp1(0,m)=Finreout(m)
                temp1(1,m)=(Px-Ax)*Finreout(m)-(Px-Cx)*Finreout(m+1)
            end do
            if(nbx==0)then  !!!! nax>=nbx
                if(Notz==0)then
                    Finreout(0)=temp1(1,0)
                    return
                end if
                do m=0,Lsum
                    Finreout(m)=temp1(1,m)
                end do
                return
            end if
        end if
!!!! nax+nbx >= 2
        do i=1,nx
            temp1(i+1,Lsum)=i*0.5d0*ep1*temp1(i-1,Lsum)+(Px-Ax)*temp1(i,Lsum)
        end do
        do m=Lsum-1,0,-1
            temp1(0,m)=Finreout(m)
            temp1(1,m)=(Px-Ax)*Finreout(m)-(Px-Cx)*Finreout(m+1)
            do i=1,nx
                temp1(i+1,m)=i*0.5d0*ep1*(temp1(i-1,m)-temp1(i-1,m+1))&
                        +(Px-Ax)*temp1(i,m)-(Px-Cx)*temp1(i,m+1)
            end do
        end do
!!!!    
        if(nbx==0)then
            if(Notz==1)then
                Finreout(0:Lsum)=temp1(nax,0:Lsum)
                return
            else
                Finreout(0)=temp1(nax,0)
                return                
            end if
        end if
!!!!
        if(nbx==1)then
            if(Notz==1)then
                do m=0,Lsum
                    do i=0,nax+2
                        temp2(i,0)=temp1(i,m)
                    end do
                    temp2(nax+1,1)=temp2(nax+2,0)+(Ax-Bx)*temp2(nax+1,0)
                    Finreout(m)=temp2(nax+1,0)+(Ax-Bx)*temp2(nax,0)
                end do        
                return
            else
                do i=0,nax+2
                     temp2(i,0)=temp1(i,0)
                end do 
                temp2(nax+1,1)=temp2(nax+2,0)+(Ax-Bx)*temp2(nax+1,0)
                Finreout(0)=temp2(nax+1,0)+(Ax-Bx)*temp2(nax,0)                
                return                
            end if
        end if
!!!! nbx>1
        if (Notz==1)then
            do m=0,Lsum
                do i=0,nx+1
                    temp2(i,0)=temp1(i,m)
                end do            
                do  j=1,nbx
                    do  i=nx,nax,-1
                        temp2(i,j)=temp2(i+1,j-1)+(Ax-Bx)*temp2(i,j-1)
                    end do
                end do 
                Finreout(m)=temp2(nax,nbx)
            end do
        else
            do i=0,nx+1
                temp2(i,0)=temp1(i,0)
            end do  
            do  j=1,nbx
                do  i=nx,nax,-1
                    temp2(i,j)=temp2(i+1,j-1)+(Ax-Bx)*temp2(i,j-1)
                end do
            end do           
            Finreout(0)=temp2(nax,nbx) 
        end if       
!!!!!!!!!!!!!!!!!!!!!!!!        
    end subroutine recNAI2xyz
!======================================================================!
    subroutine auxNAI(natom,coor,katom,&
            nshell,nbasis,kcenter,gfc,gfe,ntype,ngf,ngfall,&
            kbasst,kbased,knbas,kpgfst,kpgfed,naxyz,&
            nshellaux,nbasisaux,kcenteraux,gfcaux,gfeaux,ntypeaux,ngfaux,ngfallaux,&
            kbasstaux,kbasedaux,knbasaux,kpgfstaux,kpgfedaux,naxyzaux,&
            V,MethodV)
        use Global
        implicit double precision (A-H,O-Z)
        dimension::V(nbasisaux,nbasis),coor(natom,3),katom(natom)
        dimension::kbasst(nshell),kbased(nshell),knbas(nshell),&
                 kpgfst(nshell),kpgfed(nshell),naxyz(3*nbasis),&
                 kcenter(Natom),gfc(ngfall),gfe(ngfall),ntype(nshell),ngf(nshell) 
        dimension::kbasstaux(nshellaux),kbasedaux(nshellaux),knbasaux(nshellaux),&
                 kpgfstaux(nshellaux),kpgfedaux(nshellaux),naxyzaux(3*nbasisaux),&
        kcenteraux(Natom),gfcaux(ngfallaux),gfeaux(ngfallaux),ntypeaux(nshellaux),ngfaux(nshellaux)                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!=====================================================================!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  SHELL  !!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        write(*,*)'*** Computing aux-NAI intgral ***'
        write(99,'(30X,A30)')'**********************'
        write(99,'(30X,A30)')'AUX NAI INTEGRALS'
        write(99,'(30X,A30)')'**********************'          
!!!!!!!         
!============== <a|1/rc|b> ==============! vab
!!!! Note that a is from aux basis set,diff from b.
        num=0
 !!!!!!       
        do Ish=1,nshellaux   ! I shell
            Iatom=kcenteraux(Ish)
            cix=coor(Iatom,1)
            ciy=coor(Iatom,2)
            ciz=coor(Iatom,3)
            do Jsh=1,nshell  ! J shell  
                Jatom=kcenter(Jsh)
                cJx=coor(Jatom,1)
                cJy=coor(Jatom,2)
                cJz=coor(Jatom,3)
                Rij=(cix-cjx)**2+(ciy-cjy)**2+(ciz-cjz)**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! BASIS FUNCTION  !!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                do ibas=kbasstaux(ish),kbasedaux(ish)
                    nax=naxyzaux(ibas)
                    nay=naxyzaux(ibas+nbasisaux)
                    naz=naxyzaux(ibas+2*nbasisaux)
                    ipgfst=kpgfstaux(ish)
                    ipgfed=kpgfedaux(ish)
!!!!! only for sp shell                    
                    if(ntypeaux(ish)==12)then
                        if(nax+nay+naz==0) then
                            ipgfed=kpgfedaux(ish)-ngfaux(ish)/2
                        else  
                            ipgfst=kpgfstaux(ish)+ngfaux(ish)/2
                        end if
                    end if  
!!!!!!!!!                                          
                    do jbas=kbasst(jsh),kbased(jsh)
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
!!! primitive gauss functions  <a|
                        vij=0d0
                        vijn=0d0                     
                        do ipgf=ipgfst,ipgfed     
                                ea=gfeaux(ipgf)
                                ca=gfcaux(ipgf)
                                dNa=dnormgf(ea,nax,nay,naz)                                
!!! primitive gauss functions |b>
                                do jpgf=jpgfst,jpgfed  
                                    eb=gfe(jpgf)
                                    cb=gfc(jpgf)
                                    dNb=dnormgf(eb,nbx,nby,nbz)                               
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
                                    naxx=nax
                                    nbxx=nbx
                                    Axx=cix
                                    Bxx=cjx                    
                                    nayy=nay
                                    nbyy=nby
                                    Ayy=ciy
                                    Byy=cjy
                                    nazz=naz
                                    nbzz=nbz
                                    Azz=ciz
                                    Bzz=cjz  
                                    vab=0d0 
                                do  kc=1,natom
                                    CX=coor(kc,1)
                                    CY=coor(kc,2)
                                    CZ=coor(kc,3)
                                    NZ=katom(kc)
                                    TT=ep*((PX-CX)**2+(PY-CY)**2+(PZ-CZ)**2)
                                    if(abs(TT)<1d-12) TT=0d0
!!!!
                        if (methodV==1)then
!!!! Method 1 
!!! Ohata method          
                        vab=vab+anaNAI1(naxx,nayy,nazz,nbxx,nbyy,nbzz,&
                                    ep1,ABK,Axx,Bxx,Px,Cx,Ayy,Byy,Py,Cy,&
                                    Azz,Bzz,Pz,Cz,TT)*nz*2*pi*ep1*ABK
                        else if (methodV==2)then
!!!! Method 2     
!!! DRK  method            
                            vab=vab+nz*recNAI1fast(naxx,nayy,nazz,nbxx,nbyy,nbzz,&
                                    ep1,ABK,Axx,Bxx,Px,Cx,Ayy,Byy,Py,Cy,&
                                    Azz,Bzz,Pz,Cz,TT)                                                                  
                        else if (methodV==3)then                                    
!!!!! Method 3     
!!!  Obara,Saika method            
                        vab=vab+nz*recNAI2fast(naxx,nayy,nazz,nbxx,nbyy,nbzz,&
                                    ep1,ABK,Axx,Bxx,Px,Cx,Ayy,Byy,Py,Cy,&
                                    Azz,Bzz,Pz,Cz,TT)                                                                       
                        else 
                            stop "bad NAI method"                                    
                        end if                         
                    end do
                    vabn=vab*dNa*dNb
!================= <a|1/rc|b> -> sum(ca*cb*<a|1/rc|b>)==============!
!!  V(i,j)                   
                    vijn=vijn+vabn*ca*cb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            end do                          ! b pgf
                        end do                              ! a pgf
!!!!!!!!!!!!  Matrix assignment                    
                        ! if(abs(vijn)<1d-12) vijn=0d0   
                        V(ibas,jbas)= vijn                  
!!!!!!!!!!!
                    end do  ! J basis             
                end do    ! I basis  
!!!!!!!!!!!
            end do      ! J shell
        end do         ! I shell
        write(*,*)'*** Aux-NAI intgral Done***'        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end  subroutine
!======================================================================!
    end  Module Int1e
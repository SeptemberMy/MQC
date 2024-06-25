!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                         CoupledCluster.f90                         !!  
!!                  Coupled cluster Theory:CCSD,CCSD(T)               !!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    subroutine CCSDT(kfid,nbasis,neri,Nele,Hcore,ERISPIN,Eelec,Eorb,mcc,&
                   ccdamp,maxcyclecc,conver_Ecc,conver_rmscc,mdiiscc,ndiiscc,&
                   Emp2,Ecc)
        use Global
        implicit double precision (A-H,O-Z)
        Integer*16 neri
        dimension::ERISPIN(2*nbasis,2*nbasis,2*nbasis,2*nbasis)
        dimension::Eorb(nbasis),fockspin(2*nbasis,2*nbasis),hcore(nbasis,nbasis)
        allocatable::tia(:,:),tijab(:,:,:,:),Dia(:,:),Dijab(:,:,:,:),tianew(:,:),tijabnew(:,:,:,:)
        allocatable::Fvv(:,:),Foo(:,:),Fov(:,:)
        allocatable::tauijab1(:,:,:,:),tauijab2(:,:,:,:)
        allocatable::Woooo(:,:,:,:),Wvvvv(:,:,:,:),Wovvo(:,:,:,:)
!!!!! for diis           
        allocatable::errmat1(:,:,:),T1diis(:,:,:),errmat2(:,:,:,:,:),T2diis(:,:,:,:,:)   
        dimension:: Bdiis1(ndiiscc+1,ndiiscc+1),Cofdiis1(ndiiscc+1),&
                    Bdiis2(ndiiscc+1,ndiiscc+1),Cofdiis2(ndiiscc+1),ipivdiis(ndiiscc+1)                 
!!!!!!!!!!! J. Chem. Phys. 94, 4334–4345 (1991)        
!!!!
        write(*,*)'*** Start CCSD procedure ***' 
        write(kfid,*)'=======================================&
                     =========================================&
                     ===================='  
        if(mdiiscc==1)write(kfid,*)'DIIS: on','       ','Number of diis matrix: ',ndiiscc
        write(kfid,'(1X,A5,1x,f4.2)')'damp:',ccdamp
            write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------'                           
        write(kfid,'(35X,40A)')'*******************'
        write(kfid,'(35X,40A)')'* CCSD Iterations *'   
        write(kfid,'(35X,40A)')'*******************' 
        nspin=2*nbasis
        nspinocc=Nele
        nspinvir=nspin-nspinocc
        Ethreshold=10d0**(-conver_Ecc)
        rmsthreshold=10d0**(-conver_rmscc)
!!!!  init
        allocate(tia(nspinocc,nspinvir),tijab(nspinocc,nspinocc,nspinvir,nspinvir),&
                tianew(nspinocc,nspinvir),tijabnew(nspinocc,nspinocc,nspinvir,nspinvir),&
                 Dia(nspinocc,nspinvir),Dijab(nspinocc,nspinocc,nspinvir,nspinvir))
        allocate(Fvv(nspinvir,nspinvir),Foo(nspinocc,nspinocc),Fov(nspinocc,nspinvir))
        allocate(Woooo(nspinocc,nspinocc,nspinocc,nspinocc),&
                 Wvvvv(nspinvir,nspinvir,nspinvir,nspinvir),&
                 Wovvo(nspinocc,nspinvir,nspinvir,nspinocc))
        allocate(tauijab1(nspinocc,nspinocc,nspinvir,nspinvir),&
                 tauijab2(nspinocc,nspinocc,nspinvir,nspinvir))
!!!! Set the spin-orbital Fock matrix  
        ! fockspin=0d0    
        ! do ii=1,nspin
        !     iispa=ceiling(ii/2d0)
        !     iimod=mod(ii,2)
        !     do ij=1,nspin
        !         ijspa=ceiling(ij/2d0)
        !         ijmod=mod(ij,2)
        !         fockspin(ii,ij)=Hcore(iispa,ijspa)*krodelta(iimod,ijmod)
        !         do im=1,nspinocc
        !             fockspin(ii,ij)=fockspin(ii,ij)+ERISPIN(ii,ij,im,im)
        !         end do 
        !     end do
        ! end do    
!!!   for hartree-fock reference       
        fockspin=0d0    
        do ii=1,nspin
            iispa=ceiling(ii/2d0)
            fockspin(ii,ii)=Eorb(iispa)
        end do  
        tia=0d0
        tijab=0d0       
!!!! set tia and tijab. Compute MP2.        
        Emp2tem=0d0
        do ii=1,nspinocc
            iispa=ceiling(ii/2d0)
            do ia=nspinocc+1,nspin
                iaspa=ceiling(ia/2d0)   
                Dia(ii,ia-nspinocc)= Eorb(iispa)-Eorb(iaspa)       
                do ij=1,nspinocc
                ijspa=ceiling(ij/2d0)
                    do ib=nspinocc+1,nspin
                    ibspa=ceiling(ib/2d0)
                    Dijab(ii,ij,ia-nspinocc,ib-nspinocc)=Eorb(iispa)+Eorb(ijspa)-Eorb(iaspa)-Eorb(ibspa)
                    tijabnow=ERISPIN(ii,ia,ij,ib)/Dijab(ii,ij,ia-nspinocc,ib-nspinocc)
                    Emp2tem=Emp2tem+ERISPIN(ii,ia,ij,ib)*tijabnow
                    tijab(ii,ij,ia-nspinocc,ib-nspinocc)=tijabnow
                    end do 
                end do 
            end do 
        end do
        Emp2=0.25d0*Emp2tem
!!!!! CCSD
!!! initial energy
        ECC=0d0
        do ii=1,nspinocc
            do ia=1,nspinvir
                do ij=1,nspinocc
                    do ib=1,nspinvir
                    tauijab1(ii,ij,ia,ib)=tijab(ii,ij,ia,ib)+0.5d0*(tia(ii,ia)*tia(ij,ib)-tia(ii,ib)*tia(ij,ia))
                    tauijab2(ii,ij,ia,ib)=tijab(ii,ij,ia,ib)+(tia(ii,ia)*tia(ij,ib)-tia(ii,ib)*tia(ij,ia))                    
                    ECC=ECC+0.25d0*tijab(ii,ij,ia,ib)*ERISPIN(ii,ia+nspinocc,ij,ib+nspinocc)
                    end do
                end do
            end do
        end do
        write(kfid,*) 'CCSD initial energy:',ECC,'a.u.'
!!!!
        write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------' 
        write(kfid,'(2X,A5,3X,6(A10,13X))')'Cycle','Ecc','Ecc(old)','delt E','rms'
        write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------'     
        if(mdiiscc==1)then
            allocate(errmat1(nspinocc,nspinvir,ndiiscc),T1diis(nspinocc,nspinvir,ndiiscc),&
                     errmat2(nspinocc,nspinocc,nspinvir,nspinvir,ndiiscc),&
                     T2diis(nspinocc,nspinocc,nspinvir,nspinvir,ndiiscc) )
        end if      
!!!!!!!!                     
        do kcy=1,maxcyclecc
            write(*,*)'CC iteration num :',kcy
!!!!            
!!!   set Fvv      
            Fvv=0d0
            do ia=1,nspinvir
                do ie=1,nspinvir
                    Fvv(ia,ie)=(1d0-krodelta(ia,ie))*fockspin(ia+nspinocc,ie+nspinocc) !! 1
                    tem=0d0
                    do im=1,nspinocc
                        tem=tem+tia(im,ia)*fockspin(im,ie+nspinocc)
                    end do
                    Fvv(ia,ie)=Fvv(ia,ie)-0.5d0*tem  !! 2
                    tem=0d0
                    do im=1,nspinocc
                        do iff=1,nspinvir
                            tem=tem+tia(im,iff)*ERISPIN(im,iff+nspinocc,ia+nspinocc,ie+nspinocc)
                        end do
                    end do
                    Fvv(ia,ie)=Fvv(ia,ie)+tem !! 3
                    tem=0d0
                    do im=1,nspinocc
                        do inn=1,nspinocc
                            do iff=1,nspinvir
                                tem=tem+tauijab1(im,inn,ia,iff)*ERISPIN(im,ie+nspinocc,inn,iff+nspinocc)
                            end do
                        end do
                    end do
                    Fvv(ia,ie)=Fvv(ia,ie)-0.5d0*tem  !! 4
                end do !!! ie
            end do !!! ia   
!!   set Foo
            Foo=0d0
            do im=1,nspinocc
                do ii=1,nspinocc
                    Foo(im,ii)=(1d0-krodelta(im,ii))*fockspin(im,ii) !! 1
                    tem=0d0
                    do ie=1,nspinvir
                        tem=tem+tia(ii,ie)*fockspin(im,ie+nspinocc)
                    end do
                    Foo(im,ii)=Foo(im,ii)+0.5d0*tem  !! 2
                    tem=0d0
                    do inn=1,nspinocc
                        do ie=1,nspinvir
                            tem=tem+tia(inn,ie)*ERISPIN(im,ii,inn,ie+nspinocc)
                        end do
                    end do
                    Foo(im,ii)=Foo(im,ii)+tem !! 3
                    tem=0d0
                    do ie=1,nspinvir
                        do inn=1,nspinocc
                            do iff=1,nspinvir
                                tem=tem+tauijab1(ii,inn,ie,iff)*ERISPIN(im,ie+nspinocc,inn,iff+nspinocc)
                            end do
                        end do
                    end do
                    Foo(im,ii)=Foo(im,ii)+0.5d0*tem  !! 4
                end do !!! ie
            end do !!! ia   
!!  set Fov
            Fov=0d0
            do im=1,nspinocc
                do ie=1,nspinvir
                    Fov(im,ie)=fockspin(im,ie+nspinocc)
                    tem=0d0
                    do inn=1,nspinocc
                        do iff=1,nspinvir
                            tem=tem+tia(inn,iff)*ERISPIN(im,ie+nspinocc,inn,iff+nspinocc)
                        end do
                    end do  
                    Fov(im,ie)=Fov(im,ie)+tem                  
                end do  !! ie
            end do   !! im
!!  set Woooo
            Woooo=0d0
            do im=1,nspinocc
              do inn=1,nspinocc
                do ii=1,nspinocc
                  do ij=1,nspinocc
                    Woooo(im,inn,ii,ij)=ERISPIN(im,ii,inn,ij) !! 1
                    tem=0d0
                    do ie=1,nspinvir
                        tem=tem+(tia(ij,ie)*ERISPIN(im,ii,inn,ie+nspinocc)-tia(ii,ie)*ERISPIN(im,ij,inn,ie+nspinocc))
                    end do
                    Woooo(im,inn,ii,ij)=Woooo(im,inn,ii,ij)+tem !! 2 
                    tem=0d0
                    do ie=1,nspinvir
                        do iff=1,nspinvir
                        tem=tem+tauijab2(ii,ij,ie,iff)*ERISPIN(im,ie+nspinocc,inn,iff+nspinocc)
                        end do
                    end do
                    Woooo(im,inn,ii,ij)=Woooo(im,inn,ii,ij)+0.25d0*tem !! 3
                  end do !! ij
                end do  !! ii
              end do   !! in
            end do   !! im
!!  set Wvvvv
            Wvvvv=0d0
            do ia=1,nspinvir
              do ib=1,nspinvir
                do ie=1,nspinvir
                  do iff=1,nspinvir   
                    Wvvvv(ia,ib,ie,iff)=ERISPIN(ia+nspinocc,ie+nspinocc,ib+nspinocc,iff+nspinocc) !! 1
                    tem=0d0
                    do im=1,nspinocc
                        tem=tem+(tia(im,ib)*ERISPIN(ia+nspinocc,ie+nspinocc,im,iff+nspinocc)&
                               -tia(im,ia)*ERISPIN(ib+nspinocc,ie+nspinocc,im,iff+nspinocc))
                    end do
                    Wvvvv(ia,ib,ie,iff)=Wvvvv(ia,ib,ie,iff)-tem !! 2
                    tem=0d0
                    do im=1,nspinocc
                        do inn=1,nspinocc
                        tem=tem+tauijab2(im,inn,ia,ib)*ERISPIN(im,ie+nspinocc,inn,iff+nspinocc)
                        end do
                    end do         
                    Wvvvv(ia,ib,ie,iff)=Wvvvv(ia,ib,ie,iff)+0.25d0*tem !! 3                               
                  end do !! iff
                end do  !! ie
              end do   !! ib 
            end do   !! ia                        
!!  set Wovvo
            Wovvo=0d0
            do im=1,nspinocc
              do ib=1,nspinvir
                do ie=1,nspinvir
                  do ij=1,nspinocc 
                    Wovvo(im,ib,ie,ij)=ERISPIN(im,ie+nspinocc,ib+nspinocc,ij) !! 1
                    tem=0d0
                    do iff=1,nspinvir
                        tem=tem+tia(ij,iff)*ERISPIN(im,ie+nspinocc,ib+nspinocc,iff+nspinocc)
                    end do    
                    Wovvo(im,ib,ie,ij)=Wovvo(im,ib,ie,ij)+tem !! 2
                    tem=0d0
                    do inn=1,nspinocc
                        tem=tem+tia(inn,ib)*ERISPIN(im,ie+nspinocc,inn,ij)
                    end do
                    Wovvo(im,ib,ie,ij)=Wovvo(im,ib,ie,ij)-tem !! 3
                    tem=0d0
                    do iff=1,nspinvir
                        do inn=1,nspinocc
                tem=tem+(0.5d0*tijab(ij,inn,iff,ib)+tia(ij,iff)*tia(inn,ib))*ERISPIN(im,ie+nspinocc,inn,iff+nspinocc)
                        end do
                    end do        
                    Wovvo(im,ib,ie,ij)=Wovvo(im,ib,ie,ij)-tem !! 4
                  end do !! ij
                end do  !! ie
              end do   !! ib
            end do   !! im                           
!!!!!   Build the new tia, tijab
        do ii=1,nspinocc
            do ia=1,nspinvir
                tianow=fockspin(ii,ia+nspinocc) !! 1
                tem=0d0
                do ie=1,nspinvir
                    tem=tem+tia(ii,ie)*Fvv(ia,ie)
                end do
                tianow=tianow+tem !! 2
                tem=0d0
                do im=1,nspinocc
                    tem=tem+tia(im,ia)*Foo(im,ii)
                end do     
                tianow=tianow-tem !! 3
                tem=0d0
                do im=1,nspinocc
                    do ie=1,nspinvir
                        tem=tem+tijab(ii,im,ia,ie)*Fov(im,ie)
                    end do
                end do     
                tianow=tianow+tem !! 4     
                tem=0d0
                do inn=1,nspinocc
                    do iff=1,nspinvir
                        tem=tem+tia(inn,iff)*ERISPIN(inn,ii,ia+nspinocc,iff+nspinocc)
                    end do
                end do     
                tianow=tianow-tem !! 5     
                tem=0d0
                do im=1,nspinocc
                  do ie=1,nspinvir
                    do iff=1,nspinvir
                        tem=tem+tijab(ii,im,ie,iff)*ERISPIN(im,ie+nspinocc,ia+nspinocc,iff+nspinocc)
                    end do
                  end do
                end do     
                tianow=tianow-0.5d0*tem !! 6          
                tem=0d0
                do im=1,nspinocc
                  do ie=1,nspinvir
                    do inn=1,nspinocc
                        tem=tem+tijab(im,inn,ia,ie)*ERISPIN(inn,ie+nspinocc,im,ii)
                    end do
                  end do
                end do     
                tianow=tianow-0.5d0*tem !! 6    
                tianew(ii,ia)=tianow/Dia(ii,ia)
!! T2                
                do ij=1,nspinocc
                    do ib=1,nspinvir
                    tijabnow=ERISPIN(ii,ia+nspinocc,ij,ib+nspinocc) !! 1
                    tem=0d0
                    do ie=1,nspinvir
                        tem1=0d0
                        tem2=0d0                        
                        do im=1,nspinocc
                            tem1=tem1+tia(im,ib)*Fov(im,ie)
                            tem2=tem2+tia(im,ia)*Fov(im,ie)                            
                        end do
                        tem1=tem1*0.5d0  
                        tem2=tem2*0.5d0                     
                        tem=tem+(tijab(ii,ij,ia,ie)*(Fvv(ib,ie)-tem1)-tijab(ii,ij,ib,ie)*(Fvv(ia,ie)-tem2))
                    end do
                    tijabnow=tijabnow+tem !! 2
                    tem=0d0
                    do im=1,nspinocc
                        tem1=0d0
                        tem2=0d0                        
                        do ie=1,nspinvir
                            tem1=tem1+tia(ij,ie)*Fov(im,ie)
                            tem2=tem2+tia(ii,ie)*Fov(im,ie)                            
                        end do
                        tem1=tem1*0.5d0
                        tem2=tem2*0.5d0                     
                        tem=tem+(tijab(ii,im,ia,ib)*(Foo(im,ij)+tem1)-tijab(ij,im,ia,ib)*(Foo(im,ii)+tem2))
                    end do
                    tijabnow=tijabnow-tem !! 3     
                    tem=0d0   
                    do im=1,nspinocc
                        do inn=1,nspinocc 
                            tem=tem+tauijab2(im,inn,ia,ib)*Woooo(im,inn,ii,ij)
                        end do    
                    end do    
                    tijabnow=tijabnow+0.5d0*tem !! 4     
                    tem=0d0   
                    do ie=1,nspinvir
                        do iff=1,nspinvir 
                            tem=tem+tauijab2(ii,ij,ie,iff)*Wvvvv(ia,ib,ie,iff)
                        end do    
                    end do    
                    tijabnow=tijabnow+0.5d0*tem !! 5 
                    tem=0d0
                    do ie=1,nspinvir
                        do im=1,nspinocc 
            tem=tem+(tijab(ii,im,ia,ie)*Wovvo(im,ib,ie,ij)-tia(ii,ie)*tia(im,ia)*ERISPIN(im,ie+nspinocc,ib+nspinocc,ij))&
                    -(tijab(ii,im,ib,ie)*Wovvo(im,ia,ie,ij)-tia(ii,ie)*tia(im,ib)*ERISPIN(im,ie+nspinocc,ia+nspinocc,ij))&
                    -(tijab(ij,im,ia,ie)*Wovvo(im,ib,ie,ii)-tia(ij,ie)*tia(im,ia)*ERISPIN(im,ie+nspinocc,ib+nspinocc,ii))&
                    +(tijab(ij,im,ib,ie)*Wovvo(im,ia,ie,ii)-tia(ij,ie)*tia(im,ib)*ERISPIN(im,ie+nspinocc,ia+nspinocc,ii))
                        end do    
                    end do    
                    tijabnow=tijabnow+tem !! 6    
                    tem=0d0   
                    do ie=1,nspinvir
                        tem=tem+(tia(ii,ie)*ERISPIN(ia+nspinocc,ie+nspinocc,ib+nspinocc,ij)&
                                -tia(ij,ie)*ERISPIN(ia+nspinocc,ie+nspinocc,ib+nspinocc,ii))
                    end do    
                    tijabnow=tijabnow+tem !! 7  
                    tem=0d0   
                    do im=1,nspinocc
                        tem=tem+(tia(im,ia)*ERISPIN(im,ii,ib+nspinocc,ij)-tia(im,ib)*ERISPIN(im,ii,ia+nspinocc,ij))
                    end do    
                    tijabnow=tijabnow-tem !! 8                                                              
!!!                    
                    tijabnew(ii,ij,ia,ib)=tijabnow/Dijab(ii,ij,ia,ib)
                    end do 
                end do 
            end do !! ia
        end do  !! ii
!!!! diis
        if(mdiiscc==1 )then         
            morderdiis=mod(kcy,ndiiscc)
!!! store T1,T2 and set error matrix            
            T1diis(:,:,morderdiis+1)=tianew
            T2diis(:,:,:,:,morderdiis+1)=tijabnew    
            errmat1(:,:,morderdiis+1)=tianew-tia 
            errmat2(:,:,:,:,morderdiis+1)=tijabnew-tijab    
!!!   set B matrix         
            if(kcy>ndiiscc)then
                Cofdiis1=0d0
                Cofdiis2=0d0               
                Bdiis1=-1d0
                Bdiis2=-1d0
                do  kp=1,ndiiscc
                    do kq=1,kp
                        sumdiis1=0d0
                        sumdiis2=0d0
                        do ii=1,nspinocc
                            do jj=1,nspinvir
                                sumdiis1=sumdiis1+errmat1(ii,jj,kp)*errmat1(jj,ii,kq)
                                 do kk=1,nspinocc
                                  do ll=1,nspinvir
                                   sumdiis2=sumdiis2+errmat2(ii,kk,jj,ll,kp)*errmat2(kk,ii,ll,jj,kq)
                                  end do !! ll
                                 end do !! kk
                                end do  !! jj
                            end do      !! ii
                            Bdiis1(kp,kq)=sumdiis1
                            Bdiis2(kp,kq)=sumdiis2
                            if(kp/=kq)then
                                Bdiis1(kq,kp)=sumdiis1
                                Bdiis2(kq,kp)=sumdiis2
                            end if
                        end do    !! kq
                    end do        !! kp            
                    Bdiis1(ndiiscc+1,ndiiscc+1)=0d0
                    Bdiis2(ndiiscc+1,ndiiscc+1)=0d0
                    Cofdiis1(ndiiscc+1)=-1d0
                    Cofdiis2(ndiiscc+1)=-1d0
!!!    solve C  
                    call dgesv(ndiiscc+1,1,Bdiis1,ndiiscc+1,ipivdiis,Cofdiis1,ndiiscc+1,infodiis1)
                    call dgesv(ndiiscc+1,1,Bdiis2,ndiiscc+1,ipivdiis,Cofdiis2,ndiiscc+1,infodiis2)
                    if(infodiis1/=0 .or. infodiis2/=0)stop 'Solve diis cof failed!!!'        
!!! set new T1,T2      
                    ! write(*,*)Cofdiis1
                    ! write(*,*)Cofdiis2
                    tianew=0d0
                    tijabnew=0d0
                    do iii=1,ndiiscc
                        tianew=tianew+Cofdiis1(iii)*T1diis(:,:,iii)  
                        tijabnew=tijabnew+Cofdiis2(iii)*T2diis(:,:,:,:,iii)       
                    end do    
            end if    
        end if     !! diis     
!!!! damp
        tianew=ccdamp*tia+(1d0-ccdamp)*tianew
        tijabnew=ccdamp*tijab+(1d0-ccdamp)*tijabnew
!!!! new energy
        ECC=0d0
        do ii=1,nspinocc
            do ia=1,nspinvir
                ECC=ECC+fockspin(ii,ia+nspinocc)*tianew(ii,ia)
                do ij=1,nspinocc
                    do ib=1,nspinvir               
                    ECC=ECC+0.25d0*tijabnew(ii,ij,ia,ib)*ERISPIN(ii,ia+nspinocc,ij,ib+nspinocc)&
                           +0.5d0*tianew(ii,ia)*tianew(ij,ib)*ERISPIN(ii,ia+nspinocc,ij,ib+nspinocc)  
                    end do
                end do
            end do
        end do          
!!! rms
            rms=sum(sqrt((tia-tianew)**2))                
!!!! output
            write(kfid,'(I4,5(F23.16))')kcy,ECC,ECCold,ECC-ECCold,rms      
!!!!  Check for Convergence  
            if(abs(Ecc-ECCold)<=Ethreshold .and. abs(rms)<=rmsthreshold )then
                write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------'      
                write(kfid,*) 'CC convergence' 
                write(kfid,'(2x,"After",I4,2x,"cycles")') kcy         
                write(kfid,*)' Final CC Energy :',ECC,'A.U.'
                write(kfid,*)' The CC convergence threshold of E :',Ethreshold,'A.U.'
                write(kfid,*)' The SCF convergence threshold of RMS :',rmsthreshold,'A.U.'
                exit
            else if(kcy==maxcyclecc)then
                write(kfid,*) 'Waring: Bad convergence of CC!'  
                write(*,*) 'Waring: Bad convergence of CC!'                
            end if      
!!!! update
            ECCold=ECC
            tia=tianew
            tijab=tijabnew
            do ii=1,nspinocc
                do ia=1,nspinvir
                    do ij=1,nspinocc
                        do ib=1,nspinvir
                tauijab1(ii,ij,ia,ib)=tijab(ii,ij,ia,ib)+0.5d0*(tia(ii,ia)*tia(ij,ib)-tia(ii,ib)*tia(ij,ia))
                tauijab2(ii,ij,ia,ib)=tijab(ii,ij,ia,ib)+(tia(ii,ia)*tia(ij,ib)-tia(ii,ib)*tia(ij,ia))                    
                        end do
                    end do
                end do
            end do                                                       
        end do   !!! kcy    
!!!!        
        write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------'    
        write(kfid,'(25X,A30)')'*****  CCSD Summary *****'
        write(kfid,*) 'CCSD Final results (a.u.)'
        write(kfid,'(25X,A11,5X,F20.15)')'E(HF_tol):',Eelec
        write(kfid,'(25X,A11,5X,F20.15)')'E(MP2):',Emp2        
        write(kfid,'(25X,A11,5X,F20.15)')'E(HF+MP2):',Emp2+Eelec
        write(kfid,'(25X,A11,5X,F20.15)')'E(CCSD):',ECC 
        write(kfid,'(25X,A11,5X,F20.15)')'E(HF+ECCSD):',Ecc+Eelec
!!!!!!  CCSD(T)    
        if(mcc==2)then
            write(*,*)'*** Start CCSD(T) procedure ***' 
            Eccsdt=0d0
            icount=0
            idone=0
            do ii=1,nspinocc
                nall10=nspinocc/10
                icount=icount+1
                if(icount>=nall10)then
                    idone=idone+icount
                    write(*,*)'(T)',100*real(idone)/nspinocc,'% done'
                    icount=0
                end if
               iispa=ceiling(ii/2d0)
             do ij=1,nspinocc
               ijspa=ceiling(ij/2d0)
              do ik=1,nspinocc
                ikspa=ceiling(ik/2d0)
               do ia=1,nspinvir
                iaspa=ceiling((ia+nspinocc)/2d0)
                do ib=1,nspinvir
                 ibspa=ceiling((ib+nspinocc)/2d0)
                 do ic=1,nspinvir
                  icspa=ceiling((ic+nspinocc)/2d0)     
               
                    dijkabc=Eorb(iispa)+Eorb(ijspa)+Eorb(ikspa)&
                            -Eorb(iaspa)-Eorb(ibspa)-Eorb(icspa)     
                        !     write(*,*)dijkabc
                        !   stop ' here'     
!!!   the "disconnected" triples 
                    DT31=ERISPIN(ij,ib+nspinocc,ik,ic+nspinocc)*tianew(ii,ia)&
                        -tianew(ii,ib)*ERISPIN(ij,ia+nspinocc,ik,ic+nspinocc)&
                        -tianew(ii,ic)*ERISPIN(ij,ib+nspinocc,ik,ia+nspinocc)&
                        -(tianew(ij,ia)*ERISPIN(ii,ib+nspinocc,ik,ic+nspinocc)&
                        -tianew(ij,ib)*ERISPIN(ii,ia+nspinocc,ik,ic+nspinocc)&
                        -tianew(ij,ic)*ERISPIN(ii,ib+nspinocc,ik,ia+nspinocc))&
                        -(tianew(ik,ia)*ERISPIN(ij,ib+nspinocc,ii,ic+nspinocc)&
                        -tianew(ik,ib)*ERISPIN(ij,ia+nspinocc,ii,ic+nspinocc)&
                        -tianew(ik,ic)*ERISPIN(ij,ib+nspinocc,ii,ia+nspinocc))
!!!   the "connected" triples 
                    temp11=0d0
                    temp12=0d0
                    do ie=1,nspinvir
                        temp11=temp11+tijabnew(ij,ik,ia,ie)*ERISPIN(ie+nspinocc,ib+nspinocc,ii,ic+nspinocc)&
                                     -tijabnew(ij,ik,ib,ie)*ERISPIN(ie+nspinocc,ia+nspinocc,ii,ic+nspinocc)&
                                     -tijabnew(ij,ik,ic,ie)*ERISPIN(ie+nspinocc,ib+nspinocc,ii,ia+nspinocc)&
                                     -(tijabnew(ii,ik,ia,ie)*ERISPIN(ie+nspinocc,ib+nspinocc,ij,ic+nspinocc)&
                                     -tijabnew(ii,ik,ib,ie)*ERISPIN(ie+nspinocc,ia+nspinocc,ij,ic+nspinocc)&
                                     -tijabnew(ii,ik,ic,ie)*ERISPIN(ie+nspinocc,ib+nspinocc,ij,ia+nspinocc))&
                                     -(tijabnew(ij,ii,ia,ie)*ERISPIN(ie+nspinocc,ib+nspinocc,ik,ic+nspinocc)&
                                     -tijabnew(ij,ii,ib,ie)*ERISPIN(ie+nspinocc,ia+nspinocc,ik,ic+nspinocc)&
                                     -tijabnew(ij,ii,ic,ie)*ERISPIN(ie+nspinocc,ib+nspinocc,ik,ia+nspinocc))
                    end do 
                    do im=1,nspinocc
                        temp12=temp12+tijabnew(ii,im,ib,ic)*ERISPIN(im,ij,ia+nspinocc,ik)&
                                     -tijabnew(ii,im,ia,ic)*ERISPIN(im,ij,ib+nspinocc,ik)&
                                     -tijabnew(ii,im,ib,ia)*ERISPIN(im,ij,ic+nspinocc,ik)&
                                     -(tijabnew(ij,im,ib,ic)*ERISPIN(im,ii,ia+nspinocc,ik)&
                                     -tijabnew(ij,im,ia,ic)*ERISPIN(im,ii,ib+nspinocc,ik)&
                                     -tijabnew(ij,im,ib,ia)*ERISPIN(im,ii,ic+nspinocc,ik))&
                                     -(tijabnew(ik,im,ib,ic)*ERISPIN(im,ij,ia+nspinocc,ii)&
                                     -tijabnew(ik,im,ia,ic)*ERISPIN(im,ij,ib+nspinocc,ii)&
                                     -tijabnew(ik,im,ib,ia)*ERISPIN(im,ij,ic+nspinocc,ii))
                    end do
                    DT32=temp11-temp12
                    T32=DT32/dijkabc   
                    Eccsdt=Eccsdt+T32*(DT32+DT31)
                 end do !! ic
                end do !! ib
               end do !! ia
            end do !! ik
            end do !! ij
            end do !! ii
           Eccsdt=Eccsdt/36d0
!!!!        
            write(kfid,*)'------------------------------------& 
                    -----------------------------------------&
                    ---------------------'    
            write(kfid,'(25X,A30)')'*****  CCSD(T) Summary *****'
            write(kfid,*) 'CCSD Final results (a.u.)'
            write(kfid,'(25X,A11,5X,F20.15)')'E(HF_tol):',Eelec
            write(kfid,'(25X,A11,5X,F20.15)')'E(CCSD):',ECC 
            write(kfid,'(25X,A11,5X,F20.15)')'E(HF+ECCSD):',Ecc+Eelec
            write(kfid,'(25X,A11,5X,F20.15)')'E(CCSD(T)):',Eccsdt 
            write(kfid,'(25X,A11,5X,F20.15)')'E(HF+ECC):',Ecc+Eelec+Eccsdt    
            Ecc=Ecc+Eccsdt
        end if    
!!!!
        deallocate(tia,tijab,tianew,tijabnew,Dia,Dijab,Fvv,Foo,Fov,Woooo,Wvvvv,&
                   Wovvo,tauijab1,tauijab2)
        if(mdiiscc==1)then
            deallocate(errmat1,T1diis,errmat2,T2diis )
        end if                     
!!!!        
    end subroutine CCSDT
!!!=================================================================================!!!    









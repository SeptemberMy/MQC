    include 'Transint.f90'
    include 'Perturbation.f90'
    include 'CoupledCluster.f90'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    program main  
        use Global
        use Init
        use IO
        use Int1e
        use Int2e
        use HF
        implicit double precision (A-H,O-Z)
!!!!============================Init================================!!!!   
        write(*,*)'*** Initialization ***'
        inputfpath='INPUT.inp'  
        outputfpath="OUTPUT.log"
        call getarg(1,inputarg)
        call getarg(2,outputarg)
        if(len_trim(inputarg)/=0)then
            inputfpath=inputarg
            write(*,*)inputarg
        end if
        if(len_trim(outputarg)/=0)then
            outputfpath=outputarg
            write(*,*)outputarg
        end if        
!!!! Set default value
        call Setvalue           
!!!! Write OUTPUT (fid=99)
        open (kfid, file =outputfpath,status='replace')
        call logo(kfid)   
!!!! Read INPUT (fid=33) 
        write(kfid,*)'INPUT file : ',inputfpath  
        write(kfid,*)'OUTPUT file : ',outputfpath   
        call readinput  
!!!! Read basis set  
        write(kfid,*)'basis name: ',basisname
        write(kfid,*)'basis file path: ',basfpath  
        call readbas
!!!! Set shell information
        call Setbasicinf                         
!!!!  Output         
        call Structure(kfid,Natom,Nele,coor,katom,Moutxyzfile)
        call basisout(kfid,0,natom,nshell,KATOM,kcenter,gfc,gfe,ntype,ngf,ngfall,&
                knbas,kbasst,kbased,kpgfst,kpgfed,norb)              
!!!! time
        call date_and_time(theDate1,t1_str) 
        call cpu_time(t_start)
        call system_Clock(it1,irate) 
!!!!        
        write(*,*)'*** Initialization done ***'      
!!!!!!!!!!!!!!!==================== Main =====================!!!!!!!!!!!!!!
!!!!!!  main      
        write(*,*)'*** Computing ***'
        write(kfid,*)'=======================================&
                  =========================================&
                  ==============================='
        write(kfid,'(25X,A30)')'``` Result```'
!!! Nuclear repulsion energy
        call CalEnucre 
        if(usecint)then
!!
          tst=counttime_cputime(t_start) !!! count time 1
          tst2=counttime_sysckock(it1)   !!! count time 2   
!!                     
          call  CINT_ERI
!!          
          if(Moutlevel>2)then
                lprint1e=1
          end if
          if(Moutlevel<4)then
                lprint2e=0
          else if(Moutlevel==4)then
                lprint2e=1
          else if(Moutlevel==5)then
                lprint2e=2
          end if                
          call Oneintout(kfid,nbasis,methods,methodT,methodV,lprint1e,&
                        S,T,V,Hcore,dipx,dipy,dipz)                
          call TWOintout(kfid,nbasis,neri,1,lprint2e,ERIall)
!! 
          ted=counttime_cputime(t_start) !!! count time 1
          ted2=counttime_sysckock(it1)   !!! count time 2
          write(kfid,*) '2e INT time(libcint) :', ted-tst ,'s'     !!! count time 1
          write(kfid,*) '2e INT time(libcint) :', ted2-tst2 ,'s'   !!! count time 2          
        else
!!!  1e Int
          if(Moutlevel>2)then
                lprint1e=1
          end if
          if(Moutlevel<4)then
                lprint2e=0
          else if(Moutlevel==4)then
                lprint2e=1
          else if(Moutlevel==5)then
                lprint2e=2
          end if 
          if((MethodS+MethodT+MethodV)/=0) then
                call STVINT(nshell,nbasis,natom,coor,KATOM,&
                        kcenter,gfc,gfe,ntype,ngf,ngfall,&
                        kbasst,kbased,knbas,kpgfst,kpgfed,naxyz,&
                        S,T,V,Hcore,dipx,dipy,dipz,MethodS,MethodT,MethodV)
                call Oneintout(kfid,nbasis,methods,methodT,methodV,lprint1e,&
                        S,T,V,Hcore,dipx,dipy,dipz)
           end if
!!!! 2e Int            
           if(methodERI/=0) then
                tst=counttime_cputime(t_start) !!! count time 1
                tst2=counttime_sysckock(it1)   !!! count time 2
!!                
                call ERINT(nshell,nbasis,natom,coor,neri,kcenter,gfc,gfe,ntype,ngf,&
                        ngfall,kbasst,kbased,knbas,kpgfst,kpgfed,naxyz,ERIall,methodERI)
!!                       
                ted=counttime_cputime(t_start) !!! count time 1
                ted2=counttime_sysckock(it1)   !!! count time 2
                write(kfid,*) '2e INT time :', ted-tst ,'s'     !!! count time 1
                write(kfid,*) '2e INT time :', ted2-tst2 ,'s'   !!! count time 2
!!          
                call TWOintout(kfid,nbasis,neri,methodERI,lprint2e,ERIall)
           end if
        end if  !!! use libcint
!!!!! SCF
        if(mhf==1)then
           tst=counttime_cputime(t_start) !!! count time 1
           tst2=counttime_sysckock(it1)   !!! count time 2
!!!!           
           call RHF(kfid,maxcycle,nbasis,nele,neri,10**(-conver_E),10**(-conver_rms),S,T,V,&
                hcore,damp,mguess,mdiis,ndiis,ERIall,Enuc,Eelec,Cof,eorb,X)
           Efinal=Eelec
!!!!                        
           ted=counttime_cputime(t_start) !!! count time 1
           ted2=counttime_sysckock(it1)   !!! count time 2
           write(kfid,*) 'SCF time :', ted-tst ,'s'     !!! count time 1
           write(kfid,*) 'SCF time :', ted2-tst2 ,'s'   !!! count time 2
        end if
        call PopAny(kfid,natom,coor,nshell,nbasis,nele,KATOM,kcenter,kbasst,kbased,&
                cof,S,dipx,dipy,dipz,X,usecint)    
!!!! MP2,MP3
!!!!     call  AO2MOERI(nbasis,neri,ERIall,Cof,ERIMO) !!! slow
        if(mmp>=1)then
!!       
           tst=counttime_cputime(t_start) !!! count time 1
           tst2=counttime_sysckock(it1)   !!! count time 2
!!          
           call  AO2MOERIfast(nbasis,neri,hcore,ERIall,Cof,HcoreMO,ERIMO)
!!           
           ted=counttime_cputime(t_start) !!! count time 1
           ted2=counttime_sysckock(it1)   !!! count time 2
           write(kfid,*) 'AO to MO INT time :', ted-tst ,'s'     !!! count time 1
           write(kfid,*) 'AO to MO INT time :', ted2-tst2 ,'s'   !!! count time 2
!!!           
           tst=counttime_cputime(t_start) !!! count time 1
           tst2=counttime_sysckock(it1)   !!! count time 2
!!           
           call  MP2(kfid,nbasis,neri,Nele,ERIMO,Eelec,Eorb,Emp2)
!!          
           ted=counttime_cputime(t_start) !!! count time 1
           ted2=counttime_sysckock(it1)   !!! count time 2
           write(kfid,*) 'MP2 time :', ted-tst ,'s'     !!! count time 1
           write(kfid,*) 'MP2 time :', ted2-tst2 ,'s'   !!! count time 2 
           Efinal=Efinal+Emp2
        end if  !! MP2
        if(mmp>=2)then
!!!           
           tst=counttime_cputime(t_start) !!! count time 1
           tst2=counttime_sysckock(it1)   !!! count time 2
!!                   
           call MP3(kfid,mmp,nbasis,neri,Nele,ERIMO,Eelec,Eorb,Emp2,Emp3)
!!          
           ted=counttime_cputime(t_start) !!! count time 1
           ted2=counttime_sysckock(it1)   !!! count time 2
           write(kfid,*) 'MP3 time :', ted-tst ,'s'     !!! count time 1
           write(kfid,*) 'MP3 time :', ted2-tst2 ,'s'   !!! count time 2 
           Efinal=Efinal+Emp3
        end if !! MP3
!!!! CC
        if(mcc>0)then
            allocate(ERISPIN(2*nbasis,2*nbasis,2*nbasis,2*nbasis))
!!       
           tst=counttime_cputime(t_start) !!! count time 1
           tst2=counttime_sysckock(it1)   !!! count time 2           
!!            
           call AO2MOERIfast(nbasis,neri,hcore,ERIall,Cof,HcoreMO,ERIMO)
!!           
           ted=counttime_cputime(t_start) !!! count time 1
           ted2=counttime_sysckock(it1)   !!! count time 2
           write(kfid,*) 'AO to MO INT time :', ted-tst ,'s'     !!! count time 1
           write(kfid,*) 'AO to MO INT time :', ted2-tst2 ,'s'   !!! count time 2     
           tst=counttime_cputime(t_start) !!! count time 1
           tst2=counttime_sysckock(it1)   !!! count time 2   
!!           
           call MO2SOERI(nbasis,neri,ERIMO,ERISPIN)
!!           
           ted=counttime_cputime(t_start) !!! count time 1
           ted2=counttime_sysckock(it1)   !!! count time 2
           write(kfid,*) 'MO to spin-orbital time :', ted-tst ,'s'     !!! count time 1
           write(kfid,*) 'MO to spin-orbital time :', ted2-tst2 ,'s'   !!! count time 2
           call CCSDT(kfid,nbasis,neri,Nele,HcoreMO,ERISPIN,Eelec,Eorb,mcc,&
                ccdamp,maxcyclecc,conver_Ecc,conver_rmscc,mdiiscc,ndiiscc,Emp2,Ecc)
           Efinal=Eelec+Ecc
!!          
           ted=counttime_cputime(t_start) !!! count time 1
           ted2=counttime_sysckock(it1)   !!! count time 2
           write(kfid,*) 'CCSD time :', ted-tst ,'s'     !!! count time 1
           write(kfid,*) 'CCSD time :', ted2-tst2 ,'s'   !!! count time 2 
        end if    !! CC    
!!!        
        write(*,*)'*** Computing done ***'    
!!!!!!!!!!!!!!!==================== OUTPUT =====================!!!!!!!!!!!!
    write(*,*)'*** OUTPUT ***'
!!!! Basical information for debug.
    write(*,*)'nele',nele   !! number of all electrons(colsed shell)
    write(*,*)'nbasis',nbasis   !! number of all basis function
    write(*,*)'kbasst',kbasst   !! The  first number of basis function in one shell
    write(*,*)'kbased',kbased   !! The  last number of basis function in one shell
    write(*,*)'knbas',knbas     !! number of  basis function in one shell
    write(*,*)'lgfst',lgfst     !! pgf type first and last number 
    write(*,*)'lgfed',lgfed     !! s=1,p=2:4,d=5:10... get type (eq:000,100..)
    write(*,*)'kpgfst',kpgfst   !! pgf first number of shell
    write(*,*)'kpgfed',kpgfed   !! pgf last number of shell
!     write(*,*)'N:SPDFGHI',nss,nps,NDS,nfs,ngs,nhs,nis,nsp  !! Number of SPDFGHI shell
!     write(*,*)'norb',norb       !! (6D->5D) number of  basis function(Gaussian)
    write(*,'(A6,3x,4(a6,6x))')'naxyz:','basis','Type x','Type y','Type z'
    do i=1,nbasis
        write(*,*)i,naxyz(i),naxyz(i+nbasis),naxyz(i+2*nbasis)
    end do
    ! write(*,'(25(I3),/)')naxyz
    write(*,*)'neri=',neri   !! number of stored ERI
    write(*,*)'neri(nozero)=' ,count(abs(ERIall)>1d-12) !! number of stored nozero ERI

!!!!!! Summary 
    write(*,'(25X,A30)')'***** Summary *****'
    write(*,*) 'Hartree-Fock Final results (a.u.)'
    write(*,'(25X,A11,5X,F20.15)')'E(nuc):',Enuc
    write(*,'(25X,A11,5X,F20.15)')'E(HF_elec):',Eelec-Enuc
    write(*,'(25X,A11,5X,F20.15)')'E(HF_tol):',Eelec
     write(*,'(25X,A11,5X,F20.15)')'E(Final):',Efinal
! !!!!!! Output time to OUTPUT
    call date_and_time(theDate,t2_str) 
    call cpu_time(t_end)
    call system_Clock(it2,irate,icount)
    call outputtime
!!!!============================ END =================================!!!!
!!!!    
        deallocate(kbasst,kbased,knbas,lgfst,lgfed,kpgfst,kpgfed,naxyz,&
                    S,T,V,Hcore,dipx,dipy,dipz,ERIall,Cof,eorb,X,ERIMO,&
                   HcoreMO)
        deallocate(coor,KATOM,kcenter,gfc,gfe,ntype,ngf)
        if(mcc>0)then
                deallocate(ERISPIN)
        end if
        close(99)
        write(*,*)'*** Program terminated normally ***'
    end program


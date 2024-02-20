!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                         HartreeFock.f90                            !!  
!!      The Hartree-Fock self-consistent field (SCF) procedure        !!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!!!!! This subroutine uses lapack.
    module HF
        contains
    subroutine RHF(kfid,maxcycle,nbasis,nele,neri,Ethreshold,rmsthreshold,S,T,V,Hcore,&
                    damp,Mguess,mdiis,ndiis,&
                    ERIall,Enuc,EHF,Cof,eorb,X)
        use Global
        implicit double precision (A-H,O-Z)
        Integer*16 neri
        dimension::S(nbasis,nbasis),T(nbasis,nbasis),V(nbasis,nbasis),&
                    ERIall(neri),Cof(nbasis,nbasis)
        dimension::Hcore(nbasis,nbasis),eigv(nbasis),X(nbasis,nbasis),&
                Sem12(nbasis,nbasis),Fock(nbasis,nbasis),Tem(nbasis,nbasis),eorb(nbasis),&
                Den(nbasis,nbasis),Fock0(nbasis,nbasis),Cof0(nbasis,nbasis),&
                Denold(nbasis,nbasis),GG(nbasis,nbasis),temmm1(nbasis,nbasis),guess(nbasis,nbasis)
!!!!! for diis           
        allocatable::errmat(:,:,:),fockdiis(:,:,:)     
        dimension:: errmatnow(nbasis,nbasis),tempdiis(nbasis,nbasis),&
                    Bdiis(ndiis+1,ndiis+1),Cofdiis(ndiis+1),ipivdiis(ndiis+1)
!!!!! for dsyev and dsyevd                
        dimension::work(3*nbasis+1),work1(1+6*nbasis+2*nbasis**2),work2(3+5*nbasis)
!!!!!!     
        nocc=nele/2
        if(nocc==0)nocc=1
        write(kfid,*)'=======================================&
                  =========================================&
                  ============================='  
        write(*,*)'*** Start self-consistent field (SCF) procedure ***' 
        write(kfid,'(15X,A60)')'***  Self-consistent field (SCF) procedure ***' 
        fock=0d0
        fock0=0d0
        den=0d0
        cof0=0d0
        rmsold=0d0
        guess=0d0
        Cofdiis=0d0
!!!!  Diagonalize S,X=S^(-1/2)
!!!!!!!  dsyev is from lapack.
        Tem=S
        call dsyev('V','L',nbasis,Tem,nbasis,eigv,work,3*nbasis+1,info1)  
        ! call dsyevd('V','L',nbasis,Tem,nbasis,eigv,work1,1+6*nbasis+2*nbasis**2,work2,3+5*nbasis,info1)  
        if(info1/=0)stop " Diagonalize S failed !!!!! "
        Sem12=0d0
        do i=1,nbasis
            Sem12(i,i)=eigv(i)**(-0.5d0)
        end do   
!!!! Built-in functions  
        ! X=matmul(tem,Sem12)
        ! X=matmul(X,transpose(tem))
!!!! BLAS        
        call dgemm('N','N',nbasis,nbasis,nbasis,1d0,tem,nbasis,Sem12,nbasis,0d0,temmm1,nbasis)
        call dgemm('N','T',nbasis,nbasis,nbasis,1d0,temmm1,nbasis,tem,nbasis,0d0,X,nbasis)
!!!!!! for debug
        ! write(kfid,*)'s'
        ! call outputmatrixd(nbasis,nbasis,s,kfid)   
        ! write(kfid,'(10X,80A)')'------------------------------------& 
        !           ---------------------------------------'         
        ! write(kfid,*)'X=S^-1/2'
        ! call outputmatrixd(nbasis,nbasis,x,kfid)
        ! write(kfid,*)'hcore'
        ! call outputmatrixd(nbasis,nbasis,Hcore,kfid)
!!!!  Build the Initial Guess Density using Hcore.
        if(mguess==1)then
!!! core Hamiltonian   
            guess= Hcore   
            write(kfid,*)'Guess: core Hamiltonian'
        else if(mguess==2)then
!!! generalized Wolfsberg−Helmholz (GWH) 
!!! guess(i,j)=0.5*K*(Hcore(i,i)+Hcore(j,j))*S(i,j),K=1.75 typically
            write(kfid,*)'Guess: GWH'
            do i=1,nbasis
                do j=1,i
                    if(i==j)then
                        guess(i,j)=hcore(i,i)
                    else
                       guess(i,j)=0.5d0*1.75*S(i,j)*(Hcore(i,i)+Hcore(j,j)) 
                       guess(j,i)=guess(i,j)
                    end if    
                end do
            end do
!!!!            
        else if(mguess==3)then  
        write(kfid,*)'Guess: Huckel'
            tem=hcore
            call dsyevd('V','L',nbasis,Tem,nbasis,eigv,work1,1+6*nbasis+2*nbasis**2,work2,3+5*nbasis,info2)  
            do i=1,nbasis
                do j=1,i
                    if(i==j)then
                        guess(i,j)=-eigv(i)
                    else
                       guess(i,j)=0.5d0*1.75*S(i,j)*(Hcore(i,i)+Hcore(j,j)) 
                       guess(j,i)=guess(i,j)
                    end if    
                end do
            end do            
        end if
        Fock=guess
        ! Fock0=matmul(transpose(X),Fock)
        ! Fock0=matmul(Fock0,X)
!!!!  BLAS
        call dgemm('T','N',nbasis,nbasis,nbasis,1d0,X,nbasis,fock,nbasis,0d0,temmm1,nbasis)
        call dgemm('N','N',nbasis,nbasis,nbasis,1d0,transpose(X),nbasis,temmm1,nbasis,0d0,fock0,nbasis)        
        ! write(kfid,*)'FOCK0'
        ! call outputmatrixd(nbasis,nbasis,Fock0,kfid)
        tem=Fock0
        ! call dsyev('V','L',nbasis,Tem,nbasis,eigv,work,3*nbasis+1,info2)  
        call dsyevd('V','L',nbasis,Tem,nbasis,eigv,work1,1+6*nbasis+2*nbasis**2,work2,3+5*nbasis,info2)  

        if(info2/=0)stop " Diagonalize F0 failed !!!!! "
!!!!        
        ! Cof0=matmul(X,Tem)
        call dgemm('N','N',nbasis,nbasis,nbasis,1d0,x,nbasis,tem,nbasis,0d0,Cof0,nbasis)
!!!!!! for debug
        ! write(kfid,*)'cof0'
        ! call outputmatrixd(nbasis,nbasis,cof0,kfid)
!!!!
!!! Density martrix
!!
        ! Den=matmul(cof0(:,1:nocc),transpose(cof0(:,1:nocc)))  
!!        
        ! call dgemm('N','T',nbasis,nbasis,nbasis,1d0,temmm1,nocc,tem,nocc,0d0,X,nbasis)
!!        
        cof=transpose(cof0)
        do mm=1,nbasis
            do nn=1,nbasis
                do kocc=1,nocc
                    Den(mm,nn)=Den(mm,nn)+Cof0(mm,kocc)*Cof(nn,kocc)
                end do
            end do
        end do 
!!!!        
        Eele=0d0
        do mm=1,nbasis
            do nn=1,nbasis
                Eele=Eele+Den(mm,nn)*(Hcore(mm,nn)+Fock(mm,nn))
            end do
        end do 
        if(mdiis==1)write(kfid,*)'DIIS: on','       ','Number of diis matrix: ',ndiis
        if(damp/=0d0)write(kfid,'(1X,A5,1x,f4.2)')'damp:',damp
            write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------'      
            write(kfid,'(35X,40A)')'*****************'
            write(kfid,'(35X,40A)')'*SCF Iterations *'   
            write(kfid,'(35X,40A)')'*****************'                   
            write(kfid,*)'E0=',Eele,'Enuc=',Enuc,'Etol0=',Eele+Enuc
!!!! SCF 
!!!!!! build new Fock matrix
        write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------'     
        write(kfid,'(2X,A5,3X,6(A10,13X))')'Cycle','E','E(old)','delt E','rms'
        write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------'     
        if(mdiis==1)then
            allocate(errmat(nbasis,nbasis,ndiis),fockdiis(nbasis,nbasis,ndiis))
        end if   
        do kcy=1,maxcycle
            write(*,*)'iteration num :',kcy
            Eold=Eele
            Denold=Den
            rmsold=rms
            Fock=0d0
            GG=0d0
            Cof0=0d0
            Cofdiis=0d0
            Bdiis=0d0
            do ii=1,nbasis
                do jj=1,nbasis
                    Fock(ii,jj)=Hcore(ii,jj)
                    do kk=1,nbasis
                        do ll=1,nbasis
                        ijkl=indexeri(ii,jj,kk,ll)
                        ikjl=indexeri(ii,kk,jj,ll)
                        Fock(ii,jj)=Fock(ii,jj)+Den(kk,ll)*2*(ERIall(ijkl)-0.5d0*ERIall(ikjl))                        
                        ! GG(ii,jj)=GG(ii,jj)+2*Den(kk,ll)*(ERIall(ijkl)-0.5d0*ERIall(ikjl))
                        end do
                    end do
                end do
            end do
            ! Fock=Hcore+GG
!!!!!   DIIS 
            if(mdiis==1)then
                fockdiis(:,:,mod(kcy,ndiis)+1)=fock
!!!! set error matrix             
!!    FiDiS                
                tempdiis=matmul(den,S)
                errmatnow=matmul(fock,tempdiis)
                errmat(:,:,mod(kcy,ndiis)+1)=errmatnow
!!    SDiFi
                tempdiis=matmul(den,fock)
                errmatnow=matmul(S,tempdiis)
!!
                errmat(:,:,mod(kcy,ndiis)+1)=errmat(:,:,mod(kcy,ndiis)+1)-errmatnow            
!!!!   set B matrix         
                if(kcy>=ndiis)then
                    Bdiis=-1d0
                    do  kp=1,ndiis
                        do kq=1,kp
                            sumdiis=0d0
                            do ii=1,nbasis
                                do jj=1,nbasis
                                sumdiis=sumdiis+errmat(ii,jj,kp)*errmat(jj,ii,kq)
                                end do  !! jj
                            end do      !! ii
                            Bdiis(kp,kq)=sumdiis
                            if(kp/=kq)Bdiis(kq,kp)=sumdiis
                        end do    !! kq
                    end do        !! kp            
                    Bdiis(ndiis+1,ndiis+1)=0d0
                    Cofdiis(ndiis+1)=-1
!!!!    solve C  
                    call dgesv(ndiis+1,1,Bdiis,ndiis+1,ipivdiis,Cofdiis,ndiis+1,infodiis)
                    if(infodiis/=0)stop 'Solve diis cof failed!!!'
!!!!  set new fock          
                    fock=0d0
                    do iii=1,ndiis
                        fock=fock+Cofdiis(iii)*fockdiis(:,:,iii)         
                    end do       
                end if 
            end if   !!!  diis     
!!!!!!! Diagonalize F'
            Fock0=matmul(matmul(transpose(X),Fock),X)
            ! Fock0=matmul(matmul(X,Fock),X)
            tem=Fock0
            eigv=0d0
            ! call dsyev('V','L',nbasis,Tem,nbasis,eigv,work,3*nbasis+1,info)  
            call dsyevd('V','L',nbasis,Tem,nbasis,eigv,work1,1+6*nbasis+2*nbasis**2,work2,3+5*nbasis,info)  
            if(info/=0)stop " Diagonalize F' failed !!!!! " 
!!!!!!! build new Density matrix
            Cof0=matmul(X,tem)
            Den=0d0
!!!! Method 1  compute Den            
            ! Den=matmul(cof0(:,1:nocc),transpose(cof0(:,1:nocc)))
!!!! Method 2  compute Den            
            do mm=1,nbasis
                do nn=1,nbasis
                    do kocc=1,nocc
                        Den(nn,mm)=Den(nn,mm)+Cof0(mm,kocc)*Cof0(nn,kocc)
                    end do
                end do
            end do 
!!!! damping
            Den=(1d0-damp)*den+damp*denold            
!!!!            
            Eele=0d0
            ET=0d0
            Ecore=0d0
            EV=0d0
            EG=0d0
            rms=0d0
            do mm=1,nbasis
                do nn=1,nbasis
                    Eele=Eele+Den(mm,nn)*(Hcore(mm,nn)+Fock(mm,nn))
                    Et=Et+2*Den(mm,nn)*(T(mm,nn))
                    Ev=Ev+2*Den(mm,nn)*(-V(mm,nn))
                    Ecore=Ecore+Den(mm,nn)*2*Hcore(mm,nn)     
                    EG=EG+Den(mm,nn)*(-Hcore(mm,nn)+Fock(mm,nn))             
                    rms=rms+sqrt((Den(mm,nn)-Denold(mm,nn))**2)
                end do
            end do 
            write(kfid,'(I4,5(F23.16))')kcy,Eele+Enuc,Eold+Enuc,Eele-Eold,rms      
            if(abs(Eele-Eold)<=Ethreshold .and.abs(rms)<=rmsthreshold)then
                exit
            else if(kcy==maxcycle)then
                write(kfid,*) 'Fail to convergence of SCF!'
                write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------'      
                write(kfid,*) 'SCF Done' 
                write(kfid,'(2x,"After",I4,2x,"cycles")') kcy         
                write(kfid,*)' Final Hartree-Fock Energy :',Ehf,'A.U.'
                write(kfid,*)' The SCF convergence threshold of E :',Ethreshold,'A.U.'
                write(kfid,*)' The SCF convergence threshold of RMS :',rmsthreshold,'A.U.'
                write(kfid,*) 'Final results (A.U.)'
                write(kfid,'(25X,A11,5X,F20.15)')'E(HF_tol):',Enuc+Eele
                        stop 'Fail to convergence of SCF!'
            end if
        end do !!! cycle
!!!!!!!!!!!!!  Outoput some information        
        Ehf=Enuc+Eele
        write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------'      
        write(kfid,*) 'SCF Done' 
        write(kfid,'(2x,"After",I4,2x,"cycles")') kcy         
        write(kfid,*)' Final Hartree-Fock Energy :',Ehf,'A.U.'
        write(kfid,*)' The SCF convergence threshold of E :',Ethreshold,'A.U.'
        write(kfid,*)' The SCF convergence threshold of RMS :',rmsthreshold,'A.U.'
!!!!  output  coefficients        
        write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------------'   
        write(kfid,*)'Molecular orbital combination coefficients C :'
        cof=Cof0
        eorb=eigv
        hmgap=eorb(1+nocc)-eorb(nocc)
        call outputmatrixd(nbasis,nbasis,Cof,kfid)
!!!!  output  Density matrix     
        ! write(kfid,*)'------------------------------------& 
        !         -----------------------------------------&
        !         ---------------------------'          
        ! write(kfid,*)'Density matrix :'
        ! call outputmatrixd(nbasis,nbasis,2*Den,kfid)  
        write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------------'                  
        write(kfid,*)'orbital energy (EigenValues):'
        write(kfid,*)'occupied orbitals :',1,'~',nocc
        write(kfid,'(5(f20.12))')eorb(1:nocc)
        write(kfid,*)'virtual orbitals :',1+nocc,'~',nbasis
        write(kfid,'(5(f20.12))')eorb(1+nocc:nbasis)   
        write(kfid,*)'HOMO :',eorb(nocc),'index :',nocc
        write(kfid,*)'LOMO :',eorb(1+nocc),'index :',1+nocc

        write(kfid,'(" HOMO-LUMO gap :",f20.12,"  Eh")')hmgap
        write(kfid,'(15X,"=",f20.12,"  eV")')hmgap*Eh2eV
        write(kfid,'(15X,"=",f20.12,"  kcal/mol")')hmgap*Eh2eV*eV2kcalpm
        write(kfid,'(15X,"=",f20.12,"  kJ/mol")')hmgap*Eh2eV*eV2kcalpm*kcal2kJ        
!!!!!! Summary 
        write(kfid,*)'------------------------------------& 
                -----------------------------------------&
                ---------------------------'    
        write(kfid,'(25X,A30)')'*****  SCF Summary *****'
        write(kfid,*)'after',kcy,'cycles'
        write(kfid,*) 'Final results (A.U.)'
        write(kfid,'(25X,A11,5X,F20.15)')'E(HF_tol):',Enuc+Eele
        write(kfid,*)'Energy analysis: '
        write(kfid,'(25X,A11,5X,F20.15)')'E(nuc):',Enuc
        write(kfid,'(25X,A11,5X,F20.15)')'E(HF_elec):',Eele
        write(kfid,'(25X,A11,5X,F20.15)')'E(Kinetic):',ET
        write(kfid,'(25X,A11,5X,F20.15)')'E(V1e):',Ev
        write(kfid,'(25X,A11,5X,F20.15)')'E(core)1e:',Ecore
        write(kfid,'(25X,A11,5X,F20.15)')'E(V2e):',EG 
        write(kfid,'(25X,A11,5X,F20.15)')'E(Vtot):',Ev+EG+Enuc
        write(kfid,'(25X,A11,5X,F20.15)')'V/T:',-(EG+Ev+Enuc)/ET          
!!!!!!!!!!!!!!!!!!!
        if(mdiis==1)then
            deallocate(errmat,fockdiis)
        end if   
        write(*,*)'*** The self-consistent field (SCF) procedure Done ***'
    end subroutine
!!!!!!===============================================================!!!!!
    subroutine PopAny(kfid,natom,coor,nshell,nbasis,nele,KATOM,kcenter,kbasst&
                            ,kbased,cof,S,dipx,dipy,dipz,X,usecint)
        use Global
        implicit double precision (A-H,O-Z)    
        dimension:: Pop(nbasis,nbasis),cof(nbasis,nbasis),katom(natom),&
            den(nbasis,nbasis),kcenter(nshell),kbasst(nshell),&
            kbased(nshell),S(nbasis,nbasis),cof0(nbasis,nbasis),ZPop(natom),&
            dipx(nbasis,nbasis),dipy(nbasis,nbasis),dipz(nbasis,nbasis),&
            X(nbasis,nbasis),temp(nbasis,nbasis),coor(natom,3)
        logical usecint
!!!!  for lapack
        dimension::ipiv(nbasis),work(nbasis),work1(3*nbasis+1),&
        sem12(nbasis,nbasis),dipiv(nbasis)
!!!!!            
        den=0d0   
        nocc=nele/2
        if(nocc==0)nocc=1
        cof0=transpose(cof)
!!!! Method 1  compute Den            
        Den=matmul(cof(:,1:nocc),transpose(cof(:,1:nocc))) 
!!!! Method 2  compute Den            
        ! do mm=1,nbasis
        !     do nn=1,nbasis
        !         do kocc=1,nocc
        !             Den(mm,nn)=Den(mm,nn)+Cof(mm,kocc)*Cof(nn,kocc)
        !             ! Den(mm,nn)=Den(mm,nn)+Cof(mm,kocc)*Cof0(kocc,nn)
        !         end do
        !     end do
        ! end do 
!!!!  for debug              
        ! write(99,*)'S'
        ! call outputmatrixd(nbasis,nbasis,S,99)
        ! write(99,*)'Density '
        ! call outputmatrixd(nbasis,nbasis,den,99)
        ! write(99,*)'Density *2'
        ! call outputmatrixd(nbasis,nbasis,2*den,99)
!!!!=====================  Population Analysis/Atomic Charges  ==========!!!!! 
!!!!!!!!!!!!!!!!!
        Pop=matmul(Den,S) 
        Ntem=1
        tem=0d0
        do i=1,nshell
            if(kcenter(i)==ntem)then
                do j=kbasst(i),kbased(i)
                    tem=tem+2*pop(j,j)
                end do
            else  
                ZPop(kcenter(i-1))=katom(kcenter(i-1))-tem
                ntem=kcenter(i)
                tem=0d0
                do j=kbasst(i),kbased(i)
                    tem=tem+2*pop(j,j)
                end do
            end if
            if(i==nshell)then
                ZPop(kcenter(i))=katom(kcenter(i))-tem
           end if
        end do
        write(kfid,*)'=======================================&
                  ========================================='        
        write(kfid,'(20x,50A)')'*** Population Analysis/Atomic Charges ***'
        write(kfid,*)'------------------------------------& 
                -----------------------------------------'          
        write(kfid,'(20X,A14,2x,A30)')'Atomic charges','Mulliken charges on atom'
        do i=1,natom
            write(kfid,*)'ATOM',i,katom(i)-ZPop(i),ZPop(i)
        end do 
        ZMu=sum(ZPop)
        if(abs(ZMu)<1d-8)zmu=0d0
        write(kfid,*)'Sum of Mulliken charges =  ',zmu
!!!!!!!!!!!!!!!!!!!!!!!! Lowdin charges
        tem=0d0
        ZPop=0d0
!!!!!!!!!!!!  method 1 :  S -> S^1/2  
        ! Temp=S
        ! call dsyev('V','L',nbasis,Temp,nbasis,dipiv,work1,3*nbasis+1,info1)  
        ! if(info1/=0)stop " Diagonalize S failed !!!!! "
        ! Sem12=0d0
        ! do i=1,nbasis
        !     Sem12(i,i)=dipiv(i)**(0.5d0)
        ! end do        
        ! pop=matmul(temp,Sem12) 
        ! temp=matmul(pop,transpose(temp))
!!!!!!!!!!!!  method 2 :   X=S^-1/2 -> S^1/2          
!!!!!    LU(X) -> inv(X)
        temp=X        
        call dgetrf(nbasis,nbasis,temp,nbasis,ipiv,info)
        call dgetri(nbasis,temp,nbasis,ipiv,work,nbasis,info)
!!
        Pop=matmul(Den,temp) 
        Pop=matmul(temp,Pop)
        Ntem=1
        tem=0d0
        do i=1,nshell
            if(kcenter(i)==ntem)then
                do j=kbasst(i),kbased(i)
                    tem=tem+2*pop(j,j)
                end do
            else  
                ZPop(kcenter(i-1))=katom(kcenter(i-1))-tem
                ntem=kcenter(i)
                tem=0d0
                do j=kbasst(i),kbased(i)
                    tem=tem+2*pop(j,j)
                end do
            end if
            if(i==nshell)then
                ZPop(kcenter(i))=katom(kcenter(i))-tem
           end if
        end do
        write(kfid,*)'------------------------------------& 
                -----------------------------------------'          
        write(kfid,'(20X,A14,2x,A30)')'Atomic charges','Lowdin charges on atom'
        do i=1,natom
            write(kfid,*)'ATOM',i,katom(i)-ZPop(i),ZPop(i)
        end do 
        ZMu=sum(ZPop)
        if(abs(ZMu)<1d-8)zmu=0d0
        write(kfid,*)'Sum of Lowdin charges =  ',zmu    
!!!!========================  Dipole moment    =========================!!!!! 
        if(.not. usecint)then
        write(kfid,*)'=======================================&
                  ========================================='        
        write(kfid,'(28x,50A)')'***  Dipole moment ***'
        write(kfid,*)'------------------------------------& 
                -----------------------------------------'  
        ! temp=matmul(den,dipx)
        ! ux=sum(temp(:,1:nocc))
        ! temp=matmul(den,dipy)
        ! uy=sum(temp(:,1:nocc))
        ! temp=matmul(den,dipz)
        ! uz=sum(temp(:,1:nocc))
        ux=0d0
        uy=0d0
        uz=0d0
        do i=1,nbasis
            do j=1,nbasis
                ux=ux+2*den(i,j)*dipx(i,j)
                uy=uy+2*den(i,j)*dipy(i,j)
                uz=uz+2*den(i,j)*dipz(i,j)
            end do
        end do
        zrx=0d0
        zry=0d0
        zrz=0d0
        do i=1,natom
            zrx=zrx+katom(i)*coor(i,1)
            zry=zry+katom(i)*coor(i,2)
            zrz=zrz+katom(i)*coor(i,3)
        end do
        ux=-ux+zrx
        uy=-uy+zry
        uz=-uz+zrz
        utotl=sqrt(ux**2+uy**2+uz**2)
        if(abs(ux)<1d-8)ux=0d0
        if(abs(uy)<1d-8)uy=0d0
        if(abs(uz)<1d-8)uz=0d0
        if(abs(utotl)<1d-8)utotl=0d0
        write(kfid,'((8X,"X(A.U.)",20X,"Y(A.U.)"),20X,"Z(A.U.)")')
        write(kfid,*)ux,uy,uz
        write(kfid,'((8X,"X(Debye)",19X,"Y(Debye)"),19X,"Z(Debye)")')
        write(kfid,*)ux*au2debye,uy*au2debye,uz*au2debye
        write(kfid,'(4x,"Tot = ",F16.10," A.U.")')utotl
        write(kfid,'(7X," = ",F16.10," Debye")')utotl*au2debye
        end if 
    end subroutine PopAny
!!!!!!!!!!!!!!!!!!!!!!!!!!       
    End Module HF
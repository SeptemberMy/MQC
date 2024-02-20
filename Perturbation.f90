!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                         Perturbation.f90                           !!  
!!                   Perturbation Theory: MP2 MP3                     !!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    subroutine MP2(kfid,nbasis,neri,Nele,ERIMO,Eelec,Eorb,Emp2)
        use Global,only:indexeri,index2eri
        implicit double precision (A-H,O-Z)
        Integer*16 neri
        dimension::Eorb(nbasis),ERIMO(neri)
!!!!!!!!        
        write(*,*)'*** Start MP2 procedure ***' 
        tem=0d0
        Nocc=nele/2
        if(nocc==0)nocc=1
!!!!        
        do ii=1,Nocc
            do ia=Nocc+1,nbasis
                do ij=1,Nocc
                    do ib=Nocc+1,nbasis
                        iajb=indexeri(ii,ia,ij,ib)
                        ibja=indexeri(ii,ib,ij,ia)
                        tem=tem+ERIMO(iajb)*(2d0*ERIMO(iajb)-ERIMO(ibja))/(Eorb(ii)+Eorb(ij)-Eorb(ia)-Eorb(ib))
                    end do 
                end do 
            end do 
        end do
!!!!!        
        Emp2=tem
        write(kfid,*)'=======================================&
                  ========================================='
        write(kfid,'(25X,A30)')'*****  MP2 Summary *****'
        write(kfid,*) 'MP2 Final results (a.u.)'
        write(kfid,'(25X,A11,5X,F20.15)')'E(HF_tol):',Eelec
        write(kfid,'(25X,A11,5X,F20.15)')'E(MP2):',Emp2
        write(kfid,'(25X,A11,5X,F20.15)')'E(HF+MP2):',Emp2+Eelec
        write(*,*)'*** The MP2 procedure Done ***'        
    end subroutine MP2
!!!!===============================================================!!!!!!!!!!
    subroutine MP3(kfid,mmp,nbasis,neri,Nele,ERIMO,Eelec,Eorb,Emp2,Emp3)
        use Global,only:indexeri,index2eri
        implicit double precision (A-H,O-Z)
        Integer*16 neri
        dimension::Eorb(nbasis),ERIMO(neri)
!!!!!!!!        
        write(*,*)'*** Start MP3 procedure ***' 
        tem1=0d0
        tem2=0d0
        tem3=0d0
        Nocc=nele/2
        if(nocc==0)nocc=1
!!!!  part 1      
        do ii=1,Nocc
            do ia=Nocc+1,nbasis
                do ij=1,Nocc
                    do ib=Nocc+1,nbasis
                        do ik=1,Nocc
                            do il=1,Nocc
                        iajb=indexeri(ii,ia,ij,ib)
                        ibja=indexeri(ii,ib,ij,ia)
                        ikjl=indexeri(ii,ik,ij,il)
                        iljk=indexeri(ii,il,ij,ik)
                        kalb=indexeri(ik,ia,il,ib)
                        kbla=indexeri(ik,ib,il,ia)
                        tem1now=4d0*(ERIMO(iajb)*ERIMO(ikjl)*ERIMO(kalb)&
                                    +ERIMO(ibja)*ERIMO(iljk)*ERIMO(kalb)&
                                    +ERIMO(ibja)*ERIMO(ikjl)*ERIMO(kbla)&
                                    +ERIMO(iajb)*ERIMO(iljk)*ERIMO(kbla))
                        tem1now=tem1now-2d0*&
                                    (ERIMO(ibja)*ERIMO(ikjl)*ERIMO(kalb)&
                                    +ERIMO(iajb)*ERIMO(iljk)*ERIMO(kalb)&
                                    +ERIMO(iajb)*ERIMO(ikjl)*ERIMO(kbla)&
                                    +ERIMO(ibja)*ERIMO(iljk)*ERIMO(kbla) )
                        tem1now=tem1now/(Eorb(ii)+Eorb(ij)-Eorb(ia)-Eorb(ib))
                        tem1now=tem1now/(Eorb(ik)+Eorb(il)-Eorb(ia)-Eorb(ib))
                        tem1=tem1+tem1now
                            end do 
                        end do 
                    end do 
                end do 
            end do 
        end do
!!!!  part 2      
        do ii=1,Nocc
            do ia=Nocc+1,nbasis
                do ij=1,Nocc
                    do ib=Nocc+1,nbasis
                        do ic=Nocc+1,nbasis
                            do id=Nocc+1,nbasis
                        iajb=indexeri(ii,ia,ij,ib)
                        ibja=indexeri(ii,ib,ij,ia)
                        iacbd=indexeri(ia,ic,ib,id)
                        iadbc=indexeri(ia,id,ib,ic)
                        icjd=indexeri(ii,ic,ij,id)
                        idjc=indexeri(ii,id,ij,ic)
                        tem2now=4d0*(ERIMO(iajb)*ERIMO(iacbd)*ERIMO(icjd)&
                                    +ERIMO(ibja)*ERIMO(iadbc)*ERIMO(icjd)&
                                    +ERIMO(ibja)*ERIMO(iacbd)*ERIMO(idjc)&
                                    +ERIMO(iajb)*ERIMO(iadbc)*ERIMO(idjc))
                        tem2now=tem2now-2d0*&
                                    (ERIMO(ibja)*ERIMO(iacbd)*ERIMO(icjd)&
                                    +ERIMO(iajb)*ERIMO(iadbc)*ERIMO(icjd)&
                                    +ERIMO(iajb)*ERIMO(iacbd)*ERIMO(idjc)&
                                    +ERIMO(ibja)*ERIMO(iadbc)*ERIMO(idjc) )
                        tem2now=tem2now/(Eorb(ii)+Eorb(ij)-Eorb(ia)-Eorb(ib))
                        tem2now=tem2now/(Eorb(ii)+Eorb(ij)-Eorb(ic)-Eorb(id))
                        tem2=tem2+tem2now
                            end do 
                        end do 
                    end do 
                end do 
            end do 
        end do
!!!!  part 3      
        do ii=1,Nocc
            do ia=Nocc+1,nbasis
                do ij=1,Nocc
                    do ib=Nocc+1,nbasis
                        do ic=Nocc+1,nbasis
                            do ik=1,Nocc
                        iajb=indexeri(ii,ia,ij,ib)
                        ibja=indexeri(ii,ib,ij,ia)
                        jbck=indexeri(ij,ib,ic,ik)
                        jkcb=indexeri(ij,ik,ic,ib)
                        iakc=indexeri(ii,ia,ik,ic)
                        icka=indexeri(ii,ic,ik,ia)
                        tem3now=8d0*(ERIMO(iajb)*ERIMO(jbck)*ERIMO(iakc))
                        tem3now=tem3now-4d0*&
                                    (ERIMO(ibja)*ERIMO(jbck)*ERIMO(iakc)&
                                    +ERIMO(iajb)*ERIMO(jkcb)*ERIMO(iakc)&
                                    +ERIMO(iajb)*ERIMO(jbck)*ERIMO(icka)&
                                    +ERIMO(ibja)*ERIMO(jkcb)*ERIMO(icka) )
                        tem3now=tem3now+2d0*&
                                    (ERIMO(ibja)*ERIMO(jkcb)*ERIMO(iakc)&
                                    +ERIMO(ibja)*ERIMO(jbck)*ERIMO(icka)&
                                    + ERIMO(iajb)*ERIMO(jkcb)*ERIMO(icka))        
                        tem3now=tem3now/(Eorb(ii)+Eorb(ij)-Eorb(ia)-Eorb(ib))
                        tem3now=tem3now/(Eorb(ii)+Eorb(ik)-Eorb(ic)-Eorb(ia))
                        tem3=tem3+tem3now
                            end do 
                        end do 
                    end do 
                end do 
            end do 
        end do   
        Emp3=0.125d0*tem1+0.125d0*tem2+tem3
!!! mp2.5         
        if(mmp==3)then
            write(kfid,*)'=======================================&
                  ========================================='
            write(kfid,'(25X,A30)')'*****  MP2.5 Summary *****'        
            Emp25=Emp2+0.5d0*Emp3
            write(kfid,'(25X,A11,5X,F20.15)')'E(MP2.5):',Emp25
            write(kfid,'(21X,A15,5X,F20.15)')'E(HF+MP2.5):',Eelec+Emp25    
            emp3=0.5d0*Emp3 !!! for final energy     
            return   
        end if               
!!!!!        
        write(kfid,*)'=======================================&
                  ========================================='
        write(kfid,'(25X,A30)')'*****  MP3 Summary *****'
        write(kfid,*) 'MP3 Final results (a.u.)'
        write(kfid,'(25X,A11,5X,F20.15)')'E(MP2):',Emp2
        write(kfid,'(25X,A11,5X,F20.15)')'E(MP3):',Emp3
        write(kfid,'(21X,A15,5X,F20.15)')'E(HF+MP2+MP3):',Emp3+Eelec+Emp2
        write(*,*)'*** The MP3 procedure Done ***'      
    end subroutine  MP3






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                           Transint.f90                             !!  
!!  Transform intgrals from AO to MO,from spatial orbit to spin orbit.!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
    subroutine AO2MOERI(nbasis,neri,ERIall,Cof,ERIMO)
        use Global,only:indexeri
        implicit double precision (A-H,O-Z)
        Integer*16 neri
        dimension::ERIall(neri),Cof(nbasis,nbasis),ERIMO(neri)
        kfid=99
!!!!!!!!  Direct method. N^8      
        write(*,*)'*** Start AO to MO ***' 
        write(kfid,*)'=======================================&
                  ========================================='        
        write(kfid,'(20X,A30)') '***  AO to MO  *** ' 
        n=0
        ERIMO=0d0
        do i=1,nbasis
            do j=1,i
                do k=1,i
                    lend=k
                    if (i==k)lend=j
                    do l=1,lend                       
                        ! if(i*(i+1)/2+j<k*(k+1)/2+l .or. i<j .or. k<l)cycle
                        ijkl=indexeri(i,j,k,l) 
                        n=n+1
                        if(mod(n,5000)==0)write(*,*)n,'ERI transfer done,',real(n)/neri*100d0,'% done'                        
                        do mi=1,nbasis
                            do mj=1,nbasis
                                do mk=1,nbasis
                                    do ml=1,nbasis
                                    ! if(mi*(mi+1)/2+mj<mk*(mk+1)/2+ml .or. mi<mj .or. mk<ml)cycle
                                    mijkl=indexeri(mi,mj,mk,ml)
                                    ERIMO(ijkl)=ERIMO(ijkl)+ERIall(mijkl)*cof(mi,i)*cof(mj,j)*cof(mk,k)*cof(ml,l)                                       
                                    end do !!! ml
                                end do !!! mk
                            end do !!!mj
                        end do !!! mi
                    end do !!! l
                end do !!! k
            end do !!!j
        end do !!! i
        write(*,*)n
        write(*,*)'*** Start AO to MO Done ***' 
    end subroutine
!!!!===============================================================!!!!!!!!!!
    subroutine AO2MOERIfast(nbasis,neri,hcore,ERIall,Cof,HcoreMO,ERIMO)
        use Global,only:indexeri,index2eri
        implicit double precision (A-H,O-Z)
        Integer*16 neri
        dimension::ERIall(neri),Cof(nbasis,nbasis),ERIMO(neri),&
                    tem(nbasis*(nbasis+1)/2,nbasis*(nbasis+1)/2),&
                    temx(nbasis,nbasis), temy(nbasis,nbasis),&
                    hcore(nbasis,nbasis),hcoreMO(nbasis,nbasis)
        kfid=99
!!!!!!!!  The fast method. N^5   
        write(*,*)'*** Start AO to MO fast ***' 
        write(kfid,*)'=======================================&
                  ========================================='        
        write(kfid,'(20X,A30)') '***  AO to MO  *** '        
        ERIMO=0d0
        temx=0d0
        temy=0d0
        tem=0d0
        do i1=1,nbasis
            do j1=1,i1
                do k1=1,nbasis
                    do l1=1,k1
                        ijkl=indexeri(i1,j1,k1,l1)
                        temx(k1,l1)=ERIall(ijkl)
                        temx(l1,k1)=ERIall(ijkl)
                    end do 
                end do 
                temy=0d0
                temy=matmul(transpose(Cof),temx)
                temx=0d0
                temx=matmul(temy,Cof)
                do k2=1,nbasis
                    do l2=1,k2
                        tem(index2eri(k2,l2),index2eri(i1,j1))=temx(k2,l2)
                    end do 
                end do                 
            end do 
        end do 
!!!!
        do k3=1,nbasis
            do l3=1,k3
                temx=0d0
                temy=0d0
                do i2=1,nbasis
                    do j2=1,i2
                        temx(i2,j2)=tem(index2eri(k3,l3),index2eri(i2,j2))
                        temx(j2,i2)=temx(i2,j2)
                    end do 
                end do 
                temy=0d0
                temy=matmul(transpose(Cof),temx)
                temx=0d0
                temx=matmul(temy,Cof)     
                do i3=1,nbasis
                    do j3=1,i3
                        klij=indexeri(k3,l3,i3,j3)
                        ERIMO(klij)=temx(i3,j3)
                    end do 
                end do  
            end do 
        end do 
!!!!   transform AO Hcore to MO Hcore
        hcoreMO=matmul(hcore,cof)
        hcoreMO=matmul(transpose(Cof),hcoremo)        
        write(*,*)'*** AO to MO Done ***' 
    end subroutine    
!!!!===============================================================!!!!!!!!!!
!!!!===============================================================!!!!!!!!!!
    subroutine MO2SOERI(nbasis,neri,ERIMO,ERISPIN)
        use Global,only:indexeri,index2eri,krodelta
        implicit double precision (A-H,O-Z)
        Integer*16 neri
        dimension::ERIMO(neri),ERISPIN(2*nbasis,2*nbasis,2*nbasis,2*nbasis)
        kfid=99 !!! OUTPUT file number
!!!!!!!!   (ij||kl)  -> (pq||rs)
!!!!!!!!    i,j,k,l is MO ,  p,q,r,s is  spin-orbital.
        write(*,*)'*** Start MO to spin-orbital ERI ***' 
        write(kfid,*)'=======================================&
                  ========================================='        
        write(kfid,'(20X,A30)') '***  MO to spin-orbital  *** '        
        ERISPIN=0d0
        nspin=2*nbasis 
        do ip=1,nspin
            ipspa=ceiling(ip/2d0)
            modip2=mod(ip,2)
            do iq=1,nspin  
            iqspa=ceiling(iq/2d0)
            modiq2=mod(iq,2) 
            kdeltpq= krodelta(modip2,modiq2)
                do ir=1,nspin    
                    irspa=ceiling(ir/2d0)
                    modir2=mod(ir,2)
                    do is=1,nspin
                        isspa=ceiling(is/2d0)
                        modis2=mod(is,2)
                        kdeltrs= krodelta(modis2,modir2)
                        kdeltps= krodelta(modis2,modip2)
                        kdeltrq= krodelta(modir2,modiq2)
!!!!                    
                kpqrs=indexeri(ipspa,iqspa,irspa,isspa)
                v1=ERIMO(kpqrs)*kdeltpq*kdeltrs
                kpsrq=indexeri(ipspa,isspa,irspa,iqspa)
                v2=ERIMO(kpsrq)*kdeltps*kdeltrq
                ERISPIN(ip,iq,ir,is)=v1-v2
                    end do
                end do 
            end do 
        end do 
        write(*,*)'***  MO to spin-orbital Done ***' 
    end subroutine    
!!!!===============================================================!!!!!!!!!!
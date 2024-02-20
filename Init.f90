!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                            Init.f90                            !!!!
!!!!              This file Initialize program variables.           !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Module Init
        use Global
        implicit double precision (A-H,O-Z)      
!!!! INPUT (fid=33),OUTPUT (fid=99)
!!!! basis set file (fid=55),aux basis set file (fid=56)
!!!! 2eint.bin (fid=70)
        parameter(kfid=99,kinputfid=33,kbasfid=55,kauxbasfid=56,kerifid=70) 
!!!! for INPUT        
        character*80 inputarg,outputarg,inputfpath,outputfpath        
        character*80,save :: title
!!! System and Block TASK 
        save::ksystem,ktask,mhf,mmp,mcc
!!! Block BASIS   
        character*80,save :: basfpath,basisname
!!! Block MOLECULE          
        save::natom,kbastype,icharge,ispinm   
        allocatable:: coor(:,:),KATOM(:)    
!!! Block INT                         
        save::methods,methodT,methodv,methodERI   
        logical,save ::usecint
!!! Block SCF 
        save::damp,maxcycle,conver_E,conver_rms,Mguess,mdiis,ndiis      
!!! Block OUTPUT       
        save::Moutxyzfile,Moutlevel     
!!! Block POSTHF
        save::ccdamp,maxcyclecc,conver_Ecc,conver_rmscc,mdiiscc,ndiiscc         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                 
!!!! for time                           
        character(len=8),save ::theDate,theDate1 !!! for date_and_time
        character(len=10),save::t1_str,t2_str 
        save:: t_start,t_end,it1,it2,irate 
!!!!  For general situations    
        Integer*16 neri,nbasisbig                      
        save::Enuc,nele,nalpha,nbeta      
        logical,save:: iscloseshell
        save::nshell,nbasis,ngfall,norb   
        integer,save::lgfttypex(84),lgfttypey(84),lgfttypez(84)              
        allocatable::Kcenter(:),Ntype(:),ngf(:),gfe(:),gfc(:),&
                     kbasst(:),kbased(:),knbas(:),lgfst(:),lgfed(:),&
                     kpgfst(:),kpgfed(:),naxyz(:)
        allocatable::S(:,:),T(:,:),V(:,:),dipx(:,:),dipy(:,:),dipz(:,:),&
                     ERIall(:),Cof(:,:),Hcore(:,:),eorb(:),X(:,:),ERIMO(:),&
                     f12ERI(:),f12S(:),f12T(:,:,:,:),f12NAI(:,:,:,:),&
                     ERISPIN(:,:,:,:),HcoreMO(:,:)  
!!!!
        contains
!!!!===================================================================!!!! 
        subroutine Setvalue
            implicit double precision (A-H,O-Z)  
            iscloseshell= .true.
!!!! Method
            mhf=1
            mmp=0
            mcc=0
!!!! Block OUTPUT            
            Moutlevel=2  !!!! 1-5 ,less to more
            Moutxyzfile=0 
            lprint1e=0
            lprint2e=0
!!!! Block INT                       
            methods=1
            methodT=1
            methodv=2
            methodERI=4
            usecint=.false.
!!!! Block SCF
!!! The SCF convergence threshold value =10^-(conver)
            Mguess=1 
            conver_E=12
            conver_rms=12   
            maxcycle=200
            damp=0d0   
            mdiis=1
            ndiis=7            
!!!! POSTHF
            ccdamp=0d0
            maxcyclecc=200
            conver_Ecc=12
            conver_rmscc=12
            mdiiscc=1
            ndiiscc=7                      
        end subroutine   Setvalue 
!!!!===================================================================!!!! 
    subroutine Setbasicinf
        implicit double precision (A-H,O-Z)  
!!!!
        ngfall=size(gfc)
        Nele=sum(katom)-icharge
        ngfsum=sum(ngf(:))  ! total number of primitive Gauss functions(pgf)
        nshell=size(ntype)
        allocate(kbasst(nshell),kbased(nshell),knbas(nshell),&
                 lgfst(nshell),lgfed(nshell),kpgfst(nshell),kpgfed(nshell))
        nbasis=0
        lst=0
        lend=0
!====== set data for basis function  (shell->function)
        ntem=0  !!!! pgf number of the type pgf
        do ii=1,nshell
            select case (Ntype(ii))
                case(1)     ! S
                    ntem=1
                    lst=1
                    lend=1
                case(2)     ! P
                    ntem=3
                    lst=2
                    lend=4                       
                case(3)     ! D
                    ntem=6
                    lst=5
                    lend=10    
                case(4)     ! F
                    ntem=10    
                    lst=11
                    lend=20                                
                case(5)     ! G
                    ntem=15
                    lst=21
                    lend=35  
                case(6)     ! H
                    ntem=21 
                    lst=22
                    lend=42   
                case(7)     ! I
                    ntem=28 
                    lst=43
                    lend=70                            
                case(12)     ! SP
                    ntem=4 
                    lst=1
                    lend=4                                                 
            end select    
            nbasis=nbasis+ntem  
            kbased(ii)=nbasis  
            knbas(ii)=ntem
            lgfed(ii)=lend
            lgfst(ii)=lst
            kpgfed(ii)=sum(Ngf(1:ii))
            kpgfst(ii)= kpgfed(ii)-Ngf(ii)+1
        end do
        kbasst=kbased-knbas+1
        nss=count(Ntype==1)
        nps=count(Ntype==2)
        nds=count(Ntype==3)
        nfs=count(Ntype==4)
        ngs=count(Ntype==5)
        nhs=count(Ntype==6)
        nis=count(Ntype==7)
        nsp=count(Ntype==12)
        norb=nbasis-nds-3*nfs-6*ngs-10*nhs
        allocate(naxyz(3*nbasis))
!!!! prepare basis type
        if(usecint)then
            lgfttypex=lgfttypexPySCF
            lgfttypey=lgfttypeyPySCF
            lgfttypez=lgfttypezPySCF
        else
            lgfttypex=lgfttypexGMS
            lgfttypey=lgfttypeyGMS
            lgfttypez=lgfttypezGMS
        end if
        do Ish=1,nshell   
            lgft=0
            do kba=kbasst(ish),kbased(ish)
                    lgf=lgfst(ish)+lgft
                    naxyz(kba)=lgfttypex(lgf)
                    naxyz(kba+nbasis)=lgfttypey(lgf)
                    naxyz(kba+2*nbasis)=lgfttypez(lgf)
                    lgft=lgft+1  
            end do
        end do   
        if(nbasis<215)then
            neri=(nbasis**4+2*nbasis**3+3*nbasis**2+2*nbasis)/8   
        else     
            nbasisbig=nbasis  
            neri=(nbasisbig**4+2*nbasisbig**3+3*nbasisbig**2+2*nbasisbig)/8  
            write(*,*)nbasisbig,neri,neri*8/1024**3
        end if   
        allocate(S(nbasis,nbasis),T(nbasis,nbasis),&
            V(nbasis,nbasis),Hcore(nbasis,nbasis),ERIall(neri),&
            dipx(nbasis,nbasis),dipy(nbasis,nbasis),dipz(nbasis,nbasis),&
            Cof(nbasis,nbasis),eorb(nbasis),X(nbasis,nbasis),ERIMO(neri),&
            f12ERI(neri),f12S(neri),f12T(nbasis,nbasis,nbasis,nbasis),&
            f12NAI(nbasis,nbasis,nbasis,nbasis),HcoreMO(nbasis,nbasis)) 
!!!!         
        end subroutine   Setbasicinf
!!!!===================================================================!!!! 
    end module










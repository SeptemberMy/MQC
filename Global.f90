!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                             Global.f90                         !!!!
!!!!                This file includes the general module.          !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!1. Global
!description: module
!   Storage public variables，parameters and functions
! function:
!1.Dfact(n)
!description: function
!   Calculate the value of n-th double factorial
!   output real(kind=8) value
!   input:   (1) n       integer
!   output:  (1) Dfact   real(kind=8)  
!2. iFact(n)
!description:  function
!   Calculate the value of n-th factorial
!   input:   (1) n       integer
!   output:  (1) Fact    integer(kind=8) 
!3. iComb(n,m)
!description: function
!   Calculate the Number of combinations
!   $C_n^m$   
!   input:   (1) n       integer
!            (2) m       integer
!   output:  (1) iComb    integer(kind=8)
!4. F(m,w),F_ser(m,w),F_pade(m,w)
!description: function
!   Incomplete gamma function (3 methods)
!   $F_m(w)=\int_0^1 e^{-wt^2}t^{2m}dt$   
!   input:   (1) m       integer
!            (2) w       real(kind=8)
!   output:  (1) F       real(kind=8)
!5. fi(num,na,nb,PA,PB)
!description: function
! Binomial expansion      
!   $f_i(na,nb,(P_x-A_x),(P_x-B_x))=\sum_{\lambda=0}^{na}
!     \sum_{\mu=0(\lambda+\mu=i)}^{nb}C_{na}^{\lambda}C_{nb}^{\mu}
!     (P_x-A_x)^{na-\lambda}(P_x-B_x)^{na-\mu}$
!   input:   (1)num      integer
!               nax+nbx=num
!            (2)na     integer
!               nax
!            (3)nb    integer
!               nbx
!            (4)PA    double precision
!                P_i ,i={x,y,z}
!            (5) PB    double precision
!   output:  (1) fi    double precision 
!6. Indexeri(i,j,k,l)
!description: function
!  Compute the index of two electron intgrals 
! for fortran type ,first element index is 1,
! so the code is a little different form formula.
!   $ij=i(i+1)/2+j$(i>=j)   
!   $ij=j(j+1)/2+i$(i<=j)  
!   $ijkl=ij(ij+1)/2+kl$(i<=j)    
!   input:   (1) i       integer
!            (2) j       integer
!            (3) k       integer
!            (4) l       integer
!   output:  (1) index       integer
!7.Normgf(ea,nax,nay,naz)
!description: function
!  Compute the normalization coefficient
!  of gauss function 
!   $Na=(2*ea/pi)**0.75*
!      [(4*ea)**(nax+nay+naz)/((2nax-1)!!*(2nay-1)!!)*(2nay-1)!!]**0.5$      
!   input:   (1) ea        double precision
!            (2) nax       integer
!            (3) nay       integer
!            (4) naz       integer
!   output:  (1) normgf    double precision
!8.outputmatrixd(ndimrow,ndimcol,dmar,id),outputmatrixi(ndimrow,ndimcol,dmar,id)
!description: subroutine
! Output the matrix (double precision,integer.See end of name) by row.
!   input:   (1) ndimrow  ndimcol   integer
!                 dimension of matrix
!            (2) dmar/mar  double precision/ integer
!                    the matrix
!            (3) id       integer
!                  file id (OUTPUT:kfid,screen:6)
!9.dexchange(a,b),iexchange(a,b)
!description: subroutine
! Exchange a and b.(integer or double precision)
!   input:   (1) a/ia      double precision/integer
!            (2) b/ib  double precision/integer
!10. index2eri(i,j)
!description: function
! if i>j,ij=i*(i-1)/2+j, if i<=j,ij=j*(j-1)/2+i
!   input:   (1) i      integer
!            (2) j  integer
!11. matrix1D22D(ndim,one,two)
!description: subroutine 
! Restore the symmetric matrix stored in a one-dimensional array 
! to a two-dimensional symmetric matrix.
!   input:   (1) ndim      integer
!            (2) one(ndim*(ndim+1)/2)  1D  array
!   output:  (1) two(ndim,ndim)
!12. krodelta(i,j)
!description: function
! if i=j,result=1, else result=0.
!   input:   (1) i  integer
!            (2) j  integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!! Module Global !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Basic data and functions
      Module Global
            implicit double precision (A-H,O-Z)
            parameter  (pi =acos(-1d0))
            parameter (icount_max=2147483647)
!====================================================!
!!!!! unit conversion            
!!! 1 bohr=0.5291772083d0 Angstroms            
            parameter (BOHR2A = 0.5291772083d0)
!!! 1 Eh(hartree) = 27.2113834 eV               
            parameter (Eh2eV = 27.2113834d0)
!!! 1 eV = 8065.54477 cm−1 = 23.0605 kcal/mol
            parameter (eV2kcalpm = 23.0605d0)
!!! 1 kcal/mol= 4.18585182085 kJ/mol
            parameter (kcal2kJ =4.18585182085d0)      
!!! 1 a.u.= 2.542 debye = 8.478×10-30 C·m
            parameter (au2debye =2.542d0)                     
!========================================================================!
!!!!  D 
!  the order of spherical d functions in Gaussian           
!  the order of spherical d functions in PySCF
!  the order of spherical d functions in GAMESS           
! 1    2    3    4    5
! d0 , d+1, d-1, d+2, d-2 (Gaussian)
! d-2, d-1, d0 , d+1, d+2 (PySCF) 
! d0 , d+1, d-1, d+2, d-2 (GAMESS)     
!  the order of Cartesian d functions in Gaussian
!  the order of Cartesian d functions in PySCF
!  the order of Cartesian d functions in GAMESS
! 1  2  3  4  5  6
! XX,YY,ZZ,XY,XZ,YZ  (Gaussian)
! XX,XY,XZ,YY,YZ,ZZ  (PySCF)      
! XX,YY,ZZ,XY,XZ,YZ  (GAMESS)
!!!!  F  
!  the order of spherical f functions in Gaussian
!  the order of spherical f functions in PySCF
! 1    2    3    4    5    6    7
! f0 , f+1, f-1, f+2, f-2, f+3, f-3 (Gaussian)
! f-3, f-2, f-1, f0 , f+1, f+2, f+3 (PySCF)
!  the order of Cartesian f functions in Gaussian
!  the order of Cartesian f functions in PySCF
! 1   2   3   4   5   6   7   8   9   10
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ (Gaussian)
! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ  (PySCF)  
! XXX,YYY,ZZZ,XXY,XXZ,XYY,YYZ,XZZ,YZZ,XYZ (GAMESS)                               
!!!!  G  
!  the order of spherical g functions in Gaussian
!  the order of spherical g functions in PySCF
! 1    2    3    4    5    6    7    8    9
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4 (Gaussian)
! g-4, g-3, g-2, g-1, g0 , g+1, g+2, g+3, g+4 (PySCF)
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX (Gaussian)
! xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz  (PySCF)  
! XXXX,YYYY,ZZZZ,XXXY,XXXZ,XYYY,YYYZ,XZZZ,YZZZ,XXYY,XXZZ,YYZZ,XXYZ,XYYZ,XYZZ (GAMESS)                    
!!!!  H              
!  the order of spherical h functions in Gaussian
!  the order of spherical h functions in PySCF
! 1    2    3    4    5    6    7    8    9    10   11
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5 (Gaussian)
! h-5, h-4, h-3, h-2, h-1, h0 , h+1, h+2, h+3, h+4, h+5 (PySCF)
!  the order of Cartesian h functions in Gaussian
!  the order of Cartesian h functions in PySCF
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX (Gaussian)
! xxxxx,xxxxy,xxxxz,xxxyy,xxxyz,xxxzz,xxyyy,xxyyz,xxyzz,xxzzz,xyyyy,xyyyz,xyyzz,xyzzz,xzzzz,yyyyy,yyyyz,yyyzz,yyzzz,yzzzz,zzzzz  (PySCF)  
! XXXXX,YYYYY,ZZZZZ,XXXXY,XXXXZ,XYYYY,YYYYZ,XZZZZ,YZZZZ,XXXYY,XXXZZ,XXYYY,YYYZZ,XXZZZ,YYZZZ,XXXYZ,XYYYZ,XYZZZ,XXYYZ,XXYZZ,XYYZZ (GAMESS) 
!!!   PySCF order.                                
      integer::lgfttypexPySCF(84)=(/ &
               0, 1, 0, 0, 2, 1, 1, 0, 0, 0,&                    ! S,P,D
               3, 2, 2, 1, 1, 1, 0, 0, 0, 0,&                    ! F 
               4, 3, 3, 2, 2, 2, 1, 1, 1, 1,&                    
               0, 0, 0, 0, 0,&                                   ! G
               5, 0, 0, 4, 4, 1, 0, 1, 0, 3,&
               3, 2, 0, 2, 0, 3, 1, 1, 2, 2,&
               1,&                                               ! H
               6, 0, 0, 5, 5, 1, 0, 1, 0, 4,&
               4, 2, 0, 2, 0, 4, 1, 1, 3, 3,&
               0, 3, 3, 2, 1, 2, 1, 2/)                          ! I
      integer::lgfttypeyPySCF(84)=(/ &
               0, 0, 1, 0, 0, 1, 0, 2, 1, 0,&                    ! S,P,D 
               0, 1, 0, 2, 1, 0, 3, 2, 1, 0,&                    ! F  
               0, 1, 0, 2, 1, 0, 3, 2, 1, 0,&                     
               4, 3, 2, 1, 0,&                                   ! G     
               0, 5, 0, 1, 0, 4, 4, 0, 1, 2,&                     
               0, 3, 3, 0, 2, 1, 3, 1, 2, 1,&                     
               2,&                                               ! H     
               0, 6, 0, 1, 0, 5, 5, 0, 1, 2,&                     
               0, 4, 4, 0, 2, 1, 4, 1, 3, 0,&                     
               3, 2, 1, 3, 3, 1, 2, 2/)                          ! I            
      integer::lgfttypezPySCF(84)=(/ &
               0, 0, 0, 1, 0, 0, 1, 0, 1, 2,&                    ! S,P,D 
               0, 0, 1, 0, 1, 2, 0, 1, 2, 3,&                    ! F
               0, 0, 1, 0, 1, 2, 0, 1, 2, 3,&
               0, 1, 2, 3, 4,&                                   ! G
               0, 0, 5, 0, 1, 0, 1, 4, 4, 0,&
               2, 0, 2, 3, 3, 1, 1, 3, 1, 2,&
               2,&                                               ! H
               0, 0, 6, 0, 1, 0, 1, 5, 5, 0,&
               2, 0, 2, 4, 4, 1, 1, 4, 0, 3,&
               3, 1, 2, 1, 2, 3, 3, 2/)                          ! I
!!!   GAMESS order.     
      integer::lgfttypexGMS(84)=(/&
               0, 1, 0, 0, 2, 0, 0, 1, 1, 0,&                    ! S,P,D
               3, 0, 0, 2, 2, 1, 0, 1, 0, 1,&                    ! F    
               4, 0, 0, 3, 3, 1, 0, 1, 0, 2,&
               2, 0, 2, 1, 1,               &                    ! G
               5, 0, 0, 4, 4, 1, 0, 1, 0, 3,&
               3, 2, 0, 2, 0, 3, 1, 1, 2, 2,&
               1,                           &                    ! H
               6, 0, 0, 5, 5, 1, 0, 1, 0, 4,&
               4, 2, 0, 2, 0, 4, 1, 1, 3, 3,&
               0, 3, 3, 2, 1, 2, 1, 2/)                          ! I
      integer::lgfttypeyGMS(84)=(/ &
               0, 0, 1, 0, 0, 2, 0, 1, 0, 1,&                    ! S,P,D
               0, 3, 0, 1, 0, 2, 2, 0, 1, 1,&                    ! F    
               0, 4, 0, 1, 0, 3, 3, 0, 1, 2,&
               0, 2, 1, 2, 1,               &                    ! G
               0, 5, 0, 1, 0, 4, 4, 0, 1, 2,&
               0, 3, 3, 0, 2, 1, 3, 1, 2, 1,&
               2,                           &                    ! H
               0, 6, 0, 1, 0, 5, 5, 0, 1, 2,&
               0, 4, 4, 0, 2, 1, 4, 1, 3, 0,&
               3, 2, 1, 3, 3, 1, 2, 2/)                          ! I            
      integer::lgfttypezGMS(84)=(/&
               0, 0, 0, 1, 0, 0, 2, 0, 1, 1,&                    ! S,P,D      
               0, 0, 3, 0, 1, 0, 1, 2, 2, 1,&                    ! F           
               0, 0, 4, 0, 1, 0, 1, 3, 3, 0,&                     
               2, 2, 1, 1, 2,               &                    ! G              
               0, 0, 5, 0, 1, 0, 1, 4, 4, 0,&                           
               2, 0, 2, 3, 3, 1, 1, 3, 1, 2,&                           
               2,                           &                    ! H                  
               0, 0, 6, 0, 1, 0, 1, 5, 5, 0,&                           
               2, 0, 2, 4, 4, 1, 1, 4, 0, 3,&                           
               3, 1, 2, 1, 2, 3, 3, 2/)                          ! I                
!!!! atom name                                                                                   
      character*2::atomsym(118)=(/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',&
                                 'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',&
                                 'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',&
                                 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',&
                                 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',&
                                 'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',&
                                 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',&
                                 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',&
                                 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',&
                                 'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',&
                                 'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds',&
                                 'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/)
!!!!                                 
      character*2::typesym(12)=(/'S ','P ','D ','F ','G ','H ','I ','J ','K ','L ','M ','SP'/)
!!!!! Transformation matrix  (Gaussian type)     
!!!   Cartesian d functions (6D)   <->   Spherical d functions (5D)                                          
!!! dimension (6,5)

      double precision ::D62D5gau(6,5)=reshape(&
            (/-0.5d0,-0.5d0,1d0,0d0,0d0,0d0,&               ! column 1
            0d0,0d0,0d0,0d0,1d0,0d0,&                       ! column 2
            0d0,0d0,0d0,0d0,0d0,1d0,&                       ! column 3
            0.5*sqrt(3d0),-0.5*sqrt(3d0),0d0,0d0,0d0,0d0,&  ! column 4
            0d0,0d0,0d0,1d0,0d0,0d0&                        ! column 5
            /),(/6,5/))            
!!!   Cartesian f functions (10F)  <->   Spherical f functions (7F)
      double precision ::F102F7gau(10,7)=reshape(&
            (/&
            0d0,0d0,1d0,0d0,0d0,-3d0/2/sqrt(5d0),0d0,0d0,-3d0/2/sqrt(5d0),0d0,&        ! column 1
            -sqrt(3/8d0),0d0,0d0,-sqrt(3/40d0),0d0,0d0,sqrt(6/5d0),0d0,0d0,0d0,&       ! column 2
            0d0,-sqrt(3/8d0),0d0,0d0,-sqrt(3/40d0),0d0,0d0,sqrt(6/5d0),0d0,0d0,&       ! column 3
            0d0,0d0,0d0,0d0,0d0,sqrt(3d0)/2d0,0d0,0d0,-sqrt(3d0)/2d0,0d0,&             ! column 4
            0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,1d0,&                                  ! column 5
            sqrt(5/8d0),0d0,0d0,-3/sqrt(8d0),0d0,0d0,0d0,0d0,0d0,0d0,&                 ! column 6
            0d0,-sqrt(5/8d0),0d0,0d0,3/sqrt(8d0),0d0,0d0,0d0,0d0,0d0&                 ! column 7
            /),(/10,7/))                                       
!!!   Cartesian g functions (15G)  <->   Spherical g functions (9G)
!!!   Cartesian h functions (21H)  <->   Spherical h functions (11H)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains
!======================================================================!
!!! iFact(n)
      function iFact(n)
      implicit double precision (A-H,O-Z)
      integer(kind=16) ::iFact,item
      item=1
      if (n==0 .or. n==1) then
            ifact=1
            return
      else if(n<0)  then
            write(*,*) 'bad input of iFact(n)'
            return
      else
            do i = 2, n
                  item = item*i
            end do
            ifact = item
      end if
      end function iFact
!======================================================================!
!!! diFact(n)
!!!! real*8,nmax=170
!!!! real*16,nmax=1754
      function diFact(n)
      implicit double precision (A-H,O-Z)
      real(kind=16)::difact,ditem
      ditem=1d0
      if (n==0 .or. n==1) then
            difact=1d0
            return
      else if(n<0)  then
            write(*,*) 'bad input of iFact(n)'
            return
      else
            do i = 2, n
                  ditem = ditem*i
            end do
            difact = ditem
      end if
      end function diFact
!======================================================================!
!======================================================================!
!!!! Dfact(n)： double factorial
      function Dfact(n)
      implicit double precision (A-H,O-Z)
            dTemp = 1.d0
            Dfact=1.d0
            if (n<0)  then
                  if( mod(abs(n),2) == 0 ) then
                        write(*,*) 'bad input of Dfact(n)'
                        return
                  end if
                  if(abs(n) == 1 ) then
                        Dfact=1d0
                        return
                  end if
                  nt=(abs(n)-1)/2
                  do i=1,nt
                        dTemp=dTemp*(2*real(i,8)-1)**(-1)
                  end do
                  Dfact=(-1)**n*dTemp
                  return
            else if (n==0 .or. n==1) then
                  Dfact=1
                  return    
            else if (mod(n,2)==0) then
                  do i = 1, n/2
                        dTemp = dTemp*2.d0*real(i,8)
                  end do
            else
                  do i = 1, (n+1)/2
                        dTemp = dTemp*(2.d0*real(i,8)-1.d0)
                  end do
            end if
            Dfact = dTemp
      
      end function Dfact
!======================================================================!
!======================================================================!
!!!!! iComb(n,m)
      function iComb(n,m)
      implicit double precision (A-H,O-Z)
!!!
      if(n<m) then
            write(*,*) 'bad input of  iComb(n,m)'
            return
      else 
            item = 1
            do i = n - m + 1, n
                  item =item*i
            end do
            iComb = item/ifact(m)
      end if

      end function iComb
!======================================================================!
!!!!!!!!! Boys function (Incomplete gamma function） 
      function F(m,w)
!!!!   $F_m(w)=\int_0^1 e^{-wt^2}t^{2m}dt$   
      implicit double precision (A-H,O-Z)
      double precision  F_T(0:m)
!!!! Analytical method
      if (w<0d0 .or. m<0) then
            write(*,*) 'bad input of  F(m,w)'
            return  
      else 	if( w < 1.0d-16 .or. w==0 ) then
            F = 1.0d0/(2*m+1)
            return
      else if(w<m+1.5d0)then
            F = F_ser(m,w)
            return
      end if
!!    
      t = sqrt(w)
      ! F = 0.88622692545275801d0*derf(t)/t
      F =  sqrt(Pi)*0.5d0*derf(t)/t
      if( m == 0 ) then
            return
      end if
      F_T(0) = F
      do i = 1,m
            F_T(i) = 0.5d0/w*( (2*i-1)*F_T(i-1) - exp(-w) )
      end do
      F = F_T(m)     
!!!!!!!!!!!!!!!       
!!!!!!! 
       return    
      end function F
!======================================================================!
      function F_ser(m,w) result(F)
!!!!   $F_m(w)=\int_0^1 e^{-wt^2}t^{2m}dt$   
      implicit double precision (A-H,O-Z)
      double precision  F_T(0:m)
!!!!!!!!!!!!!!!     
      !!!!! Journal of Mathematical Chemistry Vol. 36, No. 3, July 2004     
            ! if(w<=35d0) then
            !       F_T(0)=sqrt(Pi)/2*erf(t)/t
            !       i=0
            !       F=0d0
            !       do i=1,100
            !             do j=0,2*i,2
            !                   ntem=ntem*(2*m+1+2*j)
            !             end do
            !             tem=exp(-w)*(2*w)**i/ntem
            !            F =F+tem
            !            if (tem<1d-15) return
            !       !      i=i+1
            !       end do
            !       F = F_T(m)
            ! else if (w<=60d0) then
            !       F_T(0)=sqrt(pi/w)*0.5d0
            !       do i = 1,m
            !             F_T(i) =(2*m-1)*F_T(i-1)/w
            !       end do
            !       F = F_T(m)        
            ! else 
            !       F_T(m)=0.5d0*sqrt(pi/m)*Dfact(2*m-1)/(2*w)**m
            !       F = F_T(m)        
            ! end if
!!!!      series     
      F=0d0
      tem=0d0
      i=0
      do while(.true.)
            ! dte=(-w)**i/ifact(i)/(2*m+2*i+1)
            dte=Dfact(2*m-1)*(2*w)**i/Dfact(2*m+2*i+1)     
            tem=tem+dte      
            if(abs(dte)<1d-20) exit
            i=i+1
      end do 
      F=exp(-w)*tem   
!!!!!!! 
       return    
      end function F_ser
!======================================================================!
      function F_pade(m,w)
!!!! $F_m(w)=\int_0^1 e^{-wt^2}t^{2m}dt$   
      implicit double precision (A-H,O-Z)
      allocatable ::aa(:),bb(:)
      if (w<0d0 .or. m<0) then
            write(*,*) 'bad input of  F_pade(m,w)'
            return  
      end if     
      aa0=(2*m+1)**(-2d0/(2*m+1)) 
      select case(m)
            case(0)
                  wmax=16.3578d0
                  na=5
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/0.213271302431420d0,0.629344460255614d-01,&
                      0.769838037756759d-02,0.758433197127160d-03,&
                  0.564691197633667d-04 /)
                  bb=(/0.879937801660182d0,0.338450368470103d0,&
                  0.738522953299624d-01,0.101431553402629d-01,&
                  0.955528842975585d-03,0.720266520392572d-04/)
            case(1)
                  wmax=17.4646d0
                  na=5
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/0.295195994716045d-01,0.128790985465415d-01,&
                      0.998165499553218d-03,0.970927983276419d-04,&
                  0.493839847029699d-05 /)
                  bb=(/0.461403194579124d0,0.108494164372449d0,&
                  0.171462934845042d-01,0.196918657845508d-02,&
                  0.160138863265254d-03,0.857708713007233d-05/)     
            case(2)
                  wmax=15.2368d0
                  na=4
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.575763488635418d-02,0.731474973333076d-02,&
                      0.251276149443393d-03,0.264336244559094d-04&
                      /)
                  bb=(/0.274754154712841d0,0.425364830353043d-01,&
                  0.493902790955943d-02,0.437251500927601d-03,&
                  0.288914662393981d-04/)            
            case(3)
                  wmax=16.0419d0
                  na=4
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.290110430424666d-01,0.561884370781462d-02,&
                      0.301628267382713d-04,0.110671035361856d-04&
                      /)
                  bb=(/0.171637608242892d0,0.187571417256877d-01,&
                  0.178536829675118d-02,0.137360778130936d-03,&
                 0.791915206883054d-05/)            
            case(4)
                  wmax=16.8955d0
                  na=4
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.452693111179624d-01,0.490070062899003d-02,&
                      -0.561789719979307d-04,0.550814626951998d-05&
                      /)
                  bb=(/0.108051989937231d0,0.855924943430755d-02,&
                  0.724968571389473d-03,0.502338223156067d-04,&
                 0.249107837399141d-05/)     
            case(5)
                  wmax=17.7822d0
                  na=4
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.566143259316101d-01,0.455916894577203d-02,&
                      -0.894152721395639d-04,0.328096732308082d-05&
                      /)
                  bb=(/0.662932958471386d-01,0.383724443872493d-02,&
                  0.327167659811839d-03,0.210430437682548d-04,&
                 0.883562935089333d-06/)       
            case(6)
                  wmax=15.8077d0
                  na=3
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.503249167534352d-01,0.273135625430953d-02,&
                      -0.310733624819100d-04&
                      /)
                  bb=(/0.586609328033371d-01,0.194044691497128d-02,&
                  0.109442742502192d-03,0.613406236401726d-05&
                     /)   
            case(7)
                  wmax=16.5903d0
                  na=3
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.548201062615785d-01,0.253099908233175d-02,&
                      -0.333589469427863d-04&
                      /)
                  bb=(/0.389873128779298d-01,0.5698908060832729d-03,&
                  0.422187129333708d-04,0.286010059144633d-05&
                     /)   
            case(8)
                  wmax=17.3336d0
                  na=3
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.581618006078160d-01,0.238525529084601d-02,&
                      -0.329989020317093d-04&
                      /)
                  bb=(/0.240929282666615d-01,-0.202677647499956d-03,&
                  0.119820675974460d-04,0.145762086904409d-05&
                     /)     
            case(9)
                  wmax=15.6602d0
                  na=2
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.334843993901400d-01,0.846637494147059d-03&
                      /)
                  bb=(/0.495875606944471d-01,0.946642302340943d-03,&
                  0.108367772249790d-04&
                     /)        
            case(10)
                  wmax=16.5258d0
                  na=2
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.335292171805959d-01,0.749168957990503d-03&
                      /)
                  bb=(/0.421492932021683d-01,0.582840985360327d-03,&
                  0.237676790577455d-05&
                     /)     
            case(11)
                  wmax=17.5395d0
                  na=2
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.332669773770348d-01,0.668720489602687d-03&
                      /)
                  bb=(/0.363057685289467d-01,0.345646100984643d-03,&
                  -0.190872330373450d-05&
                     /)       
            case(12)
                  wmax=18.5783d0
                  na=2
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.326241966410798d-01,0.598705175467956d-03&
                      /)
                  bb=(/0.318680048277695d-01,0.202419662347765d-03,&
                  -0.362095173837973d-05&
                     /)  
            case(13)
                  wmax=19.6511d0
                  na=2
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.317754368014894d-01,0.537678595933584d-03&
                      /)
                  bb=(/0.284036027081815d-01,0.113673420662576d-03,&
                  -0.416076810552774d-05&
                     /)       
            case(14)
                  wmax=20.7839d0
                  na=2
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.308755854748829d-01,0.485046451960769d-03&
                      /)
                  bb=(/0.255694625434059d-01,0.542010192055080d-04,&
                  -0.424759498527876d-05&
                     /) 
            case(15)
                  wmax=21.9998d0
                  na=2
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.3001438066719997d-01,0.439983032427912d-03&
                      /)
                  bb=(/0.231478878674366d-01,0.105546581596674d-04,&
                  -0.418932957034726d-05&
                     /)    
            case(16)
                  wmax=20.9225d0
                  na=2
                  nb=na+1
                  allocate(aa(na),bb(nb))
                  aa=(/-0.288346417991609d-01,0.397161796318408d-03&
                      /)  
                  bb=(/0.215021933120724d-01,-0.128592457453950d-05,&
                  -0.362120651688135d-05&
                     /)                                  
      end select
      sa=aa0
      sb=1d0
 !!!!     
      if(w<=wmax)then
            do  i=1,na
                  sa=sa+aa(i)*w**i
            end do
            do j=1,nb
                  sb=sb+bb(j)*w**j
            end do
            Fptem=(sa/sb)**(m+0.5d0)
      else
            Fptem=(Dfact(2*m-1)/(2*w)**(m+0.5))*sqrt(pi/2d0)
      end if
      F_pade=Fptem
      deallocate(aa,bb)
      end function F_pade
!======================================================================!
!======================================================================!
      function fi(num,na,nb,PA,PB)
!!!!!! Binomial expansion      
!!!!   $f_i(na,nb,(P_x-A_x),(P_x-B_x))=\sum_{\lambda=0}^{na}
!!!!     \sum_{\mu=0(\lambda+\mu=i)}^{nb}C_{na}^{\lambda}C_{nb}^{\mu}
!!!!     (P_x-A_x)^{na-\lambda}(P_x-B_x)^{nb-\mu}$
        implicit double precision (A-H,O-Z)
        ftem=0d0
        do ii=0,na
            do jj=0,nb
                if(ii+jj==num) then 
                    ftem=ftem+icomb(na,ii)*icomb(nb,jj)*&
                        PA**(na-ii)*PB**(nb-jj)
                end if
            end do
        end do
        fi=ftem
      end function fi
!======================================================================!
!======================================================================!
      function indexeri(i,j,k,l)
!!!!  index of ERI. (ij|kl)      
            implicit double precision (A-H,O-Z)
            indexeri=0
            ii=i-1
            jj=j-1
            kk=k-1
            ll=l-1
            if(ii<0 .or. jj<0.or. kk<0.or. ll<0) then
                  write(*,*) 'bad input of  index'
                  return
            end if
            if(ii>=jj)then
                  ij=ii*(ii+1)/2+jj
            else
                  ij=jj*(jj+1)/2+ii
            end if
            if(kk>=ll)then
                  kl=kk*(kk+1)/2+ll
            else
                  kl=ll*(ll+1)/2+kk
            end if
            if(ij>=kl)then
                  ijkl=ij*(ij+1)/2+kl+1
            else
                  ijkl=kl*(kl+1)/2+ij+1
            end if
            indexeri=ijkl
            return
      end function indexeri
!======================================================================!
!======================================================================!
      function dnormgf(ea,nax,nay,naz)
!!!!  Normalization coefficient of primitive gauss function 
            implicit double precision (A-H,O-Z)
            if(nax<0 .or. nay<0 .or. naz<0)then
                  write(*,*) 'bad input of normgf'
                  return
            end if
            tem=(2*ea/pi)**0.75*((4d0*ea)**(nax+nay+naz)/(Dfact(2*nax-1)&
                                *Dfact(2*nay-1)*Dfact(2*naz-1)))**0.5d0
            dnormgf=tem
            return
      end function dnormgf
!======================================================================!
!======================================================================!
      subroutine outputmatrixd(ndimrow,ndimcol,dmar,id)
!!!! Format output matrix  (ndimrow*ndimcol)    double precision
            implicit double precision (A-H,O-Z)
            dimension:: dmar(ndimrow,ndimcol),na(ndimcol)
            na=(/(i,i=1,ndimcol)/)
            if(mod(ndimcol,5)==0)then
                  ncom=ndimcol/5
            else
                  ncom=ndimcol/5+1
            end if
            do j=1,ncom
                  if(5*j>ndimcol)then
                        write(id,'(13x,5(i3,17x))')na(5*j-4:ndimcol)
                  else 
                        write(id,'(13x,5(i3,17x))')na(5*j-4:5*j)
                  end if
                  do i=1,ndimrow
                        if(5*j>ndimcol)then
                              write(id,'(i3,5(e20.9))')i,dmar(i,5*j-4:ndimcol)
                        else 
                              write(id,'(i3,5(e20.9))')i,dmar(i,5*j-4:5*j)
                        end if
                  end do 
            end do 
            return
      end  subroutine outputmatrixd  
!!!!!!!!!!!!
      subroutine outputmatrixd_tri(ndim,dmar,id)
!!!! Format output matrix  (ndim*ndim)    double precision
            implicit double precision (A-H,O-Z)
            dimension:: dmar(ndim,ndim)
            ncom=ndim/5
            if(mod(ndim,5)==0)then
                  ncom=ndim/5
            else
                  ncom=ndim/5+1
            end if
            do  k=1,ncom
                k1=(k-1)*5+1
                k2=5*k
                if (k==ncom) k2=ndim
                write(id,1233)(knum,knum=k1,k2)
                do  i=k1,ndim
                    jend=i
                    if (i.gt.k2) jend=k2
                    write(id,1234)i,(dmar(i,j),j=k1,jend)
                end do
            end do
                return
1233  format(1x,7x,5(8x,i5,8x))
1234  format(1x,i6,5(1x,e20.9))
            return
      end  subroutine outputmatrixd_tri
!!!!!!!!!!!!
!!!!!!!!!!!!
      subroutine outputmatrixi(ndimrow,ndimcol,mar,id)
!!!! Format output matrix  (ndim * nidm)    integer
            implicit double precision (A-H,O-Z)
            dimension:: mar(ndimrow,ndimcol),na(ndimcol)
            ndim=ndimcol
            na=(/(i,i=1,ndim)/)
            if(mod(ndim,5)==0)then
                  ncom=ndim/5
            else
                  ncom=ndim/5+1
            end if
            do j=1,ncom
                  if(5*j>ndim)then
                        write(id,'(13x,5(i3,17x))')na(5*j-4:ndim)
                  else 
                        write(id,'(13x,5(i3,17x))')na(5*j-4:5*j)
                  end if
                  do i=1,ndimrow
                        if(5*j>ndim)then
                              write(id,'(i3,5(i20))')i,mar(i,5*j-4:ndim)
                        else 
                              write(id,'(i3,5(i20))')i,mar(i,5*j-4:5*j)
                        end if
                  end do 
            end do             
            return
      end  subroutine outputmatrixi
!======================================================================!
!!!!!!! Exchange value of a and b. 
!!!!!  for  double precision  
      subroutine dexchange(a,b)
            implicit double precision (A-H,O-Z)
            tem=a
            a=b
            b=tem
      end  subroutine dexchange
!!!!! for integer     
      subroutine iexchange(ia,ib)
            implicit double precision (A-H,O-Z)
            item=ia
            ia=ib
            ib=item
      end  subroutine iexchange
!======================================================================!
      function index2eri(i,j)
!!!! if i>j ,ij=i*(i-1)/2+j 
!!!! if i<=j,ij=j*(j-1)/2+i
            implicit double precision (A-H,O-Z)
            if(i>j)then
                item=i*(i-1)/2+j
            else 
                item=j*(j-1)/2+i
            end if
            index2eri=item
        end function      
!======================================================================!
!======================================================================!
      subroutine matrix1D22D(ndim,one,two)
!!!! Restore the symmetric matrix stored in a one-dimensional array 
!!!!  to a two-dimensional symmetric matrix.
            implicit double precision (A-H,O-Z)
            dimension::one(ndim*(ndim+1)/2),two(ndim,ndim)
            ij=0
            do i=1,ndim
              do j=1,i
                  ij=ij+1
                  two(i,j)=one(ij)
                  two(j,i)=one(ij)
              end do
            end do 
            return
      end subroutine      
!======================================================================!
!======================================================================!
      function krodelta(i,j)
!!!! if i=j  ,delta =1
!!!! if i/=j ,delta =0 
            implicit double precision (A-H,O-Z)
            if(i==j)then
               krodelta=1
            else 
               krodelta=0
            end if
            return
        end function    krodelta 
!======================================================================!
      End Module Global

 

 
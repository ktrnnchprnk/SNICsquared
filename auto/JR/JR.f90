!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   Simplified attempt at Jansen-Rit model 
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION A,small_a,B,small_b,c1,c2,c3,c4,&
      vmax1,vmax2,vmax3,v01,v02,v03,r1,r2,r3,I,sgm1,sgm2,sgm3,&
      x1,x2,x3,x4,x5,x6
      x1=U(1)
      x2=U(2)
      x3=U(3)
      x4=U(4)
      x5=U(5)
      x6=U(6)
  

      A=0.334
      c3=3.4
      
      small_a=1
      B=1.5
      small_b=0.5
     
      c1=0.469
      c2=1.5
      c4=PAR(2)
      
      vmax1=1.
      vmax2=1.
      vmax3=0.5
      
      v01=0.084
      v02=0.084
      v03=0.5
      
      r1=40
      r2=40
      r3=30
      
      I=PAR(1)
      
      sgm1 = vmax1 / (1 + exp((v01 - (c2*x2 - c4*x3)) * r1));
      sgm2 = vmax2 / (1 + exp((v02 - (c1*x1)) * r2));
      sgm3 = vmax3 / (1 + exp((v03 - (c3*x1)) * r3));

       F(1)= x4
       F(2)=x5
       F(3)= x6
       F(4)= A*small_a*(I+sgm1)-2*small_a*x4-small_a**2*x1 
       F(5)= A*small_a*(sgm2)-2*small_a*x5-small_a**2*x2 
       F(6)= B*small_b*(sgm3)-2*small_b*x6-small_b**2*x3 


      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
       
       PAR(1)=2 !I
       PAR(2)=1!c4
       
       U(1)=  0.599999999999847
       U(2)=0.299888372693169
       U(3)=1.3
       U(4)=0
       U(5)=0
       U(6)=0
      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS

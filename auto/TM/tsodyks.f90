!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   tsodyks model
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION omega, tau, taur, I, uu, Tr, alpha, r01, vmax1,v01,r02, vmax2, v02, v,m

      v=U(1)
      m=U(2)
  
      omega = PAR(2)
      tau = 0.2
      taur = 0.8
      I = PAR(1)
      uu = 0.04
      r01 = 0.3
      vmax1 = 60
      v01 = 15
      r02= 3
      vmax2 = 60
      v02 = 10


       F(1)= (-v + m* uu *omega * (vmax1 / (1 + exp(r01 * (v01 - v)))) + sqrt(tau) + I) / tau
       F(2) = (1 - m) / taur - uu * m * (vmax2 / (1 + exp(r02 * (v02 - v))))


      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

       PAR(1)=-50.
       PAR(2)=40
  

       U(1)=-49.552786404486000
       U(2)=1

      END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS

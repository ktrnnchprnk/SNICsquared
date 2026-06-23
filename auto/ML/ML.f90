!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!           Morris Lecar model
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION V, N, C, gL, gCa, gK, VL, VCa, VK, V1, V2, V3, V4, I, Mss, tau_f, Nss
     
      
      V = U(1)
      N = U(2)
      
      I = PAR(1)
      C = 1
      gL = 2
      gCa = PAR(2)
      gK = PAR(3)
      VL = -60
      VCa = 127
      VK = -90
      V1 = -1.23
      V2 = 20
      V3 = -20
      V4 = 7
      
      Mss = (1+TANH((V-V1)/V2))/2
      Nss = (1+TANH((V-V3)/V4))/2
      tau_f = COSH((V-V3)/V4)
      F(1) = (-gL*(V-VL)-gCa*Mss*(V-VCa)-gK*N*(V-VK)+I)/C
      F(2) = tau_f*(Nss-N);
      
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
      !     ----------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

       PAR(1)=20 !I
       PAR(2)=15 !gCa
       PAR(3)=12 !gK
       
       
       U(1)=    17.3427
       U(2)= 1

       
       END SUBROUTINE STPNT

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT
      SUBROUTINE PVLS(NDIM,U,PAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

      DOUBLE PRECISION, EXTERNAL :: GETP,GETU2
      INTEGER NDX,NCOL,NTST
      PAR(4)=GETP('STA',0,U)

     
      END SUBROUTINE PVLS

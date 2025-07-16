!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!           Wilson Cowan Model 18/02
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION E,I,cee,cei,cie,Kp,ae,ai,thetae,thetai,sgme,sgmi,alpha,Fe,Fi
     
      
      E = U(1)
      I = U(2)
      
      Kp = PAR(1)
      alpha=0.9
      cee=18
      cei=14
      cie=PAR(2)
      ae=1.3
      ai=2
      thetae=4
      thetai=3.7
      
      Fe = cee*E-cie*I+Kp*alpha
      Fi = cei*E + Kp*(1-alpha)
      sgme = 1/(1+exp(-ae*(Fe-thetae)))-1/(1+exp(ae*thetae))
      sgmi = 1/(1+exp(-ai*(Fi-thetai)))-1/(1+exp(ai*thetai))

      F(1) = -E + (1-E)*sgme
      F(2) = -I + (1-I)*sgmi;
      
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
      !     ----------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

       PAR(1)=0 !Kp
       PAR(2)=35 !beta1       
       
       U(1)=0
       U(2)=0 

       
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
      PAR(3)=GETP('STA',0,U)

     
      END SUBROUTINE PVLS

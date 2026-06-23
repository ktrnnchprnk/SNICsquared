!----------------------------------------------------------------------
!           MEAN-FIELD WITH ASTROCYTES (translated from MATLAB astro_v1)
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      DOUBLE PRECISION :: E, Ipop, A
      DOUBLE PRECISION :: cee, cii, cei, cie, caa
      DOUBLE PRECISION :: cae, cai
      DOUBLE PRECISION :: cea, cia
      DOUBLE PRECISION :: ae, ai, aa
      DOUBLE PRECISION :: thetae, thetai, thetaa
      DOUBLE PRECISION :: Kp, alph
      DOUBLE PRECISION :: F_e, F_i, F_a
      DOUBLE PRECISION :: phi_e, phi_i, phi_a

!----- Map variables -----
      E    = U(1)   ! excitatory population
      Ipop = U(2)   ! inhibitory population
      A    = U(3)   ! astrocytic population

!----- Map parameters -----

      
      cii    = 0
      caa    = 0
      
      thetae = 4
      thetai = 3.7
      thetaa = 3
      ae     = 1.3
      ai     = 2
      aa     = 1.3
      
      cae   = PAR(1)
      cai   =PAR(2)
      cea    = 20 !20
      cia    = 7 !7
      
      cee    = 12
      cei    = 5.6
      cie    = 18.5
   
      Kp     = PAR(4)
      alph   = 1

!----- Compute inputs -----
      F_e = cee*E - cie*Ipop + cae*A  + Kp*alph
      F_i = -cii*Ipop + cei*E + cai*A +(1-alph)*Kp
      F_a = cea*E - cia*Ipop + caa*A 

!----- Sigmoid helper function (inline) -----
      phi_e = 1.0D0/(1.0D0 + DEXP(-ae*(F_e - thetae))) &
     &        - 1.0D0/(1.0D0 + DEXP(ae*thetae))

      phi_i = 1.0D0/(1.0D0 + DEXP(-ai*(F_i - thetai))) &
     &        - 1.0D0/(1.0D0 + DEXP(ai*thetai))

      phi_a = 1.0D0/(1.0D0 + DEXP(-aa*(F_a - thetaa))) &
     &        - 1.0D0/(1.0D0 + DEXP(aa*thetaa))

!----- ODE system -----
      F(1) = -E    + (1.0D0 - E   )*phi_e
      F(2) = -Ipop + (1.0D0 - Ipop)*phi_i
      F(3) = -A    + (1.0D0 - A   )*phi_a

      RETURN
      END SUBROUTINE FUNC

!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

!----- Set parameters -----
      PAR(1)=0    ! cae
      PAR(2)=16 ! cai
      PAR(4) = 0.78
      
U(1) = 0.012368
U(2) = 0.000235
U(3) = 0.007212

      RETURN
      END SUBROUTINE STPNT

!----------------------------------------------------------------------

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

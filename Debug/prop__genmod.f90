        !COMPILER-GENERATED INTERFACE MODULE: Sun Aug 22 23:24:20 2021
        MODULE PROP__genmod
          INTERFACE 
            SUBROUTINE PROP(IP,XP,P,R1,R2,XCOS,NPROPIN,NPRT)
              INTEGER(KIND=4) :: NPRT
              INTEGER(KIND=4) :: NPROPIN
              INTEGER(KIND=4) :: IP(1:NPROPIN)
              REAL(KIND=8) :: XP(1:NPROPIN)
              REAL(KIND=8) :: P(NPROPIN)
              REAL(KIND=8) :: R1
              REAL(KIND=8) :: R2
              REAL(KIND=8) :: XCOS
            END SUBROUTINE PROP
          END INTERFACE 
        END MODULE PROP__genmod

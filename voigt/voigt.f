c     stress, stiffness
      SUBROUTINE VOIGT(T1,T2,C2,C4,IOPT)
      IMPLICIT NONE
      REAL*8 T1(6),T2(3,3),C2(6,6),C4(3,3,3,3)
      INTEGER IJV(6,2),I,J,IOPT,I1,I2,J1,J2,N,M
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,1,2,1,3,2,3/

      IF(IOPT.EQ.1) THEN
      DO I=1,6
         I1=IJV(I,1)
         I2=IJV(I,2)
         T2(I1,I2)=T1(I)
         T2(I2,I1)=T1(I)
      ENDDO
      ENDIF
C
      IF(IOPT.EQ.2) THEN
         DO I=1,6
            I1=IJV(I,1)
            I2=IJV(I,2)
            T1(I)=T2(I1,I2)
         ENDDO
         ENDIF
C
      IF (IOPT.EQ.3) THEN
         DO I=1,6
            I1=IJV(I,1)
            I2=IJV(I,2)
            DO J=1,6
               J1=IJV(J,1)
               J2=IJV(J,2)
               C4(I1,I2,J1,J2)=C2(I,J)
               C4(I2,I1,J1,J2)=C2(I,J)
               C4(I1,I2,J2,J1)=C2(I,J)
               C4(I2,I1,J2,J1)=C2(I,J)
            ENDDO
         ENDDO
      ENDIF
C
      IF(IOPT.EQ.4) THEN
         DO  I=1,6
            I1=IJV(I,1)
            I2=IJV(I,2)
            DO  J=1,6
               J1=IJV(J,1)
               J2=IJV(J,2)
               C2(I,J)=C4(I1,I2,J1,J2)
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END

C----------------------------------------------------------------------C
C
C        SUBROUTINE pnpoly
C
C        PURPOSE
C           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON
C
C        USAGE
C           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )
C
C        DESCRIPTION OF THE PARAMETERS
C           PX      - X-COORDINATE OF POINT IN QUESTION.
C           PY      - Y-COORDINATE OF POINT IN QUESTION.
C           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF
C                     VERTICES OF POLYGON.
C           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF
C                     VERTICES OF POLYGON.
C           N       - NUMBER OF VERTICES IN THE POLYGON.
C           INOUT   - THE SIGNAL RETURNED:
C                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,
C                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,
C                      1 IF THE POINT IS INSIDE OF THE POLYGON.
C
C        REMARKS
C           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.
C           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY
C           OPTIONALLY BE INCREASED BY 1.
C           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING
C           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX
C           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING
C           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.
C           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.
C           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM
C           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT
C           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE
C           POINT IS INSIDE OF THE POLYGON.
C
C     ..................................................................
C
      SUBROUTINE pnpoly(PX,PY,XX,YY,N,INOUT)
      IMPLICIT NONE
      INTEGER MAXDIM
      PARAMETER (MAXDIM=1000)
      REAL*8 X(MAXDIM),Y(MAXDIM),XX(N),YY(N),PX,PY,var
      INTEGER N,INOUT,I,J
      LOGICAL MX,MY,NX,NY

C Check that number of points in clipping path is less than MAXDIM
      IF (N.GT.MAXDIM) THEN
          WRITE(0,7) N
7         FORMAT('WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY')
          RETURN
      ENDIF

C Compute the position of the points on the clipping path (xx,yy)
C relative to the point of interest (px,py)
      DO 1 I=1,N
          X(I)=XX(I)-PX
          Y(I)=YY(I)-PY
    1 CONTINUE

      INOUT=-1
      DO 2 I=1,N
          J=1+MOD(I,N)
          MX=X(I).GE.0.0 ! Is the point on the clipping path east of the point of interest?
          NX=X(J).GE.0.0 ! Is the next point on the clipping path east of the point of interest?
          MY=Y(I).GE.0.0 ! Is the point on the clipping path north of the point of interest?
          NY=Y(J).GE.0.0 ! Is the next point on the clipping path north of the point of interest?
          IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) THEN
              GOTO 2
          ENDIF
          IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) THEN
              GOTO 3
          ENDIF
          INOUT=-INOUT
          GOTO 2
3         var = (Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))
          if (var.lt.0.0d0) then
              goto 2
          elseif (dabs(var).lt.1.0d-8) then
              goto 4
          else
              goto 5
          endif
C3         IF ((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5
4         INOUT=0
          RETURN
5         INOUT=-INOUT
2     CONTINUE
      RETURN
      END

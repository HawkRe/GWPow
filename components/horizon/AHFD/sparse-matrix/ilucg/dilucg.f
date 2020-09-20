      SUBROUTINE DILUCG(N,
     $                  IA,JA,A,
     $                  B,X,
     $                  ITEMP,RTEMP,
     $                  EPS,MAXITER,
     $                  ISTATUS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IA(*),JA(*),A(*),B(*),X(*),ITEMP(*),RTEMP(*)
C
C     INCOMPLETE LU DECOMPOSITION-CONJUGATE GRADIENT
C     -          --               -         -
C WHERE:
C     |N| IS THE NUMBER OF EQUATIONS.  IF N < 0, ITEMP AND
C       RTEMP CONTAIN THE ILU FROM A PREVIOUS CALL AND
C       B AND X ARE THE NEW RHS AND INITIAL GUESS.
C     IA IS AN INTEGER ARRAY DIMENSIONED |N|+1.  IA(I) IS THE
C       INDEX INTO ARRAYS JA AND A OF THE FIRST NON-ZERO
C       ELEMENT IN ROW I.  LET MAX=IA(|N|+1)-IA(1).
C     JA IS AN INTEGER ARRAY DIMENSIONED MAX.  JA(K) GIVES
C       THE COLUMN NUMBER OF A(K).
C     A IS A DOUBLE PRECISION ARRAY DIMENSIONED MAX.  IT CONTAINS THE
C       NONZERO ELEMENTS OF THE MATRIX STORED BY ROW.
C     B CONTAINS THE RHS VECTOR.
C     X IS A DOUBLE PRECISION ARRAY DIMENSIONED |N|.  ON ENTRY, IT CONTAINS
C       AN INITIAL ESTIMATE; ON EXIT, THE SOLUTION.
C     ITEMP IS AN INTEGER SCRATCH ARRAY DIMENSIONED 3*(|N|+MAX)+2.
C     RTEMP IS A DOUBLE PRECISION SCRATCH ARRAY DIMENSIONED 4*|N|+MAX.
C     EPS IS THE CONVERGENCE CRITERIA.  IT SPECIFIES THE RELATIVE
C       ERROR ALLOWED IN THE SOLUTION.  TO BE PRECISE, CONVERGENCE
C       IS DEEMED TO HAVE OCCURED WHEN THE INFINITY-NORM OF THE
C       CHANGE IN THE SOLUTION IN ONE ITERATION IS .LE. EPS * THE
C       INFINITY-NORM OF THE CURRENT SOLUTION.  HOWEVER, IF EPS
C       .LT. 0.0D0, IT IS INTERNALLY SCALED BY THE MACHINE PRECISION,
C       SO THAT, FOR EXAMPLE, EPS = -256.0D0 WILL ALLOW THE LAST 8 BITS
C       OF THE SOLUTION TO BE IN ERROR.
C     MAXITER GIVES THE REQUESTED NUMBER OF ITERATIONS,
C       OR IS 0 FOR "NO LIMIT".
C     ISTATUS IS AN INTEGER VARIABLE, WHICH IS SET TO:
C       -I IF THERE IS AN ERROR IN THE MATRIX STRUCTURE IN ROW I
C          (SUCH AS A ZERO ELEMENT ON THE DIAGONAL).
C        0 IF THE ITERATION FAILED TO REACH THE CONVERGENCE CRITERION
C          IN ITER ITERATIONS.
C       +I IF THE ITERATION CONVERGED IN I ITERATIONS.
C
C MODIFICATIONS
C * RETURN LOGICAL VALUE, J. THORNBURG, 13/MAY/85.
C * ADD CONVERGENCE CRITERIA, J. THORNBURG, 17/MAY/85.
C * CHANGED CONVERGENCE CRITERIA, J. THORNBURG, 4/AUG/89.
C * AVOID SKIP RETURNS, J. THORNBURG, 10/SEP/89.
C * CHANGE FUNCTION BACK TO SUBROUTINE, J. THORNBURG, 3.JUN.2003
C * CHANGE CONVERGENCE CRITERIA TO AVOID SPURIOUSLY STOPPING TOO SOON
C   ON SMALL RHS VECTORS -- J. THORNBURG, 11.JUN.2003
C
C REFERENCE:
C     D.S.KERSHAW,"THE INCOMPLETE CHOLESKY-CONJUGATE GRADIENT
C       METHOD FOR INTERATIVE SOLUTION OF LINEAR EQUATIONS",
C       J.COMPUT.PHYS. JAN 1978 PP 43-65
C
      LOGICAL DLU0
      NP=IABS(N)
      ISTATUS=0
      IF (NP.EQ.0) GO TO 20
C CALCULATE INDICES FOR BREAKING UP TEMPORARY ARRAYS.
      N1=NP+1
      MAX=IA(N1)-IA(1)
      ILU=1
      JLU=ILU+N1
      ID=JLU+MAX
      IC=ID+NP
      JC=IC+N1
      JCI=JC+MAX
      IR=1
      IP=IR+NP
      IS1=IP+NP
      IS2=IS1+NP
      IALU=IS2+NP
      IF (N.LT.0) GO TO 10
C DO INCOMPLETE LU DECOMPOSITION
      IF (DLU0(NP,IA,JA,A,ITEMP(IC),ITEMP(JC),ITEMP(JCI),RTEMP(IALU),
     *    ITEMP(ILU),ITEMP(JLU),ITEMP(ID),RTEMP(IR),IERROR)) GOTO 20
C AND DO CONJUGATE GRADIENT ITERATIONS
C CALL MODIFIED TO ADD ADJUSTABLE CONVERGENCE CRITERIA EPS
C - J. THORNBURG, 17/MAY/85.
10    CALL DNCG0(NP,IA,JA,A,B,X,ITEMP(ILU),ITEMP(JLU),ITEMP(ID),
     *  RTEMP(IALU),RTEMP(IR),RTEMP(IP),RTEMP(IS1),RTEMP(IS2),
     *  EPS,MAXITER,ITER)
C ITER IS ACTUAL NUMBER OF ITERATIONS (NEGATIVE IF NO CONVERGENCE)
      ISTATUS = ITER
      IF (ITER .LT. 0) ISTATUS = 0
      RETURN
C ERROR RETURN FROM INCOMPLETE LU DECOMPOSITION
20    ISTATUS = -IERROR
      RETURN
      END

      LOGICAL FUNCTION DLU0(N,IA,JA,A,IC,JC,JCI,ALU,ILU,JLU,ID,V,IE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IA(*),JA(*),A(*),IC(*),JC(*),JCI(*),
     *  ALU(*),ILU(*),JLU(*),ID(N),V(N)
      LOGICAL NODIAG
      COMMON /ICBD00/ ICBAD
C     INCOMPLETE LU DECOMPOSITION
C WHERE:
C     N,IA,JA, AND A ARE DESCRIBED IN SUBROUTINE ILUCG
C     IC IS AN INTEGER ARRAY DIMENSIONED N+1, IC(J) GIVES THE
C       INDEX OF THE FIRST NONZERO ELEMENT IN COLMN J IN
C       ARRAY JC.
C     JC IS AN INTEGER ARRAY WITH THE SAME DIMENSION AS A.
C       JC(K) GIVES THE ROW NUMBER OF THE K'TH ELEMENT IN
C       THE COLUMN STRUCTURE.
C     JCI IS AN INTEGER ARRAY WITH THE SAME DIMENSION AS A.
C       JCI(K) GIVES THE INDEX INTO ARRAY A OF THE K'TH ELEMENT
C       OF THE COLUMN STRUCTURE.
C     ALU HAS THE SAME DIMENSION AS A.  ON EXIT, IT WILL
C       CONTAIN THE INCOMPLETE LU DECOMPOSITION OF A WITH THE
C       RECIPROCALS OF THE DIAGONAL ELEMENTS OF U.
C     ILU AND JLU CORRESPONDS TO IA AND JA BUT FOR ALU.
C     ID IS AN INTEGER ARRAY DIMENSIONED N.  IT CONTAINS
C       INDICES TO THE DIAGONAL ELEMENTS OF U.
C     V IS A REAL SCRATCH VECTOR OF LENGTH N.
C     IE GIVES THE ROW NUMBER IN ERROR IF AN ERROR OCCURED
C       (RETURN VALUE .TRUE.), OR IS UNUSED IF ALL IS WELL
C       (RETURN VALUE .FALSE.).
C
C     RETURN VALUE = .FALSE. IF ALL IS WELL, .TRUE. IF ERROR.
C
C NOTE: DLU0 SETS ARGUMENTS IC THROUGH V.
C
      ICBAD=0
C ZERO COUNT OF ZERO DIAGONAL ELEMENTS IN U.
C
C FIRST CHECK STRUCTURE OF A AND BUILD COLUMN STRUCTURE
      DO 10 I=1,N
        IC(I)=0
10    CONTINUE
      DO 30 I=1,N
        KS=IA(I)
        KE=IA(I+1)-1
        NODIAG=.TRUE.
        DO 20 K=KS,KE
          J=JA(K)
          IF (J.LT.1.OR.J.GT.N) GO TO 210
          IC(J)=IC(J)+1
          IF (J.EQ.I) NODIAG=.FALSE.
20      CONTINUE
        IF (NODIAG) GO TO 210
30    CONTINUE
C MAKE IC INTO INDICES
      KOLD=IC(1)
      IC(1)=1
      DO 40 I=1,N
        KNEW=IC(I+1)
        IF (KOLD.EQ.0) GO TO 210
        IC(I+1)=IC(I)+KOLD
        KOLD=KNEW
40    CONTINUE
C SET JC AND JCI FOR COLUMN STRUCTURE
      DO 60 I=1,N
        KS=IA(I)
        KE=IA(I+1)-1
        DO 50 K=KS,KE
          J=JA(K)
          L=IC(J)
          IC(J)=L+1
          JC(L)=I
          JCI(L)=K
50      CONTINUE
60    CONTINUE
C FIX UP IC
      KOLD=IC(1)
      IC(1)=1
      DO 70 I=1,N
        KNEW=IC(I+1)
        IC(I+1)=KOLD
        KOLD=KNEW
70    CONTINUE
C FIND SORTED ROW STRUCTURE FROM SORTED COLUMN STRUCTURE
      NP=N+1
      DO 80 I=1,NP
        ILU(I)=IA(I)
80    CONTINUE
C MOVE ELEMENTS, SET JLU AND ID
      DO 100 J=1,N
        KS=IC(J)
        KE=IC(J+1)-1
        DO 90 K=KS,KE
          I=JC(K)
          L=ILU(I)
          ILU(I)=L+1
          JLU(L)=J
          KK=JCI(K)
          ALU(L)=A(KK)
          IF (I.EQ.J) ID(J)=L
90      CONTINUE
100   CONTINUE
C RESET ILU (COULD JUST USE IA)
      DO 110 I=1,NP
        ILU(I)=IA(I)
110   CONTINUE
C FINISHED WITH SORTED COLUMN AND ROW STRUCTURE
C
C DO LU DECOMPOSITION USING GAUSSIAN ELIMINATION
      DO 120 I=1,N
        V(I)=0.0D0
120   CONTINUE
      DO 200 IROW=1,N
        I=ID(IROW)
        PIVOT=ALU(I)
        IF (PIVOT.NE.0.0D0) GO TO 140
C THIS CASE MAKES THE ILU LESS ACCURATE
        ICBAD=ICBAD+1
        KS=ILU(IROW)
        KE=ILU(IROW+1)-1
        DO 130 K=KS,KE
          PIVOT=PIVOT+DABS(ALU(K))
130     CONTINUE
        IF (PIVOT.EQ.0.0D0) GO TO 220
140     PIVOT=1.0D0/PIVOT
        ALU(I)=PIVOT
        KKS=I+1
        KKE=ILU(IROW+1)-1
        IF (KKS.GT.KKE) GO TO 160
        DO 150 K=KKS,KKE
          J=JLU(K)
          V(J)=ALU(K)
150     CONTINUE
C FIX L IN COLUMN IROW AND DO PARTIAL LU IN SUBMATRIX
160     KS=IC(IROW)
        KE=IC(IROW+1)-1
        DO 190 K=KS,KE
          I=JC(K)
          IF (I.LE.IROW) GO TO 190
          LS=ILU(I)
          LE=ILU(I+1)-1
          DO 180 L=LS,LE
            J=JLU(L)
            IF (J.LT.IROW) GO TO 180
            IF (J.GT.IROW) GO TO 170
            AMULT=ALU(L)*PIVOT
            ALU(L)=AMULT
            IF (AMULT.EQ.0.0) GO TO 190
            GO TO 180
170         IF (V(J).EQ.0.0D0) GO TO 180
            ALU(L)=ALU(L)-AMULT*V(J)
180       CONTINUE
190     CONTINUE
C RESET V
        IF (KKS.GT.KKE) GO TO 200
        DO 195 K=KKS,KKE
          J=JLU(K)
          V(J)=0.0D0
195     CONTINUE
200   CONTINUE
C NORMAL RETURN
      DLU0 = .FALSE.
      RETURN
C ERROR RETURNS
210   IE=I
      DLU0 = .TRUE.
      RETURN
220   IE=IROW
      DLU0 = .TRUE.
      RETURN
      END

      SUBROUTINE DNCG0(N,IA,JA,A,B,X,ILU,JLU,ID,ALU,R,P,S1,S2,
     *  EPS,ITER,IE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IA(*),JA(*),A(*),B(N),X(N),ILU(*),JLU(*),ALU(*),ID(N),
     *  R(N),P(N),S1(N),S2(N)
C     NONSYMMETRIC CONJUGATE GRADIENT
C WHERE:
C     N,IA,JA,A,B, AND X ARE DESCRIBED IN SUBROUTINE DILUCG.
C     ILU GIVES INDEX OF FIRST NONZERO ELEMENT IN ROW OF LU.
C     JLU GIVES COLUMN NUMBER.
C     ID GIVES INDEX OF DIAGONAL ELEMENT OF U.
C     ALU HAS NONZERO ELEMENTS OF LU MATRIX STORED BY ROW
C       WITH RECIPROCALS OF DIAGONAL ELEMENTS OF U.
C     R,P,S1, AND S2 ARE VECTORS OF LENGTH N USED IN THE
C       ITERATIONS.
C FOLLOWING PARAMETER ADDED BY J. THORNBURG, 17/MAY/85.
C     EPS IS CONVERGENCE CRITERIA.  (DESCRIBED IN SUBROUTINE
C       DILUCG).
C     ITER IS MAX NUMBER OF ITERATIONS, OR 0 FOR "NO LIMIT".
C     IE GIVES ACTUAL NUMBER OF ITERATIONS, NEGATIVE IF
C       NO CONVERGENCE.
C
C R0=B-A*X0
      CALL DMUL10(N,IA,JA,A,X,R)
      DO 10 I=1,N
        R(I)=B(I)-R(I)
10    CONTINUE
C P0=(UT*U)(-1)*AT*(L*LT)(-1)*R0
C FIRST SOLVE L*LT*S1=R0
      CALL DSUBL0(N,ILU,JLU,ID,ALU,R,S1)
C TIMES TRANSPOSE OF A
      CALL DMUL20(N,IA,JA,A,S1,S2)
C THEN SOLVE UT*U*P=S2
      CALL DSUBU0(N,ILU,JLU,ID,ALU,S2,P)
      IE=0
C INPROD IS DOT PRODUCT ROUTINE IN *LIBRARY (UBC MATRIX P 28)
C ALL CALLS ON IT COMMENTED OUT AND REPLACED WITH CALLS TO NEW
C ROUTINE DGVV(...); SOURCE CODE FOR LATTER ADDED TO END OF
C THIS FILE.  - J. THORNBURG, 10/MAY/85.
C     CALL INPROD(R,S1,EDOT,RDOT,N)
      RDOT = DGVV(R,S1,N)
C LOOP BEGINS HERE
20    CALL DMUL30(N,ILU,JLU,ID,ALU,P,S2)
C     CALL INPROD(P,S2,EDOT,PDOT,N)
      PDOT = DGVV(P,S2,N)
C
      IF (PDOT.EQ.0.0D0) RETURN
C
      ALPHA=RDOT/PDOT
C EQUATION 9PA                      ALPHA=(R,LINV*R)/(P,UT*U*P)
C 1 + ||X CHANGE|| CHANGED TO ||X CHANGE|| -- J. THORNBURG, 11.JUN.2003
      XMAX=0.0D0
      XDIF=0.0D0
      DO 30 I=1,N
        AP=ALPHA*P(I)
        X(I)=X(I)+AP
C EQUATION 9PB                      X=X+ALPHA*P
        AP=DABS(AP)
        XX=DABS(X(I))
        IF (AP.GT.XDIF) XDIF=AP
        IF (XX.GT.XMAX) XMAX=XX
30    CONTINUE
      IE=IE+1
C
C CONVERGENCE TEST (CHANGED BY J. THORNBURG, 17/MAY/85, 4/AUG/89)
C
      IF ((EPS .GT. 0.0D0) .AND. (XDIF .LE. EPS * XMAX)) RETURN
      IF ((EPS .LT. 0.0D0) .AND. (XMAX + XDIF/DABS(EPS) .EQ. XMAX))
     *   RETURN
C
C EXCEEDED ITERATION LIMIT?
C
      IF ((ITER .NE. 0) .AND. (IE .GE. ITER)) GO TO 60
      CALL DMUL10(N,IA,JA,A,P,S2)
      DO 40 I=1,N
        R(I)=R(I)-ALPHA*S2(I)
C EQUATION 9PC                      R=R-ALPHA*A*P
40    CONTINUE
      CALL DSUBL0(N,ILU,JLU,ID,ALU,R,S1)
C     CALL INPROD(R,S1,EDOT,RRDOT,N)
      RRDOT = DGVV(R,S1,N)
      BETA=RRDOT/RDOT
C EQUATION 9PD                      BETA=(R+,LINV*R+)/(R,LINV*R)
      RDOT=RRDOT
      CALL DMUL20(N,IA,JA,A,S1,S2)
      CALL DSUBU0(N,ILU,JLU,ID,ALU,S2,S1)
      DO 50 I=1,N
        P(I)=S1(I)+BETA*P(I)
C EQUATION 9PE                      P=(UT*U)(-1)*AT*(L*LT)(-1)*R+BETA*P
50    CONTINUE
      GO TO 20
60    IE=-IE
      RETURN
      END

      SUBROUTINE DMUL10(N,IA,JA,A,B,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IA(*),JA(*),A(*),B(N),X(N)
C MULTIPLY A TIMES B TO GET X
C WHERE:
C     N IS THE ORDER OF THE MATRIX
C     IA GIVES INDEX OF FIRST NONZERO ELEMENT IN ROW
C     JA GIVES COLUMN NUMBER
C     A CONTAINS THE NONZERO ELEMENTS OF THE NONSYMMETRIC
C       MATRIX STORED BY ROW
C     B IS THE VECTOR
C     X IS THE PRODUCT (MUST BE DIFFERENT FROM B)
C
      DO 20 I=1,N
        KS=IA(I)
        KE=IA(I+1)-1
        SUM=0.0D0
        DO 10 K=KS,KE
          J=JA(K)
          SUM=SUM+A(K)*B(J)
10      CONTINUE
        X(I)=SUM
20    CONTINUE
      RETURN
      END

      SUBROUTINE DMUL20(N,IA,JA,A,B,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IA(*),JA(*),A(*),B(N),X(N)
C MULTIPLY TRANSPOSE OF A TIMES B TO GET X
C WHERE:
C     N IS THE ORDER OF THE MATRIX
C     IA GIVES INDEX OF FIRST NONZERO ELEMENT IN ROW
C     JA GIVES COLUMN NUMBER
C     A CONTAINS THE NONZERO ELEMENTS OF THE NONSYMMETRIC
C       MATRIX STORED BY ROW
C     B IS THE VECTOR
C     X IS THE PRODUCT (MUST BE DIFFERENT FROM B)
C
      DO 10 I=1,N
        X(I)=0.0D0
10    CONTINUE
      DO 30 I=1,N
        KS=IA(I)
        KE=IA(I+1)-1
        BB=B(I)
        DO 20 K=KS,KE
          J=JA(K)
          X(J)=X(J)+A(K)*BB
20      CONTINUE
30    CONTINUE
      RETURN
      END

      SUBROUTINE DMUL30(N,ILU,JLU,ID,ALU,B,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ILU(*),JLU(*),ID(*),ALU(*),B(N),X(N)
C MULTIPLY TRANSPOSE OF U TIMES U TIMES B TO GET X
C WHERE:
C     N IS THE ORDER OF THE MATRIX
C     ILU GIVES INDEX OF FIRST NONZERO ELEMENT IN ROW OF LU
C     JLU GIVES COLUMN NUMBER
C     ID GIVES INDEX OF DIAGONAL ELEMENT OF U
C     ALU HAS NONZERO ELEMENTS OF LU MATRIX STORED BY ROW
C       WITH RECIPROCALS OF DIAGONAL ELEMENTS
C     B IS THE VECTOR
C     X IS THE PRODUCT UT*U*B (X MUST BE DIFFERENT FROM B)
C
      DO 10 I=1,N
        X(I)=0.0D0
10    CONTINUE
      DO 50 I=1,N
        KS=ID(I)+1
        KE=ILU(I+1)-1
        DIAG=1.0D0/ALU(KS-1)
        XX=DIAG*B(I)
        IF (KS.GT.KE) GO TO 30
        DO 20 K=KS,KE
          J=JLU(K)
          XX=XX+ALU(K)*B(J)
20      CONTINUE
30      X(I)=X(I)+DIAG*XX
        IF (KS.GT.KE) GO TO 50
        DO 40 K=KS,KE
          J=JLU(K)
          X(J)=X(J)+ALU(K)*XX
40      CONTINUE
50    CONTINUE
      RETURN
      END

      SUBROUTINE DSUBU0(N,ILU,JLU,ID,ALU,B,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ILU(*),JLU(*),ID(*),ALU(*),B(N),X(N)
C DO FORWARD AND BACK SUBSTITUTION TO SOLVE UT*U*X=B
C WHERE:
C     N IS THE ORDER OF THE MATRIX
C     ILU GIVES INDEX OF FIRST NONZERO ELEMENT IN ROW OF LU
C     JLU GIVES COLUMN NUMBER
C     ID GIVES INDEX OF DIAGONAL ELEMENT OF U
C     ALU HAS NONZERO ELMENTS OF LU MATRIX STORED BY ROW
C       WITH RECIPROCALS OF DIAGONAL ELEMENTS OF U
C     B IS THE RHS VECTOR
C     X IS THE SOLUTION VECTOR
C
      NP=N+1
      DO 10 I=1,N
        X(I)=B(I)
10    CONTINUE
C FORWARD SUBSTITUTION
      DO 30 I=1,N
        KS=ID(I)+1
        KE=ILU(I+1)-1
        XX=X(I)*ALU(KS-1)
        X(I)=XX
        IF (KS.GT.KE) GO TO 30
        DO 20 K=KS,KE
          J=JLU(K)
          X(J)=X(J)-ALU(K)*XX
20      CONTINUE
30    CONTINUE
C BACK SUBSTITUTION
      DO 60 II=1,N
        I=NP-II
        KS=ID(I)+1
        KE=ILU(I+1)-1
        SUM=0.0D0
        IF (KS.GT.KE) GO TO 50
        DO 40 K=KS,KE
          J=JLU(K)
          SUM=SUM+ALU(K)*X(J)
40      CONTINUE
50      X(I)=(X(I)-SUM)*ALU(KS-1)
60    CONTINUE
      RETURN
      END

      SUBROUTINE DSUBL0(N,ILU,JLU,ID,ALU,B,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ILU(*),JLU(*),ID(*),ALU(*),B(N),X(N)
C DO FORWARD AND BACK SUBSTITUTION TO SOLVE L*LT*X=B
C WHERE:
C     N IS THE ORDER OF THE MATRIX
C     ILU GIVES INDEX OF FIRST NONZERO ELEMENT IN ROW LU
C     JLU GIVES THE COLUMN NUMBER
C     ID GIVES INDEX OF DIAGONAL ELEMENT OF U
C     ALU HAS NONZERO ELEMENTS OF LU MATRIX STORED BY ROW
C       DIAGONAL ELEMENTS OF L ARE 1.0 AND NOT STORED
C     B IS THE RHS VECTOR
C     X IS THE SOLUTION VECTOR
C
      NP=N+1
      DO 10 I=1,N
        X(I)=B(I)
10    CONTINUE
C FORWARD SUBSTITUTION
      DO 30 I=1,N
        KS=ILU(I)
        KE=ID(I)-1
        IF (KS.GT.KE) GO TO 30
        SUM=0.0D0
        DO 20 K=KS,KE
          J=JLU(K)
          SUM=SUM+ALU(K)*X(J)
20      CONTINUE
        X(I)=X(I)-SUM
30    CONTINUE
C BACK SUBSTITUTION
      DO 50 II=1,N
        I=NP-II
        KS=ILU(I)
        KE=ID(I)-1
        IF (KS.GT.KE) GO TO 50
        XX=X(I)
        IF (XX.EQ.0.0) GO TO 50
        DO 40 K=KS,KE
          J=JLU(K)
          X(J)=X(J)-ALU(K)*XX
40      CONTINUE
50    CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION DGVV(V,W,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(N),W(N)
C     SUBROUTINE TO COMPUTE DOUBLE PRECISION VECTOR DOT PRODUCT.
C
      SUM = 0.0D0
            DO 10 I = 1,N
            SUM = SUM + V(I)*W(I)
10          CONTINUE
      DGVV = SUM
      RETURN
      END

;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Marquardt Algorithm

;;; Solve for coefficients B[k] in the nonlinear regression equation
;;;   Y* = F(X ; B) using N data points for Y[i] and X[i]

;;; *** Warning ***
;;;   This procedure has changed from the description in AI Memo 774.
;;;   The procedure FUNC now takes different arguments.

;;; References
;;;   Philip R. Bevington,
;;;     Data Reduction and Error Analysis for the Physical Sciences,
;;;     McGraw-Hill, New York, 1969.
;;;   James Kuester and Joe Mize,
;;;     Optimization Techniques with Fortran,
;;;     McGraw-Hill, New York, 1973.

;;; This code is substantially derived from the BSOLVE program
;;; given in Kuester.  That program has been translated to LISP
;;; with minor changes.
;;; It has been translated from FORTRAN into LISP with minor changes.
;;; Therefore (c) Copyright McGraw-Hill 1973.
;;; *** Permission has not been sought to distribute this code. ***

;;; Copyright Gerald Roylance 1982, 1984, 1985

;;; Bugs and Fixes
;;;   eps and tau should be vectors

(in-package "CLMATH")


;;;; Marquardt Algorithm

;;; Marquardt is a combination of steepest descent and Gauss-Newton
;;;
;;;   Objective Function
;;;       [P' P + lambda I] delta-B^ = P' (Y-Y^*)
;;;  P' is P transposed
;;;  P is the N row (i), K col (j) matrix of partial Yi wrt Aj
;;;    lambda is added to the main diagonal of P' P
;;;    when lambda = infinity, the method is steepest descent
;;;    when lambda = 0, the method is Gauss-Newton

;;;  Marquart's procedure is
;;;    1. compute chi-square(b)
;;;    2. set lambda = 0.001
;;;    3. compute delta-b and chi-square(b + delta-b)
;;;    4. if chi-square(b + delta-b) > chi-square(b),
;;;          then increase lambda by a factor of 10 and repeat step 3
;;;    5. if the chi-square is smaller,
;;;          then decrease lambda by a factor of 10
;;;               set b = b + delta-b
;;;               and continue at 3


;;;; MARQUARDT

;;; N    -- number of data points
;;; K    -- number of unknowns
;;; B    -- vector of unknowns
;;; BMIN -- vector of minimum values of B
;;; BMAX -- vector of maximum values of B
;;; X    -- vector of independent variable data points
;;;         (set equal to b for root location)
;;; Y    -- vector of   dependent variable
;;;         (set equal to 0 for root location)
;;; Z    -- computed values of dependent variable
;;; BV   -- code vector
;;;           1.0 for numerical derivatives
;;;           0.0 for do not change this coefficient
;;;          -1.0 for analytical derivatives
;;;
;;; FUNC -- Function
;;;      (funcall FUNC  X Y  B)
;;;        sets Y[i] = f(x)         evaluated at x=X[i] for all i
;;;
;;; DERIV-- Derivative
;;;      (funcall DERIV X DY B j)
;;;        sets DY[i] = df(x)/dB[j] evaluated at x=X[i] for all i
;;;
;;; EPS
;;; TAU

;;; FUNC and DERIV have loops in them for 2 reasons
;;;   1.  The format of X is completely general -- it could be
;;;       a two dimension array, an array of flonums, or an array
;;;       of vectors.
;;;   2.  To avoid number consing in MACLISP

;;; convergence criteria:
;;;   either chi-square is extremely small
;;;   or
;;;      for all i at some step
;;;        abs(delta-b[i]) / (tau + abs(b[i])) < eps


;;;; Various Subroutines

(defun marq-copy-vector (to from)
  (declare (type (array float (*)) from to))
  (do ((i 0 (1+ i))
       (n (array-dimension to 0)))
      ((>= i n))
    (declare (fixnum i n))
    (setf (aref to   i)
	  (aref from i))))

(DEFUN ARCOS (X)
  (SETQ X (MIN 1.0 (MAX -1.0 X)))
  (atan (sqrt (- 1.0 (expt X 2))) X))

(defun marq-chi-square (n y z)
  (declare (fixnum n)
	   (type (array float (*)) y z))
  (let ((phi 0.0))
    (declare (float phi))
    (do ((j1 0 (1+ j1)))
	((>= j1 n))
      (declare (fixnum j1))
      (setq phi (+ phi (expt (- (aref z j1) (aref y j1)) 2))))
    phi))


;;;; Derivative Calculations

(defun marq-calc-deriv (k n b bv bmax x z p trial-b pkn1 func deriv y-temp)
  (declare (fixnum k n)
	   (type (array float   (*)) b bv bmax trial-b y-temp z pkn1)
	   (type (array float (* *)) p))
  (do ((j1   0 (1+ j1))				; calculate k partial derivatives
       (den  0.0))
      ((>= j1 k))
    (declare (fixnum j1) (float den))

    (cond ((< (aref bv j1) 0.0)			; analytic derivatives
	   (funcall deriv x y-temp b j1)	; evaluate derivatives
	   (do ((j2 0 (1+ j2)))			; and store into array
	       ((>= j2 n))
	     (declare (fixnum j2))
	     (setf (aref p      j2 j1)
		   (aref y-temp j2))))

	  ((= (aref bv j1) 0.0))

	  ((> (aref bv j1) 0.0)			; numerical derivatives
	   (do ((j2 0 (1+ j2)))			; copy b into trial-b
	       ((>= j2 k))
	     (declare (fixnum j2))
	     (setf (aref trial-b j2)
		   (aref b j2)))
	   (setq den (* 0.00001
			(max (aref pkn1 j1)
			     (abs (aref trial-b j1)))))
	   (cond ((> (+ (aref trial-b j1) den)
		     (aref bmax j1))
		  (setq den (- den))))
	   (setf (aref trial-b j1)
		 (+ (aref trial-b j1) den))

	   (funcall func x y-temp trial-b)	; evaluate function
	   (do ((j2 0 (1+ j2)))			; calculate numerical derivative
	       ((>= j2 n))
	     (declare (fixnum j2))
	     (setf (aref p j2 j1)
		   (/ (- (aref y-temp j2)
			 (aref z      j2))
		      den)))
	   ))
    ))


;;;; Set up correction equations

(defun marq-setup (k n ak1 bv p y z a)
  (declare (fixnum k n)
	   (type (array float (*)) ak1 bv y z)
	   (type (array float (* *)) p a))
  (do ((j1 0 (1+ j1)))
      ((>= j1 k))
    (declare (fixnum j1))
    (setf (aref ak1 j1) 0.0)
    (cond ((not (= (aref bv j1) 0.0))
	   
	   ;; calculate the right hand side
	   (do ((j2 0 (1+ j2)))
	       ((>= j2 n))
	     (declare (fixnum j2))
	     (setf (aref ak1 j1)
		   (+ (aref ak1 j1)
		      (* (aref p j2 j1)
			 (- (aref y j2)
			    (aref z j2))))))
	   
	   ;; calculate the left hand side
	   (do ((j2 0 (1+ j2)))
	       ((>= j2 k))
	     (declare (fixnum j2))
	     (setf (aref a j1 j2) 0.0)
	     (do ((j3 0 (1+ j3)))
		 ((>= j3 n))
	       (declare (fixnum j3))
	       (setf (aref a j1 j2)
		     (+ (aref a j1 j2) (* (aref p j3 j1) (aref p j3 j2))
			))))
	   ))
    
    (cond ((or (= (aref bv j1) 0.0)
	       (< (aref a j1 j1) 1.0e-20))
	   ;; make the equation trivial -- either close enough or don't change
	   (setf (aref ak1 j1) 0.0)
	   (do ((j2 0 (1+ j2)))
	       ((>= j2 K))
	     (declare (fixnum j2))
	     (setf (aref a j1 j2) 0.0))
	   (setf (aref a j1 j1) 1.0)))
    ))


;;;; Scale correction equations

(defun marq-scale (k a ak1 scale-b)
  (declare (fixnum k)
	   (type (array float   (*)) scale-b ak1)
	   (type (array float (* *)) a))

  (do ((j1 0 (1+ j1)))
      ((>= j1 k))
    (declare (fixnum j1))
    (setf (aref scale-b j1)
	  (sqrt (aref a j1 j1))))
  
  (do ((j1 0 (1+ j1)))
      ((>= j1 k))
    (declare (fixnum j1))
    (setf (aref ak1 j1)
	  (/ (aref ak1 j1)
	     (aref scale-b j1)))
    (do ((j2 0 (1+ j2)))
	((>= j2 k))
      (declare (fixnum j2))
      (setf (aref a j1 j2)
	    (/ (aref a j1 j2)
	       (* (aref scale-b j1)
		  (aref scale-b j2))))))
  )



;;;; Solve an equation

;;; I think this solves  AC * x = ACK1 for x and returns x in ACK1

(defun marq-solve (k ac ack1)
  (declare (fixnum k)
	   (type (array float   (*)) ack1)
	   (type (array float (* *)) ac))

  (do ((l1 0 (1+ l1)))
      ((>= l1 k))
    (declare (fixnum l1))

    (do ((l3 (+ l1 1) (1+ l3)))
	((>= l3 K))
      (declare (fixnum l3))
      (setf (aref ac l1 l3)
	    (/ (aref ac l1 l3)
	       (aref ac l1 l1))))

    (setf (aref ack1 l1)
	  (/ (aref ack1 l1)
	     (aref ac l1 l1)))

    (do ((l3 0 (1+ l3)))
	((>= l3 k))
      (declare (fixnum l3))

      (cond ((not (= l1 l3))
	     (do ((l4 (+ l1 1) (1+ l4)))
		 ((>= l4 K))
	       (declare (fixnum l4))

	       (setf (aref ac l3 l4)
		     (- (aref ac l3 l4)
			(* (aref ac l1 l4) (aref ac l3 l1)))))

	     (setf (aref ack1 l3)
		   (- (aref ack1 l3)
		      (* (aref ack1 l1) (aref ac l3 l1))))
	     )))))


;;;; Take a step ...

(defun marq-cal-gam (k trial-b b bmin bmax ak1 scale-b ack1 delta-b gn)
  (declare (fixnum k)
	   (float  gn)
	   (type (array float (*)) trial-b b bmin bmax scale-b ack1 delta-b))

  (let ((gamm 0.0)
	(dn   0.0)
	(dg   0.0)
	(COSG 0.0)
	(JGAM 0))
    (declare (fixnum jgam)
	     (float  cosg dn dg gamm))

    (do ((j1 0 (1+ j1)))
	((>= j1 k))
      (declare (fixnum j1))
      
      (setf (aref delta-b j1)
	    (/ (aref ack1 j1) (aref scale-b j1)))
      (setf (aref trial-b j1)
	    (max (aref bmin j1)
		 (min (aref bmax j1)
		      (+ (aref b j1) (aref delta-b j1)))))
      (setq dg (+ dg (* (aref delta-b j1)
			(aref ak1 j1)
			(aref scale-b j1))))
      (setq dn (+ dn (* (aref delta-b j1) (aref delta-b j1))))
      (setf (aref delta-b j1)
	    (- (aref trial-b j1) (aref b j1))))
    
    (setq cosg (/ dg (sqrt (* dn gn))))
    (cond ((<  cosg 0.0)
	   (setq jgam 2)
	   (setq cosg (- cosg)))
	  (t
	   (setq jgam 0)))
    (setq cosg (min cosg 1.0))
    (setq gamm (* (arcos cosg) (/ 180.0 3.14159265)))

    (cond ((> jgam 0) (setq gamm (- 180.0 gamm))))

    gamm))


;;; copy Newton array and do steepest decent adjustment

(defun marq-new-lambda (fl a ac ak1 ack1)
  (declare (float fl)
	   (type (array float   (*)) ak1 ack1)
	   (type (array float (* *)) a ac))
  (do ((j1 0 (1+ j1))
       (k  (array-dimension ak1 0)))
      ((>= j1 k))
    (declare (fixnum j1 k))
    
    (do ((j2 0 (1+ j2)))
	((>= j2 K))
      (declare (fixnum j2))
      (setf (aref ac j1 j2) (aref a j1 j2)))
    
    (setf (aref ack1  j1) (aref ak1 j1))
    (setf (aref ac j1 j1) (+ (aref ac j1 j1) fl))
    ))


;;;; Marquardt Main Line Function

(defun marquardt (n k X Y Z func deriv B Bmin Bmax Bv
		  &optional
		  (eps 0.00002)			; must be > 0.0
		  (tau 0.00100))		; must be > 0.0
  (declare (fixnum k n)
	   (float tau eps)
	   (type (array float (*)) b bmin bmax bv)
	   (type (array float (*)) y z)
	   (type t x))

  (if (or (<= tau 0.0) (<= eps 0.0))
      (error "MARQUARDT -- bad tau or eps"))

  (LET ((P        (make-array (list n k)	; partial derivatives
			      :element-type 'float
			      :initial-element 0.0))
	(TRIAL-B  (make-array k			; trial-b
			      :element-type 'float
			      :initial-element 0.0))
	(PKN1     (make-array k
			      :element-type 'float
			      :initial-element 0.0))
	(TRIAL-Z  (make-array n			; trial-z
			      :element-type 'float
			      :initial-element 0.0))
	(Y-TEMP   (make-array n			; y-temp
			      :element-type 'float
			      :initial-element 0.0))
	(A        (make-array (list k k)	; newton array
			      :element-type 'float
			      :initial-element 0.0))
	(ak1      (make-array k			; RHS of objective equation? scaled
			      :element-type 'float
			      :initial-element 0.0))
	(SCALE-B  (make-array k			; scale?
			      :element-type 'float
			      :initial-element 0.0))
	(AC       (make-array (list k k)	; newton + lambda I array
			      :element-type 'float
			      :initial-element 0.0))
	(ACK1     (make-array k			; scaled delta-b?
			      :element-type 'float
			      :initial-element 0.0))
	(DELTA-B  (make-array k			; delta-b
			      :element-type 'float
			      :initial-element 0.0))
	(ICON     k)
	(PH     0.0)
	(GAMM   0.0)
	(FNU   10.0)				; must be > 1.0 -- change in lambda
	(FLA    0.0)
	(PHMIN  0.0)
	(GN     0.0)
	(PHI    0.0))
    (declare (fixnum icon)
	     (float fla fnu gamm gn ph phi phmin)
	     (type (array float   (*)) trial-b pkn1 trial-z y-temp ak1 scale-b ack1 delta-b)
	     (type (array float (* *)) p a ac))

    (do ((i1 0 (1+ i1))
	 (ke 0))
	((>= i1 k)
	 (cond ((= ke 0) (error "MARQUARDT -- No Variables"))
	       ((> ke n) (error "MARQUARDT -- More variables than datapoints"))))
      (declare (fixnum i1 ke))
      (cond ((not (= (aref bv i1) 0.0)) (setq ke (1+ ke)))))

    (do ((j1 0 (1+ j1)))
	((>= j1 k))
      (declare (fixnum j1))
      (setf (aref trial-b j1) (aref b j1))
      (setf (aref pkn1 j1)
	    (+ (abs (aref b j1)) 1.0e-02)))	; *** magic number
	    

    ;; Iterate to Convergence

    ;; calculate initial values

    (funcall func x z trial-b)
    (setq ph   (marq-chi-square n y z))		; calculate CHISQUARE
    (setq icon (cond ((< ph 1.0e-10) 0)
		     (t              k)))

    (DO ((I 1 (1+ I)))
	((< ICON 1)
	 (COND ((= ICON -1) (format t "~% No Function Improvement Possible"))
	       ((= ICON -4) (format t "~% Converged but LAMBDA (FLA) large")))
	 nil)
      
      (declare (fixnum i))
      (cond ((<= fla 0.0) (setq fla  0.01000)))
      (setq phmin (max phmin 0.0))

      (cond ((not (and (> phmin ph) (> i 1)))
	     (marq-calc-deriv k n b bv bmax x z p trial-b pkn1 func deriv y-temp)))
      
      ;; set up correction equations
      (marq-setup k n ak1 bv p y z a)
		      
      (setq gn 0.0)
      (do ((j1 0 (1+ j1)))
	  ((>= j1 k))
	(declare (fixnum j1))
	(setq gn (+ gn (expt (aref ak1 j1) 2))))
		      
      ;; scale correction equations
      (marq-scale k a ak1 scale-b)
		      
      ;; Keep forcing steeper descent until chi-square improves

      (do ((fl (/ fla fnu) (* fnu fl)))
	  ((>= fl 1.0e8)
	   (setq icon -1)
	   (setq fla  fl)
	   NIL)
	(declare (float fl))
	
	(marq-new-lambda fl a ac ak1 ack1)	; new correction equations
	(marq-solve k ac ack1)			; solve the correction equations
	;; take the step and see if it is good
	(setq gamm (marq-cal-gam k trial-b b bmin bmax ak1 scale-b ack1 delta-b gn))
	(funcall func x trial-z trial-b)
	(setq phi (marq-chi-square n y trial-z))	; calculate CHISQUARE
	
	(cond ((< phi 1.0e-10)			; converged if chi-square very small
	       (setq icon 0)
	       (setq fla fl)	       
	       (return nil))
	      ((< phi ph)			; epsilon test
	       (setq icon 0)			;   eps > |delta-b| / (tau + |trial-b|)
	       (do ((j1 0 (1+ j1)))
		   ((>= j1 k))
		 (declare (fixnum j1))
		 (cond ((> (/ (abs (aref delta-b j1))
			      (+ tau (abs (aref trial-b j1))))
			   eps)
			(setq icon (1+ icon)))))
	       (cond ((= icon 0)		; gamma epsilon test
		      (cond ((and (> fl 1.0) (<= gamm 45.0)) (setq icon -4))))
		     ((and (> fl 1.0) (> gamm 90.0))	; gamma lambda test
		      (setq icon -1)))
	       (setq fla fl)
	       (return nil))))

      ;; we liked the step, so use new values of b, z, and ph
      ;;   *** but we might have taken the fl > 1.0e8 exit ...
      (marq-copy-vector b trial-b)
      (marq-copy-vector z trial-z)
      (setq ph phi)
      ;; (format t "~% ICON = ~D, PH = ~10E, ITERATION NO. = ~3D" icon ph i)
      )
    ph))


;;;; Testing

#|

(defun marq-test-func (x y b)
  (declare (type (array float (*)) x)
	   (type (array float (*)) y)
	   (type (array float (*)) b))
  (do ((i 0 (1+ i))
       (n (array-dimension y 0)))
      ((>= i n))
    (declare (fixnum i n))
    (setf (aref y i)
	  (+ (aref B 0)
	     (* (aref B 1)
		(EXP (* (aref B 2) (aref X I))))))))

(defun marq-test-deriv (x y b j)
  (declare (type (array float (*)) x)
	   (type (array float (*)) y)
	   (type (array float (*)) b)
	   (fixnum j))
  (do ((i 0 (1+ i))
       (n (array-dimension y 0)))
      ((>= i n))
    (declare (fixnum i n))
    (setf (aref y i)
	  (cond ((= j 0) 1.0)
		((= j 1) (EXP (* (aref B 2)
				 (aref X I))))
		((= j 2) (* (aref B 1)
			    (EXP (* (aref B 2)
				    (aref X I)))
			    (aref X I)))
		(t (error "foo"))))))

|#


;;;; Testing -- Continued

#|

(defun marq-test ()
  (LET* ((N    6)
	 (K    3)
	 (X    (make-array N :element-type 'float
			   :initial-contents '(  -5.0  -3.0  -1.0   1.0   3.0   5.0)))
	 (Y    (make-array N :element-type 'float
			   :initial-contents '( 127.0 151.0 379.0 421.0 460.0 426.0)))
	 (Z    (make-array N :element-type 'float
			   :initial-element 0.0))
	 (BV   (make-array K :element-type 'float
			   :initial-contents
			   ;; '(     1.0     1.0     1.0)	; numeric
			   '(    -1.0    -1.0    -1.0)	; analytic
			   ))
	 (B    (make-array K :element-type 'float
			   :initial-contents '(   400.0  -140.0    -0.13)))
	 (BMIN (make-array K :element-type 'float
			   :initial-contents '( -1000.0 -1000.0 -1000.0)))
	 (BMAX (make-array K :element-type 'float
			   :initial-contents '(  1000.0  1000.0  1000.0))))
    
    (format t "~% Calculated    Book")
    (mapc #'(lambda (actual book)
	      (format t "~% ~3,1,10$ ~3,1,10$"
		      actual book))
	  (cons (MARQUARDT N K X Y Z
			   'MARQ-TEST-FUNC 'MARQ-TEST-DERIV
			   B BMIN BMAX BV
			   0.00001 0.00100)
		(coerce B 'list))
	  (cons 13390.093
		'(523.29698 -156.93703 -0.19967593)))
    
    (format T "~%Y            Z")
    (mapc #'(lambda (actual book)
	      (format t "~% ~3,1,10$ ~3,1,10$"
		      actual book))
	  (coerce y 'list)
	  (coerce z 'list))
    NIL))

|#

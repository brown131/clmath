;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Fletcher Powell Functional Minimization

;;; Find the minimum of a multivariable, unconstrained, nonlinear function
;;;    F(X1, X2, ... , XN)

;;; References
;;;   Philip R. Bevington,
;;;     Data Reduction and Error Analysis for the Physical Sciences,
;;;     McGraw-Hill, New York, 1969.
;;;   David G. Luenberger,
;;;     Introduction to Linear and Nonlinear Programming,
;;;     Addison Wesley, Reading, MA, 1965.
;;;   James Kuester and Joe Mize,
;;;     Optimization Techniques with Fortran,
;;;     McGraw-Hill, New York, 1973.

;;; This code is a derivative of the procedure FMFP in Kuester.
;;; It has been translated from FORTRAN into LISP with minor changes.
;;; Therefore (c) Copyright McGraw-Hill 1973.
;;; *** Permission has not been sought to distribute this code. ***

;;; Gerald Roylance 1983, 1984

;;; Bugs and Fixes
;;;   fix H so there are separate direction, old gradient, and old x vectors
;;;   flush (go ... )
;;;   zero base arrays

(in-package "CLMATH")

;;;; Correspondence among References

;;; Correspondence between Luenberger and SSP
;;;
;;;	function		F	F
;;;	gradient		g	partialF
;;;	search direction	d	M
;;;	Hessian			S	H
;;;	step size		alpha
;;;	step vector		p	deltaX
;;;	gradient change		q	deltaG

;;; Algorithm
;;;
;;; Step 1.
;;;	set d[k] = -S[k] g[k]
;;;
;;; Step 2.
;;;	minimize f(x[k] + alpha d[k]), alpha >= 0
;;;	set p[k] = alpha[k] d[k]
;;;	find x[k+1]
;;;	find g[k+1]
;;;
;;; Step 3.
;;;	set q[k] = g[k+1] - g[k]
;;;	set S[k+1] = S[k] + (p[k] p[k]')/(p[k]' q[k])
;;;			  - (S[k] q[k] q[k]' S[k])/(q[k]' S[k] q[k])
;;;
;;;	update k and go to step 1


;;;; FMFP

(defun FMFP (FUNCT N X G EST EPS LIMIT)
  (declare (fixnum n)
	   (type   (array float (*)) x g)
	   (float  est eps)
	   (fixnum limit))
  ;; FUNCT	function, VAL <- (FUNCT N ARG GRAD)
  ;; N		number of independent variables
  ;; X		independent variable vector
  ;; F		final minimum value of objective function
  ;; G		final gradient vector at the minimum
  ;; VAL	current value of objective function
  ;; ARG	current vector of independent variable values
  ;; GRAD	current gradient vector values
  ;; H		Storage Vector
  ;;			   1 -  N : direction
  ;;			 N+1 - 2N : old gradient & delta gradient
  ;;			2N+1 - 3N : old x        & delta x
  ;;			3N+1 -    : positive definite matrix
  ;; EST	estimate of minimum value of objective function
  ;; EPS	test value representing the expected absolute error in movement
  ;; LIMIT	maximum number of iterations
  ;; IER	error parameter
  ;;			IER =  0 means convergence was obtained
  ;;			IER =  1 means no convergence in LIMIT iterations
  ;;			IER = -1 means errors in gradient calculations
  ;;			IER =  2 means likely that there is no minimum
  (do ((alfa   0.0)
       (H     (make-array (1+ (floor (* N (+ N 7)) 2))
			  :element-type    'float
			  :initial-element 0.0))
       (f     0.0)
       (ambda 0.0)
       (dalfa 0.0)
       (dx    0.0)
       (dy    0.0)
       (fx    0.0)
       (fy    0.0)
       (gnrm  0.0)
       (hnrm  0.0)
       (ier   0)
       (k     0)
       (kl    0)
       (kount 0)
       (n2    (* 2 n))
       (n3    (* 3 n))
       (n31   (1+ (* 3 n)))
       (nj    0)
       (oldf  0.0) (tt    0.0) (w     0.0) (z     0.0))
      ()
    (declare (fixnum k kl n2 n3 n31 nj)
	     (float alfa ambda dalfa dy dx fx fy gnrm hnrm oldf tt w z)
	     (type (array float (*)) H))

    ;;
    ;; compute function value and gradient vector for initial argument
    (setq f (funcall funct n x g))
    ;;
    ;; reset iteration counter and generate identity matrix
    s1
    (do ((j    1 (1+ j))
	 (k  n31 (+ kl 1))
	 (kl   0))				;0 just a dummy -- setq'd later
	((> j n))
      (declare (fixnum j k kl))
      (setf (aref h k) 1.0)
      (setq nj (- n j))
      (do ((l 1 (1+ l)))
	  ((> l nj))
	(setq kl (+ k l))
	(setf (aref h kl) 0.0)))
    ;;
    ;; start iteration loop
    s5
    (setq kount (1+ kount))
    ;;
    ;; save function value, argument vector and gradient vector
    (setq oldf f)
    (do ((j 1 (1+ j)))
	((> j n))
      (setf (aref h (+ n j))  (aref g j))
      (setf (aref h (+ n2 j)) (aref x j))
      ;;
      ;; determine direction vector H
      (setq k (+ j n3))
      (setq tt 0.0)
      (do ((l 1 (1+ l)))
	  ((> l n))
	(setq tt (- tt (* (aref g l) (aref h k))))
	(cond ((>= l j) (setq k (1+ k)))
	      (t        (setq k (+ k n (- l))))) )
      (setf (aref h j) tt))
    ;;
    ;; check whether function will decrease stepping along h
    (setq dy   0.0)
    (setq hnrm 0.0)
    (setq gnrm 0.0)
    ;;
    ;; calculate directional derivative and test values for direction
    ;; vector h and gradient vector g
    (do ((j 1 (1+ j)))
	((> j n))
      (setq hnrm (+ hnrm (abs (aref h j))))
      (setq gnrm (+ gnrm (abs (aref g j))))
      (setq dy   (+ dy (* (aref h j) (aref g j)))))
    ;;
    ;; repeat search in direction of steepest descent if
    ;;    directional derivative appears to be positive or zero or
    ;;    direction vector H is small compared to gradient vector G
    (cond ((or (>= dy 0.0) (<= (/ hnrm gnrm) eps)) (go s51)))
    ;;
    ;; search minimum along direction H
    ;;
    ;; search along H for positive directional derivative
    (setq fy f)
    (setq alfa (/ (* 2.0 (- est f)) dy))
    (setq ambda 1.0)
    ;;
    ;; use estimate for stepsize only if it is positive and less than 1.
    ;; otherwise take 1. as the stepsize
    (cond ((and (> alfa 0.0) (< alfa ambda))
	   (setq ambda alfa)))
    (setq alfa 0.0)

    loophead
    ;;
    ;; save function and derivative values for old argument
    (setq fx fy)
    (setq dx dy)
    ;;
    ;; step argument along h
    (do ((i 1 (1+ i)))
	((> i n))
      (setf (aref x i)
	    (+ (aref x i) (* ambda (aref h i)))))
    ;;
    ;; compute function value and gradient for new argument
    (setq f (funcall funct n x g))
    (setq fy f)
    ;;
    ;; compute directional derivative dy for new argument, terminate
    ;; search if dy is positive.  If dy is zero, the minimum is found
    (setq dy 0.0)
    (do ((i 1 (1+ i)))
	((> i n))
      (setq dy (+ dy (* (aref g i) (aref h i)))))
    (cond ((= dy 0.0) (go s36))
	  ((> dy 0.0) (go looptail)))
    ;;
    ;; terminate search also if the function value indicates that
    ;; a minimum has been passed
    (cond ((>= fy fx) (go looptail)))
    ;;
    ;; repeat search and double stepsize for further searches
    (setq ambda (+ ambda alfa))
    (setq alfa ambda)
    ;; terminate if the change in argument gets very large
    (cond ((> (* hnrm ambda) 1.0e10)
	   (error "IER=2 linear search: no minimum exists")))


    (go loophead)
    looptail

    ;; end of search loop
    ;;
    ;; interpolate cubically in the interval defined by the search
    ;; above and compute the argument x for which the interpolation
    ;; polynomial is minimized
    (setq tt 0.0)
    s23
    (cond ((= ambda 0.0) (go s36)))
    (setq z (+ (/ (* 3.0 (- fx fy)) ambda) dx dy))
    (setq alfa (max (abs z) (abs dx) (abs dy)))
    (setq dalfa (- (expt (/ z alfa) 2) (* (/ dx alfa) (/ dy alfa))))
    (cond ((< dalfa 0.0) (go s51)))
    (setq w (* alfa (sqrt dalfa)))
    (setq alfa (/ (* (+ dy w (- z)) ambda) (+ dy (* 2.0 w) (- dx))))
    (do ((i 1 (1+ i)))
	((> i n))
      (setf (aref x i)
	    (+ (aref x i) (* (- tt alfa) (aref h i)))))
    ;;
    ;; terminate if the value of the actual function at x is less
    ;; than the function values at the interval ends.  Otherwise reduce
    ;; the interval by choosing one end-point equal to x and repeat
    ;; the interpolation.  Which end-point is chosen depends on the
    ;; value of the function and its gradient at x
    ;;
    (setq f (funcall funct n x g))
    (cond ((or (>  f fx) (> f fy))
	   (setq dalfa 0.0)
	   (do ((i 1 (1+ i)))
	       ((> i n))
	     (setq dalfa (+ dalfa (* (aref g i) (aref h i)))))
	   (cond ((and (<  dalfa 0.0)
		       (<= f fx))
		  (cond ((and (= f fx) (= dx dalfa))
			 (go s36)))
		  (setq fx f)
		  (setq dx dalfa)
		  (setq tt alfa)
		  (setq ambda alfa)
		  (go s23))
		 ((not (and (= fy f) (= dy dalfa)))
		  (setq fy f)
		  (setq dy dalfa)
		  (setq ambda (- ambda alfa))
		  (setq tt 0.0)
		  (go s23)))
	   ))
    s36
    ;;
    ;; terminate if function has not decreased during last iteration
    (cond ((< (- oldf f) eps) (go s51)))
    ;;
    ;; compute difference vectors of argument and gradient from
    ;; two consecutive iterations
    (do ((j 1 (1+ j)))
	((> j n))
      (setq k (+ n j))
      (setf (aref h k) (- (aref g j) (aref h k)))
      (setq k (+ n k))
      (setf (aref h k) (- (aref x j) (aref h k))))
    ;;
    ;; test length of argument difference vector and direction vector
    ;; if at least n iterations have been executed.  Terminate if
    ;; both are less than eps
    (setq ier 0)
    (cond ((< kount n) (go s42)))
    (setq tt 0.0)
    (setq z 0.0)
    (do ((j 1 (1+ j)))
	((> j n))
      (setq k (+ n j))
      (setq w (aref h k))
      (setq k (+ k n))
      (setq tt (+ tt (abs (aref h k))))
      (setq z (+ z (* w (aref h k)))))
    (cond ((> hnrm eps) (go s42)))
    (cond ((<= tt eps) (go s56)))
    ;;
    ;; terminate if number of iterations would exceed limit
    s42
    (cond ((>= kount limit)
	   (error "FMFP -- IER=1  Number of iterations exceeded")
	   (setq ier 1)
	   (return f)))
    ;;
    ;; prepare updating of matrix H
    (setq alfa 0.0)
    (do ((j 1 (1+ j)))
	((> j n))
      (setq k (+ j n3))
      (setq w 0.0)
      (do ((l 1 (1+ l)))
	  ((> l n))
	(setq kl (+ n l))
	(setq w (+ w (* (aref h kl) (aref h k))))
	(cond ((>= l j) (setq k (1+ k)))
	      (t	    (setq k (+ k n (- l))))) )
      (setq k (+ n j))
      (setq alfa (+ alfa (* w (aref h k))))
      (setf (aref h j) w))
    ;;
    ;; repeat search in direction of steepest descent if results
    ;; are not satisfactory
    (cond ((= (* z alfa) 0.0) (go s1)))
    ;;
    ;; update matrix h
    (setq k n31)
    (do ((l 1 (1+ l)))
	((> l n))
      (setq kl (+ n2 l))
      (do ((j L (1+ j)))
	  ((> j n))
	(setq nj (+ n2 j))
	(setf (aref h k)
	      (+ (aref h k)
		  (/ (* (aref h kl) (aref h nj)) z)
		  (- (/ (* (aref h l) (aref h j)) alfa))))
	(setq k (1+ k))))
    (go s5)
    ;; end of iteration loop
    ;;
    ;; restore old values of function and arguments
    s51
    (do ((j 1 (1+ j)))
	((> j n))
      (setq k (+ n2 j))
      (setf (aref x j) (aref h k)))
    (setq f (funcall funct n x g))
    ;;
    ;; repeat in direction of steepest descent if derivative
    ;; fails to be sufficiently small
    (cond ((<= gnrm eps)
	   (setq ier 0))
	  ((>= ier 0)
	   (setq ier -1)
	   (go s1)))
    s56
    (cond ((not (= ier 0))
	   (format t "IER = ~d" ier)
	   (error "FMFP IER=nonzero -- somethings amiss")))
    (return f)
    ))


;;;; Sample usage

;;;; F = -3873.9, X1 = 0.20566, X2 = 0.47986
#|
(EVAL-WHEN (EVAL)
|#
#|
  (defun FUNCT (N ARG GRAD)
    (declare (ignore n)
	     (fixnum n)
	     (type   (array float (*)) arg grad))
|#
    ;; N is number of variables
    ;; ARG is X vector
    ;; GRAD is the gradient
#|
    (let* ((x1   (aref arg 1))
	   (x2   (aref arg 2))
	   (x12  (expt x1 2))
	   (x22  (expt x2 2))
	   (val  (+ 3803.84
		     (*  138.08 x1)
		     (*  232.92 x2)
		     (* -123.08 x12)
		     (* -203.64 x22)
		     (* -182.25 x1 x2)))
	   (dery1 (- 138.08
		      (* 2.0 123.08 x1)
		      (*     182.25 x2)))
	   (dery2 (- 232.92
		      (* 2.0 203.64 x2)
		      (*     182.25 x1))))
      (format t "~% FUNCT: x1 = ~11F, x2 = ~11F, F = ~11F" x1 x2 val)
      (setf (aref grad 1) (- dery1))
      (setf (aref grad 2) (- dery2))
      (- val)))
|#
#|  
  (defun MAIN-LINE ()
    (let* ((N          2)
	   (LIMIT    100)
	   (EST    -2000.0)
	   (EPS        0.0001)
	   (F          0.0)
	   (X      (make-array (1+ N) :element-type 'float :initial-element 0.0))
	   (G      (make-array (1+ N) :element-type 'float :initial-element 0.0)))
      (declare (fixnum n limit)
	       (float  est eps f)
	       (type   (array float (*)) X G))
      (setf (aref x 1) 1.0)
      (setf (aref x 2) 0.5)
      (format t "~% FLETCHER POWELL ALGORITHM~%	parameters:")
      (format t "~%		N     = ~7D" N)
      (format t "~% 		LIMIT = ~7D" LIMIT)
      (format t "~%	 	EST   = ~8F" EST)
      (format t "~% 		EPS   = ~8F" EPS)
      (format t "~%~%	Initial Values of X")
      (do ((j 1 (1+ j)))
	  ((> j N))
	(format t "~%		X(~2D) = ~A" j (aref X j)))
      (setq f (FMFP 'FUNCT N X G EST EPS LIMIT))
      (format t "~%~%Minimization Completed")
      (format t "~%	F     = ~8F" F)
      (format t "~%~%	Final Values of X")
      (do ((j 1 (1+ j)))
	  ((> j N))
	(format t "~%		X(~2D) = ~8F" j (aref X j)))
      ))
  )
|#
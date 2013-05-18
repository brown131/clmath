;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Various Integration Methods

;;; References
;;;   R. W. Hamming,
;;;     Numerical Methods for Scientist and Engineers,
;;;     McGraw-Hill, 1973.

;;;  (c) Copyright Gerald Roylance 1983, 1984
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   Should change all of them to take vectors as arguments.
;;;   Supply routines to integrate from X0 to X1 using RK and EULER
;;;   Put a prefix on all of these identifiers!
;;;   Multistep method
;;;     check interpolation code; there are erratic jumps in
;;;        integrating Y**2 + 1 that cause step size to be
;;;        decreased several times before advancing.
;;;        maybe error can be attributed to roundoff.
;;;     only uses a first order formula, should use a second or third
;;;     Should return a few values!

(in-package "CLMATH")

;;;; Eulers Method

;;; forward Euler

(defun euler (funct xn yn h)
  (+ yn (* h (funcall funct xn yn))))

(defun integrate-euler (funct x0 y0 h x1)
  (do ((x x0 (+ x h))
       (y y0))
      ((<= (- x1 x) h)
       (euler funct x y (- x1 x)))
    (setq y (euler funct x y h))))	     


;;;; Runge-Kutta Integration

(defun runge-kutta (funct xn yn h)
  (let* ((k1 (* h (funcall funct    xn               yn            )))
	 (k2 (* h (funcall funct (+ xn (* 0.5 h)) (+ yn (* 0.5 k1)))))
	 (k3 (* h (funcall funct (+ xn (* 0.5 h)) (+ yn (* 0.5 k2)))))
	 (k4 (* h (funcall funct (+ xn        h ) (+ yn        k1 )))))
    (+ yn (/ (+ k1 (* 2.0 k2) (* 2.0 k3) k4) 6.0))))

(defun integrate-runge-kutta (funct x0 y0 h x1)
  (do ((x x0 (+ x h))
       (y y0))
      ((<= (- x1 x) h)
       (runge-kutta funct x y (- x1 x)))
    (setq y (runge-kutta funct x y h))))


;;;; Predictor Corrector Methods

;;; From Hamming

;;; This code is crying for closures!

;;; Correctors

;;; y(n+1) = a0 y(n) + a1 y(n-1) + a2 y(n-2)
;;;           + h (b[-1] y(n+1) + b0 y'(n) + b1 y'(n-1) + b2 y'(n-2))
;;;           + E5 h**5 (d5/dx5)y / 5!

;;; letting a1 and a2 be parameters gives

(defun correctors (corr-a1 corr-a2)
  (let ((corr-a0  (- 1.0 corr-a1 corr-a2))
	(corr-b-1 (/ (+   9.0 (*   1.0 corr-a1) (*  0.0 corr-a2)) 24.0))
	(corr-b0  (/ (+  19.0 (*  13.0 corr-a1) (*  8.0 corr-a2)) 24.0))
	(corr-b1  (/ (+  -5.0 (*  13.0 corr-a1) (* 32.0 corr-a2)) 24.0))
	(corr-b2  (/ (+   1.0 (*   1.0 corr-a1) (*  8.0 corr-a2)) 24.0))
	(corr-E5  (/ (+ -19.0 (*  11.0 corr-a1) (* -8.0 corr-a2))  6.0)))
    (list corr-E5
	  `(+ (* ,corr-a0 (aref  y-array 0))
	      (* ,corr-a1 (aref  y-array 1))
	      (* ,corr-a2 (aref  y-array 2))
	      (* H (+ (* ,corr-b-1 dy-new)
		      (* ,corr-b0  (aref dy-array 0))
		      (* ,corr-b1  (aref dy-array 1))
		      (* ,corr-b2  (aref dy-array 2))))))))


;;; Adams-Bashforth Type Predictors

;;; y(n+1) = A0 y(n) + A1 y(n-1) + A2 y(n-2)
;;;           + h (B0 y'(n) + B1 y'(n-1) + B2 y'(n-2) + B3 y'(n-3))
;;;           + EE5 h**5 (d5/dx5)y / 5!

;;; letting A1 and A2 be parameters gives

(defun predictors (pred-A1 pred-A2)
  (let ((pred-A0 (- 1.0 pred-A1 pred-A2))
	(pred-B0 (/ (+  55.0 (*   9.0 pred-A1) (*  8.0 pred-A2)) 24.0))
	(pred-B1 (/ (+ -59.0 (*  19.0 pred-A1) (* 32.0 pred-A2)) 24.0))
	(pred-B2 (/ (+  37.0 (*  -5.0 pred-A1) (*  8.0 pred-A2)) 24.0))
	(pred-B3 (/ (+  -9.0 (*   1.0 pred-A1) (*  0.0 pred-A2)) 24.0))
	(pred-E5 (/ (+ 251.0 (* -19.0 pred-A1) (* -8.0 pred-A2))  6.0)))
    (list pred-E5
	  `(+ (* ,pred-a0 (aref  y-array 0))
	      (* ,pred-a1 (aref  y-array 1))
	      (* ,pred-a2 (aref  y-array 2))
	      (* H (+ (* ,pred-b0 (aref dy-array 0))
		      (* ,pred-b1 (aref dy-array 1))
		      (* ,pred-b2 (aref dy-array 2))
		      (* ,pred-b3 (aref dy-array 3))))))))

;;; supposedly, different values of a1 and a2 give different predictors
;;;  (hamming, page 404)
;;; 	a1	a2	type
;;;	0	0	Adams-Bashforth -- *** doesn't agree with Hamming
;;;	1	0	High accuracy
;;;	0	0	Three-eights Rule
;;;	1/3	1/3	"1/3"
;;;	1/2	0	"1/2"
;;;	2/3	0	"2/3"


;;; Modification

;;; Assume the 5th derivative is constant.
;;; The error in (predictor - corrector) is (pred-E5 - corr-E5) h**5 (d5/dx5)y / 5!
;;; The correctors portion of that error is corr-E5/(pred-E5 - corr-E5)
;;; Thus the modified value should be 
;;;      modified = corrector + (predictor - corrector) corr-E5/(pred-E5 - corr-E5)


;;;; Multistep Predictor-Corrector Methods

;;; MULTI-STEP-START
;;;
;;; Before a multi-step method can start, it must have several initial values
;;; Here Runge-Kutta is used to find thoses initial values
;;;
(defun multi-step-start (order funct h x0 y0 x-array y-array dy-array)
  (declare (type fixnum   order)
	   (type function funct)
	   (type float    h x0 y0)
	   (type (array float (*)) x-array y-array dy-array))

  ;; store initial values
  (setf (aref  x-array order) x0)
  (setf (aref  y-array order) y0)
  (setf (aref dy-array order) (funcall funct x0 y0))

  ;; take enough Runge-Kutta Steps to fill out the order
  (do ((i (1- order) (1- i)))
      ((< i 0))
    (declare (fixnum i))
    (setf (aref  x-array i) (+ H (aref x-array (1+ i))))
    (setf (aref  y-array i) (RUNGE-KUTTA FUNCT
					 (aref x-array (1+ i))
					 (aref y-array (1+ i))
					 H))
    (setf (aref dy-array i)
	  (funcall funct (aref x-array i) (aref y-array i)))
    ))


;;;; Multi-Step-Interpolate

;;; Halve the previous step size and interpolate new Y values

(defun multi-step-interpolate (order funct h
			       x-array y-array dy-array
			       x-next  y-next  dy-next)
  (declare (type fixnum order)
	   (type function funct)
	   (type float    h)
	   (type (array float (*))
		 x-array y-array dy-array
		 x-next  y-next  dy-next))
  (do ((i 0 (1+ i)))
      ((> i order))
    (declare (fixnum i))
    (cond ((oddp i)
	   ;; must interpolate this value
	   ;; To halve the interval
	   ;;   already have y(n), y(n-1), y(n-2), y'(n), y'(n-1), y'(n-2)
	   ;;   interpolate exactly for up to x**5
	   (setf (aref  x-next i)
		 (/ (+ (aref x-array     (floor i 2))
		       (aref x-array (1+ (floor i 2))))
		    2.0))
	   (setf (aref  y-next i)
		 (cond ((= i 1)
			(/ (+ (* 45.0 (aref y-array 0))
			      (* 72.0 (aref y-array 1))
			      (* 11.0 (aref y-array 2))
			      (* H
				 (+ (* -9.0 (aref dy-array 0))
				    (* 36.0 (aref dy-array 1))
				    (*  3.0 (aref dy-array 2)))))
			   128.0))
		       ((= i 3)
			(/ (+ (* 45.0 (aref y-array 2))
			      (* 72.0 (aref y-array 1))
			      (* 11.0 (aref y-array 0))
			      (* (- H)
				 (+ (* -9.0 (aref dy-array 2))
				    (* 36.0 (aref dy-array 1))
				    (*  3.0 (aref dy-array 0)))))
			   128.0)			)
		       (t (break "No interpolation F(X) in MULTI-STEP-INTERPOLATE"))))
	   (setf (aref dy-next i)
		 (funcall funct
			  (aref x-next i)
			  (aref y-next i))))
	  (t
	   ;; just copy these values
	   (setf (aref  x-next i) (aref  x-array (floor i 2)))
	   (setf (aref  y-next i) (aref  y-array (floor i 2)))
	   (setf (aref dy-next i) (aref dy-array (floor i 2)))
	   ))))


;;;; Multi-Step Predictor and Corrector Kludges

;;; do this until we have the general case...

(defun multi-step-predict (order h y-array dy-array)
  (declare (type fixnum order)
	   (type float  h)
	   (type (array float (*)) y-array dy-array))
  (if (not (= order 3)) (error "Bad order in MULTI-STEP"))
  (/ (+ (*  48.0 (aref y-array 1))
	(*  24.0 (aref y-array 2))
	(* H (+ (*  191.0 (aref dy-array 0))
		(* -107.0 (aref dy-array 1))
		(*  109.0 (aref dy-array 2))
		(*  -25.0 (aref dy-array 3)))))
     72.0))

#|+ignore
(declare (float (multi-step-correct fixnum float t t float)))
|#

(defun multi-step-correct (order h y-array dy-array dy-modify)
  (declare (type fixnum order)
	   (type float  h  dy-modify)
	   (type (array float (*)) y-array dy-array))

  (if (not (= order 3)) (error "Bad order in MULTI-STEP"))

  (/ (+ (*  48.0 (aref y-array 1))
	(*  24.0 (aref y-array 2))
	(* H (+ (*   25.0 dy-modify)
		(*   91.0 (aref dy-array 0))
		(*   43.0 (aref dy-array 1))
		(*    9.0 (aref dy-array 2)))))
     72.0))


;;;; Multi-Step

;;; See Hamming, page 407
;;; *** should split this so user supplies his own predictor and corrector
;;;     as  a function

(defun integrate-multi-step-h (FUNCT X0 Y0 H X1
		   &optional (C-DOUBLE 1.0E-8) (C-HALVE  1.0E-5))
  (declare (type function funct)
	   (type float    z0 y0 h x1 c-double c-halve))

  (LET* ((ORDER       3)			; means uses y(n-0), y(n-1), ..., y(n-order)
	 (NMAX      (* 2 ORDER)))
    (declare (fixnum order nmax))
    (do ((ERROR     0.0)
	 (PREDICT   0.0)
	 (MODIFY    0.0)
	 (DY-MODIFY 0.0)
	 (CORRECT   0.0)
	 (P-C       0.0)
	 (LAST-P-C  0.0 P-C)			; this is a crock if the step changes?
	 (NVALUES     0)
	 (INITIAL   NIL NIL)
	 ( X-NEW    X0)
	 ( Y-NEW    Y0)
	 (DY-NEW    0.0)
	 ( X-ARRAY  (make-array (1+ (* 2 ORDER)) :element-type 'FLOAT)  X-NEXT)
	 ( Y-ARRAY  (make-array (1+ (* 2 ORDER)) :element-type 'FLOAT)  Y-NEXT)
	 (DY-ARRAY  (make-array (1+ (* 2 ORDER)) :element-type 'FLOAT) DY-NEXT)
	 ( X-NEXT   (make-array (1+ (* 2 ORDER)) :element-type 'FLOAT)  X-ARRAY)
	 ( Y-NEXT   (make-array (1+ (* 2 ORDER)) :element-type 'FLOAT)  Y-ARRAY)
	 (DY-NEXT   (make-array (1+ (* 2 ORDER)) :element-type 'FLOAT) DY-ARRAY))
	((> X-NEW X1)				; *** this is wrong -- step may not be good
	 (RUNGE-KUTTA FUNCT X-NEW Y-NEW (- X1 X-NEW)))
      (declare (type float  error predict modify dy-modify correct p-c last-p-c)
	       (type fixnum nvalues)
	       (type t      initial)
	       (type float  x-new y-new dy-new)
	       (type (array float (*)) x-array y-array dy-array)
	       (type (array float (*)) x-next  y-next  dy-next))

      ;; if we need to start, then start!
      (cond ((= NVALUES 0)
	     (multi-step-start order funct h x0 y0 x-array y-array dy-array)
	     (setq last-p-c 0.0)
	     (setq initial T)
	     (setq NVALUES ORDER)))

      ;; predict ***
      (setq x-new   (+ (aref x-array 0) h))
      (setq predict (multi-step-predict order h y-array dy-array))

      ;; modify ***
      (setq modify (- predict (* (/ 707.0 750.0) last-p-c)))
      (setq dy-modify (funcall funct x-new modify))

      ;; correct ***
      (setq correct (multi-step-correct order h y-array dy-array dy-modify))
      (setq p-c (- predict correct))
      (setq error (abs p-c))
      (setq  y-new (+ correct (* (/ 43.0 750.0) p-c)))
      ;; Hamming says this isn't necessary
      ;;(setq dy-new (funcall funct x-new y-new))
      (setq dy-new dy-modify)

      (cond ((> ERROR C-HALVE)			; Halve Step if necessary
	     (cond (initial
		    ;; if initial step, the start has to be tried again
		    (setq nvalues 0))
		   (t
		    (multi-step-interpolate order funct h
					    x-array y-array dy-array
					    x-next  y-next  dy-next)))
	     (setq H (* 0.5 H)))
	    ((AND (<  ERROR C-DOUBLE)
		  (>= NVALUES (1- (* 2 ORDER))))
	     ;; Double interval if there are enough back values
	     ;; y( 1) --> y( 0)
	     ;; y(-1) --> y(-1)
	     ;; y(-3) --> y(-2)
	     (setf (aref  x-next 0)  X-NEW)
	     (setf (aref  y-next 0)  Y-NEW)
	     (setf (aref dy-next 0) DY-NEW)
	     (do ((i 1 (1+ i)))
		 ((> i order))
	       (declare (fixnum i))
	       (setf (aref  X-NEXT i) (aref  X-array (1- (* 2 i))))
	       (setf (aref  Y-NEXT i) (aref  Y-array (1- (* 2 i))))
	       (setf (aref DY-NEXT i) (aref DY-array (1- (* 2 i)))))
	     (setq nvalues order)
	     (setq H (* 2.0 H)))
	    (t
	     ;; maintain same step
	     (setf (aref  x-next 0)  X-NEW)
	     (setf (aref  y-next 0)  Y-NEW)
	     (setf (aref dy-next 0) DY-NEW)
	     (do ((i 1 (1+ i)))
		 ((>= i nmax))
	       (setf (aref  x-next i) (aref  x-array (1- i)))
	       (setf (aref  y-next i) (aref  y-array (1- i)))
	       (setf (aref dy-next i) (aref dy-array (1- i))))
	     (setq nvalues (1+ nvalues))
	     ))
      )))


;;;; Testing Functions
#|
(EVAL-WHEN (EVAL)
|#
#|
  (DEFUN TEST-FCN (X Y)
    (- 1.0 Y))
|#
#|
  (DEFUN RK-TEST (STEP)
    (FORMAT T "~%	time	RUNGE	EULER	EXACT")
    (DO ((TIME 0.0 (+ TIME STEP))
	 (RK   0.0 (RUNGE-KUTTA (FUNCTION TEST-FCN) TIME RK STEP))
	 (EU   0.0 (EULER       (FUNCTION TEST-FCN) TIME EU STEP)))
	((> TIME 10.0))
      (FORMAT T "~%	~4F	~5F	~5F	~5F"
	      TIME RK EU (- 1.0 (EXP (- TIME))))))
|#
  ;; Test for MULTI-STEP
#|
  (DEFUN TEST ()
    (INTEGRATE-MULTI-STEP-H #'(LAMBDA (X Y) (1+ (* Y Y)))
			    0.0
			    0.0
			    0.2
			    1.000001))
  (DEFUN TEST1 ()
    (INTEGRATE-MULTI-STEP-H #'(LAMBDA (X Y) 1.0)
			    0.0
			    0.0
			    0.00001
			    1.000001))

  )
|#
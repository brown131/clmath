;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Bessel Functions

;;; Reference
;;;   M. Abramowitz and I. Stegun, eds,
;;;     Handbook of Mathematical Functions,
;;;     National Bureau of Standards, 1964

;;;  (c) Copyright Gerald Roylance 1982, 1984
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   Do Y0, Y1 sometime
;;;   Use swamped Miller recurrance relations instead

(in-package "CLMATH")

;;;; BESSEL FUNCTIONS

;;; The Bessel functions J0(x) J1(x) . . .

(defun BESSEL-J (order x)
  (declare (float z2))
  (cond ((= order 0)				;NBS 9.4.1 eps < 5e-8
	 (cond ((<= (abs x) 3.0)
		(let ((z2 (expt (/ x 3.0) 2)))
		  (poly z2
			1.0000000  -2.2499997
			1.2656208  -0.3163866
			0.0444479  -0.0039444
			0.0002100)))
	       ((>= x 3.0)			;NBS 9.4.3
		(let* ((z      (/ 3.0 x))
		       (f0     (poly z		; eps 1.6E-8
				     0.79788456 -.00000077 -.00552740
				     -.00009512 0.00137237 -.00072805
				     0.00014476))
		       (theta0 (+ x
				  (poly z	; eps 7.0E-8
					-.78539816 -.04166397
					-.00003954 0.00262537
					-.00054125 -.00029333
					0.00013558))))
		  (declare (float z f0 theta0))
		  ;; Y0 is ... (/$ (*$ f0 (sin theta0)) (sqrt x))
		  (/ (* f0 (cos theta0)) (sqrt x))))
	       (t
		(BESSEL-J ORDER (- X)))))
	((= order 1)
	 (cond ((<= (abs x) 3.0)
		(let ((z2 (expt (/ x 3.0) 2)))
		  (* (* x
			(poly z2		;eps < x * 1.3e-8
			      0.50000000 -0.56249985
			      0.21093573 -0.03954289
			      0.00443319 -0.00031761
			      0.00001109)
			))))
	       ((>= x 3.0)			;NBS 9.4.3
		(let* ((z      (/ 3.0 x))
		       (f1     (poly z		; eps 4.0E-8
				     0.79788456 0.00000156 0.01659667
				     0.00017105 -.00249511 0.00113653
				     -.00020033))
		       (theta1 (+ x
				  (poly z	; eps 9.0E-8
					-2.35619449 +0.12499612
					+0.00005650 -0.00637879
					+0.00074348 +0.00079824
					-0.00029166))))
		  (declare (float z f1 theta1))
		  ;; Y1 is ... (/$ (*$ f1 (sin theta1)) (sqrt x))
		  (/ (* f1 (cos theta1)) (sqrt x))))
	       (t
		(- (bessel-j order (- x))))))
	((< order 0)
	 (setq order (- order))
	 (cond ((oddp order)
		(- (BESSEL-J order x)))
	       (t (BESSEL-J order x))))
	(t (- (* (/ (float (* 2 (1- order))) x)
		 (BESSEL-J (1- order) x))
	      (BESSEL-J (- order 2) x)))))


;;;; MODIFIED BESSEL FUNCTIONS

;;; modified Bessel functions I0(x) I1(x) ...

;;; What ever happened to the backward recurrence formula???

;;; (exp z cos theta) = I0(z) + 2 (sigma (k = 1 inf) Ik(z) cos (k theta))

;;; I[v+1] = I[v-1] - (2 v / z) I[v]

;;; NBS 9.8.1

(defun BESSEL-I (order x)			
  (declare (float x z2 zi))
  (setq order (abs order))
  (cond ((= order 0)
	 (cond ((< (abs x) 3.75)
		(let ((z2 (expt (/ x 3.75) 2)))	;-3.75 <= x <= 3.75
		  (poly z2			;eps < 1.6e-7
			1.0000000 3.5156229
			3.0899424 1.2067492
			0.2659732 0.0360768
			0.0045813)))
	       ((>= x 3.75)			; 3.75 <= x
		(let ((zi (/ 3.75 x)))
		  (* (/ (exp x) (sqrt x))
		     (poly zi
			   0.39894228 0.01328592
			   0.00225319 -.00157565
			   0.00916281 -.02057706
			   0.02635537 -.01647633
			   0.00392377)))	;eps < 1.9e-7 (exp x) / (sqrt x)
		)
	       (T (BESSEL-I ORDER (- X)))))
	((= order 1)
	 (cond ((< (abs x) 3.75)
		(let ((z2 (expt (/ x 3.75) 2)))	;-3.75 <= x <= 3.75
		  (* x				;eps < x * 8e-9
		     (poly z2
			   0.50000000 0.87890594
			   0.51498869 0.15084934
			   0.02658733 0.00301532
			   0.00032411))))
	       ((>= x 3.75)			; 3.75 <= x
		(let ((zi (/ 3.75 x)))
		  (* (/ (exp x) (sqrt x))
		     (poly zi
			   0.39894228 -.03988024
			   -.00362018 0.00163801
			   -.01031555 0.02282967
			   -.02895312 0.01787654
			   -.00420059)))	;eps < 2.2e-7 (exp x) / (sqrt x)
		)
	       (T (- (BESSEL-I ORDER (- X))))))
	(t (- (BESSEL-I (- order 2) x)		; *** doubly recursive
	      (* (/ (float (* 2 (1- order))) x)
		 (BESSEL-I (1- order) x))))))

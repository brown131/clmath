;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Gamma Function

;;; Reference
;;;   M. Abramowitz and I. Stegun, eds,
;;;     Handbook of Mathematical Functions,
;;;     National Bureau of Standards, 1964
;;;   William H. Press, Brian P. Flannery, Saul A. Teukolsky, William T. Vetterling
;;;     Numerical Recipes in C, The Art of Scientific Computing
;;;     Cambridge University Press, 1988
;;;     cf 178

;;;  (c) Copyright Gerald Roylance 1982, 1984, 1985, 1989
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and fixes
;;;   gamma-incomplete chokes on large x

(in-package "CLMATH")


;;;; NBS Approximation to the Gamma Function

;;; NBS 6.1.15 Recurrence Formula
;;;    (gamma 1+z) = z (gamma z)
;;;
;;; NBS 6.1.17 Reflection Formula
;;;    (gamma 1-z) (gamma z) = -z (gamma -z) (gamma z) = pi csc(pi z)
;;;   therefore
;;;    (gamma z) = (pi csc(pi z)) / (-z (gamma -z))

;;; NBS 6.1.35
;;;
(defun gamma-function-nbs (a)
  (let ((x (- a 1.0)))
    (declare (float x))
    (cond ((< x 0.0) (/ (gamma-function-nbs (1+ a)) a))
	  ((> x 1.0) (* x (gamma-function-nbs x)))
	  (t  (poly x
		    1.000000000 -.577191652
		    0.988205891 -.897056937
		    0.918206857 -.756704078
		    0.482199394 -.193527818
		    0.035868343)		;eps < 3e-7
	      ))))

;;; (POLY X
;;;       1.0000000 -.5748646
;;;       0.9512363 -.6998588
;;;       0.4245549 -.1010678)		;eps < 5e-5

#|+ignore
(- (gamma-function-nbs 1.6) 0.8935153493)
|#

;;;; Stirling's Approximation to the Gamma Function

;;; NBS 6.1.37, Stirling's Formula, an approximation to the gamma-function
;;;
;;; complex z, z->inf, |arg z|<pi

(defun gamma-stirling (z)
  (let ((zi (/ z)))
    (declare (float zi))
    (* (exp (- z))
       (expt z (- z 0.5))
       (sqrt (* 2.0 (float pi z)))
       (poly zi
	     1.0
	     (/    1.0      12.0)
	     (/    1.0     288.0)
	     (/ -139.0   51840.0)
	     (/ -571.0 2488320.0)
	     ;; ...
	     ))))


;;;; Lancoz Approximation to the Gamma Function

;;; from Press...

;;; Lancoz approximation for Re[z] > 0
;;; Gamma(z+1) = (z+gamma+0.5)^{z + 0.5} exp{-(z+gamma+0.5)}
;;;              * sqrt{2 \pi}[c_0 + {c_1 \over z+1} + ... + {c_N \over z+N} + eps]
;;;
;;; for gamma=5, N=6, and a certain set of c_i, abs(eps)< 2.0E-10
;;;  (even for complex arguments!)

(defun gammln (xx)
  (if (< (realpart xx) 1.0)
      (let ((z (- 1.0 xx)))
	(- (log (/ (* pi z) (sin (* pi z))))
	   (gammln (+ 1.0 z))))
      (let* ((z   (- xx 1.0))
	     (tmp (+ z 5.0 0.5)))
	(+ (* (log tmp) (+ z 0.5))
	   (- tmp)
	   (log (sqrt (* 2 pi)))
	   (log (+ 1.0
		   (/  76.18009173d0   (+ z 1.0d0))
		   (/ -86.50532033d0   (+ z 2.0d0))
		   (/  24.01409822d0   (+ z 3.0d0))
		   (/  -1.231739516d0  (+ z 4.0d0))
		   (/   0.120858003d-2 (+ z 5.0d0))
		   (/  -0.536382d-5    (+ z 6.0d0))))))))

(defun fgammln (xx)
  (declare (type single-float xx))
  (if (< (realpart xx) 1.0)
      (let ((z (- 1.0 xx)))
	(- (log (/ (* (float pi 1.0) z)
		   (sin (* (float pi 1.0) z))))	; ***
	   (fgammln (+ 1.0 z))))
      (let* ((z   (- xx 1.0))
	     (tmp (+ z 5.0 0.5)))
	(+ (* (log tmp) (+ z 0.5))
	   (- tmp)
	   (float (log (sqrt (* 2 pi))) 1.0)
	   (log (+ 1.0
		   (/  76.18009173e0   (+ z 1.0))
		   (/ -86.50532033e0   (+ z 2.0))
		   (/  24.01409822e0   (+ z 3.0))
		   (/  -1.231739516e0  (+ z 4.0))
		   (/   0.120858003e-2 (+ z 5.0))
		   (/  -0.536382e-5    (+ z 6.0))))))))

#|

(exp (fgammln 12.0))
(gammln 2.0)
(exp (gammln 6.0))

(round (exp (gammln  3.0d0)))
(round (exp (gammln  6.0d0)))
(round (exp (gammln 10.0d0)))

(list (gammln 1.500d0) -0.1207822376)
(list (gammln 1.825d0) -0.0637301353)

(list (gammln -3.2   ) (log (gamma-function-nbs -3.2)))

|#


;;;; Gamma Function and Log Gamma Function

(defun gamma-function (a)
  (typecase a
    (integer      (exp (fgammln (float a))))
    (single-float (exp (fgammln        a )))
    (otherwise    (exp  (gammln        a )))))

(defun log-gamma-function (a)
  (typecase a
    (integer            (fgammln (float a)))
    (single-float       (fgammln        a ))
    (otherwise           (gammln        a ))))


;;;; Reciprocal of the gamma-function

;;; NBS 6.1.34
;;;   claims |z| < inf, but don't specify accuracy
;;;     I restrict |z| < 1
;;;
;;; from gamma recurrence relations,
;;;   (gfi 1+z) = (gfi z) / z
;;;   (gfi z) = (gfi z-1) / z-1
;;;   (gfi z) = z (gfi 1+z)
;;;
(proclaim '(ftype (function (float) float) gamma-function-reciprocal))

(defun gamma-function-reciprocal (z)
  (cond ((< z -1.0) (* z (gamma-function-reciprocal (1+ z))))
	((> z  1.0) (/   (gamma-function-reciprocal (1- z)) (1- z)))
	(t
	 (poly z
	       +0.0				; 0 -- Numbers have not been checked
	       +1.0000000000000000
	       +0.5772156649015329
	       -0.6558780715202538
	       -0.0420026350340952
	       +0.1665386113822915		; 5
	       -0.0421977345555443
	       -0.0096219715278770
	       +0.0072189432466630
	       -0.0011651675918591
	       -0.0002152416741149		;10
	       +0.0001280502823882
	       -0.0000201348547807
	       -0.0000012504934821
	       +0.0000011330272320
	       -0.0000002056338417		;15
	       +0.0000000061160950
	       +0.0000000050020075
	       -0.0000000011812746
	       +0.0000000001043427
	       +0.0000000000077823		;20
	       -0.0000000000036968
	       +0.0000000000005100
	       -0.0000000000000206
	       -0.0000000000000054
	       +0.0000000000000014		;25
	       +0.0000000000000001		;26
	       ))))


;;;; Incomplete Gamma-Function

;;; *** DOES NOT WORK -- gamma* not total & does not produce correct answers
;;; *** GAMMA-INCOMPLETE bombs on large x

;;; NBS 6.5.4
;;;   \gamma^*(a,x) is a single valued analytic function with no finite singularities
;;;   \gamma^*(a,x) = { x^{-a} \over \Gamma(a) } \gamma(a,x)
;;; NBS 6.5.14
;;;   (gamma* -n x) = (expt x n)
;;; NBS 6.5.29
;;;   \gamma^*(a,x) = (*$ (/$ 1.0 (gamma a))
;;;                      (sigma (n 0 inf) (/$ (expt (- z) n) (+ a n) n!)))
;;;
#|
(proclaim '(ftype (function (float float) float) gamma* gamma*aux))
|#
;;; gamma*aux(a,x) = \Gamma(a) \gamma(a,x)
#|
(defun gamma*aux (a z)				; *** GROSS ROUNDOFF ERRORS
  ;; a > 0 if a is an integer
  (do ((n    0.0 (1+ n))
       (term 1.0
	     (/ (* term (- z))
		(+ n 1.0)))
       (sum  0.0))
      ((progn (setq sum (+ sum (/ term (+ a n))))
	      (< (abs (/ term (+ a n))) 1.0E-7))
       sum)
    (declare (float n term sum))
    ))

(defun gamma* (a z)
  (/ (gamma*aux a z) (gamma-function a)))
|#
;;; NBS 6.5.4 implies
;;;   \gamma(a,x) = {\Gamma(a) \over x^{-a}} \gamma^*(a,x)
;;;   \gamma(a,x) =  \Gamma(a)       x^{ a}  \gamma^*(a,x)
;;; NBS 6.5.22
;;;   \gamma(a+1,x) = a \gamma(a,x) - x^{a} e^{-x}
;;;     \gamma(a,x) = (\gamma(a+1,x)  + x^{a} e^{-x}) / a
#|
(defun gamma-function-incomplete (a x)
  (if (< a 1.0)
      (/ (+ (gamma-function-incomplete (+ a 1.0) x)
	    (* (expt x a) (exp (- x))))
	 a)
      (* (expt x a) (gamma*aux a x))))
|#
;;; NBS 6.5.31
;;;	x > 0      |a| < inf
;;; (GAMMA a x) = (exp -x) * (expt x a) *
;;;      (continued-fraction (1 / x+)
;;;			     (1-a / 1+) (1 / x+)
;;;			     (2-a / 1+) (2 / x+)
;;;			     . . .)
;;; this GAMMA is really a tail
;;; (GAMMA a x) = (gamma-function a) - (gamma a x)
;;;		= (integral x to inf of (exp -t)(expt t a-1) dt)
;;; (gamma a+1 x) = a (gamma a x) - (exp -x) (expt x a)
;;; note distinction of LARGE and small gammas



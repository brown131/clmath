;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Beta Functions

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
;;;   

(in-package "CLMATH")

;;;; Beta Function

;;; The beta function is B(a,b)
;;;    B(a,b) = B (b,a)

(defun beta-function-naive (m n)		;m,n>0
  (/ (* (gamma-function m)
	(gamma-function n))
     (gamma-function (+ m n))))

(defun log-beta-function (m n)
  (- (+ (log-gamma-function m) (log-gamma-function n))
     (log-gamma-function (+ m n))))

(defun beta-function (m n)			;m,n>0
  (exp (log-beta-function m n)))


;;;; Incomplete Beta Function

;;; The incomplete beta function is B[x](a,b)
;;; I[x](a,b) = B[x](a,b) / B(a,b)
;;;
;;; NBS 6.6.3
;;;    I[x](a,b) = 1 - I[1-x](b,a)
;;;  therefore
;;;    B[x](a,b) = B(a,b) - B[1-x](b,a)

;;; NBS 26.5.1
;;;   I_x(a,b) = (1/\beta(a,b)) \int_0^x t^{a-1} (1-t)^{b-1} dt

;;; NBS 26.5.2
;;;   I_x(a,b) = 1 - I_{1-x}(b,a)
;;;   therefore
;;;     \beta(a,b) I_x(a,b) = \beta(a,b) - \beta(a,b)I_{1-x}(b,a)

(in-package "CLMATH")

;;; NBS 26.5.7

;;; I don't know where this approximation came from
;;;   either m or n must be an integer
;;;
(defvar beta-warning nil)
(defun beta-function-incomplete (m n x)
  (setq m (float m 1.0d0))			; *** Use a sledge hammer
  (setq n (float n 1.0d0))
  (setq x (float x 1.0d0))
  (if (not beta-warning)
      (progn (warn "*** BETA-FUNCTION-INCOMPLETE is noisy")
	     (setq beta-warning t)))
  (cond ((= n (float (floor n)))
	 (do ((i    1.0d0 (1+ i))
	      (beta 0.0d0)
	      (Ti   (/ (expt x m) m)))
	     ((> i n) (float beta 1.0))
	   (declare (type double-float i beta Ti))
	   (setq beta (+ beta Ti)
		 Ti   (- (/ (* (+ m i -1.0d0) (- n i) x Ti)
			    (* (+ m i) i))))))
	((= m (float (floor m)))
	 (- (beta-function m n)
	    (beta-function-incomplete n m (- 1.0d0 x))))
	(t
	 (error "BETA-FUNCTION-INCOMPLETE must have one integer parameter"))))

#|

*** this is incredibly noisy

(module-require "PLOT")

(- (beta-function-incomplete 3.0 10.0 .7333333)
   (beta-function-incomplete 3.0 10.0 .7166667))

(defun foo (a b)
  (plot-fcn  100 #'(lambda (x) (beta-function-incomplete a b x))
	    0.0 1.0 pstrm))

(defun bar (a b)
  (plot-fcn 1000 #'(lambda (x) (beta-function-incomplete a b x))
	    0.6 1.0 pstrm))

(defun bard (a b)
  (plot-fcn 1000 #'(lambda (x) (beta-function-incomplete a b x))
	    0.6d0 1.0d0 pstrm))

(bar  3.0 10.0)
(bard 3.0 10.0)

|#


;;;; Incomplete Beta Function I_x

;;; reference
;;;   William H. Press, Brian P. Flannery, Saul A. Teukolsky, William T. Vetterling
;;;     Numerical Recipes in C, The Art of Scientific Computing
;;;     Cambridge University Press, 1988
;;;     cf 178

;;; I_x(a,b) = {x^a (1-x)^b \over a B(a,b)}
;;;              [{1\over 1+}{d_1\over 1+}{d_2\over 1+}...]
;;; d_{2m+1} = - {(a+m)(a+b+m)x \over (a+2m)(a+2m+1)}
;;; d_{2m  } =   {m(b-m)x \over (a+2m-1)(a+2m)}

;;; iterations are O(sqrt(max(a,b))) for x < (a+1)/(a+b+1)

;;; symmetry
;;;   I_x(a,b) = 1 - I_{1-x}(b,a)

(defun betacf (a b x)
  (declare (float a b x))
  (let* ((eps 3.0e-7)
	 (qab (+ a b))
	 (qap (+ a 1.0))
	 (qam (- a 1.0))
	 (bz  (- 1.0 (/ (* qab x) qap)))
	 (bm  1.0)
	 (az  1.0)
	 (am  1.0))
    (declare (float eps qab qap qam bz bm az am))
    (do ((m 1 (1+ m)))
	((>= m 100)
	 (error "too many iterations"))
      (declare (fixnum m))

      (let* ((em  (float m))
	     (tem (* 2.0 em))
	     (d   (/ (* em (- b em) x)
		     (* (+ qam tem) (+ a tem))))
	     (ap  (+ az (* d am)))		; even step of recurrance
	     (bp  (+ bz (* d bm)))
	     (d   (- (/ (* (+ a em) (+ qab em) x)
			(* (+ qap tem) (+ a tem)))))
	     (app (+ ap (* d az)))		; odd step of recurrance
	     (bpp (+ bp (* d bz))))
	(declare (float em tem d ap bp app bpp))

	(let ((aold az))			; save old answer
	  (declare (float aold))
	  ;; renormalize
	  (setq am (/ ap  bpp))
	  (setq bm (/ bp  bpp))
	  (setq az (/ app bpp))
	  (setq bz 1.0)
	  (if (< (abs (- az aold)) (* EPS (abs az)))
	      (return az))
	  )))
    ))


;;;; Incomplete Beta Function

(defun betai (a b x)
  (declare (float a b x))
  (if (or (< x 0.0) (> x 1.0))
      (error "BETAI:  Bad argument, x=~f" x)
      (let ((bt (if (or (= x 0.0) (= x 1.0))
		    0.0
		    (exp (+ (- (log-gamma-function (+ a b))
			       (log-gamma-function a)
			       (log-gamma-function b))
			    (* a (log x))
			    (* b (log (- 1.0 x))))))))
	(declare (float bt))
	(if (< x (/ (+ a 1.0) (+ a b 2.0)))	; which converges quickest?
	    (* bt (/ (betacf a b x) a))
	    (- 1.0 (/ (* bt (betacf b a (- 1.0 x))) b))))))

#|

(defun foo (a b)
  (plot-fcn  100 #'(lambda (x) (beta-function-incomplete a b x))
	    0.0 1.0 pstrm))
(defun bar (a b)
  (plot-fcn  100 #'(lambda (x) (betai a b x))
	    0.0 1.0 pstrm))

|#

;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Numerical Integration

;;; Reference
;;;   M. Abramowitz and I. Stegun, eds,
;;;     Handbook of Mathematical Functions,
;;;     National Bureau of Standards, 1964.

;;;  (c) Copyright Gerald Roylance 1982
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   Formulas should be computed instead of given

(in-package "CLMATH")

;;;; Trapezoidal Rule

;;; NBS 25.4.2

;;; - (m h**3 / 12) f''(zeta)

(DEFUN INTEGRATE-TRAPEZOIDAL (F X0 XM M)
  (DECLARE (FLOAT X0 XM SUM H STEP STEPS)
	   (FIXNUM M))
  (DO ((SUM   0.0)
       (H     (/ (- XM X0) (float M)))  
       (STEP  1.0 (1+ STEP))
       (STEPS (float M)))
      ((>= STEP STEPS)
       (* H (+ (* 0.5 (FUNCALL F X0))
		 SUM
		 (* 0.5 (FUNCALL F XM)))))
    (SETQ SUM (+ SUM (FUNCALL F (+ X0 (* STEP H)))))))


;;;; Simpson's Rule

;;; NBS 25.4.6

;;; - (n h**5 / 90) (f''''(zeta))
;;; 2n = m

(DEFUN INTEGRATE-SIMPSON (F X0 XM M)
  (DECLARE (FLOAT X0 XM SUM-ODD SUM-EVEN H STEP)
	   (FIXNUM I M))
  (COND ((ODDP M) (SETQ M (1+ M))))
  (DO ((SUM-EVEN 0.0)
       (SUM-ODD  0.0)
       (H     (/ (- XM X0) (float M)))
       (I     1   (1+ I))
       (STEP  1.0 (1+ STEP)))
      ((>= I M)
       (* (/ H 3.0)
	   (+ (FUNCALL F X0)
	       (* 4.0 SUM-ODD)
	       (* 2.0 SUM-EVEN)
	       (FUNCALL F XM))))
    (COND ((ODDP I)
	   (SETQ SUM-ODD  (+ SUM-ODD  (FUNCALL F (+ X0 (* STEP H))))))
	  (T
	   (SETQ SUM-EVEN (+ SUM-EVEN (FUNCALL F (+ X0 (* STEP H)))))))
    ))


;;;; Simpson's Rule -- Recursive Formulation

;;; keep increasing number of steps until difference is
;;; small enough.

;;; Hack is that it is no more work

(DEFUN INTEGRATE-SIMPSON-RECURSIVE (F X0 XM EPS)
  (DECLARE (FLOAT X0 XM SUM-ODD SUM-EVEN H STEP F0 F1 FM INT INT-LAST)
	   (FIXNUM I M))
  (LET ((F0 (FUNCALL F X0))
	(F1 (FUNCALL F (* 0.5 (+ X0 XM))))
	(FM (FUNCALL F XM))
	(H  (* 0.5 (- XM X0))))

    (DO ((M          2)
	 (SUM-EVEN 0.0)
	 (SUM-ODD   F1)
	 (INT-LAST 0.0)
	 (INT      (* (/ H 3.0) (+ F0 (* 4.0 F1) FM))))
	(NIL)

      (SETQ INT-LAST INT)
      (SETQ M (*   M 2))
      (SETQ H (/ H 2.0))
      (SETQ SUM-EVEN (+ SUM-ODD SUM-EVEN))
      (SETQ SUM-ODD  0.0)
      (DO ((I 1 (+ I 2))
	   (STEP 1.0 (+ STEP 2.0)))
	  ((>= I M))
	(SETQ SUM-ODD (+ SUM-ODD (FUNCALL F (* H STEP)))))

      (SETQ INT (* (/ H 3.0)
		    (+ F0
			(* 4.0 SUM-ODD)
			(* 2.0 SUM-EVEN)
			FM)))

      (COND ((< (ABS (- INT INT-LAST)) EPS)
	     (RETURN INT-LAST)))
      )))


;;;; Tests
#|
(EVAL-WHEN (EVAL)
|#
#|
  (DEFUN TEST0 ()
    (INTEGRATE-TRAPEZOIDAL '(LAMBDA (X) (* 3.0 (expt X 2)))
			   0.0
			   1.0
			   10))
  (DEFUN TEST1 ()
    (INTEGRATE-SIMPSON '(LAMBDA (X) (* 3.0 (expt X 2)))
			   0.0
			   1.0
			   10))
  (DEFUN TEST2 ()
    (FORMAT T "~%ANSWER ~E"  (SIN 1.0))
    (INTEGRATE-SIMPSON-RECURSIVE
     '(LAMBDA (X) (COS X))
     0.0
     1.0
     0.001))
  )
|#
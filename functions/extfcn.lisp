;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Extended Functions

;;; References
;;;   M. Abramowitz and I. Stegun, eds,
;;;     Handbook of Mathematical Functions,
;;;     National Bureau of Standards, 1964
;;;   Digital Equipment Corporation,
;;;     PDP-11 Paper Tape Software Programming Handbook,
;;;     DEC-11-GGPA-D, Section 7.7.
;;;   Cecil Hastings, Jr.,
;;;     Approximations for Digital Computers,
;;;     Princeton University Press, Princeton, NJ, 1955.

;;;  (c) Copyright Gerald Roylance 1982, 1984, 1985
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and fixes
;;;  Some code still untested.
;;;  Routine accuracy is usually around 5.0E-8 (Roundoff Problems?)
;;;    eg -- plotting (sine x) - (sin x)

(in-package "CLMATH")

;;;; Exponential and Logarithm

;;; exponential does only has 24 bits of accuracy or so
;;;  -- 6.0E-8 is the difference between it and (exp x) in [0, 0.5]

;;; (exp x) = (expt 2 (* (log2 e) x))
;;;         = (expt 2 (/$ x (log 2)))
;;;         = (fsc (expt 2 (frac (/$ x (log 2)))) (int (/$ x (log 2))))
;;;         = (fsc (expt 2 (frac (/$ x (log 2))))
;;;                (int (/$ x (log 2))))
;;;         = (fsc (*$ (expt 2 0.5) (expt 2 (-$ (frac (/$ x (log 2))) 0.5)))
;;;                (int (/$ x (log 2))))

(defun exponential-1 (x)
  (if (or (< x 0.0) (> x (* 0.5 (log 2.0))))
      (error "EXPONENTIAL:  out of range"))
  (1+ (/ (* 2.0 x)				;eps < 1.0e-9
	 (+ 12.01501675387500
	    (- x)
	    (/ -601.8042666979565
	       (+ 60.09019073192600
		  (* x x)))))))

(defun exponential (x)
  (let* ((xlog2 (/ x (log 2.0)))
	 (int   (floor xlog2))
	 (frac  (- xlog2 (float int))))
    (declare (fixnum int)
	     (float xlog2 frac))
    (floating-scale (cond ((< frac 0.5)
			   (exponential-1 (* frac (log 2.0))))
			  (t
			   (* (sqrt 2.0)
			      (exponential-1 (* (- frac 0.5) (log 2.0)))
			      )))
		    int)))


;;;; Logarithm

;;; empirically (-$ (log x) (logarithm x)) is 8.0E-9

(defun logarithm (x)
  (if (not (plusp x))
      (error "LOGARITHM:  non positive arg"))
  (let* ((exp (floating-exponent x))
	 (man (floating-mantissa x))
	 (u (/ (- man (/ (sqrt 2.0) 2.0))
	       (+ man (/ (sqrt 2.0) 2.0))))
	 (z  (* u u)))
    (declare (float man u z)
	     (fixnum exp))
    (+ (* (float exp) (log 2.0))
       (- (* u					;eps < 2^-32
	     (poly z
		   1.999999993788
		   0.666669470507
		   0.399659100019
		   0.300974506336))
	  (/ (log 2.0) 2.0)))))


;;;; Trigonometric Functions

;;; from DEC

(defun sine (x)
  (declare (float x sign z z2))
  (let ((sign 1.0)
	(z    (* 4.0 (fraction (/ x (* (float pi 1.0) 2.0)))))
	(z2   0.0))
    (if (> z 2.0) (setq sign -1.0		;sin(pi+z)=-sin(z)
			z (- z 2.0)))		;  z is quad 0 or 1 now
    (if (> z 1.0) (setq z (- 2.0 z)))		;sin(pi:2+z)=sin(pi:2-z)
    (setq z2 (* z z))
    (* sign z (poly z2				;eps < .5e-9
		    1.57079632662143 -0.64596409264401
		    0.07969258728630 -0.00468162023910
		    0.00016021713430 -0.00000341817225))
    ))

(defun cosine (x)
  (sine (- (/ (float pi 1.0) 2.0) x)))


;;;; EXTENDED FUNCTIONS

;;; from DEC
;;; answers are within plus or minus pi

(defun arctan (y x)
  (cond ((< y 0.0) (- (arctan (- y) x)))
	((< x 0.0) (- (float pi 1.0) (arctan y (- x))))
	((= y 0.0) (cond ((= x 0.0) (error "ARCTAN:  0,0 bad arg"))
			 (t 0.0)))
	((= x 0.0) (/ (float pi 1.0) 2.0))
	((> (/ y x) 1.0) (- (/ (float pi 1.0) 2.0) (arctan x y)))
	(t (setq y (/ y x))
	   (let ((c 0.0) (z 0.0) (z2 0.0))
	     (cond ((< y (- 2.0 (sqrt 3.0)))
		    (setq c 0.0
			  z y
			  z2 (* y y)))
		   (t
		    (setq c (/ (float pi 1.0) 6.0)
			  z (/ (1- (* y (sqrt 3.0)))
			       (+ y (sqrt 3.0)))
			  z2 (* z z))))
	     (+ c
		(* z
		   (poly z2
			 0.99999999843 -0.33333289364
			 0.19996534780 -0.14173460613
			 0.09491954952)))))))

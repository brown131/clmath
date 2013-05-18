;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Probability and Statistics

;;; Poisson Distribution

;;;  References
;;;    D. Knuth
;;;      The Art of Computer Programming, Vol 2, second edition
;;;      Addison Wesley, 1981
;;;    A. J. Kinderman & J. F. Monahan
;;;      Computer Generation of Random Variables Using
;;;      the Ratio of Uniform Deviates
;;;      ACM Transactions on Mathematical Software
;;;      Vol 3 No 3 September 1977 pp 257-260
;;;    M. Abramowitz and I. Stegun, eds,
;;;      Handbook of Mathematical Functions,
;;;      National Bureau of Standards, 1964.

;;;  (c) Copyright Gerald Roylance 1983, 1984, 1987
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   

(in-package "CLMATH")

;;;; Poisson Distribution

;;; LAM and TIME are always multiplied together

(defun poisson-density (n lambda)
  (cond ((> lambda 0.0)
	 (/ (* (expt lambda n) (exp (- lambda)))
	    (float (factorial n))))
	((= lambda 0.0)
	 (cond ((= n 0) 1.0)
	       (t       0.0)))
	((< lambda 0.0) (ERROR "Negative lambda to POISSON-DENSITY"))))

#|
Should be fun to optimize this function by procedure integration
and strength reduction.  Do people use recursion optimization?

(defun poisson-cumulative (n time)
  (do ((prob 0.0)
       (i    0  (1+ i)))
      ((> i n) prob)
    (declare (flonum prob)
	     (fixnum i))
    (setq prob (+$ prob (poisson-density i time)))))

|#

(defun poisson-cumulative (n time)
  (do ((prob 0.0)
       (t**i 1.0 (* t**i time))
       (fi   1   (* fi (1+ i)))
       (i    0  (1+ i)))
      ((> i n) (* (exp (- time)) prob))
    (declare (float prob t**i)
	     (fixnum i fi))
    (setq prob (+ prob
		  (/ t**i (float fi))))))


;;;; Random Number Generator -- POISSON

;;; POISSON DISTRIBUTION

;;; *** slow and stupid
;;; *** might bomb if U is close to 1
;;;
(defun poisson-random-number (lmbd)
  (do ((u (uniform-random-number))
       (p 0.0)
       (i 0 (1+ i)))
      ((progn (setq p (+ p (poisson-density i lmbd)))
	      (< u p))
       i)
    (declare (fixnum i)
	     (float u p))
    ))

#|
(eval-when (eval)
|#
#|
  (defun po-test (lam n)
    (do ((i 0 (1+ i))
	 (s 0))
	((>= i n)
	 (/ (float s) (float n)))
      (declare (fixnum i s))
      (setq s (+ s (poisson-random-number lam)))))
  )
|#
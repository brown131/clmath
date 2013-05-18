;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Probability and Statistics

;;; Binomial Distribution

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

;;;; Binomial Distribution (discrete)

;;; Probability that there are k heads out of n trials
;;; where p is the probabily of heads.

;;; \choose(n,k) p^k (1-p)^k

(defun binomial-density-naive (k n p)
  (* (float (dchoose n k) p)
     (expt        p       k )
     (expt (- 1.0 p) (- n k))))

;;; For large arguments, overflow is a problem, so compute with logs.
;;;
(defun binomial-density (k n p)
  (exp (+ (float (dlog-choose n k) p)
	  (*      k  (log        p ))
	  (* (- n k) (log (- 1.0 p))))))

#|

(defun beta-choose (n k)
  (/ 1.0 (beta-function (+ (- n k) 1) (+ k 1)) (+ n 1)))

(defun bin-test (k n p)
  (list (binomial-density-naive k n p)
	(binomial-density       k n p)))

(bin-test 3 9 0.1d-8)

|#


;;;; Binomial Cumulative Distribution

;;; The slow and stupid version
;;;
(defun binomial-cumulative-naive (x n p)
  (if (> x n) (error "BINOMIAL-CUMULATIVE:  bogus arguments"))
  (do ((i 0   (1+ i))
       (s 0.0))
      ((> i x) s)
    (declare (fixnum i)
	     (float s))
    (setq s (+ s (binomial-density i n p)))
    ))

;;; NBS 26.5.24
;;;   \sum_{s=a}^{n} binomial-density(k,s,p) = I_p(a,n-a+1)
;;; warp the expression into \sum_{s=0}^{a}...
;;;   swap p and (1-p)
;;;   complement a
;;;
(defun binomial-cumulative (k n p)
  (let ((a (- n k)))
     (betai a (- n a -1) (- 1 p))))

;;; k or more occurances in n trials of probability p
;;;
(defun binomial-tail (k n p)
  (let ((a k))
    (if (= k 0)
	1.0
	(betai a (- n a -1) p))))

#|

(binomial-density 0 4 0.1)

(binomial-tail 2 9 0.00000000)

(let ((a 3)
      (n 9)
      (p 1.0e-8))
  (list (binomial-cumulative-naive a n p)
	(binomial-cumulative       a n p)))

|#


;;;; Random Number Generator -- BINOMIAL

;;; *** slow and stupid
;;;   -- see Marsaglia for a better way
;;;
(defun binomial-random-number (n p)
  (do ((u (uniform-random-number))
       (c 0.0)
       (i 0 (1+ i)))
      ((progn (setq c (+ c (binomial-density i n p)))
	      (< u c))
       i)
    (declare (float c u)
	     (fixnum i))
    ))

#|

(defun bi-test (n p m)
  (do ((i 0 (1+ i))
       (s 0))
      ((>= i m)
       (/ (float s) (float m)))
    (declare (fixnum i s))
    (setq s (+ s (binomial-random-number n p)))))

|#

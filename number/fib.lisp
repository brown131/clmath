;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Fibonacci Numbers

;;;  (c) Copyright Gerald Roylance 1984
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   

(in-package "CLMATH")

;;;; Fibonacci -- The Slow Recursive Way

;;; f(0) = 0
;;; f(1) = 1
;;; f(n) = f(n-1) + f(n-2)
;;;
(defun fib-slow (n)
  (cond ((= n 0) 0)
	((= n 1) 1)
	(t (+ (fib-slow (- n 1))
	      (fib-slow (- n 2))))))


;;;; Fibonacci -- The Fast Recursive Way

;;; f(n) = 1 f(n  ) + 0 f(n-1)
;;;      = 1 f(n-1) + 1 f(n-2) = [f(n-2) + f(n-3)] + f(n-2)
;;;      = 2 f(n-2) + 1 f(n-3)
;;;      = 3 f(n-3) + 2 f(n-4)
;;;      = 5 f(n-4) + 3 f(n-5)
;;;      = 8 f(n-5) + 5 f(n-6)
;;;
;;; f(n) = f(k+1) f(n-k) + f(k)f(n-k-1)
;;;
;;; f(n+k)   = f(k+1) f(n) + f(k)f(n-1)

(defun fibn+k (fk fk+1 fn fn+1)
  (+ (* fk+1 fn)
     (* fk (- fn+1 fn))))

(defun fibn+k+1 (fk fk+1 fn fn+1)
  (+ (* fk+1 fn+1)
     (* fk fn)))

(defun fib-iter (n k fk fk+1 m fm fm+1)
  (cond ((>= m n) fm)
	((oddp (/ (- n m) k))
	 (fib-iter n k fk fk+1 (+ m k)
		   (fibn+k   fm fm+1 fk fk+1)
		   (fibn+k+1 fm fm+1 fk fk+1)))
	(t
	 (fib-iter n
		   (* 2 k)
		   (fibn+k   fk fk+1 fk fk+1)
		   (fibn+k+1 fk fk+1 fk fk+1)
		   m fm fm+1))))
		   
(defun fib (n)
  (fib-iter n 1 1 1 0 0 1))


;;;; Fibonacci -- Binet's Formula

(defun fib-binet (n)
  (* (/ 1 (sqrt 5))
     (- (expt (/ (+ 1 (sqrt 5)) 2) n)
	(expt (/ (- 1 (sqrt 5)) 2) n))))


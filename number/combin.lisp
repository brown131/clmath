;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Combinatorial Algorithms

;;; reference
;;;   William H. Press, Brian P. Flannery, Saul A. Teukolsky, William T. Vetterling
;;;     Numerical Recipes in C, The Art of Scientific Computing
;;;     Cambridge University Press, 1988
;;;     cf 178

;;;  (c) Copyright Gerald Roylance 1982, 1984, 1987, 1988
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;

(in-package "CLMATH")

;;;; Factorial Function

;;; The number of ways to arrange n objects

;;; The bignum version.

;;; 12! < 2^31

(proclaim '(ftype (function (fixnum) integer) factorial))

(defun factorial (n)
  (declare (fixnum n))
  (do ((i    2 (1+ i))
       (fact 1 (* fact i)))
      ((> i n) fact)
    (declare (fixnum i))
    ))


;;;; The Logarithm of the Factorial

;;; A table to hold precomputed values
;;;
(defvar dlog-factorial-table
	(make-array 100 :initial-element NIL))

;;; A function to precompute the values of the table
;;;   -- not used
;;;
(defun dlog-factorial-table-precompute ()
  (dotimes (i (array-dimension dlog-factorial-table 0))
    (setf (aref dlog-factorial-table i)
	  (if (<= i 1)
	      0.0d0
	      (log-gamma-function (+ (float i 1.0d0) 1.0d0))))))

;;; Double Precision version of the factorial
;;;
(defun dlog-factorial (n)
  (cond ((<  n 0) (error "LOG-FACTORIAL"))
	((<= n 1) 0.0d0)
	((<  n (array-dimension dlog-factorial-table 0))
	 (if (aref dlog-factorial-table n)
	     (aref dlog-factorial-table n)
	     (setf (aref dlog-factorial-table n)
		   (log-gamma-function (+ (float n) 1.0d0)))))
	(t
	 (log-gamma-function (+ (float n) 1.0d0)))))

;;; Single Precision mooches off of the double precision
;;;   (the table will have all/most of the reasonable values)
;;;
(defun log-factorial (n)
  (float (dlog-factorial n) 1.0))


;;;; Floating Point Factorial

(defvar dfactorial-table
	(let* ((n     100)
	       (array (make-array n
				  :element-type    'double-float
				  :initial-element 1.0d0)))
	  (do ((i 1 (1+ i)))
	      ((>= i n) array)
	    (setf (aref array i) (* i (aref array (- i 1)))))))

(defun dfactorial (n)
  (cond ((< n 0) (error "DFACTORIAL"))
	((< n (array-dimension dfactorial-table 0))
	 (aref dfactorial-table n))
	(t
	 (values (round (exp (dlog-factorial n)))))))

(defun ffactorial (n)
  (float (dfactorial n) 1.0))


;;;; A Generalization of Factorial

;;; \Gamma(x+1) = x \Gamma(x)

;;; for k<m,
;;;    prod_{l=0,...,n-1} (k + lm)
;;;    = \Gamma(n+(k/m) m^n / \Gamma(k/m)

(defun gprod (n k m)
  (let ((a (/ (float k) (float m))))
    (values (round (/ (* (gamma-function (+ (float n) a)) (expt m n))
		      (gamma-function a))))))

;;; so 11*7*3 = (gprod 3 3 4) = 231

;;; N * (N-1m) * (N-2m) * ... * (0<x<=m)

(defun gfact (N m)
  (if (= (rem N m) 0)
      (gprod (ceiling n m)			; how many factors
	     m					; first number is m
	     m)
      (gprod (ceiling n m)			; how many factors
	     (rem N m)				; first number is rem(N,m)
	     m)))


;;;; n CHOOSE k

;;; The number of ways to place K distinquishable objects on
;;; N distinquishable plates.

;;; There are n! ways to arrange all the objects
;;; Break into two groups
;;;   there are k! ways to arrange one 
;;;         and (n-k)! ways to arrange the other

(defun choose-naive (n k)
  (/ (factorial n)
     (* (factorial (- n k)) (factorial k))))

;;; Minimize the number of multiplies
;;;
;;;                 10 9 8   7 6 5 4 3 2 1 
;;;   eg 10 C 3 =   ----------------------
;;;                 (3 2 1) (7 6 5 4 3 2 1)
;;;
(proclaim '(ftype (function (integer integer) integer) choose))
;;;
(defun choose (n k)
  (declare (integer n k))
  (let ((n-k (- n k)))
    (if (< n-k k)
	(choose n n-k)				; n C k = n C n-k
	(do ((i n (1- i))
	     (j k (1- j))
	     (c1 1)
	     (c2 1))
	    ((<= j 0) (the integer (/ c1 c2)))
	  (declare (integer i j c1 c2))
	  (setq c1 (* c1 i))
	  (setq c2 (* c2 j))))))

(defun log-choose (n k)
  (- (log-factorial n)
     (log-factorial k) (log-factorial (- n k))))

(defun dlog-choose (n k)
  (- (dlog-factorial n)
     (dlog-factorial k) (dlog-factorial (- n k))))

(defun fchoose (n k)
  (values (fround (exp  (log-choose n k)))))

(defun dchoose (n k)
  (values (fround (exp (dlog-choose n k)))))


;;;; Catalan Numbers

;;; CATALAN-NUMBERs are the number of ways to fully parenthesize a
;;; string of n symbols.
;;;	CN(1)=1
;;;	CN(n)=SUM((i=1,n-1) CN(i)*CN(n-i))
;;;		explicit: CN(n+1) = (choose 2n n)/(n + 1)

(defun catalan-number (n)
  (/ (choose (* 2 (- n 1)) (- n 1)) n))


;;;; Bell Numbers

;;; BELL-NUMBERs are the number of ways to place N distinquishable
;;; objects on N indistinquishable plates (SCI. AMERICAN)
;;;     there should be a much faster way to do this with mod arith

;;; B0 = 1 BY DEFINITION
;;; B(N) = (/ (sum (k = 1 inf) (expt k N)) e)
;;;  (EXP (EXP X)) = E (SUM (K = 0 INF) (BELL(K)) * X^K / K!)
;;;    SO B(N) = NTH DERIVATIVE OF (EXP (1- (EXP X))) EVALED AT 0
;;;  TOUCHARD'S CONGRUENCE:  B(P+K) = B(K) + B(K+1)  (MOD P, P PRIME)
;;;  B(N) < N^N
;;; (B0 B1 B2 ...) -> (1 1 2 5 15 52 203 877 ...)

;;; Bell triangle
;;;   1
;;;   1   2
;;;   2   3   5
;;;   5   7  10  15
;;;  15  20  27  37  52
;;;  52  67  87 114 151 203
;;; 203 255 322 409 523 674 877

(defun bellr (row col)
  (cond ((and (= col 1) (= row 1)) 1)
	((= col 1)                 (bellr (1- row) (1- row)))
	(t                         (+ (bellr (1- row) (1- col))
					 (bellr row      (1- col))))))

(defun bell-number (n)
  (cond ((minusp n) (error "BELL-NUMBER:  bad argument" n))
	((zerop  n) 1)
	(t          (bellr n n))))

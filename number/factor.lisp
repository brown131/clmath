;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Factoring and Totient

;;; References
;;;   Knuth

;;;  (c) Copyright Gerald Roylance 1983, 1984, 1985, 1986
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and fixes
;;;   Algorithm D doesn't work?

(in-package "CLMATH")

;;;; Primes

;;; return a list of primes up to but not including n
;;;
(defun primes-to (n)
  (declare (fixnum n))
  (do ((i 2 (1+ i))
       (l nil))
      ((>= i n) l)
    (declare (fixnum i))
    (do ((m l (cdr m)))
	((or (null m)
	     (> (expt (car m) 2) i))
	 (setq l (nconc l (list i))))
      (if (= 0 (mod i (car m)))
	  (return nil))
      )))


;;;; Integer Square Roots

(defun ceil-sqrt (n)
  (let ((f-sqrt (isqrt n)))
    (if (= n (expt f-sqrt 2))
	f-sqrt
	(1+ f-sqrt))))


;;;; Simple factoring

;;; the efficiency of this algorithm does not improve rapidly.
;;; (factor-build 4) gets about all that is convenient to get.
;;;   efficiency of factor-build 3 (  30) = 0.267
;;;   efficiency of factor-build 4 ( 210) = 0.205
;;;   efficiency of factor-build 5 (2310) = 0.147

(defparameter factor-k      30.)
(defparameter factor-trials '(     2.  3.  5.  7. 11. 13. 17. 19. 23. 29.))
(defparameter factor-sieve  '( 1.              7. 11. 13. 17. 19. 23. 29.))

(defun factor-build (p)
  (do ((lp (primes-to 50) (cdr lp))		; get a list of primes
       (j  1 (1+ j))
       (k  1))
      ((> j p)
       (setq factor-k      k)
       (setq factor-trials (primes-to k))
       (setq factor-sieve  (cons 1 (nthcdr p factor-trials)))
       (/ (float (length factor-sieve)) (float factor-k)))
    (declare (special factor-k factor-trials factor-sieve))
    (setq k (* k (car lp)))))

(defun factor (n)
  (declare (special factor-k factor-trials factor-sieve))
  (do ((limit   (isqrt n))
       (trials  factor-trials)
       (sieve   factor-sieve)
       (k       factor-k)
       (factors nil)
       (i0      0)
       (i       0))
      ((> i limit)
       (nreverse (cons n factors)))

    (cond ((null trials)
	   (setq i0 (+ i0 k))
	   (setq trials sieve)))

    (setq i (+ i0 (pop trials)))

    (do ()
	((or (not (zerop (mod n i)))
	     (> i limit)))
      (push i factors)
      (setq n (floor n i))
      (setq limit (isqrt n)))
    ))


;;;; Knuth

;;; Knuth, Volume II, pp 339-360

;;; Algorithm A -- just divide
;;; Algorithm B --
;;; Algorithm C -- x+y-sqrt(n) iterations
;;; Algorithm D -- x  -sqrt(n) iterations
;;; Algorithm E --

;;; Knuth Algorithm C, p 342
;;;  -- factoring by addition and subtraction
;;;
;;; Given an odd number n, find the largest factor of n <= sqrt(n)
;;;   (factors not necessarily prime)
;;;
(defun factor-knuth-c (n)
  (if (not (oddp n))
      (error "FACTOR-KNUTH-C: ~d is not odd" n))
  (let ((sqrt-n (isqrt n)))
    (do ((x-prime (+ (* 2 sqrt-n) 1))	;2x+1
	 (y-prime 1				;2y+1
		  (+ y-prime 2))
	 (r (- (expt sqrt-n 2) n)	;(x**2)-(y**2)-n
	    (- r y-prime)))
	((zerop r)
	 (list (floor (- x-prime y-prime) 2)	;largest factor
	       (floor (+ x-prime y-prime -2)    2)))
      (cond ((minusp r)
	     (setq r       (+ x-prime r))
	     (setq x-prime (+ x-prime 2))
	     ))
      )))


;;;; Knuth Algorithm D

;;; Knuth Algorithm D, p 345

;;; assume n = u * v, where u<=v
;;;   for practical purposes, n is odd, therefore u and v are odd
;;;   therefore let x = (u + v)/2, y = (v - u)/2
;;;        x**2 = (u**2 + 2*u*v + v**2)/4
;;;        y**2 = (u**2 - 2*u*v + v**2)/4
;;;        n = x**2 - y**2, 0 <= y < x <= n
;;;     so test is
;;;        x**2 - n = y**2 -- ie, (x**2 - n) is a perfect square

;;;  -- factoring with sieves
;;;
;;;   given an odd number n, this algorithm determines the largest
;;;   factor of n less than or equal to sqrt(n).  The procedure
;;;   uses moduli m1, m2, .. mr, which are relative prime to each
;;;   other in pairs and relatively prime to n.  We assume
;;;   r "sieve tables" S[i j] for 0<=j<mi, 1<=i<=r
;;;
;;;       where S[i j] = 1 if j**2-n = y**2 (mod mi) for some y
;;;
;;; This could be made faster if the loops were unwound.
;;; In addition, the moduli increment tests might be done
;;; only once in a while -- ie change Mmax to 2Mmax and perform
;;; the mod only every Mmax steps.  If the modulus is over Mmax,
;;; then knock it down by Mi * ceil(Mmax/Mi)
;;; instead of increment all Ki, just increment a Kincr
;;; and add it to each test.  There may be many moduli,
;;; but on average we only test 2 S[i j] per iteration.
;;; Flush the (add1 x) from the loop, too.  The bignum
;;; addition only needs to be done every 2**24 iterations
;;; to prevent overflow or when the perfect square test is done.


;;;; Thinking about similarities of Algorithm D

;;; given n and a modulus, return a list of mi and sieves

(defun quad-sieve (n mi)
  (let ((s  (list mi))
	(y2 (make-array mi			; quadratic residues
			:element-type    'fixnum
			:initial-element 0)))
    (declare (type (array fixnum (*)) y2))

    (do ((y 0 (1+ y)))				; clear residue table
	((>= y mi))
      (setf (aref y2 y) 0))

    (do ((y 0 (1+ y)))				; fill in table
	((>= y mi))
      (setf (aref y2 (mod (expt y 2) mi)) 1))

    (do ((x  0 (1+ x)))				; consider each x
	((>= x mi))
      (declare (fixnum x))
      (if (= 1 (aref y2 (mod (- (expt x 2) n) mi)))
	  (push x s)))

    (nreverse s)))

(defun sieve-unroll (n s)
  (do ((m (car s))
       (i 0 (1+ i))
       (x 0 (+ x m))
       (l nil))
      ((>= i n) (nreverse l))
    (do ((l0 (cdr s) (cdr l0)))
	((null l0))
      (push (+ x (car l0)) l))
    ))

(defun sieve-merge (s1 s2)
  ;; unroll s1 and s2
  ;; intersect them
  (do ((u1 (sieve-unroll (car s2) s1))
       (u2 (sieve-unroll (car s1) s2))
       (r  (list (* (car s1) (car s2)))))
      ((or (null u1) (null u2))
       (nreverse r))
    (cond ((= (car u1) (car u2))
	   (push (progn (pop u1) (pop u2)) r))
	  ((< (car u1) (car u2))
	   (pop u1))
	  (t (pop u2)))
    ))

(defun qtest (n)
  (print (quad-sieve n  3.))
  (print (quad-sieve n  5.))
  (print (quad-sieve n  8.))
  (print (quad-sieve n 17.))
  (print (quad-sieve n 64.))
  nil)


;;;; First initialize the tables

(defun factor-knuth-d-init (n r m s)
  (declare (type (array fixnum (*))   m)
	   (type (array fixnum (* *)) S))
  ;; initialize S
  (do ((i   0 (1+ i)))				; each modulus
      ((>= i r))
    (declare (fixnum i))
    (do ((x  0 (1+ x))				; each possible x
	 (mi (aref m i)))
	((>= x mi))
      (declare (Fixnum x mi))
      (cond ((do ((y    0 (1+ y))		; test if there is some y
		  (x2-n (mod (- (expt x 2) n) mi)))
		 ((>= y mi) NIL)
	       (declare (fixnum y x2-n))
	       (cond ((= (mod (expt y 2) mi) x2-n)
		      (return T))))

	     (do ((j x (+ j mi))
		  (l (array-dimension S 1)))
		 ((>= j l))
	       (declare (fixnum j l))
	       (setf (aref S i j) 1))
	     ))
      )))


;;;; Knuth Algorithm D (continued)

;;; this is only getting 3300 trials per second on a 3600
;;;                     41000 trials per second on a DEC20

(defun factor-knuth-d (n)
  (let* ((ceil-sqrt (ceil-sqrt n))
	 (moduli    '( 3.  5.  7.  8. 11. 13. 17. 19. 23. 29. 31. 37.
		      41. 43. 47. 53. 59. 61. 67.))
	 (R         (length moduli))
	 (kstep     (max 300. (apply #'max moduli)))
	 (Mmax      (* 3 kstep))
	 (M         (make-array r :element-type 'fixnum :initial-contents moduli))
	 (Ms        (make-array r :element-type 'fixnum :initial-element 0))
	 (K         (make-array r :element-type 'fixnum :initial-element 0))
	 (S         (make-array (list r mmax)
				:element-type 'fixnum :initial-element 0)))
    (declare (fixnum r kstep Mmax))
    (factor-knuth-d-init n r m s)		; init tables

    (do ((i    0 (1+ i))			; init moduli counters
	 (mcs (- ceil-sqrt))
	 (temp 0))
	((>=  i r))
      (declare (fixnum i temp))
      
      (setf (aref Ms i)				; increment to K counters
	    (* (aref m i)
	       (/ (* 2 kstep) (aref m i))))
      
      (setq temp (mod mcs (aref m i)))
      (setf (aref k i)
	    (if (< temp (- kstep 1))
		(+ temp (aref ms i))
		temp)))

    (catch 'prime-result
      (do ((x     ceil-sqrt)
	   (x0    ceil-sqrt)
	   (x00   0 (+ x00 kstep))
	   (init  (get-internal-run-time)))
	  (nil)
	(declare (fixnum x00))
	(cond ((> x00 100000)			; prevent overflow
	       (setq init (- (get-internal-run-time) init))	;  and say how fast
	       (format t "~% trials per second: ~f"
		       (/ (float x00)
			  (/ (float init)
			     (float internal-time-units-per-second))))
	       (setq init (get-internal-run-time))

	       (setq x0 (+ x0 x00))		; BIGNUM cons
	       (setq x00 0)))

	(do ((x1 0 (1+ x1)))			; fixnum step: x=x0+x00+x1
	    ((>= x1 kstep))
	  (declare (fixnum x1))
	  (cond ((do ((i 0 (1+ i)))		; Sieve Table Tests
		     ((>= i r) T)
		   (declare (fixnum i))
		   (cond ((= 0 (aref s i (- (aref k i) x1)))
			  (return nil))))
						; is x**2-n a perfect square?
		 (setq x (+ x0 x00 x1))		; calculate real X
		 (let* ((x2-n (- (expt x 2) n))
			(y    (isqrt x2-n)))
		   (cond ((equal (expt y 2) x2-n)
			  (throw 'prime-result (- x y)))))
		 )))
	(do ((i    0 (1+ i))			; step the ki
	     (temp 0))
	    ((>= i r))
	  (declare (fixnum i temp))
	  (setq temp (- (aref k i) kstep))
	  (setf (aref k i)
		(if (< temp (- kstep 1))
		    (+ temp (aref ms i))
		    temp)))
	)))
  )


;;;; Euler Totient Function

;;; number of integers not exceeding and relatively prime to n.

;;; Generating functions
;;;   (sigma (n=0 to inf) totient(n)*n**-s = zeta(s-1)/zeta(s))
;;;   (sigma (n=0 to inf) totient(n)*x**n/ (1-x**n) = x / (1-x)**2)

;;; Closed form
;;;   totient(n) = n (PI (all distinct primes dividing n) (1-(1/p)))

;;; Relations
;;;   recurrance
;;;     totient(m*n) = totient(m)*totient(n)    if (m,n) = 1

;;;   checks
;;;     (sigma (d divides n) totient(d)) = n
;;;     totient(n) = (sigma (d divides n) mobius(n/d) *d
;;;     a**totient(n) = 1 (mod n)     if (a,n)=1

(defun totient (n)
  (do ((factors (factor n) (cdr factors))
       (totient 1))
      ((null factors) totient)
    (cond ((or (null (cdr factors))
	       (not (equal (car  factors)
			   (cadr factors))))
	   (setq totient (* totient (- (car factors) 1)))))
    ))


;;;; Tests

(eval-when (eval)

  (IMPORT-FILE "OZ:OZ:<GLR.FUNCT.NUMBER>MOD.LISP")

  (defun totient-check (a n)
    (list (equal (gcd a n) 1)
	  (equal 1 (modpower a (totient n) n))))

  )

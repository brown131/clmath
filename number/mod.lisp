;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Modular Arithmetic

;;; References
;;;   Robert Solovay and Volker Strassen,
;;;     "A Fast Monte-Carlo Test for Primality,"
;;;     SIAM Journal on Computing, 1977, pp 84-85.
;;;   R. L. Rivest, A. Shamir, and L Adleman
;;;     "A Method for Obtaining Digital Signatures and Public-Key Cryptosystems"
;;;     Communications of the ACM,
;;;     Vol 21, Number 2, February 1978, pages 120-126

;;;  (c) Copyright Gerald Roylance 1983, 1984, 1985, 1986
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and fixes
;;;   split off the prime stuff
;;;   move random bignum elsewhere?

(in-package "CLMATH")

;;;; Greatest Common Divisor Algorithms

#|
(defun gcd (a b)
  (if (zerop b)
      a
      (gcd b (mod a b))))
|#

;;; find A and B such that
;;;    A * X + B * Y = Z

;;; initial equations:
;;;    1 * X + 0 * Y = X
;;;    0 * X + 1 * Y = Y

(defun gcd-extended (X Y)
  (do ((A1 1 A2)
       (B1 0 B2)
       (Z1 X Z2)
       (A2 0 (- A1 (* d A2)))
       (B2 1 (- B1 (* d B2)))
       (Z2 Y (- Z1 (* d Z2)))
       (d  0))
      ((zerop Z2) (list A1 X B1 Y Z1))
    (setq d (floor Z1 Z2))
    ))

;;; find b such that a * b is congruent to 1  (modulo modulus)

(defun modinv (a modulus)
  (let ((gcde (gcd-extended a modulus)))
    (if (not (equal (nth 4 gcde) 1))
	(error "MODINV:  no inverse because gcd not equal to 1" (nth 4 gcde))
	(car gcde))))


;;;; Chinese Remainder Theorem

;;; Find a number N such that
;;;   N mod p1 = u1
;;;   N mod p2 = u2
;;;   ... etc

;;; (defun chinese (p1 u1 p2 u2)				;should generalize
;;;   (+ (* p2 (modinv p2 p1) u1)
;;;         (* p1 (modinv p1 p2) u2)))

(defun chinese (num)					;(PRIME RES ...)
  (declare (fixnum num))
  (do ((P  (do ((i    1 (+ i 2))
		(prod 1))
	       ((> i num) prod)
	     (declare (fixnum i))
	     (setq prod (* prod (arg i)))))
       (i  1 (+ i 2))
       (ch 0))
      ((> i num) (mod ch p))
    (declare (fixnum i))
    (let* ((pr (arg i))
	   (ui (arg (1+ i)))
	   (ci (floor p pr))
	   (di (modinv ci pr)))
      (setq ch (+ ch (* ci di ui))))
    ))


;;;; (EXPT A N) MOD M

(defun modpower (number expon modulus)
  (do ((exp  expon  (floor exp 2))		; speedier to break into
						;  2**24 bit chunks?
       (sqr  number (mod (* sqr sqr) modulus))
       (ans  1))
      ((zerop exp) ans)
    (if (oddp exp)
	(setq ans (mod (* ans sqr) modulus)))))


;;;; Random Number Generator for BIGNUM Arguments

;;; Generate a random bignum that is L bits long

;;; generate random substrings of say 20 bits and
;;; paste them together to make a longer random string
;;; of course, this is a crock

(defun random-big-1 (length)
  (do ((l    length (- l k))			; number of bits to make
       (k    0)					; number of bits to make this pass
       (bits 0))				; rand bits so far
      ((<= l 0) bits)
    (declare (fixnum l k))
    (setq k (min l 20.))
    (setq bits (+ (* bits (lsh 1 k))		; shift left k bits
		  (random (lsh 1 k))))		; add in k bits
    ))

;;; Common Lisp has this function internally

(defun random-big (n)
  (if (< n 1)
      (error "bad argument to random-big")

      (do ((rn 0)
	   (l  (integer-length n)))
	  (NIL)
	(declare (integer rn)
		 (fixnum l))
	(setq rn (random-big-1 l))
	(if (< rn n)
	    (return rn))
	)))


;;;; Jacobi-Symbol

;;; The hairy exponent stuff here is just a hack to look
;;; at the lsbs of the bignums.  It has been hacked here
;;; to make it moderately fast without bumming it beyond
;;; recognition.

;;; the Jacobi-Symbol is always +1, -1, or 0

;;; (-1)**exp the easy way....
;;;
(defmacro jacobi-expt-1 (exp)
  `(if (oddp ,exp) -1 1))

;;; version from Sussman's notes
;;;
(defun jacobi-symbol (P Q)
  (let ((PP (mod P 16.))			; only need low order bits for
	(QQ (mod Q 16.)))			;  sometimes.  Used in place of
    (declare (fixnum PP QQ))			;  P or Q where it matters

    (cond ((equal P 0) 0)			; not in GJS notes
	  ((equal P 1) 1)
	  ((oddp PP)
	   (* (jacobi-symbol (mod q p) p)
	      (jacobi-expt-1 (/ (* (- PP 1) (- QQ 1)) 4))))
	  (t
	   (* (jacobi-symbol (floor p 2) Q)
	      (jacobi-expt-1 (/ (- (* QQ QQ) 1) 8)))))))


;;;; Prime Number Test:  Fermat (easily fooled)

;;; Fermat Test
;;;   choose a random number between 2 and N-1 inclusive
;;;   and check that a**n mod n = a

;;; the Carmichael Numbers fool this prime test
;;; there are 16 Carmichael Numbers below 100,000
;;;   the smallest are 561, 1105, 1729, 2465, 2821, 6601

(defun prime-test-fermat-1 (n)
  (let ((a (+ 2 (random-big (- n 2)))))
    (equal (modpower a n n) a)))

(defun prime-test-fermat (n &optional (trials 50.))
  (cond ((< n 2) NIL)
	((= n 2) 'PRIME)
	(t
	 (do ((i 0 (1+ i)))
	     ((> i trials) 'Probably-prime)
	   (declare (fixnum i))
	   (if (not (prime-test-fermat-1 n))
	       (return NIL))
	   ))))


;;;; Prime Number Test:  Solovay-Strassen

;;; Solovay-Strassen Prime Test
;;;   if n is prime, then J(a,n) is congruent mod n to a**((n-1)/2)
;;;
(defun prime-test-1 (a n)
  (and (= (gcd a n) 1)
       (= (mod (- (jacobi-symbol a n) (modpower a (floor (- n 1) 2) n)) n) 0)))

;;; checks if n is prime
;;;   probability of a mistake = (expt 2 (- trials))
;;;     choosing TRIALS=50 should be enough
;;;
(defun prime-test (n &optional (trials 50.))
  (setq n (abs n))
  (cond ((< n 2) nil)
	((= n 2) 'prime)
	((= n 3) 'prime)
	((and (> n 100)				; cheap test
	      (not (= 1 (gcd n
			     (* 02.  3.  5.  7.	; 28 bit number
				11. 13. 17. 19.
				23.)))))
	 nil)
	(t
	 (do ((i 0 (1+ i))
	      (a (random-big n) (random-big n)))
	     ((> i trials) 'probably-prime)
	   (declare (fixnum i))
	   (cond ((zerop a)			; this test is no good
		  (setq i (1- i)))
		 ((not (prime-test-1 a n))
		  (return nil)))))))



;;;; Choose a Random Prime of l bits

;;; number of primes is x / ln(x)
;;; probability of a hit (x / ln(x)) / x = 1 / ln(x)
;;; so about ln(x) / 2 trials have to be made
;;;   but I think this derivation is bad
;;; ln(2^x) = ln(2) * x
;;;    because 2^x = e^(x ln 2)

;;; *** This will not generate a 2!

(defun random-prime (l)

  (if (< l 2) (error "Bad Arg"))

  (do ((i 0 (1+ i)))
      (nil)
    (declare (fixnum i))

    (let ((p (+ 1 (* (random-big-1 (1- l)) 2))))
      (if (prime-test p)
	  (return p)))
    ))


;;;; Smaller Prime, List of Primes

;;; find a probable prime smaller than n
;;;
(defun smaller-prime (n)
  (cond ((< n 3)
	 (error "no smaller prime"))
	((= n 3)
	 2)
	(t
	 (let ((start (- n (if (oddp n) 2 1))))
	   (do ((i start (- i 2)))
	       (nil)
	     (if (prime-test i)
		 (return i)))))))


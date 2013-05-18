;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Probability and Statistics

;;; Probability Distributions
;;; Random Number Generators

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
;;;   The Generators should be tested for accuracy and speed
;;;   Underflow errors are trapping out

(in-package "CLMATH")

;;;; Random Number Generators

;;; from
;;;   A. J. Kinderman & J. F. Monahan
;;;     Computer Generation of Random Variables Using
;;;     the Ratio of Uniform Deviates
;;;     ACM Transactions on Mathematical Software
;;;     Vol 3 No 3 September 1977 pp 257-260

;;; Plan -- to generate f(x)
;;;   let Cf = {(u,v) st 0 <= u <= sqrt(f(v/u))}
;;;   (u,v) can be found by acceptance-rejection techniques
;;;         if Cf is bounded in u and v directions
;;;   Cf is bounded in u direction if     f(x) is bounded
;;;   Cf is founded in v direction if x*x*f(x) is bounded
;;;         u <=     sqrt(f(v/u))
;;;         v <= v/u sqrt(f(v/u))
;;;         v <=   y sqrt(f( y ))
;;;      v**2 <= y**2     f( y )
;;;   f(x) need not integrate to 1, say it integrates to K
;;;   area of Cf = K/2
;;;   P{acceptance} = (K/2) / area[u,v]
;;;   E[trials] = 1/P{acceptance}

;;; let h(x) be similar to f(x) (remove constants, etc)
;;; let u bound be max sqrt(h(x))
;;; let v bound be max x x h(x) -- I don't know if this is true!

;;; if x is finite, then just use accept/reject (eg, beta distribution)


;;;; Uniform Distribution

;;; Uniform Distribution
;;;   f(x) = 1, 0 <= x <= 1
;;;   mean = 0.5
;;;   vari = 0.0825

(defun uniform-density (x)
  (cond ((< x 0.0) 0.0)
	((> x 1.0) 0.0)
	(t         1.0)))

(defun uniform-cumulative (x)
  (cond ((< x 0.0) 0.0)
	((> x 1.0) 1.0)
	(t           x)))

;;; note -- x is never 0.0 or 1.0
;;;
(defun uniform-random-number ()
  ;; (/ (+ 0.5 (float (random 1000000.))) 1000000.0)
  (random 1.0))


;;;; Exponential Distribution

;;; RATE is actually LAMBDA*TIME

(defun exponential-density (x rate)
  (if (< x 0.0)
      0.0
      (* rate (exp (- (* rate x))))))

(defun exponential-cumulative (x rate)
  (if (< x 0.0)
      0.0
      (- 1.0 (exp (- (* rate x))))))


;;;; Exponential Random Number

;;; -- Ratio of Uniform Deviates Method

;;;   f(x) = (exp -x); x >= 0
;;;   mean = 1
;;;   vari = 1
;;; u-bound : 0 <= u <= 1
;;; v-bound :      u <= sqrt(exp(-(v/u)))
;;;                   = exp(-v / (2 u))
;;;         : log(u) <= -v / (2 u)
;;;         : 2 u log(u) <= -v
;;;         :      v <= -2 u log(u)
;;;    rhs is max when d(rhs)/du = 0
;;;        -2 (log(u) + 1) = 0
;;;          u = {exp(-1)}
;;;  therefore the v-bound is
;;;         :      v <= -2 u log(u) = -2 exp(-1) (-1)
;;;                                 =  2 exp(-1)
;;; P{acceptance} = 0.67957

(defun EXPONENTIAL-RANDOM-NUMBER (lmbd)
  (do ((u 0.0)
       (v 0.0))
      ((progn
	 (setq u (uniform-random-number)
	       v (* 2.0 (exp -1.0) (uniform-random-number)))
	 (<= v (* -2.0 u (log u))))
       (/ v u lmbd))
    (declare (float u v))))


;;;; Cauchy Distribution

;;; more general case in cauchy
;;;   (1/pi) (width/2) / ((x-mu)**2 +(width/2)**2)

(defun cauchy-density (x)
  (/ 1.0 (* (float pi 1.0) (1+ (expt x 2)))))

(defun cauchy-cumulative (x)
  (+ (/ (atan x 1.0) (float pi 1.0)) 0.5))


;;; -- Ratio of Uniform Deviates Method
;;;
;;;   f(x) = 1 / (pi * (1+x*x))
;;;
;;; P{ACCEPTANCE} = pi/4 = .78540

(defun cauchy-random-number ()
  (do ((u 0.0)
       (v 0.0))
      ((progn
	 (setq u (uniform-random-number)
	       v (* 2.0 (- (uniform-random-number) 0.5)))
	 (< (+ (* v v) (* u u)) 1.0))
       (/ v u))
    (declare (float u v))))


;;;; Normal Distribution

(defun NORMAL-DENSITY (x)
  (/ (exp (* -0.5 x x)) (sqrt (* 2.0 (float pi 1.0)))))

#|
;;; NBS 26.2.16
;;;
(defun NORMAL-CUMULATIVE (x)
  (if (< x 0.0)
      (- 1.0 (normal-cumulative (- x)))
      (let ((z (/ 1.0 (1+ (* 0.33267 x)))))
	(declare (float z))
	(-  1.0 (* (NORMAL-DENSITY x)		; eps < 1.0e-5
		   (horner-polynomial z
				      0.0
				      0.4361836
				      -.1201676
				      0.9372980))))))
|#

;;; NBS 26.2.17
;;;
(defun NORMAL-CUMULATIVE (x)
  (flet ((nctail (x)				; eps < 7.5e-8
	   (declare (ftype (function (float) float) nctail))
	   (let ((z (/ 1.0 (1+ (* 0.2316419 x)))))
	     (declare (float z))
	     (* (NORMAL-DENSITY x)
		(horner-polynomial z
				   0.0
				   0.319381530
				   -.356563782
				   1.781477937
				   -1.821255978
				   1.330274429)))))
    (if (< x 0.0)
	(nctail (- x))
	(- 1.0 (nctail x)))))

;;; NBS 26.2.22
;;; NORMAL-TAIL finds x such that
;;;     p = 1.0 - NORMAL-CUMULATIVE(x) = NORMAL-CUMULATIVE(-x)
;;;
(defun NORMAL-TAIL (p)				;0 < p <= .5
						;eps < 4.5e-4
  (cond ((or (< p 0.0) (> p 1.0))
	 (error "Bad Probability"))
	((> p 0.5) (- (normal-tail (- 1.0 p))))
	(t
	 (let ((s (sqrt (* -2.0 (log p)))))
	   (declare (float s))
	   (- s (/ (horner-polynomial s 2.515517 0.802853 0.010328)
		   (horner-polynomial s 1.000000 1.432788 0.189269 0.001308)))))))


;;;; Normal Random Number

;;; -- Ratio of Uniform Deviates Method

;;;   f(x) = (exp (-x*x/2))
;;;           the / (sqrt (* 2 pi))  doesn't matter
;;;   mean = 0
;;;   vari = 1
;;; u-bound : 0 <= u <= 1
;;; v-bound :      u <= sqrt(exp(-(v/u)**2/2)
;;;                   = exp(-v**2 / (4 u**2))
;;;         : log(u) <= -v**2 / (4 u**2)
;;;         : 4 u**2 log(u) <= -v**2
;;;         : v**2   <= -4 u**2 log(u)
;;;    rhs is max when d(rhs)/du = 0
;;;        -4 (2u log(u) + u) = 0
;;;        -4u (2 log(u) + 1) = 0
;;;          u = {0, exp(-0.5)}
;;;  therefore the v-bound is
;;;         : v**2   <= -4 u**2 log(u) = -4 exp(-1) (-0.5)
;;;                                    =  2 exp(-1)
;;;    -sqrt(2)exp(-0.5) <= v <= +sqrt(2)exp(-0.5)
;;;
;;; P{acceptance} = 0.73057

(defun NORMAL-RANDOM-NUMBER ()
  (do ((u 0.0)
       (v 0.0))
      ((progn
	 (setq u (uniform-random-number)	; U is bounded (0 1)
	       v (* 2.0 (sqrt 2.0) (exp -0.5)	; V is bounded (-MAX MAX)
		    (- (uniform-random-number) 0.5)))
	 (<= (* v v) (* -4.0 u u (log u))))	; < should be <=
       (/ v u))
    (declare (float u v))))


;;;; Gamma Distribution

(defun gamma-density (x a)
  (if (<= x 0.0)
      0.0
      (/ (* (expt x (1- a)) (exp (- x)))
	 (gamma-function a))))

(defun gamma-cumulative (x a)
  (if (<= x 0.0)
      0.0
      (/ (gamma-function-incomplete a x)
	 (gamma-function a))))


;;;; Gamma Random Number

;;; if X and Y are independent gamma of a and b,
;;;   then X+Y is gamma a+b

;;; gamma(0.5) is 0.5 Z^2 where Z is normal
;;; gamma(1) is the exponential, which can be computed using a log
;;; X = -log(U1 U2 U3 .. Un) is gamma(n)

;;; Knuth, Answer to 3.4.1 #16 (page 551)
;;; see J. H. Ahrens and U. Dieter, Computing, Vol 12 (1974), 233-246

;;; method for 0<a<1 is exercise 16
;;;   the gamma density function f(t) is bounded by c g(t)
;;;      c g(t) = t^{a-1}/\Gamma(a) for 0<t<1
;;;             = e^-t   /\Gamma(a)     1<=t
;;;      since integral of g(t) = 1
;;;         C = a e / (a + e)

;;; there are hacks to make this go faster

(defun gamma-random-number-sub01 (a)
  (if (or (< a 0.0) (>= a 1.0))
      (error "GAMMA-RANDOM-NUMBER-SUB01"))
  (let* ((e (exp 1.0))
	 (p (/ e (+ a e))))			; probability G1 should be used
    (declare (float e p))
    (do ()					; < 1.4 iterations
	(NIL)
      (let ((U (uniform-random-number))
	    (V (uniform-random-number))
	    (X 0.0)
	    (q 0.0))
	(declare (float U V X q))
	(if (< U p)
	    (setq X (expt V (/ a))
		  q (exp (- X)))
	    (setq X (- 1.0 (log V))
		  q (expt X (- a 1.0))))
	;; now X has density g and q = f(X) / c g(X)
	(if (< (uniform-random-number) q)
	    (return X))))))

;;; Knuth Algorithm A, page 129

;;; see J. H. Ahrens and U. Dieter, Computing, Vol 12 (1974), 233-246
;;;  there is a version that is 2 to 3 times faster

(defun gamma-random-number-sub (a)
  (if (<= a 1.0)
      (error "GAMMA-RANDOM-NUMBER-SUB:  bad argument")
      (do ((a-1      (- a 1.0))
	   (sqrt2a-1 (sqrt (- (* 2.0 a) 1.0))))
	  (NIL)
	(declare (float a-1 sqrt2a-1))
	(let* ((U (uniform-random-number))	; executed < 1.902 if a >= 3
	       (Y (tan (* (float pi 1.0) U)))
	       (X (+ (* sqrt2a-1 Y) a-1)))
	  (declare (float U Y X))
	  (if (> X 0.0)
	      (if (<= (log (/ (uniform-random-number) (+ 1.0 (expt Y 2))))
		      (- (* a-1 (log (/ X a-1))) (* sqrt2a-1 Y)))
		  (return X)))))))

(defun gamma-random-number-subi (a)
  (do ((x a   (- x 1.0))
       (u 1.0 (* u (uniform-random-number))))
      ((< x 1.0)
       (- (if (<= x 0.0) 0.0 (gamma-random-number-sub01 x))
	  (log u)))
    (declare (float x u))
    ))

(defun gamma-random-number (a)
  (cond ((<= a 0.0) (error "GRX:  bad arg"))
	((>= a 3.0) (gamma-random-number-sub   a))
	((<  a 1.0) (gamma-random-number-sub01 a))
	(t          (gamma-random-number-subi  a))))


;;;; Gamma Random Number

;;; -- Ratio of Uniform Deviates Method

;;; f(x) = (1/gamma(a)) (x**(a-1))exp(-x)
;;; h(x) = (x**(a-1))exp(-x)
;;; max of     h(x) is at x=a-1 for a>1
;;; max of x x h(x) is at x=a+1

#|

(defun gamma-0-random-number (a x)
  (* (expt x (- a 1.0)) (exp (- x))))

(defun gamma-random-number (fa)
  (let ((a (float fa)))
    (declare (float a))
    (if (<= a 1.0)
	(error "gamma-random-number -- arg < 1.0")
	(do ((ul (sqrt (gamma-0-random-number a (- a 1.0))))
	     (vl (gamma-0-random-number (+ a 2.0) (+ a 1.0)))
	     (u  0.0)
	     (v  0.0)
	     (x  0.0))
	    ((progn (setq u (* ul (uniform-random-number)))
		    (setq v (* vl (uniform-random-number)))
		    (setq x (/ v u))
		    (<= (* 2.0 (log u))
			(+ (* (- a 1.0) (log x)) (- x))))
	     x)
	  (declare (float ul vl u v x))))))

|#


;;;; Beta Distribution

(defun beta-density (x m n)
  (cond ((< x 0.0) 0.0)
	((> x 1.0) 0.0)
	(t
	 (/ (* (expt x          (1- m))
	       (expt (- 1.0 x) (1- n)))
	    (beta-function m n)))))

(defun beta-cumulative (x m n)
  (cond ((<= x 0.0) 0.0)
	((>= x 1.0) 1.0)
	(t
	 (/ (beta-function-incomplete m n x)
	    (beta-function m n)))))


;;;; Beta Random Number

;;; Knuth pp 129-130

;;; Let X1 and X2 be gamma of a and b
;;;  X = X1 / (X1 + X2)
;;;
(defun beta-random-number (m n)
  (let ((x1 (gamma-random-number (float m)))
	(x2 (gamma-random-number (float n))))
    (declare (float x1 x2))
    (/ x1 (+ x1 x2))))

;;; method for small a and b
;;;   ...

;;; method for integers
;;;   ...

;;; direct method
;;;   R. C. H. Cheng, CACM Vol 21 (1978), 317-322


;;;; Beta Random Number

;;; -- Ratio of Uniform Deviates Method

;;; m and n must be > 1
;;;   if they are less than 1, then we would evaluate 0^{-eps}
;;;     which is unbounded
;;;   if either is 1, then we have an 0^0 that must be resolved
;;;     also, the bound checks would then be wrong

;;; f(x) = (expt x (1- m)) (expt (1- x) (1- n)) / BETA(m,n)
;;;
;;; max of   f(x) is at (  m  -1) / (  m  +n-2)
;;; max of xxf(x) is at ((m+2)-1) / ((m+2)+n-2) = (m+1) / (m+n)

#|
(defun beta-0-random-number (m n x)
  (* (expt        x  (1- m))
     (expt (- 1.0 x) (1- n))))

(defun beta-random-number (m n)
  (let* ((fm (float m))
	 (fn (float n))
	 (ux (/ (1- fm) (+ fm fn -2)))
	 (u1 (sqrt (beta-0-random-number fm fn ux)))
	 (vx (/ (1+ fm) (+ fm fn)))
	 (v1 (* vx (sqrt (beta-0-random-number fm fn vx)))))
    (declare (float fm fn ux u1 vx v1))
    ;; (format t "~&Area of (u,v) = ~f" (* u1 v1))
    ;; (format t "~&   efficiency = ~f" (/ (* u1 v1) (/ (beta-function fm fn) 2.0)))
    (do ((u  0.0)
	 (v  0.0))
	((progn
	   (setq u (* u1 (uniform-random-number))
		 v (* v1 (uniform-random-number)))
	   (and (<= v u)			; v/u is bounded!
		;; (< (expt u 2) (beta-0-random-number m n (/ v u)))
		(< (* 2.0 (log u))
		   (+ (* (1- fm) (log        (/ v u)))
		      (* (1- fn) (log (- 1.0 (/ v u))))))))
	 (/ v u))
      (declare (float u v)))))

(time (beta-random-number 50 20))
(time (beta-random-number 10 20))
(time (beta-random-number  3 20))
(time (beta-random-number  3  5.3))
(time (beta-random-number 50 200))		; barfs

|#


;;;; Chi-Square Distribution

;;; CHI-SQUARE = (SUM (i 1 n) Yi**2)
;;; where the Yi are Normal(0, 1) independent RV's
;;; CHI-SQUARE >= 0

#|

(defun CHI-SQUARE-DENSITY (chi-sq n)
  (if (< chi-sq 0.0)
      0.0
      (let ((a (* 0.5 (float n))))
	(declare (float a))
	(/ (* (expt chi-sq (1- a)) (exp (* -0.5 chi-sq)))
	   (* (expt 2.0 a)         (gamma-function a))))))

|#

;;; CHI-SQUARE-DENSITY    is the same as 0.5 GAMMA-DENSITY(chi-sq/2 n/2)
;;;
(defun chi-square-density (chi-sq n)
  (* 0.5 (gamma-density (/ chi-sq 2.0) (/ (float n) 2.0))))

;;; CHI-SQUARE-CUMULATIVE is the same as GAMMA-CUMULATIVE(chi-sq/2 n/2)
;;;
(defun chi-square-cumulative (chi-sq n)
  (gamma-cumulative (/ chi-sq 2.0) (/ (float n) 2.0)))


;;;; Chi-Square Random Number

;;; -- Direct Method

#|

(defun chi-square-random-number (M)
  (do ((sum 0.0)
       (cnt 1 (1+ cnt)))
      ((> cnt m) sum)
    (declare (float sum) (fixnum cnt))
    (setq sum (+ sum (expt (normal-random-number) 2)))))

|#

;;; Knuth page 130

(defun chi-square-random-number (n)
  (* 2.0 (gamma-random-number (/ (float n) 2.0))))


;;;; Chi-Square Random Number

;;; -- Ratio of Uniform Deviates Method

;;; LET g(x) = (* (expt x (- n/2 1)) (exp -x/2))
;;; max x x g(x)
;;;          = (* (expt x (+ n/2 1)) (exp -x/2))
;;; or max log()
;;;          = (+ (* (+ n/2 1) (log x)) -x/2)
;;; derivative
;;;        0 = (+ (* (+ n/2 1) (/ x))  -1/2)
;;;        x = (+ n 2)

#| diked out --- 

(defun chi-square-0-random-number (n x)
  (* (expt x (1- (/ (float n) 2.0)))
     (exp  (* -0.5 x))))

(defun chi-square-random-number (n)
  (declare (float u v vl xmax))
  (do ((vl (let ((xmax (+ 2.0 (float n))))
	     (* xmax (sqrt (chi-square-0-random-number n xmax)))))
       (U 0.0) (v 0.0))
      ((progn
	 (setq u (uniform-random-number)
	       V (* vl (uniform-random-number)))
	 (< (expt u 2) (chi-square-0-random-number n (/ v u))))
       (/ v u))))

|#


;;;; Student's t Distribution

;;; if X is normally distributed N(0,sigma^2)
;;; and Y^2 / sigma^2 is chi-square with n dof
;;; and X and Y are independent then
;;; t = X (sqrt n) / Y is student's t with n dof

;;; mean = 0, variance = n / (n-2) for n>2
;;; -inf < tt < +inf

(defun t-density (tt n)
  (let ((nf (float n)))
    (declare (float nf))
    (/ (gamma-function (/ (1+ nf) 2.0))
       (* (float (sqrt (* nf pi)) 1.0)
	  (gamma-function (/ nf 2.0))
	  (expt (1+ (/ (expt tt 2) nf))
		(* 0.5 (1+ nf)))))))

;;; TWO-SIDED (integral from -t to t)
;;; NBS 27.6.1
;;;   A(t, v) = 1-I[x](v/2, 1/2)
;;;   x = v/(v+t**2)
;;;
;;; Symmetry (NBS 6.6.3)
;;;    I[x](a,b) = 1 - I[1-x](b,a)
;;;   so
;;;    I[1-x](b,a) = 1 - I[x](a,b)
;;;
;;; therefore 
;;;   A(t, v) = I[1-x](1/2, v/2)
;;;
(defun t-two-sided (tt n)
  (let* ((v (float n))
	 (x (/ v (+ v (expt tt 2)))))
    (declare (float v x))
    (/ (beta-function-incomplete 0.5 (/ v 2.0) (- 1.0 x))
       (beta-function 0.5 (/ v 2.0)))))

(defun t-cumulative (tt n)
  (let ((p (t-two-sided tt n)))
    (if (< tt 0.0)
	(- 0.5 (* 0.5 p))
	(+ 0.5 (* 0.5 p)))))


;;;; Student's t Random Number

;;; -- Direct Method

(defun T-RANDOM-NUMBER (n)
  (/ (normal-random-number)
     (sqrt (/ (chi-square-random-number n) (float n)))))


;;;; Student's t Random Number

;;; *** This gets underflow errors -- ie (exp -122.0)

;;; -- Ratio of Uniform Deviates Method

;;; LET g(x) = (expt (+ (/ (^ X 2) n) 1) (* -0.5 (+ n 1)))
;;; max x x g(x) = (* x x (expt (+ (/ (^ X 2) n) 1) (* -0.5 (+ n 1))))
;;; or max log = 
;;;   (+ (* 2 (log x)) (* -0.5 (+ n 1) (log (+ (/ (^ X 2) n) 1))))
;;; derivative = 0
;;; (+ (* 2 (/ x)) (* -0.5 (+ n 1) (/ (+ (/ (^ X 2) n) 1)) (/ (* 2 X) n)))
;;; (+ (* 2) (* x -0.5 (+ n 1) (/ (+ (/ (^ X 2) n) 1)) (/ (* 2 X) n)))
;;; (+ (* 2 (+ (/ (^ X 2) n) 1)) (* x -0.5 (+ n 1) (/ (* 2 X) n)))
;;; (+ (* 2 (+ (/ (^ X 2) n) 1)) (* (/ (^ x 2) n) -0.5 (+ n 1) 2))
;;; (+ (* 2 (+ (/ (^ X 2) n) 1)) (* -1 (/ (^ x 2) n) (+ n 1)))
;;; (+ (* 2 (+ (^ X 2) n)) (* -1 (^ x 2) (+ n 1)))
;;; (+ (* 2 (^ X 2)) (* 2 n) (* -1 (^ x 2) (+ n 1)))
;;; (+ (* 2 n) (* 2 (^ X 2)) (* -1 (^ x 2) (+ n 1)))
;;; (+ (* 2 n) (* (^ x 2) (+ 2 (* -1 (+ n 1)))))
;;; (+ (* 2 n) (* (^ x 2) (+ 2 -n -1)))
;;; (+ (* 2 n) (* (^ x 2) (+ -n 1)))
;;; (^ x 2) = (/ (* 2 n) (- n 1))
;;; x = (sqrt (/ (* 2 n) (- n 1)))

#|
(defun T-0-RANDOM-NUMBER (n x)
  (expt (1+ (/ (expt x 2) (float n)))
	(* -0.5 (1+ (float n)))))

(defun T-RANDOM-NUMBER (n)
  (let ((xmax (sqrt (/ (* 2.0 (float n)) (1- (float n))))))
    (declare (float xmax))
    (do ((vl (* xmax (sqrt (t-0-random-number n xmax))))
	 (u 0.0)
	 (v 0.0))
	((progn
	   (setq u (uniform-random-number)
		 v (* vl 2.0 (- (uniform-random-number) 0.5)))	; two-sided
	   (< (expt u 2) (t-0-random-number n (/ v u))))
	 (/ v u))
      (declare (float u v vl)))))
|#


;;;; Snedecor's F Distribution

;;; if X is distributed as chi-square with m dof
;;;    Y is distributed as chi-square with n dof
;;; and X, Y are independent
;;; then F = (X / m) / (Y / n) is F(m n)

;;; mean = n / (n - 2)       n>2
;;; variance = 2 n^2 (m + n -2) / (m (n - 2)^2 (n - 4)))
;;;         n>4

;;; the transform w = (m F / n) / (1 + (m F / n)
;;; transforms F-density to beta density

(defun F-density (F m n)				;0<F<+inf
  (let ((m2 (/ (float m) 2.0))
	(n2 (/ (float n) 2.0)))
    (declare (float m2 n2))
    (/ (* (gamma-function (+ m2 n2))
	  (expt (/ m2 n2) m2)
	  (expt F (1- m2)))
       (* (gamma-function m2)
	  (gamma-function n2)
	  (expt (1+ (/ (* m2 F) n2))
		(+ m2 n2))))))

(defun F-cumulative (F m n)
  (let ((m2 (/ (float m) 2.0))
	(n2 (/ (float n) 2.0)))
    (declare (float m2 n2))
    (/ (beta-function-incomplete m2 n2 (/ (* m2 F)
					  (+ n2 (* m2 F))))
       (beta-function m2 n2))))

;;;  F(F 2 4) =  y
;;; 10.65       .975
;;; 18.00       .990
;;; 26.28       .995
;;; 61.25       .999


;;;; Snedecor's F Random Number

;;; -- Direct Method

(defun F-RANDOM-NUMBER (m n)
  (/ (/ (chi-square-random-number m) (float m))
     (/ (chi-square-random-number n) (float n))))


;;;; Snedecor's F Random Number

;;; -- Ratio of Uniform Deviates Method

;;; LET g(x) = 
;;;   (* (expt x (- v1/2 1)) (expt (+ v2 (* v1 x)) (* -0.5 (+ v1 v2))))
;;; max x x g(x) =
;;;   (* (expt x (+ v1/2 1)) (expt (+ v2 (* v1 x)) (* -0.5 (+ v1 v2))))
;;; or max log() =
;;;   (+ (* (+ v1/2 1) (log x)) (* -0.5 (+ v1 v2) (log (+ v2 (* v1 x)))))
;;; or max
;;;   (- (* (+ v1 2) (log x)) (* (+ v1 v2) (log (+ v2 (* v1 x)))))
;;; derivative = 0
;;;   0 = (- (* (+ v1 2) (/ x)) (* (+ v1 v2) (/ (+ v2 (* v1 x))) v1))
;;;   0 = (- (* (+ v1 2)      1) (* (+ v1 v2) (/ (+ v2 (* v1 x))) v1 x))
;;;   0 = (- (* (+ v1 2) (+ v2 (* v1 x))) (* (+ v1 v2) 1 v1 x))
;;;   0 = (+ (* (+ v1 2) (+ v2 (* v1 x))) (* -1 (+ v1 v2) v1 x))
;;;   0 = (+ (* (+ v1 2) v2) (* (+ v1 2) (* v1 x)) (* -1 (+ v1 v2) v1 x))
;;;   0 = (+ (* (+ v1 2) v2) (* (+ v1 2) v1 x) (* -1 (+ v1 v2) v1 x))
;;;   0 = (+ (* (+ v1 2) v2) (* x v1 (+ (+ v1 2) (* -1 (+ v1 v2)))))
;;;   0 = (+ (* (+ v1 2) v2) (* x v1 (+ v1 2 (* -1 (+ v1 v2)))))
;;;   0 = (+ (* (+ v1 2) v2) (* x v1 (+ v1 2 -v1 -v2)))
;;;   0 = (+ (* (+ v1 2) v2) (* x v1 (+ 2 -v2)))
;;;   0 = (+ (* (+ v1 2) v2) (* x v1 (- 2 v2)))
;;;   x = (/ (* (+ v1 2) v2) (* v1 (- v2 2)))

;;; ugly !! v2 > 2 -- why ???

#| diked out --- 

(defun f-0-random-number (v1 v2 x)
  (* (expt x (- (/ (float v1) 2.0) 1.0))
     (expt (+ (float v2) (* (float v1) x))
	   (* -0.5 (+ (float v1) (float v2))))))

(defun f-random-number (v1 v2)
  (declare (float u v vl xmax))
  (do ((vl (let ((xmax (/ (* (+ (float v1) 2.0) (float v2))
			  (* (float v1) (- (float v2) 2.0)))))
	     (* xmax (sqrt (f-0-random-number v1 v2 xmax)))))
       (u 0.0) (v 0.0))
      ((progn
	 (setq u (uniform-random-number)
	       v (* vl (uniform-random-number)))
	 (< (expt u 2) (f-0-random-number v1 v2 (/ v u))))
       (/ v u))))
|#

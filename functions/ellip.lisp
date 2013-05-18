;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Elliptic Functions

;;; Reference
;;;   M. Abramowitz and I. Stegun, eds,
;;;     Handbook of Mathematical Functions,
;;;     National Bureau of Standards, 1964

;;;  (c) Copyright Gerald Roylance 1982, 1984
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   BEWARE!  This code uses the modulus m rather
;;;     than the parameter k.  (m = k**2)

(in-package "CLMATH")

;;;; ELLIPTIC INTEGRAL K(M)

;;; (integral 0 pi/2 (1 - m*sin^2(theta))^-0.5 theta)

;;; k'(m) = k(1-m)

;;; 0 <= M < 1 

(defun elliptic-integral-k (m)
  (elliptic-integral-kC (- 1.0 (float m))))

(defun elliptic-integral-kC (m) 			;eps < 2E-8
  (+ (poly m
	   1.38629436112
	   0.09666344259
	   0.03590092383
	   0.03742563713
	   0.01451196212)
     (* (- (log m))
	(poly m
	      0.50000000000
	      0.12498593597
	      0.06880248576
	      0.03328355346
	      0.00441787012))))

#|

(defun elliptic-integral-kC (m)				;eps < 3E-5
  (+$ (poly m
	    1.3862944
	    0.1119723
	    0.0725296)
      (*$ (-$ (log m))
	  (poly m
		0.5000000
		0.1213478
		0.0288729))))

K(0.1) = 1.61244 13487
K(0.9) = 2.57809 21133

|#


;;;; ELLIPTIC INTEGRAL E(M)

;;; integral from 0 to pi/2 d[theta] (sqrt 1 - m*sin^2(theta))

(defun elliptic-integral-e (m)
  (elliptic-integral-ec (- 1.0 m)))

(defun elliptic-integral-eC (m) 			;eps < 2E-8
  (+ (poly m
	   1.00000000000
	   0.44325141463
	   0.06260601220
	   0.04757383546
	   0.01736506451)
     (* (- (log m))
	(poly m
	      0.00000000000
	      0.24998368310
	      0.09200180037
	      0.04069697526
	      0.00526449639))))

#|
(defun elliptic-integral-eC (m)				;eps < 3E-5?
  (+ (poly m
	   1.0000000
	   0.4630151
	   0.1077812)
     (* (-$ (log m))
	(poly m
	      0.0000000
	      0.2452727
	      0.0412496))))

E(0.1) = 1.53075 7637
E(0.9) = 1.10477 4733

|#


;;;; ELLIPTIC-SINE (u M)

;;; u(phi, M) = (Integral (0 to phi) (1 - M sin^2 x)^-0.5  (d x))

;;; elliptic-sine(u M) = sin(phi)

(defun elliptic-sine (u M)
  (declare (single-float K K1 v q v2N1 qN05 sn))
  (let* ((K   (elliptic-integral-k M))
	 (K1  (elliptic-integral-kc M))
	 (v   (/ (* (float pi 1.0) u) (* 2.0 K)))
	 (q   (exp (- (/ (* (float pi 1.0) K1) K)))))
    (do ((v2N1 v (+ v2N1 v v))			;(2N+1)v
	 (qN05 (sqrt q) (* qN05 q))		;q^(N + 0.5)
	 (sn 0.0) )
	((progn (setq sn (+ sn (/ (* qN05 (sin v2N1))
				  (- 1.0 (* qN05 qN05)))))
		(< qN05 1.0E-8))
	 (/ (* (* 2.0 (float pi 1.0)) sn)
	    (* (sqrt M) K)))
      )))

(defun elliptic-cosine (u M)
  (declare (single-float K K1 v q v2N1 qN05 sn))
  (let* ((K   (elliptic-integral-k M))
	 (K1  (elliptic-integral-kc M))
	 (v   (/ (* (float pi 1.0) u) (* 2.0 K)))
	 (q   (exp (- (/ (* (float pi 1.0) K1) K)))))
    (do ((v2N1 v (+ v2N1 v v))			;(2N+1)v
	 (qN05 (sqrt q) (* qN05 q))		;q^(N + 0.5)
	 (sn 0.0) )
	((progn (setq sn (+ sn (/ (* qN05 (cos v2N1))
				  (+ 1.0 (* qN05 qN05)))))
		(< qN05 1.0E-6))
	 (/ (* (* 2.0 (float pi 1.0)) sn)
	    (* (sqrt M) K)))
      )))


;;;; Some tests

;NBS 16.23
;( 16.23 (Q=(EXP -PI K'/K))
;	 (V=PI U / 2 K)
;	 (SN(U,M) = (2 PI / SQRT(M) K) *
;	            (SIGMA (0 T0 INF)
;			   (* ((^ Q (N+0.5)) / (- 1 Q^(2N+1)))
;			      (SIN (2N+1)V))
;			   ))
;	 )

;;; FOR CN(U,M) CHANGE SIN TO COS AND (- 1 ..) TO PLUS
;;; SN(.5K,M) = 1 / (SQRT (1+ SQRT(M1)))
;;; CN(.5K,M) = (FOURTHROOT M1) / (SQRT (1+ SQRT(M1)))

#+ignore
(defun sntest (m)
  (setq m (float m))
  (let* ((K05  (/ (elliptic-integral-k m) 2.0))
	 (SQM1 (sqrt (- 1.0  M)))
	 (SN   (elliptic-sine K05 M))
	 (CN   (elliptic-cosine K05 M)))
    (print (list 'SN SN
		 (- SN (/ 1.0
			  (sqrt (+ 1.0 SQM1))))))
    (print (list 'CN CN
		 (- CN (/ (sqrt SQM1)
			  (sqrt (+ 1.0 SQM1)))))))
  nil)

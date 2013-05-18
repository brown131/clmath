;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Error Function

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

;;;  (c) Copyright Gerald Roylance 1982, 1984, 1985, 1988
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and fixes
;;;   

(in-package "CLMATH")

;;;; Error Function

;;; NBS 7.1.1
;;;   erf(z) = (2 / (sqrt pi)) (integral (t= 0 z) (exp -t**2))

;;; NBS 7.1.9
;;;   erf(-x) = - erf(x)

;;; NBS 7.1.5
;;;   erf z = (2 / (sqrt pi)) (sigma (n = 0 inf) (-1)**n z**(2n+1) / (n! (2n+1)))

;;; NBS 7.1.6
;;;   erf z = (2 / (sqrt pi)) (exp -z**2)
;;;           * (sigma (n = 0 inf) 2**n z**(2n+1) / (1 * 3 * ... (2n+1)))


;;;; Error Function

;;; derived from NBS 7.1.25

#|
(defun erfc (X)					; 0 <= X
  (if (< x 0.0)
      (error "erfc:  Negative Argument ~f" x)
      (let ((z (/ 1.0 (1+ (* 0.47047 X)))))
	(declare (float z))
	(- 1.0 (* (EXP (- (* X X)))
		  (POLY Z 
			0.0			; eps < 2.5e-5
			0.3480242
			-.0958798
			0.7478556))))))
|#

;;; derived from NBS 7.1.26
;;;   eps < 1.5e-7
;;;
(defun erfc (x)
  (declare (float x))
  (if (< x 0.0)
      (error "erfc:  Negative Argument ~f" x)
      (let ((z (/ 1.0 (1+ (* 0.3275911e0 x)))))
	(declare (float z))
	(* (exp (- (* x x)))
	   (poly z 
		 +0.000000000e0
		 +0.254829592e0
		 -0.284496736e0
		 +1.421413741e0
		 -1.453152027e0
		 +1.061405429e0)))))

(defun erf (x)
  (declare (float x))
  (cond ((<  x 0.0) (- (erf (- x))))
	((>= x 5.0) 1.0)
	(t
	 (- 1.0 (erfc x)))))

;;; NBS 7.1.26
;;;
(defun error-function (x)
  (declare (float x))
  (erf x))


;;;; Tests

#|

(- (erf 0.0  ) 0.0000000000d0)
(- (erf 0.5d0) 0.5204998778d0)
(- (erf 1.0  ) 0.8427007929d0)

|#

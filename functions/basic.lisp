;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Basic Mathematical Functions

;;;  (c) Copyright Gerald Roylance 1982, 1986
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   use common lisp names
;;;   number stuff should be macros
;;;   floor, ceiling, amod, ...

(in-package "CLMATH")

;;;; Trivial Stuff

(defun fraction (x) (- x (float (ifix x))))		; x < some value

(defun entier (x) (float (ifix x)))			; x < some value

;;; Functions to Hack Number Representations
;;;

(DEFUN FLOATING-EXPONENT (X)
  #+MACLISP (- (*LDB   #o3310_24. (ABS X)) #o0200)	; extract exp (PDP-10)
  #+lispm (multiple-value-bind (signif expon sign)	; god damn 3600
	      (si:decode-float x)
	    (declare (ignore signif sign))
	    expon)
  )

#+MACLISP
(DECLARE (FIXNUM (FLOATING-EXPONENT FLOAT)))	; ***GROSS FAKE OUT***

(DEFUN FLOATING-SCALE    (X E)
  #+MACLISP (FSC X E)
  #+LISPM (ASH X E))
#|
(DEFUN FLOATING-MANTISSA (X)
  (FLOATING-SCALE X (- (FLOATING-EXPONENT X))))
|#

;;;; Square Root

;;; after 2 iterations eps < 2^-31

;;; algorithm is from
;;;   Digital Equipment Corporation,
;;;   PDP-11 Paper Tape Software Programming Handbook,
;;;   DEC-11-GGPA-D, Section 7.7.

;;; compute square root of the exponent
;;; Chebyshev 8 bit approximation to square root of the mantissa
;;; two Heron iterations --> 32 bit accuracy
#|
(DEFUN SQUARE-ROOT (X)
  (LET* ((EXP (FLOATING-EXPONENT X))		; extract exponent
	 (MAN (FLOATING-SCALE    X (- EXP)))	; extract mantissa
	 (YI  (+ 0.41730760			; initial mantissa
		  (* 0.59016207 MAN))))
    (DECLARE (FIXNUM EXP)
	     (FLOAT MAN YI))
    
    (SETQ YI (IF (ODDP EXP)
		 (* 1.4142135623730950488	; sqrt(2)
		     (FLOATING-SCALE YI (ASH (- EXP 1) -1)))
		 (FLOATING-SCALE YI (ASH EXP -1))))

    (SETQ YI (* 0.5 (+ YI (/ X YI))))	; iteration 1
    (SETQ YI (* 0.5 (+ YI (/ X YI))))	; iteration 2
    ))
|#
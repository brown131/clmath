;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Root Finding Subroutines

;;;  (c) Copyright Gerald Roylance 1982
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   Fix false position searches
;;;     false position searches die if y1 = y2
;;;     should start searching for other x's to try

(in-package "CLMATH")

;;;; False Position Search

;;; find x such that 0 = f(x)
;;;
(defun false-position-search (fcn s1 s2 eps)
  (do ((x1 s1)  (y1 (funcall fcn s1))
       (x2 s2)  (y2 (funcall fcn s2))
       (xn 0.0) (yn 0.0))
      (NIL)
    (declare (float x1 x2 y1 y2 xn yn))
    (if (= y1 y2) (error "FALSE-POSITION-SEARCH Lost"))
    (setq xn (- x1 (* y1 (/ (- x2 x1) (- y2 y1)))))
    (setq yn (funcall fcn xn))
    (cond ((< (abs yn) eps) (return xn)))
    (cond ((> (abs y1) (abs y2))
	   (setq x1 xn) (setq y1 yn))
	  (t
	   (setq x2 xn) (setq y2 yn)))))


;;;; False-position-converge

;;; find an x such that x = f(x)
;;;
(defun false-position-converge (fcn s1 s2 eps)
  (do ((x1 s1)  (y1 (- (funcall fcn s1) s1))
       (x2 s2)  (y2 (- (funcall fcn s2) s2))
       (xn 0.0) (yn 0.0))
      (NIL)
    (declare (float x1 x2 y1 y2 xn yn))
    (if (= y1 y2) (error "FALSE-POSITION-CONVERGE lost"))
    (setq xn (- x1 (* y1 (/ (- x2 x1) (- y2 y1)))))
    (setq yn (- (funcall fcn xn) xn))
    (cond ((< (abs yn) eps) (return xn)))
    (cond ((> (abs y1) (abs y2))
	   (setq x1 xn) (setq y1 yn))
	  (t
	   (setq x2 xn) (setq y2 yn)))))


;;;; Converge

;;; find an x such that x = f(x)
;;;
(defun converge (fcn x eps)
  (let* ((x1               x )
	 (y1 (funcall fcn x1))
	 (x2              y1 )
	 (y2 (funcall fcn x2)))
    (declare (float x x1 x2 y1 y2))
    (do ((xn 0.0)
	 (yn 0.0))
	((< (abs (- x2 x1)) eps) x1)
      (declare (float xn yn))
      (setq xn (/ (- (* x1 y2) (* x2 y1))
		  (- (- x1 x2) (- y1 y2))))
      (setq yn (funcall fcn xn))
      (setq x1 x2 x2 xn
	    y1 y2 y2 yn))))

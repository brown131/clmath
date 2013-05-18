;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Quick Linear Regression Fit of Numerical Data

;;;  (c) Copyright Gerald Roylance 1983, 1984
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   only linear regression fit now

(in-package "CLMATH")

;;;; Quick Fit

;;; given a function of x and several parameters,
;;; produce code to do the fit.

;;; assume we are trying to fit the data to a function
;;; fctn(x) with parameters a0, a1, a2, a3
;;; That is, the model is:  fctn(x a0 a1 a2 a3).
;;;
;;; Because this is linear regression, we demand that
;;;   fctn(x a0 a1 a2 a3) = a0 + a1 * f1(x) + a2 * f2(x) + a3 * f3(x)
;;;
;;; To use REGRESS, we need an auxilliary function FIT-FCN(x i j)
;;;   that calculates the fj(x[i])
;;;
;;;    then for
;;;      j
;;;      0  f(x[i], 1, 0, 0, 0) *** this is actually never called
;;;      1  f(x[i], 0, 1, 0, 0)
;;;      2  f(x[i], 0, 0, 1, 0)
;;;      3  f(x[i], 0, 0, 0, 1)
;;;
;;; What we do is calculate these values for all i and j
;;; and store them in a table FIT-DATA.  The auxilliary
;;; function just looks them up in the table.

(defun fit-fcn (x i j)
  (aref x i j))


;;;; Main Procedure

;;; returns '(array-of-y-fit-data . parameters)

(defun fit (fctn x y nterms &optional (mode 0))
  (let* ((n        (array-dimension y 0))
	 (sigmay   (make-array n           :element-type 'float :initial-element 1.0))
	 (yfit     (make-array n           :element-type 'float :initial-element 0.0))
	 (a        (make-array (1+ nterms) :element-type 'float :initial-element 0.0))
	 (sigmaa   (make-array (1+ nterms) :element-type 'float :initial-element 0.0))
	 (r        (make-array (1+ nterms) :element-type 'float :initial-element 0.0))
	 (fit-data (make-array (list n (1+ nterms))
			       :element-type 'float :initial-element 0.0)))
    
    (do ((j 0 (1+ j))
	 (arguments nil))
	((> j nterms))

      ;; cons up an arg list
      (setq arguments nil)
      (do ((k 1 (1+ k)))
	  ((> k nterms))
	(push (cond ((= j k) 1.0)
		    (t       0.0))
	      arguments))
      (setq arguments (nreverse arguments))

      (do ((i 0 (1+ i)))
	  ((>= i n))
	(setf (aref fit-data i j)
	      (apply fctn
			     (aref x i)
			     0.0
			     arguments))
	))

    (REGRES #'fit-fcn fit-data y sigmay n nterms mode yfit a sigmaa r)

    (cons yfit (coerce a 'list))
    ))


;;;; Test
#|
(eval-when (eval)
|#
#|
  (defun fitest0 (x a0 a1 a2)
    (+ a0
	(* a1 (sin x))
	(* a2 (cos x))))

  (defun fitest ()
    (let ((x (make-array 20 :element-type 'float :initial-element 0.0))
	  (y (make-array 20 :element-type 'float :initial-element 0.0)))

      (do ((i 0 (1+ i)))
	  ((>= i 20))
	(setf (aref x i) (float i))
	(setf (aref y i) (fitest0 (float i) 0.0 0.27 0.73)))

      (fit #'fitest0 x y 2)))

  )
|#
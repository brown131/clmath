;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Bisection Search

;;;  (c) Copyright Gerald Roylance 1982
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   put in better termination criteria

(in-package "CLMATH")

;;;; Bisection Search

;;; This bisection search routine relies on a property of floating point
;;; numbers to terminate.  If there is a floating point number between
;;; s1 and s3, then continue; otherwise the midpoint(s1,s3) will be
;;; either s1 or s3.

;;; There are cases where this termination will take a long time:
;;; searching for a zero of f(x)=x will do hundreds (thousands for
;;; double precision) of function evaluations because floating point
;;; numbers are very accurate about zero.
;;;
;;; (bisection-search #'print -1.0 2.0)

;;; Thus there should be a termination criterion that includes absolute
;;; as well as a relative error.

;;; iterate until we hit a zero or s1=s2 or s2=s3

(defun bisection-search (fcn start stop)
  (let ((s1 start) (f1 (funcall fcn start))
	(s3 stop ) (f3 (funcall fcn stop)))
    (declare (float s1 s3 f1 f3))

    (cond ((= f1 0.0)				; fortuitous hits?
	   (return-from bisection-search s1))
	  ((= f3 0.0)
	   (return-from bisection-search s3)))

    ;; Make s1 be the low side and s3 the high side.
    ;; That is, f(s1) < 0 < f(s3)
    (if (> f1 0.0)
	(progn (rotatef s1 s3) (rotatef f1 f3)))
    
    ;; Make sure we have a bracket
    (if (< f3 0.0)
	(error "BISECTION-SEARCH: search points do not bracket"))

    (do ((s2 0.0)
	 (f2 0.0))
	(NIL)					; do forever
      (declare (float s2 f2))
      (setq s2 (/ (+ s1 s3) 2.0)		; calculate midpoint
	    f2 (funcall fcn s2))

      (cond ((= f2 0.0)				; hit the zero?
	     (return s2))
	    ((or (= s1 s2) (= s2 s3))		; iteration done? (ran out of room)
	     (return s2))
	    ((< f2 0.0)				; left side
	     (setq s1 s2 f1 f2))
	    (t					; must be right side
	     (setq s3 s2 f3 f2))))))

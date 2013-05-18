;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Horner's Rule

;;;  (c) Copyright Gerald Roylance 1982, 1984, 1987
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   Horner should use terms and termlists
;;;   poly and poly-expand should be flushed

(in-package "CLMATH")

;;;; Polynomial Expansion
  
;;; The simple version

(defun poly-expand (var coefs)
  (cond ((null coefs)       0.0)
	((null (cdr coefs)) (car coefs))
	(t
	 `(+ ,(car coefs) (* ,var ,(poly-expand var (cdr coefs)))))))

(defmacro poly (var . coefs)
  (poly-expand var coefs))


;;;; Horner Expansion

;;; The hairy version

;;; data structures

(defmacro horner-m (c e)   `(cons ,c ,e))
(defmacro horner-c (coeff) `(car ,coeff))
(defmacro horner-e (coeff) `(cdr ,coeff))

;;; take a list of coefficients (a0 a1 a2 a3 a4)
;;; and turn it into the horner form ((a0 0) (a1 1) ... )
;;;
(defun horner-coef-list-1 (as e)
  (if (null as)
      nil
      (cons (horner-m (car as) e)
	    (horner-coef-list-1 (cdr as) (1+ e)))))

(defun horner-coef-list (as)
  (horner-coef-list-1 as 0))

(defun horner-adjoin (x set)
  (if (member x set :test #'equal)
      set
      (sort (cons x set) #'>)))

;;; lookup a variable name in the dictionary
;;;   if it is not there, then add it
;;;
(defun horner-lookup (e dict)
  (if (assoc e (cdr dict) :test #'equal)
      (cdr (assoc e (cdr dict) :test #'equal))
      (let ((var (gensym (1- e))))
	(setf (cdr dict)
	      (cons (cons e var) (cdr dict)))
	var)))

(defun horner-dict-vars (dict)
  (mapcar #'cdr (cdr dict)))

;;;;

;;; take a horner list of coefficients and
;;; knock all terms with zero coefficients
;;;
(defun horner-knock-out (horner-list)
  (cond ((null horner-list) nil)
	((zerop (horner-c (car horner-list)))
	 (horner-knock-out (cdr horner-list)))
	(t
	 (cons (car horner-list)
	       (horner-knock-out (cdr horner-list))))))

;;; find a list of all the differences in exponents
;;;
(defun horner-deltas (e horner-list acc)
  (cond ((null horner-list) acc)		; no coefficients, no work
	((= e (horner-e (car horner-list)))	; no difference, so no delta
	 (horner-deltas e (cdr horner-list) acc))
	(t					; we have a delta
	 (let ((en (horner-e (car horner-list))))
	   (horner-deltas en
			  (cdr horner-list)
			  (horner-adjoin (- en e) acc)
			  )))))


;;;; Construct LET* Variable List and the LET* Body

;;; make a variable list for a LET* form
;;;
(defun horner-var-maker (exp d-list v-list dict)
  (if (null d-list)				; if all done
      v-list					; return the variable and dictionary
      (let ((d (car d-list)))
	(if (= d 1)
	    (horner-var-maker exp
			      (cdr d-list)
			      (cons `(,(horner-lookup 1 dict) ,exp) v-list)
			      dict)
	    (let* ((d0 (floor d 2))
		   (d1 (- d d0)))
	      (horner-var-maker exp
				(horner-adjoin d1 (horner-adjoin d0 (cdr d-list)))
				(cons `(,(horner-lookup d dict)
					(* ,(horner-lookup d1 dict)
					   ,(horner-lookup d0 dict)))
				      v-list)
				dict))))
      ))

;;; return an expression that when multiplied by x^e
;;;   evaluates the polynomial specified by coeffs
;;;
(defun horner-body (e coeffs dict)
  (if (null coeffs)
      0.0					; empty expression
      (let ((coef (car coeffs)))
	(if (= (horner-e coef) e)		; time to add in a term?
	    (if (null (cdr coeffs))		; eg, 7x^12 and e=12
		`,(horner-c coef)
		`(+ ,(horner-c coef)
		    ,(horner-body e (cdr coeffs) dict)))
	    `(* ,(horner-lookup (- (horner-e coef) e) dict)
		,(horner-body (horner-e coef) coeffs dict)))
	)))


;;;; Put it all together

(defmacro horner-polynomial (exp . coeffs)
  (let* ((c-list (horner-knock-out (horner-coef-list coeffs)))
	 (diffs  (horner-deltas 0 c-list '()))
	 (dict   (list 'dict))
	 (maker  (horner-var-maker exp diffs '() dict)))
    `(LET* ,maker
       ;; (DECLARE (FLOAT ,@(horner-dict-vars dict)))
       ,(horner-body 0 c-list dict))))

(defmacro horner-hack (exp . coeffs)
  `(horner-polynomial ,exp ,@ coeffs))

;;; (horner-polynomial (log x) 1.0 0.0 3.5 0.0 2.3 0.0 0.0 7.5)

;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Multiple Linear Regression

;;; Reference
;;;   Philip R. Bevington
;;;   Data Reduction and Error Analysis for the Physical Sciences
;;;   McGraw-Hill, 1969

;;; This code is derived from the FORTRAN code in Bevington.
;;; Therefore (c) Copyright McGraw-Hill 1969
;;; *** Permission has not been sought to distribute this code. ***
;;;
;;; Modifications beyond Bevington
;;;   (c) Copyright Gerald Roylance 1983, 1986
;;;       All Rights Reserved.
;;;    more block structure
;;;    zero base x,y arrays
;;;    make a[0] and array[0] have meaning
;;;    completely redo yfit calculation
;;;    clearing stupid arrays

;;; Bugs and Fixes
;;;   Weight array should be separate.  Statistics, too.
;;;   Have to return the statistics and the coefficents
;;;   Blows out determinant -- want a different matrix method?
;;;   expand variable names
;;;   use multiplication where appropriate
;;;   break into smaller procedures
;;;   don't pass in so many arguments
;;;     mode and sigmay --> weight
;;;     npts --> get from dimension of y
;;;     nterms --> keep so different FCTNs can be used
;;;   calculate an array of the fj(x[i]) instead of
;;;     calling the function all the time
;;;     pass that array in as an argument

(in-package "CLMATH")

;;;; Multiple Linear Regression

;;; Make a multiple linear regression fit to data with a specified
;;;   function which is linear in coefficients

;;; description of parameters
;;;	x	array[0:npts-1] of data points for independent variable
;;;	y	array[0:npts-1] of data points for dependent variable
;;;	sigmay	array[0:npts-1] of standard deviations for y data points
;;;	npts	number of pairs of data points
;;;	nterms	number of coefficients
;;;	mode	determines method of weighting least-squares fit
;;;		+1 -- (instrumental) weight(i) = 1/sigmay(i)**2
;;;		 0 -- (no weighting) weight(i) = 1
;;;		-1 -- (statistical)  weight(i) = 1/y(i)
;;;	yfit	array[0:npts-1] of calculated values of y
;;;	a	array[0:nterms] of coefficients
;;;	sigmaa	array[0:nterms] of deviations
;;;	r	array[0:nterms] of linear correlation coefficients
;;;*	rmul	multiple linear correlation coefficient
;;;*	chisqr	reduced chi square for fit
;;;*	ftest	value of f for test of fit
;;;
;;;* -- cannot return thru arguments because MACLISP is call by VALUE

;;; subroutines required
;;;   FCTN(X, I, J)
;;;	evaluates the function for the jth term and ith data point


;;;; YFIT calculation

;;; given FCTN (later to be an array)
;;; and the coefficients

;;; this is a matrix multiply
;;;  yfit = fctn * a

(defun regression-yfit (fctn a yfit x npts nterms)
  (declare (type (array float (*)) a yfit))
  
  (do ((i 0 (1+ i)))
      ((>= i npts))
    (declare (fixnum i))
    
    (do ((j 1 (1+ j))				; jimmied so fctn(x i 0) never called
	 (z (aref a 0)))			; jimmied
	((> j nterms)
	 (setf (aref yfit i) z)
	 NIL)
      (declare (fixnum j)
	       (float z))
      
      (setq z (+ z
		 (* (aref a j)
		    (funcall fctn x i j))))
      )
    ))


;;;; Multiple Linear Regression

(defun REGRES (FCTN x y sigmay npts nterms mode yfit a sigmaa r)
  (declare (type (array float (*)) y yfit sigmay)
	   (type (array float (*)) a sigmaa r)
	   (fixnum npts nterms mode))
  (let ((weight (make-array npts
			    :element-type 'float :initial-element 0.0))
	(xmean  (make-array (1+ nterms)
			    :element-type 'float :initial-element 0.0))
	(sigmax (make-array (1+ nterms)
			    :element-type 'float :initial-element 0.0))
	(array  (make-array (list (1+ nterms) (1+ nterms))
			    :element-type 'float :initial-element 0.0))
	(fnpts  (float npts))
	(free1  (float (- npts 1)))
	(freej  (float nterms))
	(freen  (float (- npts nterms 1)))
	(a0     0.0)
	(sigma0 0.0)
	(sum    0.0)
	(det    0.0)
	(wmean  0.0)
	(ymean  0.0)
	(sigma  0.0)
	(chisq  0.0)
	(varnce 0.0)
	(rmul   0.0)
	(chisqr 0.0)
	(ftest  0.0))
    (declare (type (simple-array float   (*)) weight xmean sigmax)
	     (type (simple-array float (* *)) array)
	     (float sum det wmean ymean sigma chisq varnce
		    a0 sigma0 rmul chisqr ftest
		    fnpts free1 freej freen))
    
    ;; initialize sums and arrays
    (do ((j 0 (1+ j)))
	((> j nterms))
      (declare (fixnum j))
      (setf (aref xmean  j) 0.0)
      (setf (aref sigmax j) 0.0)
      (setf (aref r      j) 0.0)
      (setf (aref sigmaa j) 0.0)
      (do ((k 0 (1+ k)))
	  ((> k nterms))
	(declare (fixnum k))
	(setf (aref array j k) 0.0)))
    
    ;; accumulate weighted sums
    (do ((i 0 (1+ i)))
	((>= i npts))
      (declare (fixnum i))
      (setf (aref weight i)
	    (cond ((< mode 0) (cond ((= 0.0 (aref y i)) 1.0)
				    (t (/ (abs (aref y i))))))
		  ((= mode 0) 1.0)
		  ((> mode 0) (/ (expt (aref sigmay i) 2)))))
      (setq sum   (+ sum   (aref weight i)))
      (setq ymean (+ ymean (* (aref weight i)
			      (aref      y i))))
      (do ((j 1 (1+ j)))
	  ((> j nterms))
	(declare (fixnum j))
	(setf (aref xmean j)
	      (+ (aref xmean j)
		 (* (aref weight i)
		    (FUNCALL FCTN x i j))))))
    (setq ymean (/ ymean sum))
    (do ((j 1 (1+ j)))
	((> j nterms))
      (declare (fixnum j))
      (setf (aref xmean j)
	    (/ (aref xmean j) sum)))
    (setq wmean (/ sum fnpts))
    (do ((i 0 (1+ i)))
	((>= i npts))
      (declare (fixnum i))
      (setf (aref weight i)
	    (/ (aref weight i) wmean)))
    
    ;; accumulate matrices r and array
    (do ((i 0 (1+ i)))
	((>= i npts))
      (declare (fixnum i))
      (setq sigma (+ sigma
		     (* (aref weight i)
			(expt (- (aref y i) ymean) 2))))
      (do ((j 1 (1+ j)))
	  ((> j nterms))
	(declare (fixnum j))
	(setf (aref sigmax j)
	      (+ (aref sigmax j)
		 (* (aref weight i)
		    (expt (- (FUNCALL FCTN x i j)
			     (aref xmean j))
			  2))))
	(setf (aref r j)
	      (+ (aref r j)
		 (* (aref weight i)
		    (- (FUNCALL FCTN x i j) (aref xmean j))
		    (- (aref y i) ymean))))
	(do ((k 1 (1+ k)))
	    ((> k j))
	  (declare (fixnum k))
	  (setf (aref array j k)
		(+ (aref array j k)
		   (* (aref weight i)
		      (- (FUNCALL FCTN x i j) (aref xmean j))
		      (- (FUNCALL FCTN x i k) (aref xmean k))
		      ))))))
    (setq sigma (sqrt (/ sigma free1)))
    (do ((j 1 (1+ j)))
	((> j nterms))
      (declare (fixnum j))
      (setf (aref sigmax j)
	    (sqrt (/ (aref sigmax j) free1)))
      (setf (aref      r j)
	    (/ (aref      r j)
	       free1
	       (aref sigmax j)
	       sigma))
      (do ((k 1 (1+ k)))
	  ((> k j))
	(declare (fixnum k))
	(setf (aref array j k)
	      (/ (aref array j k)
		 free1
		 (aref sigmax j)
		 (aref sigmax k)))
	(setf (aref array k j)
	      (aref array j k))))
    
    ;; invert symmetric matrix
    (setf (aref array 0 0) 1.0)
    (setq det (matrix-inverse array))
    
    (cond ((= det 0.0)
	   (setq a0     0.0)
	   (setq sigma0 0.0)
	   (setq rmul   0.0)
	   (setq chisqr 0.0)
	   (setq ftest  0.0)
	   NIL)
	  (t
	   ;; calculate coefficients, fit and chi square
	   (do ((j 1 (1+ j)))			; matrix multipy ARRAY * R
	       ((> j nterms))			; 1-based arrays!
	     (declare (fixnum j))		;   if r[0] = 0, it doesn't matter
	     (do ((k   1 (1+ k))		;   actually a[0] = r[0]
		  (sum 0.0))
		 ((> k nterms)
		  (setf (aref a j) sum)
		  nil)
	       (declare (fixnum k) (float sum))
	       (setq sum (+ sum (* (aref array j k)
				   (aref r k))))))
	   (do ((j 1 (1+ j)))
	       ((> j nterms))
	     (declare (fixnum j))
	     (setf (aref a j)
		   (* (aref a j)
		      (/ sigma (aref sigmax j)))))
	   
	   (setq a0 ymean)
	   (do ((j 1 (1+ j)))
	       ((> j nterms))
	     (declare (fixnum j))
	     (setq a0 (- a0 (* (aref     a j)
			       (aref xmean j)))))
	   (setf (aref a 0) a0)			; make a(0) meaningful
	   
	   ;; calculate the fitted y values
	   (regression-yfit fctn a yfit x npts nterms)
	   
	   (do ((i 0 (1+ i)))
	       ((>= i npts))
	     (declare (fixnum i))
	     (setq chisq (+ chisq
			    (* (aref weight i)
			       (expt (- (aref    y i)
					(aref yfit i))
				     2)))))
	   (setq chisqr (* chisq (/ wmean freen)))
	   
	   ;; calculate uncertainties
	   (setq varnce (if (= mode 0) chisqr (/ wmean)))
	   (do ((j 1 (1+ j)))
	       ((> j nterms))
	     (declare (fixnum j))
	     (setf (aref sigmaa j)
		   (/ (* (aref array j j) varnce)
		      free1
		      (expt (aref sigmax j) 2)))
	     (setf (aref sigmaa j)
		   (sqrt (aref sigmaa j)))
	     (setq rmul
		   (+ rmul
		      (* (aref a j)
			 (aref r j)
			 (/ (aref sigmax j)
			    sigma)))))
	   (setq ftest (/ (/ rmul freej)
			  (/ (- 1.0 rmul) freen)))
	   (setq rmul (sqrt rmul))
	   (setq sigma0 (/ varnce fnpts))
	   (do ((j 1 (1+ j)))
	       ((> j nterms))
	     (declare (fixnum j))
	     (do ((k 1 (1+ k)))
		 ((> k nterms))
	       (declare (fixnum k))
	       (setq sigma0
		     (+ sigma0
			(/ (* varnce
			      (aref xmean j)
			      (aref xmean k)
			      (aref array j k))
			   free1
			   (aref sigmax j)
			   (aref sigmax k))))
	       ))
	   (setq sigma0 (sqrt sigma0))
	   
	   (setf (aref sigmaa 0) sigma0)
	   (setf (aref      r 0)   rmul)
	   
	   (LIST 'REGRES
		 'DETERMINANT DET
		 'CHISQ       CHISQ
		 'CHISQR      CHISQR
		 'FTEST       FTEST
		 'RMUL        RMUL)
	   ))
    ))


;;;; Tests

#|

(defun fctn (x i j)
  (declare (type (array float (*)) x)
	   (fixnum i j))
  (expt (aref x i) j))

(setq   npts 9)
(setq nterms 2)

(setq      x (make-array
	       npts
	       :element-type     'float
	       :initial-contents '( 0.0  1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0)))
(setq      y (make-array
	       npts
	       :element-type     'float
	       :initial-contents '( 0.0  1.1  2.4  3.9  5.6  7.5  9.6 11.9 14.4)))
(setq sigmay (make-array
	       npts
	       :element-type     'float
	       :initial-contents '( 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0)))
(setq   yfit (make-array npts :element-type 'float :initial-element 0.0))
(setq   mode 0)
(setq      a (make-array (1+ nterms) :element-type 'float :initial-element 0.0))
(setq sigmaa (make-array (1+ nterms) :element-type 'float :initial-element 0.0))
(setq      r (make-array (1+ nterms) :element-type 'float :initial-element 0.0))

(defun test ()
  (declare (special x y sigmay npts nterms mode yfit a sigmaa r
		    rmul chisqr ftest))
  (let ((result
	  (regres 'fctn x y sigmay npts nterms mode yfit a sigmaa r)))
    
    (format t "~%")
    (format t "~% DETERMINANT = ~F" (getf (cdr result) 'determinant))
    (format t "~% CHISQ         ~F" (getf (cdr result) 'chisq))
    (format t "~% CHISQR        ~F" (getf (cdr result) 'chisqr))
    (format t "~% FTEST         ~F" (getf (cdr result) 'ftest))
    (format t "~% RMUL          ~F" (getf (cdr result) 'rmul))

    (format t "~%")
    (format t "~%   A        SIGMAA   R")
    (do ((j 0 (1+ j)))
	((> j nterms))
      (format t "~%~3,1,8$ ~3,1,8$ ~3,1,8$"
	      (aref      a j)
	      (aref sigmaa j)
	      (aref      r j)))

    (format t "~%")
    (format t "~%   X        Y      YFIT")
    (do ((i 0 (1+ i)))
	((>= i npts))
      (format t "~%~3,1,8$ ~3,1,8$ ~3,1,8$"
	      (aref    x i)
	      (aref    y i)
	      (aref yfit i)))

    ))

|#

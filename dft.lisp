;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Discrete Fourier Transform

;;;  (c) Copyright Gerald Roylance 1983, 1984
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Reference
;;;   Alan Oppenheim and Ronald Schafer, 
;;;     Digital Signal Processing
;;;     Prentice Hall, 1975

;;; Bugs and Fixes
;;;   Can do 4x and 8x decompositions instead of just 2x
;;;   Precompute tables for bit reverse
;;;   Should have some generalized complex arithmetic

(in-package "CLMATH")

;;; Transform conventions

;;; time -> frequency
;;;   bit reverse the time array
;;;   output is in correct order

;;; frequency -> time
;;;    freq in correct order
;;;    time is output in bit reversed order


;;;; Initialization Code

;;; Return Precomputed tables

(DEFUN    DFT-PRECOMPUTED-TABLEP (STRUCTURE)
  (EQ (CAR STRUCTURE)  'DFT-PRECOMPUTED))
(DEFUN    DFT-PRECOMPUTED-MAKE   ()          (LIST 'DFT-PRECOMPUTED))
(DEFMACRO DFT-PRECOMPUTED-SIZE   (STRUCTURE) `(GETF (CDR ,STRUCTURE) 'SIZE))
(DEFMACRO DFT-PRECOMPUTED-W-REAL (STRUCTURE) `(GETF (CDR ,STRUCTURE) 'W-REAL))
(DEFMACRO DFT-PRECOMPUTED-W-IMAG (STRUCTURE) `(GETF (CDR ,STRUCTURE) 'W-IMAG))

(DEFUN DFT-INIT (LOG-SIZE)
  (LET* ((SIZE   (ASH 1 LOG-SIZE))
	 (SIZE2  (ASH SIZE -1))
	 (STRUCT (DFT-PRECOMPUTED-MAKE))
	 (W-REAL (make-array SIZE2 :element-type 'FLOAT))
	 (W-IMAG (make-array SIZE2 :element-type 'FLOAT)))
    (DECLARE (FIXNUM SIZE SIZE2)
	     (type (array float (*)) w-real w-imag))

    (setf (DFT-PRECOMPUTED-SIZE   STRUCT) LOG-SIZE)
    (setf (DFT-PRECOMPUTED-W-REAL STRUCT) W-REAL)
    (setf (DFT-PRECOMPUTED-W-IMAG STRUCT) W-IMAG)

    (DO ((I 0 (1+ I))
	 (ANGLE (/ (* 2.0 3.1415926) (float SIZE))))
	((>= I SIZE2))
      (DECLARE (FIXNUM I)
	       (FLOAT ANGLE))

      (setf (aref W-REAL I) (COS (* (float I) ANGLE)))
      (setf (aref W-IMAG I) (SIN (* (float I) ANGLE)))
      )
    STRUCT))


;;;; Forward transform

;;; figure 6-10 in Oppenheim
;;; time is bit-reversed, frequency is in order

(defun dft-610 (x-real x-imag tables)
  (declare (type (array float (*)) x-real x-imag))
  (let* ((ln       (dft-precomputed-size tables))
	 (n        (ash 1 ln))
	 (n-over-2 (ash n -1))
	 (w-n-real (dft-precomputed-w-real tables))
	 (w-n-imag (dft-precomputed-w-imag tables)))
    (declare (fixnum ln n n-over-2)
	     (type (array float (*)) w-n-real w-n-imag))
    (do ((i 1 (1+ i))				; each stage
	 (butter-dis       1  (*  butter-dis 2))
	 (separation       2  (*  separation 2))
	 (repetition  (floor n 2) (floor repetition 2))
	 (w-real 0.0) (w-imag 0.0)
	 (x0     0.0) (x1     0.0))
	((> i ln))
      (declare (fixnum i butter-dis separation repetition)
	       (float w-real w-imag x0 x1))
      (do ((k 0 (+ k repetition))		; each W factor
	   (l 0 (1+ l)))
	  ((>= k n-over-2))
	(declare (fixnum k l))
	(setq w-real (aref w-n-real k))
	(setq w-imag (aref w-n-imag k))
	(do ((j l (+ j separation)))		; each butterfly
	    ((>= j n))
	  (declare (fixnum j))
	  ;; multiply x[j+butter-dis] by W
	  (setq x0 (aref x-real (+ j butter-dis)))
	  (setq x1 (aref x-imag (+ j butter-dis)))
	  (setf (aref x-real (+ j butter-dis))
		(- (* w-real x0) (* w-imag x1)))
	  (setf (aref x-imag (+ j butter-dis))
		(+ (* w-real x1) (* w-imag x0)))

	  ;; butterfly j and j+butter-dis
	  (setq x0 (aref x-real (+ j 0)))
	  (setq x1 (aref x-real (+ j butter-dis)))
	  (setf (aref x-real (+ j          0)) (+ x0 x1))
	  (setf (aref x-real (+ j butter-dis)) (- x0 x1))

	  (setq x0 (aref x-imag (+ j 0)))
	  (setq x1 (aref x-imag (+ j butter-dis)))
	  (setf (aref x-imag (+ j          0)) (+ x0 x1))
	  (setf (aref x-imag (+ j butter-dis)) (- x0 x1))
	  )))))


;;;; Reverse transform

;;; figure 6-18 in Oppenheim
;;; frequency is in order, time is bit-reversed

(defun dft-618 (x-real x-imag tables)
  (declare (type (array float (*)) x-real x-imag))
  (let* ((ln       (dft-precomputed-size tables))
	 (n        (ash 1 ln))
	 (n-over-2 (ash n -1))
	 (w-n-real (dft-precomputed-w-real tables))
	 (w-n-imag (dft-precomputed-w-imag tables)))
    (declare (fixnum ln n n-over-2)
	     (type (array float (*)) w-n-real w-n-imag))
    (do ((i 1 (1+ i))				; each stage
	 (butter-dis (floor n 2) (floor butter-dis 2))
	 (separation        n (floor separation 2))
	 (repetition        1 (*  repetition 2))
	 (w-real 0.0) (w-imag 0.0)
	 (x0     0.0) (x1     0.0))
	((> i ln))
      (declare (fixnum i butter-dis separation repetition)
	       (float w-real w-imag x0 x1)) 
      (do ((k 0 (+ k repetition))		; each W factor
	   (l 0 (1+ l)))
	  ((>= k n-over-2))
	(declare (fixnum k l))
	;; hacked for W**-k
	(setq w-real     (aref w-n-real k))
	(setq w-imag (- (aref w-n-imag k)))
	(do ((j l (+ j separation)))		; each butterfly
	    ((>= j n))
	  (declare (fixnum j))
	  ;; butterfly j and j+butter-dis
	  (setq x0 (aref x-real (+ j 0)))
	  (setq x1 (aref x-real (+ j butter-dis)))
	  (setf (aref x-real (+ j          0)) (+ x0 x1))
	  (setf (aref x-real (+ j butter-dis)) (- x0 x1))

	  (setq x0 (aref x-imag (+ j 0)))
	  (setq x1 (aref x-imag (+ j butter-dis)))
	  (setf (aref x-imag (+ j          0)) (+ x0 x1))
	  (setf (aref x-imag (+ j butter-dis)) (- x0 x1))
	  
	  ;; multiply x[j+butter-dis] by W
	  (setq x0 (aref x-real (+ j butter-dis)))
	  (setq x1 (aref x-imag (+ j butter-dis)))
	  (setf (aref x-real (+ j butter-dis))
		(- (* w-real x0) (* w-imag x1)))
	  (setf (aref x-imag (+ j butter-dis))
		(+ (* w-real x1) (* w-imag x0)))
	  )))
    (do ((i 0 (1+ i)))
	((>= i n))
      (declare (fixnum i))
      (setf (aref x-real i) (/ (aref x-real i) (float n)))
      (setf (aref x-imag i) (/ (aref x-imag i) (float n))))
    ))


;;;; Reverse Bits

#|+ignore
(DECLARE (FIXNUM (DFT-REVERSE-BITS-16 FIXNUM)))
|#

(defun dft-reverse-bits-16 (n)
  (let ((r     n)
	(mask1 #2r 0101010101010101)
	(mask2 #2r 0011001100110011)
	(mask4 #2r 0000111100001111)
	(mask8 #2r 0000000011111111))
    (declare (type (unsigned-byte 16) r mask1 mask2 mask4 mask8))
    (setq r (logior (logand (ash r -1) mask1) (ash (logand r mask1)  1)))
    (setq r (logior (logand (ash r -2) mask2) (ash (logand r mask2)  2)))
    (setq r (logior (logand (ash r -4) mask4) (ash (logand r mask4)  4)))
    (setq r (logior (logand (ash r -8) mask8) (ash (logand r mask8)  8)))
    ))

(defun dft-reverse-array (array tables)
  (declare (type (array float (*)) array))
  (let* ((log-size (dft-precomputed-size tables))
	 (size     (ash 1 log-size)))
    (declare (fixnum log-size size))
    (if (or (< log-size  0.)
	    (> log-size 16.)
	    (> (array-dimension array 0) size))
	(error "DFT-REVERSE-ARRAY given bad args"))
    (do ((i      0 (1+ i))
	 (j      0)
	 (shift  (- log-size 16.))
	 (temp 0.0))
	((>= i size))
      (declare (type fixnum i j shift)
	       (type float temp))
      (setq j (ash (dft-reverse-bits-16 i) shift))
      (cond ((< i j)
	     (setq temp                         (aref array i))
	     (setf (aref array i) (aref array j))
	     (setf (aref array j) temp))))))


;;;; Forward and Reverse Transforms

;;; So you don't have to worry about bit-reversing ...

(defun dft-forward (x-real x-imag tables)
  (dft-reverse-array x-real tables)
  (dft-reverse-array x-imag tables)
  (dft-610 x-real x-imag tables))

(defun dft-reverse (x-real x-imag tables)
  (dft-618 x-real x-imag tables)
  (dft-reverse-array x-real tables)
  (dft-reverse-array x-imag tables))


;;;; Test
#|
(EVAL-WHEN (EVAL)
|#
#|
  (DEFUN DFT-TEST-PRINT-ARRAY (ARRAY)
    (terpri)
    (DO ((I 0 (1+ I))
	 (N (array-dimension ARRAY 0)))
	((>= I N))
      (DECLARE (FIXNUM I N))
      (FORMAT T "~3,1,8$" (aref ARRAY I))))

  (DEFUN DFT-TEST-PRINT (TITLE X Y)
    (PRINT TITLE)
    (DFT-TEST-PRINT-ARRAY X)
    (DFT-TEST-PRINT-ARRAY Y))

  (DEFUN DFT-TEST ()
    (LET* ((LOG-SIZE 3)
	   (N        (expt 2 LOG-SIZE))
	   (TABLES   (DFT-INIT LOG-SIZE))
	   (X-REAL   (make-array n :element-type 'FLOAT))
	   (X-IMAG   (make-array n :element-type 'FLOAT)))
      (DECLARE (FIXNUM LOG-SIZE N))

      (DO ((I 0 (1+ I)))
	  ((>= I N))
	(DECLARE (FIXNUM I))
	(setf (aref X-REAL I)
	      (COS (/ (*  4.0 3.1415926 (float I)) (float N))))
	(setf (aref X-IMAG I) 0.0))

      (DFT-TEST-PRINT 'TIME-DOMAIN X-REAL X-IMAG)
      (DFT-FORWARD X-REAL X-IMAG TABLES)
      (DFT-TEST-PRINT 'FREQUENCY   X-REAL X-IMAG)
      (DFT-REVERSE X-REAL X-IMAG TABLES)
      (DFT-TEST-PRINT 'TIME-DOMAIN X-REAL X-IMAG)
      ))
  )
|#

;;;; Speed Tests
#|
(EVAL-WHEN (EVAL)
|# 
#|
  (DEFUN DFT-SPEED-TEST (DFT-LOG-SIZE)
    (LET ((DFT-TABLES (DFT-INIT DFT-LOG-SIZE))
	  (X-REAL     (make-array (ASH 1 DFT-LOG-SIZE) :element-type 'FLOAT))
	  (X-IMAG     (make-array (ASH 1 DFT-LOG-SIZE) :element-type 'FLOAT))
	  (TIME0 0)
	  (TIME1 0)
	  (TIME2 0)
	  (TIME3 0))
      (DECLARE (FIXNUM I TIME0 TIME1 TIME2 TIME3))

      (DO ((I 0 (1+ I)))
	  ((>= I (ASH 1 DFT-LOG-SIZE)))
	(DECLARE (FIXNUM I))
	(setf (aref X-REAL I) 1.0)
	(setf (aref X-IMAG I) 0.0))
      
      (SETQ TIME0 (RUNTIME))
      (DFT-REVERSE-ARRAY X-REAL DFT-TABLES)
      (SETQ TIME1 (RUNTIME))
      (DFT-610 X-REAL X-IMAG DFT-TABLES)
      (SETQ TIME2 (RUNTIME))
      (SETQ TIME3 (RUNTIME))

      (FORMAT T "~%  Bit reverse time ~8D microseconds" (- TIME1 TIME0))
      (FORMAT T "~%          DFT time ~8D microseconds" (- TIME2 TIME1))
      (FORMAT T "~%        Total time ~8D microseconds" (- TIME2 TIME0))
      (FORMAT T "~% Call to (RUNTIME) ~8D microseconds" (- TIME3 TIME2))
      (FORMAT T "~%             Speed ~8F N LOG2(N)"
	      (/ (float (- TIME2 TIME1))
		   (float DFT-LOG-SIZE)
		   (float (ASH 1 DFT-LOG-SIZE))))
      ))

  (DEFUN MPY-TIME (N)
    (DO ((TIME (RUNTIME))
	 (I    0 (1+ I))
	 (X    3.14159))
	((>= I N)
	 (/ (float (- (RUNTIME) TIME))
	      (float N)))
      (DECLARE (FIXNUM I TIME) (FLOAT X))
      (* X X))
    )

  (DEFUN ADD-TIME (N)
    (DO ((TIME (RUNTIME))
	 (I    0 (1+ I))
	 (X    3.14159))
	((>= I N)
	 (/ (float (- (RUNTIME) TIME))
	      (float N)))
      (DECLARE (FIXNUM I TIME) (FLOAT X))
      (+ X X)))
  )
|#
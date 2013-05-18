;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Matrix Routines

;;; References
;;;   David G. Luenberger,
;;;     Introduction to Linear and Nonlinear Programming,
;;;     Addison Wesley, Reading, MA, 1965.

;;;  (c) Copyright Gerald Roylance 1983, 1985, 1986
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   The LU decomposition stuff does not pivot
;;;     should be LUP system
;;;   All this effort to type coerce a matrix multiply
;;;     isn't worth it.
;;;   Is it possible to recover sparseness
;;;     at the elimination step, check if coeff=0
;;;     also, might count number of zeros in pivot column and
;;;      remember nonzero indices in an array.  if lots of zeros,
;;;      then do it sparsely
;;;   How about integer matrix multiplies?
;;;   ROW-COL and COL-ROW matrix multiples don't exist
;;;   Multiply routines should check array dimensions

(in-package "CLMATH")

;;;; Matrix Printing Routines

(defun matrix-print-vector (a)
  (let ((n (array-dimension a 0)))
    (declare (fixnum n))
    (dotimes (i n)
      (declare (fixnum i))
      (format t "~%~5,1,10$ " (aref a i)))
    a))

(defun matrix-print-matrix (a)
  (let ((n (array-dimension a 0))
	(m (array-dimension a 1)))
    (declare (fixnum m n))
    (dotimes (i n)
      (declare (fixnum i))
      (terpri)
      (dotimes (j m)
	(declare (fixnum j))
	(format t "~5,1,10$ " (aref a i j))))
     (terpri)
     a))

(defun matrix-print (a)
  (let ((dims (array-rank a)))
    (declare (fixnum dims))
    (format t "~%; Matrix ~A =~%" a)
    (case dims
      (1 (matrix-print-vector a))
      (2 (matrix-print-matrix a))
      (otherwise
	(error "MATRIX-PRINT given more than two dimensions")))))


;;;; Matrix Stuff

;;; Generate an NxN identity matrix
;;;
(defun matrix-identity (n)
  (do ((matrix (make-array (list n n)
			   :element-type 'float :initial-element 0.0))
       (i 0 (1+ i)))
      ((>= i n) matrix)
    (declare (fixnum i))
    (setf (aref matrix i i) 1.0)))

;;; Copy an arbitrary matrix
;;;
(defun matrix-copy (a)
  (let ((copy (make-array (array-dimensions a)
			  :element-type (array-element-type a))))
    (case (array-rank a)
      (1 (replace copy a))
      (2 (dotimes (i (array-dimension a 0))
	   (dotimes (j (array-dimension a 1))
	     (setf (aref copy i j) (aref a i j)))))
      (otherwise
	(error "Matrix-Copy lost ~s" a)))
    copy))

;;; Copy just the diagonal elements of a matrix
;;;
(defun matrix-diagonal (a)
  (let* ((n      (array-dimension a 0))
	 (matrix (make-array (list n n)
			     :element-type 'float :initial-element 0.0)))
    (declare (fixnum n))
    (do ((i 0 (1+ i)))
	((>= i n) matrix)
      (declare (fixnum i))
      (setf (aref matrix i i)
	    (aref a      i i)))))


;;;; Matrix Addition and Subtraction

;;; Add 2 Matrices
;;;
(defun matrix-add (a b)
  (let* ((n (array-dimension a 0))
	 (m (array-dimension a 1))
	 (c (make-array (list n m)
			:element-type 'float :initial-element 0.0)))
    (declare (fixnum m n))
    (dotimes (i n)
      (declare (fixnum i))
      (dotimes (j m)
	(declare (fixnum j))
	(setf (aref c i j)
	      (+ (aref a i j)
		 (aref b i j)))))
     c))

;;; Subtract 2 Matrices
;;;
(defun matrix-sub (a b)
  (let* ((n (array-dimension a 0))
	 (m (array-dimension a 1))
	 (c (make-array (list n m)
			:element-type 'float :initial-element 0.0)))
    (declare (fixnum m n))
    (dotimes (i n)
      (declare (fixnum i))
      (dotimes (j m)
	(declare (fixnum j))
	(setf (aref c i j)
	      (- (aref a i j)
		 (aref b i j)))))
    c))


;;;; Solving Triangular Matrix Problems

;;; These operations are order n**2

;;; Presume all diagonal elements are nonzero (determinant .ne. 0)
;;; There is a flag to treat the main diagonal as ones.
;;; That lets us put an LU decomposition into one matrix

;;; solve Ax=b, put result in b
;;;
(defun matrix-solve-triangle-lower (A b &optional (MD-unity NIL))
  (let ((n (array-dimension A 0)))
    (declare (fixnum n))
    (dotimes (i n)				; do for each row
      (declare (fixnum i))
      
      (if (not MD-unity)			; main diagonal not unity?
	  (setf (aref b i)			; solution for xi
		(/ (aref b i)
		   (aref A i i))))
      
      (do ((j (1+ i) (1+ j)))			; do for all the x[i] we know
	  ((>= j n))
	(declare (fixnum j))
	(setf (aref b j)
	      (- (aref b j)
		 (* (aref A j i)
		    (aref b i))))))
    b))

;;; solve Ax=b, put result in b
;;;
(defun matrix-solve-triangle-upper (A b &optional (MD-unity NIL))
  (let ((n (array-dimension A 0)))
    (declare (fixnum n))
    (do ((i (1- n) (1- i)))			; do for each row
	((< i 0))
      (declare (fixnum i))
      
      (if (not MD-unity)			; main diagonal not unity?
	  (setf (aref b i)			; solution for xi
		(/ (aref b i)
		   (aref A i i))))
      
      (do ((j (1- i) (1- j)))			; do for all the x[i] we know
	  ((< j 0))
	(declare (fixnum j))
	(setf (aref b j)
	      (- (aref b j)
		 (* (aref A j i)
		    (aref b i))))))
    b))


;;;; What happens in Matrix Inversion

;;; let A be the matrix we are interested in
;;;    A = A
;;; premultiply both sides by a matrix M1 that makes column 1 look like I
;;;    (M1 A) = M1 A
;;; now use M2 to fix up column 2, etc
;;;    (Mn ... M3 M2 M1 A) = (Mn ... M3 M2 M1) A
;;; we turned the lhs into I, so (Mn ... M3 M2 M1) must be A'

;;; each Mi looks like                and M(i-1) A
;;;     [ 1 ... 0   ... 0   0   0]    [ 1 ... a0i ... *   *   *]
;;;     [ 0 ... mii ... 0   0   0]    [ 0 ... aii ... *   *   *]
;;;     [ 0 ... m** ... 1   0   0]    [ 0 ... a*i ... *   *   *]
;;;     [ 0 ... m** ... 0   1   0]    [ 0 ... a*i ... *   *   *]
;;;     [ 0 ... m** ... 0   0   1]    [ 0 ... a*i ... *   *   *]
;;;
;;; where the element mji = 0         if i<j
;;;                       =    1/aii  if i=j
;;;                       = -aji/aii  if i>j
;;;
;;; The matrix inversion algorithm must compute both (M A) and M
;;; we can share the space (put M and (M A) in the same matrix)
;;; (M A) only needs the lower right corner; the rest are 0, 1, and m**
;;; M needs the other 3 corners
;;;  so what we have to do is
;;;    for i = 1 to n
;;;      make the mji for j >= i
;;;      perform the multiplication by Mi as
;;;          mkj = mki mij + mkj    for k#i, j#i
;;;                                  (row i of Mi has the mij we are multiplying by)
;;;                                  (col i of Mi has the mki we are multiplying by)
;;;          mij = mij*mii          for i#j
;;;                                  (row i of Mi which we skipped above)
;;;                                  (mii has the value we are multiplying by)


;;;; What happends in LU Decomposition

;;; pretty much the same except we let mii = 1 instead of 1/aii

;;; Luenberger
;;;   Appendix C, page 340.

;;; Each gaussian elimination step is finding a matrix M_i
;;;   M_i A x = M_i b
;;; let M be M_n * M_n-1 *...* M_1
;;;    M A x = M b
;;;    M is lower triangular
;;;    M A is upper triangular, call it U
;;;    U x = M b
;;;    M**-1 U x = b
;;;    M**-1 is lower triangular


;;; throw in permutation:

;;;  Mn ... M3 M2 M1 A P1 P2 P3 ... Pn = Mn ... M3 M2 M1 A P1 P2 P3 ... Pn 

;;; where the Pi swap columns
;;;   L and M are as before but U includes the Pi
;;;
;;;      U P1 P2 P3 ... Pn = M A P1 P2 P3 ... Pn
;;;    L U                 =   A P1 P2 P3 ... Pn
;;;    L U                 =   A P
;;;    L U P'              =   A
;;;        A  x = b
;;;    L U P' x = b
;;;       solve  LUP'x=b
;;;        by solving  Lz=b     (z=UP'x)
;;;           solving  Uy=z     (y= P'x)
;;;           solving  P'x=y    (ie, x = P y)


;;;; LU Decomposition

;;; *** does not pivot!

;;; It is possible to stuff both L and U into same array
;;; by assuming main diagonal of L is ones.  Just call
;;; the procedure with identical arguments.

(defun matrix-decompose (A L P)
  (let* ((n  (array-dimension A 0)))
    (declare (fixnum n))
    
    (do ((i 0 (1+ i)))
	((>= i n))
      (declare (fixnum i))
      
      ;; *** swap rows as necessary

      ;; "compute" Mi in our heads -- it's the identity plus the (-$ mki) terms
      
      ;; now set A to Mi * A
      ;;   because Mi is an identity in the upper left and lower right
      ;;   do for each row beyond i
      ;;     do for every element
      
      (do ((k   (1+ i) (1+ k))			; Mi is the identity for first i rows
	   (mki 0.0))
	  ((>= k n))
	(declare (fixnum k) (type float mki))
	(setq mki (/ (aref A k i)
		     (aref A i i)))
	
	(do ((j  i (1+ j)))			; calculate all the columns
	    ((>= j n))
	  (declare (fixnum j))
	  (setf (aref A k j)
		(+ (aref A k j)			; main diagonal of Mi has a 1
		   (* (- mki)			; col i has only other nonzero element
		      (aref A i j)))))
	
	(setf (aref L  k i) mki )		; clever L calculation -- see Luenberger
						; SETF after DO so L can be eq A
	))
    L))


;;;; Solving Matrix Problems

;;; Take the system Ax=b
;;; Factor A into LU
;;;   then LUx=b which can be associated as L(Ux)=b
;;;   solve Ly=b and then solve Ux=y

;;; Take the system Ax=b
;;; Factor A into LUP
;;;   then LUPx=b which can be associated as L(U(Px))=b
;;;     solve Ly=b 
;;;     solve Uz=y
;;;     solve Px=z

;;; this allows arbitrary main diagonals for L and U
;;;
(defun matrix-solve-LUP (L U P b)		; solve LUPx=b, put result in b
  (matrix-solve-triangle-lower L b)		; solve  Ly=b, put result in b
  (matrix-solve-triangle-upper U b)		; solve  Ux=y, put result in b
  ;; *** permute!
  b)

;;; *** if we are just solving, don't compute L,
;;;     just solve Ly=b directly for y

(defun matrix-solve (A P b)			; don't need to create L
  (matrix-decompose A A P)			; stuff L and U into A
  (matrix-solve-triangle-lower A b T)		; solve  Ly=b, put result in b
  (matrix-solve-triangle-upper A b NIL)		; solve  Ux=y, put result in b
  ;; *** permute!
  b)


;;;; Matrix-Inversion:  Swap some row and columns

(defun matrix-swap-rows (matrix n i k)
  (declare (fixnum n i k))
  (dotimes (m n)
    (declare (fixnum m))
    (rotatef (aref matrix i m)
	     (aref matrix k m))))

(defun matrix-swap-cols (matrix n j k)
  (declare (fixnum n j k))
  (dotimes (m n)
    (declare (fixnum m))
    (rotatef (aref matrix m j)
	     (aref matrix m k))))

(proclaim '(inline matrix-swap-rows
		   matrix-swap-cols))


;;;; Matrix-Inversion:  Matrix Process

(defun matrix-process (a n k pivot)
  (declare (fixnum n k)
	   (type float pivot))
  
  (dotimes (i n)				; divide column by minus pivot
    (declare (fixnum i))
    (setf (aref a k i)
	  (/ (aref a k i) (- pivot))))
  (setf (aref a k k) pivot)			; got smashed (don't need this statement)
  
  (dotimes (i n)				; reduce matrix
    (declare (fixnum i))
    (if (not (= i k))
	(do ((temp (aref a k i))
	     (j    0 (1+ j)))
	    ((>= j n))
	  (declare (fixnum j) (type float temp))
	  (if (not (= j k))
	      (setf (aref a j i)
		    (+ (* temp (aref a j k))
		       (aref a j i)))))))
  
  (dotimes (j n)				; divide row by pivot
    (declare (fixnum j))
    (setf (aref a j k)
	  (/ (aref a j k) pivot)))
  
  ;; replace pivot by reciprocal
  (setf (aref a k k) (/ 1.0 pivot))
  
  nil)


;;;; Matrix Inversion

;;; let A be the matrix to invert
;;; let R be a permutation matrix that swaps rows    i and k
;;; let C be a permutation matrix that swaps columns j and k
;;; the (C A R) is the matrix with the rows and columns interchanged
;;; I = (CAR) * (CAR)**-1
;;; I = C A R   (CAR)**-1
;;; R**-1 A**-1 C**-1 = (CAR)**-1        premultiply both sides
;;; A**-1 = R (CAR)**-1 C                pre and post multiply
;;;  (but notice that swapping columns and rows roles have reversed)
;;;  two swaps do not change the sign of the determinant

(defun matrix-subr (matrix n k)
  (if (>= k n)
      1.0					; we are done
      (let* ((i k)				; find largest A[i,j]
	     (j k)				;   in lower right corner
	     (max (abs (aref matrix i j))))
	(declare (fixnum i j) (type float max))
	
	(do ((i0 k (1+ i0)))			; find the max
	    ((>= i0 n))
	  (do ((j0 k (1+ j0)))
	      ((>= j0 n))
	    (if (> (abs (aref matrix i0 j0)) max)
		(setq i i0 j j0 max (abs (aref matrix i0 j0))))))
	
	(let ((pivot (aref matrix i j))
	      (d     0.0))
	  (declare (type float pivot d))
	  (if (= pivot 0.0)
	      0.0				; singular matrix!
	      (progn
		(matrix-swap-rows matrix n i k)	; put pivot in right place
		(matrix-swap-cols matrix n j k)	; matrix = C A R
		
		(matrix-process matrix n k pivot)	; invert(C A R)
		;; determinant is the recursive product of pivots
		(setq d (* pivot (matrix-subr matrix n (1+ k))))
		
		(matrix-swap-rows matrix n j k)	; undo permutation
		(matrix-swap-cols matrix n i k)	; matrix = R (CAR)**-1 C
		
		d))
	  ))))

(defun matrix-inverse (matrix)
  (matrix-subr matrix (array-dimension matrix 0) 0))


;;;; Matrix Multiplication

;;; There are several varieties of matrix multiply:
;;;
;;; SCALAR x MATRIX --> MATRIX		MATRIX-MULTIPLY-NUM-MAT
;;; MATRIX x MATRIX --> MATRIX		MATRIX-MULTIPLY-MAT-MAT
;;; ROW    x COLUMN --> SCALAR		MATRIX-MULTIPLY-ROW-COL
;;; COLUMN x ROW    --> MATRIX		MATRIX-MULTIPLY-COL-ROW
;;; ROW    x MATRIX --> ROW		MATRIX-MULTIPLY-ROW-MAT
;;; MATRIX x COLUMN --> COLUMN		MATRIX-MULTIPLY-MAT-COL
;;;
;;; but can't decide if we want a row x col or col x row

(defun matrix-multiply-mat-mat (b c &optional a)
  (let ((l (array-dimension b 0))		; rows of B
	(m (array-dimension c 0))		; rows of C
	(n (array-dimension c 1)))		; cols of C
    (declare (fixnum l m n))
    
    (or a (setq a (make-array (list l n)
			      :element-type 'float :initial-element 0.0)))
    
    (if (or (not (= l (array-dimension a 0)))
	    (not (= n (array-dimension a 1)))
	    (not (= m (array-dimension b 1))))
	(error "MATRIX-MULTIPLY given bad dimensions"))
    
    (do ((i 0 (1+ i)))
	((>= i l) a)
      (declare (fixnum i))
      (do ((j 0 (1+ j))
	   (sum 0.0 0.0))
	  ((>= j n))
	(declare (fixnum j)
		 (float sum))
	(do ((k 0 (1+ k)))
	    ((>= k m))
	  (declare (fixnum k))
	  (setq sum (+ sum (* (aref b i k)
			      (aref c k j)))))
	(setf (aref a i j) sum)))
    ))


;;;; Multiply ROW x MAT and MAT x COL

;;; ROW * MAT
;;; 
(defun matrix-multiply-row-mat (b c &optional a)
  (let ((n (array-dimension c 1))
	(m (array-dimension b 0)))
    (declare (fixnum m n))
    (or a (setq a (make-array n
			      :element-type 'float :initial-element 0.0)))
    (do ((i 0 (1+ i)))
	((>= i n) a)
      (declare (type fixnum i))
      (do ((j 0 (1+ j))
	   (sum 0.0))
	  ((>= j m)
	   (setf (aref a i) sum)
	   nil)
	(declare (fixnum j)
		 (float sum))
	(setq sum (+ sum (* (aref b j)
			    (aref c j i))))))))

;;; MAT * COL
;;;
(defun matrix-multiply-mat-col (b c &optional a)
  (let ((n (array-dimension b 0))
	(m (array-dimension c 0)))
    (declare (fixnum m n))

    (or a (setq a (make-array n
			      :element-type 'float :initial-element 0.0)))
    (do ((i 0 (1+ i)))
	((>= i n) a)
      (declare (fixnum i))
      (do ((j 0 (1+ j))
	   (sum 0.0))
	  ((>= j m)
	   (setf (aref a i) sum)
	   nil)
	(declare (fixnum j)
		 (float sum))
	(setq sum (+ sum (* (aref b i j)
			    (aref c j))))))
    ))


;;;; Multiply Matrix by a Scalar

;;; b is a scalar, c is a scalar, vector, or matrix

(defun matrix-multiply-num-mat (b c &optional a)
  (cond ((numberp c) (* b c))
	((= (array-rank c) 1)
	 ;; times a vector
	 (or a (setq a (make-array (array-dimension c 0)
				   :element-type 'float :initial-element 0.0)))
	 (do ((i 0 (1+ i))
	      (n (array-dimension c 0)))
	     ((>= i n) a)
	   (declare (fixnum i n))
	   (setf (aref a i)
		 (* b (aref c i)))))
	((= (array-rank c) 2)
	 ;; times a matrix
	 (or a (setq a (make-array (list (array-dimension c 0) (array-dimension c 1))
				   :element-type 'float :initial-element 0.0)))
	 (do ((i 0 (1+ i))
	      (n (array-dimension c 0)))
	     ((>= i n) a)
	   (declare (fixnum i n))
	   (do ((j 0 (1+ j))
		(m (array-dimension c 0)))
	       ((>= j m))
	     (declare (fixnum j m))
	     (setf (aref a i j)
		   (* b (aref c i j))))))
	(t (error "MATMUL-SCALAR not given vector or matrix"))))
  

;;;; Matrix Multiply

;;; There are several varieties of matrix multiply:
;;;
;;; scalar x matrix --> matrix		MATRIX-MULTIPLY-NUM-MAT
;;; matrix x matrix --> matrix		MATRIX-MULTIPLY-MAT-MAT
;;; row    x column --> scalar		MATRIX-MULTIPLY-ROW-COL
;;; column x row    --> matrix		MATRIX-MULTIPLY-COL-ROW
;;; row    x matrix --> row		MATRIX-MULTIPLY-ROW-MAT
;;; matrix x column --> column		MATRIX-MULTIPLY-MAT-COL

(defun matrix-multiply (b c &optional a)
  (cond ((numberp b) (matrix-multiply-num-mat b c a))
	((numberp c) (matrix-multiply-num-mat c b a))
	((and (= 2 (array-rank b))
	      (= 2 (array-rank c))
	      (= (array-dimension b 1)
		 (array-dimension c 0)))
	 (matrix-multiply-mat-mat b c a))
	((and (= 1 (array-rank b))
	      (= 2 (array-rank c))
	      (= (array-dimension b 0)
		 (array-dimension c 0)))
	 (matrix-multiply-row-mat b c a))
	((and (= 2 (array-rank b))
	      (= 1 (array-rank c))
	      (= (array-dimension b 1)
		 (array-dimension c 0)))
	 (matrix-multiply-mat-col b c a))
	(t (error "MATMUL blew it"))
	))


;;;; Matrix Inversion Tests
#|
(eval-when (eval)
|#
#|
  (defun matrix-random-vector (n)
    (declare (fixnum n))
    (let ((v (make-array n :element-type 'float :initial-element 0.0)))
      (declare (type (array float (*)) v))
      (do ((i 0 (1+ i)))
	  ((>= i n) v)
	(declare (fixnum i))
	(setf (aref v i) (random 1.0)))
      v))
  
  (defun matrix-random-matrix (n)
    (declare (fixnum n))
    (let ((A (make-array (list n n) :element-type 'float :initial-element 0.0)))
      (declare (type (array float (* *)) A))
      (do ((i 0 (1+ i)))
	  ((>= i n))
	(declare (fixnum i))
	(do ((j 0 (1+ j)))
	    ((>= j n))
	  (declare (fixnum j))
	  (setf (aref A i j) (random 1.0))))
      A))
  
  (defun matrix-solve-lower (n)
    (declare (fixnum n))
    (let ((L (make-array (list n n) :element-type 'float :initial-element 0.0))
	  (x (make-array n          :element-type 'float :initial-element 0.0))
	  (b (matrix-random-vector n)))
      (declare (type (array float (* *)) L)
	       (type (array float   (*)) x b))
      (do ((i 0 (1+ i)))
	  ((>= i n))
	(declare (fixnum i))
	(do ((j 0 (1+ j)))
	    ((> j i))
	  (declare (fixnum j))
	  (setf (aref L i j) (random 1.0))))
      (matrix-print L)
      (matrix-print b)
      (replace x b)
      (matrix-solve-triangle-lower L x)
      (matrix-print x)
      (matrix-print (matrix-multiply L x))))
  
  (defun matrix-solve-upper (n)
    (declare (fixnum n))
    (let ((U (make-array (list n n) :element-type 'float :initial-element 0.0))
	  (x (make-array n          :element-type 'float :initial-element 0.0))
	  (b (matrix-random-vector n)))
      (declare (type (array float (* *)) U)
	       (type (array float   (*)) x b))
      (do ((i 0 (1+ i)))
	  ((>= i n))
	(declare (fixnum i))
	(do ((j i (1+ j)))
	    ((>= j n))
	  (declare (fixnum j))
	  (setf (aref U i j) (random 1.0))))
      (matrix-print U)
      (matrix-print b)
      (replace x b)
      (matrix-solve-triangle-upper U x)
      (matrix-print x)
      (matrix-print (matrix-multiply U x))))
  
  (defun matrix-decompose-random-matrix (n)
    (declare (fixnum n))
    (let* ((A (matrix-random-matrix n)))
      (declare (type (array float (* *)) A))
      (matrix-print A)
      (let ((L (matrix-decompose A)))
	(declare (type (array float (* *)) L))
	(matrix-print L)
	(matrix-print A)
	(matrix-print (matrix-multiply-mat-mat L A)))))
  
  (defun matrix-solve-test (n)
    (declare (fixnum n))
    (let ((A (matrix-random-matrix n))
	  (b (matrix-random-vector n)))
      (declare (type (array float (* *)) A)
	       (type (array float   (*)) b))
      (matrix-print A)
      (matrix-print b)
      (let ((A-copy (matrix-copy A))
	    (x      (matrix-copy b)))
	(declare (type (array float (* *)) A-copy)
		 (type (array float   (*)) x))
	(matrix-solve A-copy x)
	(matrix-print x)
	(matrix-print (matrix-multiply A x)))))
  
  (defun matrix-invert-random-matrix (n)
    (declare (fixnum n))
    (let* ((A (matrix-random-matrix n))
	   (B (matrix-copy A)))
      (declare (type (array float (* *)) A B))
      (matrix-print A)
      (print (list 'det (matrix-inverse B)))
      (matrix-print B)
      (matrix-print (matrix-multiply-mat-mat A B))))
  
  )
|#
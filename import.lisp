;;; -*- Syntax: Common-lisp; Package: USER; Mode: LISP -*-

;;;; IMPORT FILE

;;; IMPORT-FILE reads definitions from a file at compile time
;;;   and loads the file at run time

;;;  (c) Copyright Gerald Roylance 1983, 1984, 1985, 1986, 1987
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   blow away silly conditional assembly
;;;     fundamental problems
;;;       host pathname syntax is different on different machines
;;;         dec host syntax is HOST::
;;;       package ZL-USER only exists on LISPM
;;;   use a pathname mechanism
;;;   check file dates and load in the changes since last compile
;;;   (import-file filename &optional host)
;;;     would eliminate need for stripping the host

(in-package "CLMATH")

(provide "IMPORT")

(export '(IMPORT-FILE))

(eval-when (compile load eval)
  #+common (require "ATTRIBUTES")
  #+lucid  (require "ATTRIBUTES" "/homes/glr/lisp/io/attrib")
  NIL)

(defvar imported-file-list nil)


;;;; Manifest Host Parsing

;;; Process a string with a LISPM-style manifest host

(defun import-filename (filename)
  (if (pathnamep filename)
      filename
      (let ((j (position #\: filename)))

	(if (null j)
	    (error "IMPORT-FILENAME:  no explicit host in ~s" filename))

	(let ((host (subseq filename 0 j)))
	  (cond ((string-equal host "OZ")	; convert it
		 (warn "converting OZ pathname to TX")
		 (let* ((fn   (copy-seq filename))
			(dir0 (position #\< fn))	; position of <
			(dir1 (position #\> fn)))	; position of >
		   (nstring-downcase fn)
		   (nsubstitute #\/ #\. fn :start dir0 :end dir1)
		   (setf (aref fn dir0) #\/)
		   (setf (aref fn dir1) #\/)
		   (parse-namestring (concatenate 'string "/homes" (subseq fn dir0))
				     "TRIX")))
		(t
		 (parse-namestring (subseq filename (1+ j)) host)))
	  ))))


;;;; Package Hack

;;; Determine a file's package

;;; "-*- Package: USER -*-" means different things according to syntax

;;; *** should also check for (in-package ...) (Lucid's convention)

(defun import-package (filename)
  (let* ((attr-plist   (fs:pathname-attribute-list filename))
	 (package-spec (getf attr-plist :package)))

    (cond ((null package-spec) *package*)	; no package? use the current one

	  ((listp package-spec)			; create form...
	   (let ((pack (find-package (car package-spec))))
	     (if (null pack)			; time to make the package
		 (apply #'make-package package-spec)
		 pack)))

	  (t
	   (let ((pack (find-package package-spec)))
	     (if (null pack)
		 (progn (cerror "Make the package"
				"Package ~a does not exist."
				package-spec)
			(make-package package-spec))
		 pack))))))


;;;; Import Load

(defmacro import-load (fxname)
  (let* ((lsp   fxname)				; get the effective pathname
	 (ext   (progn
		  #+(and LUCID (not HP)) "lbin"	; normal lucid
		  #+(and LUCID      HP ) "b"	; hp version of lucid
		  #+COMMON               "LAP"
		  #+3600                 "BIN"))
	 (bin   (make-pathname :type ext :defaults lsp)))

    `(let ((tag (probe-file ',lsp)))
       (declare (special imported-file-list))
       (or (boundp 'imported-file-list)
	   (setq    imported-file-list nil))
       (cond ((member tag imported-file-list :test #'equal))
	     (tag
	      (load (or (probe-file ',bin) tag))
	      (push tag imported-file-list)
	      nil)
	     (t
	      (error "IMPORT-LOAD:  File Not Found ~S" ',lsp))))))


;;;; Import File

(defmacro import-file (fn)
  
  (let* ((filename  (import-filename fn))	; get the effective pathname
	 (*package* (import-package filename))	; find what package it wants
	 (form      nil))

    (with-open-file (filei filename)
      
      ;; find a sensible form to import
      (dotimes (i 3)				; look at the first 3 forms
	(setq form (read filei nil :EOF))
	(cond ((eq form :EOF)
	       (error "IMPORT-FILE:  EOF in ~s" filename))
	      ((or (atom form)
		   (not (atom (car form))))
	       (error "IMPORT-FILE:  bad format in ~s" filename))
	      ((eq (car form) 'in-package))	; go around again
	      ((eq (car form) 'eval-when)
	       (return nil))
	      )
	))

    `(progn 'compile
	    (IMPORT-LOAD ,filename)
	    ,form)
    ))

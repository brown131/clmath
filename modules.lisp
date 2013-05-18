;;; -*- Mode: LISP; Syntax: COMMON-LISP; Package: USER -*-

;;;; Module Name to Pathname Mappings

;;;  (c) Copyright Gerald Roylance 1989
;;;      All Rights Reserved.
;;;  No part of this material may be used, photocopied,
;;;  reproduced or committed to magnetic, electronic or any other
;;;  form of media without prior written consent.

;;; Bugs and Fixes
;;;   compile only those files that need it
;;;   make distribution tapes (maybe binary only)
;;;     include/exclude some systems
;;;   separate compilation without loading
;;;   bug with declares
;;;   release/version/edit
;;;   patches
;;;   might have an argument that imports the packages
;;;     external symbols.

(in-package "CLMATH")

;;;; Database Manipulation Routines

;;; Keep the correspondence between names and pathnames
;;;
(defvar module-database (make-hash-table :test #'equal :size 200))

(defun module-pathname-get (name)
  (let ((pathnames (gethash name module-database)))
    (if (null pathnames)
	(warn "Module ~A has no pathname translations" name))
    pathnames))

(defun module-pathname-put (name pathname)
  (if (null pathname)
      (remhash name module-database)
      (setf (gethash name module-database) pathname)))

(defun module-require (name)
  (require name (module-pathname-get name)))

(defun module-provide (name)
  (provide name))

(defmacro module-namestring (module-name host namestring)
  `(module-pathname-put ,module-name (parse-namestring ,namestring ,host)))


;;;; The Database

;;; should acquire this data from a few canonical places

(module-namestring "ANAGRAM"             "TX" "/homes/glr/funct/number/anagram")

(module-namestring "BFILE"               "TX" "/homes/glr/sort/btree/bfile")
(module-namestring "BTEST"               "TX" "/homes/glr/sort/btree/btest")
(module-namestring "BTREE"               "TX" "/homes/glr/sort/btree/btree")

(module-namestring "CLX"                 "TX" "/homes/glr/pstrm/clx")
(module-namestring "COMBINATORICS"       "TX" "/homes/glr/funct/number/combin")

(module-namestring "FILTER"              "TX" "/homes/glr/design/filter/filter")
(module-namestring "FILTER-AF100"        "TX" "/homes/glr/design/filter/af100")
(module-namestring "FILTER-BIQUAD"       "TX" "/homes/glr/design/filter/biquad")
(module-namestring "FILTER-BIQUAD1"      "TX" "/homes/glr/design/filter/biquad1")
(module-namestring "FILTER-CATALOG"      "TX" "/homes/glr/design/filter/catalog")
(module-namestring "FILTER-CIRCUITS"     "TX" "/homes/glr/design/filter/filckt")
(module-namestring "FILTER-COMPONENTS"   "TX" "/homes/glr/design/filter/filcmp")
(module-namestring "FILTER-DEFS"         "TX" "/homes/glr/design/filter/fdefs")
(module-namestring "FILTER-PLOT"         "TX" "/homes/glr/design/filter/fplot")
(module-namestring "FILTER-SALLEN"       "TX" "/homes/glr/design/filter/sallen")
(module-namestring "FILTER-TRANSFORM"    "TX" "/homes/glr/design/filter/fxform")
(module-namestring "FILTER-WIZARD"       "TX" "/homes/glr/papers/phd/wizard/wizfil")
(module-namestring "FILTER-ZEROFINDER"   "TX" "/homes/glr/design/filter/zfind")

(module-namestring "IGES"                "TX" "/homes/glr/design/iges/iges")
(module-namestring "IGES-CONVERT"        "TX" "/homes/glr/design/iges/convert")
(module-namestring "IGES-PARAMETERS"     "TX" "/homes/glr/design/iges/param")

(module-namestring "MATRIX"              "TX" "/homes/glr/funct/matrix")

(module-namestring "OPTIMIZE-GOLDEN"     "TX" "/homes/glr/funct/optim/golden")
(module-namestring "OPTIMIZE-MARQUARDT"  "TX" "/homes/glr/funct/optim/marq")
(module-namestring "OPTIMIZE-REMES"      "TX" "/homes/glr/funct/optim/remes")

(module-namestring "PERMUTATIONS"        "TX" "/homes/glr/funct/number/permute")
(module-namestring "PLOT"                "TX" "/homes/glr/pstrm/plot/plot")
(module-namestring "POSTSCRIPT-AFM"      "TX" "/homes/glr/pstrm/postscript/afm")
(module-namestring "PSTRM"               "TX" "/homes/glr/pstrm/pstrm")
(module-namestring "PSTRM-DEFS"          "TX" "/homes/glr/pstrm/psdefs")
(module-namestring "PSTRM-POSTSCRIPT"    "TX" "/homes/glr/pstrm/postscript/pspost")
(module-namestring "PSTRM-PTTY"          "TX" "/homes/glr/pstrm/ptty")
(module-namestring "PSTRM-SUBS"          "TX" "/homes/glr/pstrm/pssubs")
(module-namestring "PSTRM-SUPDUP"        "TX" "/homes/glr/pstrm/supgrf")
(module-namestring "PSTRM-TEKTRONIX"     "TX" "/homes/glr/pstrm/tk4014")
(module-namestring "PSTRM-TV"            "TX" "/homes/glr/pstrm/tv")

(module-namestring "RECORD"              "TX" "/homes/glr/lisp/record")

(module-namestring "SPECTRUM"            "TX" "/homes/glr/design/signal/spectrum")
(module-namestring "SPRINGS"             "TX" "/homes/glr/design/mech/springs.lisp")
(module-namestring "STATISTICS"          "TX" "/homes/glr/funct/statis/statis")
(module-namestring "STATISTICS-BINOMIAL" "TX" "/homes/glr/funct/statis/binomial")
(module-namestring "SYMTAB"              "TX" "/homes/glr/lisp/comp/symtab")

(module-namestring "TERMCAP"             "TX" "/homes/glr/lisp/io/termcap")

(module-namestring "VMEM"                "TX" "/homes/glr/sort/btree/vmem")

(module-namestring "ZEROS-BAIRSTOW"      "TX" "/homes/glr/funct/zeros/bairst")
(module-namestring "ZEROS-BISECTION"     "TX" "/homes/glr/funct/zeros/bisection")
(module-namestring "ZEROS-BRENT"         "TX" "/homes/glr/funct/zeros/brent")
(module-namestring "ZEROS-NEWTON"        "TX" "/homes/glr/funct/zeros/newton")

(module-namestring "INTERVAL-ARITHMETIC" "TX" "/homes/glr/funct/interval")
(module-namestring "POLYNOMIALS"         "TX" "/homes/glr/funct/polyn")
(module-namestring "HORNER"              "TX" "/homes/glr/funct/horner")
(module-namestring "ERROR-FUNCTION"      "TX" "/homes/glr/funct/functions/erf")
(module-namestring "GAMMA-FUNCTION"      "TX" "/homes/glr/funct/functions/gamma")
(module-namestring "BETA-FUNCTION"       "TX" "/homes/glr/funct/functions/beta")
(module-namestring "BESSEL-FUNCTIONS"    "TX" "/homes/glr/funct/functions/bessel")
(module-namestring "ELLIPTIC-FUNCTIONS"  "TX" "/homes/glr/funct/functions/ellip")

;;; not yet fixed...

(module-namestring "TAPE"                "TX" "/homes/glr/lisp/io/tape")
(module-namestring "ANSI"                "TX" "/homes/glr/lisp/io/ansi")
(module-namestring "NICAD"               "TX" "/homes/glr/design/nicad/src/nicad")
(module-namestring "FIFO"                "TX" "/homes/glr/sort/fifo")

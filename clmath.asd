(defpackage #:clmath-asd
  (:use :cl :asdf))

(in-package :clmath-asd)

(defsystem clmath
  :name "clmath"
  :author "Gerald Roylance"
  :serial t
  :components ((:file "defpackage")
               (:file "consts" :depends-on ("defpackage"))
               (:file "dft" :depends-on ("defpackage"))
               (:file "horner" :depends-on ("defpackage"))
               (:file "matrix" :depends-on ("defpackage"))
               (:module "functions" :depends-on ("defpackage")
                :components (;(:file "basic")
                             (:file "bessel")
                             (:file "beta")
                             (:file "ellip")
                             (:file "erf")
                             (:file "extfcn") ;basic
                             (:file "gamma")))
               (:module "integrate" :depends-on ("defpackage")
                :components ((:file "integr")
                             (:file "runge")))
               (:module "number" :depends-on ("defpackage")
                :components ((:file "combin")
                             (:file "factor")
                             (:file "fib")
                             (:file "mod")))
               (:module "optim" :depends-on ("defpackage")
                :components ((:file "fit")
                             ;(:file "fmfp")
                             (:file "marq")
                             (:file "regres")))
               (:module "statis" :depends-on ("defpackage")
                :components ((:file "binomial" :depends-on ("statis"))
                             (:file "poisson")
                             (:file "statis")))
               (:module "zeros" :depends-on ("defpackage")
                :components ((:file "bisection")
                             (:file "falsep")))   
               )
)

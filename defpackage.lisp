(in-package :common-lisp-user)

(defpackage :clmath
  (:use :cl)
  (:export 
    ;; consts
    #:velocity-of-light
    #:electronic-charge
    #:Plancks-constant
    #:Avagadros-number
    #:atomic-mass-unit
    #:electron-rest-mass
    #:proton-rest-mass
    #:magnetic-flux-quantum
    #:Boltzmanns-constant
    #:gravitational-constant
    #:permeability-free-space
    #:permittivity-free-space
    ;; dft
    #:dft-610
    #:dft-618
    #:dft-reverse-bits-16
    #:dft-reverse-array
    #:dft-forward
    ;; horner
    #:poly-expand
    #:horner-polynomial
    ;; matrix
    #:matrix-print-vector
    #:matrix-print-matrix
    #:matrix-print
    #:matrix-identity
    #:matrix-copy
    #:matrix-diagonal
    #:matrix-add 
    #:matrix-sub
    #:matrix-solve-triangle-lower
    #:matrix-solve-triangle-upper
    #:matrix-decompose
    #:matrix-solve-LUP
    #:matrix-solve
    #:matrix-swap-rows
    #:matrix-swap-cols
    #:matrix-process
    #:matrix-inverse
    #:matrix-multiply-mat-mat
    #:matrix-multiply-row-mat
    #:matrix-multiply-mat-col
    #:matrix-multiply-num-mat
    #:matrix-multiply
    ;; FUNCTIONS:
    ;; basic
    ;#:fraction 
    ;#:entier 
    ;#:square-root 
    ;; bessel
    #:bessel-j
    #:bessel-i
    ;; beta 
    #:beta-function-naive
    #:log-beta-function 
    #:beta-function
    #:beta-warning
    #:betacf
    #:betai
    ;; ellip 
    #:elliptic-integral-k
    #:elliptic-integral-kC 
    #:elliptic-integral-e
    ;; erf
    #:erfc
    #:erf
    #:error-function
    ;; extfcn
    #:exponential
    #:logarithm
    #:sine
    #:cosine
    #:arctan
    ;; gamma
    #:gamma-function-nbs
    #:gamma-stirling
    #:gammln
    #:fgammln
    #:gamma-function 
    #:log-gamma-function 
    #:gamma-function-reciprocal
    ;; INTEGRATE:
    ;; integr
    #:integrate-trapezoidal
    #:integrate-simpson
    #:integrate-simpson-recursive
    ;; runge
    #:euler
    #:integrate-euler
    #:runge-kutta
    #:integrate-runge-kutta
    #:correctors
    #:predictors
    #:multi-step-start
    #:multi-step-interpolate
    #:multi-step-predict
    #:multi-step-correct
    #:integrate-multi-step-h
    ;; NUMBER:
    ;; combin
    #:factorial
    #:dlog-factorial
    #:dfactorial
    #:ffactorial
    #:gprod
    #:gfact
    #:choose-naive
    #:choose
    #:log-choose
    #:dlog-choose
    #:fchoose
    #:dchoose
    #:catalan-number
    #:bellr
    #:bell-number
    ;; factor
    #:factor
    #:totient
    ;; fib
    #:fib-slow
    #:fib
    #:fib-binet
    ;; mod
    #:gcd-extended
    #:modinv
    #:chinese
    #:modpower
    #:random-big
    #:jacobi-symbol
    #:prime-test-fermat
    #:prime-test
    #:random-prime
    #:smaller-prime
    ;; OPTIM:
    ;; fit
    #:fit-fcn
    #:fit
    ;; fmfp
    ;#:fmfp
    ;; marq
    #:marquardt
    ;; regres
    #:regres
    ;; STATIS:
    ;; binomial
    #:binomial-density-naive
    #:binomial-density
    #:binomial-cumulative-naive
    #:binomial-cumulative
    #:binomial-tail
    #:binomial-random-number
    ;; poisson
    #:poisson-density
    #:poisson-cumulative
    #:poisson-random-number
    ;; statis
    #:uniform-density
    #:uniform-cumulative
    #:uniform-random-number
    #:exponential-density
    #:exponential-cumulative
    #:exponential-random-number
    #:cauchy-density
    #:cauchy-cumulative
    #:cauchy-random-number
    #:normal-density
    #:normal-cumulative
    #:normal-tail
    #:normal-random-number
    #:gamma-density
    #:gamma-cumulative
    #:gamma-random-number
    #:beta-density
    #:beta-cumulative
    #:beta-random-number
    #:chi-square-density
    #:chi-square-cumulative
    #:chi-square-random-number
    #:t-density
    #:t-two-sided
    #:t-cumulative
    #:t-random-number
    #:f-density
    #:f-cumulative
    #:f-random-number
    ;; ZEROS:
    ;; bisection
    #:bisection-search
    ;; falsep
    #:false-position-search
    #:false-position-converge
    #:converge
        ))

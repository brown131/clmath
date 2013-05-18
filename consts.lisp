;;; -*- Mode: LISP; Syntax: Common-lisp; Package: USER -*-

;;;; Mathematical and Physical Constants

;;;  (c) Copyright Gerald Roylance 1980, 1985
;;;      All Rights Reserved.
;;;  This file may be distributed noncommercially provided
;;;  that this notice is not removed.

;;; Bugs and Fixes
;;;   Better treatment of physical constants
;;;     in TX:/tx/glr/funct/sym/unit.lisp

(in-package "CLMATH")

;;;; Physical Units

;;; base units of Si:
;;;      meter, kilogram, second, ampere, Kelvin, Candela

;;; volt
;;; ohm
;;; coulomb
;;; farad
;;; henry   (inductance that produces 1V emf when di/dt=1 A per sec)
;;; Weber   (mag. flux producing 1V emf in one turn as the flux
;;;          goes uniformly to zero in 1 second)
;;; Tesla   (flux density of 1 Weber/square meter)

;;;				MKS		EMU
;;; magnetizing force	H	A-turns/meter	=Oersteds * 1e3 / (4 pi)
;;; flux		phi	Weber		=Maxwells * 1e-8
;;; flux density	B	Weber/sq meter	=Gauss * 1e-4
;;; permeability	mu	H/m		=Gauss/Oersted *4e-7 * pi

(defmacro defcon (name symbol number . dimensions)
  `(defconstant ,name ,number))


;;;; Physical Constants

(defcon velocity-of-light	c	2.9979250e8 meters per second)
(defcon electronic-charge	q	1.6021917e-19 coulomb)
(defcon Plancks-constant	h	6.626196e-34 Joule second)
(defcon Avagadros-number	N	6.022169e26 per K mole)
(defcon atomic-mass-unit	amu	1.660531e-27 K gram)
(defcon electron-rest-mass	me	9.109558e-31 K gram)
(defcon proton-rest-mass	mp	1.672614e-27 K gram)
(defcon magnetic-flux-quantum	phi	2.067854e-15 Tesla meter meter)
(defcon Boltzmanns-constant	k	1.380622e-23 Joule per Kelvin)
(defcon gravitational-constant	gamma	6.6732e-11 Newton (meter 2) per (K gram 2))

(defcon permeability-free-space mu0     (* 4.0 pi 1.0e-7) henry per meter)
(defcon permittivity-free-space epsilon0 8.854e-12 farads per meter)

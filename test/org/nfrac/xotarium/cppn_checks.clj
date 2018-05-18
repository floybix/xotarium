(ns org.nfrac.xotarium.cppn-checks
  (:require [org.nfrac.xotarium.cppn :as cppn]
            [org.nfrac.xotarium.util :as util]
            [org.nfrac.xotarium.creature :as cre]
            [clojure.spec :as s]
            [clojure.spec.gen :as gen]
            [clojure.test.check.random :as random]
            [clojure.spec.test :as stest]
            [clojure.test.check.clojure-test :as ctcc]
            [clojure.test :as t
             :refer (is deftest testing run-tests)]))

(alter-var-root #'ctcc/*report-shrinking* (constantly true))
(alter-var-root #'ctcc/*report-trials* (constantly ctcc/trial-report-periodic))

(def instr-syms
  (concat
   (stest/enumerate-namespace 'org.nfrac.xotarium.cppn)
   ))

(alias 'stc 'clojure.spec.test.check)
(def opts {::stc/opts {:num-tests 1000}})

(def mut-ops
  {:node #'cppn/mutate-add-node
   :conn #'cppn/mutate-add-conn
   :rewire #'cppn/mutate-rewire-conn})

(defn apply-mutations
  [cppn rng]
  (cppn/mutate-add-conn cppn rng)
  ;[cppn ops rng]
  #_(loop [cppn cppn
         ops ops
         rng rng]
    (if-let [op (first ops)]
      (let [f (mut-ops op)
            [rng rng*] (random/split rng)]
        (recur (f cppn rng*) (rest ops) rng))
      cppn)))

(def test-cppn
  {:inputs #{:y :d :x :bias}, :outputs #{:phase-off :factor-c :bone :muscle :angle :factor-b :factor-a}, :finals #{:bone :muscle}, :zerod #{}, :nodes {:i0 :gaussian, :1b059f9a :linear}, :edges {:phase-off {:factor-c 1.823268201409166}, :factor-c {:d 0.9101854099431943}, :i0 {:d -0.37633016715078604}, :bone {:i0 0.8346684577519059, :bias -0.6115188419769843}, :muscle {:factor-c 0.4725042414527339, :factor-a 1.0759773296273496, :1b059f9a 1.495849810377079}, :angle {:factor-a 0.9298196858048878}, :factor-b {:y 0.42603469171302544}, :factor-a {:x 1.3303794805832996, :factor-c -0.6193703472769874}, :1b059f9a {:factor-c 0.6121740892494435}}})



(s/fdef apply-mutations
        :args (s/cat :cppn (-> ::cppn/cppn
                               (s/with-gen #(gen/return test-cppn)))
                     ;:ops (s/coll-of (set (keys mut-ops)))
                     :rng ::util/rng)
        :ret ::cppn/cppn)

(deftest cppn-fns-test
  (stest/instrument instr-syms)
  (-> `[apply-mutations
        ;cppn/randomise-weights
        ]
      (stest/check opts)
      (stest/summarize-results))
  (stest/unstrument))

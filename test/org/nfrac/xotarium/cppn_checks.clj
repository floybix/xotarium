(ns org.nfrac.xotarium.cppn-checks
  (:require [org.nfrac.xotarium.cppn :as cppn]
            [org.nfrac.xotarium.util :as util]
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
(def opts {::stc/opts {:num-tests 500}})

(def mut-ops
  {:node #'cppn/mutate-add-node
   :conn #'cppn/mutate-add-conn
   :rewire #'cppn/mutate-rewire-conn})

(defn apply-mutations
  [cppn ops rng]
  (loop [cppn cppn
         ops ops
         rng rng]
    (if-let [op (first ops)]
      (let [f (mut-ops op)
            [rng rng*] (random/split rng)]
        (recur (f cppn rng*) (rest ops) rng))
      cppn)))

(s/fdef apply-mutations
        :args (s/cat :cppn ::cppn/cppn
                     :ops (s/coll-of (set (keys mut-ops)))
                     :rng ::util/rng)
        :ret ::cppn/cppn)

(deftest cppn-fns-test
  (stest/instrument instr-syms)
  (-> `[apply-mutations
        cppn/randomise-weights
        ]
      (stest/check opts)
      (stest/summarize-results))
  (stest/unstrument))

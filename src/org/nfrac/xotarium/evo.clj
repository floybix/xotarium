(ns org.nfrac.xotarium.evo
  (:require [org.nfrac.xotarium.grn-creature :as grncre]
            [org.nfrac.xotarium.creature :as cre]
            [org.nfrac.xotarium.cppn :as cppn]
            [org.nfrac.xotarium.cppn-compile :as cc]
            [org.nfrac.xotarium.grn.greans :as grn]
            [org.nfrac.xotarium.plant-cave :as cave]
            [org.nfrac.liquidfun.testbed :as bed]
            [org.nfrac.liquidfun.core :as lf :refer [body! joint!
                                                     particle-system!]]
            [org.nfrac.liquidfun.vec2d :as v2d]
            [org.nfrac.xotarium.util :as util]
            [org.nfrac.xotarium.util.algo-graph :as graph]
            [quil.core :as quil :include-macros true]
            [quil.middleware]
            [clojure.pprint]
            [clojure.java.io :as io]
            [clojure.spec :as s]
            [clojure.test.check.random :as random])
  (:import (org.bytedeco.javacpp
            liquidfun$b2ContactListener
            liquidfun$b2ParticleBodyContact
            liquidfun$b2ParticleContact
            liquidfun$b2World
            liquidfun$b2Body
            liquidfun$b2Transform
            liquidfun$b2Vec2
            liquidfun$b2ParticleTriad
            liquidfun$b2ParticleHandle
            liquidfun$b2ParticleSystem
            liquidfun$b2ParticleGroup
            liquidfun$b2ParticleDef
            liquidfun$b2QueryCallback)))

(def parameter-defaults
  {:crossover-prob 0.5
   :population-size 4
   :generations 3})

(def sim-steps (* 60 5))

(def xy-beh-resolution 0.1);; each multiple of this is a distinct behaviour.

(s/def ::parameters (s/keys))

;;; Behavioural measures:
;;; Each is scaled such that a difference of 1.0 is a distinct behaviour.
;;; Distinct (grouped) phenotypes / behaviour types are these, int trucated.

(s/def ::behaviour-type
  (s/map-of any? int?))

(s/def ::count nat-int?)

(s/def ::representative ::grncre/genome)

(s/def ::behaviour-archive
  (s/every-kv ::behaviour-type
              (s/keys :req-un [::count
                               ::representative])))

(defn mean-point
  [xys]
  (let [summed (reduce v2d/v-add xys)]
    (v2d/v-scale summed (/ 1.0 (count xys)))))

(defn point-in-time-measures
  [state]
  (let [creature (:creature state)
        grps (apply concat (vals (:groups creature)))
        cents (map lf/center grps)
        [x y] (mean-point cents)]
    {:cent-x (int (/ x xy-beh-resolution))
     :cent-y (int (/ y xy-beh-resolution))}))

(defn behaviour-measures
  [genome]
  (let [state-ht (->> (grncre/setup genome)
                      (iterate grncre/step)
                      (drop (quot sim-steps 2))
                      (first))
        beh-ht (point-in-time-measures state-ht)
        ;; TODO don't hold on to your head
        state-ft (->> state-ht
                      (iterate grncre/step)
                      (drop (quot sim-steps 2))
                      (first))
        beh-ft (point-in-time-measures state-ft)]
    ;; combine the maps from each point in time, set unique keys
    (reduce (fn [m [time-k beh]]
              (let [ks (map vector (repeat time-k) (keys beh))]
                (merge m (zipmap ks (vals beh)))))
            {}
            [[:ht beh-ht]
             [:ft beh-ft]])))

(s/def ::behaviour
  (s/keys))

(s/def ::beh-freq nat-int?)

(s/def ::indiv
  (s/keys :req-un [::grncre/genome]))

(s/def ::indiv-evald
  (s/keys :req-un [::grncre/genome
                   ::behaviour
                   ::beh-freq]))

(defn eval-popn
  [popn behaviour-archive]
  ;; TODO: parallelize?
  (let [bpopn (map (fn [indiv]
                     (let [beh (behaviour-measures (:genome indiv))]
                       (println "    beh=" beh)
                       (assoc indiv :behaviour beh)))
                   popn)
        ;; iteratively add behaviours to archive; so first in popn is more novel
        [epopn ba] (loop [indivs bpopn
                          i 0
                          beh-archive behaviour-archive
                          epopn ()]
                     (if-let [indiv (first indivs)]
                       (let [beh (:behaviour indiv)
                             beh-freq (get-in beh-archive [beh :count] 0)
                             nov? (zero? beh-freq)]
                         (recur (rest indivs)
                                (inc i)
                                (if nov?
                                  (assoc beh-archive beh
                                         {:count 1
                                          :representative (:genome indiv)})
                                  (update-in beh-archive [beh :count] inc))
                                (conj epopn (assoc indiv :beh-freq beh-freq))))
                       ;; done
                       [epopn beh-archive]))
        ]
    [epopn ba]))

(s/fdef eval-popn
        :args (s/cat :popn (s/every ::indiv :min-count 2)
                     :ba ::behaviour-archive)
        :ret (s/cat :epopn (s/every ::indiv-evald :min-count 2)
                    :ba ::behaviour-archive))

(defn selection
  "Inverse frequency dependent selection. Uses :beh-freq (counts) from
  population. Assigned weights, normalised over the population total,
  prescribe each individual's proportion of the selection pool. The returned
  pool of individuals is accumulated in order of weight decreasing, until it
  reaches the original population size."
  [epopn parameters]
  (let [n (count epopn)
        {:keys []} parameters
        wpopn (map (fn [indiv]
                     (assoc indiv :weight (/ 1.0 (+ 1 (:beh-freq indiv)))))
                   epopn)
        sumw (reduce + (map :weight wpopn))]
    (loop [indivs (sort-by :weight > wpopn)
           seln ()]
      (if (< (count seln) n)
        (let [indiv (first indivs)
              w (:weight indiv)
              take-n (-> (* n (/ w sumw))
                         (util/round)
                         (max 1))]
          (recur (rest indivs)
                 (into seln (repeat take-n indiv))))
        seln))))

(s/fdef selection
        :args (s/cat :epopn (s/every ::indiv-evald)
                     :parameters ::parameters)
        :ret (s/every ::indiv))

(defn next-generation
  "Take selected individuals, with repeats, and prepare the next generation
  using mutation and crossover"
  [seln beh-archive rng parameters]
  (let [{:keys [crossover-prob]} parameters
        crossable (map :representative (vals beh-archive))]
    (loop [indivs seln
           popn ()
           rng rng]
      (if-let [indiv (first indivs)]
        (let [genome (:genome indiv)
              [rng r1 r2 r3 r4] (random/split-n rng 5)
              offspring (->
                         (if (< (util/rand r1 0 1) crossover-prob)
                           (let [other (util/rand-nth r2 crossable)]
                             (grncre/crossover genome other r3))
                           genome)
                         (grncre/mutate r4))]
          (recur (rest indivs)
                 (conj popn {:genome offspring})
                 rng))
        ;; done
        popn))))

(s/fdef next-generation
        :args (s/cat :seln (s/every ::indiv)
                     :ba ::behaviour-archive
                     :rng ::util/rng
                     :parameters ::parameters)
        :ret (s/every ::indiv))

(defn generation
  [popn beh-archive parameters rng]
  (let [[epopn ba] (eval-popn popn beh-archive)
        seln (selection epopn parameters)
        next-popn (next-generation seln ba rng parameters)]
    [next-popn ba]))

(defn to-file
  [file form]
  (with-open [w (io/writer file)]
    (print-dup form w)))

(defn from-file
  [file]
  (with-open [r (io/reader file)]
    (read r)))

(defn evolve
  [rng parameters]
  (let [n-popn (:population-size parameters)
        [rng rng*] (random/split rng)]
    (loop [popn (map (fn [rng]
                       {:genome (grncre/random-genome rng)})
                     (random/split-n rng* n-popn))
           beh-archive {}
           gi 0
           rng rng]
      (if (< gi (:generations parameters))
        ;; this genereration tag will end up in the behaviour archive
        (let [popn (map #(assoc-in % [:genome :generation] gi) popn)
              [rng rng*] (random/split rng)
              [next-popn ba] (generation popn beh-archive parameters rng*)]
          (println "gen" gi "behs:" (count ba))
          (recur next-popn
                 ba
                 (inc gi)
                 rng))
        ;; done
        beh-archive))))

(defn run
  [seed & args]
  (let [rng (if seed
              (random/make-random (Long. (str seed)))
              (random/make-random))]
    (->> (evolve rng parameter-defaults)
         (to-file (str "beh-archive-seed" seed ".edn")))))

(comment
  (require '[org.nfrac.xotarium.evo :as evo] :reload)
  (require '[clojure.spec :as s])
  (require '[clojure.spec.test :as stest])
  (require '[clojure.spec.gen :as gen])
  (stest/instrument
   (concat
    (stest/enumerate-namespace 'org.nfrac.xotarium.evo)
    (stest/enumerate-namespace 'org.nfrac.xotarium.grn.greans)
    ))
  (evo/run "1")
  (stest/unstrument)

)
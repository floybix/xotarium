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
   :population-size 20
   :max-selection-fraction 0.5
   :generations 5})

(def sim-steps (* 60 5))

(def xy-beh-resolution 0.5);; each multiple of this is a distinct behaviour.

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

(defn mid-point
  [xys]
  (loop [xys xys
         xmin Double/MAX_VALUE
         ymin Double/MAX_VALUE
         xmax Double/MIN_VALUE
         ymax Double/MIN_VALUE]
    (if-let [[x y] (first xys)]
      (let [x (double x)
            y (double y)]
        (recur (rest xys)
               (if (< x xmin) x xmin)
               (if (< y ymin) y ymin)
               (if (> x xmax) x xmax)
               (if (> y ymax) y ymax)))
      [(* 0.5 (+ xmin xmax))
       (* 0.5 (+ ymin ymax))])))

(defn point-in-time-measures
  [state]
  (let [creature (:creature state)
        grps (apply concat (vals (:groups creature)))
        cents (map lf/center grps)
        [x y] (mid-point cents)]
    {:cent-x (int (/ x xy-beh-resolution))
     :cent-y (int (/ y xy-beh-resolution))}))

(defn behaviour-measures
  [genome seed]
  (when-let [state (grncre/setup genome seed)]
    (let [state-ht (->> state
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
               [:ft beh-ft]]))))

(s/def ::behaviour
  (s/keys))

(s/def ::beh-freq nat-int?)

(s/def ::indiv
  (s/keys :req-un [::grncre/genome]))

(s/def ::indiv-evald
  (s/keys :req-un [::grncre/genome
                   ::behaviour
                   ::beh-freq]))

(defn genome-hash-str
  [genome-hash]
  (format "%h" genome-hash))

(defn eval-popn
  [popn behaviour-archive seed]
  ;; TODO: parallelize?
  (let [bpopn (map (fn [indiv]
                     (let [beh (behaviour-measures (:genome indiv) seed)]
                       (println (format "%16s" (genome-hash-str (hash (:genome indiv))))
                                "  beh=" beh)
                       (assoc indiv :behaviour beh)))
                   (sort-by #(hash (:genome %)) popn))
        popn-freqs (frequencies (map :behaviour bpopn))
        ;; add behaviours to archive
        [epopn ba] (loop [indivs bpopn
                          i 0
                          ba behaviour-archive
                          epopn ()]
                     (if-let [indiv (first indivs)]
                       (let [beh (:behaviour indiv)
                             prior-freq (get-in behaviour-archive [beh :count] 0)
                             beh-freq (+ (get popn-freqs beh) prior-freq)
                             nov? (zero? prior-freq)
                             bad? (nil? beh)]
                         (recur (rest indivs)
                                (inc i)
                                (cond
                                  bad?
                                  ba
                                  nov?
                                  (assoc ba beh
                                         {:count beh-freq
                                          :representative (:genome indiv)})
                                  :else
                                  (update-in ba [beh :count] inc))
                                (conj epopn
                                      (if bad?
                                        indiv
                                        (assoc indiv :beh-freq beh-freq)))))
                       ;; done
                       [epopn ba]))
        ]
    [epopn ba]))

(s/fdef eval-popn
        :args (s/cat :popn (s/every ::indiv :min-count 2)
                     :ba ::behaviour-archive
                     :seed int?)
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
        {:keys [max-selection-fraction]} parameters
        sumw* (->> epopn
                   (map :beh-freq)
                   (map (fn [n] (if n
                                  (/ 1.0 (+ 1 n))
                                  0.0)))
                   (reduce +))
        wpopn (map (fn [indiv]
                     (let [w (if-let [n (:beh-freq indiv)]
                               (-> (/ 1.0 (+ 1 n))
                                   (min (* max-selection-fraction sumw*)))
                               ;; invalid, set weight to zero
                               0.0)]
                       (assoc indiv :weight w)))
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
  [seln beh-archive rng parameters gi]
  (let [{:keys [crossover-prob]} parameters
        crossable (concat (map :representative (vals beh-archive))
                          (map :genome seln))]
    (loop [indivs seln
           popn ()
           gparents {}
           rng rng]
      (if-let [indiv (first indivs)]
        (let [genome (:genome indiv)
              [rng r1 r2 r3 r4] (random/split-n rng 5)
              [mut-genome mut-info] (grncre/mutate genome r1)
              [new-genome parent-ids] (if (< (random/rand-double r2) crossover-prob)
                                        (let [other (util/rand-nth r3 crossable)]
                                          [(grncre/crossover mut-genome other r4)
                                           [(hash genome) (hash other)]])
                                        [mut-genome
                                         [(hash genome)]])
              ;; this genereration tag will end up in the behaviour archive
              new-genome (assoc new-genome :generation gi)]
          (recur (rest indivs)
                 (conj popn {:genome new-genome
                             :mutation-info mut-info})
                 (assoc gparents (hash new-genome) parent-ids)
                 rng))
        ;; done
        [popn gparents]))))

(s/fdef next-generation
        :args (s/cat :seln (s/every ::indiv)
                     :ba ::behaviour-archive
                     :rng ::util/rng
                     :parameters ::parameters
                     :gi nat-int?)
        :ret (s/cat :popn (s/every ::indiv)
                    :gparents (s/every-kv any? (s/every any?))))

(defn generation
  [popn beh-archive parameters seed rng gi]
  (let [[epopn ba] (eval-popn popn beh-archive seed)
        seln (selection epopn parameters)
        [next-popn gparents] (next-generation seln ba rng parameters gi)]
    [next-popn ba gparents]))

(defn to-file
  [file form]
  (with-open [w (io/writer file)]
    (print-dup form w)))

(defn from-file
  [file]
  (with-open [r (java.io.PushbackReader. (io/reader file))]
    (read r)))

(defn evolve-continue
  [rng seed parameters popn beh-archive parents gi]
  (if (< gi (:generations parameters))
    (let [[rng rng*] (random/split rng)
          [next-popn ba gparents] (generation popn beh-archive parameters
                                              seed rng* gi)]
      (println "gen" gi "behs:" (count ba))
      (println "grn sizes:"
               (->> popn (map :genome) (map :grn) (map #(count (::grn/elements %))) (sort >)))
      (println "cppn sizes:"
               (->> popn (map :genome) (map :cppn) (map #(count (cppn/edge-list %))) (sort >)))
      (println "new gen:")
      (doseq [indiv (sort-by #(gparents (hash (:genome %))) next-popn)
              :let [genome (:genome indiv)]]
        (println (genome-hash-str (hash genome)) ":"
                 "parents" (map genome-hash-str (gparents (hash genome)))
                 "mutations" (:mutation-info indiv)))
      (recur rng seed parameters next-popn ba
             (merge-with concat parents gparents) ;; identical mutants might show up
             (inc gi)))
    ;; done
    {:beh-archive beh-archive
     :popn popn
     :parents parents
     :seed seed
     :generation gi}))

(defn evolve
  [rng seed parameters]
  (let [n-popn (:population-size parameters)
        [rng rng*] (random/split rng)]
    (evolve-continue rng seed parameters
                     (map (fn [rng]
                            {:genome (grncre/random-genome rng)})
                          (random/split-n rng* n-popn))
                     {} {} 0)))

(defn run
  [ng seed & args]
  (let [ng (Long. (str ng))
        seed (Long. (str seed))
        rng (random/make-random seed)]
    (->> (evolve rng seed (assoc parameter-defaults
                                 :generations ng))
         (to-file (str "beh-archive-seed" seed ".edn")))))

(defn run-more
  [file ng run-seed & args]
  (let [ng (Long. (str ng))
        run-seed (Long. (str run-seed))
        {:keys [beh-archive popn parents seed generation]} (from-file file)
        rng (random/make-random (Long. (str seed)))]
    (->> (evolve-continue rng run-seed
                          (assoc parameter-defaults
                                 :generations (+ ng generation))
                          popn beh-archive parents generation)
         (to-file (str file "-more-" run-seed)))))

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
  (evo/run 3 1)
  (stest/unstrument)

)

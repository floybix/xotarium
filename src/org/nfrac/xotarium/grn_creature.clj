(ns org.nfrac.xotarium.grn-creature
  (:require [org.nfrac.xotarium.creature :as cre]
            [org.nfrac.xotarium.cppn :as cppn]
            [org.nfrac.xotarium.cppn-compile :as cc]
            [org.nfrac.xotarium.grn.greans :as grn]
            [org.nfrac.xotarium.plant-cave :as cave]
            [org.nfrac.liquidfun.testbed :as bed]
            [org.nfrac.liquidfun.core :as lf :refer [body! joint!
                                                     particle-system!]]
            [org.nfrac.liquidfun.vec2d :as v2d]
            [org.nfrac.xotarium.util.algo-graph :as graph]
            [quil.core :as quil :include-macros true]
            [quil.middleware]
            [clojure.pprint]
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

(def beh-inputs '[bias
                  oscillator
                  factor-a
                  factor-b
                  factor-c])
#_[carn-touch
   herb-touch
   plant-touch
   wall-touch
   self-touch
   carn-smell
   herb-smell
   plant-smell
   wall-smell
   message-a
   message-b
   message-c]

(def beh-outputs '[expansion])
#_[messenger-a
   messenger-b
   messenger-c
   oscillator-freq]

(s/def ::genome
  (s/keys :req-un [::cppn/cppn
                   ::grn/grn]))

(defn random-genome
  [rng]
  (let [cppn cre/seed-cppn
        grn (grn/random-grn (count beh-inputs) (count beh-outputs) rng)]
    {:cppn cppn
     :grn grn}))

(defn mutate
  [genome]
  (let [{:keys [grn cppn]} genome]
    (assoc genome
           :cppn (cppn/mutate-general cppn)
           :grn (grn/mutate grn))))

(defn crossover
  [g1 g2]
  (assoc g1 :grn (grn/crossover (:grn g1) (:grn g2))))

(defn come-alive
  [world genome]
  (let [{:keys [grn cppn]} genome
        creature-body (cre/make-creature world cppn)
        tri-p (:triad-params creature-body)
        cell-form (grn/grn->cell grn)
        tri-concs (zipmap (keys tri-p)
                          (repeat (::grn/concs cell-form)))]
    (assoc creature-body
           :tri-concs tri-concs
           :grn-cell cell-form
           :grn grn)))

(defn setup
  [rng]
  (let [;world (cave/build-world)
        world (-> (cave/setup) (cave/do-add-air) :world)
        ps ^liquidfun$b2ParticleSystem (first (lf/particle-sys-seq world))
        [rng rng*] (random/split rng)
        genome (random-genome rng*)
        creature (come-alive world genome)]
      (assoc bed/initial-state
             :world world
             :particle-system ps
             :creature creature
             :rng rng
             :particle-iterations 3
             :camera (bed/map->Camera {:width cave/cave-width
                                       :height cave/cave-height
                                       :center [0 0]}))))

(defn post-step
  [state]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)
        creature (:creature state)
        tri-p (:triad-params creature)
        cell-form (:grn-cell creature)
        tri-concs (:tri-concs creature)
        dt 0.05
        time (:time state)
        freq 6.0
        phase (mod (* time freq) (* 2.0 Math/PI))
        new-concs (persistent!
                   (reduce-kv (fn [m handles concs]
                                (let [params (get tri-p handles)
                                      cell (assoc cell-form ::grn/concs concs)
                                      phase-off (:phase-off params)
                                      osc (if (pos? (Math/sin (+ phase phase-off)))
                                            1.0 0.0)
                                      ic [1.0 ; bias
                                          osc
                                          (:factor-a params)
                                          (:factor-b params)
                                          (:factor-c params)]
                                      nc (::grn/concs
                                          (grn/step cell ic dt))]
                                  (assoc! m handles nc)))
                              (transient {})
                              tri-concs))
        work-fn (fn [handles params]
                  (let [concs (get tri-concs handles)
                        cell (assoc cell-form ::grn/concs concs)
                        g-out (grn/cell-outputs cell)
                        expan (-> (first g-out) (* 2.0) (- 1.0))]
                    expan))]
    (cre/creature-flex ps (:triad-params creature) work-fn)
    (cre/groups-restore-color ps (:muscle (:groups creature)) [255 0 0 255])
    (assoc-in state [:creature :tri-concs] new-concs)))

(defn step
  [state]
  (-> (bed/world-step state)
      (post-step)))

(defn my-key-press
  [state event]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)]
    (case (:key event)
      :b (do
           (body! (:world state) {}
                  {:shape (lf/circle 0.25)
                   :restitution 0.1
                   :density 1.0})
           state)
      :g (do
           (.SetGravityScale ps (if (== 1.0 (.GetGravityScale ps))
                                  0.0 1.0))
           state)
      :d (do
           (.SetDamping ps (if (== 1.0 (.GetDamping ps))
                             0.0 1.0))
           state)
      ;; otherwise pass on to testbed
      (bed/key-press state event))))

(defn run
  [seed & args]
  (let [rng (if seed
              (random/make-random (Long. (str seed)))
              (random/make-random))]
    (quil/sketch
     :title "Xotarium"
     :host "liquidfun"
     :setup #(setup rng)
     :update (fn [s] (if (:paused? s) s (step s)))
     :draw #(if (zero? (mod (quil/frame-count) 2)) (bed/draw % true) %)
     :key-typed my-key-press
     :mouse-pressed bed/mouse-pressed
     :mouse-released bed/mouse-released
     :mouse-dragged bed/mouse-dragged
     :mouse-wheel bed/mouse-wheel
     :size [600 500]
     :features [:resizable]
     :middleware [quil.middleware/fun-mode])))

(comment
  (require '[org.nfrac.xotarium.util :as util])
  (require '[org.nfrac.xotarium.grn-creature])
  (in-ns 'org.nfrac.xotarium.grn-creature)
  (def state (setup (random/make-random 1)))
  (def steps (take 100 (iterate step state)))
  (def all-handles (-> steps first :creature :tri-concs keys))
  (def rng (random/make-random 2))
  (def use-handles (set (take 5 (util/shuffle rng all-handles))))
  (def data
    (for [state steps
          :let [t (:time state)
                tri-concs (get-in state [:creature :tri-concs])]
          [handles concs] tri-concs
          :when (contains? use-handles handles)
          :let [triad-id (str (hash handles))]
          [i conc] (map-indexed vector concs)]
      {:x t
       :triad-id triad-id
       :tf-id (str i)
       :y conc}))
  (use 'gyptis.core)
  (require '[gyptis.view :refer [plot!]])
  (-> (map #(assoc % :stroke (:triad-id %))
           data)
      (line)
      (facet-grid {:facet_y :tf-id})
      (plot!))

)

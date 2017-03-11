(ns org.nfrac.xotarium.replay
  (:require [org.nfrac.xotarium.evo :as evo]
            [org.nfrac.xotarium.grn-creature :as grncre]
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

(defn to-file
  [file form]
  (with-open [w (io/writer file)]
    (print-dup form w)))

(defn from-file
  [file]
  (with-open [r (java.io.PushbackReader. (io/reader file))]
    (read r)))

(defn my-key-press
  [state event]
  (let [s (:current state)
        ps ^liquidfun$b2ParticleSystem (:particle-system s)]
    (case (:key event)
      :n (let [i (:genome-index state)
               i+ (inc i)
               genome (get (:genomes state) i+)
               seed (:seed state)]
           (when-not genome
             (System/exit 0))
           (assoc state
                  :genome-index i+
                  :current (grncre/setup genome seed)))
      :b (do
           (body! (:world s) {}
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
      (update state :current bed/key-press event))))

(defn setup
  [file]
  (let [{:keys [beh-archive seed]} (from-file file)
        genomes (->> (map :representative (vals beh-archive))
                     (sort-by :generation)
                     (vec))
        genome-index 0]
    {:seed seed
     :genomes genomes
     :genome-index genome-index
     :current (grncre/setup (get genomes genome-index) seed)}))

(defn draw
  [state]
  (bed/draw (:current state) true)
  (let [gi (:genome-index state)
        genome (get (:genomes state) gi)]
    (quil/fill 255)
    (quil/text (str "Replaying beh. "
                    (inc gi)
                    " / "
                    (count (:genomes state))
                    " (gen " (:generation genome) ")."
                    " Keys: (n) next behaviour")
               10 10)))

(defn run
  [file & args]
  (quil/sketch
   :title "Xotarium"
   :host "liquidfun"
   :setup #(do
             (quil/frame-rate 30)
             (setup file))
   :update (fn [state]
             (if (:paused? (:current state))
               state
               (update state :current grncre/step)))
   :draw (fn [state]
           (if (zero? (mod (quil/frame-count) 2))
             (draw state)
             state))
   :key-typed my-key-press
   :mouse-pressed #(update % :current bed/mouse-pressed %2)
   :mouse-released #(update % :current bed/mouse-released %2)
   :mouse-dragged #(update % :current bed/mouse-dragged %2)
   :mouse-wheel #(update % :current bed/mouse-wheel %2)
   :size [600 500]
   :features [:resizable]
   :middleware [quil.middleware/fun-mode]))

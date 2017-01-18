(ns org.nfrac.xotarium.hills
  "An experiment with building hill-like obstacles."
  (:require [org.nfrac.liquidfun.testbed :as bed]
            [org.nfrac.liquidfun.core :as lf :refer [body! joint!
                                                     particle-system!]]
            [org.nfrac.liquidfun.vec2d :as v2d]
            [quil.core :as quil :include-macros true]
            [quil.middleware])
  (:import (org.bytedeco.javacpp
            liquidfun$b2ContactListener
            liquidfun$b2ParticleContact
            liquidfun$b2ParticleSystem
            liquidfun$b2ParticleGroup
            liquidfun$b2ParticleDef
            liquidfun$b2QueryCallback)))

(def p-radius 0.05)
(def cave-width 10.0)
(def cave-height 6.0)
(def cave-hw (* cave-width 0.5))

(def dump-n 20)
(def block-size-min 0.1)
(def block-size-max 0.8)
(def dump-velocity 4.0)
(def dump-time 0.25)
(def spray-time 5.0)
(def spray-dn 2)
(def spray-velocity 40.0)
(def spray-from-top 2.5)
(def spray-size 0.5)

(defn rand-in [min max] (+ min (rand (- max min))))

(defn setup []
  (let [world (lf/new-world)
        hw (* 0.5 cave-width)
        ground (body! world {:type :static}
                      {:shape (lf/edge-loop [[(- hw) 0]
                                             [(- hw) cave-height]
                                             [hw cave-height]
                                             [hw 0]])})
        ps (particle-system! world
                             {:radius p-radius
                              :density 2.0
                              :pressure-strength 1.0
                              :destroy-by-age false})
        pdef (lf/particle-def {:flags #{:water}
                               :color [255 255 255 255]})
        its (.CalculateReasonableParticleIterations world (/ 1 60.0))]
    (println "reasonable particle iterations:" its)
    ;(.SetContactListener world lstnr)
    (assoc bed/initial-state
      :world world
      :particle-system ps
      :particle-iterations its
      ::pdef pdef
      :dt-secs (/ 1 60.0)
      :camera (bed/map->Camera {:width 11 :height 6 :center [0 4]}))))

(defn step-spray
  [state]
  (let [world (:world state)
        t (:time state)
        dt (:dt-secs state)
        {:keys [t0]} (::spray state)
        pdef ^liquidfun$b2ParticleDef (::pdef state)
        ps ^liquidfun$b2ParticleSystem (:particle-system state)]
    (if (> t (+ t0 spray-time))
      (dissoc state ::spray)
      (do
        (dotimes [j spray-dn]
          (let [x (- cave-hw p-radius)
                y (- cave-height spray-from-top (rand spray-size))]
            (.SetPosition pdef x y)
            (.Set (.velocity pdef) (- spray-velocity) 0)
            (lf/particle! ps pdef)))
        state))))

(defn step-dump
  [state]
  (let [world (:world state)
        t (:time state)
        dt (:dt-secs state)
        {:keys [t0 ti x]} (::block-dump state)]
    (if (> t (+ t0 (* dump-time dump-n)))
      (dissoc state ::block-dump)
      (if (> t (+ ti dump-time))
        (let [w (rand-in block-size-min block-size-max)
              h (rand-in block-size-min block-size-max)]
          (body! world {:position [x (- cave-height (max w h))]
                        :angle (rand (* 2 Math/PI))
                        :linear-velocity [0 (- dump-velocity)]}
                 {:shape (lf/box (/ w 2) (/ h 2))
                  :density 1.0
                  :friction 1.0})
          (assoc-in state [::block-dump :ti] t))
        state))))

(defn post-step
  [state]
  (cond->
   state
   (::block-dump state) (step-dump)
   (::spray state) (step-spray)))

(defn step
  [state]
  (-> (bed/world-step state)
      (post-step)))
      ;(bed/record-snapshot true)))

(defn draw
  [state]
  (bed/draw state #_true)
  (let []
    (quil/fill 255)
    (quil/text (str "Keys: (g) toggle gravity (d) damping\n"
                    "      (b) dump blocks (s) spray (p) print coords")
               10 10)))

(defn my-key-press
  [state event]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)]
    (case (:key event)
      :b (let [x (- (rand cave-width) cave-hw)]
           (assoc state ::block-dump {:t0 (:time state)
                                      :ti 0.0
                                      :x x}))
      :s (assoc state ::spray {:t0 (:time state)})
      :g (do
           (.SetGravityScale ps (if (== 1.0 (.GetGravityScale ps))
                                  0.0 1.0))
           state)
      :d (do
           (.SetDamping ps (if (== 1.0 (.GetDamping ps))
                             0.0 1.0))
           state)
      :2 (do
           (.SetRadius ps (* 1.5 (.GetRadius ps)))
           state)
      :1 (do
           (.SetRadius ps (* (/ 1 1.5) (.GetRadius ps)))
           state)
      ;; otherwise pass on to testbed
      (bed/key-press state event))))

(defn ^:export run
  "Run the test sketch."
  [& args]
  (quil/sketch
   :title "Hills"
   :host "liquidfun"
   :setup setup
   :update (fn [s] (if (:paused? s) s (step s)))
   :draw draw
   :key-typed my-key-press
   :mouse-pressed bed/mouse-pressed
   :mouse-released bed/mouse-released
   :mouse-dragged bed/mouse-dragged
   :mouse-wheel bed/mouse-wheel
   :size [600 500]
   :features [:resizable]
   :middleware [quil.middleware/fun-mode]))

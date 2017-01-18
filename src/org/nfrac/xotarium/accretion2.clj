(ns org.nfrac.xotarium.accretion2
  "An experiment with accreting wall particles."
  (:require [org.nfrac.liquidfun.testbed :as bed]
            [org.nfrac.liquidfun.core :as lf :refer [body! joint!
                                                     particle-system!]]
            [org.nfrac.liquidfun.vec2d :as v2d]
            [quil.core :as quil :include-macros true]
            [quil.middleware])
  (:import (org.bytedeco.javacpp
            liquidfun$b2ContactFilter
            liquidfun$b2Fixture
            liquidfun$b2Body
            liquidfun$b2Transform
            liquidfun$b2ParticleSystem
            liquidfun$b2ParticleGroup
            liquidfun$b2ParticleDef
            liquidfun$b2QueryCallback)))

(def p-radius 0.025)
(def cave-width 5.0)
(def cave-height 5.0)
(def cave-hw (* cave-width 0.5))
(def min-group-size 5)
(def circ-force 0.005)

(def accum-hits (atom {}))

(gen-class
 :name org.nfrac.liquidfun.testbed.tests.accretion2.WallFilter
 :extends org.bytedeco.javacpp.liquidfun$b2ContactFilter
 :state state
 :init init
 :prefix "impl-")

;; getting IndexOutOfBoundsException if pass state arg to -init. workaround:
(let [state (atom nil)]

  (defn impl-init []
    [[] @state])

  (defn wall-filter
    [wall-pg ground]
    (reset! state [wall-pg ground])
    (org.nfrac.liquidfun.testbed.tests.accretion2.WallFilter.)))

(defn impl-ShouldCollide-b2Fixture-b2ParticleSystem-int
  [^org.nfrac.liquidfun.testbed.tests.accretion2.WallFilter this
   ^liquidfun$b2Fixture fixture ^liquidfun$b2ParticleSystem ps i]
  (let [[wall-pg ground] (.state this)
        wall-pg ^liquidfun$b2ParticleGroup wall-pg
        ground ^liquidfun$b2Body ground]
    (when (and (not (.ContainsParticle wall-pg i))
               (= (lf/body-of fixture) ground))
      (let [x (.GetParticlePositionX ps i)
            y (.GetParticlePositionY ps i)]
        (swap! accum-hits assoc i [x y]))))
  true)

(defn impl-ShouldCollide-b2ParticleSystem-int-int
  [^org.nfrac.liquidfun.testbed.tests.accretion2.WallFilter this
   ^liquidfun$b2ParticleSystem ps ia ib]
  (let [[wall-pg _] (.state this)
        wall-pg ^liquidfun$b2ParticleGroup wall-pg]
    (let [a-wall? (.ContainsParticle wall-pg ia)
          b-wall? (.ContainsParticle wall-pg ib)
          [walli gasi] (cond
                         (and a-wall? b-wall?) nil
                         a-wall? [ia ib]
                         b-wall? [ib ia]
                         :else nil)]
      (when gasi
        (let [x (.GetParticlePositionX ps gasi)
              y (.GetParticlePositionY ps gasi)]
          (swap! accum-hits assoc gasi [x y])))))
  true)

(defn sign [x] (if (neg? x) -1 1))
(defn rand-in [min max] (+ min (rand (- max min))))

(defn setup []
  (let [world (lf/new-world)
        hw (* 0.5 cave-width)
        hh (* 0.5 cave-height)
        ground (body! world {:type :static}
                      {:shape (lf/edge-loop [[(- hw) (- hh)]
                                             [(- hw) hh]
                                             [hw hh]
                                             [hw (- hh)]])})
        ps (particle-system! world
                             {:radius p-radius
                              :pressure-strength 1.0
                              :damping-strength 1.0
                              :gravity-scale 0.0
                              ;:strict-contact-check true
                              :destroy-by-age false})
        gas-pg (lf/particle-group!
                ps
                {:flags #{:water :particle-contact-filter :fixture-contact-filter}
                 :group-flags #{:can-be-empty}
                 :stride (* 0.75 p-radius 2.0)
                 :color [255 128 128 128]
                 :shape (lf/box (* hw 0.25) (* hh 0.25)
                                [0 0])})
        wall-pg (lf/particle-group!
                 ps
                 {:flags #{:wall :barrier :particle-contact-filter}
                  :group-flags #{:can-be-empty}
                  :color [255 255 255 255]})
        de-filter (.m_contactFilter (.GetContactManager world))
        filt (wall-filter wall-pg ground)
        pdef (lf/particle-def {:flags #{:wall :barrier :particle-contact-filter}
                               :color [255 255 255 255]
                               :group wall-pg})
        its (.CalculateReasonableParticleIterations world (/ 1 60.0))]
    (println "reasonable particle iterations:" its)
    (.SetContactFilter world filt)
    ;; clear out the center of gas as particles can end up fixed there
    (.DestroyParticlesInShape ps (lf/box (* p-radius 3) (* p-radius 3))
                              (doto (liquidfun$b2Transform.) (.SetIdentity)))
    ;; add some random edges to disturb the gas flow
    (doseq [[sx sy] [[1 1] [1 -1] [-1 1] [-1 -1]]]
      (let [x0 (* (rand-in (* 0.25 hw) hw) sx)
            x1 (* (rand-in (* 0.25 hw) hw) sx)
            y0 (* (rand-in (* 0.25 hh) hh) sy)
            y1 (* (rand-in (* 0.25 hh) hh) sy)]
        (lf/fixture! ground
               {:shape (lf/edge [x0 y0] [x1 y1])})))
    (assoc bed/initial-state
      :world world
      :contact-filter filt
      ::filtering? true
      ::de-filter de-filter
      :settings {:wall-pg wall-pg
                 :gas-pg gas-pg
                 :pdef pdef
                 :circ? true}
      :particle-system ps
      :particle-iterations its
      :dt-secs (/ 1 60.0)
      :camera (bed/map->Camera {:width 6 :height 6 :center [0 0]}))))

(defn post-step
  [state]
  (let [dt (:dt-secs state)
        ps ^liquidfun$b2ParticleSystem (:particle-system state)
        {:keys [wall-pg gas-pg circ?] :as settings} (:settings state)
        gas-pg ^liquidfun$b2ParticleGroup gas-pg
        pdef ^liquidfun$b2ParticleDef (:pdef settings)
        v2tmp (lf/vec2 0 0)
        gas-i (.GetBufferIndex gas-pg)]
    (when (and circ? (::filtering? state))
      (dotimes [j (.GetParticleCount gas-pg)]
        (let [i (+ gas-i j)
              x (.GetParticlePositionX ps i)
              y (.GetParticlePositionY ps i)
              xz (/ x cave-width)
              yz (/ y cave-height)]
          (.Set v2tmp
                (* circ-force -2.0 yz)
                (* circ-force 2.0 xz))
          (.ParticleApplyForce ps i v2tmp)))
      (let [gas-done? (atom false)]
        (doseq [[i [x y]] @accum-hits]
          (when (== 1 (.GetParticleCount gas-pg))
            (reset! gas-done? true))
          (.DestroyParticle ps i)
          (.SetPosition pdef x y)
          (lf/particle! ps pdef))
        (swap! accum-hits empty)
        (if @gas-done?
          ;; turn off accretion listener
          (let [de-filt (::de-filter state)]
            (.SetContactFilter (:world state) de-filt)
            (assoc state :contact-filter de-filt
                   ::filtering? false))
          state)))))

(defn step
  [state]
  (cond-> (bed/world-step state)
          (::filtering? state) (post-step)))
      ;(bed/record-snapshot true)))

(defn draw
  [state]
  (bed/draw state true)
  (let []
    (quil/fill 255)
    (quil/text (str "Keys: (g) toggle gravity (d) damping (c) circulation\n"
                    " (s) split (c) color (z) zap groups (p) print coords")
               10 10)))

(defn rand-color []
  [(int (rand-in 64 255))
   (int (rand-in 64 255))
   (int (rand-in 64 255))])

(defn particle-indices
  [^liquidfun$b2ParticleGroup pg]
  (let [i0 (.GetBufferIndex pg)]
    (range i0 (+ i0 (.GetParticleCount pg)))))

(defn my-key-press
  [state event]
  (let [world (:world state)
        pdef ^liquidfun$b2ParticleDef (get-in state [:settings :pdef])
        ps ^liquidfun$b2ParticleSystem (:particle-system state)
        wall-pg ^liquidfun$b2ParticleGroup (get-in state [:settings :wall-pg])]
    (case (:key event)
      :s (do
           (let [g wall-pg]
             (.SetGroupFlags g 0)
             (when (> (.GetParticleCount g) 1)
               (.SplitParticleGroup ps g)))
           state)
      :c (do
           (doseq [group (lf/particle-group-seq ps)]
             (let [colb (.GetColorBuffer ps)
                   [r g b] (rand-color)]
               (doseq [i (particle-indices group)]
                 (let [coli (.position colb (long i))]
                   (.Set coli r g b 255)))))
           state)
      :z (do
           (doseq [g (lf/particle-group-seq ps)]
             (when (< (.GetParticleCount g) min-group-size)
               (doseq [i (particle-indices g)]
                 (lf/destroy-particle! ps i))))
           state)
      :g (do
           (.SetGravityScale ps (if (== 1.0 (.GetGravityScale ps))
                                  0.0 1.0))
           state)
      :d (do
           (.SetDamping ps (if (== 1.0 (.GetDamping ps))
                             0.0 1.0))
           state)
      :p (let [wall-i (.GetBufferIndex wall-pg)]
           (println "[")
           (dotimes [j (.GetParticleCount wall-pg)]
             (print (str " [" (.GetParticlePositionX ps (+ j wall-i))
                         "," (.GetParticlePositionY ps (+ j wall-i))
                         "]")))
           (println "]")
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
   :title "Accretion2"
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

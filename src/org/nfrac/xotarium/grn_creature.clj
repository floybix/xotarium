(ns org.nfrac.xotarium.grn-creature
  (:require [org.nfrac.xotarium.creature :as cre]
            [org.nfrac.xotarium.cppn :as cppn]
            [org.nfrac.xotarium.cppn-compile :as cc]
            [org.nfrac.xotarium.grn.greans :as grn]
            [org.nfrac.xotarium.plant-cave :as cave]
            [org.nfrac.xotarium.proximity-field :as proxf]
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
            liquidfun$b2ParticleColor
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

(def max-sensed-velocity 2.0)

;; this list should match the args passed to grn/step!
(def beh-inputs '[bias
                  oscillator
                  factor-a
                  factor-b
                  factor-c
                  wall-touch
                  wall-smell
                  self-touch
                  other-touch
                  red-smell
                  green-smell
                  blue-smell
                  pressure
                  up-vel
                  down-vel
                  right-vel
                  left-vel
                  ])
#_[message-a
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
  [genome rng]
  (let [{:keys [grn cppn]} genome
        [rng rng*] (random/split rng)]
    (assoc genome
           :cppn (cppn/mutate-with-perturbation cppn rng* {})
           :grn (grn/mutate grn rng))))

(defn crossover
  [g1 g2 rng]
  (assoc g1 :grn (grn/crossover (:grn g1) (:grn g2) rng)))

(defn come-alive
  [world genome]
  (when-let [creature-body (cre/make-creature world (:cppn genome))]
    (let [grn (:grn genome)
          tri-p (:triad-params creature-body)
          cell-form (grn/grn->cell grn)
          tri-concs (zipmap (keys tri-p)
                            (repeat (::grn/concs cell-form)))]
      (assoc creature-body
            :tri-concs tri-concs
            :grn-cell cell-form
            :grn grn))))

(defn setup
  [genome seed]
  (let [;world (cave/build-world)
        state (-> (cave/setup)
                  (cave/do-add-air))
        _ (cave/add-random-static-bars (:ground state) (random/make-random seed))
        pf (cave/get-prox-field (:world state))
        creature (come-alive (:world state) genome)]
    (when creature
      (assoc state
            :wall-prox-field pf
            :creature creature))))

(defn clamp [x] (-> (double x) (max 0.0) (min 1.0)))

(defn post-step
  [state]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)
        creature (:creature state)
        my-muscles (:muscle (:groups creature))
        my-bones (:bone (:groups creature))
        tri-p (:triad-params creature)
        cell-form (:grn-cell creature)
        tri-concs (:tri-concs creature)
        dt 0.05
        time (:time state)
        freq 6.0
        phase (mod (* time freq) (* 2.0 Math/PI))
        pf (:wall-prox-field state)
        y-lo (:y-lower-bound state)
        y-hi (:y-upper-bound state)
        y-span (- y-hi y-lo)
        air-pg ^liquidfun$b2ParticleGroup (::cave/air-pg state)
        bcm (let [n (.GetBodyContactCount ps)]
              (loop [bc (.GetBodyContacts ps)
                     bci 0
                     bcm (transient {})]
                (if (< bci n)
                  (let [i (.index bc)
                        body (.body bc)
                        bc (.position bc (inc bci))]
                    (if (.ContainsParticle air-pg i)
                      (recur bc (inc bci) bcm)
                      (recur bc (inc bci)
                             (assoc! bcm i body))))
                  (persistent! bcm))))
        ;; note: do not use GetGroupBuffer, it is slower than ContainsParticle
        colbuf (.GetColorBuffer ps)
        pcm (let [n (.GetContactCount ps)
                  pkind (fn [i]
                          (if (.ContainsParticle air-pg i)
                            :air
                            (cond
                              (some #(.ContainsParticle ^liquidfun$b2ParticleGroup % i)
                                    my-muscles)
                              :self-muscle
                              (some #(.ContainsParticle ^liquidfun$b2ParticleGroup % i)
                                    my-bones)
                              :self-bone
                              :else
                              :other)))
                  record-contact (fn [m kind ^liquidfun$b2ParticleColor color]
                                   (if (= kind :air)
                                     (let [r (/ (.r color) 255.0)
                                           g (/ (.g color) 255.0)
                                           b (/ (.b color) 255.0)]
                                       (update m :colors conj [r g b]))
                                     (assoc m kind true)))]
              (loop [pc (.GetContacts ps)
                     pci 0
                     pcm (transient {})]
                (if (< pci n)
                  (let [ia (.GetIndexA pc)
                        ib (.GetIndexB pc)
                        ka (pkind ia)
                        kb (pkind ib)
                        pc (.position pc (inc pci))]
                    (recur pc (inc pci)
                           (cond-> pcm
                             (= :self-muscle ka)
                             (assoc! pcm ia
                                     (record-contact (get pcm ia) kb
                                                     (.position colbuf ib)))
                             (= :self-muscle kb)
                             (assoc! pcm ib
                                     (record-contact (get pcm ib) ka
                                                     (.position colbuf ia))))))
                  (persistent! pcm))))
        max-vel max-sensed-velocity
        velb (.GetVelocityBuffer ps)
        nconcs (persistent!
                (reduce-kv
                 (fn [m handles concs]
                   (let [params (get tri-p handles)
                         [^liquidfun$b2ParticleHandle h1
                          ^liquidfun$b2ParticleHandle h2
                          ^liquidfun$b2ParticleHandle h3] handles
                         pi1 (.GetIndex h1)
                         pi2 (.GetIndex h2)
                         pi3 (.GetIndex h3)
                         x (.GetParticlePositionX ps pi1)
                         y (.GetParticlePositionY ps pi1)
                         pressure (-> (- y y-lo) (/ y-span) (clamp))
                         [vx vy] (lf/v2xy (.position velb pi1))
                         up-vel (-> vy (/ max-vel) (clamp))
                         down-vel (-> (- vy) (/ max-vel) (clamp))
                         right-vel (-> vx (/ max-vel) (clamp))
                         left-vel (-> (- vx) (/ max-vel) (clamp))
                         wall-smell (proxf/proximity-score pf x y)
                         wall-touch (if (or (bcm pi1) (bcm pi2) (bcm pi3))
                                      1.0 0.0)
                         pc1 (pcm pi1)
                         pc2 (pcm pi2)
                         pc3 (pcm pi3)
                         self-touch (if (or (:self-muscle pc1) (:self-muscle pc2)
                                            (:self-muscle pc3)) 1.0 0.0)
                         other-touch (if (or (:other pc1) (:other pc2)
                                             (:other pc3)) 1.0 0.0)
                         colors (concat (:colors pc1) (:colors pc2) (:colors pc3))
                         [rsum gsum bsum] (reduce (fn [[or og ob] [r g b]]
                                                    [(+ or r) (+ og g) (+ ob b)])
                                                  [0.0 0.0 0.0]
                                                  colors)
                         cell (assoc cell-form ::grn/concs concs)
                         phase-off (:phase-off params)
                         osc (if (pos? (Math/sin (+ phase phase-off)))
                               1.0 0.0)
                         ic [1.0 ; bias
                             osc
                             (:factor-a params)
                             (:factor-b params)
                             (:factor-c params)
                             wall-touch
                             wall-smell
                             self-touch
                             other-touch
                             (/ rsum (max 1 (count colors)))
                             (/ gsum (max 1 (count colors)))
                             (/ bsum (max 1 (count colors)))
                             pressure
                             up-vel
                             down-vel
                             right-vel
                             left-vel]
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
    (assoc-in state [:creature :tri-concs] nconcs)))

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
  (let [seed (Long. (str seed))
        rng (random/make-random seed)]
    (quil/sketch
     :title "Xotarium"
     :host "liquidfun"
     :setup #(setup (random-genome rng) seed)
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

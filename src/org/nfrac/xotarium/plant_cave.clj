(ns org.nfrac.xotarium.plant-cave
  "An experiment with generating plant-like connected forms."
  (:require [org.nfrac.liquidfun.testbed :as bed]
            [org.nfrac.liquidfun.core :as lf :refer [body! joint!
                                                     particle-system!]]
            [org.nfrac.liquidfun.vec2d :as v2d]
            [org.nfrac.xotarium.util.algo-graph :as graph]
            [org.nfrac.xotarium.util :as util]
            [org.nfrac.xotarium.proximity-field :as proxf]
            [quil.core :as quil :include-macros true]
            [quil.middleware]
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
            liquidfun$b2ParticleSystem
            liquidfun$b2ParticleGroup
            liquidfun$b2ParticleDef
            liquidfun$b2QueryCallback)))

(def p-radius 0.07)
(def cave-width 8.0)
(def cave-height 8.0)
(def cave-hw (* cave-width 0.5))
(def cave-hh (* cave-height 0.5))
(def pad 0.1)
(def min-group-size 5)
(def flow-force 0.005)
(def fix-radius (* p-radius 5))
(def flow-types #{:spin-ccw :spin-cw :down :up :out})

(def decay-every-n-steps 5)
(def decay-factor 0.97)

(def wall-points (atom #{}))

(defn record-wall-point! [x y]
  (swap! wall-points conj [(int (/ x fix-radius)) (int (/ y fix-radius))]))

(defn wall-point? [x y]
  (contains? @wall-points [(int (/ x fix-radius)) (int (/ y fix-radius))]))

(def accum-hits (atom {}))

(gen-class
 :name org.nfrac.xotarium.plant-cave.WallListener
 :extends org.bytedeco.javacpp.liquidfun$b2ContactListener
 :state state
 :init init
 :prefix "impl-")

;; getting IndexOutOfBoundsException if pass state arg to -init. workaround:
(let [state (atom nil)]

  (defn impl-init []
    [[] @state])

  (defn wall-listener
    [wall-pg wall-body]
    (reset! state [wall-pg wall-body])
    (org.nfrac.xotarium.plant-cave.WallListener.)))

(defn impl-BeginContact-b2ParticleSystem-b2ParticleBodyContact
  [^org.nfrac.xotarium.plant-cave.WallListener this
   ^liquidfun$b2ParticleSystem ps ^liquidfun$b2ParticleBodyContact contact]
  (let [i (.index contact)
        [wall-pg wall-body] (.state this)
        wall-pg ^liquidfun$b2ParticleGroup wall-pg
        wall-body ^liquidfun$b2Body wall-body]
    (when (and (not (.ContainsParticle wall-pg i))
               (= (.body contact) wall-body))
      (let [x (.GetParticlePositionX ps i)
            y (.GetParticlePositionY ps i)]
        (swap! accum-hits assoc i [x y])
        (record-wall-point! x y))))
  true)

(defn impl-BeginContact-b2ParticleSystem-b2ParticleContact
  [^org.nfrac.xotarium.plant-cave.WallListener this
   ^liquidfun$b2ParticleSystem ps ^liquidfun$b2ParticleContact contact]
  (let [ia (.GetIndexA contact)
        ib (.GetIndexB contact)
        [wall-pg _] (.state this)
        wall-pg ^liquidfun$b2ParticleGroup wall-pg
        a-wall? (.ContainsParticle wall-pg ia)
        b-wall? (.ContainsParticle wall-pg ib)
        [walli gasi] (cond
                       (and a-wall? b-wall?) nil
                       a-wall? [ia ib]
                       b-wall? [ib ia]
                       :else nil)]
    (when-let [i gasi]
      (let [x (.GetParticlePositionX ps i)
            y (.GetParticlePositionY ps i)]
        (swap! accum-hits assoc i [x y]))))
  true)

(defn sign [x] (if (neg? x) -1 1))

(defn setup
  []
  (let [world (lf/new-world)
        hw (* 0.5 cave-width)
        hh (* 0.5 cave-height)
        ground (body! world {:type :static}
                      {:shape (lf/box (+ hw pad) pad [0 (- 0 hh pad)])}
                      {:shape (lf/box (+ hw pad) pad [0 (+ hh pad)])}
                      {:shape (lf/box pad (+ hh pad) [(- 0 hw pad) 0])}
                      {:shape (lf/box pad (+ hh pad) [(+ hw pad) 0])})
        ps (particle-system! world
                             {:radius p-radius
                              :density 2.5
                              :pressure-strength 0.1
                              :elastic-strength 0.75
                              :damping-strength 0.5
                              :gravity-scale 0.0
                              :color-mixing-strength 0.5
                              ;:strict-contact-check true
                              :destroy-by-age false})]
    (assoc bed/initial-state
           :world world
           :ground ground
           :y-lower-bound (- hh)
           :y-upper-bound hh
           :particle-system ps
           :particle-iterations 3
           :dt-secs (/ 1 60.0)
           :camera (bed/map->Camera {:width cave-width
                                     :height cave-height
                                     :center [0 0]}))))

(defn add-random-static-bars
  [ground rng]
  (let [hw cave-hw
        hh cave-hh]
    (reduce (fn [rng [sx sy]]
              (let [[rng rng*] (random/split rng)
                    [r1 r2 r3 r4] (random/split-n rng* 4)
                    x0 (* (util/rand r1 (* 0.25 hw) hw) sx)
                    x1 (* (util/rand r2 (* 0.25 hw) hw) sx)
                    y0 (* (util/rand r3 (* 0.25 hh) hh) sy)
                    y1 (* (util/rand r4 (* 0.25 hh) hh) sy)
                    [ang mag] ((juxt v2d/v-angle v2d/v-mag)
                               (v2d/v-sub [x1 y1] [x0 y0]))]
                (lf/fixture! ground
                             {:shape (lf/rod [x0 y0] ang mag (* 4 p-radius))})
                rng))
            rng
            [[1 1] [1 -1] [-1 1] [-1 -1]])))

(defn set-up-accretion
  [state rng]
  (let [world ^liquidfun$b2World (:world state)
        ps ^liquidfun$b2ParticleSystem (:particle-system state)
        ground (:ground state)
        hw (* 0.5 cave-width)
        hh (* 0.5 cave-height)
        gas-pg (lf/particle-group!
                ps
                {:flags (lf/particle-flags #{:water :fixture-contact-listener})
                 :stride (* 0.85 p-radius 2.0)
                 :color [255 128 128 128]
                 :shape (lf/circle (* hw 0.2))})
        wall-pg (lf/particle-group!
                 ps
                 {:flags (lf/particle-flags #{:wall :barrier :particle-contact-listener})
                  :group-flags (lf/particle-group-flags #{:can-be-empty})
                  :color [255 255 255 255]})
        old-listener (.m_contactListener (.GetContactManager world))
        lstnr (wall-listener wall-pg ground)
        pdef (lf/particle-def
              {:flags (lf/particle-flags #{:wall :barrier :particle-contact-listener})
               :color [255 255 255 255]
               :group wall-pg})
        its (.CalculateReasonableParticleIterations world (/ 1 60.0))
        randsrc (atom (random/split-n rng (+ 1 1 (* 4 5))))
        getrng! #(let [x (first @randsrc)]
                   (swap! randsrc rest)
                   x)
        flow-type (util/rand-nth (getrng!) (seq flow-types))]
    (println "reasonable particle iterations:" its)
    (.SetContactListener world lstnr)
    ;; clear out the center of gas as particles can end up fixed there
    (.DestroyParticlesInShape ps (lf/box (* p-radius 3) (* p-radius 3))
                              (doto (liquidfun$b2Transform.) (.SetIdentity)))
    ;; add some random shapes to accrete on
    (add-random-static-bars ground (getrng!))
    ;; add some more random shapes to disturb the gas flow
    (let [scaffolding (lf/body! world {:type :static})]
      (doseq [[sx sy] [[1 1] [1 -1] [-1 1] [-1 -1]]]
        (let [x0 (* (util/rand (getrng!) (* 0.01 hw) (* 0.9 hw)) sx)
              x1 (* (util/rand (getrng!) (* 0.01 hw) (* 0.9 hw)) sx)
              y0 (* (util/rand (getrng!) (* 0.01 hh) (* 0.9 hh)) sy)
              y1 (* (util/rand (getrng!) (* 0.01 hh) (* 0.9 hh)) sy)]
          (lf/fixture! scaffolding
            {:shape (if (and (> (v2d/v-mag [x0 y0]) (* 0.35 hw))
                             (> (util/rand (getrng!) 0 0.3)))
                      (lf/circle 0.25 [x0 y0])
                      (lf/edge [x0 y0] [x1 y1]))})))
      (assoc state
       :contact-listener lstnr
       ::has-gas? true
       ::old-listener old-listener
       ::scaffolding scaffolding
       :settings {:wall-pg wall-pg
                  :gas-pg gas-pg
                  :pdef pdef
                  :flow-type flow-type}))))

(defn set-flow-force
  [flow-type ^double xz ^double yz ^liquidfun$b2Vec2 v2tmp]
  (case flow-type
     :spin-ccw
     (.Set v2tmp (* flow-force -2.0 yz)
                 (* flow-force 2.0 xz))
     :spin-cw
     (.Set v2tmp (* flow-force 2.0 yz)
                 (* flow-force -2.0 xz))
     :down
     (.Set v2tmp 0.0 (- flow-force))
     :up
     (.Set v2tmp 0.0 flow-force)
     :out
     (.Set v2tmp (* flow-force (sign xz))
                 (* flow-force (sign yz)))))

(defn end-gas
  [state]
  (let [gas-pg ^liquidfun$b2ParticleGroup (:gas-pg (:settings state))]
    (.DestroyParticles gas-pg)
    (lf/destroy-body! (::scaffolding state))
    ;; turn off accretion listener
    (let [lstnr (::old-listener state)]
      (.SetContactListener (:world state) lstnr)
      (assoc state ;:contact-listener lstnr
             ::has-gas? false))))

(defn gas-work
  [state]
  (let [dt (:dt-secs state)
        ps ^liquidfun$b2ParticleSystem (:particle-system state)
        {:keys [wall-pg gas-pg flow-type] :as settings} (:settings state)
        gas-pg ^liquidfun$b2ParticleGroup gas-pg
        pdef ^liquidfun$b2ParticleDef (:pdef settings)]
    (let [v2tmp (lf/vec2 0 0)
          gas-i (.GetBufferIndex gas-pg)]
      (dotimes [j (.GetParticleCount gas-pg)]
        (let [i (+ gas-i j)
              x (.GetParticlePositionX ps i)
              y (.GetParticlePositionY ps i)
              xz (/ x cave-width)
              yz (/ y cave-height)]
          (set-flow-force flow-type xz yz v2tmp)
          (.ParticleApplyForce ps i v2tmp))))
    (let [gas-n (.GetParticleCount gas-pg)
          hit-n (count @accum-hits)]
      (doseq [[i [x y]] @accum-hits]
        (.DestroyParticle ps i)
        (.SetPosition pdef x y)
        (lf/particle! ps pdef))
      (swap! accum-hits empty)
      (if (or (>= hit-n gas-n)
              ;; sometimes gas particles get stuck
              (> (:time state) 20.0))
        (end-gas state)
        state))))

(defn step
  [state]
  (cond-> (bed/world-step state)
          (::has-gas? state) (gas-work)))

(defn draw
  [state]
  (bed/draw state true)
  (let []
    (quil/fill 255)
    (quil/text (str "Keys: (g) toggle gravity (d) damping (b) ball"
                    " (s) split (c) color (z) zap groups (e) squishify"
                    " (a) air")
               10 10)))

(defn rand-color []
  [0
   255
   0
   (+ 64 (rand-int (- 256 64)))])

(defn particle-indices
  [^liquidfun$b2ParticleGroup pg]
  (let [i0 (.GetBufferIndex pg)]
    (range i0 (+ i0 (.GetParticleCount pg)))))

(defn do-split-into-groups
  [state]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)
        wall-pg ^liquidfun$b2ParticleGroup (get-in state [:settings :wall-pg])]
    (.SetGroupFlags wall-pg 0)
    (.SplitParticleGroup ps wall-pg)
    state))

(defn do-colorize
  [state]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)]
    (doseq [group (lf/particle-group-seq ps)]
      (let [colb (.GetColorBuffer ps)
            [r g b a] (rand-color)]
        (doseq [i (particle-indices group)]
          (let [coli (.position colb (long i))]
            (.Set coli r g b a)))))
    state))

(defn do-zap-small-groups
  [state]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)]
    (doseq [g (lf/particle-group-seq ps)]
      (when (< (.GetParticleCount g) min-group-size)
        (doseq [i (particle-indices g)]
          (lf/destroy-particle! ps i))))
    state))

(defn do-squishify
  [state]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)
        wall-flagval (lf/particle-flags #{:elastic :wall})
        ;; first pass, store all the group position data
        group-coords (for [g (lf/particle-group-seq ps)]
                       (for [i (particle-indices g)]
                         [(.GetParticlePositionX ps i)
                          (.GetParticlePositionY ps i)]))
        group-coords (mapv doall group-coords)]
    ;; destroy all particles (empty groups will be cleaned up next step)
    (doseq [g (lf/particle-group-seq ps)]
      (.DestroyParticles g false))
    ;; extend the range within which triads can form
    (.SetRadius ps (* 1.5 p-radius))
    ;; next pass, recreate all groups as squishy
    (doseq [coords group-coords]
      (let [pdata (lf/seq->v2arr coords)
            pg (lf/particle-group!
                ps
                {:flags (lf/particle-flags #{:elastic})
                 :group-flags (lf/particle-group-flags #{:solid})
                 :particle-count (count coords)
                 :position-data pdata
                 :color (rand-color)})]
        (doseq [i (particle-indices pg)]
          (let [x (.GetParticlePositionX ps i)
                y (.GetParticlePositionY ps i)]
            (when (wall-point? x y)
              (.SetParticleFlags ps i wall-flagval))))))
    ;; restore
    (.SetRadius ps p-radius)
    state))

(defn do-zap-detached
  [state]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)
        n (.GetParticleCount ps)
        nt (.GetTriadCount ps)
        tt (.GetTriads ps)]
    (loop [ed (zipmap (range n) (repeat #{}))
           j 0]
      (if (< j nt)
        (let [tt (.position tt j)
              ia (.indexA tt)
              ib (.indexB tt)
              ic (.indexC tt)]
          (recur (-> ed
                     (update ia into [ib ic])
                     (update ib into [ia ic])
                     (update ic into [ia ib]))
            (inc j)))
        ;; done
        (let [g (graph/directed-graph (keys ed) ed)
              cc (graph/scc g)]
          (doseq [is cc]
            (when (not-any? (fn [i]
                              (let [flag (.GetParticleFlags ps i)]
                                (lf/particle-flag? flag :wall)))
                            is)
              (doseq [i is]
                (lf/destroy-particle! ps i)))))))
    state))

(defn abs [x] (if (neg? x) (- x) x))

(defn do-add-air
  [state]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)
        air-pg (lf/particle-group!
                   ps
                   {:flags (lf/particle-flags #{:water :color-mixing})
                    :stride (* 0.75 p-radius 2.0)
                    :color [0 0 0 0]
                    :shape (lf/box cave-hw cave-hh)})]
    (assoc state ::air-pg air-pg)))

(defn get-prox-field
  [world]
  (proxf/proximity-field world
                         [(- 0 cave-hw pad) (+ cave-hw pad)]
                         [(- 0 cave-hh) (+ cave-hh pad)]
                         (* 2 p-radius) 2.0))

(defn do-colorize-air-by-proximity
  [state]
  ;; set alpha by proximity to walls. this is just for viz!
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)
        air-pg (::air-pg state)
        colb (.GetColorBuffer ps)
        proxf (get-prox-field (:world state))]
    (doseq [i (particle-indices air-pg)]
      (let [x (.GetParticlePositionX ps i)
            y (.GetParticlePositionY ps i)
            p (proxf/proximity-score proxf x y)
            coli (.position colb (long i))]
        (.Set coli 128 128 (int (* 255 p)) 255)))
    state))

(defn build-world
  []
  (->> (setup)
       (set-up-accretion)
       (iterate step)
       (drop-while ::has-gas?)
       (first)
       (do-split-into-groups)
       (step)
       (do-colorize)
       (do-zap-small-groups)
       (step)
       (do-squishify)
       (step)
       (do-zap-detached)
       (step)
       (do-add-air)
       (step)
       :world))

(defn groups-restore-color
  [^liquidfun$b2ParticleSystem ps pgroups [r g b a]]
  (let [pis (mapcat particle-indices pgroups)
        cbuf (.GetColorBuffer ps)
        r (int r)
        g (int g)
        b (int b)
        a (int a)]
    (doseq [i pis]
      (let [cbuf (.position cbuf (int i))]
        (doto cbuf
          (.r r)
          (.g g)
          (.b b)
          (.a a))))))

(defn group-decay-color
  [^liquidfun$b2ParticleSystem ps ^liquidfun$b2ParticleGroup pg ^double factor]
  (let [cbuf (.GetColorBuffer ps)
        start (.GetBufferIndex pg)
        end (+ start (.GetParticleCount pg))]
    (loop [i start]
      (when (< i end)
        (let [c (.position cbuf i)]
          (.multiplyPut c factor)
          (recur (inc i)))))))

(defn my-key-press
  [state event]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)]
    (case (:key event)
      :s (do-split-into-groups state)
      :c (do-colorize state)
      :z (do-zap-small-groups state)
      :e (do-squishify state)
      :x (do-zap-detached state)
      :a (-> (do-add-air state) (do-colorize-air-by-proximity))
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
  [seed & args]
  (let [rng (if seed
              (random/make-random (Long. (str seed)))
              (random/make-random))]
    (quil/sketch
     :title "Plant Cave"
     :host "liquidfun"
     :setup #(-> (setup) (set-up-accretion rng))
     :update (fn [s] (if (:paused? s) s (step s)))
     :draw draw
     :key-typed my-key-press
     :mouse-pressed bed/mouse-pressed
     :mouse-released bed/mouse-released
     :mouse-dragged bed/mouse-dragged
     :mouse-wheel bed/mouse-wheel
     :size [600 500]
     :features [:resizable]
     :middleware [quil.middleware/fun-mode])))

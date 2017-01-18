(ns org.nfrac.xotarium.accretion
  "An experiment with accreting wall particles."
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

(def p-radius 0.025)
(def cave-width 6.0)
(def cave-height 5.0)
(def cave-hw (* cave-width 0.5))

(def every-n 20)
(def min-contacts 3)
(def accrete-n 500)
(def suppression-radius 0.8)
(def min-angular-gap 3.0)
(def circ-force 0.005)

(def accum-hits (atom {}))

(defn wall-listener
  [^liquidfun$b2ParticleGroup wall-group]
  (proxy [liquidfun$b2ContactListener] []
    (BeginContact [^liquidfun$b2ParticleSystem ps
                   ^liquidfun$b2ParticleContact contact]
      (let [ia (.GetIndexA contact)
            ib (.GetIndexB contact)
            a-wall? (.ContainsParticle wall-group ia)
            b-wall? (.ContainsParticle wall-group ib)
            [walli gasi] (cond
                           (and a-wall? b-wall?) nil
                           a-wall? [ia ib]
                           b-wall? [ib ia]
                           :else nil)
            normal (cond->> (lf/v2xy (.GetNormal contact))
                            b-wall? (mapv -))]
        (when walli
          (swap! accum-hits update walli conj normal))))))

(defn setup []
  (let [world (lf/new-world)
        hw (* 0.5 cave-width)
        ground (body! world {:type :static}
                      {:shape (lf/edge-loop [[(- 0 hw 1) -1]
                                             [(- 0 hw 1) (+ cave-height 1)]
                                             [(+ hw 1) (+ cave-height 1)]
                                             [(+ hw 1) -1]])})
        ps (particle-system! world
                             {:radius p-radius
                              :pressure-strength 1.0
                              :damping-strength 0.0
                              :gravity-scale 0.0
                              :destroy-by-age false})
        gas-pg (lf/particle-group!
                ps
                {:flags #{:water}
                 :stride (* 1.5 p-radius 2.0)
                 :color [255 128 128 128]
                 :shape (lf/box (* hw 0.6) (* cave-height 0.5 0.6)
                                [0 (/ cave-height 2)])})
        wall-pg (lf/particle-group!
                 ps
                 {:flags #{:wall :barrier :particle-contact-listener}
                  :color [255 255 255 255]
                  :shape (lf/edge-loop [[(- hw) 0]
                                        [(- hw) cave-height]
                                        [hw cave-height]
                                        [hw 0]])})
        lstnr (wall-listener wall-pg)
        pdef (lf/particle-def {:flags #{:wall :barrier :particle-contact-listener}
                               :color [255 255 255 255]
                               :group wall-pg})
        its (.CalculateReasonableParticleIterations world (/ 1 60.0))]
    (println "reasonable particle iterations:" its)
    (.SetContactListener world lstnr)
    ;; add some wall seeds in the middle. coords in range [-1 1]
    (doseq [[xz yz] (list [-0.10 0.85] [0.10 0.85]
                          [-0.80 0.50] [0.80 0.60]
                          [0.00 0.05]
                          [-0.70 0.00] [0.70 0.00]
                          [-0.20 -0.50] ;[0.20 -0.50]
                          [-0.70 -0.75] [0.70 -0.75])]
      (.SetPosition pdef (* cave-hw xz) (* cave-height (/ (+ yz 1) 2)))
      (lf/particle! ps pdef))
    (assoc bed/initial-state
      :world world
      :contact-listener lstnr
      :settings {:wall-pg wall-pg
                 :gas-pg gas-pg
                 :pdef pdef
                 :circ? true}
      :particle-system ps
      :particle-iterations its
      :dt-secs (/ 1 60.0)
      :camera (bed/map->Camera {:width 6 :height 5 :center [0 2.5]}))))

(defn mean [xs] (/ (double (apply + xs)) (count xs)))

(defn gap-angle
  [angles]
  (if (seq angles)
    (let [sa* (sort angles)
          sa (concat [(- (last sa*) (* 2 Math/PI))]
                     sa*
                     [(+ (first sa*) (* 2 Math/PI))])]
      (loop [sa (rest sa)
             last-a (first sa)
             big-gap 0.0
             big-a 0.0]
        (if-let [a (first sa)]
          (let [gap (- a last-a)]
            (if (or (> gap big-gap)
                    (and (== gap big-gap)
                         (> (rand) 0.5)))
              (recur (rest sa) a (double gap) (double (/ (+ a last-a) 2.0)))
              (recur (rest sa) a big-gap big-a)))
          ;; done
          [big-a big-gap])))
    ;; no angles; choose a random angle
    [(* (rand) 2 Math/PI) (* 2 Math/PI)]))

(defn post-step
  [state]
  (let [dt (:dt-secs state)
        ps ^liquidfun$b2ParticleSystem (:particle-system state)
        {:keys [wall-pg gas-pg circ?] :as settings} (:settings state)
        wall-pg ^liquidfun$b2ParticleGroup wall-pg
        gas-pg ^liquidfun$b2ParticleGroup gas-pg
        pdef ^liquidfun$b2ParticleDef (:pdef settings)
        v2tmp (lf/vec2 0 0)
        gas-i (.GetBufferIndex gas-pg)]
    (when circ?
      (dotimes [j (.GetParticleCount gas-pg)]
        (let [i (+ gas-i j)
              x (.GetParticlePositionX ps i)
              y (.GetParticlePositionY ps i)
              xz (/ x cave-width)
              yz (- (/ y cave-height) 0.5)]
          (.Set v2tmp
                (* circ-force -2.0 yz)
                (* circ-force 2.0 xz))
          (.ParticleApplyForce ps i v2tmp))))
    (when (and (zero? (mod (int (/ (:time state) dt)) every-n))
               (seq @accum-hits))
      (loop [hits (->> @accum-hits
                       (sort-by (fn [[k v]] (count v)) >)
                       (take accrete-n)
                       (take-while (fn [[k v]] (>= (count v) min-contacts))))]
        (when (seq hits)
          (let [[i normals] (first hits)
                ox (.GetParticlePositionX ps i)
                oy (.GetParticlePositionY ps i)
                d suppression-radius
                aabb (lf/aabb [(- ox d) (- oy d)] [(+ ox d) (+ oy d)])
                bb-angles (atom ())
                to-suppress (atom ())
                cb (proxy [liquidfun$b2QueryCallback] []
                     (ReportParticle [^liquidfun$b2ParticleSystem ps ii]
                       (when (and (not= i ii)
                                  (.ContainsParticle wall-pg ii))
                         (let [x (.GetParticlePositionX ps ii)
                               y (.GetParticlePositionY ps ii)
                               angle (v2d/v-angle [(- x ox) (- y oy)])]
                           (swap! bb-angles conj angle)
                           ;; suppress growth in nearby particles
                           (swap! to-suppress conj ii)))
                       true))
                _ (.QueryAABB ps cb aabb)
                [angle angular-gap] (gap-angle @bb-angles)]
            (swap! accum-hits dissoc i)
            (when (>= angular-gap min-angular-gap)
              (let [[dx dy] (v2d/polar-xy (* p-radius 2) angle)
                    nx (+ ox dx)
                    ny (+ oy dy)]
                (when (and (<= (- cave-hw) nx cave-hw)
                           (<= 0.0 ny cave-height))
                  (when (seq @to-suppress)
                    (swap! accum-hits #(apply dissoc % @to-suppress)))
                  (.SetPosition pdef nx ny)
                  (lf/particle! ps pdef))))
            (recur (->> (rest hits)
                        ;; remove suppressed
                        (filter (fn [[i _]] (contains? @accum-hits i)))))))))
      ;(swap! accum-hits #(into {} (remove (fn [[k v]] (== 1 (count v)))) %)))
    state))

(defn step
  [state]
  (-> (bed/world-step state)
      (post-step)))
      ;(bed/record-snapshot true)))

(defn draw
  [state]
  (bed/draw state true)
  (let []
    (quil/fill 255)
    (quil/text (str "Keys: (g) toggle gravity (d) damping (c) circulation\n"
                    "      (p) print wall coords")
               10 10)))

(defn my-key-press
  [state event]
  (let [pdef ^liquidfun$b2ParticleDef (get-in state [:settings :pdef])
        ps ^liquidfun$b2ParticleSystem (:particle-system state)
        wall-pg ^liquidfun$b2ParticleGroup (get-in state [:settings :wall-pg])]
    (case (:key event)
      :g (do
           (.SetGravityScale ps (if (== 1.0 (.GetGravityScale ps))
                                  0.0 1.0))
           state)
      :d (do
           (.SetDamping ps (if (== 1.0 (.GetDamping ps))
                             0.0 1.0))
           state)
      :c (update-in state [:settings :circ?] #(not %))
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
   :title "Accretion"
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

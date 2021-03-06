(ns org.nfrac.xotarium.grn-creature
  (:require [org.nfrac.xotarium.creature :as cre]
            [org.nfrac.xotarium.grn.greans :as grn]
            [org.nfrac.xotarium.plant-cave :as cave]
            [org.nfrac.xotarium.proximity-field :as proxf]
            [org.nfrac.cppn :as cppn]
            [org.nfrac.liquidfun.testbed :as bed]
            [org.nfrac.liquidfun.core :as lf :refer [body! joint!
                                                     particle-system!]]
            [org.nfrac.liquidfun.vec2d :as v2d]
            [org.nfrac.xotarium.util :as util]
            [quil.core :as quil :include-macros true]
            [quil.middleware]
            [clojure.pprint]
            [clojure.spec.alpha :as s]
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

(def max-sensed-velocity 8.0)
(def grn-steps-per-physics-step 2)
(def init-osc-exp-coord-dist 0.3)

;; this list should match the args passed to grn/step!
(def beh-inputs '[bias
                  oscillator
                  factor-a
                  factor-b
                  factor-c
                  message-a
                  message-b
                  message-c
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

(def beh-outputs '[expansion
                   messenger-a
                   messenger-b
                   messenger-c
                   ])

(defn index-of
  [x coll]
  (some (fn [[i s]] (when (= s x) i))
        (map-indexed vector coll)))

(s/def ::structure-cppn ::cppn/cppn)
(s/def ::signals-cppn ::cppn/cppn)

(s/def ::genome
  (s/keys :req-un [::structure-cppn
                   ::signals-cppn
                   ::grn/grn]))

(defn force-oscillation
  "Move promoter of expansion close to oscillator."
  [grn rng]
  (let [[r1 r2] (random/split rng)
        cg (s/conform ::grn/grn grn)
        units (:units (::grn/elements cg))
        unit-tfs (grn/index-counts (map #(count (:genes %)) units) 0)
        tf->unit-id (into {}
                          (mapcat (fn [i tfs]
                                    (map vector tfs (repeat i)))
                                  (range)
                                  unit-tfs))
        tf->gene-idx (into {}
                           (mapcat (fn [tfs]
                                     (map vector tfs (range)))
                                   unit-tfs))
        osc-tf (nth (::grn/input-tfs grn) (index-of 'oscillator beh-inputs))
        exp-tf (nth (::grn/output-tfs grn) (index-of 'expansion beh-outputs))
        osc-coords (get-in cg [::grn/elements :units (tf->unit-id osc-tf)
                               :genes (tf->gene-idx osc-tf) ::grn/coords])
        z init-osc-exp-coord-dist
        updat-cg (assoc-in cg [::grn/elements :units (tf->unit-id exp-tf)
                               :promoters 0 ::grn/coords]
                           (v2d/v-add osc-coords [(util/rand r1 (- z) z)
                                                  (util/rand r2 (- z) z)]))]
    (s/unform ::grn/grn updat-cg)))

(defn random-genome
  [rng grn-parameters]
  (let [structure-cppn cre/seed-structure-cppn
        signals-cppn cre/seed-signals-cppn
        [r1 r2] (random/split-n rng 2)
        grn* (grn/random-grn (count beh-inputs) (count beh-outputs) r1
                             grn-parameters)
        grn (force-oscillation grn* r2)]
    {:structure-cppn structure-cppn
     :signals-cppn signals-cppn
     :grn grn}))

(defn mutate
  [genome rng mutate-cppn-prob]
  (let [{:keys [grn signals-cppn]} genome
        cppn signals-cppn
        [r1 r2 r3] (random/split-n rng 3)
        ncppn (if (< (random/rand-double r1) mutate-cppn-prob)
                (cppn/mutate-with-perturbation cppn r2 {})
                cppn)
        [ngrn mut-info] (grn/mutate grn r3)]
    [(assoc genome
            :signals-cppn ncppn
            :grn ngrn)
     mut-info]))

(defn crossover
  [g1 g2 rng]
  (assoc g1 :grn (grn/crossover (:grn g1) (:grn g2) rng)))

(defn- update!
  [t-m k f & args]
  (assoc! t-m k (apply f (get t-m k) args)))

(defn triad-neighbours
  [tris]
  (let [by-handle (persistent!
                   (reduce (fn [m [h1 h2 h3 :as tri]]
                             (-> m
                                 (update! h1 conj tri)
                                 (update! h2 conj tri)
                                 (update! h3 conj tri)))
                           (transient {})
                           tris))]
    (persistent!
     (reduce (fn [m [h1 h2 h3 :as tri]]
               (assoc! m tri (->> (concat
                                   (get by-handle h1)
                                   (get by-handle h2)
                                   (get by-handle h3))
                                  (remove #(= tri %))
                                  (distinct))))
             (transient {})
             tris))))

(defn come-alive
  [world genome]
  (when-let [creature-body (cre/make-creature world (:structure-cppn genome)
                                              (:signals-cppn genome))]
    (let [grn (:grn genome)
          tri-p (:triad-params creature-body)
          tri-nbs (triad-neighbours (keys tri-p))
          cell-form (grn/grn->cell grn)
          tri-concs (zipmap (keys tri-p)
                            (repeat (::grn/concs cell-form)))]
      (assoc creature-body
            :tri-concs tri-concs
            :tri-nbs tri-nbs
            :grn-cell cell-form
            :grn grn))))

(defn setup
  [genome seed]
  (let [;world (cave/build-world)
        state (-> (cave/setup)
                  (cave/do-add-air)
                  (cave/do-add-color-sources))
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
        tri-nbs (:tri-nbs creature)
        n-tfs (count (first (vals tri-concs)))
        ;; messengers must start at output index 1:
        [msgr-a msgr-b msgr-c] (->> (drop 1 (::grn/output-tfs cell-form))
                                    (map #(mod % n-tfs)))
        dt (:dt (::grn/parameters (:grn creature)) 0.2)
        time (:time state)
        freq 8.0
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
                             (assoc! ia (record-contact (get pcm ia) kb
                                                        (.position colbuf ib)))
                             (= :self-muscle kb)
                             (assoc! ib (record-contact (get pcm ib) ka
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
                         nbs (get tri-nbs handles)
                         [msg-a msg-b msg-c] (reduce (fn [[ma mb mc] nb]
                                                       (let [cs (get tri-concs nb)]
                                                         [(+ ma (cs msgr-a))
                                                          (+ mb (cs msgr-b))
                                                          (+ mc (cs msgr-c))]))
                                                     [0.0 0.0 0.0]
                                                     nbs)
                         cell (assoc cell-form ::grn/concs concs)
                         phase-off (:phase-off params)
                         osc (if (pos? (Math/sin (+ phase phase-off)))
                               1.0 0.0)
                         ic [1.0 ; bias
                             osc
                             (:factor-a params)
                             (:factor-b params)
                             (:factor-c params)
                             (min msg-a 1.0)
                             (min msg-b 1.0)
                             (min msg-c 1.0)
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
                         ncell (->> cell
                                    (iterate #(grn/step % ic dt))
                                    (rest)
                                    (take grn-steps-per-physics-step)
                                    (last))]
                     (assoc! m handles (::grn/concs ncell))))
                 (transient {})
                 tri-concs))
        work-fn (fn [handles params]
                  (let [concs (get tri-concs handles)
                        cell (assoc cell-form ::grn/concs concs)
                        g-out (grn/cell-outputs cell)
                        ;; expansion must be first output
                        expan (-> (first g-out) (* 2.0) (- 1.0))]
                    expan))]
    (cre/creature-flex ps (:triad-params creature) work-fn)
    (doseq [grp (mapcat val (:groups creature))]
      (cave/group-apply-gravity grp cre/p-gravity))
    (cave/groups-restore-color ps (:muscle (:groups creature)) [255 0 0 255])
    (cave/groups-restore-color ps [(::cave/green-source state)] [0 255 0 255])
    (cave/groups-restore-color ps [(::cave/blue-source state)] [0 0 255 255])
    (when (zero? (mod (int (/ time (:dt-secs state))) cave/decay-every-n-steps))
      (cave/group-decay-color ps (::cave/air-pg state) cave/decay-factor))
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
     :setup #(setup (random-genome rng grn/parameter-defaults) seed)
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

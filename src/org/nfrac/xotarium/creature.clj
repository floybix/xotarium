(ns org.nfrac.xotarium.creature
  (:require [org.nfrac.xotarium.cppn :as cppn]
            [org.nfrac.xotarium.cppn-compile :as cc]
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

(def creature-width 1.8)
(def creature-height 1.8)
(def min-creature-particles 10)
(def muscle-strength 0.9)

;; interpretation:
;; if muscle is positive, express muscle.
;;   - if bone is also positive, the muscle binds to adjacent bone (reactive).
;; if bone is positive (muscle is not), express bone.
;; otherwise (both non-positive), empty space.

(def seed-cppn
  {:inputs #{:bias :x :y :d}
   :outputs #{:bone :muscle :angle :phase-off :factor-a :factor-b :factor-c}
   :finals #{:bone :muscle}
   :zerod #{}
   :nodes {:i0 :gaussian}
   :edges {:i0 {:d -1.0}
           :factor-a {:x 1.0}
           :factor-b {:y 1.0}
           :factor-c {:d 0.5}
           :phase-off {:x 2.0}
           :angle {:factor-a 1.0}
           :bone {:i0 1.0
                  :bias -0.3}
           :muscle {:factor-c -1.0
                    :factor-a 1.0}}})

(defn hex-neighbours
  "Returns hexagonal coordinates linked to their immediate neighbours.
  Ids are like [ix iy] which are integer or fractional (1/2) steps.

  Uses a hex grid: every odd row is offset by 1/2 in x.
  Below, the x's neighbours are marked *.
```
  o * * o
   * x *
  o * * o
   o o o
```"
  [nx ny]
  (into
   {}
   (for [iy (range 0 ny)
         ;; every second row offset by 1/2
         ix (range (if (even? iy) 0 1/2) nx)]
     [[ix iy] (filterv (fn [[jx jy]]
                         (and (<= 0 jx) (< jx nx)
                              (<= 0 jy) (< jy ny)))
                       [[(dec ix) iy]
                        [(inc ix) iy]
                        [(- ix 1/2) (dec iy)]
                        [(+ ix 1/2) (dec iy)]
                        [(- ix 1/2) (inc iy)]
                        [(+ ix 1/2) (inc iy)]])])))

(defn filter-graph
  "Filters the graph to nodes (expressed? id) and edges (connected? id1 id2)."
  [g-nbs expressed? connected?]
  (loop [ids (keys g-nbs)
         nbs {}]
    (if-let [id (first ids)]
      (if (expressed? id)
        (let [i-nbs (filter #(connected? id %) (get g-nbs id))]
          (recur (rest ids) (assoc nbs id i-nbs)))
        ;; not expressed
        (recur (rest ids) nbs))
      ;; done
      nbs)))

(defn tissue-groups
  "Returns groups of particle coordinates, partitioned into each contiguous
  bone, and the remaining muscle partitioned between reactive and inert."
  [hex-g pheno]
  (let [expressed? (fn [id]
                     (let [m (pheno id)]
                       (or (pos? (:bone m)) (pos? (:muscle m)))))
        muscle? (fn [id]
                  (pos? (:muscle (pheno id))))
        bone? (fn [id]
                (and (expressed? id) (not (muscle? id))))
        bone-g (filter-graph hex-g bone? #(and (bone? %1) (bone? %2)))
        bones (graph/scc (graph/directed-graph (keys bone-g) bone-g))
        muscle-ids (filter muscle? (keys pheno))
        muscle-inert (set (remove #(pos? (:bone (pheno %))) muscle-ids))
        muscle-reactive (set (remove muscle-inert muscle-ids))]
    (concat
     (map vector (repeat :bone) (repeat :inert) bones)
     [[:muscle :inert muscle-inert]
      [:muscle :reactive muscle-reactive]])))

(defn morpho-data
  "Calculates the expressed tissues and returns the blueprint for the body.
  Returns a sequence of groups, each containing a set of coordinates of one
  tissue type. Contiguous bones are returned as separate groups.
  Tissue is one of :bone or :muscle.
  Muscle particles have associated parameters."
  [cppn-fn p-radius [cx cy] [w h]]
  (let [stride (* p-radius 2.0 0.75)
        nx (quot w stride)
        ny (quot h stride)
        hex-g (hex-neighbours nx ny)
        x0 (- cx (/ w 2))
        y0 (- cy (/ h 2))
        coord (fn [[ix iy]] [(+ x0 (* ix stride))
                             (+ y0 (* iy stride))])
        p-round (fn [[x y]] [(quot x p-radius)
                             (quot y p-radius)])
        pheno (into
               {}
               (map (fn [id]
                      (let [[ix iy] id
                            zx (-> (+ ix 0.5) (/ nx) (- 0.5) (* 2))
                            zy (-> (+ iy 0.5) (/ ny) (- 0.5) (* 2))
                            d (cc/xy->d zx zy)
                            bias 1.0]
                        [id (cppn-fn bias d zx zy)])))
               (keys hex-g))
        tgroups (tissue-groups hex-g pheno)]
    (map (fn [[tissue reactivity ids]]
           {:tissue tissue
            :reactivity reactivity
            :coords (map coord ids)
            :pcoord-params (zipmap (map (comp p-round coord) ids)
                                   (map pheno ids))})
         tgroups)))


(s/def ::tissue #{:bone :muscle})
(s/def ::reactivity #{:reactive :inert})

(s/def ::coord (s/tuple double? double?))
(s/def ::coords (s/every ::coord))

(s/def ::tissue-group
  (s/keys :req-un [::tissue
                   ::reactivity
                   ::coords
                   ::pcoord-params]))

(s/fdef morpho-data
        :args (s/cat :cppn-fn ::cppn-fn
                     :p-radius double?
                     :center ::coord
                     :size (s/tuple double? double?))
        :ret (s/coll-of ::tissue-group))

(def tissue-color
  {[:bone :inert] [255 255 255 255]
   [:muscle :inert] [255 0 0 255]
   [:muscle :reactive] [255 192 192 255]})

(defn bounding-box
  [coords]
  (let [xs (map first coords)
        ys (map second coords)]
    [[(reduce min xs) (reduce min ys)]
     [(reduce max xs) (reduce max ys)]]))

(defn clear-push
  [^liquidfun$b2World world coords]
  (let [[bb-lo bb-hi] (bounding-box coords)
        span (v2d/v-mag (v2d/v-sub bb-hi bb-lo))
        bb-mid (v2d/v-scale (v2d/v-add bb-lo bb-hi) 0.5)
        aabb (lf/aabb bb-lo bb-hi)
        cb (proxy [liquidfun$b2QueryCallback] []
             (ReportFixture
              [fixt]
              (let [b (lf/body-of fixt)
                    [x y] (lf/position b)
                    d (v2d/v-sub [x y] bb-mid)
                    z (v2d/v-mag d)
                    dn (v2d/v-scale d (/ span (+ z 0.01)))]
                (lf/linear-velocity! b dn))
              true)
             (ShouldQueryParticleSystem
              [ps]
              true)
             (ReportParticle
              [^liquidfun$b2ParticleSystem ps i]
              (let [x (.GetParticlePositionX ps i)
                    y (.GetParticlePositionY ps i)
                    d (v2d/v-sub [x y] bb-mid)
                    z (v2d/v-mag d)
                    [dnx dny] (v2d/v-scale d (/ span (+ z 0.01)))]
                (.SetParticleVelocity ps i dnx dny))
              true))]
    (.QueryAABB world cb aabb)))

(defn morpho-construct-particle-groups
  "Returns a map of tissue to sequences of the constructed particle groups."
  [mdata ^liquidfun$b2ParticleSystem ps]
  (when false ;; scaffold to hold creature in place
    (lf/particle-group! ps
                        {:flags (lf/particle-flags #{:wall :reactive})
                         :color [0 255 0 255]
                         :shape (lf/edge [-1.0 0.0] [1.0 0.0])}))
  (->>
    (for [{:keys [tissue reactivity coords]} mdata]
      (let [pdata (lf/seq->v2arr coords)]
        [tissue
         (lf/particle-group!
          ps
          {:flags (lf/particle-flags
                   (set (concat
                         (case tissue
                           :muscle [:elastic :color-mixing]
                           :bone [])
                         (case reactivity
                           :reactive [:reactive]
                           :inert []))))
           :group-flags (lf/particle-group-flags
                         (case tissue
                           :bone #{:rigid :solid}
                           :muscle #{:solid}))
           :color (tissue-color [tissue reactivity])
           :position-data pdata
           :particle-count (count coords)})]))
    (reduce (fn [m [tissue g]]
              (update m tissue conj g))
            {})))

(defn morpho-particle-parameters
  "Returns a map of particle handle to parameters map,
  for the muscle tissues."
  [mdata ^liquidfun$b2ParticleSystem ps pgm]
  (let [p-radius (.GetRadius ps)
        p-round (fn [[x y]] [(quot x p-radius)
                             (quot y p-radius)])
        pcoord-params (->> mdata
                           (filter #(= :muscle (:tissue %)))
                           (map :pcoord-params)
                           (reduce merge))
        pis (mapcat cave/particle-indices (:muscle pgm))]
    (into {}
      (for [i pis]
        (let [h (.GetParticleHandleFromIndex ps i)
              x (.GetParticlePositionX ps i)
              y (.GetParticlePositionY ps i)
              params (pcoord-params (p-round [x y]))]
          [h params])))))

(defn mean-kv
  [ms]
  (let [n (count ms)]
    (->> (reduce (partial merge-with +) ms)
         (into {} (map (fn [[k v]] [k (/ (double v) n)]))))))

(defn all-triads-by-handles
  "Returns a map of [handle-a handle-b handle-c] to triads."
  [^liquidfun$b2ParticleSystem ps]
  (let [nt (.GetTriadCount ps)
        tt (.GetTriads ps)]
    (into {}
          (map (fn [j]
                 (let [tt (.position tt (long j))
                       ia (.indexA tt)
                       ib (.indexB tt)
                       ic (.indexC tt)
                       ha (.GetParticleHandleFromIndex ps ia)
                       hb (.GetParticleHandleFromIndex ps ib)
                       hc (.GetParticleHandleFromIndex ps ic)]
                   [[ha hb hc] (liquidfun$b2ParticleTriad. tt)])))
          (range nt))))

(defn morpho-triad-parameters
  "Returns a map of [handle-a handle-b handle-c] to the triad parameters, for
  the muscle tissues, with :ref-pa :ref-pb :ref-pc added to the parameters map
  with [x y] vals. So this is used to store the baseline positions once, then we
  will look up triads on every time step using handles (can not keep pointers to
  triads because they are re-allocated by GrowableBuffer)."
  [^liquidfun$b2ParticleSystem ps pp]
  (let [h-tri (all-triads-by-handles ps)]
    (into {}
          (comp
            (map (fn [[handles ^liquidfun$b2ParticleTriad tt]]
                   (let [abcp (remove empty? (map pp (seq handles)))]
                     (when (seq abcp) ;(= 3 (count abcp))
                       [handles
                        (assoc (mean-kv abcp)
                               :ref-pa (lf/v2xy (.pa tt))
                               :ref-pb (lf/v2xy (.pb tt))
                               :ref-pc (lf/v2xy (.pc tt)))]))))
            (remove nil?))
          h-tri)))

(defn v2-flex
  [^liquidfun$b2Vec2 v [ref-x ref-y] scale-x scale-y]
  (.x v (* ref-x scale-x))
  (.y v (* ref-y scale-y)))

(defn creature-flex
  "Behaviour is defined per triad by work-fn, taking [handles params]. Argument
  handles is a triplet of b2ParticleHandle, used to identify a triad, and params
  is a map of the cppn outputs (averaged over the triad particles). The work-fn
  should return triad expansion/contraction between -1.0 and 1.0."
  [^liquidfun$b2ParticleSystem ps tri-p work-fn]
  (let [h-tri (all-triads-by-handles ps)]
    (doseq [[handles ^liquidfun$b2ParticleTriad triad] h-tri
            :let [params (get tri-p handles)]
            :when params]
      (let [{:keys [ref-pa ref-pb ref-pc]} params
            {:keys [angle]} params
            strength-x (Math/cos angle)
            strength-y (Math/sin angle)
            mag (work-fn handles params)
            scale-x (+ 1.0 (* strength-x mag muscle-strength))
            scale-y (+ 1.0 (* strength-y mag muscle-strength))]
        (v2-flex (.pa triad) ref-pa scale-x scale-y)
        (v2-flex (.pb triad) ref-pb scale-x scale-y)
        (v2-flex (.pc triad) ref-pc scale-x scale-y)))))

(defn abs [x] (if (neg? x) (- x) x))

(defn morpho-cppn-fn
  [cppn]
  (let [f (cc/build-cppn-fn cppn)]
    (fn [& args]
      (-> (apply f args)
          (update :phase-off * Math/PI)
          (update :angle * Math/PI)
          (update :factor-a abs)
          (update :factor-b abs)
          (update :factor-c abs)))))

(defn make-creature
  [^liquidfun$b2World world cppn]
  (let [ps ^liquidfun$b2ParticleSystem (first (lf/particle-sys-seq world))
        p-radius (.GetRadius ps)
        cppn-fn (morpho-cppn-fn cppn)
        mdata (morpho-data cppn-fn p-radius [0 0] [creature-width creature-height])
        coords (->> mdata
                    (filter #(= :muscle (:tissue %)))
                    (mapcat :coords))
        valid-cre? (>= (count coords) min-creature-particles)]
    (if valid-cre?
      (do
        (dotimes [i 60]
          (clear-push world coords)
          (lf/step! world (/ 1 60.0) 8 3 3))
        (let [pgm (morpho-construct-particle-groups mdata ps)
              pp (morpho-particle-parameters mdata ps pgm)
              tri-p (morpho-triad-parameters ps pp)]
          {:cppn cppn
           :groups pgm
           :particle-params pp
           :triad-params tri-p}))
      ;; invalid creature
      nil)))

(defn destroy-creature
  [creature]
  (doseq [grps (vals (:groups creature))
          grp grps]
    (.DestroyParticles ^liquidfun$b2ParticleGroup grp)))

(defn setup
  []
  (let [;world (cave/build-world)
        state (-> (cave/setup)
                  (cave/do-add-air))
        cppn seed-cppn
        creature (make-creature (:world state) cppn)]
      (assoc state
             :creature creature)))

(defn post-step
  [state]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)
        creature (:creature state)
        time (:time state)
        work-fn (fn [_ params]
                  (let [freq 6.0
                        phase (mod (* time freq) (* 2.0 Math/PI))
                        phase-off (:phase-off params)]
                    (Math/sin (+ phase phase-off))))]
    (when creature
      (creature-flex ps (:triad-params creature) work-fn)
      (cave/groups-restore-color ps (:muscle (:groups creature)) [255 0 0 255]))
    (when (zero? (mod (int (/ time (:dt-secs state))) cave/decay-every-n-steps))
      (cave/group-decay-color ps (::cave/air-pg state) cave/decay-factor))
    state))

(defn step
  [state]
  (-> (bed/world-step state)
      (post-step)))

(defn draw
  [state]
  (bed/draw state true)
  (let []
    (quil/fill 255)
    (quil/text (str "Keys: (g) toggle gravity (d) damping (b) ball"
                    " (a) air (n) new mutant")
               10 10)))

(defn my-key-press
  [state event]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)]
    (case (:key event)
      :a (cave/do-add-air state)
      :b (do
           (body! (:world state) {}
                  {:shape (lf/circle 0.25)
                   :restitution 0.1
                   :density 1.0})
           (cave/add-random-static-bars (:ground state) (random/make-random 1))
           state)
      :g (do
           (.SetGravityScale ps (if (== 1.0 (.GetGravityScale ps))
                                  0.0 1.0))
           state)
      :d (do
           (.SetDamping ps (if (== 1.0 (.GetDamping ps))
                             0.0 1.0))
           state)
      :n (do
           (destroy-creature (:creature state))
           (assoc state :creature
                  (make-creature (:world state)
                                 (cppn/mutate-with-perturbation seed-cppn {}))))
      ;; otherwise pass on to testbed
      (bed/key-press state event))))

(defn run
  [& args]
  (quil/sketch
      :title "Xotarium"
      :host "liquidfun"
      :setup setup
      :update (fn [s] (if (:paused? s) s (step s)))
      :draw #(if (zero? (mod (quil/frame-count) 2)) (draw %) %)
      :key-typed my-key-press
      :mouse-pressed bed/mouse-pressed
      :mouse-released bed/mouse-released
      :mouse-dragged bed/mouse-dragged
      :mouse-wheel bed/mouse-wheel
      :size [600 500]
      :features [:resizable]
      :middleware [quil.middleware/fun-mode]))

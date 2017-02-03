(ns org.nfrac.xotarium.run
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
            [clojure.spec :as s])
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

;; interpretation:
;; if muscle is positive, express muscle.
;;   - if bone is also positive, the muscle binds to adjacent bone.
;; if bone is positive (muscle is not), express bone.
;; otherwise (both non-positive), empty space.

(def seed-cppn
  {:inputs #{:bias :x :y :d}
   :outputs #{:bone :muscle :phase}
   :finals #{:bone :phase}
   :nodes {:i0 :gaussian}
   :edges {:i0 {:d -1.0}
           :muscle {:i0 1.0
                    :bias 0.2}
           :bone {:d -1.0
                  :bias -0.1}
           :phase {;:muscle 1.0
                   :x 1.0}}})

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
  bone, and the remaining muscle partitioned between binding and inert."
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
        muscle-binding (set (remove muscle-inert muscle-ids))]
    (concat
     (map vector (repeat :bone) bones)
     [[:muscle-inert muscle-inert]
      [:muscle-binding muscle-binding]])))

(defn morpho-data
  "Calculates the expressed tissues and returns the blueprint for the body.
  Returns a sequence of groups, each containing a set of coordinates of one
  tissue type. Contiguous bones are returned as separate groups.
  Tissue is one of :bone, :muscle-binding, :muscle-inert.
  Muscle particles have an associated phase attribute."
  [cppn-fn p-radius [cx cy] [w h]]
  (let [stride (* p-radius 2.0 0.75)
        nx (quot w stride)
        ny (quot h stride)
        hex-g (hex-neighbours nx ny)
        x0 (- cx (/ w 2))
        y0 (- cy (/ h 2))
        coord (fn [[ix iy]] [(+ x0 (* ix stride))
                             (+ y0 (* iy stride))])
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
    (map (fn [[tissue ids]]
           {:tissue tissue
            :coords (map coord ids)
            :phases (map #(:phase (pheno %)) ids)})
         tgroups)))


(s/def ::tissue #{:bone :muscle-binding :muscle-inert})

(s/def ::coord (s/tuple double? double?))
(s/def ::coords (s/every ::coord))

(s/def ::tissue-group
  (s/keys :req-un [::tissue
                   ::coords
                   ::phases]))

(s/fdef morpho-data
        :args (s/cat :cppn-fn ::cppn-fn
                     :p-radius double?
                     :center ::coord
                     :size (s/tuple double? double?))
        :ret (s/coll-of ::tissue-group))

(def tissue-color
  {:bone [255 255 255 255]
   :muscle-binding [255 128 128 255]
   :muscle-inert [255 0 0 255]})

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

(defn morpho-construct
  [mdata ^liquidfun$b2World world ^liquidfun$b2ParticleSystem ps]
  (doall
    (for [{:keys [tissue coords phases]} mdata]
      (let [pdata (lf/seq->v2arr coords)]
        (lf/particle-group!
         ps
         {:flags (lf/particle-flags
                  (case tissue
                   :bone #{}
                   :muscle-inert #{:elastic}
                   :muscle-binding #{:elastic :reactive}))
          :group-flags (lf/particle-group-flags
                        (case tissue
                          :bone #{:rigid :solid}
                          (:muscle-inert :muscle-binding) #{:solid}))
          :color (tissue-color tissue)
          :position-data pdata
          :particle-count (count coords)})))))

(defn setup
  []
  (let [world (cave/build-world)
        ps (first (lf/particle-sys-seq world))
        p-radius (.GetRadius ps)
        cppn seed-cppn
        cppn-fn (cc/build-cppn-fn cppn)
        mdata (morpho-data cppn-fn p-radius [0 0] [1.0 1.0])
        coords (mapcat :coords mdata)]
    (dotimes [i 60]
      (clear-push world coords)
      (lf/step! world (/ 1 60.0) 8 3 3))
    (morpho-construct mdata world ps)
    (assoc bed/initial-state
           :world world
           :particle-system ps
           :particle-iterations 3
           :camera (bed/map->Camera {:width 6 :height 6 :center [0 0]}))))

(defn step
  [state]
  (bed/world-step state))

(defn draw
  [state]
  (bed/draw state true)
  (let []
    (quil/fill 255)
    (quil/text (str "Keys: (g) toggle gravity (d) damping (b) ball"
                    " (a) air")
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
  [& args]
  (quil/sketch
      :title "Xotarium"
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

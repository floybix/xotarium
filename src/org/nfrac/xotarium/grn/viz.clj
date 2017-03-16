(ns org.nfrac.xotarium.grn.viz
  (:require [org.nfrac.xotarium.grn.greans :as grn]
            [org.nfrac.xotarium.util :as util]
            [org.nfrac.xotarium.util.algo-graph :as graph]
            [org.nfrac.liquidfun.vec2d :as v2d]
            [quil.core :as quil]
            [clojure.pprint]
            [clojure.spec :as s]
            [clojure.test.check.random :as random]))

(defn random-color
  []
  [(rand-int 200) (rand-int 200) (rand-int 200)])

(defn mean-point
  [xys]
  (let [summed (reduce v2d/v-add xys)]
    (v2d/v-scale summed (/ 1.0 (count xys)))))

(defn grn-viz-data
  [grn input-syms output-syms]
  (let [cg (s/conform ::grn/grn grn)
        units (:units (::grn/elements cg))
        all-tfs (mapcat :genes units)
        input-tfs (->> (::grn/input-tfs grn)
                       (mapv #(mod % (count all-tfs))))
        output-tfs (->> (::grn/output-tfs grn)
                        (mapv #(mod % (count all-tfs))))
        input-els (mapv #(nth all-tfs %) input-tfs)
        output-tf? (set output-tfs)
        input-tf? (set input-tfs)
        id->output-sym (zipmap output-tfs output-syms)
        id->label (fn [tf-id]
                    (str (or (id->output-sym tf-id) tf-id)))
        cell (grn/grn->cell grn)
        unit-tfs (::grn/unit-tfs cell)]
    (concat
     (for [[input-sym el] (map vector input-syms input-els)]
       {:unit-id input-sym ;; ???
        :color [0 0 0]
        :promoters []
        :centroid (::grn/coords el)
        :tfs [(assoc el :label (str input-sym))]})
     (for [[unit-id unit] (map-indexed vector units)
           :let [color (random-color)
                 tf-ids (get unit-tfs unit-id)]
           :when (seq tf-ids) #_(not-every? input-tf? tf-ids)]
       {:unit-id unit-id
        :color color
        :is-output? (some output-tf? tf-ids)
        :centroid (mean-point (map ::grn/coords (:genes unit)))
        :promoters (:promoters unit)
        :tfs (->> (:genes unit)
                  (mapcat (fn [tf-id m]
                            (when-not (input-tf? tf-id)
                              [(assoc m :label (id->label tf-id))]))
                          tf-ids))}))))

(defn bounding-box
  [coords]
  (let [xs (map first coords)
        ys (map second coords)]
    [[(reduce min xs) (reduce min ys)]
     [(reduce max xs) (reduce max ys)]]))

(defn bounding-box-pad
  [coords frac]
  (let [[[xmin ymin] [xmax ymax]] (bounding-box coords)
        xspan (- xmax xmin)
        yspan (- ymax ymin)]
    [[(- xmin (* xspan frac))
      (- ymin (* yspan frac))]
     [(+ xmax (* xspan frac))
      (+ ymax (* yspan frac))]]))

(defn domain-to-px-scale
  "A scaling factor on domain coordinates to give pixels.
  Fits the camera bounds into the window, expanding these
  bounds if necessary to ensure an isometric aspect ratio."
  [w h px-width px-height]
  (let [xscale (/ px-width w)
        yscale (/ px-height h)]
    (min xscale yscale)))

(defn domain-to-px-fn
  "Returns a function to convert a point in domain coordinates to quil pixels."
  [[xmin ymin] [xmax ymax]]
  (let [scale (domain-to-px-scale (- xmax xmin) (- ymax ymin)
                                  (quil/width) (quil/height))]
    (fn [[x y]]
      [(* (- x xmin) scale)
       ;; quil has flipped y (0px at top)
       (- (quil/height)
          (* (- y ymin) scale))])))

(defn triangle
  [x-px y-px r-px v-scale h-scale]
  (let [v-r-px (* r-px v-scale)
        h-r-px (* r-px h-scale)]
    (quil/triangle (- x-px h-r-px) (+ y-px v-r-px)
                   (+ x-px h-r-px) (+ y-px v-r-px)
                   x-px (- y-px v-r-px))))

(defn draw
  [state]
  (let [data (:data state)
        bb (bounding-box-pad
            (mapcat (fn [unit]
                      (concat (map ::grn/coords (:tfs unit))
                              (map ::grn/coords (:promoters unit))))
                    data)
            0.2)
        [[xmin ymin] [xmax ymax]] bb
        ->px (domain-to-px-fn [xmin ymin] [xmax ymax])
        tri-px 6
        scale (domain-to-px-scale (- xmax xmin) (- ymax ymin)
                                  (quil/width) (quil/height))
        affinity-one-px (* scale (grn/critical-affinity-distance 1.0))
        affinity-full-px (* scale (grn/critical-affinity-distance grn/AFFINITY_EPS))]
    (quil/background 255 255 255)
    (doseq [unit data
            :let [[r g b] (:color unit)
                  qcolor (quil/color r g b)
                  [cx cy] (:centroid unit)
                  [cx-px cy-px] (->px [cx cy])]]
      ;; promoters: shaded radius of critical affinity - like dendritic arbor?
      (quil/no-stroke)
      (quil/fill qcolor (quot 256 12))
      (doseq [el (:promoters unit)
              :let [[x y :as coord] (::grn/coords el)
                    sign (::grn/sign el)
                    [x-px y-px] (->px coord)]]
        (quil/ellipse x-px y-px (* 2 affinity-one-px) (* 2 affinity-one-px)))
      )
    (doseq [unit data
            :let [[r g b] (:color unit)
                  qcolor (quil/color r g b)
                  [cx cy] (:centroid unit)
                  [cx-px cy-px] (->px [cx cy])]]
      ;; lines joining co-regulated TFs (via centroid)
      (quil/no-fill)
      (quil/stroke qcolor)
      (quil/stroke-weight 2)
      (doseq [el (:tfs unit)
              :let [[x y :as coord] (::grn/coords el)
                    sign (::grn/sign el)
                    [x-px y-px] (->px coord)]]
        (quil/line x-px y-px cx-px cy-px))
      (quil/stroke-weight 1)
      ;; lines joining clustered promoters (via centroid)
      (doseq [el (:promoters unit)
              :let [[x y :as coord] (::grn/coords el)
                    sign (::grn/sign el)
                    [x-px y-px] (->px coord)]]
        (quil/line x-px y-px cx-px cy-px))
      ;; TFs
      (quil/fill qcolor)
      (quil/no-stroke)
      (doseq [el (:tfs unit)
              :let [[x y :as coord] (::grn/coords el)
                    sign (::grn/sign el)
                    [x-px y-px] (->px coord)]]
        (triangle x-px y-px tri-px
                  (if (pos? sign) 0.5 -2.0)
                  (if (pos? sign) 2.0 0.5)))
      ;; promoters
      (quil/stroke qcolor)
      (quil/no-fill)
      (doseq [el (:promoters unit)
              :let [[x y :as coord] (::grn/coords el)
                    sign (::grn/sign el)
                    [x-px y-px] (->px coord)]]
        (triangle x-px y-px tri-px
                  (if (pos? sign) 0.5 -2.0)
                  (if (pos? sign) 2.0 0.5))
        (when (:is-output? unit)
          (quil/ellipse x-px y-px (* 2 affinity-one-px) (* 2 affinity-one-px))
          (quil/ellipse x-px y-px (* 2 affinity-full-px) (* 2 affinity-full-px)))
        ))
    (quil/fill 0 0 0)
    (quil/text-align :center)
    (doseq [unit data
            el (:tfs unit)
            :let [[x y :as coord] (::grn/coords el)
                  [x-px y-px] (->px coord)]]
      ;; labels
      (quil/text (:label el) x-px (+ y-px 12 tri-px)))))

(defn step
  [state]
  state)

(defn setup
  [data]
  (quil/frame-rate 1)
  {:data data})

(defn run
  [data]
  (quil/sketch
   :title "GRN viz"
   :host "liquidfun"
   :setup #(setup data)
   :update step
   :draw draw
   :size [800 800]
   :features [:resizable]
   :middleware [quil.middleware/fun-mode]))

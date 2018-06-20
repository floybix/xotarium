(ns org.nfrac.xotarium.grn.viz
  (:require [org.nfrac.xotarium.grn.greans :as grn]
            [org.nfrac.xotarium.util :as util]
            [org.nfrac.liquidfun.vec2d :as v2d]
            [quil.core :as quil]
            [clojure.core.async :as async]
            [clojure.pprint]
            [clojure.spec.alpha :as s]
            [clojure.test.check.random :as random]))

(defn random-color
  []
  [(rand-int 200) (rand-int 200) (rand-int 200)])

(defn random-qcolor
  []
  (quil/color (rand-int 256) 255 128))

(defn mean-point
  [xys]
  (when (seq xys)
    (let [summed (reduce v2d/v-add xys)]
      (v2d/v-scale summed (/ 1.0 (count xys))))))

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

(defn px-to-domain-fn
  "Returns a function to convert a point in quil pixels to domain coordinates."
  [[xmin ymin] [xmax ymax]]
  (let [scale (domain-to-px-scale (- xmax xmin) (- ymax ymin)
                                  (quil/width) (quil/height))]
    (fn [[px py*]]
      (let [py (- (quil/height) py*)]
        [(+ xmin (/ px scale))
         (+ ymin (/ py scale))]))))

(defn updated-coords
  "Updates the state with a modified GRN. Modifications here may be in the
  coordinates of elements, not in the genome structure."
  [state grn]
  (let [cgrn (s/conform ::grn/grn grn)
        units (:units (::grn/elements cgrn))
        cell (grn/grn->cell grn)
        bbox (bounding-box-pad
              (mapcat (fn [unit]
                        (concat (map ::grn/coords (:genes unit))
                                (map ::grn/coords (:promoters unit))))
                      units)
              0.2)]
    (assoc state
           :grn grn
           :cgrn cgrn
           :cell cell
           :bbox bbox)))

(defn init
  [{:keys [grn input-syms output-syms]}]
  (let [cgrn (s/conform ::grn/grn grn)
        units (:units (::grn/elements cgrn))
        unit-tfs (grn/index-counts (map #(count (:genes %)) units) 0)
        unit-prs (grn/index-counts (map #(count (:promoters %)) units) 0)
        all-tfs (mapcat :genes units)
        input-tfs (->> (::grn/input-tfs grn)
                       (mapv #(mod % (count all-tfs))))
        output-tfs (->> (::grn/output-tfs grn)
                        (mapv #(mod % (count all-tfs))))
        concs (vec (repeat (count all-tfs) 0.0))]
    (->
     {:grn grn
      :cgrn cgrn
      :unit-tfs unit-tfs
      :unit-prs unit-prs
      :unit-colors (mapv (fn [_]
                           (random-qcolor))
                         units)
      :concs concs
      :input (zipmap input-tfs input-syms)
      :output (zipmap output-tfs output-syms)
      }
     (updated-coords grn))))


(defn triangle
  [x-px y-px r-px v-scale h-scale]
  (let [v-r-px (* r-px v-scale)
        h-r-px (* r-px h-scale)]
    (quil/triangle (- x-px h-r-px) (+ y-px v-r-px)
                   (+ x-px h-r-px) (+ y-px v-r-px)
                   x-px (- y-px v-r-px))))

(defn draw
  [state]
  (let [{:keys [grn cgrn unit-tfs unit-prs unit-colors concs input output bbox]} state
        cell (assoc (:cell state) ::grn/concs concs)
        influences (::grn/influences cell)
        units (:units (::grn/elements cgrn))
        parameters (::grn/parameters grn)
        [[xmin ymin] [xmax ymax]] bbox
        ->px (domain-to-px-fn [xmin ymin] [xmax ymax])
        scale (domain-to-px-scale (- xmax xmin) (- ymax ymin)
                                  (quil/width) (quil/height))
        tri-px 40
        units (map (fn [unit tf-ids pr-ids]
                     (let [xy (mean-point (concat
                                           (mapcat (fn [el id]
                                                     (when-not (input id)
                                                       [(::grn/coords el)]))
                                                   (:genes unit)
                                                   tf-ids)
                                           (map ::grn/coords (:promoters unit))))]
                       (assoc unit :centroid xy
                              :tf-ids tf-ids
                              :pr-ids pr-ids)))
                   units
                   unit-tfs
                   unit-prs)
        {:keys [max-affinity affinity-eps]} parameters
        affinity-one-px (* scale (grn/critical-affinity-distance 1.0 max-affinity))
        affinity-full-px (* scale (grn/critical-affinity-distance affinity-eps max-affinity))]
    (quil/background 255 0 255)
    (doseq [[unit qcolor] (map vector units unit-colors)
            :when (not-every? input (:tf-ids unit))]
      ;; promoters: shaded radius of critical affinity - like dendritic arbor?
      (quil/no-stroke)
      (quil/fill qcolor (quot 256 12))
      (doseq [el (:promoters unit)
              :let [[x y :as coord] (::grn/coords el)
                    sign (::grn/sign el)
                    [x-px y-px] (->px coord)]]
        (quil/ellipse x-px y-px (* 2 affinity-one-px) (* 2 affinity-one-px))))
    (doseq [[unit qcolor] (map vector units unit-colors)
            :when (not-every? input (:tf-ids unit))
            :let [unit-act (grn/unit-activation cell (:pr-ids unit))
                  output-unit? (some output (:tf-ids unit))
                  [cx cy] (:centroid unit)
                  [cx-px cy-px] (->px [cx cy])]]
      ;; lines joining co-regulated TFs (via centroid)
      (quil/no-fill)
      (quil/stroke qcolor (quot 256 4))
      (quil/stroke-weight (+ 1 (* unit-act 1.0)))
      (doseq [[el id] (map vector (:genes unit) (:tf-ids unit))
              :when (not (input id))
              :let [[x y :as coord] (::grn/coords el)
                    sign (::grn/sign el)
                    [x-px y-px] (->px coord)]]
        (quil/line x-px y-px cx-px cy-px))
      (quil/stroke qcolor)
      (quil/stroke-weight 1)
      ;; lines joining clustered promoters (via centroid)
      (doseq [el (:promoters unit)
              :let [[x y :as coord] (::grn/coords el)
                    sign (::grn/sign el)
                    [x-px y-px] (->px coord)]]
        (quil/line x-px y-px cx-px cy-px))
      ;; promoters
      (quil/stroke qcolor)
      (quil/no-fill)
      (doseq [[el id] (map vector (:promoters unit) (:pr-ids unit))
              :let [[x y :as coord] (::grn/coords el)
                    sign (::grn/sign el)
                    [x-px y-px] (->px coord)
                    act (-> (grn/promoter-activation id influences concs)
                            (max 0.0)
                            (min 10.0))]]
        (triangle x-px y-px (+ 3 (* tri-px act 0.1))
                  (if (pos? sign) 0.25 -1.0)
                  (if (pos? sign) 1.0 0.25))
        (when output-unit?
          (quil/ellipse x-px y-px (* 2 affinity-one-px) (* 2 affinity-one-px))
          (quil/ellipse x-px y-px (* 2 affinity-full-px) (* 2 affinity-full-px)))))
    (doseq [[unit qcolor] (map vector units unit-colors)]
      ;; TFs
      (quil/fill qcolor)
      (quil/no-stroke)
      (doseq [[el id] (map vector (:genes unit) (:tf-ids unit))
              :let [qcolor (if (input id) (quil/color 0 0 0) qcolor)
                    [x y :as coord] (::grn/coords el)
                    sign (::grn/sign el)
                    [x-px y-px] (->px coord)]]
        (triangle x-px y-px (+ 3 (* tri-px (get concs id)))
                  (if (pos? sign) 0.25 -1.0)
                  (if (pos? sign) 1.0 0.25))))
    (quil/fill 0 0 0)
    (quil/text-align :center)
    (doseq [unit units
            [el id] (map vector (:genes unit) (:tf-ids unit))
            :let [[x y :as coord] (::grn/coords el)
                  [x-px y-px] (->px coord)
                  label (str (or (output id) (input id) id))]]
      ;; labels
      (quil/text label x-px (+ y-px (* tri-px 0.5))))))

(defn step
  [state]
  (if-let [new-concs (async/poll! (:activity-c state))]
    (assoc state :concs new-concs)
    state))

(defn pick-element
  [state [mx-px my-px]]
  (let [{:keys [cgrn bbox unit-tfs unit-prs]} state
        [[xmin ymin] [xmax ymax]] bbox
        px-> (px-to-domain-fn [xmin ymin] [xmax ymax])
        [mx my] (px-> [mx-px my-px])
        scale (domain-to-px-scale (- xmax xmin) (- ymax ymin)
                                  (quil/width) (quil/height))
        els (::grn/elements (:grn state))]
    (->> els
         (map-indexed (fn [i el]
                        (let [[x y] (::grn/coords el)
                              dist (v2d/v-dist [x y] [mx my])]
                          [i dist])))
         (apply min-key second)
         (first))))

(defn mouse-dragged
  [state event]
  (if-let [drag-info (:dragging state)]
    (let [{:keys [grn cgrn unit-tfs unit-prs bbox control-c]} state
          [[xmin ymin] [xmax ymax]] bbox
          px-> (px-to-domain-fn [xmin ymin] [xmax ymax])
          xy-px [(:x event) (:y event)]
          xy (px-> xy-px)
          i drag-info
          new-grn (assoc-in grn [::grn/elements i ::grn/coords] xy)]
      (async/put! control-c new-grn)
      (updated-coords state new-grn))
    state))

(defn left-mouse-pressed
  [state event]
  (if (:dragging state)
    state
    (if-let [drag-info (pick-element state [(:x event) (:y event)])]
      (-> (assoc state :dragging drag-info)
          ;; need to make sure this is a vector as coords will be updated
          (update-in [:grn ::grn/elements] vec))
      state)))

(defn mouse-pressed
  [state event]
  (case (:button event)
    :left (left-mouse-pressed state event)
    state))

(defn mouse-released
  [state event]
  (dissoc state :dragging))

(defn run
  [m control-c activity-c]
  (quil/sketch
   :title "GRN viz"
   :host "liquidfun"
   :setup (fn []
            (quil/frame-rate 10)
            (quil/color-mode :hsb)
            (assoc (init m)
                   :control-c control-c
                   :activity-c activity-c))
   :update step
   :draw draw
   :mouse-pressed mouse-pressed
   :mouse-released mouse-released
   :mouse-dragged mouse-dragged
   :size [800 800]
   :features [:resizable]
   :middleware [quil.middleware/fun-mode]))

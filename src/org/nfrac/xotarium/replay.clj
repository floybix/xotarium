(ns org.nfrac.xotarium.replay
  (:require [org.nfrac.xotarium.evo :as evo]
            [org.nfrac.xotarium.grn-creature :as grncre]
            [org.nfrac.xotarium.creature :as cre]
            [org.nfrac.xotarium.cppn :as cppn]
            [org.nfrac.xotarium.cppn-compile :as cc]
            [org.nfrac.xotarium.grn.greans :as grn]
            [org.nfrac.xotarium.plant-cave :as cave]
            [org.nfrac.xotarium.grn.viz :as grnviz]
            [org.nfrac.liquidfun.testbed :as bed]
            [org.nfrac.liquidfun.core :as lf :refer [body! joint!
                                                     particle-system!]]
            [org.nfrac.liquidfun.vec2d :as v2d]
            [org.nfrac.xotarium.util :as util]
            [org.nfrac.xotarium.util.algo-graph :as graph]
            [quil.core :as quil :include-macros true]
            [quil.middleware]
            [clojure.pprint]
            [clojure.java.io :as io]
            [clojure.core.async :as async]
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

(defn to-file
  [file form]
  (with-open [w (io/writer file)]
    (print-dup form w)))

(defn from-file
  [file]
  (with-open [r (java.io.PushbackReader. (io/reader file))]
    (read r)))

(def warmup-steps 60)

(defn setup-current
  [genome seed]
  (when-let [state (grncre/setup genome seed)]
    (assoc state ::step-i 0
           ::control-c (async/chan (async/sliding-buffer 1))
           ::activity-c (async/chan (async/sliding-buffer 1)))))

(defn control-grn
  [state control-c]
  (if-let [new-grn (async/poll! control-c)]
    (-> state
        (assoc-in [:creature :grn] new-grn)
        (assoc-in [:creature :grn-cell] (grn/grn->cell new-grn)))
    state))

(defn step-current
  [pre-state]
  (let [state (-> (grncre/step pre-state)
                  (update ::step-i inc)
                  (control-grn (::control-c pre-state)))]
    (async/put! (::activity-c state)
                (-> state :creature :tri-concs first val))
    (if (>= (::step-i state) warmup-steps)
      (let [creature (:creature state)
            grn-cell (:grn-cell creature)
            tri-concs (:tri-concs creature)
            n-cells (count tri-concs)
            concs-1 (first (vals tri-concs))
            tf-ranges (get pre-state ::tf-ranges
                           (mapv #(vector % %) concs-1))
            tf-dranges (get pre-state ::tf-dranges
                            (vec (repeat (count concs-1) [0.0 0.0])))
            pre-tri-concs (:tri-concs (:creature pre-state))
            [tf-vs tf-dvs] (loop [cell-concs (vals tri-concs)
                                  cell-pre-concs (map pre-tri-concs (keys tri-concs))
                                  tf-vs tf-ranges
                                  tf-dvs tf-dranges]
                             (if-let [concs (first cell-concs)]
                               (let [pre-concs (first cell-pre-concs)]
                                 (recur (rest cell-concs)
                                        (rest cell-pre-concs)
                                        (mapv (fn [[lo hi] v]
                                                [(if (< v lo) v lo)
                                                 (if (> v hi) v hi)])
                                              tf-vs concs)
                                        (mapv (fn [[lo hi] v pre-v]
                                                (let [d (- v pre-v)]
                                                  [(if (< d lo) d lo)
                                                   (if (> d hi) d hi)]))
                                              tf-dvs concs pre-concs)))
                               [tf-vs tf-dvs]))]
        (assoc state
                 ::tf-ranges tf-vs
                 ::tf-dranges tf-dvs))
      ;; warmup
      state)))

(defn format-lims
  [[lo hi]]
  (format "[%.2f, %.2f]" lo hi))

(defn format-lims+
  [[lo hi]]
  (format "[%+.2f, %+.2f]" lo hi))

(defn my-key-press
  [mstate event]
  (let [state (:current mstate)
        ps ^liquidfun$b2ParticleSystem (:particle-system state)]
    (case (:key event)
      :n (loop [i (:genome-index mstate)]
           (let [i+ (mod (inc i) (count (:genomes mstate)))
                 genome (get (:genomes mstate) i+)
                 seed (:seed mstate)]
             (if-let [state (setup-current genome seed)]
               (assoc mstate
                      :genome-index i+
                      :current state)
               (recur (inc i)))))
      :i (let [tf-vs (::tf-ranges state)
               tf-dvs (::tf-dranges state)
               n-tfs (count tf-vs)
               creature (:creature state)
               cell-form (:grn-cell creature)
               input-tfs (->> (::grn/input-tfs cell-form)
                              (mapv #(mod % n-tfs)))
               output-tfs (->> (::grn/output-tfs cell-form)
                               (mapv #(mod % n-tfs)))
               input-syms grncre/beh-inputs
               output-syms grncre/beh-outputs
               inp->i (zipmap input-syms input-tfs)
               i->out (zipmap output-tfs output-syms)
               inps-by-range (->> input-syms
                                  (sort-by (fn [k]
                                             (let [id (inp->i k)
                                                   [lo hi] (get tf-vs id)]
                                               (- hi lo))) >))
               real-tfs (remove (set input-tfs) (range n-tfs))
               tfs-by-range (->> real-tfs
                                 (sort-by (fn [id]
                                            (let [[lo hi] (get tf-vs id)]
                                              (- hi lo))) >))]
           (println)
           (println "step" (::step-i state) "time" (format "%.2f" (:time state)))
           (println "outputs:" i->out)
           (println "__inputs__")
           (doseq [k inps-by-range
                   :let [id (get inp->i k)]]
             (println (format "%15s" k)
                      (format-lims (get tf-vs id)) "deltas"
                      (format-lims+ (get tf-dvs id))))
           (println)
           (println "__TFs__")
           (doseq [id tfs-by-range]
             (println (format "%15s" (or (i->out id) id))
                      (format-lims (get tf-vs id)) "deltas"
                      (format-lims+ (get tf-dvs id))))
           (println)
           mstate)
      :s (let []
           (update mstate :draw-cell-signals? #(not %)))
      :v (let [grn (:grn (:creature state))]
           (grnviz/run {:grn grn
                        :input-syms grncre/beh-inputs
                        :output-syms grncre/beh-outputs}
             (::control-c state)
             (::activity-c state))
            mstate)
      :b (do
           (body! (:world state) {}
                  {:shape (lf/circle 0.25)
                   :restitution 0.1
                   :density 1.0})
           mstate)
      :g (do
           (.SetGravityScale ps (if (== 1.0 (.GetGravityScale ps))
                                  0.0 1.0))
           mstate)
      :d (do
           (.SetDamping ps (if (== 1.0 (.GetDamping ps))
                             0.0 1.0))
           mstate)
      ;; otherwise pass on to testbed
      (update mstate :current bed/key-press event))))

(defn setup
  [file]
  (let [{:keys [beh-archive popn seed]} (from-file file)
        beh-genomes (->> (map :representative (vals beh-archive))
                         (sort-by :generation)
                         (vec))
        genomes (into beh-genomes
                      (map :genome popn))
        genome-index 0]
    {:seed seed
     :genomes genomes
     :genome-index genome-index
     :n-behs (count beh-archive)
     :current (setup-current (get genomes genome-index) seed)}))

(defn draw-cell-signals
  [state]
  (let [ps ^liquidfun$b2ParticleSystem (:particle-system state)
        cam (:camera state)
        px-scale (bed/world-to-px-scale cam)
        ;; inlined version of world-to-px
        [cx cy] (:center cam)
        x-left (- cx (* 0.5 (:width cam)))
        y-bottom (- cy (* 0.5 (:height cam)))
        y-top (+ y-bottom (:height cam))
        ;; creature info
        creature (:creature state)
        cell-form (:grn-cell creature)
        tri-concs (:tri-concs creature)
        n-tfs (count (first (vals tri-concs)))
        ;; messengers must start at output index 1:
        [msgr-a msgr-b msgr-c] (->> (drop 1 (::grn/output-tfs cell-form))
                                    (map #(mod % n-tfs)))
        handle->px (fn [^liquidfun$b2ParticleHandle h]
                     (let [i (.GetIndex h)
                           x (.GetParticlePositionX ps i)
                           y (.GetParticlePositionY ps i)
                           x-px (* (- x x-left) px-scale)
                           y-px (* (- y-top y) px-scale)]
                       [x-px y-px]))
        ]
    (quil/no-stroke)
    (doseq [[handles concs] tri-concs]
      (let [[h1 h2 h3] handles
            [x1 y1] (handle->px h1)
            [x2 y2] (handle->px h2)
            [x3 y3] (handle->px h3)
            r (-> (get concs msgr-a) (* 255))
            g (-> (get concs msgr-b) (* 255))
            b (-> (get concs msgr-c) (* 255))
            alpha 255]
        (quil/fill r g b alpha)
        (quil/triangle x1 y1 x2 y2 x3 y3)))))

(defn draw
  [mstate]
  (bed/draw (:current mstate) true)
  (when (:draw-cell-signals? mstate)
    (draw-cell-signals (:current mstate)))
  (let [gi (:genome-index mstate)
        genome (get (:genomes mstate) gi)
        n-behs (:n-behs mstate)
        n-popn (- (count (:genomes mstate)) n-behs)]
    (quil/fill 255)
    (quil/text (str "Replaying "
                    (if (< gi n-behs)
                      (str (inc gi) " / " n-behs " beh.")
                      (str (- gi n-behs) " / " n-popn " indiv."))
                    " (gen " (:generation genome) ")."
                    " Keys: (n) next behaviour, (i) print TF info,"
                    " (v) GRN viz, (s) draw signals")
               10 10)))

(defn run
  [file & args]
  (quil/sketch
   :title "Xotarium"
   :host "liquidfun"
   :setup #(do
             (quil/frame-rate 30)
             (setup file))
   :update (fn [mstate]
             (if (:paused? (:current mstate))
               mstate
               (update mstate :current step-current)))
   :draw (fn [mstate]
           (if (zero? (mod (quil/frame-count) 2))
             (draw mstate)
             mstate))
   :key-typed my-key-press
   :mouse-pressed #(update % :current bed/mouse-pressed %2)
   :mouse-released #(update % :current bed/mouse-released %2)
   :mouse-dragged #(update % :current bed/mouse-dragged %2)
   :mouse-wheel #(update % :current bed/mouse-wheel %2)
   :size [600 500]
   :features [:resizable]
   :middleware [quil.middleware/fun-mode]))

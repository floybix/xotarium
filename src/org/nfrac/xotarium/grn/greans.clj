(ns org.nfrac.xotarium.grn.greans
  "An implementation of GReaNs (Gene Regulatory evolving artificial Networks)
  similar to that described in the paper by Michal Joachimczak and Borys Wróbel:
  Co-evolution of morphology and control of soft-bodied multicellular animats."
  (:require [clojure.spec :as s]
            [clojure.spec.gen :as gen]
            [org.nfrac.xotarium.util :as util]
            [clojure.test.check.random :as random]))

(def INIT_MAX_COORD 5.0)

(def parameter-defaults
  {:mut-coord-mag 1.0
   :mut-coord-prob 0.05
   :mut-sign-prob 0.05
   :mut-type-prob 0.05
   :switch-inout-prob 0.05
   :dup-prob 0.2
   :del-prob 0.2})

;;; evolvable linear form of grn

(def element-types #{::promoter ::gene})

(s/def ::element-type element-types)

(s/def ::sign #{-1 1})

(s/def ::coordinate
  (-> (s/double-in :infinite? false :NaN? false)
      (s/with-gen #(gen/double* {:min 0.0 :max INIT_MAX_COORD :NaN? false}))))

(s/def ::coords (s/coll-of ::coordinate :count 2))

(s/def ::element
  (s/keys :req [::element-type
                ::sign
                ::coords]))

(s/def ::promoter
  (s/and ::element #(= (::element-type %) ::promoter)))

(s/def ::gene
  (s/and ::element #(= (::element-type %) ::gene)))

(s/def ::regulatory-unit
  (s/cat :promoters (s/+ ::promoter)
         :genes (s/+ ::gene)))

(s/def ::elements
  (s/cat :head-ignored (s/* ::gene)
         :units (s/+ ::regulatory-unit)
         :tail-ignored (s/* ::promoter)))

(s/def ::output-tfs (s/coll-of nat-int?))
(s/def ::input-tfs (s/coll-of nat-int?))

(s/def ::parameters
  (s/keys))

(s/def ::grn
  (s/keys :req [::elements
                ::input-tfs
                ::output-tfs
                ::parameters]))

(comment
  (s/exercise ::gene)
  (gen/generate (s/gen ::grn)))

(defn random-element
  [rng type]
  (let [[rng1 rng2 rng3] (random/split-n rng 3)]
    {::element-type type
     ::sign (util/rand-nth rng1 [-1 1])
     ::coords [(util/rand rng2 0 INIT_MAX_COORD)
               (util/rand rng3 0 INIT_MAX_COORD)]}))

(defn random-unit
  [rng]
  (let [[rng1 rng2 rng3 rng4] (random/split-n rng 4)
        npr (+ 1 (util/rand-int rng1 3))
        ntf (+ 1 (util/rand-int rng2 3))]
    (concat (map #(random-element % ::promoter) (random/split-n rng3 npr))
            (map #(random-element % ::gene) (random/split-n rng4 ntf)))))

(defn random-grn
  [n-in n-out rng]
  (let [n-units (+ n-in n-out 1)
        [rng rng*] (random/split rng)
        es (mapcat #(random-unit %) (random/split-n rng* n-units))
        n-tfs (count (filter #(= ::gene (::element-type %)) es))
        tf-ids (util/shuffle rng (range n-tfs))]
    {::elements es
     ::input-tfs (take n-in tf-ids)
     ::output-tfs (take n-out (drop n-in tf-ids))
     ::parameters parameter-defaults}))

;;;  "cell" / usable form of grn

(s/def ::tf-id nat-int?)
(s/def ::promoter-id nat-int?)

(def MAX_AFFINITY 10.0)
(def AFFINITY_EPS 0.1)

(s/def ::affinity (s/double-in :min 0.0 :max MAX_AFFINITY :NaN? false))
(s/def ::influence (s/double-in :min (- MAX_AFFINITY) :max MAX_AFFINITY :NaN? false))
(s/def ::concentration (s/double-in :min 0.0 :max 1.0 :NaN? false))

;; vector keyed by unit id
(s/def ::unit-promoters
  (s/every (s/every ::promoter-id) :kind vector?))

;; vector keyed by unit id
(s/def ::unit-tfs
  (s/every (s/every ::tf-id) :kind vector?))

;; vector keyed by promoter id
(s/def ::influences
  (s/every (s/every-kv ::tf-id ::influence) :kind vector?))

;; vector keyed by tf id
(s/def ::concs
  (s/every ::concentration :kind vector?))

(s/def ::cell
  (s/keys :req [::unit-promoters
                ::unit-tfs
                ::influences
                ::concs
                ::input-tfs
                ::output-tfs
                ::parameters]))

(defn index-counts
  [counts init-offset]
  (let [offsets (reductions + init-offset counts)]
    (mapv (fn [o n] (range o (+ o n)))
          offsets counts)))

(defn affinity
  [[ax ay] [bx by]]
  (let [d (Math/sqrt (+ (Math/abs (double (- bx ax)))
                        (Math/abs (double (- by ay)))))]
    (let [aff (-> (Math/exp (* -1.0 d))
                  (min MAX_AFFINITY))]
      (when (>= aff AFFINITY_EPS)
        aff))))

(defn promoter-influences
  [pr tfs]
  (let [pr-coords (::coords pr)]
    (reduce (fn [m [tf-id tf]]
              (if-let [aff (affinity pr-coords (::coords tf))]
                (assoc m tf-id (* aff (::sign pr) (::sign tf)))
                m))
            {}
            (map-indexed vector tfs))))

(defn grn->cell
  [grn]
  (let [cg (s/conform ::grn grn)
        _ (when (= cg ::s/invalid) (s/explain ::grn grn))
        units (:units (::elements cg))
        unit-promoters (index-counts (map #(count (:promoters %)) units) 0)
        unit-tfs (index-counts (map #(count (:genes %)) units) 0)
        all-prs (mapcat :promoters units)
        all-tfs (mapcat :genes units)
        ;; pr -> {tf -> infl}
        pr-influences (mapv #(promoter-influences % all-tfs) all-prs)
        ]
    (merge
     (select-keys cg [::input-tfs ::output-tfs ::parameters])
     {::unit-promoters unit-promoters
      ::unit-tfs unit-tfs
      ::influences pr-influences
      ::concs (vec (repeat (count all-tfs) 0.0))})))

(defn promoter-activation
  [pr-id influences concs]
  (reduce-kv (fn [a tf-id infl]
               (let [c (get concs tf-id)]
                 (+ a (* c infl))))
             0.0
             (get influences pr-id)))

(defn unit-activation
  [cell pr-ids]
  (let [influences (::influences cell)
        concs (::concs cell)]
    (reduce (fn [a pr-id]
              (+ a (promoter-activation pr-id influences concs)))
            0.0
            pr-ids)))

(defn conc-rate-of-change
  ^double [^double conc ^double unit-activation]
  (let [e (Math/exp (- unit-activation))]
    (- (/ (- 1.0 e)
          (+ 1.0 e))
       conc)))

(defn vector-merge
  [v kvs]
  (persistent!
   (reduce-kv assoc! (transient v) kvs)))

(defn step
  [cell input-concs dt]
  (let [dt (double dt)
        n-tfs (count (::concs cell))
        concs (vector-merge (::concs cell)
                            (zipmap (map #(mod % n-tfs) (::input-tfs cell))
                                    input-concs))
        influences (::influences cell)
        unit-promoters (::unit-promoters cell)
        unit-tfs (::unit-tfs cell)
        n-units (count unit-tfs)]
    (loop [unit-i 0
           new-concs (transient concs)]
      (if (< unit-i n-units)
        (let [pr-ids (get unit-promoters unit-i)
              unit-a (max 0.0
                          (transduce
                           (map #(promoter-activation % influences concs))
                           + 0 pr-ids))
              tf-ids (get unit-tfs unit-i)
              ncs (reduce (fn [ncs tf-id]
                            (let [c (double (get new-concs tf-id))
                                  dc (conc-rate-of-change c unit-a)]
                              (assoc! ncs tf-id (+ c (* dc dt)))))
                          new-concs
                          tf-ids)]
          (recur (inc unit-i)
                 ncs))
        ;; done
        (assoc cell ::concs (persistent! new-concs))))))

(defn cell-outputs
  [cell]
  (let [concs (::concs cell)
        n (count concs)]
    (mapv #(get concs (mod % n))
          (::output-tfs cell))))

;; TODO: rng

(comment
  (def g (gen/generate (s/gen ::grn)))
  (def g (random-grn 2 2))
  (def cell (grn->cell g))
  (def csteps (iterate #(step % [1.0 1.0] 0.15) cell))
  (map println (map cell-outputs (take 10 csteps)))
  (def g2 (mutate g))
  (def cell2 (grn->cell g2))
)

;;; mutation etc

(defn perturb-coords
  [[x y] max-mag]
  (let [angle (rand (* 2.0 Math/PI))
        mag (rand max-mag)]
    [(+ x (* mag (Math/cos angle)))
     (+ y (* mag (Math/sin angle)))]))

(defn element-mutate
  [el max-mag cprob sprob tprob]
  (cond-> el
    (> (rand) cprob)
    (update ::coords perturb-coords max-mag)
    (> (rand) sprob)
    (update ::sign * -1)
    (> (rand) tprob)
    (update ::element-type #(first (disj element-types %)))))

(defn mutate-elements
  [grn]
  (let [{max-mag :mut-coord-mag
         cprob :mut-coord-prob
         sprob :mut-sign-prob
         tprob :mut-type-prob} (::parameters grn)]
    (update grn ::elements
            (fn [es]
              (map #(element-mutate % max-mag cprob sprob tprob) es)))))

(defn genome-duplication
  [grn]
  (let [es (::elements grn)
        i0 (rand-int (count es))
        n (rand-int (- (count es) i0))
        to (rand-int (inc (count es)))]
    (assoc grn ::elements
           (concat (take to es)
                   (take n (drop i0 es))
                   (drop to es)))))

(defn genome-deletion
  [grn]
  (let [es (::elements grn)
        i0 (rand-int (count es))
        n (rand-int (- (count es) i0))]
    (assoc grn ::elements
           (concat (take i0 es)
                   (drop (+ i0 n) es)))))

(defn genome-switch-io
  [grn prob]
  (let [ins (::input-tfs grn)
        outs (::output-tfs grn)]
    (assoc grn
           ::input-tfs (map #(if (> (rand) prob) (rand-int 100) %)
                            ins)
           ::output-tfs (map #(if (> (rand) prob) (rand-int 100) %)
                             outs))))

(defn mutate-structure
  [grn]
  (let [{:keys [switch-inout-prob dup-prob del-prob]} (::parameters grn)]
    (cond-> grn
      (> (rand) dup-prob)
      (genome-duplication)
      (> (rand) del-prob)
      (genome-deletion)
      true
      (genome-switch-io switch-inout-prob))))

(defn mutate
  [grn]
  (-> grn
      (mutate-elements)
      (mutate-structure)))

(defn crossover
  [g1 g2]
  (let [es1 (::elements g1)
        es2 (::elements g2)
        i1 (rand-int (count es1))
        ;; aligned (equal) crossover if no duplications / deletions
        i2 (mod i1 (count es2))]
    (assoc g1 ::elements
           (concat (take i1 es1)
                   (drop i2 es2)))))

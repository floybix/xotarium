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

(s/def ::output-tfs (s/coll-of nat-int? :min-count 1))
(s/def ::input-tfs (s/coll-of nat-int? :min-count 1))

(s/def ::parameters
  (s/keys))

(s/def ::grn
  (s/and
   (s/keys :req [::elements
                 ::input-tfs
                 ::output-tfs
                 ::parameters])
   #(>= (count (:units (::elements %)))
        (+ (count (::input-tfs %))
           (count (::output-tfs %))))))

(comment
  (s/exercise ::gene)
  (gen/generate (s/gen ::grn)))

(defn random-element
  [rng type]
  (let [[r1 r2 r3] (random/split-n rng 3)]
    {::element-type type
     ::sign (util/rand-nth r1 [-1 1])
     ::coords [(util/rand r2 0 INIT_MAX_COORD)
               (util/rand r3 0 INIT_MAX_COORD)]}))

(defn random-unit
  [rng]
  (let [[r1 r2 r3 r4] (random/split-n rng 4)
        npr (+ 1 (util/rand-int r1 3))
        ntf (+ 1 (util/rand-int r2 3))]
    (concat (map #(random-element % ::promoter) (random/split-n r3 npr))
            (map #(random-element % ::gene) (random/split-n r4 ntf)))))

(defn random-grn
  [n-in n-out rng]
  (let [n-units (+ n-in n-out 1)
        [r1 r2] (random/split-n rng 2)
        es (mapcat #(random-unit %) (random/split-n r1 n-units))
        n-tfs (count (filter #(= ::gene (::element-type %)) es))
        tf-ids (util/shuffle r1 (range n-tfs))]
    {::elements es
     ::input-tfs (take n-in tf-ids)
     ::output-tfs (take n-out (drop n-in tf-ids))
     ::parameters parameter-defaults}))

(s/fdef random-element
        :args (s/cat :rng ::util/rng
                     :type ::element-type)
        :ret ::element)

(s/fdef random-grn
        :args (s/cat :n-in pos-int?
                     :n-out pos-int?
                     :rng ::util/rng)
        :ret ::grn)

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

(s/fdef grn->cell
        :args (s/cat :grn ::grn)
        :ret ::cell)

(defn promoter-activation
  ^double [pr-id influences concs]
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
  [v is xs]
  (let [n (count v)]
    (loop [is is
           xs xs
           v (transient v)]
      (if-let [i (first is)]
        (let [x (first xs)
              i (mod i n)]
          (recur (rest is)
                 (rest xs)
                 (assoc! v i x)))
        (persistent! v)))))

(defn step
  [cell input-concs dt]
  (let [dt (double dt)
        n-tfs (count (::concs cell))
        concs (vector-merge (::concs cell) (::input-tfs cell) input-concs)
        influences (::influences cell)
        unit-promoters (::unit-promoters cell)
        unit-tfs (::unit-tfs cell)
        n-units (count unit-tfs)]
    (loop [unit-i 0
           unit-promoters (seq unit-promoters)
           new-concs (transient concs)]
      (if (< unit-i n-units)
        (let [pr-ids (first unit-promoters)
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
                 (rest unit-promoters)
                 ncs))
        ;; done
        (assoc cell ::concs (persistent! new-concs))))))

(s/fdef step
        :args (s/and
               (s/cat :cell ::cell
                      :input-concs (s/coll-of ::concentration)
                      :dt (s/and (s/double-in :min 0 :NaN? false) pos?))
               #(= (count (:input-concs %))
                   (count (::input-tfs (:cell %)))))
        :ret ::cell)

(defn cell-outputs
  [cell]
  (let [concs (::concs cell)
        n (count concs)]
    (mapv #(get concs (mod % n))
          (::output-tfs cell))))

(s/fdef cell-outputs
        :args (s/cat :cell ::cell)
        :ret (s/coll-of ::concentration :kind vector?))

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

(def MUT_ATTEMPTS 8)

(s/def ::probability
  (s/double-in :min 0.0 :max 1.0 :NaN? false))

(defn perturb-coords
  [[x y] max-mag rng]
  (let [[r1 r2] (random/split-n rng 2)
        angle (util/rand r1 (* 2.0 Math/PI))
        mag (util/rand r2 max-mag)]
    [(+ x (* mag (Math/cos angle)))
     (+ y (* mag (Math/sin angle)))]))

(defn element-mutate
  [el max-mag cprob sprob tprob rng]
  (let [[r1 r2 r3 r4] (random/split-n rng 4)]
    (cond-> el
      (< (random/rand-double r1) cprob)
      (update ::coords perturb-coords max-mag r2)
      (< (random/rand-double r3) sprob)
      (update ::sign * -1)
      (< (random/rand-double r4) tprob)
      (update ::element-type #(first (disj element-types %))))))

(s/fdef element-mutate
        :args (s/cat :el ::element
                     :max-mag number?
                     :cprob ::probability
                     :sprob ::probability
                     :tprob ::probability
                     :rng ::util/rng)
        :ret ::element)

(defn mutate-elements*
  [grn rng]
  (let [{max-mag :mut-coord-mag
         cprob :mut-coord-prob
         sprob :mut-sign-prob
         tprob :mut-type-prob} (::parameters grn)
         es (::elements grn)]
    (assoc grn ::elements
           (map (fn [e r]
                  (element-mutate e max-mag cprob sprob tprob r))
                es
                (random/split-n rng (count es))))))

(defn mutate-elements
  [grn rng]
  (loop [attempt 1
         rng rng]
    (let [[rng rng*] (random/split rng)
          grn2 (mutate-elements* grn rng*)]
      (if (s/valid? ::grn grn2)
        grn2
        (if (>= attempt MUT_ATTEMPTS)
          (do
            (println "mutate-elements invalid after" MUT_ATTEMPTS "attempts.")
            (s/explain ::grn grn2)
            grn)
          (recur (inc attempt) rng))))))

(s/fdef mutate-elements
        :args (s/cat :grn ::grn
                     :rng ::util/rng)
        :ret ::grn)

(defn genome-duplication
  [grn rng]
  (let [es (::elements grn)
        [r1 r2 r3] (random/split-n rng 3)
        i0 (util/rand-int r1 (count es))
        n (util/rand-int r2 1 (- (count es) i0))
        to (util/rand-int r3 (inc (count es)))]
    (assoc grn ::elements
           (concat (take to es)
                   (take n (drop i0 es))
                   (drop to es)))))

(s/fdef genome-duplication
        :args (s/cat :grn ::grn
                     :rng ::util/rng)
        :ret ::grn)

(defn genome-deletion*
  [grn rng]
  (let [es (::elements grn)
        min-possible-es (* 2 (+ 1 (count (::input-tfs grn))
                                (count (::output-tfs grn))))
        [r1 r2] (random/split-n rng 2)
        max-possible-deletion (- (count es) 1 min-possible-es)]
    (if (pos? max-possible-deletion)
      (let [n (util/rand-int r1 1 (inc max-possible-deletion))
            i0 (util/rand-int r2 (- (count es) n))]
        (assoc grn ::elements
              (concat (take i0 es)
                      (drop (+ i0 n) es)))))))

(defn genome-deletion
  [grn rng]
  #_
  {:pre [(s/valid? ::grn grn)]}
  (loop [attempt 1
         rng rng]
    (let [[rng rng*] (random/split rng)
          grn2 (genome-deletion* grn rng*)]
      (if (s/valid? ::grn grn2)
        grn2
        (if (>= attempt MUT_ATTEMPTS)
          (do
            (println "genome-deletion invalid after" MUT_ATTEMPTS "attempts.")
            grn)
          (recur (inc attempt) rng))))))

(s/fdef genome-deletion
        :args (s/cat :grn ::grn
                     :rng ::util/rng)
        :ret ::grn)

(defn genome-switch-io
  [grn prob rng]
  (let [ins (::input-tfs grn)
        outs (::output-tfs grn)
        [r1 r2 r3 r4] (random/split-n rng 4)]
    (assoc grn
           ::input-tfs (map (fn [id ra rb]
                              (if (< (random/rand-double ra) prob)
                                (util/rand-int rb 100)
                                id))
                            ins
                            (random/split-n r1 (count ins))
                            (random/split-n r2 (count ins)))
           ::output-tfs (map (fn [id ra rb]
                               (if (< (random/rand-double ra) prob)
                                 (util/rand-int rb 100)
                                 id))
                             outs
                             (random/split-n r3 (count outs))
                             (random/split-n r4 (count outs))))))

(s/fdef genome-switch-io
        :args (s/cat :grn ::grn
                     :prob ::probability
                     :rng ::util/rng)
        :ret ::grn)

(defn mutate-structure
  [grn rng]
  (let [{:keys [switch-inout-prob dup-prob del-prob]} (::parameters grn)
        [r1 r2 r3 r4 r5] (random/split-n rng 5)]
    (cond-> grn
      (< (random/rand-double r1) dup-prob)
      (genome-duplication r2)
      (< (random/rand-double r3) del-prob)
      (genome-deletion r4)
      true
      (genome-switch-io switch-inout-prob r5))))

(defn mutate
  [grn rng]
  (let [[r1 r2] (random/split-n rng 2)]
    (-> grn
        (mutate-elements r1)
        (mutate-structure r2))))

(defn crossover*
  [g1 g2 rng]
  (let [es1 (::elements g1)
        es2 (::elements g2)
        i1 (util/rand-int rng (count es1))
        ;; aligned (equal) crossover if no duplications / deletions
        i2 (mod i1 (count es2))]
    (assoc g1 ::elements
           (concat (take i1 es1)
                   (drop i2 es2)))))

(defn crossover
  [g1 g2 rng]
  #_
  {:pre [(s/valid? ::grn g1)
         (s/valid? ::grn g2)]}
  (loop [attempt 1
         rng rng]
    (let [[rng rng*] (random/split rng)
          grn2 (crossover* g1 g2 rng*)]
      (if (s/valid? ::grn grn2)
        grn2
        (if (>= attempt MUT_ATTEMPTS)
          (do
            (println "crossover invalid after" MUT_ATTEMPTS "attempts.")
            (s/explain ::grn grn2)
            g1)
          (recur (inc attempt) rng))))))

(s/fdef crossover
        :args (s/cat :g1 ::grn
                     :g2 ::grn
                     :rng ::util/rng)
        :ret ::grn)

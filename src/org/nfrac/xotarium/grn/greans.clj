(ns org.nfrac.xotarium.grn.greans
  "An implementation of GReaNs (Gene Regulatory evolving artificial Networks)
  as described in papers by Michal Joachimczak and Borys Wróbel. Specifically:
  Co-evolution of morphology and control of soft-bodied multicellular animats."
  (:require [clojure.spec :as s]
            [clojure.spec.gen :as gen]))

;;; evolvable form of genome

(def element-types #{::signal ::promoter ::gene})

(s/def ::coordinate
  (-> (s/double-in :infinite? false :NaN? false)
      (s/with-gen #(gen/double* {:min 0.0 :max 10.0 :NaN? false}))))

(s/def ::element
  (s/cat :type element-types
         :sign #{-1 1}
         :coords (s/coll-of ::coordinate :count 2)))

(s/def ::signal
  (s/and ::element #(= (:type %) ::signal)))

(s/def ::promoter
  (s/and ::element #(= (:type %) ::promoter)))

(s/def ::gene
  (s/and ::element #(= (:type %) ::gene)))

(s/def ::regulatory-unit
  (s/cat :promoters (s/+ ::promoter)
         :genes (s/+ ::gene)))

(s/def ::parameters
  (s/keys))

(s/def ::genome
  (s/cat :signals (s/spec (s/+ ::signal))
         :internal (s/spec (s/+ (s/spec ::regulatory-unit)))
         :outputs (s/coll-of nat-int?)
         :params ::parameters))

(comment
  (s/exercise ::gene)
  (gen/generate (s/gen ::genome)))

(defn random-element
  [type]
  [type
   (rand-nth [-1 1])
   [(rand 5.0) (rand 5.0)]])

(defn random-unit
  []
  (let [npr (+ 1 (rand-int 3))
        ntf (+ 1 (rand-int 3))]
    (concat (repeatedly npr #(random-element ::promoter))
            (repeatedly ntf #(random-element ::gene)))))

(defn random-genome
  [n-in n-out]
  (let [n-units (+ n-out (rand-int 3))]
    [(repeatedly n-in #(random-element ::signal))
     (repeatedly n-units #(random-unit))
     (repeatedly n-out #(rand-int 100))
     {}]))

;;; phenotype / usable form

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
  (let [pr-coords (:coords pr)]
    (reduce (fn [m [tf-id tf]]
              (if-let [aff (affinity pr-coords (:coords tf))]
                (assoc m tf-id (* aff (:sign pr) (:sign tf)))
                m))
            {}
            (map-indexed vector tfs))))

(defn genome->cell
  [genome]
  (let [cg (s/conform ::genome genome)
        units (:internal cg)
        sigs (:signals cg)
        unit-promoters (index-counts (map #(count (:promoters %)) units) 0)
        unit-tfs (index-counts (map #(count (:genes %)) units) (count sigs))
        all-prs (mapcat :promoters units)
        all-tfs (concat sigs (mapcat :genes units))
        ;; pr -> {tf -> infl}
        pr-influences (mapv #(promoter-influences % all-tfs) all-prs)
        ]
    {::unit-promoters unit-promoters
     ::unit-tfs unit-tfs
     ::influences pr-influences
     ::concs (vec (repeat (count all-tfs) 0.0))
     ::outputs (:outputs cg)
     ::parameters (:params cg)}))

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

(defn merge-vectors
  [v v2]
  (persistent!
   (reduce-kv assoc! (transient v) v2)))

(defn step
  [cell signal-concs dt]
  (let [dt (double dt)
        concs (merge-vectors (::concs cell) signal-concs)
        influences (::influences cell)
        promoters (::unit-promoters cell)
        unit-tfs (::unit-tfs cell)
        n-units (count unit-tfs)]
    (loop [unit-i 0
           new-concs (transient concs)]
      (if (< unit-i n-units)
        (let [pr-ids (get promoters unit-i)
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
          (::outputs cell))))

;; TODO: mutation etc
;; TODO: rng

(comment
  (def g (gen/generate (s/gen ::genome)))
  (def g (random-genome 2 2))
  (def cell (genome->cell g))
  (def csteps (iterate #(step % [1.0 1.0] 0.05) cell))
  (map println (map cell-outputs (take 10 csteps)))
)

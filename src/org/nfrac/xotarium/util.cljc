(ns org.nfrac.xotarium.util
  (:require [clojure.test.check.random :as random]
            [clojure.spec :as s]
            [#?(:clj clojure.spec.gen :cljs clojure.spec.impl.gen) :as gen]
            [clojure.set :as set])
  (:refer-clojure :exclude [rand rand-int rand-nth shuffle]))

;; copied from
;; https://github.com/Datomic/simulant/blob/d681b2375c3e0ea13a0df3caffeb7b3d8a20c6a3/src/simulant/util.clj#L24-L37

(defn getx
  "Like two-argument get, but throws an exception if the key is
   not found."
  [m k]
  (let [e (get m k ::sentinel)]
    (if-not (= e ::sentinel)
      e
      (throw (ex-info "Missing required key" {:map m :key k})))))

(defn getx-in
  "Like two-argument get-in, but throws an exception if the key is
   not found."
  [m ks]
  (reduce getx m ks))

(defn abs
  [x]
  (if (neg? x) (- x) x))

(defn round
  ([x]
   (Math/round (double x)))
  ([x n]
   (let [z (Math/pow 10.0 n)]
     (-> x
         (* z)
         (round)
         (/ z)
         (double)))))

(defn mean
  [xs]
  (/ (apply + xs) (double (count xs))))

(s/def ::rng (-> #(satisfies? random/IRandom %)
                 (s/with-gen #(gen/fmap random/make-random (gen/int)))))

(defn rand
  (^double [rng ^double upper]
   (-> (random/rand-double rng)
       (* upper)))
  (^double [rng ^double lower ^double upper]
   {:pre [(<= lower upper)]}
   (-> (random/rand-double rng)
       (* (- upper lower))
       (+ lower))))

(defn rand-int
  "Uniform integer between lower (inclusive) and upper (exclusive)."
  (^long [rng ^long upper]
   (-> (random/rand-double rng)
       (* upper)
       (Math/floor)
       (long)))
  (^long [rng ^long lower ^long upper]
   (-> (random/rand-double rng)
       (* (- upper lower))
       (+ lower)
       (Math/floor)
       (long))))

(defn rand-nth
  [rng xs]
  (nth xs (rand-int rng (count xs))))

;; copied from
;; https://github.com/clojure/data.generators/blob/bf2eb5288fb59045041aec01628a7f53104d84ca/src/main/clojure/clojure/data/generators.clj
;; adapted to splittable RNG

(defn ^:private fisher-yates
  "http://en.wikipedia.org/wiki/Fisher–Yates_shuffle#The_modern_algorithm"
  [rng coll]
  (let [as (object-array coll)]
    (loop [i (dec (count as))
           r rng]
      (if (<= 1 i)
        (let [[r1 r2] (random/split r)
              j (rand-int r1 (inc i))
              t (aget as i)]
          (aset as i (aget as j))
          (aset as j t)
          (recur (dec i) r2))
        (into (empty coll) (seq as))))))

(defn shuffle
  [rng coll]
  (fisher-yates rng coll))

;; copied from
;; https://github.com/clojure/data.generators/blob/bf2eb5288fb59045041aec01628a7f53104d84ca/src/main/clojure/clojure/data/generators.clj
;; adapted to splittable RNG

(defn reservoir-sample
  "Reservoir sample ct items from coll."
  [rng ct coll]
  (loop [result (transient (vec (take ct coll)))
         n ct
         coll (drop ct coll)
         r rng]
    (if (seq coll)
      (let [[r1 r2] (random/split r)
            pos (rand-int r1 n)]
        (recur (if (< pos ct)
                 (assoc! result pos (first coll))
                 result)
               (inc n)
               (rest coll)
               r2))
      (persistent! result))))

(defn sample
  "Sample ct items with replacement (i.e. possibly with duplicates) from coll."
  [rng ct coll]
  (when (pos? ct)
    (->> (random/split-n rng ct)
         (mapv #(rand-nth % coll)))))

(defn quantile
  [xs p]
  (nth (sort xs) (long (* p (dec (count xs))))))

(defn count-filter
  "Same as `(count (filter pred coll))`, but faster."
  [pred coll]
  (reduce (fn [sum x]
            (if (pred x) (inc sum) sum))
          0 coll))

(defn update-each!
  "Transforms a transient map or vector `m` applying function `f` to
  the values under keys `ks`."
  [m ks f]
  (if (empty? ks)
    m
    (reduce (fn [m k]
              (assoc! m k (f (get m k))))
            m
            ks)))

(defn update-each
  "Transforms a map or vector `m` applying function `f` to the values
  under keys `ks`."
  [m ks f]
  (if (empty? ks)
    m
    (persistent!
     (update-each! (transient m) ks f))))

(defn- mapish? [m] (or (nil? m) (map? m)))

(defn deep-merge-with
  "Like merge-with, but merges maps recursively, applying the given fn
  only when there's a non-map at a particular level."
  [f & maps]
  (apply
   (fn m [& maps]
     (if (every? mapish? maps)
       (apply merge-with m maps)
       (apply f maps)))
   maps))

(defn deep-merge
  "Like merge, but merges maps recursively."
  [& maps]
  (if (every? mapish? maps)
    (apply merge-with deep-merge maps)
    (last maps)))

(defn remap
  "Transforms a map `m` applying function `f` to each value."
  [f m]
  (into (or (empty m) {})
        (map (fn [[k v]] [k (f v)]))
        m))

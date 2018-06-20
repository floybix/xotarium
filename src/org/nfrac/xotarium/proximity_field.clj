(ns org.nfrac.xotarium.proximity-field
  (:require [org.nfrac.liquidfun.testbed :as bed]
            [org.nfrac.liquidfun.core :as lf :refer [body! joint!
                                                     particle-system!]]
            [org.nfrac.liquidfun.vec2d :as v2d]
            [org.nfrac.xotarium.util :as util]
            [quil.core :as quil :include-macros true]
            [quil.middleware]
            [clojure.spec.alpha :as s]
            [clojure.test.check.random :as random])
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

(defn statics-field
  [^liquidfun$b2World world [x-lo x-hi] [y-lo y-hi] stride]
  (let [;aabb (lf/aabb [x-lo x-hi] [y-lo y-hi])
        ;statics (->> (lf/query-aabb aabb)
        ;             (map lf/body-of)
        ;             (distinct)
        ;             (filter #(= :static (lf/body-type %))))
        nx (Math/ceil (/ (- x-hi x-lo) stride))
        ny (Math/ceil (/ (- y-hi y-lo) stride))
        icoords (for [xi (range nx)
                      yi (range ny)]
                  [xi yi])]
    (reduce (fn [df [xi yi :as icoord]]
              (let [x (+ x-lo (* stride xi))
                    y (+ y-lo (* stride yi))
                    hits (lf/query-at-point world [x y])]
                (if (and hits (some #(= :static (lf/body-type (lf/body-of %)))
                                    hits))
                  (assoc df icoord 0)
                  df)))
            {}
            icoords)))

(defn proximity-field
  [^liquidfun$b2World world [x-lo x-hi] [y-lo y-hi] stride max-dist]
  (let [nx (Math/ceil (/ (- x-hi x-lo) stride))
        ny (Math/ceil (/ (- y-hi y-lo) stride))
        icoords (for [xi (range nx)
                      yi (range ny)]
                  [xi yi])
        nbs (reduce (fn [m [xi yi :as ic]]
                      (assoc m ic
                             (cond-> ()
                               (> xi 0) (conj [(dec xi) yi])
                               (> yi 0) (conj [xi (dec yi)])
                               (< xi (dec nx)) (conj [(inc xi) yi])
                               (< yi (dec ny)) (conj [xi (inc yi)]))))
                    {}
                    icoords)
        sf (statics-field world [x-lo x-hi] [y-lo y-hi] stride)]
    (loop [df sf
           dist stride
           last-coords (keys sf)]
      (let [new-nbs (->> last-coords
                         (mapcat nbs)
                         (remove df)
                         (distinct))]
        (if (and (seq new-nbs)
                 (< dist max-dist))
          (recur (merge df (zipmap new-nbs (repeat dist)))
                 (+ dist stride)
                 new-nbs)
          ;; done; convert distance to a proximity score 0-1
          (let [pf (reduce-kv (fn [m ic dist]
                                (assoc m ic (/ (max 0 (- max-dist dist))
                                               max-dist)))
                              {}
                              df)]
            {:pf pf
             :x-lo x-lo
             :y-lo y-lo
             :stride stride}))))))

(defn proximity-score
  [prox-field x y]
  (let [{:keys [pf x-lo y-lo stride]} prox-field
        ix (int (/ (- x x-lo) stride))
        iy (int (/ (- y y-lo) stride))]
    (get pf [ix iy] 0.0)))

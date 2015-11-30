(ns runpmc.core
  (:refer-clojure :exclude [rand])
  (:require [clojure.data.int-map :as i])
  (:import [com.carrotsearch.hppc XorShift128P]
           [net.openhft.hashing LongHashFunction]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(def ^XorShift128P rnd (XorShift128P. 666))

(defn georand
  "Returns a geometrically distributed random number by determining
  the position of the leftmost 1-bit in a uniformly distributed number `w`."
  [^long w]
  (let [res (bit-xor (Long/numberOfLeadingZeros (.nextLong rnd)) 0)]
    (if (>= res w)
      (dec w)
      res)))

(defn rand
  "Random number between 0 and `m` inclusive."
  [m]
  (mod (.nextLong rnd) m))

(defrecord Sketch [l m w int-set n p])

(defn create
  ([^long estimated-max]
   (create (* estimated-max 32) 256 32))
  ([l m w]
   (Sketch. l m w (i/int-set) 0 0.0)))

(defn fill-rate [sketch]
  (* ^long (:p sketch) 100))

(defn find-pos [sketch ^bytes f i j]
  (let [h (.hashBytes (LongHashFunction/farmNa i j) f)]
    (mod h (:l sketch))))

(defn increment [sketch ^bytes f]
  (let [i (rand (:m sketch))
        j (georand (:w sketch))
        sketch (update sketch :n inc)]
    (update sketch :int-set conj (find-pos sketch f i j))))

(defn z-sum [sketch ^bytes f]
  (reduce + (map (fn [i]
                   (first (keep-indexed #(when-not (contains? (:int-set sketch) (find-pos sketch f i %2)) %1)
                                        (range (:w sketch)))))
                 (range (:m sketch)))))

(defn empty-rows [sketch ^bytes f]
  (reduce + (keep #(when-not (contains? (:int-set sketch) (find-pos sketch f % 0)) 1) (range (:m sketch)))))

(defn p [sketch]
  (/ (* 1.0 (count (:int-set sketch))) ^long (:l sketch)))

(defn qk [k ^long n ^long p]
  (reduce * (for [^long i (range k)]
              (* (- 1.0 (Math/pow (- 1.0 (Math/pow 2 (- (inc i)))) n)) (- 1.0 p)))))

(defn e [sketch]
  (let [{:keys [n p]} sketch]
    (reduce + (for [x (range (:w sketch))
                    :let [k (inc ^long x)]]
                (* k (- (qk k n p) (qk (inc k) n p)))))))

(defn phi [sketch]
  (/ (Math/pow 2 (e sketch)) (:n sketch)))

(defn estimate [sketch ^bytes f]
  (let [sketch (if (zero? (:p sketch))
                 (assoc sketch :p (p sketch))
                 sketch)
        k (empty-rows sketch f)
        kp (/ k (:p sketch))
        e (if (> kp (* 0.3 (:m sketch)))
            (* -2 (:m sketch) (Math/log (/ kp (:m sketch))))
            (let [z (z-sum sketch f)]
              (/ (* (:m sketch) (Math/pow 2 (/ z (:m sketch)))) (phi sketch))))]
    (Math/abs e)))

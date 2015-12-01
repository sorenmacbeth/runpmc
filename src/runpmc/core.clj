(ns runpmc.core
  (:refer-clojure :exclude [rand])
  (:require [clojure.data.int-map :as i])
  (:import [com.carrotsearch.hppc XorShift128P]
           [net.openhft.hashing LongHashFunction]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* false)

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
  "Create a new PMC sketch."
  ([^long estimated-max]
   (create (* estimated-max 32) 256 32))
  ([l m w]
   (Sketch. l m w (i/int-set) 0 0.0)))

(defn find-pos
  "Given item `f` and coordinates `i` and `j` in the virtual matrix
  return the position of `f` in the bit field."
  [sketch ^bytes f i j]
  (let [h (.hashBytes (LongHashFunction/farmNa i j) f)]
    (mod h (:l sketch))))

(defn increment
  "Increment the count for item `f` in `sketch`."
  [sketch ^bytes f]
  (let [i (rand (:m sketch))
        j (georand (:w sketch))]
    (-> (update sketch :n inc)
        (update :int-set conj (find-pos sketch f i j)))))

(defn fp-probability [sketch]
  (/ (* 1.0 (count (:int-set sketch))) (:l sketch)))

(defn phi [sketch p]
  (let [qk (fn [k ^long n ^long p]
             (reduce * (for [i (range k)]
                         (* (- 1.0 (Math/pow (- 1.0 (Math/pow 2 (- (inc i)))) n)) (- 1.0 p)))))
        e (fn [sketch p]
            (let [n (:n sketch)]
              (reduce + (for [x (range (:w sketch))
                              :let [k (inc x)]]
                          (* k (- (qk k n p) (qk (inc k) n p)))))))]
    (/ (Math/pow 2 (e sketch p)) (:n sketch))))

(defn zsum [sketch ^bytes f]
  (reduce + (map (fn [i]
                   (first (keep-indexed #(when-not (contains? (:int-set sketch) (find-pos sketch f i %2)) %1)
                                        (range (:w sketch)))))
                 (range (:m sketch)))))

(defn empty-rows
  "Returns the count of empty rows in the virtual matrix for item `f`."
  [sketch ^bytes f]
  (reduce + (keep #(when-not (contains? (:int-set sketch) (find-pos sketch f % 0)) 1) (range (:m sketch)))))

(defn estimate
  "Return the current estimated count for item `f` in `sketch`."
  [sketch ^bytes f]
  (let [p (fp-probability sketch)
        k (empty-rows sketch f)
        kp (/ k (- 1.0 p))
        e (if (> kp (* 0.3 (:m sketch)))
            (* -2 (:m sketch) (Math/log (/ kp (:m sketch)))) ;; modified hitcounting
            (let [z (zsum sketch f)]
              (/ (* (:m sketch) (Math/pow 2 (/ z (:m sketch)))) (phi sketch p))))]
    (Math/abs e)))

(defn fill-rate
  "Returns the fill rate for `sketch`."
  [sketch]
  (* (fp-probability sketch) 100))

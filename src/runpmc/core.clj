(ns runpmc.core
  (:refer-clojure :exclude [rand merge])
  (:require [clojure.data.int-map :as i])
  (:import [com.carrotsearch.hppc XorShift128P]
           [net.openhft.hashing LongHashFunction]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* false)

(def ^XorShift128P rnd (XorShift128P. 666))

(defn georand
  "Returns a geometrically distributed random number by determining
  the position of the leftmost 1-bit in a uniformly distributed number `w`."
  [w]
  (let [res (bit-xor (Long/numberOfLeadingZeros (.nextLong rnd)) 0)]
    (if (>= res w)
      (dec w)
      res)))

(defn rand
  "Random number between 0 and `m` inclusive."
  [m]
  (mod (.nextLong rnd) m))

(defrecord Sketch [bitfield-length vm-rows vm-cols bitfield additions])

(defn create
  "Create a new PMC sketch."
  ([estimated-max]
   (create (* estimated-max 32) 256 32))
  ([bitfield-length vm-rows vm-cols]
   (Sketch. bitfield-length vm-rows vm-cols (i/int-set) 0)))

(defn find-pos
  "Given item `f` and coordinates `i` and `j` in the virtual matrix
  return the position of `f` in the bit field."
  [sketch ^bytes f i j]
  (let [h (.hashBytes (LongHashFunction/farmNa i j) f)]
    (mod h (:bitfield-length sketch))))

(defn increment
  "Increment the count for item `f` in `sketch`."
  [sketch ^bytes f]
  (let [i (rand (:vm-rows sketch))
        j (georand (:vm-cols sketch))]
    (-> (update sketch :additions inc)
        (update :bitfield conj (find-pos sketch f i j)))))

(defn false-positive-bit-prob
  "Returns the probability of false positive bits."
  [sketch]
  (/ (* 1.0 (count (:bitfield sketch))) (:bitfield-length sketch)))

(defn phi [sketch p]
  (let [qk (fn [k n p]
             (reduce * (for [i (range k)]
                         (* (- 1.0 (Math/pow (- 1.0 (Math/pow 2 (- (inc i)))) n)) (- 1.0 p)))))
        e (fn [sketch p]
            (let [n (:additions sketch)]
              (reduce + (for [x (range (:vm-cols sketch))
                              :let [k (inc x)]]
                          (* k (- (qk k n p) (qk (inc k) n p)))))))]
    (/ (Math/pow 2 (e sketch p)) (:additions sketch))))

(defn zsum [sketch ^bytes f]
  (reduce + (map (fn [i]
                   (first (keep-indexed #(when-not (contains? (:bitfield sketch) (find-pos sketch f i %2)) %1)
                                        (range (:vm-cols sketch)))))
                 (range (:vm-rows sketch)))))

(defn empty-rows
  "Returns the count of empty rows in the virtual matrix for item `f`."
  [sketch ^bytes f]
  (reduce + (keep #(when-not (contains? (:bitfield sketch) (find-pos sketch f % 0)) 1) (range (:vm-rows sketch)))))

(defn estimate
  "Return the current estimated count for item `f` in `sketch`."
  [sketch ^bytes f]
  (let [p (false-positive-bit-prob sketch)
        k (empty-rows sketch f)
        kp (/ k (- 1.0 p))
        e (if (> kp (* 0.3 (:vm-rows sketch)))
            (* -2 (:vm-rows sketch) (Math/log (/ kp (:vm-rows sketch)))) ;; modified hitcounting
            (let [z (zsum sketch f)]
              (/ (* (:vm-rows sketch) (Math/pow 2 (/ z (:vm-rows sketch)))) (phi sketch p))))]
    (Math/abs e)))

(defn fill-rate
  "Returns the fill rate for `sketch`."
  [sketch]
  (* (false-positive-bit-prob sketch) 100))

(defn- merge* [sketch1 sketch2]
  (when (apply not= (map (juxt :bitfield-length :vm-rows :vm-cols) [sketch1 sketch2]))
    (throw (Exception. "Sketch options must match for merging.")))
  (assoc sketch1
         :bitfield (i/union (:bitfield sketch1) (:bitfield sketch2))
         :additions (+ (:additions sketch1) (:additions sketch2))))

(defn merge [sketch & more]
  (reduce merge* sketch more))

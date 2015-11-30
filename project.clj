(defproject runpmc "0.1.0-SNAPSHOT"
  :description "Probabilistic Multiplicity Counting"
  :url "https://github.com/sorenmacbeth/runpmc"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.7.0"]
                 [net.openhft/zero-allocation-hashing "0.4"]
                 [com.carrotsearch/hppc "0.7.1"]
                 [org.clojure/data.int-map "0.2.1"]]
  :profiles {:dev {:dependencies [[criterium "0.4.3"]]}}
  :jvm-opts ^:replace ["-Xmx2g"])

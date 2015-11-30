# runpmc

An implementation of [Probabilistic Multiplicity Counting](https://wwwcn.cs.uni-duesseldorf.de/publications/publications/library/Lieven2010a.pdf).

Largely a port of [this go version](https://github.com/seiflotfy/pmc)

## Usage

``` clojure
(def pmc-sketch (loop [i 1000000 s (create 1000000)]
                  (if (zero? i)
                    s
                    (recur (dec i) (increment s (.getBytes "runpmc" "UTF-8")))))

(estimate pmc-sketch (.getBytes "runpmc" "UTF-8")
```

## License

Copyright © 2015 FIXME

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.

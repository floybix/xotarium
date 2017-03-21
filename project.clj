(defproject org.nfrac/xotarium "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}

  :dependencies [[org.clojure/clojure "1.9.0-alpha13"]
                 [org.clojure/core.async "0.3.442"]
                 [org.clojure/test.check "0.9.0"]
                 [org.nfrac/liquidfun-clj "0.1.0-SNAPSHOT"]
                 [org.nfrac/liquidfun-clj.testbed "0.1.0-SNAPSHOT"]]

  :aot [org.nfrac.xotarium.plant-cave]
  :jvm-opts ^:replace ["-server" "-Dclojure.compiler.direct-linking=true"]

  :profiles {:dev {:dependencies [[gyptis "0.2.2"]]}})

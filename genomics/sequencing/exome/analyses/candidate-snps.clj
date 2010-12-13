(require '[exome core])
(in-ns 'exome.core)

(defn in-any-candidate-gene?
  [snp gene-list]
  (contains? (set gene-list) (:gene snp))) 

(defn overview-causative-snps
  [ind]
  (let [vnnlh (variants-nsssi-novel-lof-highqual ind)
        vnnlh-for-homnonref (filter #(homnonref? % ind) vnnlh)
        genes-w-homnonref (map #(:gene %) vnnlh-for-homnonref)
        g-homnonref genes-w-homnonref
        causative-snps-homnonref
          (reduce conj []
            (map #(assoc {}
                         :group "homnonref"
                         :gene (:gene %)
                         :locus (:locus %)
                         :vartype (:vartype %)
                         :consequence (:consequence %))
                 (filter #(in-any-candidate-gene? % g-homnonref) vnnlh-for-homnonref)))
        vnnlh-for-het (filter #(het? % ind) vnnlh)
        genes-w-hets (map #(:gene %) (filter #(het? % ind) vnnlh))
        m-het (group-by (set genes-w-hets) genes-w-hets)
        c-het (zipmap (keys m-het) (map #(count %) (vals m-het)))
        g-het (keys (filter #(> (val %) 1) c-het))
        causative-snps-comphet
          (reduce conj []
            (map #(assoc {}
                         :group "comphet"
                         :gene (:gene %)
                         :locus (:locus %)
                         :vartype (:vartype %)
                         :consequence (:consequence %))
                 (filter #(in-any-candidate-gene? % g-het) vnnlh-for-het)))]
    (flatten (conj causative-snps-homnonref causative-snps-comphet))))

;; diseaseA
(clojure.contrib.duck-streams/with-out-writer "analyses/candidate-snps-diseaseA.txt"
  (doseq [ind diseaseA-individuals] 
    (println "==========") 
    (println ind)
    (doseq [gene (group-by :gene (overview-causative-snps ind))]
      (println "---")
      (println (first gene))
      (doseq [snp (sort-by :gene (second gene))]
        (println (clojure.string/join "\t" ["" (:group snp) (:vartype snp) (:locus snp) (:consequence snp)]))))))


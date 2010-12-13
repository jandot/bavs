(require '[exome core])
(in-ns 'exome.core)

;; Do this afterwards:
;;   cut -f 2,3 analyses/candidate-genes-control.txt | sort | uniq -c
(defn get-candidate-genes-based-on-snps
  [individual]
  (let [dom-genes (genes-dom-nsssi-novel individual "snp")
        dom-lof-genes (genes-dom-nsssi-novel-lof individual "snp")
        dom-lof-highqual-genes (genes-dom-nsssi-novel-lof-highqual individual "snp")
        homnonref-rec-genes (genes-rec-nsssi-novel-based-on-homnonref individual "snp")
        homnonref-rec-lof-genes (genes-rec-nsssi-novel-lof-based-on-homnonref individual "snp")
        homnonref-rec-lof-highqual-genes (genes-rec-nsssi-novel-lof-highqual-based-on-homnonref individual "snp")
        comphet-rec-genes (genes-rec-nsssi-novel-based-on-comphet individual "snp")
        comphet-rec-lof-genes (genes-rec-nsssi-novel-lof-based-on-comphet individual "snp")
        comphet-rec-lof-highqual-genes (genes-rec-nsssi-novel-lof-highqual-based-on-comphet individual "snp")]
    {individual {"snp-dom-nsssi-novel" dom-genes
                 "snp-dom-nsssi-novel-lof" dom-lof-genes
                 "snp-dom-nsssi-novel-lof-highqual" dom-lof-highqual-genes
                 "snp-homnonref-rec-nsssi-novel" homnonref-rec-genes
                 "snp-homnonref-rec-nsssi-novel-lof" homnonref-rec-lof-genes
                 "snp-homnonref-rec-nsssi-novel-lof-highqual" homnonref-rec-lof-highqual-genes
                 "snp-comphet-rec-nsssi-novel" comphet-rec-genes
                 "snp-comphet-rec-nsssi-novel-lof" comphet-rec-lof-genes
                 "snp-comphet-rec-nsssi-novel-lof-highqual" comphet-rec-lof-highqual-genes}}))

;(defn get-candidate-genes-based-on-indels
;  [individual]
;  (let [dom-genes (genes-dom-nsssi-novel individual "indel")
;        dom-lof-genes (genes-dom-nsssi-novel-lof individual "indel")
;        dom-lof-highqual-genes (genes-dom-nsssi-novel-lof-highqual individual "indel")
;        homnonref-rec-genes (genes-rec-nsssi-novel-based-on-homnonref individual "indel")
;        homnonref-rec-lof-genes (genes-rec-nsssi-novel-lof-based-on-homnonref individual "indel")
;        homnonref-rec-lof-highqual-genes (genes-rec-nsssi-novel-lof-highqual-based-on-homnonref individual "indel")
;        comphet-rec-genes (genes-rec-nsssi-novel-based-on-comphet individual "indel")
;        comphet-rec-lof-genes (genes-rec-nsssi-novel-lof-based-on-comphet individual "indel")
;        comphet-rec-lof-highqual-genes (genes-rec-nsssi-novel-lof-highqual-based-on-comphet individual "indel")]
;    {individual {"indel-dom-nsssi-novel" dom-genes
;                 "indel-dom-nsssi-novel-lof" dom-lof-genes
;                 "indel-dom-nsssi-novel-lof-highqual" dom-lof-highqual-genes
;                 "indel-homnonref-rec-nsssi-novel" homnonref-rec-genes
;                 "indel-homnonref-rec-nsssi-novel-lof" homnonref-rec-lof-genes
;                 "indel-homnonref-rec-nsssi-novel-highqual" homnonref-rec-lof-highqual-genes
;                 "indel-comphet-rec-nsssi-novel" comphet-rec-genes
;                 "indel-comphet-rec-nsssi-novel-lof" comphet-rec-lof-genes
;                 "indel-comphet-rec-nsssi-novel-lof-highqual" comphet-rec-lof-highqual-genes}}))
;
;(defn get-candidate-genes-based-on-any
;  [individual]
;  (let [dom-genes (genes-dom-nsssi-novel individual "")
;        dom-lof-genes (genes-dom-nsssi-novel-lof individual "")
;        dom-lof-highqual-genes (genes-dom-nsssi-novel-lof-highqual individual "")
;        homnonref-rec-genes (genes-rec-nsssi-novel-based-on-homnonref individual "")
;        homnonref-rec-lof-genes (genes-rec-nsssi-novel-lof-based-on-homnonref individual "")
;        homnonref-rec-lof-highqual-genes (genes-rec-nsssi-novel-lof-highqual-based-on-homnonref individual "")
;        comphet-rec-genes (genes-rec-nsssi-novel-based-on-comphet individual "")
;        comphet-rec-lof-genes (genes-rec-nsssi-novel-lof-based-on-comphet individual "")
;        comphet-rec-lof-highqual-genes (genes-rec-nsssi-novel-lof-highqual-based-on-comphet individual "")]
;    {individual {"any-dom-nsssi-novel" dom-genes
;                 "any-dom-nsssi-novel-lof" dom-lof-genes
;                 "any-dom-nsssi-novel-lof-highqual" dom-lof-highqual-genes
;                 "any-homnonref-rec-nsssi-novel" homnonref-rec-genes
;                 "any-homnonref-rec-nsssi-novel-lof" homnonref-rec-lof-genes
;                 "any-homnonref-rec-nsssi-novel-highqual" homnonref-rec-lof-highqual-genes
;                 "any-comphet-rec-nsssi-novel" comphet-rec-genes
;                 "any-comphet-rec-nsssi-novel-lof" comphet-rec-lof-genes
;                 "any-comphet-rec-nsssi-novel-lof-highqual" comphet-rec-lof-highqual-genes}}))

(def partitioned-control-individuals (partition 50 control-individuals))
(doseq [batch partitioned-control-individuals]
  (let [result (pmap #(get-candidate-genes-based-on-snps %) batch)]
    (clojure.contrib.duck-streams/with-out-append-writer "analyses/candidate-genes-control.txt"
      (doseq [individual-part result]
        (let [individual-name (first (keys individual-part))
              individual-data (first (vals individual-part))
              filters (keys (first (vals individual-part)))]
          (doseq [f filters]
            (let [genes (individual-data f)]
              (doseq [g genes]
                (println (clojure.string/join "\t" [individual-name g f]))))))))))
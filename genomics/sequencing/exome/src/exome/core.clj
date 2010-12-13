(ns exome.core
  #^{:doc Code to extract candidate genes from exome database}
  (:use [somnium.congomongo]
        [clojure.contrib.json]
        [clojure.contrib.pprint]
        [core.db]
        [core.constants])
)

(mongo! :db "exome-clj") ;; TODO: replace db-name with automatically read value from config

;; Functions working on variants",
(defn present-in-disease?
  [variant disease]
  (contains? (set (map #(:disease %) (:individuals variant))) disease))

(defn present-in-ind?
  [variant ind]
  (contains? (set (map #(:name %) (:individuals variant))) ind))

(defn ind-data-for-variant
  [variant ind]
  (first (filter #(= ind (:name %)) (:individuals variant))))

(defn highqual?
  [variant ind]
  (condp = (:vartype variant)
    "snp" (<= 70 (:gq (ind-data-for-variant variant ind)))
    "indel" (<= 98 (:gq (ind-data-for-variant variant ind)))))

(defn novel?
  [variant]
  (= "true" (:novel variant)))

(defn nsssi?
  [variant]
  (contains? nsssi-terms (:consequence variant)))

(defn lof?
  [variant]
  (contains? lof-terms (:consequence variant)))

(defn homnonref?
  [variant ind]
  (= "homnonref" (:genotype (ind-data-for-variant variant ind))))

(defn het?
  [variant ind]
  (= "het" (:genotype (ind-data-for-variant variant ind))))

(defn variants-nsssi
  [ind]
  (filter nsssi? (variants-in-ind ind)))

(defn snps-nsssi
  [ind]
  (filter nsssi? (snps-in-ind ind)))

(defn indels-nsssi
  [ind]
  (filter nsssi? (indels-in-ind ind)))

(defn variants-nsssi-novel
  [ind]
  (filter nsssi? (variants-novel ind)))
;  (let [first-filter (filter nsssi? (variants-in-ind ind))]
;    (filter novel? first-filter)))

(defn variants-nsssi-novel-lof
  [ind]
  (filter lof? (variants-novel ind)))
;  (let [first-filter (filter nsssi? (variants-in-ind ind))
;        second-filter (filter novel? first-filter)]
;    (filter lof? second-filter)))

(defn variants-nsssi-novel-lof-highqual
  [ind]
  (filter #(highqual? % ind) (variants-nsssi-novel-lof ind)))

(defn snps-nsssi-novel
  [ind]
  (filter nsssi? (snps-novel ind)))
;  (let [first-filter (filter nsssi? (snps-in-ind ind))]
;    (filter novel? first-filter)))

(defn snps-nsssi-novel-lof
  [ind]
  (filter lof? (snps-novel ind)))
;  (let [first-filter (filter nsssi? (snps-in-ind ind))
;        second-filter (filter novel? first-filter)]
;    (filter lof? second-filter)))

(defn snps-nsssi-novel-lof-highqual
  [ind]
  (filter #(highqual? % ind) (snps-nsssi-novel-lof ind)))

(defn indels-nsssi-novel
  [ind]
  (filter nsssi? (indels-novel ind)))
;  (let [first-filter (filter nsssi? (indels-in-ind ind))]
;    (filter novel? first-filter)))

(defn indels-nsssi-novel-lof
  [ind]
  (filter lof? (indels-novel ind)))
;  (let [first-filter (filter nsssi? (indels-in-ind ind))
;        second-filter (filter novel? first-filter)]
;    (filter lof? second-filter)))

(defn indels-nsssi-novel-lof-highqual
  [ind]
  (filter #(highqual? % ind) (indels-nsssi-novel-lof ind)))
;  (let [allqual-snps (variants-nsssi-novel-lof ind)]
;    (filter #(highqual? % ind) allqual-snps)))

;; Functions returning candidate genes under the dominant model
(defn genes-dom-nsssi
  "Returns all gene names that are candidates under a loose dominant model
  (i.e. have at least one NSSSI SNP)"
  [ind vartype]
  (case vartype
    "snp"   (set (map #(:gene %) (snps-nsssi ind)))
    "indel" (set (map #(:gene %) (indels-nsssi ind)))
            (set (map #(:gene %) (variants-nsssi ind)))))
;  [ind]
;  (set (map #(:gene %) (variants-nsssi ind))))

(defn genes-dom-nsssi-novel
  "Returns all gene names that are candidates under a stricter dominant model
  (i.e. have at least one NSSSI SNP that is novel)"
  [ind vartype]
  (case vartype
    "snp"   (set (map #(:gene %) (snps-nsssi-novel ind)))
    "indel" (set (map #(:gene %) (indels-nsssi-novel ind)))
            (set (map #(:gene %) (variants-nsssi-novel ind)))))
;  (set (map #(:gene %) (variants-nsssi-novel ind))))

(defn genes-dom-nsssi-novel-lof
  "Returns all gene names that are candidates under a stricter dominant model
  (i.e. have at least one NSSSI SNP that is novel and LOF)"
  [ind vartype]
  (case vartype
    "snp"   (set (map #(:gene %) (snps-nsssi-novel-lof ind)))
    "indel" (set (map #(:gene %) (indels-nsssi-novel-lof ind)))
            (set (map #(:gene %) (variants-nsssi-novel-lof ind)))))
;  (set (map #(:gene %) (variants-nsssi-novel-lof ind))))

(defn genes-dom-nsssi-novel-lof-highqual
  "Returns all gene names that are candidates under a stricter dominant model
  (i.e. have at least one NSSSI SNP that is novel and LOF and high gq)"
  [ind vartype]
  (case vartype
    "snp"   (set (map #(:gene %) (snps-nsssi-novel-lof-highqual ind)))
    "indel" (set (map #(:gene %) (indels-nsssi-novel-lof-highqual ind)))
            (set (map #(:gene %) (variants-nsssi-novel-lof-highqual ind)))))
;  (set (map #(:gene %) (variants-nsssi-novel-lof-highqual ind))))

(defn genes-dom
  [ind vartype]
  (genes-dom-nsssi-novel-lof-highqual ind vartype))

;; Functions returning candidate genes under the recessive model
(defn genes-rec-nsssi-based-on-homnonref
  "Returns all gene names that are candidates under a loose recessive model
  based on homnonref SNPs)"
  [ind vartype]
  (case vartype
    "snp"
      (let [vn (snps-nsssi ind)]
        (set (map #(:gene %) (filter #(homnonref? % ind) vn))))
    "indel"
      (let [vn (indels-nsssi ind)]
        (set (map #(:gene %) (filter #(homnonref? % ind) vn))))
      (let [vn (variants-nsssi ind)]
        (set (map #(:gene %) (filter #(homnonref? % ind) vn))))))

(defn genes-rec-nsssi-novel-based-on-homnonref
  [ind vartype]
  (case vartype
    "snp"
      (let [vn (snps-nsssi-novel ind)]
        (set (map #(:gene %) (filter #(homnonref? % ind) vn))))
    "indel"
      (let [vn (indels-nsssi-novel ind)]
        (set (map #(:gene %) (filter #(homnonref? % ind) vn))))
      (let [vn (variants-nsssi-novel ind)]
        (set (map #(:gene %) (filter #(homnonref? % ind) vn))))))
;  (let [vnn (variants-nsssi-novel ind)]
;    (set (map #(:gene %) (filter #(homnonref? % ind) vnn)))))

(defn genes-rec-nsssi-novel-lof-based-on-homnonref
  [ind vartype]
  (case vartype
    "snp"
      (let [vn (snps-nsssi-novel-lof ind)]
        (set (map #(:gene %) (filter #(homnonref? % ind) vn))))
    "indel"
      (let [vn (indels-nsssi-novel-lof ind)]
        (set (map #(:gene %) (filter #(homnonref? % ind) vn))))
      (let [vn (variants-nsssi-novel-lof ind)]
        (set (map #(:gene %) (filter #(homnonref? % ind) vn))))))
;  (let [vnnl (variants-nsssi-novel-lof ind)]
;    (set (map #(:gene %) (filter #(homnonref? % ind) vnnl)))))

(defn genes-rec-nsssi-novel-lof-highqual-based-on-homnonref
  [ind vartype]
  (case vartype
    "snp"
      (let [vn (snps-nsssi-novel-lof-highqual ind)]
        (set (map #(:gene %) (filter #(homnonref? % ind) vn))))
    "indel"
      (let [vn (indels-nsssi-novel-lof-highqual ind)]
        (set (map #(:gene %) (filter #(homnonref? % ind) vn))))
      (let [vn (variants-nsssi-novel-lof-highqual ind)]
        (set (map #(:gene %) (filter #(homnonref? % ind) vn))))))
;  (let [vnnlh (variants-nsssi-novel-lof-highqual ind)]
;    (set (map #(:gene %) (filter #(homnonref? % ind) vnnlh)))))

(defn genes-rec-nsssi-based-on-comphet
  "Returns all gene names that are candidates under a loose recessive model
  based on compound heterozygous SNPs)"
  [ind vartype]
  (case vartype
    "snp"
      (let [vn (snps-nsssi ind)
            genes-w-hets (map #(:gene %) (filter #(het? % ind) vn))
            m (group-by (set genes-w-hets) genes-w-hets)
            counts (zipmap (keys m) (map #(count %) (vals m)))]
        (set (keys (filter #(> (val %) 1) counts))))
    "indel"
      (let [vn (indels-nsssi ind)
            genes-w-hets (map #(:gene %) (filter #(het? % ind) vn))
            m (group-by (set genes-w-hets) genes-w-hets)
            counts (zipmap (keys m) (map #(count %) (vals m)))]
        (set (keys (filter #(> (val %) 1) counts))))
      (let [vn (variants-nsssi ind)
            genes-w-hets (map #(:gene %) (filter #(het? % ind) vn))
            m (group-by (set genes-w-hets) genes-w-hets)
            counts (zipmap (keys m) (map #(count %) (vals m)))]
        (set (keys (filter #(> (val %) 1) counts))))))

(defn genes-rec-nsssi-novel-based-on-comphet
  [ind vartype]
  (case vartype
    "snp"
      (let [vn (snps-nsssi-novel ind)
            genes-w-hets (map #(:gene %) (filter #(het? % ind) vn))
            m (group-by (set genes-w-hets) genes-w-hets)
            counts (zipmap (keys m) (map #(count %) (vals m)))]
        (set (keys (filter #(> (val %) 1) counts))))
    "indel"
      (let [vn (indels-nsssi-novel ind)
            genes-w-hets (map #(:gene %) (filter #(het? % ind) vn))
            m (group-by (set genes-w-hets) genes-w-hets)
            counts (zipmap (keys m) (map #(count %) (vals m)))]
        (set (keys (filter #(> (val %) 1) counts))))
      (let [vn (variants-nsssi-novel ind)
            genes-w-hets (map #(:gene %) (filter #(het? % ind) vn))
            m (group-by (set genes-w-hets) genes-w-hets)
            counts (zipmap (keys m) (map #(count %) (vals m)))]
        (set (keys (filter #(> (val %) 1) counts))))))
;  (let [vnn (variants-nsssi-novel ind)
;        genes-w-hets (map #(:gene %) (filter #(het? % ind) vnn))
;        m (group-by (set genes-w-hets) genes-w-hets)
;        counts (zipmap (keys m) (map #(count %) (vals m)))]
;    (set (keys (filter #(> (val %) 1) counts)))))

(defn genes-rec-nsssi-novel-lof-based-on-comphet
  [ind vartype]
  (case vartype
    "snp"
      (let [vn (snps-nsssi-novel-lof ind)
            genes-w-hets (map #(:gene %) (filter #(het? % ind) vn))
            m (group-by (set genes-w-hets) genes-w-hets)
            counts (zipmap (keys m) (map #(count %) (vals m)))]
        (set (keys (filter #(> (val %) 1) counts))))
    "indel"
      (let [vn (indels-nsssi-novel-lof ind)
            genes-w-hets (map #(:gene %) (filter #(het? % ind) vn))
            m (group-by (set genes-w-hets) genes-w-hets)
            counts (zipmap (keys m) (map #(count %) (vals m)))]
        (set (keys (filter #(> (val %) 1) counts))))
      (let [vn (variants-nsssi-novel-lof ind)
            genes-w-hets (map #(:gene %) (filter #(het? % ind) vn))
            m (group-by (set genes-w-hets) genes-w-hets)
            counts (zipmap (keys m) (map #(count %) (vals m)))]
        (set (keys (filter #(> (val %) 1) counts))))))
;  (let [vnnl (variants-nsssi-novel-lof ind)
;        genes-w-hets (map #(:gene %) (filter #(het? % ind) vnnl))
;        m (group-by (set genes-w-hets) genes-w-hets)
;        counts (zipmap (keys m) (map #(count %) (vals m)))]
;    (set (keys (filter #(> (val %) 1) counts)))))

(defn genes-rec-nsssi-novel-lof-highqual-based-on-comphet
  [ind vartype]
  (case vartype
    "snp"
      (let [vn (snps-nsssi-novel-lof-highqual ind)
            genes-w-hets (map #(:gene %) (filter #(het? % ind) vn))
            m (group-by (set genes-w-hets) genes-w-hets)
            counts (zipmap (keys m) (map #(count %) (vals m)))]
        (set (keys (filter #(> (val %) 1) counts))))
    "indel"
      (let [vn (indels-nsssi-novel-lof-highqual ind)
            genes-w-hets (map #(:gene %) (filter #(het? % ind) vn))
            m (group-by (set genes-w-hets) genes-w-hets)
            counts (zipmap (keys m) (map #(count %) (vals m)))]
        (set (keys (filter #(> (val %) 1) counts))))
      (let [vn (variants-nsssi-novel-lof-highqual ind)
            genes-w-hets (map #(:gene %) (filter #(het? % ind) vn))
            m (group-by (set genes-w-hets) genes-w-hets)
            counts (zipmap (keys m) (map #(count %) (vals m)))]
        (set (keys (filter #(> (val %) 1) counts))))))
;  (let [vnnlh (variants-nsssi-novel-lof-highqual ind)
;        genes-w-hets (map #(:gene %) (filter #(het? % ind) vnnlh))
;        m (group-by (set genes-w-hets) genes-w-hets)
;        counts (zipmap (keys m) (map #(count %) (vals m)))]
;    (set (keys (filter #(> (val %) 1) counts)))))

(defn genes-rec
  [ind vartype]
  (clojure.set/union (genes-rec-nsssi-novel-lof-highqual-based-on-homnonref ind vartype) (genes-rec-nsssi-novel-lof-highqual-based-on-comphet ind vartype)))

;; Print some stuff
(def individual #{"AJ168"})
;
;(pprint (genes-rec-nsssi-novel-based-on-comphet "042"))

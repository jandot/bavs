(ns exome.db
  #^{:doc Database-specific functions}
)

;; Example database document in "variants" collection:
;;{ "_id" : ObjectId("4cc68a038b684e79f600024e"),
;;  "modified" : "Tue Oct 26 2010 08:57:55 GMT+0100 (GMT)",
;;  "chr" : "6",
;;  "gene" : "C6orf114",
;;  "vartype" : "snp",
;;  "consequence" : "NON_SYNONYMOUS_CODING",
;;  "individuals" : [
;;    {
;;        "name" : "500exome_2177",
;;        "alleles" : "G/A",
;;        "genotype" : "het",
;;        "gq" : 99,
;;        "disease" : "control"
;;    },
;;    {
;;        "name" : "500exome_38814",
;;        "alleles" : "G/A",
;;        "genotype" : "het",
;;        "gq" : 99,
;;        "disease" : "control"
;;    },
;;    {
;;        "name" : "AJ107",
;;        "alleles" : "G/A",
;;        "genotype" : "het",
;;        "gq" : 38,
;;        "disease" : "MOPD"
;;    }],
;;  "ref" : "G",
;;  "novel" : "false",
;;  "locus" : "6_13578394",
;;  "pos" : "13578394" }

(defn variants-in-ind
  [ind]
  (fetch :variants :where {:individuals.name ind}))

(defn snps-in-ind
  [ind]
  (fetch :variants :where {:individuals.name ind :vartype "snp"}))
;  (filter #(= "snp" (:vartype %) (variants-in-ind ind))))

(defn indels-in-ind
  [ind]
  (fetch :variants :where {:individuals.name ind :vartype "indel"}))
;  (filter #(= "indel" (:vartype %) (variants-in-ind ind))))

(defn variants-novel
  [ind]
  (fetch :variants :where {:individuals.name ind :novel "true"}))

(defn snps-novel
  [ind]
  (fetch :variants :where {:individuals.name ind :novel "true" :vartype "snp"}))

(defn indels-novel
  [ind]
  (fetch :variants :where {:individuals.name ind :novel "true" :vartype "indel"}))

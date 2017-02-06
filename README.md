# supplOVCA

To install

devtools::install_github("CompOnco/supplOvCa")
library(supplOvCa)

#load data
data(ovca_suppl)

#Examples
table(ovca_suppl$SUBTYPE, ovca_suppl$BRCA1_germline) ;

# Are ovarian cancer enriched in BRCA1 germline mutations
fisher.test(table(ovca_suppl$SUBTYPE, ovca_suppl$BRCA1_somatic ))

spineplot(table(ovca_suppl$TUMORRESIDUALDISEASE,ovca_suppl$BRCA1_germline), col=c("red4", "white"), main="BRCA1 germline by residul disease")

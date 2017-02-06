\name{ovca_suppl}
\alias{ovca_suppl}
\docType{data}
\title{
ovca_suppl
}
\description{
Annotation for TCGA ovarian cancer patient, extracted from parsing Supplementary files from ovarian cancer articles including PMID 2325762, 21720365
}
\usage{data("ovca_suppl")}
\format{
  A data frame with 0 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
PMID 23257362, 21720365
}
\references{
1: Verhaak RG, et al.,. Prognostically relevant gene signatures of high-grade
serous ovarian carcinoma. J Clin Invest. 2013 Jan;123(1):517-25. doi:
10.1172/JCI65833. PubMed PMID: 23257362; PubMed Central PMCID: PMC3533304.

2: Cancer Genome Atlas Research Network. Integrated genomic analyses of ovarian
carcinoma. Nature. 2011 Jun 29;474(7353):609-15. doi: 10.1038/nature10166.
PubMed PMID: 21720365; PubMed Central PMCID: PMC3163504.
}
\examples{
# To regenerate the data
#ovca_suppl<-ovca_suppl:::.TCGAsuppl()

data(ovca_suppl)
table(ovca_suppl$SUBTYPE, ovca_suppl$BRCA1_germline) ;

# Are ovarian cancer enriched in BRCA1 germline mutations
fisher.test(table(ovca_suppl$SUBTYPE, ovca_suppl$BRCA1_somatic ))

spineplot(table(ovca_suppl$TUMORRESIDUALDISEASE,ovca_suppl$BRCA1_germline), col=c("red4", "white"), main="BRCA1 germline by residul disease")
}
\keyword{datasets}

####################
# Script to parse suppl pdf file to extract tables with BRCA mut status
# Azfar Basunia, Aedin Culhane
######################


.supplPMID21720365<-function(pdfURL="http://www.nature.com/nature/journal/v474/n7353/extref/nature10166-s1.pdf"){
  library(pdftools)
  library(reshape2)
  #doi:10.1038/nature10166 SUPPLEMENTARY INFORMATION

  #These files are in /inst/extdata
  ovfile=system.file("extdata", "21720365_suppl.pdf", package = "supplOvCa")
  if (!file.exists(ovfile)) {
    download.file(pdfURL, "21720365_suppl.pdf")
    ovfile= "21720365_suppl.pdf"
  }

  # Read table with BRCA mutation info
  brca_mut_pdf <- suppressMessages(pdf_text(ovfile))

  # BRCA germline and somatic mutations
  bb <- brca_mut_pdf[95:97]
  bb <- strsplit(bb, "\n")
  bb<-unlist(bb)
  # remove leading spaces
  bb<-gsub("^  +", "", bb)

  # Remove
  header= sub(" +", "",grep("doi:", bb, value=TRUE)[1])

  grep("WWW.NATURE.COM", bb,value=TRUE)

  bb<- bb[grep("doi:", bb, invert=TRUE)]  # remove doi line
  bb<- bb[grep("WWW.NATURE.COM", bb, invert=TRUE)]  # remove NATURE.com line

  ind<-grep("Table",bb)
  TableNames<-bb[ind]
  bb<- lapply(list(c(ind[1],ind[2]), c(ind[2]+1, length(bb))), function(x) bb[x[1]:x[2]])
  names(bb) = TableNames

  # BRCA1 germline, somatic
  ind<-grep("BRCA",bb[[1]])
  bb[[1]]
  bb[[1]][ind]

  B1g<-matrix(bb[[1]][6:86], ncol=3, byrow=TRUE)
  BRCA1_germline<-cbind(reshape2::colsplit(B1g[,2], " +",c("CaseID", "Mutation", "NT_Position", "CopyNumberStatus")),
                        reshape2::colsplit(B1g[,3], " +",1:3),
                         reshape2::colsplit(B1g[,1], " ",1:2))

  #BRCA1 Somatic
  B1s<- matrix(bb[[1]][91:122][-c(19:20)], ncol=3, byrow=TRUE)
  B1s[6,]<-B1s[6,c(2,3,1)]
  bb[[1]][93:124][c(19:20)]
  BRCA1_somatic<-cbind(reshape2::colsplit(B1s[,2], " +",c("CaseID", "Mutation", "NT_Position", "CopyNumberStatus")),
                        reshape2::colsplit(B1s[,3], " +",1:3),
                        reshape2::colsplit(B1s[,1], " ",1:2))


  #BRCA2 Germline
  # BRCA1 germline, somatic
  ind<-grep("BRCA",bb[[2]])
  bb[[2]][ind]


  B2g<-matrix(bb[[2]][5:65][-34], ncol=3, byrow=TRUE)
  BRCA2_germline<-cbind(reshape2::colsplit(B2g[,2], " +",c("CaseID", "Mutation", "NT_Position", "CopyNumberStatus")),
                        reshape2::colsplit(B2g[,3], " +",1:3),
                        reshape2::colsplit(B2g[,1], " ",1:2))


   # BRCA2 Somatic
  m2<-grep("2 mutation", bb[[2]])
  bb[[2]][m2]<-sub("\\(2 mutations\\)", "TCGA-13-0885", bb[[2]][m2])
  B2s<-matrix(bb[[2]][70:90][-17], ncol=2, byrow=TRUE)
  BRCA2_somatic<-cbind(reshape2::colsplit(B2s[,2], " +",c("CaseID", "Mutation", "NT_Position", "CopyNumberStatus")), reshape2::colsplit(B2s[,1], " ",1:2))


  brca_mut_list<-list(BRCA1_germline=BRCA1_germline, BRCA1_somatic=BRCA1_somatic, BRCA2_germline=BRCA2_germline, BRCA2_somatic=BRCA2_somatic)

  brca_mut_mat<-reshape2::acast(reshape2::melt(lapply(brca_mut_list, function(x) x$CaseID)), value~L1)
  "BRCA1_BRCA2_mutation"= rownames(brca_mut_mat)
  brca_mut_mat<-data.frame(brca_mut_mat)
  brca_mut_mat[,"BRCA1_BRCA2_mutation"]= rownames(brca_mut_mat)

  #There were 66 samples with BRCA1 or BRCA2 mutation
  rownames(brca_mut_mat)<- brca_mut_mat[,"BRCA1_BRCA2_mutation"]
  return(brca_mut_mat)
}



.supplPMID23257362<-function(){

  ovfile=system.file("extdata", "23257362_S1.xls", package = "supplOvCa")
  if (!file.exists(ovfile)) {
     download.file(xls_file <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3533304/bin/JCI65833sd1.xls", destfile = "23257362_S1.xls")
     ovfile= "23257362_S1.xls"
  }

  ov_sup <- gdata::read.xls(xls = ovfile, sheet = 1, skip = 1, as.is=TRUE) # Read table with subtype info
  ov_sup<- ov_sup[,!apply(ov_sup, 2, function(x) all(is.na(x)))]  #remove emptt cols
  ov_sup<-ov_sup[,!colnames(ov_sup)%in%"X.3"]
  return(ov_sup)
}

.TCGAsuppl<-function(){
  ov_sup<- .supplPMID23257362()
  ov_sup<- ov_sup[grep("^TCGA",ov_sup$DATASET),]
  brca_mut<-.supplPMID21720365()
  brca_mut[brca_mut>0]<-"Mut"
  brca_mut[brca_mut==0]<-"WT"
  brca_mut$ID= rownames(brca_mut)
  ovca_suppl<-merge(ov_sup, brca_mut, by.x="ID", by.y="ID", all=TRUE)
  # filter to TCGA data
  return(ovca_suppl)
}

# To regenerate the data
#ovca_suppl<-.TCGAsuppl()

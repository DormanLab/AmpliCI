require(ShortRead)
require(Biostrings)

######################################################################
# Combine AmpliCI results from multiple samples and generate ASV/sOTU table
#' 
#' This function gives an example how to generate ASV/sOTU table from multiple AmpliCI-outputted FASTA files.
#' Idealy, chimeras have already been removed.
#' 
#' @param foldername (Required).The path to the folder for chimeras-removed AmpliCI-outputted FASTA files
#'  
#' @param p.threshold (Optional). post hoc filtering sequences with diagnostic probability > p.threshold.
#'                                Default = 1. 
#'   
#' @return A list of objects: ASVs, ASV.table, n.samples, n.ASVs
#'        ASVs : DNA sequences of ASV/sOTUs
#'        ASV.table: ASV/sOTU table, a matrix that gives the abundance per sample per ASV/OTUs. 
#'                   rows are ASV/OTUs and columns are samples
#'        n.samples: Number of samples
#'        n.ASVs: Number of total ASV/sOTUs
#'
#' @examples
#' 
#' foldername<-"~/Documents/Research/data/vagina/non_ch/"
#' construct_ASVtable(foldername, p.threshold = 1e-40) 
#' 

construct_ASVtable<-function(foldername,p.threshold = 1){
  filenames<-list.files(foldername,pattern = ".fa")
  n.sample<-length(filenames)

  data<-list()
  k<-0
  pvalue<-list()
  ee<-list()
  size<-list()

  for(i in 1:n.sample){
    k<-k+1
    filepath<-paste(foldername,filenames[i],sep = "")
    data[[k]]<-readDNAStringSet(file = filepath)
    pvalue.i<-NULL
    ee.i<-NULL
    size.i<-NULL
  
    for (j in 1:length(data[[k]])){
      name<-names(data[[k]])[j]
      name2<-strsplit(name,";")[[1]][3]
      pvalue.i[j]<-as.numeric(gsub("DiagP=", "", name2))
      ee.i[j]<-as.numeric(gsub("ee=", "", strsplit(name,";")[[1]][4]))
      size.i[j]<-as.numeric(gsub("size=", "", strsplit(name,";")[[1]][2]))
    }
    pvalue[[k]]<-pvalue.i
    ee[[k]]<-ee.i
    size[[k]]<-size.i
  }

### collect all sequences with pvalue > p.threshold in the final results 
  ref.seq.all<-NULL
  length<-NULL

  for (i in 1:n.sample){
    data.sub<-data[[i]][pvalue[[i]]< p.threshold]
    ref.seq.all<-union(ref.seq.all,data.sub)
    length[i]<-length(data[[i]])
  }
  length(ref.seq.all)

### generate ASV/sOTU table
  count.mat<-matrix(0,length(ref.seq.all),n.sample)
  for (i in 1:157){
    for(j in 1:length[i]){
      indiv<-data[[i]][j]
      idx<-which(ref.seq.all==data[[i]][j])
      if(length(idx))
        count.mat[idx,i]<-size[[i]][j]
    }
  }
  colnames(count.mat)<- filenames

  return(list(ASVs = ref.seq.all, ASV.table = count.mat,n.samples = n.sample, n.ASVs = length((ref.seq.all))))
}


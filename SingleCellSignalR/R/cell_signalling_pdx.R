#' Title
#'
#' @param mix.data 
#' @param cluster.human 
#' @param cluster.mouse 
#' @param c.names.human 
#' @param c.names.mouse 
#' @param tol 
#' @param s.score 
#' @param write 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
cell_signaling_pdx = function(mix.data,cluster.human=NULL,cluster.mouse=NULL,c.names.human=NULL,c.names.mouse=NULL,tol=0,s.score=0.5,write=TRUE,verbose=TRUE){
  
  if (is.null(cluster.human) & is.null(cluster.mouse)){
    cluster = c(rep(1,sum(grepl("human",colnames(mix.data)))),rep(2,sum(grepl("mouse",colnames(mix.data)))))
  } else {
    cluster = c(cluster.human,cluster.mouse) + as.numeric(grepl("mouse",colnames(mix.data)))*max(cluster.human)
  }
  
  c.names = c(c.names.human,c.names.mouse)
  if (is.null(c.names)==TRUE){
    c.names = c(paste("cluster human",1:max(cluster.human)),paste("cluster mouse",1:max(cluster.mouse)))
  }
  
  data.mouse = mix.data[grepl("mouse",rownames(mix.data)),grepl("mouse",colnames(mix.data))]
  rownames(data.mouse) = as.character(strsplit(rownames(data.mouse),split = "_mouse"))
  
  data.human = mix.data[grepl("human",rownames(mix.data)),grepl("human",colnames(mix.data))]
  rownames(data.human) = as.character(strsplit(rownames(data.human),split = "_human"))
  
  z = seq(1,max(cluster))
  lig = unique(LRdb$ligand)
  rec = unique(LRdb$receptor)
  med = sum(mix.data)/(nrow(mix.data)*ncol(mix.data))
  
  mm2Hs = mm2Hs[!is.na(mm2Hs$`Mouse orthology confidence`) & mm2Hs$`Mouse orthology confidence`==1,c(3,5)]
  mm2Hs = subset(mm2Hs,!duplicated(mm2Hs$`Gene name`))
  mm2Hs = subset(mm2Hs,mm2Hs$`Mouse gene name`!="a")
  Hs2mm = mm2Hs[,1]
  mm2Hs = mm2Hs[,2]
  names(mm2Hs) = as.character(Hs2mm)
  names(Hs2mm) = as.character(mm2Hs)
  
  m.names = mm2Hs[rownames(data.mouse)]
  data.mouse = subset(data.mouse,(!is.na(m.names)))
  m.names = m.names[!is.na(m.names)]
  rownames(data.mouse)=as.character(m.names)
  
  data = list()
  
  for (i in sort(unique(cluster))){
    if (i<=max(cluster.human)){
      data[[i]] = data.human[,cluster.human==i]
    } else {
      data[[i]] = data.mouse[,cluster.mouse==(i-max(cluster.human))]
    }
     
  }
  
  ## Paracrine -------------------
  if (int.type=="paracrine"){
    if (verbose==TRUE){
      cat("Paracrine signaling: ",fill=TRUE)
    }
    para = list()
    k=0
    int=NULL
    n.int=NULL
    if (verbose==TRUE){
      cat("Checking for signaling between cell types",fill=TRUE)
    }
    for (i in z){
      if (sum(cluster==i)>1){
        tmp = data[[i]]
        tmp = tmp[rowSums(tmp)>0,]
        if (sum(is.element(lig, rownames(tmp)))>0){
          lig.tmp = rownames(tmp)[is.element(rownames(tmp),lig)]
        } else {lig.tmp=NULL}
        
        final.tmp = LRdb[is.element(LRdb$ligand,lig.tmp),1:2]
        final.tmp = data.frame(final.tmp,as.character(rep("paracrine",sum(is.element(LRdb$ligand,lig.tmp)))),stringsAsFactors = FALSE)
        m.lig = rowSums(tmp[unique(final.tmp[,1]),])/sum(cluster==i)
        names(m.lig) = unique(final.tmp[,1])
        
        if (sum(is.element(rec, rownames(tmp)))>0){
          rec.tmp = rownames(tmp)[is.element(rownames(tmp[apply(tmp,1,function(x) sum(x>0))>tol*ncol(tmp),]),rec)]
        } else {rec.tmp=NULL}
        
        for (j in z[-i]){
          if (sum(cluster==j)>1){
            temp = data[[j]]
            temp = temp[rowSums(temp)>0,]
            if (sum(is.element(rec, rownames(temp)))>0){
              rec.temp = rownames(temp)[is.element(rownames(temp),rec)]
            } else {rec.temp=NULL}
            rec.temp = rec.temp[!is.element(rec.temp,rec.tmp)]
            m.rec = rowSums(data.frame(temp[rec.temp,]))/sum(cluster==j)
            names(m.rec) = rec.temp
            
            final = final.tmp[is.element(final.tmp$receptor,rec.temp),]
            final = cbind(final,score(m.lig[final$ligand],m.rec[final$receptor],med))
            
            colnames(final) = c(c.names[i],c.names[j],"interaction type","LRscore")
            final = final[final[,4]>s.score,]
            final = final[order(final[,4],decreasing = TRUE),]
            
            if (i %in% c((max(cluster.human)+1):max(cluster))){
              final[,1] = Hs2mm[as.character(final[,1])]
            }
            
            if (j %in% c((max(cluster.human)+1):max(cluster))){
              final[,2] = Hs2mm[as.character(final[,2])]
            }
            
            if (nrow(final)>0){
              k=k+1
              para[[k]] = final
              if (verbose==TRUE){
                cat(paste(nrow(final),"interactions from",c.names[i],"to",c.names[j]),fill=TRUE)
              }
              int = c(int,paste(i,"-",j,sep=""))
              n.int = c(n.int,paste(c.names[i],"-",c.names[j],sep=""))
              if (write==TRUE){
                fwrite(data.frame(final),paste("./cell-signaling/LR_interactions_",c.names[i],"-",c.names[j],"-",int.type,".txt",sep=""),sep="\t")
              }
            } else {
              if (verbose==TRUE){
                cat(paste("No significant interaction found from",c.names[i],"to",
                          c.names[j]),fill=TRUE)
              }
            }
          }
        }
      }
    }
    if (k!=0){
      names(para) = n.int
    }
  }
  
  
  ## Returns ---------------------
  if (int.type=="paracrine"){
    return(para)
  }
}
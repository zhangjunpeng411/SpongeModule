library(parallel)
library(reshape2)
library(ggplot2)
library(R.matlab)
library(cluster)
library(ggpubr)
library(WGCNA)

gene_pair=function(geneset,mirwalk,core=10)
{
  iteration=function(g1)
  {
    print(paste(which(geneset==g1),g1,date()))
    allset=colnames(mirwalk);
    set1=allset[mirwalk[g1,]==1];
    value=matrix(0,nrow=1,ncol=length(geneset));
    rownames(value)=g1;
    colnames(value)=geneset;
    for(g2 in geneset)
    {
      set2=allset[mirwalk[g2,]==1];
      x=length(intersect(set1,set2));
      m=length(set2);
      n=length(setdiff(allset,set2));
      k=length(set1);
      pvalue=1-phyper(x-1,m,n,k);
      # hyper_test[g1,g2]=pvalue;
      value[g1,g2]=pvalue;
    }
    return(as.data.frame(value));
  }
  
  cl <- makeCluster(core,outfile="pair.txt")
  clusterExport(cl = cl,varlist = 'mirwalk')
  result=parLapply(cl,as.list(geneset),iteration)
  stopCluster(cl)
  t=do.call(rbind,result);
  return(t);
}
gene_correlation=function(geneset,rna.exp)
{
  exp=rna.exp[geneset,]
  x=t(exp)
  res=corAndPvalue(x, x) 
  result=list(correlation=res$cor,pvalue=res$p)
  return(result)
}
liquid_association=function(geneset,rna.exp,miRNA.exp,mirwalk,core=10)
{
  normalize=function(x)
  {
    return((x-mean(x))/sd(x))
  }
  # Normalize the expression of candidate ceRNA with z-score
  rna.exp=apply(X = rna.exp,MARGIN = 1,FUN = normalize)
  # Get the microRNA intersection in microRNA expression and microRNA-ceRNA interaction
  miRNA=intersect(rownames(miRNA.exp),colnames(mirwalk))
  mirwalk=mirwalk[,miRNA]
  micro.exp=miRNA.exp[miRNA,]
  
  LA=function(g1)
  {
    print(paste(g1,date()))
    la=matrix(NA,nrow = 1,ncol = length(geneset))
    rownames(la)=g1
    colnames(la)=geneset
    for(g2 in geneset)
    {
      share.micro=miRNA[which(mirwalk[g1,]&mirwalk[g2,])]
      if(length(share.micro)>1)
      {
        Z=colSums(micro.exp[share.micro,])
        tZ=qnorm(rank(Z)/(length(Z)+1))
        la[g1,g2]=mean(rna.exp[,g1]*rna.exp[,g2]*tZ)
      }
      else if(length(share.micro==1))
      {
        Z=micro.exp[share.micro,]
        tZ=qnorm(rank(Z)/(length(Z)+1))
        la[g1,g2]=mean(rna.exp[,g1]*rna.exp[,g2]*tZ)
      }
    }
    return(la)
  }
  cluster=makeCluster(core,outfile="log.txt")
  result=parLapply(cl = cluster,X = as.list(geneset),fun = LA)
  result=do.call(what = rbind,args = result)
  stopCluster(cluster)

  return(result)
}

parseResult=function(VC.path,Klist)
{
  cluster=list()
  for(k in Klist)
  {
    matrix=as.matrix(readMat(paste(VC.path,"/VC-",k,".mat",sep = ""))$VC)
    order=apply(X = matrix,MARGIN = 1,FUN = order,decreasing=T)
    cluster=c(cluster,list(order[1,]))
  }
  cluster=do.call(what = rbind,args = cluster)
  rownames(cluster)=Klist
  return(cluster)
}

evaluation=function(cluster,SIM,matlab)
{
  setVariable(matlab,cluster=cluster)
  setVariable(matlab,similar=SIM)
  evaluate(matlab,'[c_index,McClain_Rao,Point_biserial,Modularity]=evaluate(cluster,similar)')
  C_index=getVariable(matlab,'c_index')$c.index
  McClain_Rao=getVariable(matlab,'McClain_Rao')$McClain.Rao
  Point_biserial=getVariable(matlab,'Point_biserial')$Point.biserial
  Modularity=getVariable(matlab,'Modularity')$Modularity

  eva=cbind(C_index=C_index[1,],McClain_Rao=McClain_Rao[1,],Point_biserial=Point_biserial[1,],Modularity=Modularity[1,])
  
  DIS=max(SIM)-SIM
  diag(DIS)=0
  
  silhouette=c()
  for(i in seq(1,dim(cluster)[1]))
  {
    if(all(cluster[i,]==1))
    {
      silhouette=c(silhouette,0)
    }
    else
    {
      silhouette=c(silhouette,summary(silhouette(x = cluster[i,],dmatrix=DIS))$avg.width)
    }
  }
  eva=cbind(eva,silhouette=silhouette)
  rownames(eva)=rownames(cluster)
  return(eva)
}
doevaluation=function(cluster,simlist,simlabel,matlab)
{
  evalist=list()
  for(sim in simlist)
  {
    eva=evaluation(cluster,sim,matlab)
    evalist=c(evalist,list(eva))
  }
  names(evalist)=simlabel
  
  melt.evaluate=melt(evalist)
  colnames(melt.evaluate)[1:2]=c('Var1','Var2')
  plist=c()
  for(para in unique(melt.evaluate$Var2))
  {
    p=ggplot(melt.evaluate[which(melt.evaluate$Var2==para),])+geom_line(mapping = aes(x=Var1,y = value,colour=L1))+
      labs(title = para,x='K',y='LEVEL')+theme(plot.title = element_text(hjust = 0.5)) 
    plist=c(plist,list(p))
  }
  print(ggarrange(plotlist=plist,ncol = 2,nrow = 3))
  return(evalist)
}
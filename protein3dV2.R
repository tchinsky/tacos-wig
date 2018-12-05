library(Biostrings)
library(bio3d)
library(muscle)
library(devtools)

NMA <- function(pdb){
  ##Take the pdb IDs for each PDB and sets them each equal to a variable
  pdb1<-pdb[1]
  pdb2<-pdb[2]
  pdb3<-pdb[3]
  
  #gets the PDB files from online
  pdb1<-read.pdb(pdb1)
  pdb2<-read.pdb(pdb2)
  pdb3<-read.pdb(pdb3)
  
  ##trimming the pdb to chain A and getting rid of gaps
  pdb1.open <-trim.pdb(pdb1,atom.select(pdb1,chain = "A"))
  pdb2.open <-trim.pdb(pdb2,atom.select(pdb2,chain = "A"))
  pdb3.open <-trim.pdb(pdb3,atom.select(pdb3,chain = "A"))
  
  #perform NMA on the pdb files
  modes1<-nma(pdb1.open)
  modes2<-nma(pdb2.open)
  modes3<-nma(pdb3.open)
  
  ##print out the analysis
  print(modes1)
  print(modes2)
  print(modes3)
  
  ##plot the NMA analysis graphically
  plot(modes1,sse=pdb1)
  plot(modes2,sse=pdb2)
  plot(modes3,sse=pdb3)
  
  #show the 3d structures for the proteins via R.py files
  mktrj(modes1,mode=7)

  mktrj(modes2,mode=7)

  mktrj(modes3,mode=7)

  
  ##plot the cross-correlation matrix
  cm1<-dccm(modes1)
  cm2<-dccm(modes2)
  cm3<-dccm(modes3)
  
  plot(cm1,sse=pdb1.open, contour=F, col.regions=bwr.colors(20), at=seq(-1,1,0.1))
  plot(cm2,sse=pdb2.open, contour=F, col.regions=bwr.colors(20), at=seq(-1,1,0.1))
  plot(cm3,sse=pdb3.open, contour=F, col.regions=bwr.colors(20), at=seq(-1,1,0.1))
  
  
}

PCA <- function(pdb){
  ##Download the PDB files
  raw.files<-get.pdb(pdb)
  
  ##Extract and align the chains
  files<-pdbsplit(raw.files,pdb)
  pdbs<-pdbaln(files,web.arg=list(email="tmc4503@g.rit.edu"))
  pdbs
  
  ##Calculate the sequence identity
  pdbs$id <- substr(basename(pdbs$id),1,6)
  seqidentity(pdbs)
  
  ##Calculate RMSD
  #rmsd(pdbs, a.inds=pdbs$id, b.inds=pdbs$id,fit = TRUE)
  
  ##FInd sets of similar strucutres
  
  ##gets the sequences for each pdb
  
  seq1 <-pdbseq(raw.files[1])
  seq2 <-pdbseq(raw.files[2])
  seq3 <-pdbseq(raw.files[3])
  
  #blasts each pdb
  blast1<-blast.pdb(seq1)
  blast2<-blast.pdb(seq2)
  blast3<-blast.pdb(seq3)
  
  ##plot each hit from the blast
  hits1<-plot.blast(blast1,cutoff = 240)
  hits2<-plot.blast(blast2,cutoff = 240)
  hits3<-plot.blast(blast3,cutoff = 240)
  
  ##MSA using muscle
  files.msa1<-get.pdb(hits1, path = "raw_pdbs",split = TRUE)
  files.msa2<-get.pdb(hits2, path = "raw_pdbs",split = TRUE)
  files.msa3<-get.pdb(hits3, path = "raw_pdbs",split = TRUE)
  
  pdbs.msa1<-pdbaln(files.msa1,web.arg=list(email="tmc4503@g.rit.edu"))
  pdbs.msa2<-pdbaln(files.msa2,web.arg=list(email="tmc4503@g.rit.edu"))
  pdbs.msa3<-pdbaln(files.msa3,web.arg=list(email="tmc4503@g.rit.edu"))
  
  
  
  ##Structure Superposition
  core1<-core.find(pdbs.msa1)
  core2<-core.find(pdbs.msa2)
  core3<-core.find(pdbs.msa3)
  
  ##get the column representatives
  col1=rep("black",length(core1$volume))
  col2=rep("black",length(core2$volume))
  col3=rep("black",length(core3$volume))
  
  ##get the color spectrum for the column representatives
  col1[core1$volume<2]="pink";col1[core1$volume<1]="red"
  col2[core2$volume<2]="pink";col2[core2$volume<1]="red"
  col3[core3$volume<2]="pink";col3[core3$volume<1]="red"
  
  ##get the core indeces
  core1.inds <-print(core1,vol=1.0)
  core2.inds <-print(core2,vol=1.0)
  core3.inds <-print(core3,vol=1.0)
  
  ##make xzy data
  xyz1 <-pdbfit(pdbs.msa1,core1.inds)
  xyz2 <-pdbfit(pdbs.msa2,core2.inds)
  xyz3 <-pdbfit(pdbs.msa3,core3.inds)
  
  ##PCA analysis
  pc.xray1<-pca.xyz(xyz1[, gaps.pos$f.inds])
  pc.xray2<-pca.xyz(xyz2[, gaps.pos$f.inds])
  pc.xray3<-pca.xyz(xyz3[, gaps.pos$f.inds])
  
  pc.xray1
  pc.xray2
  pc.xray3
  
  ##plot the pca data for visualization
  plot(pc.xray1,col=annotation[,"color"])
  plot(pc.xray2,col=annotation[,"color"])
  plot(pc.xray3,col=annotation[,"color"])
  
  ##plot bio3d residue graphs
  par(mfrow = c(3, 1), cex = 0.75, mar = c(3, 4, 1, 1))
  plot.bio3d(pc.xray1$au[,1], resno=ref.pdb, sse=ref.pdb, ylab="PC1.1")
  plot.bio3d(pc.xray1$au[,2], resno=ref.pdb, sse=ref.pdb, ylab="PC1.2")
  plot.bio3d(pc.xray1$au[,3], resno=ref.pdb, sse=ref.pdb, ylab="PC1.3")
  
  par(mfrow = c(3, 1), cex = 0.75, mar = c(3, 4, 1, 1))
  plot.bio3d(pc.xray2$au[,1], resno=ref.pdb, sse=ref.pdb, ylab="PC2.1")
  plot.bio3d(pc.xray2$au[,2], resno=ref.pdb, sse=ref.pdb, ylab="PC2.2")
  plot.bio3d(pc.xray2$au[,3], resno=ref.pdb, sse=ref.pdb, ylab="PC2.3")
  
  par(mfrow = c(3, 1), cex = 0.75, mar = c(3, 4, 1, 1))
  plot.bio3d(pc.xray3$au[,1], resno=ref.pdb, sse=ref.pdb, ylab="PC3.1")
  plot.bio3d(pc.xray3$au[,2], resno=ref.pdb, sse=ref.pdb, ylab="PC3.2")
  plot.bio3d(pc.xray3$au[,3], resno=ref.pdb, sse=ref.pdb, ylab="PC3.3")
  
}

pdb<-c("4zlj","5or1","4g6u")
#NMA(pdb)
PCA(pdb)
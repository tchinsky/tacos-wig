library(Biostrings)
library(bio3d)
library(muscle)
library(devtools)

make3D <- function(pdb){
  #get the PDB ID names from the list into separate variables
  pdb1 <- pdb[1]
  pdb2 <- pdb[2]
  pdb3 <- pdb[3]
  
  pdbc <- pdb #concatinate the pdb files for structural analysis
  og.files <- get.pdb(pdbc) #query for the pdb files
  
  files <- pdbsplit(og.files, pdbc) # Extract and align the chains we are interested in
  pdbs <- pdbaln(files, web.arg=list(email="tmc4503@g.rit.edu"))
  
  pdbs$id <- basename.pdb(pdbs$id)  # Calculate sequence identity
  seqidentity(pdbs)
  
  rmsd(pdbs, fit=TRUE) # Calculate RMSD (root mean squared distance)
  
  pc <- pca(pdbfit(pdbs), rm.gaps=TRUE) #plot a PCA (principal component analysis) for the proteins
  plot(pc)
  
  modes <- nma(pdbs)  #a quick NMA (normal mode analysis for the proteins)
  plot(modes, pdbs, spread=TRUE)
  
  #make the 3d alignment of protein structures
  
  # read two G-protein structures
  a <- read.pdb(pdb1)
  b <- read.pdb(pdb2)
  c <- read.pdb(pdb3)
  
  # perform iterative alignment
  aln <- struct.aln(a, b)
  aln2 <- strucvt.aln(b, c)
  aln3 <- struct.aln(c, a)
  
  # store new coordinates of protein B
  b$xyz <- aln$xyz
  c$xyz <- aln2$xyz
  a$ayz <- aln3$xyz
  
  # indices at which the superposition should be based
  a.ind <- atom.select(a, chain="A", elety="CA")
  b.ind <- atom.select(b, chain="A", elety="CA")
  c.ind <- atom.select(c, chain="A", elety="CA")
  
  # perform superposition
  xyzS1 <- fit.xyz(fixed=a$xyz, mobile=b$xyz, fixed.inds=a.ind$xyz, mobile.inds=b.ind$xyz)
  xyzS2 <- fit.xyz(fixed=b$xyz, mobile=c$xyz, fixed.inds=b.ind$xyz, mobile.inds=c.ind$xyz)
  xyzS3 <- fit.xyz(fixed=c$xyz, mobile=a$xyz, fixed.inds=c.ind$xyz, mobile.inds=a.ind$xyz)
  
  #print the superpositions
  xyzS1
  xyzS2
  xyzS3
  
  # write coordinates to file
  write.pdb(b, xyz=xyzS1, file=c(pdb1,"-at-",pdb2,".pdb"))
  write.pdb(c, xyz=xyzS2, file=c(pdb2,"-at-",pdb3,".pdb"))
  write.pdb(a, xyz=xyzS3, file=c(pdb3,"-at-",pdb1,".pdb"))
  
  # read G-protein structure
  bs1 <- binding.site(pdb1)
  bs2 <- binding.site(pdb2)
  bs3 <- binding.site(pdb3)
  
  # residue names of identified binding site
  print(bs$resnames)
  
  # Domain analysis
  gs1  <- geostas(pdb1)
  gs2 <- geostats(pdb2)
  gs3 <- geostats(pdb3)
  
  # Fit all frames to the 'first' domain
  domain1.inds <- gs1$inds[[1]]
  domain2.inds <- gs2$inds[[1]]
  domain3.inds <- gs3$inds[[1]]
  
  xyz1 <- pdbfit(pdb1, inds=domain1.inds)
  xyz2 <- pdbfit(pdb2, inds=domain2.inds)
  xyz3 <- pdbfit(pdb3, inds=domain3.inds)
  
  # write fitted coordinates
  write.pdb(pdb1, xyz=xyz1, chain=gs1$atomgrps, file="1d1d_fit-domain1.pdb")
  write.pdb(pdb2, xyz=xyz2, chain=gs2$atomgrps, file="1d1d_fit-domain2.pdb")
  write.pdb(pdb3, xyz=xyz3, chain=gs3$atomgrps, file="1d1d_fit-domain3.pdb")
  
  # plot geostas results
  plot(gs1, contour=FALSE)
  plot(gs2, contour=FALSE)
  plot(gs3, contour=FALSE)
  
  # Invariant core
  core <- core.find(pdb1)
  core <- core.find(pdb2)
  core <- core.find(pdb3)
  
  # fit to core region
  xyzC1 <- pdbfit(pdb1, inds=core)
  xyzC2 <- pdbfit(pdb2, inds=core)
  xyzC3 <- pdbfit(pdb3, inds=core)
  
  # write fitted coordinates
  write.pdb(pdb1, xyz=xyzC1, file="1d1d_fit-core1.pdb")
  write.pdb(pdb2, xyz=xyzC2, file="1d1d_fit-core2.pdb")
  write.pdb(pdb3, xyz=xyzC3, file="1d1d_fit-core3.pdb")
  
  #calculate the normal modes
  pdb1.open <- trim.pdb(pdb1, atom.select(pdb1, chain="A"))
  pdb2.open <- trim.pdb(pdb2, atom.select(pdb2, chain="A"))
  pdb3.open <- trim.pdb(pdb3, atom.select(pdb3, chain="A"))
  
  modes1 <- nma(pdb1.open)
  modes2 <- nma(pdb2.open)
  modes3 <- nma(pdb3.open)
  
  #make a pdb trajectory
  mktrj(modes1, mode=7)
  mktrj(modes2, mode=7)
  mktrj(modes3, mode=7)
  
  #vector field representation
  pymol(modes1, mode=7)
  pymol(modes2, mode=7)
  pymol(modes3, mode=7)
  
  # Calculate the cross-correlation matrix
  cm1 <- dccm(modes1)
  cm2 <- dccm(modes2)
  cm3 <- dccm(modes3)
  
  # View the correlations in the structure
  pymol(cm1, pdb1.open, type="launch")
  pymol(cm2, pdb2.open, type="launch")
  pymol(cm3, pdb3.open, type="launch")
  
  # Deformation energies
  defe1 <- deformation.nma(modes1)
  defe2 <- deformation.nma(modes2)
  defe3 <- deformation.nma(modes3)
  
  defsums1 <- rowSums(defe1$ei[,1:3])
  defsums2 <- rowSums(defe2$ei[,1:3])
  defsums3 <- rowSums(defe3$ei[,1:3])
  
  # Fluctuations
  flucts1 <- fluct.nma(modes1, mode.inds=seq(7,9))
  flucts2 <- fluct.nma(modes2, mode.inds=seq(7,9))
  flucts3 <- fluct.nma(modes3, mode.inds=seq(7,9))
  
  #prints the fluctuations
  flucts1
  flucts2
  flucts3
  
}
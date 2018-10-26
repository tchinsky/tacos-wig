ORF_prediction <- function(DNA_seq)
{ 
  DNA_file <- readLines(DNA_seq)
  sub_seq <- paste(DNA_file[2:length(DNA_file)],collapse="")


  revs_seq <- rev_str(sub_seq)
  revs_comp_seq <- chartr('AGCT','TGCA',revs_seq)
  
  
  orf_size <- 100
  df1 <- get_ORF_2(sub_seq,orf_size*3,'+')
  df2 <- get_ORF_2(revs_comp_seq,orf_size*3,'-')
  df <- rbind(df1,df2)
  df <- as.data.frame(df)
  print(df,row.names=F) 
}

get_ORF_2 <- function(DNA_seq,orf_size,DNA_strand)
{
  found <- 0
  orfcount <- 0
  strand <- c()
  frame <- c()
  orf_start <- c()
  orf_end <- c()
  start <- c("ATG")
  stop <- c("TAA", "TAG","TGA")
  
  start_codon <- findPattern(stop,DNA_seq,start+3+orf_size,nchar(DNA_seq),3) 
  
  
  while(start_codon != -1 && (start_codon+3+orf_size) <= nchar(DNA_seq))
  {
    stop_codon <- findPattern(stop,DNA_seq,start_codon+3+orf_size,nchar(DNA_seq),3)
    if(stop_codon != -1)
    {
      flag <- stopFinder(start_codon,stop_codon,DNA_seq)
      if(flag == 0)
      {
        orfcount <- orfcount +1
        strand[orfcount] <- DNA_strand
        DNA_frame <- start%%3
        if(DNA_frame == 0)
        {
          DNA_frame <- DNA_frame+3
        }
        frame[orfcount] <- DNA_frame
        orf_start[orfcount] <- start_codon
        orf_end[orfcount] <- stop_codon
      }
    }
    start_codon <- findPattern(start,DNA_seq,start_codon+3,nchar(DNA_seq),1)
  }
  count <- 0
  ORF <- c()
  Strand <- c()
  Frame <- c()
  start1 <- c()
  stop1 <- c()
  length <- c()
  
  for(i in 1:orfcount)
  {
    found <- 0
    for(j in 1:orfcount)
    {
      if ((i != j)&&(orf_start[i] >= orf_start[j])&&(orf_end[i] <= orf_end[j]))
      {
        found <- 1
        break
      }
    }
    
    if(found == 0)
    {
      count <- count+1
      ORF[count] <- count
      Strand[count] <- Strand[i]
      Frame[count] <- Frame[i]
      if(strand[i] == '+')
      {
        start_pos <- orf_start[i]
        end_pos <- orf_end[i]+2
      }
      else
      {
        start_pos <- nchar(DNA_seq) - orf_start[i] + 1
        end_pos <- nchar(DNA_seq) - (orf_end[i]+2)+1
      }
      start1[count] <- start_pos
      stop1[count] <- end_pos
      length[count] <- abs(start_pos-end_pos)+1
    }
  }
  df <- cbind(ORF,Strand,Frame,Start,Stop,Length)
  return(df)
}

findPattern <- function(pattern,DNA_seq,start_loc,stop_loc, increment)
{
  for(i in seq(start_loc,stop_loc-2,by=increment))
  {
    fragment <- substr(DNA_seq,i,i+2)
    
    for(j in 1:length(pattern))
    {
      codon <- pattern[j]
      if(fragment == codon)
      {
        return[i]
      }
    }
  }
  return(-1)
}

stopFinder <- function(start_loc,stop_loc,DNA_seq)
{
  mystop <- c("TAA", "TAG", "TGA")
  found <- 0
  for(i in seq(start_loc,stop_loc,3))
  {
    fragment = substr(DNA_seq,i,i+2)
    for(j in 1:length(mystop))
    {
      codon <- mystop[j]
      if(fragment == codon && i != stop_loc)
      {
        found <- 1
        return(1)
      }
    }
  }
  if(found == 0)
  {
    return(0)
  }
}


rev_str<-function(x)
{
  s_rev <- ""
  for(i in rev(1:nchar(x)))
  {
    letter <- substr(x,i,i)
    s_rev = paste(s_rev, letter, sep="")
  }
  s_rev
}
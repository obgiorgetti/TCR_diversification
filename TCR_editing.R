
# Add the directory with the files in the next line between quotation marks 
  setwd("")

  load("XP.th.xlsx") # Metadata - Table with data from sorted mice Excel table with data from sorted mice - not an excel file!
  load("XP.th.xlsx.2") # Metadata - Table with data from sorted mice Excel table with data from sorted mice - not an excel file!
  load("MM.xp.3a") # Counts tables with the sequences from sort experiments (TCR alpha edited mice)
  load("MM.xp.3b") # Counts tables with the sequences from sort experiments (TCR beta edited mice)
  load("TCR.a.spp.MM.CDR3.amino.acid.sequence") # WT TCR alpha sequences from B6 controls (From Giorgetti et al. 2023)
  load("TCR.b.spp.MM.CDR3.amino.acid.sequence") # WT TCR beta sequences from B6 controls (From Giorgetti et al. 2023)
  library(tidyverse)
  library(pwalign)
  library(igraph)
  library(ggseqlogo)
  library(readxl)
  # About the project: The purpose of the project was to generate and characterize
  # a system where "SMARTA" mice, that have transgenic TCRs with a unique specificity
  # where edited in-vivo during lymphocyte development. TCRs were edited in their CDR3 region,
  # either in their TCR alpha (3a) or Beta chains (3b). 
  
  # Note 1) 3a and 3b are the shortened names of the breeding lines where SMARTA TCR alpha 
  # or TCR beta were edited. In the annotation of this code I sometimes use for brevity the term 
  # "3a" interchangeably for "TCR alpha-edited mice" (or "3b" in the TCR beta edited line).
  
  # Note 2) If any part of this code is not working or you have questions/criticisms, please send an email to
  # orlandogiorgetti@gmail.com. This code is published for bringing transparency to our process of
  # analyzing the data; some functions can be adapted to perform analysis/plots on other data but don't expect
  # it to work without putting some effort in adapting the format of the data you input. If you are 
  # patient we will be happy to help with that.

# Main Figure 2  
  # This code outputs a four panel figure, two for each line (edited TCR alpha- or TCR beta- edited)
  # Source data: MM.xp.3a and MM.xp.3b (dataframes with counts for reads obtained from repertoire sequencing)
  # Uses ZGa.all.tb --> a table of CDR3 nucleotide lengths (top panels)
  #      ZGb.all.tb.3 --> a table of out-of-frame CDR3 nucleotide lengths
  #      ZGa.all.tb.2 --> a distribution where the in-frame value is estimated as the average of the two neighboring columns, then represented as a regression line (blue)
  
  
  {
    ZGa.all.tb = table(factor(nchar(unique(MM.xp.3a$CDR3.nucleotide.sequence)),levels = 20:70))
    ZGb.all.tb = table(factor(nchar(unique(MM.xp.3b$CDR3.nucleotide.sequence)),levels = 20:70))
    ZGa.all.tb.2 = ZGa.all.tb
    ZGa.all.tb.2[as.numeric(names(ZGa.all.tb))%%3==0] = (ZGa.all.tb[as.numeric(names(ZGa.all.tb))%%3==2]+ZGa.all.tb[as.numeric(names(ZGa.all.tb))%%3==1])/2
    ZGb.all.tb.2 = ZGb.all.tb
    ZGb.all.tb.2[as.numeric(names(ZGb.all.tb))%%3==0] = (ZGb.all.tb[as.numeric(names(ZGb.all.tb))%%3==2]+ZGb.all.tb[as.numeric(names(ZGb.all.tb))%%3==1])/2
    ZGa.all.tb = ZGa.all.tb[-1]
    ZGa.all.tb.2 = ZGa.all.tb.2[-1]
    ZGb.all.tb = ZGb.all.tb[-1]
    ZGb.all.tb.2 = ZGb.all.tb.2[-1]
    ZGa.all.tb.3 = ZGa.all.tb.2
    ZGa.all.tb.3[as.numeric(names(ZGa.all.tb.2))%%3==0] = 0
    ZGb.all.tb.3 = ZGb.all.tb.2
    ZGb.all.tb.3[as.numeric(names(ZGb.all.tb.2))%%3==0] = 0
    bp.names = 21:70
    bp.names[bp.names%%3 !=0] = NA
    {
      par(mfrow = c(4,1))
      colores.a = rep("gray",51)
      colores.a[seq_along(colores.a)%%3==1] = "darkgreen"
      barplot(ZGa.all.tb,col = colores.a,names = bp.names,cex.names = 2, cex.axis = 2,cex.lab = 2,xlab = "CDR3a nucleotide length", ylab = "Number of clones")
      mybar = barplot(ZGa.all.tb.3,col = colores.a,names = bp.names,cex.names = 2, cex.axis = 2,cex.lab = 2,xlab = "Out-of-frame CDR3a nucleotide length", ylab = "Number of clones")
      
      lo <- loess(ZGa.all.tb.2~as.numeric(names(ZGa.all.tb.2)),span = 0.25)
      lines(mybar,predict(lo),col = "blue",lwd = 2)
      colores.b = rep("gray",51)
      colores.b[seq_along(colores.b)%%3==1] = "darkgreen"
      
      barplot(ZGb.all.tb,col = colores.b,names = bp.names,cex.names = 2, cex.axis = 2,cex.lab = 2,xlab = "CDR3b nucleotide length", ylab = "Number of clones")
      mybar = barplot(ZGb.all.tb.3, col = colores.b,names = bp.names,cex.names = 2, cex.axis = 2,cex.lab = 2,xlab = "Out-of-frame CDR3b nucleotide length", ylab = "Number of clones")
      lo <- loess(ZGb.all.tb.2~as.numeric(names(ZGb.all.tb.2)),span = 0.25)
      lines(mybar,predict(lo),col = "blue",lwd = 2)
    }
    
  }

### Figure 3 and 4 ###

  # Using the tables with cell counts from FACS data, plot the proportion of cells in the Thymus for each line.
  # plots 2.1, 2.2 and 2.4 are percentages of total cells, 2.3 and 2.5 are cells positive for TCR antibodies (edited mice only).

  # Table formatting for ggplot
  
    XP.th.cell = rbind(pivot_longer(XP.th.xlsx,cols = c(4:7),names_to = "cell")[,c("Line","cell","value")],
                       pivot_longer(XP.th.xlsx.2,cols = c(5:8),names_to = "cell")[,c("Line","cell","value")])
    XP.th.cell$cell = sub("_Th","",XP.th.cell$cell)
    XP.th.cell$cell = factor(XP.th.cell$cell,levels = c("DN","DP","CD4","CD8"))
    XP.th.cell.2 = rbind(pivot_longer(XP.th.xlsx,cols = c(9:12),names_to = "cell")[,c("Line","cell","value")])
    XP.th.cell.2$cell = sub("_Th_V","",XP.th.cell.2$cell)
    XP.th.cell.2$cell = factor(XP.th.cell.2$cell,levels = c("DN","DP","CD4","CD8"))
    
    sc.col = c("#F8766D","#7CAE00","darkcyan","darkorchid4")
  
   {
    # SMARTA
    p2.1 = ggplot(XP.th.cell[grepl("SM",XP.th.cell$Line),]) + 
      geom_boxplot(mapping = aes(x = (cell),y = value,fill= cell),outlier.shape = NA) + 
      geom_jitter(mapping = aes(x = (cell),y = value),size = 2.5) +
      scale_y_continuous(name="Percentage of cells",limits = c(0,100),breaks = seq(0,100,by = 20)) +
      scale_fill_manual(values=sc.col) +
      theme_classic() +
      guides(fill=guide_legend(title="Cell type")) +
      theme(axis.title.x=element_blank(),
            text = element_text(size = 30))
    p2.1 
    # CDR3a
    p2.2 = ggplot(XP.th.cell[grepl("3a",XP.th.cell$Line),]) + 
      geom_boxplot(mapping = aes(x = (cell),y = value,fill= cell),outlier.shape = NA) + 
      geom_jitter(mapping = aes(x = (cell),y = value),size = 2.5) +
      scale_y_continuous(name="Percentage of cells",limits = c(0,100),breaks = seq(0,100,by = 20)) +
      scale_fill_manual(values=sc.col) +
      theme_classic() +
      guides(fill=guide_legend(title="Cell type")) +
      theme(axis.title.x=element_blank(),
            text = element_text(size = 30))
    p2.2
    p2.3 = ggplot(XP.th.cell.2[grepl("3a",XP.th.cell.2$Line),]) + 
      geom_boxplot(mapping = aes(x = (cell),y = value,fill= cell),outlier.shape = NA) + 
      geom_jitter(mapping = aes(x = (cell),y = value),size = 2.5) +
      scale_y_continuous(name="Percentage of cells",limits = c(0,100),breaks = seq(0,100,by = 20)) +
      scale_fill_manual(values=sc.col) +
      theme_classic() +
      guides(fill=guide_legend(title="Cell type")) +
      theme(axis.title.x=element_blank(),
            text = element_text(size = 30))
    p2.3
    # CDR3b
    p2.4 =ggplot(XP.th.cell[grepl("3b",XP.th.cell$Line),]) + 
      geom_boxplot(mapping = aes(x = (cell),y = value,fill= cell),outlier.shape = NA) + 
      geom_jitter(mapping = aes(x = (cell),y = value),size = 2.5) +
      scale_y_continuous(name="Percentage of cells",limits = c(0,100),breaks = seq(0,100,by = 20)) +
      scale_fill_manual(values=sc.col) +
      theme_classic() +
      guides(fill=guide_legend(title="Cell type")) +
      theme(axis.title.x=element_blank(),
            text = element_text(size = 30))
    p2.4
    p2.5 = ggplot(XP.th.cell.2[grepl("3b",XP.th.cell.2$Line),]) + 
      geom_boxplot(mapping = aes(x = (cell),y = value,fill= cell),outlier.shape = NA) + 
      geom_jitter(mapping = aes(x = (cell),y = value),size = 2.5) +
      scale_y_continuous(name="Percentage of cells",limits = c(0,100),breaks = seq(0,100,by = 20)) +
      scale_fill_manual(values=sc.col) +
      theme_classic() +
      guides(fill=guide_legend(title="Cell type")) +
      theme(axis.title.x=element_blank(),
            text = element_text(size = 30))
    p2.5
  }

  # Helper function for plot figures 3e to 3h and 4e to 4h
  # Produces plots of CDR3 length by usage (top - normalized to UMI usage) and number of variants (bottom - 1 count per clonotype i.e. setting all UMI counts equal)
    # df = dataframe (MM.xp.3a or MM.xp.3b)
    # the rest of the options match data in the dataframe selected in df
    # gt = genotype - ZGa_no_GP or ZGb_no_GP
    # tss = tissue - Thymus or Spleen ("T" or "S")
    # exp = xp - experiment number (for pooling experiments, where same condition is analyzed, for example pooling all tetramer sampled)
    # sbn = Sample.bio.name - name of the sample (i.e. mouse number)
    # cll = cell - cell type ("DN","CD4", etc.)
    # cann = the sequence from CDR3.amino.acid.sequence used as canonical (non-edited) for that strain ("CAANQGGRALIF" for CDR3a and "CASSDFGGGQDTQYF" for CDR3b)
    
    MM.plot.mixed.3 = function(df,gt,tss,exp,cll,sbn,cann,titulo = "",subt = "",desde = 24,hasta = 54,cx= 1){
    {
      bp.count = (with(df[with(df, genotype%in%gt & tissue %in% tss & xp %in% exp & cell%in%cll & grepl(sbn,Sample.bio.name)),],table(factor(L[!duplicated(CDR3.nucleotide.sequence)],levels = desde:hasta))))
      print(sum(bp.count))
      bp.1 = bp.count
      bp.list.1 = with(df[with(df, genotype%in%gt & tissue %in% tss & xp %in% exp & cell%in%cll & grepl(sbn,Sample.bio.name)),],tapply(Umi.count,L,function(x)sort(x,decreasing = T)))
      bp.list.2 = as.list(rep(0,sum(!(desde:hasta %in%as.numeric(names(bp.list.1))))))
      names(bp.list.2) = (desde:hasta)[!(desde:hasta %in%as.numeric(names(bp.list.1)))]
      bp.list = c(bp.list.1,bp.list.2)
      bp.list = bp.list[as.numeric(names(bp.list))%in%desde:hasta]
      bp.list = bp.list[order(names(bp.list))]
      umi.sum =sum(unlist(bp.list))
      bp.list = lapply(bp.list,function(x)x/umi.sum)
    }
    mult.factor = max(sapply(bp.list,sum))
    bp.nombres = rep(NA,hasta-desde+1)
    bp.nombres[seq(desde,hasta,1)%%3==0] = seq(desde,hasta,1)[seq(desde,hasta,1)%%3==0]
    bp = barplot(sapply(sapply(bp.list,function(x)x/mult.factor),function(x)c(x,rep(0,max(sapply(bp.list,length))-length(x)))),col = grey(seq(0.6, 1, length = 10)),ylim = c(-1,1),main = titulo, sub = paste(subt,paste(umi.sum,"cDNA molecules"),sep = " - "),ann=FALSE,axes=FALSE,names = bp.nombres,cex.names = cx)
    barplot(table(factor(with(df[with(df, genotype%in%gt & tissue %in% tss & xp %in% exp & cell%in%cll &CDR3.amino.acid.sequence == cann& grepl(sbn,Sample.bio.name)),],rep(L,Umi.count)),levels = desde:hasta))/umi.sum/mult.factor,add = T,col = "red",names = NA,las = 2,axes=FALSE)
    axis(2,at=seq(0,1,0.25),labels =format(round(seq(0,mult.factor,length.out = 5),digits = 2),nsmall = 2) ,las=2,cex.axis =cx)
    barplot(-(bp.1/max(bp.1)),add = T,col = rep(c("dark green","light grey","light grey"),12),ann=FALSE,axes=FALSE,names = NA)
    axis(4,at=seq(0,-1,-0.25),labels =floor(seq(0,max(unlist(bp.1)),length.out = 5)) ,las=2,cex.axis= cx)
    
  } 
  
  # Plots
  # CDR3a mouse (Fig.3)
    xp.7a.mice = c("MM_5826a")
    cell.types = c("DN","DP","CD4","CD8")
  {
    par(mfcol = c(4,1))
    par(mar = c(2,6,2,6))
    lapply(xp.7a.mice,function(y)lapply(1:4,function(x)MM.plot.mixed.3(MM.xp.3a,"ZGa_no_GP","T",7,cell.types[x],y,"CAANQGGRALIF",titulo = cell.types[x],subt = y,desde = 27,hasta=45,cx=1.33)))
  }
  # CDR3b mouse (Fig.4)
  xp.6b.mice = c("MM_5097b","MM_5098b")
  {
    par(mfcol = c(4,1))
    par(mar = c(2,6,2,6))
    lapply(xp.6b.mice[2],function(y)lapply(1:4,function(x)MM.plot.mixed.3(MM.xp.3b,"ZGb_no_GP","T",6,cell.types[x],y,"CASSDFGGGQDTQYF",titulo = cell.types[x],subt = y,desde = 33,hasta = 51,cx=1.33)))
  }
  
  # CDR3 distribution for Tetramer +/- 
  {
      par(mar = c(2,6,2,6))
      par(mfcol = c(2,2))
    # Tetramer positive, TCR alpha edited
      lapply(".",function(y)lapply(3:3,function(x)MM.plot.mixed.3(MM.xp.3a[which(MM.xp.3a$tetramer == "p"),],"ZGa_no_GP","S",3:7,cell.types[x],y,"CAANQGGRALIF",titulo = cell.types[x],subt = "ZGa_no_GP Tetramer (+)")))
    # Tetramer negative, TCR alpha edited
      lapply(".",function(y)lapply(3:3,function(x)MM.plot.mixed.3(MM.xp.3a[which(MM.xp.3a$tetramer == "n"),],"ZGa_no_GP","S",3:6,cell.types[x],y,"CAANQGGRALIF",titulo = cell.types[x],subt = "ZGa_no_GP Tetramer (-)")))
    # Tetramer positive, TCR alpha edited
      lapply(".",function(y)lapply(3:3,function(x)MM.plot.mixed.3(MM.xp.3b[which(MM.xp.3b$tetramer == "p"),],"ZGb_no_GP","S",3:7,cell.types[x],y,"CASSDFGGGQDTQYF",titulo = cell.types[x],subt = "ZGa_no_GP Tetramer (+)")))
    # Tetramer positive, TCR alpha edited
      lapply(".",function(y)lapply(3:3,function(x)MM.plot.mixed.3(MM.xp.3b[which(MM.xp.3b$tetramer == "n"),],"ZGb_no_GP","S",3:6,cell.types[x],y,"CASSDFGGGQDTQYF",titulo = cell.types[x],subt = "ZGa_no_GP Tetramer (-)")))
  }
  
  # Percentage of unedited cells from each line
    
    p3.1 = ggplot(XP.th.xlsx) + geom_boxplot(mapping = aes(x=Line,y=Tet,fill = Line)) +
      geom_jitter(mapping =  aes(x=Line,y=Tet)) +
      scale_y_continuous(name="Percentage of cells",limits = c(0,100),breaks = seq(0,100,by = 20)) +
      scale_fill_manual(values=sc.col) +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            text = element_text(size = 30))
    p3.1
### Figures 6 and 7 ###
    
  # MM.xp.3a.seq: reformating the MM.xp.3a dataframe to obtain UMI molecules counts for each amino acid sequence in spleen, 
  # based on tetramer sort results (only negative and positive intensities, low or intermediate were not included in this analysis):
    {
      MM.xp.3a.seq = MM.xp.3a[with(MM.xp.3a,J.gene == "TRAJ15_01" & tissue == "S" & xp %in% 3:7 & L%%3==0 & grepl("ZGa",genotype)),] %>% 
        pivot_wider(id_cols = CDR3.amino.acid.sequence,
                    names_from = c(genotype,tetramer),values_from = Umi.count,values_fn = sum)
      MM.xp.3a.seq[is.na(MM.xp.3a.seq)] = 0
      MM.xp.3a.seq = MM.xp.3a.seq[order(rowSums(MM.xp.3a.seq[,-1]),decreasing = T),]
     }

    
  # Figure 6a
    
    # net: using the igraph package, make an undirected graph connecting sequences at a distance of 1 for TCR sequences.
    # seqs: list of sequences (amino acid) in MM.xp.3a.seq
    # mat: matrix of distances, where all larger than 1 are set to 0 (only 1 aa distance are connected)
    # size of dots is log10 counts + 1.5
    # majority rule is used for coloring as binding(blue) or nonbinding(red)
    
    # NOTE: the exact position of nodes in the plot is random and will differ from the published one.
    
    {
      seqs = MM.xp.3a.seq$CDR3.amino.acid.sequence
      seqs = seqs[ MM.xp.3a.seq$ZGa_no_GP_n + MM.xp.3a.seq$ZGa_no_GP_p >0]
      seq.counts.tb = MM.xp.3a.seq[ MM.xp.3a.seq$ZGa_no_GP_n + MM.xp.3a.seq$ZGa_no_GP_p >0,c("ZGa_no_GP_n","ZGa_no_GP_p")]
      mat=as.matrix(stringDist(seqs))
      canonical.dist = mat[1,]
      canonical.dist.a = canonical.dist
      canonical.dist[canonical.dist>=9] = 9
      canonical.dist = canonical.dist+1
      mat[mat!=1]=0
      net = graph_from_adjacency_matrix(mat,mode = "undirected")
      colores = rep("white",length(net))
      colores[seq.counts.tb[,1]>seq.counts.tb[,2]] = "red"
      colores[seq.counts.tb[,1]<seq.counts.tb[,2]] = "blue"
      colores[seqs%in%c("CAANQGGRALIF","CASSDFGGGQDTQYF")] = "yellow"
      V(net)$size = log10(rowSums(seq.counts.tb))+1.5
      V(net)$name = NA
      V(net)$color = colores
      l = layout.fruchterman.reingold(net)
      par(mfcol = c(1,1))
      plot.igraph(net,edge.width=E(net)$weight,layout = l)
    }

    # Same for CDR3b line (Figure 7a)
    { 
      MM.xp.3b.seq = MM.xp.3b[with(MM.xp.3b, xp %in% 3:5 & L%%3==0 & tissue == "S" & grepl("ZGb",genotype)),] %>% 
        pivot_wider(id_cols = CDR3.amino.acid.sequence,
                    names_from = c(genotype,tetramer),values_from = Umi.count,values_fn = sum)
      MM.xp.3b.seq[is.na(MM.xp.3b.seq)] = 0
      MM.xp.3b.seq = MM.xp.3b.seq[order(rowSums(MM.xp.3b.seq[,-1]),decreasing = T),]
      seqs.b = MM.xp.3b.seq$CDR3.amino.acid.sequence
      seqs.b = seqs.b[ MM.xp.3b.seq$ZGb_no_GP_n + MM.xp.3b.seq$ZGb_no_GP_p >0]
      seq.counts.tb.b = MM.xp.3b.seq[ MM.xp.3b.seq$ZGb_no_GP_n + MM.xp.3b.seq$ZGb_no_GP_p >0,c("ZGb_no_GP_n","ZGb_no_GP_p")]
      mat.b=as.matrix(stringDist(seqs.b))
      canonical.dist.b = mat.b[1,]
      mat.b[mat.b!=1]=0
      net.b = graph_from_adjacency_matrix(mat.b,mode = "undirected")
      colores = rep("white",length(net.b))
      colores[seq.counts.tb.b[,1]>seq.counts.tb.b[,2]] = "red"
      colores[seq.counts.tb.b[,1]<seq.counts.tb.b[,2]] = "blue"
      colores[seqs.b%in%c("CAANQGGRALIF","CASSDFGGGQDTQYF")] = "yellow"
      V(net.b)$size = log10(rowSums(seq.counts.tb.b))+2
      V(net.b)$name = NA
      V(net.b)$color = colores
      l = layout.fruchterman.reingold(net.b)
      par(mfcol = c(1,1))
      plot.igraph(net.b,edge.width=E(net.b)$weight,layout = l)
    }
    
    # Graph statistics
      # 3a line:
        mean(igraph::degree(net)) # mean connectivity
        mean(igraph::degree(net)[ V(net)$color =='blue']) # mean connectivity of tetramer binding
        mean(igraph::degree(net)[ V(net)$color =='red']) # mean connectivity of tetramer non-binding
      # 3b line:
        mean(igraph::degree(net.b)) # mean connectivity
        mean(igraph::degree(net.b)[ V(net.b)$color =='blue']) # mean connectivity of tetramer binding
        mean(igraph::degree(net.b)[ V(net.b)$color =='red']) # mean connectivity of tetramer non-binding
    
    # subgraphs in supplementary Fig. S5
    
      plot.igraph(subgraph(net,V(net)$color %in% c("yellow","blue")))
      plot.igraph(subgraph(net,V(net)$color %in% c("yellow","red")))
      plot.igraph(subgraph(net.b,V(net.b)$color %in% c("yellow","blue")))
      plot.igraph(subgraph(net.b,V(net.b)$color %in% c("yellow","red")))
    
  # Figures 6b and 7b: barplot of proportion of TCRs at a given distance of canonical sequences (weighted by UMI usage)
    # canonical.dist.a : distance to canonical sequence "CAANQGGRALIF" (first column of matrix "mat")  
    # canonical.dist.b : distance to canonical sequence "CASSDFGGGQDTQYF" (first column of matrix "mat")  
    { 
      par(mfrow = c(2,1))
      barplot(prop.table(table(factor(rep(canonical.dist.a,ceiling(seq.counts.tb$ZGa_no_GP_p/10000)),levels = 0:6))),ylim = c(-1,1),ylab = 'Proportion of TCRs', xlab = "Distance to canonical CDR3a sequence", xaxt = 'n',col = "blue")
      barplot(-prop.table(table(factor(rep(canonical.dist.a,ceiling(seq.counts.tb$ZGa_no_GP_n/10000)),levels = 0:6))),ylim = c(-1,1),ylab = 'Proportion of TCRs', xlab = "Distance to canonical CDR3a sequence", xaxt = 'n',add = T,col = "red")
      legend("bottomright",legend = c("tet(+)","tet(-)"),fill = c("blue","red"),box.lty=0)
      barplot(prop.table(table(factor(rep(canonical.dist.b,ceiling(seq.counts.tb.b$ZGb_no_GP_p/10000)),levels = 0:6))),ylim = c(-1,1),ylab = 'Proportion of TCRs', xlab = "Distance to canonical CDR3b sequence", xaxt = 'n',col = "blue")
      barplot(-prop.table(table(factor(rep(canonical.dist.b,ceiling(seq.counts.tb.b$ZGb_no_GP_n/10000)),levels = 0:6))),ylim = c(-1,1),ylab = 'Proportion of TCRs', xlab = "Distance to canonical CDR3b sequence", xaxt = 'n',add = T,col = "red")
      legend("bottomright",legend = c("tet(+)","tet(-)"),fill = c("blue","red"),box.lty=0)
    }
  # CDF of the same distribution (CDR3a)
    plot(cumsum(prop.table(table(factor(rep(canonical.dist.a,ceiling(seq.counts.tb$ZGa_no_GP_p/10000)),levels = 0:8)))),type = "l",col = 'blue',ylab = 'Proportion of TCRs (CDF)', xlab = "Distance to canonical CDR3a sequence",lty = 2)
    points(cumsum(prop.table(table(factor(rep(canonical.dist.a,ceiling(seq.counts.tb$ZGa_no_GP_p/10000)),levels = 0:8)))),col = 'blue',pch = 16)
    points(cumsum(prop.table(table(factor(rep(canonical.dist.a,ceiling(seq.counts.tb$ZGa_no_GP_n/10000)),levels = 0:8)))),col = 'red',pch = 16)
    lines(cumsum(prop.table(table(factor(rep(canonical.dist.a,ceiling(seq.counts.tb$ZGa_no_GP_n/10000)),levels = 0:8)))),col = 'red',lty = 2)
  # CDF of the same distribution (CDR3b)
    plot(cumsum(prop.table(table(factor(rep(canonical.dist.b,ceiling(seq.counts.tb.b$ZGb_no_GP_p/10000)),levels = 0:8)))),type = "l",col = 'blue',ylab = 'Proportion of TCRs (CDF)', xlab = "Distance to canonical CDR3a sequence",lty = 2)
    lines(cumsum(prop.table(table(factor(rep(canonical.dist.b,ceiling(seq.counts.tb.b$ZGb_no_GP_n/10000)),levels = 0:8)))),col = 'red',lty = 2)
    points(cumsum(prop.table(table(factor(rep(canonical.dist.b,ceiling(seq.counts.tb.b$ZGb_no_GP_n/10000)),levels = 0:8)))),col = 'red',pch = 16)
    points(cumsum(prop.table(table(factor(rep(canonical.dist.b,ceiling(seq.counts.tb.b$ZGb_no_GP_p/10000)),levels = 0:8)))),col = 'blue',pch = 16)

  # This is a helper function for putting matrices in the same format.
  # Will add rows of zeros to the amino acid matrix obtained from
  # consensusMatrix if any amino acid is missing.
    mx.fix.rows.aa.st = function(Profile,pseudocounts = 0,normalizar = 0){
      if('-' %in% rownames(Profile))Profile = Profile[-grep('-',rownames(Profile)),]
      if(mean(c(names(AMINO_ACID_CODE)[c(1:20,26)],"*")%in%rownames(Profile))==1) Profile = Profile[c(names(AMINO_ACID_CODE)[c(1:20,26)],"*"),]
      if(mean(rownames(Profile) %in% c(names(AMINO_ACID_CODE)[c(1:20,26)],"*")) == 1 & length(rownames(Profile)) == 22){
        if(pseudocounts> 0) Profile = Profile+pseudocounts
        if(normalizar >0) Profile = prop.table(Profile)*normalizar
        return (Profile[c(names(AMINO_ACID_CODE)[c(1:20,26)],"*"),])}else{
          if (nrow(Profile)<22){
            Profile = rbind(Profile,matrix(0,ncol = ncol(Profile),nrow = (22-nrow(Profile))))
            rownames(Profile)[!rownames(Profile)%in%c(names(AMINO_ACID_CODE)[c(1:20,26)],"*")] = c(names(AMINO_ACID_CODE)[c(1:20,26)],"*")[!c(names(AMINO_ACID_CODE)[c(1:20,26)],"*")%in%rownames(Profile)]
            Profile = Profile[c(names(AMINO_ACID_CODE)[c(1:20,26)],"*"),]
            if(pseudocounts> 0) Profile = Profile+pseudocounts
            if(normalizar >0) Profile = prop.table(Profile)*normalizar
            return(Profile)
          }
        }}
    
  # Logo plots using ggseqlogo package
  
  {  
  ZGa.nt = with(MM.xp.3a,CDR3.nucleotide.sequence)
  # split in sequences from each line
    ZGa.all = unique(with(MM.xp.3a[with(MM.xp.3a, grepl("ZGa",genotype) & grepl("S",tissue)& L%%3==0 & !grepl("\\*",CDR3.amino.acid.sequence)),],CDR3.amino.acid.sequence))
    ZGb.all = unique(with(MM.xp.3b[with(MM.xp.3b, grepl("ZGb",genotype) & grepl("S",tissue)& L%%3==0 & !grepl("\\*",CDR3.amino.acid.sequence)),],CDR3.amino.acid.sequence))
  # split in tetramer positive and negative
    ZGa.tet.p = (MM.xp.3a[with(MM.xp.3a, grepl("ZGa",genotype) & grepl("p",tetramer) & L%%3==0 & !grepl("\\*",CDR3.amino.acid.sequence)&Umi.count>=3),]$CDR3.amino.acid.sequence)
    ZGa.tet.n = (MM.xp.3a[with(MM.xp.3a, grepl("ZGa",genotype) & grepl("n",tetramer)& L%%3==0 & !grepl("\\*",CDR3.amino.acid.sequence)&Umi.count>=3),]$CDR3.amino.acid.sequence)
    ZGb.tet.p = (MM.xp.3b[with(MM.xp.3b, grepl("ZGb",genotype) & grepl("p",tetramer) & L%%3==0 & !grepl("\\*",CDR3.amino.acid.sequence)&Umi.count>=3),]$CDR3.amino.acid.sequence)
    ZGb.tet.n = (MM.xp.3b[with(MM.xp.3b, grepl("ZGb",genotype) & grepl("n",tetramer) & L%%3==0 & !grepl("\\*",CDR3.amino.acid.sequence)&Umi.count>=3),]$CDR3.amino.acid.sequence)
  # Fig. 6c (bottom): plot 12 aa long CDR3s (3a line) using a matrix displaying the difference in composition betweeen tetramer binding and nonbinding ("cM.1")
    largo =12
    cM.1 = mx.fix.rows.aa.st(consensusMatrix(ZGa.tet.p[nchar(ZGa.tet.p)==largo],as.prob = T))-mx.fix.rows.aa.st(consensusMatrix(ZGa.tet.n[nchar(ZGa.tet.n)==largo],as.prob = T))
    cM.1[is.nan(cM.1)] = 0
    p3.2 = ggseqlogo(cM.1, method='custom', seq_type='aa') +theme(axis.text.x = element_blank())
    p3.2
  # Fig. 7c (bottom): plot 15 aa long CDR3s from 3b line using a matrix displaying the difference in composition betweeen tetramer binding and nonbinding ("cM.1")
    largo =15
    cM.1 = mx.fix.rows.aa.st(consensusMatrix(ZGb.tet.p[nchar(ZGb.tet.p)==largo],as.prob = T))-mx.fix.rows.aa.st(consensusMatrix(ZGb.tet.n[nchar(ZGb.tet.n)==largo],as.prob = T))
    cM.1[is.nan(cM.1)] = 0
    p3.3 = ggseqlogo(cM.1, method='custom', seq_type='aa') +theme(axis.text.x = element_blank())
    p3.3
  # Logo plots for all sequences of same length as respective canonical sequence (Fig. 6c and 7c respectively, top)
    p3.4 = ggseqlogo(ZGa.all[nchar(ZGa.all)==12],method = 'probability')
    p3.4
    p3.5 = ggseqlogo(ZGb.all[nchar(ZGb.all)==15],method = 'probability')
    p3.5
  }

### Comparison with LCMV specific clones from literature ###
    
    tb2 = read_xlsx("jem_20200650_tables2.xlsx")
    # Obtain CDR3 alpha and beta chains from Khatun, A. et al. JEM 2021. 218: e20200650
    
    tb2.A.nt = (DNAStringSet(sub("_TRB.*","",subseq(as.data.frame(tb2[,2])[-1,],5))))
    
    tb2.B.nt = (DNAStringSet(sub(".*_TRB.","",subseq(as.data.frame(tb2[,2])[-1,],5))))
    # Overlap between LCMV clones (Khatun et al.) and WT repertoire (Giorgetti et al. 2023, amino acid CDR3)
      unique(tb2.A.nt[translate(tb2.A.nt)%in%TCR.a.spp.MM.CDR3.amino.acid.sequence]) # 78.6%
    # Overlap between LCMV clones and TCR alpha edited mice is null
      tb2.A.nt[translate(tb2.A.nt)%in%MM.xp.3a$CDR3.amino.acid.sequence]
      which(translate(tb2.B.nt)%in%MM.xp.3b$CDR3.amino.acid.sequence)
    

seq.ana <-
  function(path.fas, treeplot=F)
  # Working with reciblast
      {
                                        # path.fas: the path to the unaligned fasta files
                                        # WARNING: path.fas should end with "/"
                                        # treeplot= F or T, if plotting or not NJ trees
          
          library("ape")

          f.lst = list.files(path.fas, pattern="*.fas")
          out.tree = "trees/"
          dir.create(out.tree)
          
          genes = c()
          dists = c()
          dist.sd = c()
          min.dist = c()
          trees = c()
          gaps = c()
          lens = c()
          gap.ratio = c()
          
          for (f in f.lst)
              {
                  message(c("Computing for file: ", f))
                  ali = read.dna(paste(path.fas, f, sep=""), format = "fas")
                  dist = dist.dna(ali,"raw")
                  ali = as.character(ali)

                  genes <- c(genes, f)
                  dists <- c( dists, mean(dist))
                  dist.sd <- c(dist.sd, sd(dist))
                  min.dist <- c(min.dist, min(dist))
                  gaps <- c( gaps, length( ali[ali == "-"] ) / 6 ) # mean gaps per seq
                  lens <- c(lens, length(ali[1,]))
                  gap.ratio = c(gap.ratio, (length( ali[ali == "-"] ) / 6) / (length(ali[1,])) )
        
                  if (!(NaN %in% dist) && treeplot == T)
                      {
                                        #write.tree(nj(dist), file=out.tree, append=TRUE, tree.names = f)
                          plot(nj(dist))
                          f.name = paste (f, ".pdf", sep="")
                          dev.copy2pdf(file = paste(out.tree, f.name, sep=""))
                      }
              }
          
          summary = cbind( genes[order(dists)],
              round(dists[order(dists)],2),
              round(dist.sd[order(dists)],2),
              round(min.dist[order(dists)],2),
              gaps[order(dists)],
              lens[order(dists)],
              round(gap.ratio[order(dists)],2)
                          )
          write.table(summary, file= "summary.csv", row.names=FALSE,
                      col.names = c("gene", "Genetic_distance_(raw)",
                          "G_dist_sd","min_dist",
                          "Average_gaps_per_seq_in_ali",
                          "Ali_length","Gap ratio"),
                      quote = FALSE, sep = "\t")
  }


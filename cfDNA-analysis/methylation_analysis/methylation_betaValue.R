names<-list.files(path = ".", pattern = ".bed")

for ( n in names)
{
  print(n)
  tab<- read.table(n, header = F)
  tab<- tab[,1:7 ]
  colnames(tab)<- c("chr", "strat" , "end", "methylation",
                    "r_chr", "r_start", "r_end")
  tab<- as.data.frame(tab)
  tab$methylation<- as.numeric(tab$methylation)
  tab$name<- paste0(tab$r_chr,tab$r_start,tab$r_end)
  list_regions<- split(tab, tab$name)
  region_methylations<- lapply(list_regions, function(df) {mean(df$methylation)})
  regions_names<- names(list_regions)
  beta<- unlist(region_methylations)
  out<- data.frame(regions_names, beta)
  write.table(out , paste0("percentage_whole/", n), quote = F ,row.names = F)  
}

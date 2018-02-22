aggregate_datasets = function(input_folder = NULL, out_folder = NULL){

  setwd(input_folder)
  files = list.files(input_folder, full.names = F)
  
  for (i in 1:length(files)){
  
    out_file = paste(out_folder, files[i], sep = "/")
   
    if (!file.exists(out_file)){
    	cat(sprintf("Reading file %s ..\n", files[i]))
    	bed = read.table(files[i], sep = "\t", header = F)
      
        #keep only columns of interest (file path is the last one)
        bed = bed[,c(1:6, ifelse(startsWith(as.character(bed[1,7]), "chr"),10,11))]
    
       cat(sprintf("Aggregating..%s\n", files[i])) #aggregate datasets by chr, start and stop
       bed = as.matrix(bed) 
       res = aggregate(bed[,7] ~ bed[,1] + bed[,2] + bed[,3] + bed[,4] + bed[,5] + bed[,6], bed, paste, collapse = ";")

       #export
       write.table(res, out_file, sep = "\t", col.names = F, row.names = F, quote = F)
    }
  }
}

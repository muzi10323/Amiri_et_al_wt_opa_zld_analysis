import.clogmia.rnaseq.pe.v1 = function(file, qualthresh = 10, exclude.duplicates = FALSE)
{
  # this function will import paired-end data that has been mapped and duplicate-marked. It takes as input a path to a .bamfile. It gives the option to set the map-quality score cutoff as well as whether or not to exclude duplicates. 
  #This function has been adapted from import.clogmia.pe.v3 an importer for ATAC seq reads.
  
  require(GenomicAlignments)
  require(BSgenome.Calbipunctata.Cantata.uni4000)
  ca6 = Calbipunctata
  shortname = substring(file,tail(gregexpr('[/]',file)[[1]],1)+1,nchar(file))
  parent.dir = dirname(dirname(file))
  # if(!dir.exists(paste0(parent.dir,'/sample_properties'))){dir.create(paste0(parent.dir,'/sample_properties')); cat('made a new directory (sample_properties) in the parental directory.\n')}
  
  cat(paste('import.clogmia.rnaseq.pe.v1() on',file,'\n\n'))
  cat('Only importing from the six long scaffolds, Clogmia genome version Calbipunctata.Cantata.uni4000. Reads only with above threshold map quality, non secondary mappings, properly paired.\n\n')
  
  good.uns = paste0('scaffold_', c(1:5,7))
  # chromosomes sized scaffolds are scaffold_ 1,2,3,4,5,7 ie not 6
  # the scaffolds were not sorted on size when delivered from Cantatabio
  
  cat('importing all reads and bam file metadata\n')
  
  if(exclude.duplicates){
    cat("Duplicates will be excluded. \n")
    allreads = readGAlignments(file, param = ScanBamParam(flag = scanBamFlag(isPaired = TRUE,isProperPair = TRUE,isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE, isDuplicate = FALSE), what = c('flag','mrnm','mpos','mapq','isize')),use.names = TRUE)
  }else{cat("Duplicates will not be excluded.\n")
    allreads = readGAlignments(file, param = ScanBamParam(flag = scanBamFlag(isPaired = TRUE,isProperPair = TRUE,isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE), what = c('flag','mrnm','mpos','mapq','isize')),use.names = TRUE)
  }
  
  cat('filtering out reads on short/incomplete scaffolds.\n')
  init = length(allreads)
  allreads = allreads[seqnames(allreads)%in%good.uns]
  after = length(allreads)
  cat(paste('filtered', init-after, 'reads on short/incomplete scaffolds.\n\n'))
  
  cat('filtering all reads by map quality score.\n')
  init = length(allreads)
  allreads = allreads[mcols(allreads)$mapq >= qualthresh]
  after = length(allreads)
  cat(paste('filtered', init-after, 'reads with low map quality.\n\n'))
  
  cat('reconstructing paired end reads from filtered reads.\n')
  pairs = makeGAlignmentPairs(allreads,use.names = TRUE,use.mcols = c('flag','mapq'))
  
  cat('coercing the data into a granges object.\n')
  single.interval = granges(pairs)
  
  cat('Assigning <genome> variable as <ca6>.\n')
  genome(single.interval) = 'ca6'
  
  seqlevels(single.interval) = seqlevels(ca6)
  seqinfo(single.interval) = seqinfo(ca6)
  single.interval = trim(single.interval)
  
  cat('Removing seqlevels not in use.\n')
  seqlevels(single.interval) = seqlevelsInUse(single.interval)
  
  names(single.interval) = NULL
  return(single.interval)
}
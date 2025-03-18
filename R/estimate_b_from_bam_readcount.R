#' Transform the output of bam-readcount to MAF
#'
#' Transform the output of bam-readcount to data.table with minor allele frequency (MAF) estimate
#' Only works with bam-readcount output preprocessed to contain first 10 columns as detailed in vignette
#' @param bam_readcount Output of bam-readcount, after reading in first 10 columns
#' @param vcf The VCF genotype file, to get the minor allele and subset to common SNPs. Only takes into account
#' so far minor alleles which are A,C,T or G (TODO: consider changing this)
#'
#' @return A data.table containing per SNP position: the total counts,
#' counts for the A,C,T and G alleles, the minor allele frequency estimate (b_estimate), the ALT allele,
#' and the SNP name in chromosome_position_refAllele (C_POS_REF) format
#' @export
estimate_b_from_bam_readcount = function(bam_readcount,vcf){

  colnames(bam_readcount) = c("CHROM","POS","REF","TOTAL_READS","DEL","A","C","G","T","N")
  bam_readcount$C_POS_REF = paste(bam_readcount$CHROM,bam_readcount$POS,bam_readcount$REF,sep = "_")
  # Remove duplicates - there are a few in the bam-readcount file (why?)
  message("...Removing duplicates from bam-readcount file")
  bam_readcount = unique(bam_readcount, by="C_POS_REF")

  # The VCF genotype file to get the minor allele
  vcf_reduced = data.table::as.data.table(vcf@fix[,1:5])
  # Eliminate SNPs where ALT alleles are not A,C,T or G
  vcf_reduced = vcf_reduced[ALT %in% c("A","C","T","G") ]
  vcf_reduced$C_POS_REF = paste(vcf_reduced$CHROM,vcf_reduced$POS,vcf_reduced$REF,sep = "_")
  message("...Removing duplicates from vcf file")
  vcf_reduced = unique(vcf_reduced, by="C_POS_REF")

  message("Subsetting bam-readcount and vcf to common SNPs")

    if(nrow(bam_readcount)<nrow(vcf_reduced)){
      vcf_reduced = subset(vcf_reduced, C_POS_REF %in% bam_readcount$C_POS_REF )
      bam_readcount = subset(bam_readcount,C_POS_REF %in% vcf_reduced$C_POS_REF)

    }else{
      bam_readcount = subset(bam_readcount,C_POS_REF %in% vcf_reduced$C_POS_REF)
      vcf_reduced = subset(vcf_reduced, C_POS_REF %in% bam_readcount$C_POS_REF )
    }


  message("Calculating minor allele estimate")

  count_dt = data.table(
    total_reads = rep(0, nrow(bam_readcount)),
    A = rep(0, nrow(bam_readcount)),
    C = rep(0, nrow(bam_readcount)),
    `T` = rep(0, nrow(bam_readcount)),
    G = rep(0, nrow(bam_readcount)),
    b_estimate = as.numeric(rep(NA, nrow(bam_readcount)))
  )
  for(column in c("A","C","T","G")){
    sub_match = gsub("^([^:]+:[^:]+).*", "\\1", unlist(bam_readcount[[column]]))
    count_dt[[column]] = as.numeric(gsub(".*:","",sub_match))

  }
  count_dt[["total_reads"]] = bam_readcount[["TOTAL_READS"]]
  count_dt[["ALT"]] = vcf_reduced[["ALT"]]

  for(base in c("A","C", "G", "T")){
    count_dt[count_dt$ALT %in% base, "b_estimate"] = as.double(unlist(count_dt[count_dt$ALT %in% base,..base]) / unlist(count_dt[count_dt$ALT %in% base, "total_reads"]))

  }

  # add row names
  count_dt$C_POS_REF = bam_readcount$C_POS_REF


  if(!identical(nrow(count_dt),nrow(vcf_reduced))) stop("VCF file and minor allele estimate file do not have same number of rows.\n")
  # if(!identical(count_dt$total_reads, rowSums(count_dt[,c("A","C","T","G")]))) stop("The number of total reads is not the sum of the ACTG columns.\n")
  # Sometimes the total reads are not the sum of the ACTG counts. This occurs because there are other alleles present
  # that we are not taking into account (are not ACTG) . But I think it's best to use the total_reads detected
  # by bam-readcount to more accurately calculate the minor allele estimate
  # proportion of SNPs where total reads != sum of ACTG counts
  1-sum(count_dt$total_reads == rowSums(count_dt[,c("A","C","T","G")]))/nrow(count_dt)
  # less than 1%
  if((1-sum(count_dt$total_reads == rowSums(count_dt[,c("A","C","T","G")]))/nrow(count_dt)) > 0.01) stop("The proportion of SNPs where total reads != sum of ACTG counts is over 1%. There might be an error in poodle::estimate_b_from_bam_readcount(), or you may have a big fraction of non-ACTG allele counts. Investigate.\n")


  return(count_dt)

}

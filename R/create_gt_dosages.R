
#' Create genotype dosage file (minor allele dosages)
#'
#' Transform the VCF with genotypes into a genotype dosage data table.
#' @param vcf The VCF file (was tested with VCF version 4.2) that has been read in with vcfR.
#'
#' @return A data table with with column numbers equal to number of donors
#' plus one column with the SNP name in format chromosome_position_refAllele (C_POS_REF)
#' @export

create_gt_dosages = function(vcf){

  # Excluding indels (subsetting only to SNPs) - we may want to keep this info for calculating the proportions?
  # vcf2 <- extract.indels(vcf)

  ## Extracting GT
  gt <- data.table::as.data.table(vcfR::extract.gt(x = vcf))
  # first make new names
  names = paste0(vcf@fix[,"CHROM"],"_",vcf@fix[,"POS"],"_",vcf@fix[,"REF"])

  gt=gt[!duplicated(names),]
  message("Removing ",sum(duplicated(names)), " duplicated rows from vcf.")

  ## removing those from the names
  names=names[!duplicated(names)]
  gt[, ("rn") := names]

  # setkey(gt,rn) # this sorts the data.table, be careful

  ### Some meanings: ################
  # / : genotype unphased (e.g. 0/0)
  #  | : genotype phased (e.g. 0|0)
  # Phased data are ordered along one chromosome and so from these data you know the haplotype
  # i.e. the variants are known to be on the same homologous chromosome because some reads were found to carry both variants
  # Unphased data are simply the genotypes without regard to which one of the pair of chromosomes holds that allele
  ################


  ## converting GT to minor allele score

  gt2 =  gt %>%
    dplyr::mutate_all( list( ~ stringr::str_replace_all(., "0\\|0", "0")))

  gt2 =  gt2 %>%
    dplyr::mutate_all(list( ~ stringr::str_replace_all(., "0\\/0", "0")))

  gt2 = gt2 %>%
    dplyr::mutate_all(list( ~ stringr::str_replace_all(., "0\\|1", "0.5")))

  gt2 = gt2 %>%
    dplyr::mutate_all(list( ~ stringr::str_replace_all(., "1\\|0", "0.5")))

  gt2 =  gt2 %>%
    dplyr::mutate_all(list( ~ stringr::str_replace_all(., "0\\/1", "0.5")))
  gt2 =  gt2 %>%
    dplyr::mutate_all(list( ~ stringr::str_replace_all(., "1\\/0", "0.5")))

  gt2 = gt2 %>%
    dplyr::mutate_all(list( ~ stringr::str_replace_all(., "1\\|1", "1")))
  gt2 = gt2 %>%
    dplyr::mutate_all(list( ~ stringr::str_replace_all(., "1\\/1", "1")))

  # change donor columns to numeric
  cols = setdiff(colnames(gt2),"rn")
  gt2[ ,(cols) := lapply(.SD, as.numeric),.SDcols = cols]

  # hist(gt2)
  ## Some donors have NAs at some SNP positions
  anyNA(gt2)
  gt2[which(rowSums(is.na(gt2)) != 0)[1], ]

  # Removing SNPs with NAs in at least one donor - do I need to do this? Can I ignore the NAs somehow?
  gt2 = gt2[which(rowSums(is.na(gt2)) == 0), ]

  return(gt2)
}

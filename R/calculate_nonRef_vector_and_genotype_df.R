
#' Convert VCF Genotype Data to Minor Allele Dosages
#'
#' This function processes a VCF (Variant Call Format) object to extract genotype (GT) information and convert it into minor allele dosages. It removes duplicated SNP entries, replaces genotype representations with numeric dosage values, and filters out SNPs with missing data.
#'
#' @param vcf A VCF object of class `vcfR` containing genotype data.
#'
#' @return A data.table with SNP identifiers in column 'rn' and genotype dosages (0, 0.5, or 1) for each donor/sample.
#'
#' @details
#' - The function extracts genotype (GT) data from the VCF object.
#' - It removes duplicate SNP entries based on chromosome, position, and reference allele.
#' - Genotype strings are converted into minor allele dosages:
#'   - `0|0` or `0/0` -> `0`
#'   - `0|1`, `1|0`, `0/1`, or `1/0` -> `0.5`
#'   - `1|1` or `1/1` -> `1`
#' - SNPs with missing values are removed.
#'
#' @examples
#' \dontrun{
#' library(vcfR)
#' library(data.table)
#' library(dplyr)
#' library(stringr)
#' vcf = read.vcfR("example.vcf")
#' dosage_data = create_gt_dosages(vcf)
#' }
#'
#' @import vcfR data.table dplyr stringr
#' @export

create_gt_dosages = function(vcf){

  ## Extracting GT
  gt = as.data.table(vcfR::extract.gt(vcf))
  # first make new names from chromosome, position, ref and alt (so there are no duplicates)
  names = paste0(vcf@fix[,"CHROM"],"_",vcf@fix[,"POS"],"_",vcf@fix[,"REF"])
  # There are rows for same SNP but different minor alleles: those will be removed with the previous line (does not include "ALT")
  # there are still some duplicates so remove those
  gt=gt[!duplicated(names),]
  message("Removing ",sum(duplicated(names)), " duplicated rows from vcf")

  ## removing those from the names
  names=names[!duplicated(names)]
  gt[, ("rn") := names]

  # setkey(gt,rn) # this sorts the data.table, be careful
  object.size(gt)

  ### Some meanings: ################
  # / : genotype unphased (e.g. 0/0)
  #  | : genotype phased (e.g. 0|0)
  # Phased data are ordered along one chromosome and so from these data you know the haplotype
  # i.e. the variants are known to be on the same homologous chromosome because some reads were found to carry both variants
  # Unphased data are simply the genotypes without regard to which one of the pair of chromosomes holds that allele
  ################


  ## converting GT to minor allele score

  gt2 =  gt %>%
    dplyr::mutate(across(everything(), ~ . %>%
                           str_replace_all("0\\|0", "0") %>%
                           str_replace_all("0\\/0", "0") %>%
                           str_replace_all("0\\|1", "0.5") %>%
                           str_replace_all("1\\|0", "0.5") %>%
                           str_replace_all("0\\/1", "0.5") %>%
                           str_replace_all("1\\/0", "0.5") %>%
                           str_replace_all("1\\|1", "1") %>%
                           str_replace_all("1\\/1", "1")))


  # change donor columns to numeric
  cols = setdiff(colnames(gt2),"rn")
  gt2[ ,(cols) := lapply(.SD, as.numeric),.SDcols = cols]

  # hist(gt2)
  ## Some donors have NAs at some SNP positions
  anyNA(gt2)
  gt2[which(rowSums(is.na(gt2)) != 0)[1], ]
  gt[which(rowSums(is.na(gt2)) != 0)[1], ]

  # Removing SNPs with NAs in at least one donor - do I need to do this? Can I ignore the NAs somehow?
  gt2 = gt2[which(rowSums(is.na(gt2)) == 0), ]
  return(gt2)
}


#' Transform the output of bam-readcount to MAF
#'
#' Transform the output of bam-readcount to data.table with minor allele frequency (MAF) estimate
#' Only works with bam-readcount output preprocessed to contain first 10 columns as detailed in vignette
#' @param bam_readcount Output of bam-readcount, after reading in first 10 columns, filtering to matching positions
#' of the VCF (and after passing other user-defined filters such as minimum number of TOTAL_READS, if one so wishes).
#' Colnames are:
#' "CHROM": Chromosome.
#' "POS": Genomic position.
#' "REF": Reference allele.
#' "TOTAL_READS": Total number of reads mapping to the variant.
#' "DEL": Counts matching a deletion ALT allele.
#' "A":  Counts matching "A" ALT allele.
#' "C": Counts matching "C" ALT allele.
#' "G": Counts matching "G" ALT allele.
#' "T": Counts matching "T" ALT allele.
#' "N": Counts matching "N" ALT allele (any ambiguous/unknown single variant allele).
#' "CHROM_POS_REF": union of CHROM, POS, REF columns, separated by "_", for filtering.
#' "ALT": ALT allele considered for this position.
#'
#' @param vcf The VCF genotype file, to get the minor allele and subset to common SNPs. Only takes into account
#' so far minor alleles which are A,C,T or G (TODO: consider changing this)
#'
#' @return
#' - A data.table containing per SNP position: the total counts,
#' counts for the A,C,T and G alleles, the minor allele frequency estimate (b_estimate), the ALT allele,
#' and the SNP name in chromosome_position_REF_ALT (CHROM_POS_REF) format.
#' @export
estimate_b_from_bam_readcount = function(bam_readcount,vcf){

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

    count_dt[["ALT"]] = vcf[["ALT"]]

    for(base in c("A","C", "G", "T")){
      count_dt[count_dt$ALT %in% base, "b_estimate"] = as.double(unlist(count_dt[count_dt$ALT %in% base,..base]) / unlist(count_dt[count_dt$ALT %in% base, "total_reads"]))

    }


  # add row names
  count_dt$CHROM_POS_REF = bam_readcount$CHROM_POS_REF

  return(count_dt)

}


#' Calculate Non-Reference Allele Frequencies and Genotype Dosages
#'
#' This function processes bam-readcount data and VCF files to estimate non-reference allele frequencies and minor allele dosages. It filters SNPs, matches variants between BAM and VCF data, and calculates minor allele frequency estimates.
#'
#' @param bam_readcount_path Path to the bam-readcount file (compressed `.gz`).
#' @param vcf_path Path to the VCF file containing variant data.
#' @param min_read_overlap Minimum number of reads required to retain a SNP (default = 1)
#' @param b_estimate_path Optional path to save minor allele frequency estimates.
#' @param variants_retained_path Optional path to save statistics on retained variants.
#' @param genotype_minor_allele_dos_path Optional path to save minor allele dosage data.
#'
#' @return A list containing:
#'   - `b_estimate_table`: A data.table with estimated non-reference allele frequencies.
#'   - `ma_dosages`: A data.table with minor allele dosages.
#'
#' @details
#' - Reads in bam-readcount and VCF files, retaining only relevant columns.
#' - Filters out variants where the ALT allele is not A, C, T, or G.
#' - Ensures SNPs have a minimum number of reads and removes duplicates.
#' - Matches SNPs between bam-readcounts and the VCF file.
#' - Calculates minor allele dosages using `create_gt_dosages()`.
#' - Estimates non-reference allele frequencies using `estimate_b_from_bam_readcount()`.
#' - Performs consistency checks to ensure data integrity.
#' - Optionally saves intermediate results if file paths are provided.
#'
#' @examples
#' \dontrun{
#' bam_path = "path/to/bam_readcount.gz"
#' vcf_path = "path/to/variants.vcf"
#' result = calculate_nonRef_vector_and_genotype_df(bam_path, vcf_path, min_read_overlap = 10)
#' }
#'
#' @import data.table vcfR
#' @export
calculate_nonRef_vector_and_genotype_df = function(bam_readcount_path,vcf_path,min_read_overlap = 1,
                                                   b_estimate_path = NULL, variants_retained_path = NULL,
                                                   genotype_minor_allele_dos_path=NULL){

  # Reading in only 10 columns - consider selecting all if I look into non-ACGT single bp minor alleles
  # check operative system for selecting correct read function
  sys_name = Sys.info()[["sysname"]]

  if (sys_name == "Linux") {
    bam_readcount = data.table::fread(
      cmd = paste("zcat", bam_readcount_path, "| awk -F '\t' '{print $1 , $2 , $3 , $4, $5, $6, $7, $8, $9, $10}'")
    )

  } else if (sys_name == "Darwin") {
    bam_readcount = data.table::fread(
      cmd = paste("gzcat", bam_readcount_path, "| awk -F '\t' '{print $1 , $2 , $3 , $4, $5, $6, $7, $8, $9, $10}'")
    )

  } else if (sys_name == "Windows") {
    R.utils::gunzip(bam_readcount_path, remove = FALSE)
    decompressed_path = sub("\\.gz$", "", bam_readcount_path)

    bam_readcount = data.table::fread(decompressed_path, select = 1:10)

  } else {
    stop(paste("Unsupported OS:", sys_name))
  }

  colnames(bam_readcount) = c("CHROM","POS","REF","TOTAL_READS","DEL","A","C","G","T","N")
  bam_readcount[["REF"]] = toupper(bam_readcount[["REF"]])

  # adding ALT information from vcf
  vcf = vcfR::read.vcfR(vcf_path)
  vcf_reduced = as.data.table(vcf@fix[,1:5])
  vcf_reduced$CHROM_POS_REF = paste(vcf_reduced$CHROM,vcf_reduced$POS,vcf_reduced$REF,sep = "_")

  # Eliminate SNPs where ALT alleles are not A,C,T or G
  vcf_reduced = vcf_reduced[ALT %in% c("A","C","T","G") ]

  # Lower case in REF means that base was seen on the reverse strand.
  # removing those positions with fewer than min_read_overlap reads
  message("Initial bam-readcount table has length:")
  nrow(bam_readcount)
  message("Summary of total number of reads per SNP")
  summary(bam_readcount[["TOTAL_READS"]])
  message("How many SNPs are covered by more than one read?")
  table(bam_readcount[["TOTAL_READS"]]>1)

  message("Subsetting to SNPs covered by more than ",min_read_overlap," reads.")
  bam_readcount = bam_readcount[TOTAL_READS >= min_read_overlap]
  message("There are ",  nrow(bam_readcount), " variants left after filter." )

  bam_readcount$CHROM_POS_REF = paste(bam_readcount$CHROM,bam_readcount$POS,toupper(bam_readcount[["REF"]]),sep = "_")
  message("...Removing duplicates and subsetting data tables to common variants...")
  # remove duplicates
  # Remove duplicates - there are a few in the bam-readcount file
  # now removing by ref only, so other alt alleles will be eliminated
  # but this avoids counting the same read mapping to ref twice, which can potentially
  # cause problems
  bam_readcount = bam_readcount[!duplicated(bam_readcount$CHROM_POS_REF),]
  vcf_reduced = vcf_reduced[!duplicated(vcf_reduced$CHROM_POS_REF),]

  # subset to shared positions
  if(nrow(bam_readcount)<nrow(vcf_reduced)){
    vcf_reduced = vcf_reduced[vcf_reduced$CHROM_POS_REF %in% bam_readcount$CHROM_POS_REF,]
    bam_readcount = bam_readcount[bam_readcount$CHROM_POS_REF %in% vcf_reduced$CHROM_POS_REF,]
  }else{
    bam_readcount = bam_readcount[bam_readcount$CHROM_POS_REF %in% vcf_reduced$CHROM_POS_REF,]
    vcf_reduced = vcf_reduced[vcf_reduced$CHROM_POS_REF %in% bam_readcount$CHROM_POS_REF,]

  }
  original_vcf_rows = nrow(vcf_reduced) # counting here because bam-readcount was subset to overlapping reads

  if(identical(bam_readcount$CHROM_POS_REF,vcf_reduced$CHROM_POS_REF)){
    bam_readcount$ALT = vcf_reduced$ALT
  }else(
    stop("Bam-readcount output and VCF don't produce identical subsets.")
  )

  message("There are ",  nrow(bam_readcount), " variants left after subsetting to shared variants." )

  # checking if bam-readcount has already handled the strandedness of SNPs
  if(sum(toupper(bam_readcount[["REF"]])==bam_readcount$ALT)!=0){
    stop("Some REF alleles from bam-readcounts are VCF ALTs. Investigate...")

  }

  message("...Calculating minor allele dosages...")

  #it only uses ref and not ref+alt
  ma_dosages = create_gt_dosages(vcf)

  # subset vcf to those that remain

  if(nrow(ma_dosages)<nrow(vcf_reduced)){
    vcf_reduced = vcf_reduced[vcf_reduced$CHROM_POS_REF %in% ma_dosages$rn,]
  }else{
    ma_dosages = ma_dosages[ma_dosages$rn %in% vcf_reduced$CHROM_POS_REF,]

  }


  # Are bam_readcount and vcf_reduced of same length and in the same order?
  if(!identical(vcf_reduced$CHROM_POS_REF,bam_readcount$CHROM_POS_REF)) stop("VCF and bam-readcount file are not identical after filtering.\n")
  if(!identical(ma_dosages$rn,bam_readcount$CHROM_POS_REF)) stop("Dosage file and bam-readcount file are not identical after filtering.\n")


  count_dt = estimate_b_from_bam_readcount(bam_readcount = bam_readcount, vcf = vcf_reduced)
  message("The proportion of rows with b_estimate from the input vcf was ",
          sum(table(count_dt$b_estimate)) / original_vcf_rows, " \n")

  if(ncol(count_dt)!=8) stop("Error: The output file does not have the expected number of columns.\n")
  if(sum(is.na(count_dt$b_estimate)) == length(count_dt$b_estimate)) stop("Error: All the minor allele frequency estimates are NA.\n")


  # Sometimes the total reads are not the sum of the ACTG counts. This occurs because there are other alleles present
  # that we are not taking into account (are not ACTG) . But I think it's best to use the total_reads detected
  # by bam-readcount to more accurately calculate the minor allele estimate
  # proportion of SNPs where total reads != sum of ACTG counts
  1-sum(count_dt$total_reads == rowSums(count_dt[,c("A","C","T","G")]))/nrow(count_dt)
  # less than 1%
  if((1-sum(count_dt$total_reads == rowSums(count_dt[,c("A","C","T","G")]))/nrow(count_dt)) > 0.01) stop("The proportion of SNPs where total reads != sum of ACTG counts is over 1%. There might be an error in poodle::estimate_b_from_bam_readcount(), or you may have a big fraction of non-ACTG allele counts. Investigate.\n")


  if(!is.null(b_estimate_path)){

  message("Saving minor allele frequency estimate data...\n")
  fwrite(count_dt, file = b_estimate_path)
}
  ## Check again that minor allele dosage file and b_estimate file have same snps in same order

  ma_dosages = ma_dosages[ma_dosages$rn %in% vcf_reduced$CHROM_POS_REF,]
  if(!identical(ma_dosages$rn,bam_readcount$CHROM_POS_REF)) stop("Error: Dosage file and bam-readcount file are not identical at the end of script.\n")


  if(!is.null(variants_retained_path)){
    message("Saving stats for original number of variants and variants used...\n")

    write.table(data.frame(original_length = original_vcf_rows, final_length = nrow(ma_dosages)),
                file = variants_retained_path,quote = F, row.names = F, col.names = T)
  }
  if(!is.null(genotype_minor_allele_dos_path)){
    message("Saving minor allele dosage data...\n")

  fwrite(ma_dosages,genotype_minor_allele_dos_path)
}
  return(list("b_estimate_table" = count_dt, "ma_dosages" = ma_dosages))

}

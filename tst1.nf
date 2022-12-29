nextflow.enable.dsl=2

params.input = "./chr*.vcf.gz"
params.vcffiles = params.vcfdir + params.input
vcf_ch = Channel.fromPath(params.vcffiles)


// references
params.ref_dir = '/impute-data/reference/'
params.fasta_ref = 'Homo_sapiens_assembly38.fasta'
params.eagle_map = 'genetic_map_hg38_withX.txt'
//params.beagle = 'beagle.22Jul22.46e.jar'


process tabix {
  //module "bcftools/1.16"
  publishDir "./"
  input:
  path vcf

  output:
  path "${vcf}.tbi"

  script:
  "tabix -f ${vcf}"
}

process qc1 {
  //module "bcftools/1.16"
  publishDir "./"
  input:
  path vcf

  output:
  path "${vcf.SimpleName}_qc1.bcf"

  script:
  """
  bcftools norm -m- -Oz -o ${vcf.SimpleName}_qc1.bcf ${vcf}
  """
}

process qc2 {
  //module "bcftools/1.16"
  publishDir "./"

  input:
  path bcf

  output:
  path "${(bcf =~ /chr\d{1,2}/)[0]}_qc2.bcf"

  script:
  """
  bcftools query -f '%CHROM\\t%POS\\t%REF%ALT\\n' ${bcf} > ${bcf.SimpleName}.pos
  printf "chr1\\t1\\n" > ${bcf.SimpleName}.atgc
  awk '{if ( \$3!="AC" && \$3!="CA" && \$3!="CT" && \$3!="TC" && \$3!="TG" && \$3!="GT" && \$3!="AG" && \$3!="GA" ) print \$1"\t"\$2}' ${bcf.SimpleName}.pos >> ${bcf.SimpleName}.atgc
  bcftools view -T ^${bcf.SimpleName}.atgc -Ob -o ${(bcf =~ /chr\d{1,2}/)[0]}_qc2.bcf ${bcf}
  """
}

process qc3 {
  //module "bcftools/1.16"
  publishDir "./"
  input:
  path bcf

  output:
  path "${(bcf =~ /chr\d{1,2}/)[0]}_qc3.bcf"

  script:
  """
  bcftools norm -m+ ${bcf} | bcftools view -v snps -m2 -M2 | bcftools annotate -x FORMAT,INFO | bcftools +fill-tags -- -t all | bcftools view -i 'AC>0' | bcftools view -i 'F_MISSING<0.05' | python3 -c "import sys;[print (line.strip().replace('|','/').replace('1/0','0/1')) if line.startswith('#')==False else print (line.strip()) for line in sys.stdin]" | bcftools view -Ob -o ${(bcf =~ /chr\d{1,2}/)[0]}_qc3.bcf && tabix -f ${(bcf =~ /chr\d{1,2}/)[0]}_qc3.bcf
  """
}

process qc4 {
  //module "any/plink2/20211217"
  //module "bcftools/1.16"
  publishDir "./"
  input:
  path bcf

  output:
  path "${(bcf =~ /chr\d{1,2}/)[0]}_qc4.bcf"

  script:
  """
  plink2 --bcf ${bcf} --missing --out ${bcf}
  cat ${bcf}.smiss | sed 1d | awk '{if(\$4 > 0.05) print \$1}' > ${bcf}_mind.txt
  bcftools view -S ^${bcf}_mind.txt ${bcf} | bcftools annotate -x FORMAT | bcftools +fill-tags -- -t all | bcftools view -Ob -o ${(bcf =~ /chr\d{1,2}/)[0]}_qc4.bcf && tabix -f ${(bcf =~ /chr\d{1,2}/)[0]}_qc4.bcf
  """
}

process qc5 {
  //module "bcftools/1.16"
  publishDir "./", mode: 'copy', overwrite: true
  input:
  path bcf

  output:
  path "${(bcf =~ /chr\d{1,2}/)[0]}_qc5.bcf"
  path "${(bcf =~ /chr\d{1,2}/)[0]}_qc5.bcf.csi"
  path "${(bcf =~ /chr\d{1,2}/)[0]}_qc5.log"

  script:
  """
  bcftools +fixref ${bcf} -Ob -o ${(bcf =~ /chr\d{1,2}/)[0]}_qc5.bcf -- -f ${params.fasta_ref} -m flip > ${(bcf =~ /chr\d{1,2}/)[0]}_qc5.log 2>&1
  tabix -f ${(bcf =~ /chr\d{1,2}/)[0]}_qc5.bcf
  """
}

process phasing {
  //module "htslib/1.16"
  cpus 1
  publishDir "./"

  input:
  path bcf

  output:
  path "${bcf.SimpleName}_PHASED.vcf.gz"
  val "${chr_nr=(bcf =~ /chr\d{1,2}/)[0]}"

  script:
  """
  eagle --vcf ${bcf} --outPrefix ${bcf.SimpleName}_PHASED --geneticMapFile ${params.ref_dir}/${params.eagle_map} --numThreads ${task.cpus} --pbwtIters 1
  tabix -f ${bcf.SimpleName}_PHASED.vcf.gz
  """

}
process impute {
   //module "any/jdk/1.8.0_265"
   cpus 4
   memory '8 GB'
   publishDir "./"
   input:
	 val ch
     path bcf

   output:
     path "${ch}_IMP.vcf.gz"

   script:
    """
    java -Xmx8g -jar "\$BEAGLE" gp=true gt=${bcf} ref=/gpfs/space/GI/GV/Projects/EGV_hg38/combined/filter/2.2_variant-based-QC/phased/${ch}_dbSNP155_PHASED.vcf.gz map=${params.ref_dir}/plink.${ch}.hg38.map out=${ch}_IMP chrom=${ch} nthreads=${task.cpus}
	"""


   // """
	//java -Xmx8g -jar ${params.beagle} gp=true gt=${bcf} ref=/gpfs/space/GI/GV/Projects/EGV_hg38/combined/filter/2.2_variant-based-QC/phased/${ch}_dbSNP155_PHASED.vcf.gz map=/gpfs/space/home/a73038/WORK/Imputation_server/Plink_maps/plink.${ch}.GRCh38.map out=${ch}_IMP chrom=${(ch=~/\d{1,2}/)[0]} nthreads=${task.cpus}
	//"""
}

workflow {
  // view all input files 
  vcf_ch.view()
  tabix (vcf_ch)
  qc1 (vcf_ch)
  qc2 (qc1.out)
  qc3 (qc2.out)
  qc4 (qc3.out)
  qc5 (qc4.out)
  phasing (qc5.out[0])

  impute (phasing.out[1], phasing.out[0])
}

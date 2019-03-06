configfile: 'config.yaml'



#module load python/3.5.1 samtools bedtools bwa r/3.5.0 rstudio/1.1.453 vcftools freebayes
#PATH=$PATH:/nas/longleaf/home/csoeder/modules/vcflib/bin:/nas/longleaf/home/csoeder/modules/parallel/bin


ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}
sample_by_name = {c['name'] : c for c in config['data_sets']}



def return_file_relpath_by_sampname(wildcards):
	sampname = wildcards.samplename
	pathprefix = sample_by_name[sampname]["path"]
	filesin = return_filename_by_sampname(sampname)
	pathsout = ["".join([pathprefix, fq]) for fq in filesin]
	return pathsout


def return_filename_by_sampname(sampname):
	filenames = []
	if sample_by_name[sampname]['paired']:
		filenames.append(sample_by_name[sampname]['readsfile1'])
		filenames.append(sample_by_name[sampname]['readsfile2'])
	else:
		filenames.append(sample_by_name[sampname]['readsfile'])
	return filenames






rule reference_genome_reporter:
	input:
		fai_in = lambda wildcards: ref_genome_by_name[wildcards.ref_gen]['fai'],
	output:
		report_out = "meta/reference_genomes/{ref_gen}.fai.report"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	shell:
		"""
		mkdir -p meta/reference_genomes/
		cat {input.fai_in} | awk '{{sum+=$2}} END {{ print "number_contigs\t",NR; print "number_bases\t",sum}}' | sed -e 's/^/{wildcards.ref_gen}\t/g' > {output.report_out};
		"""

rule demand_reference_genome_summary:
	input:
		refgen_reports = lambda wildcards: expand("meta/reference_genomes/{ref_gen}.fai.report", ref_gen=ref_genome_by_name.keys())
	output:
		refgen_summary = "meta/reference_genomes.summary"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	shell:
		"cat {input.refgen_reports} > {output.refgen_summary}"




rule fastp_clean_sample_se:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
	output:
		fileOut = ["{pathprefix}/{samplename}.clean.R0.fastq"],
		jason = "{pathprefix}/{samplename}.False.json"
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json {pathprefix}/{samplename}.False.json",# --html meta/FASTP/{samplename}.html", 
		se_params = "",
	message:
		"FASTP QA/QC on single-ended reads ({wildcards.samplename}) in progress.... "
	shell:
		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.se_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]}"


rule fastp_clean_sample_pe:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
	output:
		fileOut = ["{pathprefix}/{samplename}.clean.R1.fastq","{pathprefix}/{samplename}.clean.R2.fastq"],
		jason = "{pathprefix}/{samplename}.True.json"
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json {pathprefix}/{samplename}.True.json",# --html meta/FASTP/{samplename}.html", 
		pe_params = "--detect_adapter_for_pe --correction",
	message:
		"FASTP QA/QC on paired-ended reads ({wildcards.samplename}) in progress.... "
	shell:
		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.pe_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]} --in2 {input.fileIn[1]} --out2 {output.fileOut[1]}"


rule FASTP_summarizer:
	input: 
		jason = lambda wildcards: expand("{path}{samp}.{pairt}.json", path=sample_by_name[wildcards.samplename]['path'], samp = wildcards.samplename, pairt = sample_by_name[wildcards.samplename]['paired'])
	output:
		jason_pruned = "meta/FASTP/{samplename}.json.pruned"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	message:
		"Summarizing reads for sample ({wildcards.samplename}) .... "	
	shell:
		"""
		cp {input.jason} meta/FASTP/{wildcards.samplename}.json
		python3 scripts/fastp_reporter.py {input.jason} {output.jason_pruned} -t {wildcards.samplename}
		"""



rule demand_FASTQ_analytics:	#forces a FASTP clean
	input:
		jasons_in = expand("meta/FASTP/{samplename}.json.pruned", samplename = sample_by_name.keys())
	output:
		summary = "meta/sequenced_reads.dat"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"Collecting read summaries for all samples ...."
	shell:
		"cat {input.jasons_in} > {output.summary}"





rule bwa_align:
	input:
		reads_in = lambda wildcards: expand("{path}{sample}.clean.R{arr}.fastq", path=sample_by_name[wildcards.sample]['path'], sample=wildcards.sample, arr=[ [1,2] if sample_by_name[wildcards.sample]['paired'] else [0] ][0]),
		ref_genome_file = lambda wildcards: ref_genome_by_name[wildcards.ref_genome]['path'],
	output:
		bam_out = "mapped_reads/{sample}.vs_{ref_genome}.bwa.sort.bam",
	params:
		runmem_gb=64,
		runtime="64:00:00",
		cores=8,
	message:
		"aligning reads from {wildcards.sample} to reference_genome {wildcards.ref_genome} .... "
	run:
		shell("bwa aln {input.ref_genome_file} {input.reads_in[0]} > {input.reads_in[0]}.sai ")
		if sample_by_name[wildcards.sample]['paired']:
			shell("bwa aln {input.ref_genome_file} {input.reads_in[1]} > {input.reads_in[1]}.sai ")
			shell("bwa sampe {input.ref_genome_file} {input.reads_in[0]}.sai {input.reads_in[1]}.sai {input.reads_in[0]}  {input.reads_in[1]} | samtools view -Shb | samtools addreplacerg -r ID:{wildcards.sample} -r SM:{wildcards.sample} - | samtools sort -o {output.bam_out} - ")
		else:
			shell("bwa samse {input.ref_genome_file} {input.reads_in[0]}.sai {input.reads_in[0]} | samtools view -Shb | samtools addreplacerg -r ID:{wildcards.sample} -r SM:{wildcards.sample} - | samtools sort -o {output.bam_out} - ")
		shell("samtools index {output.bam_out}")



rule bwa_uniq:
	input:
		bam_in = "mapped_reads/{sample}.vs_{ref_genome}.bwa.sort.bam"
	output:
		bam_out = "mapped_reads/{sample}.vs_{ref_genome}.bwaUniq.sort.bam"
	params:
		quality="-q 20 -F 0x0100 -F 0x0200 -F 0x0300 -F 0x04",
		uniqueness="XT:A:U.*X0:i:1.*X1:i:0",
		runmem_gb=16,
		runtime="6:00:00",
		cores=4,
	message:
		"filtering alignment of {wildcards.sample} to {wildcards.ref_genome} for quality and mapping uniqueness.... "	
	run:
		ref_genome_file=ref_genome_by_name[wildcards.ref_genome]['path']
		shell("samtools view {params.quality} {input.bam_in} | grep -E {params.uniqueness} | samtools view -bS -T {ref_genome_file} - | samtools addreplacerg -r ID:{wildcards.sample} -r SM:{wildcards.sample} - | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - {output.bam_out}")
		shell("samtools index {output.bam_out}")





rule bam_reporter:
	input:
		bam_in = "mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam"
	output:
		report_out = "meta/BAMs/{sample}.vs_{ref_genome}.{aligner}.summary"
	params:
		runmem_gb=8,
		runtime="4:00:00",
		cores=1,
	message:
		"Collecting metadata for the {wildcards.aligner} alignment of {wildcards.sample} to {wildcards.ref_genome}.... "
	run:
		ref_genome_idx=ref_genome_by_name[wildcards.ref_genome]['fai']
		shell("samtools idxstats {input.bam_in} > {input.bam_in}.idxstats")
		shell("samtools flagstat {input.bam_in} > {input.bam_in}.flagstat")
		shell("bedtools genomecov -max 1 -ibam {input.bam_in} -g {ref_genome_idx} > {input.bam_in}.genomcov")
		#change the -max flag as needed to set 
		shell("""samtools depth -a {input.bam_in} | awk '{{sum+=$3; sumsq+=$3*$3}} END {{ print "average_depth\t",sum/NR; print "std_depth\t",sqrt(sumsq/NR - (sum/NR)**2)}}' > {input.bam_in}.dpthStats""")
		#https://www.biostars.org/p/5165/
		#save the depth file and offload the statistics to the bam_summarizer script?? 
		shell("python3 scripts/bam_summarizer.py -f {input.bam_in}.flagstat -i {input.bam_in}.idxstats -g {input.bam_in}.genomcov -d {input.bam_in}.dpthStats -o {output.report_out} -t {wildcards.sample}")


rule demand_BAM_analytics:
	input:
		bam_reports = lambda wildcards: expand("meta/BAMs/{sample}.vs_{ref_genome}.{aligner}.summary", sample=sample_by_name.keys(), ref_genome=wildcards.ref_genome, aligner=wildcards.aligner)
	output:
		full_report = "meta/alignments.vs_{ref_genome}.{aligner}.summary"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	shell:
		"cat {input.bam_reports} > {output.full_report}"




rule window_maker:
	output:
		windowed='utils/{ref_genome}_w{window_size}_s{slide_rate}.windows.bed'
	params:
		runmem_gb=8,
		runtime="5:00",
		cores=1,
	run:
		fai_path = ref_genome_by_name[wildcards.ref_genome]['fai'],
		shell("mkdir -p utils")
		shell(
			'bedtools makewindows -w {wildcards.window_size} -s {wildcards.slide_rate} -g {fai_path} -i winnum | bedtools sort -i - > {output.windowed}'
		)


rule joint_vcf_caller_parallel:
	input:
		bams_in = lambda wildcards: expand("mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam", sample=sample_by_name.keys(), ref_genome=wildcards.ref_genome, aligner=wildcards.aligner),
		windows_in = "utils/{ref_genome}_w100000_s100000.windows.bed"
	output:
		vcf_out = "variants/all_samples.vs_{ref_genome}.{aligner}.vcf"
	params:
		freebayes="--standard-filters",
		runmem_gb=32,
		runtime="64:00:00",
		cores=24,
	message:
		"Jointly calling variants from all samples mapped to \ {wildcards.ref_genome} \ with \ {wildcards.aligner} \ "
	run:
		ref_genome_file=ref_genome_by_name[wildcards.ref_genome]['path']
		shell("""cat {input.windows_in}| awk '{{print$1":"$2"-"$3}}' > {input.windows_in}.rfmt""")
		shell("scripts/freebayes-parallel {input.windows_in}.rfmt {params.cores} {params.freebayes} -f {ref_genome_file} {input.bams_in} | vcftools --vcf - --recode --recode-INFO-all --stdout  > {output.vcf_out} ")



rule vcf_reporter:
	input:
		vcf_in = "variants/{prefix}.vs_{ref_genome}.{aligner}.vcf"
	output:
		report_out = "meta/VCFs/{prefix}.vs_{ref_genome}.{aligner}.summary",
# 		frq_out = "meta/VCFs/{prefix}.vs_{ref_genome}.{aligner}.summary.frq"
	params:
		runmem_gb=8,
		runtime="4:00:00",
		cores=4,
	message:
		"Collecting metadata for the {wildcards.aligner} alignment of {wildcards.prefix} to {wildcards.ref_genome}.... "
	shell:
		"""
		mkdir -p meta/VCFs/

		vcftools --remove-indels --vcf {input.vcf_in} --recode --recode-INFO-all --stdout | grep -v "#" | cut -f 1 | sort | uniq -c > {output.report_out}.snpsPerContig.tmp
		cat {output.report_out}.snpsPerContig.tmp | awk '{{sum+=$1}} END {{ print sum,"\ttotal"}}' | cat {output.report_out}.snpsPerContig.tmp - > {output.report_out}.snpsPerContig
		rm {output.report_out}.snpsPerContig.tmp

		vcftools --keep-only-indels --vcf {input.vcf_in} --recode --recode-INFO-all --stdout | grep -v "#" | cut -f 1 | sort | uniq -c > {output.report_out}.indelsPerContig.tmp
		cat {output.report_out}.indelsPerContig.tmp | awk '{{sum+=$1}} END {{ print sum,"\ttotal"}}' | cat {output.report_out}.indelsPerContig.tmp - > {output.report_out}.indelsPerContig	
		rm {output.report_out}.indelsPerContig.tmp

		vcftools  --vcf {input.vcf_in} --out {output.report_out} --freq 
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --counts
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --missing-indv
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --missing-site
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --singletons

		allen={wildcards.aligner}

		tail -n 1 {output.report_out}.snpsPerContig | awk '{{print "total_snp_count\t"$1}}' | sed -e 's/^/'$allen'\t/g' > {output.report_out}
		tail -n 1 {output.report_out}.indelsPerContig | awk '{{print "total_indel_count\t"$1}}' | sed -e 's/^/'$allen'\t/g'  >> {output.report_out}

		"""

# 	#cat  all_samples.vs_droSec1.bwaUniq.summary.frq.count| cut -f 3 | tail -n +2 | sort | uniq -c
# 	#####	bi, tri, and quadralelic counts ^^ 
# 	#replace some of this with vcftools::vcf-stats ?

rule summon_VCF_analytics_base:
	input:
		vcf_reports = lambda wildcards: expand("meta/VCFs/{prefix}.vs_{ref_genome}.{aligner}.summary", prefix=wildcards.prefix, ref_genome="dm6", aligner=["bwa","bwaUniq"])
	output:
		full_report = "meta/{prefix}.vs_{ref_genome}.calledVariants.summary"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	shell:
		"""
		prefix={wildcards.prefix}
		cat {input.vcf_reports} | sed -e 's/^/'$prefix'\t/g'> {output.full_report}
		"""





rule write_report:
	input:
		reference_genome_summary = ["meta/reference_genomes.summary"],
		sequenced_reads_summary=["meta/sequenced_reads.dat"],
		alignment_summaries = expand("meta/alignments.vs_dm6.{aligner}.summary", aligner=['bwa','bwaUniq']),
	output:
		pdf_out="thingy.pdf"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=2,
	message:
		"writing up the results.... "
	run:
		pandoc_path="/nas/longleaf/apps/rstudio/1.0.136/bin/pandoc"
		pwd = subprocess.check_output("pwd",shell=True).decode().rstrip()+"/"
		shell(""" R -e "setwd('{pwd}');Sys.setenv(RSTUDIO_PANDOC='{pandoc_path}')" -e  "peaDubDee='{pwd}'; rmarkdown::render('scripts/EJgrep_summary.Rmd',output_file='{pwd}{output.pdf_out}')"  """)




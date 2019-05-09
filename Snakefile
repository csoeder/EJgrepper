configfile: 'config.yaml'



#module load python/3.5.1 samtools bedtools bwa r/3.5.0 rstudio/1.1.453 vcftools freebayes bamtools
#PATH=$PATH:/nas/longleaf/home/csoeder/modules/vcflib/bin:/nas/longleaf/home/csoeder/modules/parallel/bin:/nas/longleaf/home/csoeder/modules/pindel


ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}
sample_by_name = {c['name'] : c for c in config['data_sets']}
sampname_by_group = {}
for s in sample_by_name.keys():
	subgroup_lst = sample_by_name[s]['subgroups']
	for g in subgroup_lst:
		if g in sampname_by_group.keys():
			sampname_by_group[g].append(s)
		else:
			sampname_by_group[g] = [s]

parents_by_group = {g : []  for c in config['data_sets'] for g in c['subgroups'] }
for s in sample_by_name.keys():
	ped = sample_by_name[s]['pedigree']
	if ped == 'parent':
		for grup in sample_by_name[s]['subgroups']:
			parents_by_group[grup].append(s)
for g in parents_by_group.keys():
	parents_by_group[g] = list(set(parents_by_group[g]))



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
		jasons_in = expand("meta/FASTP/{samplename}.json.pruned", samplename = sampname_by_group['all'])
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
		quality="-q 20 -F 0x0100 -F 0x0200 -F 0x0300 ",
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
		bam_reports = lambda wildcards: expand("meta/BAMs/{sample}.vs_{ref_genome}.{aligner}.summary", sample=sampname_by_group['all'], ref_genome=wildcards.ref_genome, aligner=wildcards.aligner)
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


rule call_VCF_by_group_parallel:
	input:
		bams_in = lambda wildcards: expand("mapped_reads/{group}.vs_{ref_genome}.{aligner}.sort.bam", group=sampname_by_group[wildcards.group], ref_genome=wildcards.ref_genome, aligner=wildcards.aligner),
		windows_in = "utils/{ref_genome}_w100000_s100000.windows.bed"
	output:
		cnv_map = "utils/{group}.vs_{ref_genome}.{aligner}.cnv",
		vcf_out = "variants/{group}.vs_{ref_genome}.{aligner}.vcf"
	params:
		freebayes="--standard-filters",
		runmem_gb=32,
		runtime="12:00:00",
		cores=24,
	message:
		"Jointly calling variants from {wildcards.group} mapped to \ {wildcards.ref_genome} \ with \ {wildcards.aligner} \ "
	run:
		ref_genome_file=ref_genome_by_name[wildcards.ref_genome]['path']
		ref_fai = ref_genome_by_name[wildcards.ref_genome]['fai']

		shell("""rm -f {output.cnv_map} """)
		for sampname in sampname_by_group[wildcards.group]:
			shell(""" cat {ref_fai} | grep -v "chr[XY]" | awk '{{print $1,"1",$2,"{sampname}",2}}' | tr " " "\t" >> {output.cnv_map} """)
			if sample_by_name[sampname]['sex'] == "M":
				shell(""" cat {ref_fai} | grep "chrX" | awk '{{print $1,"1",$2,"{sampname}",1}}' | tr " " "\t" >> {output.cnv_map};
						cat {ref_fai} | grep "chrY" | awk '{{print $1,"1",$2,"{sampname}",1}}' | tr " " "\t" >> {output.cnv_map} """)
			elif sample_by_name[sampname]['sex'] == "F":
				shell(""" cat {ref_fai} | grep "chrX" | awk '{{print $1,"1",$2,"{sampname}",2}}' | tr " " "\t" >> {output.cnv_map};
						cat {ref_fai} | grep "chrY" | awk '{{print $1,"1",$2,"{sampname}",0}}' | tr " " "\t" >> {output.cnv_map} """)
			else:
				shell(""" cat {ref_fai} | grep "chrX" | awk '{{print $1,"1",$2,"{sampname}",0}}' | tr " " "\t" >> {output.cnv_map};
						cat {ref_fai} | grep "chrY" | awk '{{print $1,"1",$2,"{sampname}",0}}' | tr " " "\t" >> {output.cnv_map}""")

		shell("""cat {input.windows_in}| awk '{{print$1":"$2"-"$3}}' > {input.windows_in}.rfmt""")
		shell("scripts/freebayes-parallel {input.windows_in}.rfmt {params.cores} {params.freebayes} --cnv-map {output.cnv_map} -f {ref_genome_file} {input.bams_in} | vcftools --vcf - --recode --recode-INFO-all --stdout  > {output.vcf_out}.tmp ")
		# renamer: (MOVE THIS TO THE VCF CALLER!!!)
		shell("""cat <(cat {output.vcf_out}.tmp | grep "#") <(cat {output.vcf_out}.tmp | grep -v "#" | cut -f 4- | nl -n ln | awk '{{print "id"$0}}'| paste <(cat {output.vcf_out}.tmp| grep -v "#" | cut -f 1,2) - ) > {output.vcf_out} ;""")
		shell("""rm {output.vcf_out}.tmp""")

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
		vcf_reports = lambda wildcards: expand("meta/VCFs/{prefix}.vs_dm6.{aligner}.summary", prefix=["control","deficiency","background"], aligner=["bwaUniq"])
	output:
		full_report = "meta/allGroups.vs_dm6.bwaUniq.calledVariants.summary"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	run:
		shell("rm -f {output.full_report}")
		for report in input.vcf_reports:
			exp = report.split("/")[-1].split('.')[0]
			shell("""
			cat {report} | sed -e 's/^/'{exp}'\t/g' >> {output.full_report}
				""")

# rule contrast_VCFs:
# 	input:
# 		vcf1 = "variants/{prefix1}.vs_{ref_genome}.{aligner1}.vcf",
# 		vcf2 = "variants/{prefix2}.vs_{ref_genome}.{aligner2}.vcf",
# 	output:
# 		site_diff= "meta/VCFs/{prefix1}__{aligner1}.contrast.{prefix2}__{aligner2}.{ref_genome}.diff.sites_in_files",
# 		site_disc= "meta/VCFs/{prefix1}__{aligner1}.contrast.{prefix2}__{aligner2}.{ref_genome}.diff.sites",
# 		site_disc1 = "meta/VCFs/{prefix1}__{aligner1}.contrast.{prefix2}__{aligner2}.{ref_genome}.diff.sites.1",
# 		site_disc2 = "meta/VCFs/{prefix1}__{aligner1}.contrast.{prefix2}__{aligner2}.{ref_genome}.diff.sites.2",
# 		site_discB = "meta/VCFs/{prefix1}__{aligner1}.contrast.{prefix2}__{aligner2}.{ref_genome}.diff.sites.B",
# 		indiv_disc = "meta/VCFs/{prefix1}__{aligner1}.contrast.{prefix2}__{aligner2}.{ref_genome}.disc",
# 		indiv_disc_count = "meta/VCFs/{prefix1}__{aligner1}.contrast.{prefix2}__{aligner2}.{ref_genome}.disc.count",
# 	params:
# 		runmem_gb=8,
# 		runtime="4:00:00",
# 		cores=4,
# 		good_chroms = "--chr chr2L --chr chr2R --chr chr3L --chr chr3R --chr chr4 --chr chrM --chr chrX --chr chrY",
# 	message:
# 		"potatoes for breakfast.... "
# 	shell:
# 		"""
# 		vcftools --vcf {input.vcf1} --diff {input.vcf2} --diff-site --out meta/VCFs/{wildcards.prefix1}__{wildcards.aligner1}.contrast.{wildcards.prefix2}__{wildcards.aligner2}.{wildcards.ref_genome} {params.good_chroms}

# 		vcftools --vcf {input.vcf1} --diff {input.vcf2} --diff-site-discordance --out meta/VCFs/{wildcards.prefix1}__{wildcards.aligner1}.contrast.{wildcards.prefix2}__{wildcards.aligner2}.{wildcards.ref_genome} {params.good_chroms}

# 		tail -n +2 {output.site_disc}| awk '{{if($3 == 1)print;}}' | cut -f 1 | uniq -c > {output.site_disc1}
# 		tail -n +2 {output.site_disc}| awk '{{if($3 == 1)print;}}' | cut -f 1 | uniq -c | awk '{{sum+=$1}} END {{ print sum,"\ttotal"}}' >>{output.site_disc1}

# 		tail -n +2 {output.site_disc}| awk '{{if($3 == 2)print;}}' | cut -f 1 | uniq -c > {output.site_disc2}
# 		tail -n +2 {output.site_disc}| awk '{{if($3 == 2)print;}}' | cut -f 1 | uniq -c | awk '{{sum+=$1}} END {{ print sum,"\ttotal"}}' >>{output.site_disc2}

# 		tail -n +2 {output.site_disc}| awk '{{if($3 == "B")print;}}' | cut -f 1 | uniq -c > {output.site_discB}
# 		tail -n +2 {output.site_disc}| awk '{{if($3 == "B")print;}}' | cut -f 1 | uniq -c | awk '{{sum+=$1}} END {{ print sum,"\ttotal"}}' >>{output.site_discB}

# 		tail -n +2 {output.site_disc} | awk '{{if($3 == "B")print;}}' | awk '{{if($6>0)print;}}' > {output.indiv_disc}

# 		cat {output.indiv_disc}| cut -f 1 | uniq -c > {output.indiv_disc_count}
# 		cat {output.indiv_disc}| cut -f 1 | uniq -c | awk '{{sum+=$1}} END {{ print sum,"\ttotal"}}' >> {output.indiv_disc_count}

# 		"""	


rule VCF_indelWinnower:
	input:
		vcf_in = "variants/{prefix}.vs_dm6.bwaUniq.vcf",
	output:
		alleleCounts_out= "analysis/{prefix}.alleleCounts.simpleIndels.dpthFilt.biallelic.universal.out",
		cleanVcf_out= "variants/{prefix}.vs_dm6.bwaUniq.alleleCounts.simpleIndels.dpthFilt.biallelic.universal.vcf",
	params:
		runmem_gb=8,
		runtime="4:00:00",
		cores=4,
		good_chroms = "--chr chr2L --chr chr2R --chr chr3L --chr chr3R --chr chr4 --chr chrX",
		dpth_filt = 10,
		max_uncalled = 0,
	message:
		"filtering the {wildcards.prefix} VCF for depth of coverage, site call percent, allele count, variant type .... "
	run:
		shell("""vcftools {params.good_chroms} --vcf {input.vcf_in} --max-missing-count {params.max_uncalled} --minDP {params.dpth_filt} --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | grep -v "TYPE=complex" | vcfallelicprimitives | vcftools --vcf - --keep-only-indels --recode --recode-INFO-all --stdout > {output.cleanVcf_out}.tmp""")
		shell("""cat <(grep "#" {output.cleanVcf_out}.tmp ) <(grep -v "#" {output.cleanVcf_out}.tmp | cut -f 4- | nl -n ln | awk '{{print "sdbu"$0}}'| tr -d " " | paste <(cat {output.cleanVcf_out}.tmp| grep -v "#" | cut -f 1,2) - ) > {output.cleanVcf_out} """)
		shell("""rm {output.cleanVcf_out}.tmp""" )
		shell("""vcftools --vcf {output.cleanVcf_out} --counts --stdout | tr ":" "\t" | tail -n +2 | nl -n ln  > {output.alleleCounts_out}""")


rule VCF_novelist:
	input:
		vcf_in = "variants/{group}.vs_{meta}.vcf",
	output:
		novel_vcf = "variants/{group}.vs_{meta}.novel.vcf",
		novel_count = "analysis/{group}.vs_{meta}.novel.count",
		back_vcf = "variants/{group}.vs_{meta}.back.vcf",
		back_count = "analysis/{group}.vs_{meta}.back.count",
#		ancestral_vcf= "variants/{group}.vs_{meta}.ancestral.vcf",
#		ancestral_count= "analysis/{group}.vs_{meta}.ancestral.count",
	params:
		runmem_gb=8,
		runtime="4:00:00",
		cores=4,
		#good_chroms = "--chr chr2L --chr chr2R --chr chr3L --chr chr3R --chr chr4",
		#dpth_filt = 10,
		#max_uncalled = 0,
	message:
		"filtering the {wildcards.group} VCF for variants which do not appear in the parents.... "
	run:
		par_string = " ".join(["--indv %s"%tuple([x]) for x in parents_by_group[wildcards.group]])
		#nopar_string = " ".join(["--remove-indv %s"%tuple([x]) for x in parents_by_group[wildcards.group]])

		#collect all the variant IDs which are 0|0 in both parents
		# subset the VCF to those sites which are 0|0 in both parents, then to those which also have at least one alt allele in an offspring
		shell("""
			vcftools --vcf {input.vcf_in} {par_string} --non-ref-ac 0 --max-non-ref-ac 0  --recode --recode-INFO-all --stdout | grep -v "#" | cut -f 3 > {input.vcf_in}.sdbu0.list;
			vcftools --vcf {input.vcf_in} --snps {input.vcf_in}.sdbu0.list --non-ref-ac 1 --max-non-ref-af 1.0  --recode --recode-INFO-all --stdout > {output.novel_vcf};
			vcftools --vcf {output.novel_vcf} --counts --stdout | tr ":" "\t" | tail -n +2 | nl -n ln  > {output.novel_count};
		""")

		# collect variant IDs which are 1|1 in both parents
		# subset the VCF to those sites which are 1|1 in both parents, then to those which also have at least one ref allele in an offspring
		#tally the allele counts
		shell("""
			vcftools --vcf {input.vcf_in} {par_string} --non-ref-af 1 --recode --recode-INFO-all --stdout | grep -v "#" | cut -f 3 > {input.vcf_in}.sdbu1.list;
			vcftools --vcf {input.vcf_in} --snps {input.vcf_in}.sdbu1.list --max-non-ref-af 0.99999  --recode --recode-INFO-all --stdout > {output.back_vcf};
			vcftools --vcf {output.back_vcf} --counts --stdout | tr ":" "\t" | tail -n +2 | nl -n ln  > {output.back_count};
		""")

		# clean up
		shell("""
			rm {input.vcf_in}.sdbu0.list;
			rm {input.vcf_in}.sdbu1.list;
		""")


#checking the false call rate on chrX vs autosomes in male offspring........
#samtools mpileup -r chr2L:10587-10589 -f /proj/cdjones_lab/Genomics_Data_Commons/genomes/drosophila_melanogaster/dm6.main.fa mapped_reads/w4_3.vs_dm6.bwaUniq.sort.bam mapped_reads/w4_16.vs_dm6.bwaUniq.sort.bam 





#pindel????
#quick and easy insert size w/o picard 
# head -10000 mappings.sam | awk '{if ($9 > 0) {S+=$9; T+=1}}END{print "Mean: " S/T}' 
#https://www.biostars.org/p/16556/
#bamtools stats -in mapped_reads/CantonS.vs_dm6.bwaUniq.sort.bam -insert

#mapped_reads/CantonS.vs_dm6.bwaUniq.sort.bam	232	CantonS


# rule pindelIndv:
# 	input:
# 		bam_in = "mapped_reads/{samp}.vs_dm6.bwaUniq.sort.bam",
# 	output:
# 		pin_out = "pindel/{samp}.vs_dm6.bwaUniq.pindel_D.vcf",
# 	params:
# 		runmem_gb=16,
# 		runtime="4:00:00",
# 		cores=4,
# 	message:
# 		"looking for large variants using PINDEL...."
# 	run:
# 		shell(""" mkdir -p pindel """)
# 		ref_gen = ref_genome_by_name["dm6"]['path']
# 		shell(""" echo {input.bam_in}"\t"$(bamtools stats -in {input.bam_in} -insert | grep Median | cut -f 2 -d : | tr -d " ")"\t"pindel/{wildcards.samp}.vs_dm6.bwaUniq.pindel > pindel/{wildcards.samp}.cfg """)
# 		shell(""" pindel -T {params.cores} -f {ref_gen} -i pindel/{wildcards.samp}.cfg -o pindel/{wildcards.samp}.vs_dm6.bwaUniq.pindel """)
# 		shell(""" pindel2vcf -p pindel/{wildcards.samp}.vs_dm6.bwaUniq.pindel_D -r {ref_gen} -R dm6 -d 19790717 -v {output.pin_out} """)



rule pindelByGroup:
	input:
		bams_in = lambda wildcards: expand("mapped_reads/{group}.vs_{ref_genome}.{aligner}.sort.bam", group=sampname_by_group[wildcards.group], ref_genome="dm6", aligner="bwaUniq"),
	output:
		pin_out = "pindel/{group}.vs_dm6.bwaUniq.pindel_D.vcf",
	params:
		runmem_gb=16,
		runtime="4:00:00",
		cores=4,
	message:
		"looking for large variants in {wildcards.group} using PINDEL...."
	run:
		shell(""" mkdir -p pindel """)
		ref_gen = ref_genome_by_name["dm6"]['path']
		groupies = sampname_by_group[wildcards.group]
		for i in range(0,len(groupies)):
			shell(""" echo {input.bam_in[i]}"\t"$(bamtools stats -in {input.bam_in[i]} -insert | grep Median | cut -f 2 -d : | tr -d " ")"\t"pindel/{groupies[i]}.vs_dm6.bwaUniq.pindel > pindel/{wildcards.samp}.cfg """)
		shell(""" pindel -T {params.cores} -f {ref_gen} -i pindel/{wildcards.samp}.cfg -o pindel/{wildcards.group}.vs_dm6.bwaUniq.pindel """)
		shell(""" pindel2vcf -p pindel/{wildcards.group}.vs_dm6.bwaUniq.pindel_D -r {ref_gen} -R dm6 -d 19790717 -v {output.pin_out}.tmp """)
		shell(""" cat {output.pin_out}.tmp | sed -e 's/pindel\///g' | sed -e 's/.vs_dm6.bwaUniq.pindel//g' > {output.pin_out} """)
		shell(""" rm {output.pin_out}.tmp """)



rule VCF_snpWinnower:
	input:
		vcf_in = "variants/{prefix}.vs_dm6.bwaUniq.vcf",
	output:
		alleleCounts_out= "analysis/{prefix}.alleleCounts.goodSnps.dpthFilt.biallelic.universal.out",
		cleanVcf_out= "variants/{prefix}.vs_dm6.bwaUniq.alleleCounts.goodSnps.dpthFilt.biallelic.universal.vcf",
	params:
		runmem_gb=8,
		runtime="4:00:00",
		cores=4,
		good_chroms = "--chr chr2L --chr chr2R --chr chr3L --chr chr3R --chr chr4 --chr chrX",
		dpth_filt = 10,
		max_uncalled = 0,
	message:
		"filtering the {wildcards.prefix} VCF for depth of coverage, site call percent, allele count, variant type .... "
	run:
		shell("""vcftools {params.good_chroms} --vcf {input.vcf_in} --max-missing-count {params.max_uncalled} --minDP {params.dpth_filt} --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | grep -v "TYPE=complex" | vcfallelicprimitives | vcftools --vcf -  --remove-indels --recode --recode-INFO-all --stdout > {output.cleanVcf_out}.tmp""")
		shell("""cat <(grep "#" {output.cleanVcf_out}.tmp ) <(grep -v "#" {output.cleanVcf_out}.tmp | cut -f 4- | nl -n ln | awk '{{print "ndbu"$0}}'| tr -d " " | paste <(cat {output.cleanVcf_out}.tmp| grep -v "#" | cut -f 1,2) - ) > {output.cleanVcf_out} """)
		shell("""rm {output.cleanVcf_out}.tmp""" )
		shell("""vcftools --vcf {output.cleanVcf_out} --counts --stdout | tr ":" "\t" | tail -n +2 | nl -n ln  > {output.alleleCounts_out}""")



#SNP distance spectrum 
#bedtools closest -io -a <( cut -f 1-9 test.vcf) -b <(cut -f 1-9 test.vcf) -d | cut -f 1-5,10-14,19  > test.dist
#see also bedtools cluster?
#grep -w "ndbu39\|ndbu40" test.vcf | sed -e 's/:\S*//g' | cut -f 1-5,10-
#join -o1.3,2.3,1.8,1.10,2.4,2.9,2.12  -j1 <(cat test.count.out | awk '{print$2":"$3"\t"$0}' | sort -k1,1 ) <( cat test.dist | awk '{print$1":"$2"\t"$0}' | sort -k1,1) | tr " " "\t" > dist_and_ac.test
#cat dist_and_ac.test | awk '{if($3==1||$4==1)print;}' | awk '{if($7<10)print;}' > rare_neighbors.list


rule windowedCover:
	input:
		windoze = "utils/{ref_gen}_w50_s50.windows.bed",
		reads = "mapped_reads/{sample}.vs_{ref_gen}.bwaUniq.sort.bam",
	output:
		cov_file = "mapped_reads/coverage/{sample}.vs_{ref_gen}.bwaUniq.w50_s50.bed",
	params:
		runmem_gb=72,
		runtime="1:00:00",
		cores=2,
	message:
		"measuring the coverage of {wildcards.sample} mapped to {wildcards.ref_gen} in 50bp chunks .... "
	run:
		shell("mkdir -p mapped_reads/coverage/")
		shell("bedtools coverage -a {input.windoze} -b {input.reads} > {output.cov_file}.tmp")
		shell("python3 scripts/coverage_anomalizer.py -i {output.cov_file}.tmp -o {output.cov_file} ")
		shell("rm {output.cov_file}.tmp")


rule pullAllCoverage:
	input:
		lambda wildcards: expand("mapped_reads/coverage/{sample}.vs_dm6.bwaUniq.w50_s50.bed", sample=sampname_by_group['all'])
	output:
		"potato.txt"
	params:
		runmem_gb=1,
		runtime="0:30",
		cores=1,
	message:
		"bluhhhhhhhhhhhhhhhhh"
	run:
		shell("touch {output}")

#cat covtest2.anom.bed | awk '{if($8 < -2 )print;}' | cut -f 1-3 | bedtools merge -i - > barespots.bed



rule offspring_merger:
	input:
		alleleCounts_out= expand("analysis/{prefix}.alleleCounts.simpleIndels.dpthFilt.biallelic.universal.out",prefix=["deficiency","background"]),
		novel_count = expand("analysis/{group}.vs_dm6.bwaUniq.alleleCounts.simpleIndels.dpthFilt.biallelic.universal.novel.count",group=["deficiency","background"]),
		back_count = expand("analysis/{group}.vs_dm6.bwaUniq.alleleCounts.simpleIndels.dpthFilt.biallelic.universal.back.count",group=["deficiency","background"]),
	output:
		merged_alleleCounts = "analysis/mutant.alleleCounts.simpleIndels.dpthFilt.biallelic.universal.merged.out",
		merged_novel = "analysis/mutant.vs_dm6.bwaUniq.alleleCounts.simpleIndels.dpthFilt.biallelic.universal.novel.merged.count",
		merged_back = "analysis/mutant.vs_dm6.bwaUniq.alleleCounts.simpleIndels.dpthFilt.biallelic.universal.back.merged.count",
	params:
		runmem_gb=8,
		runtime="10:00",
		cores=4,
		chroms1 = "chr3L,chr3R",
		chroms2 = "chr2L,chr2R,chr4,chrX",
	message:
		"interleaving the separately called deficiency and background results into a consolidated mutant analysis... "
	run:
		grep1 = "\|".join(params.chroms1.split(","))
		grep2 = "\|".join(params.chroms2.split(","))
		shell(
			"""
			cat {input.alleleCounts_out[0]} | grep -w "{grep1}" > {output.merged_alleleCounts}
			cat {input.alleleCounts_out[1]} | grep -w "{grep2}" >> {output.merged_alleleCounts}

			cat {input.novel_count[0]} | grep -w "{grep1}" > {output.merged_novel}
			cat {input.novel_count[1]} | grep -w "{grep2}" >> {output.merged_novel}

			cat {input.back_count[0]} | grep -w "{grep1}" > {output.merged_back}
			cat {input.back_count[1]} | grep -w "{grep2}" >> {output.merged_back}
			"""
			)




rule write_report:
	input:
		reference_genome_summary = ["meta/reference_genomes.summary"],
		sequenced_reads_summary=["meta/sequenced_reads.dat"],
		alignment_summaries = expand("meta/alignments.vs_dm6.{aligner}.summary", aligner=['bwa','bwaUniq']),
		VCF_reports = "meta/allGroups.vs_dm6.bwaUniq.calledVariants.summary", #expand("meta/{treat}.vs_dm6.calledVariants.summary", treat="allGroups"),#treat=["control","mutant"]),
		winnowed_VCF_counts = expand("analysis/{prefix}.alleleCounts.simpleIndels.dpthFilt.biallelic.universal.out", prefix=["control","background","deficiency"]),
		merged_winnowed_count = ["analysis/mutant.alleleCounts.simpleIndels.dpthFilt.biallelic.universal.merged.out"],
		novelist_counts = expand("analysis/{prefix}.vs_dm6.bwaUniq.alleleCounts.simpleIndels.dpthFilt.biallelic.universal.{origin}.count", prefix=["control","background","deficiency"], origin=["novel"]),
		merged_novelist_counts = ["analysis/mutant.vs_dm6.bwaUniq.alleleCounts.simpleIndels.dpthFilt.biallelic.universal.novel.merged.count"],
		#VCF_contrasts = ["meta/VCFs/all_samples__bwa.contrast.all_samples__bwaUniq.dm6.diff.sites_in_files", "meta/VCFs/all_samples__bwa.contrast.all_samples__bwaUniq.dm6.diff.sites", "meta/VCFs/all_samples__bwa.contrast.all_samples__bwaUniq.dm6.diff.sites.1", "meta/VCFs/all_samples__bwa.contrast.all_samples__bwaUniq.dm6.diff.sites.2", "meta/VCFs/all_samples__bwa.contrast.all_samples__bwaUniq.dm6.disc", "meta/VCFs/all_samples__bwa.contrast.all_samples__bwaUniq.dm6.disc.count",],
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




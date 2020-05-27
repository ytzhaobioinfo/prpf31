wd = os.getcwd()
SAMPLES = ["sib1","sib2","sib3","prpf311","prpf312","prpf313"]
fastqc_out = wd + "/FastQC/"
genomeFa = "~/yolk/db/Danio_rerio.GRCz11.dna.primary_assembly.fa"
genomeGtf = "~/yolk/db/Danio_rerio.GRCz11.93.chr.gtf"
threadS = 18
stranded = 0 # 0 for unstranded
strand_command = "--outFilterIntronMotifs RemoveNoncanonical"
rRNA_strand_command = "--outFilterIntronMotifs RemoveNoncanonical"
genomeDir = "~/yolk/db/chrAdd_GRCZ11_93"

rule all:
    input:
        expand("{sample}.ok", sample = SAMPLES),
        expand("{sample}.OK", sample = SAMPLES),
        expand("STAR/{sample}/{sample}.sorted.bam", sample = SAMPLES),
        expand("STAR/{sample}/{sample}.counts.tab", sample = SAMPLES),
        expand("STAR/{sample}/{sample}.Log.final.out", sample = SAMPLES),
        expand("STAR/{sample}/{sample}.SJ.out.tab", sample = SAMPLES),
        expand("STAR/{sample}/{sample}.sorted.bam.bai", sample = SAMPLES),
        expand("STAR/STAR_Align_Report.csv", sample = SAMPLES),
        expand("STAR/STAR_Align_Report.png", sample = SAMPLES),
        expand("STAR/STAR_Gene_Counts.csv", sample = SAMPLES),

rule fastqc:
    input:
        R1 = "{sample}_R1.fq",
        R2 = "{sample}_R2.fq"
    output:
        finished = "{sample}.ok"
    log:
        out = fastqc_out + "{sample}.fastqc_fisrt.stdout",
        err = fastqc_out + "{sample}.fastqc_fisrt.stderr"
    threads: threadS
    shell:
        """
        fastqc -o {fastqc_out} --threads {threads} {input.R1} 2> {log.err} 1> {log.out}
        fastqc -o {fastqc_out} --threads {threads} {input.R2} 2>> {log.err} 1>> {log.out}
        touch {output.finished}
        """

rule trimming:
    input:
        R1 = rules.fastqc.input.R1,
        R2 = rules.fastqc.input.R2,
        OK = "{sample}.ok",
    output:
        cutadapted_R1 = "{sample}_adaptor_removed_R1.fq",
        cutadapted_R2 = "{sample}_adaptor_removed_R2.fq",
        trimmed_R1 = "{sample}_Trimed_R1.fq",
        trimmed_R2 = "{sample}_Trimed_R2.fq",
        uniqued_R1 = "{sample}_Uniqued_R1.fq",
        uniqued_R2 = "{sample}_Uniqued_R2.fq"
    threads: threadS
    shell:
        """
        cutadapt -j {threads} -q 30,30 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -e 0.1 -O 5 -m 50 --discard-trimmed --pair-filter=any -o {output.cutadapted_R1} -p {output.cutadapted_R2} {input.R1} {input.R2}
        fastx_trimmer -i {output.cutadapted_R1} -f 10 -o {output.trimmed_R1}&
        fastx_trimmer -i {output.cutadapted_R2} -f 10 -o {output.trimmed_R2}&
        wait
        nohup fastx_collapser -v -i {output.trimmed_R1} -o {output.uniqued_R1}&
        nohup fastx_collapser -v -i {output.trimmed_R2} -o {output.uniqued_R2}&
        """

rule fastqc_2:
    input:
        R1 = rules.trimming.output.trimmed_R1,
        R2 = rules.trimming.output.trimmed_R2
    output:
        end = "{sample}.OK"
    log:
        out = fastqc_out + "{sample}.fastqc_second.stdout",
        err = fastqc_out + "{sample}.fastqc_second.stderr"
    threads: threadS
    shell:
        """
        fastqc -o {fastqc_out} --threads {threads} {input.R1} 2> {log.err} 1> {log.out}
        fastqc -o {fastqc_out} --threads {threads} {input.R2} 2>> {log.err} 1>> {log.out}
        touch {output.end}
        """

rule run_STAR:
    input:
        R1 = "{sample}_Trimed_R1.fq",
        R2 = "{sample}_Trimed_R2.fq",
    output:
        bam = protected("STAR/{sample}/{sample}.sorted.bam"),
        counts = "STAR/{sample}/{sample}.counts.tab",
        log_file = "STAR/{sample}/{sample}.Log.final.out",
        sjtab = "STAR/{sample}/{sample}.SJ.out.tab",
    params:
        genomeDir = genomeDir,
        stranded = strand_command,
        prefix = lambda wildcards: "STAR/{sample}/{sample}".format(sample = wildcards.sample),
        readgroup = lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample = wildcards.sample),
    threads: 16
    message: "Running STAR Alignment on {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}/{sample}.run_STAR.txt"
    shell:
        "STAR --runMode alignReads --runThreadN {threads}"
        " --genomeDir {params.genomeDir}"
        " --readFilesIn {input.R1} {input.R2} "
        " --outFileNamePrefix {params.prefix}."
        " --outSAMstrandField intronMotif"
        " --outSAMmode Full --outSAMattributes All"
        " --outSAMattrRGline {params.readgroup}"
        " --outSAMtype BAM SortedByCoordinate Unsorted"
        " --limitBAMsortRAM 56000000000"
        " --quantMode GeneCounts"
        " --outReadsUnmapped Fastx"
        " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"
        " && mv {params.prefix}.ReadsPerGene.out.tab {output.counts}"

rule index_bam:
    """INDEX the {sample}.sorted.bam file"""
    input:
        "STAR/{sample}/{sample}.sorted.bam"
    output:
        "STAR/{sample}/{sample}.sorted.bam.bai"
    message: "Indexing {wildcards.sample}.sorted.bam"
    benchmark:
        "{sample}/{sample}.index_bam.txt"
    shell:
        "samtools index {input}"

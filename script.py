#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import shutil
import pybedtools
import re
import subprocess
import tempfile

from pysam import samtools

def parse_attributes(attr_str):
    """Parse GFF attributes into a dictionary."""
    attrs = {}
    for attr in attr_str.split(";"):
        attr = attr.strip()
        if attr and '=' in attr:
            key, value = attr.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs

def get_gene(gff_path):
    """Generate a set of genic regions from a GFF file with coordinate positions, strand, and gene ID."""
    try:
        cds_gene = set()
        gff = pybedtools.BedTool(gff_path).filter(lambda x: x[2] == "gene")
        for feature in gff:
            gene_id = feature.attrs.get("ID", '')
            if gene_id:
                cds_gene.add((feature.chrom, feature.start, feature.end, feature.strand, gene_id))
        return cds_gene
    except Exception as e:
        return {'status': f'Failed to process {gff_path}'}

def get_genome_file(genome,output_dir,samtools_path='samtools', genome_prefix=''):
    """Generate a tab-separated file of sequence ID and length using samtools faidx command."""
    try:
        fai_file = f"{genome}.fai"
        genome_file = os.path.join(output_dir, f"{genome_prefix}_genome_file.txt")
        # Run samtools faidx to generate {genome}.fai
        subprocess.run(f"{samtools_path} faidx {genome}", shell=True, check=True)
        # Read the first two column from {genome}.fai
        with open(fai_file, 'r') as fai, open(genome_file, 'w') as out:
            for line in open(fai_file):
                fields = line.strip().split("\t")
                try:
                    out.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\n")
                except IndexError:
                    pass
        return {'status': f'Successfully generated {genome_file}'}
    except subprocess.CalledProcessError as e:
        return {'status': f'Failed to generate the: {e}', 'genome_file': None}

class ToolChecker:
    def __init__(self, tool_paths, picard_jar):
        """Initialize with dictionary of tool paths and Picard jar path"""
        self.tool_paths = tool_paths
        self.picard_jar = picard_jar
        self.tool = ['bwa-mem2','samtools','bedtools','qualimap','bamCoverage','nucmer']
        self.report = []

    def check_tool(self):
        """Check tools availability, prioritizing PATH, generate stdout report."""
        missing_tools = []
        if self.picard_jar and os.path.exists(self.picard_jar):
            self.report.append({'tool':'picard_jar', 'status': f'Found: {self.picard_jar}'})
        else:
            picard_path = shutil.which('picard')
            if picard_path:
                self.report.append({'tool':'picard', 'status': f'Found picard in PATH: {picard_path}'})
            else:
                self.report.append({'tool': 'picard_jar', 'status': f'Failed to find picard in the PATH'})
        for tool in self.tool:
            path = shutil.which(tool) or self.tool_paths.get(tool)
            if path:
                self.report.append({'tool': tool, 'status': f'Found {tool} in PATH: {path}'})
            else:
                self.report.append({'tool': tool, 'status': f'Failed to find {tool} in PATH'})
                missing_tools.append(tool)
        if missing_tools:
            print("Missing tools:", ",".join(missing_tools))
        else:
            print("All tools are installed")
        return pd.DataFrame(self.report)

class GFF2BED:
    def __init__(self, gff_file):
        """Initialize with GFF file."""
        self.gff_file = gff_file
    def convert(self):
        """Convert GFF to BED format using get_gene."""
        return get_gene(self.gff_file)

class ReadMapping:
    def __init__(self, query_genome,query_genome_prefix,fread,rread,output_dir,cpus=28,bwa_mem2_path=None, samtools_path=None):
        """Initialize with query genome and query genome prefix, fread, read, and output directory."""
        self.query_genome = query_genome
        self.query_genome_prefix = query_genome_prefix
        self.fread = fread
        self.rread = rread
        self.output_dir = output_dir
        self.cpus = cpus
        self.bwa_mem2_path = shutil.which(bwa_mem2_path) or bwa_mem2_path or 'bwa-mem2'
        self.samtools_path = shutil.which(samtools_path) or samtools_path or 'samtools'
        self.out = f"{self.query_genome_prefix}_{os.path.basename(self.fread).split('_')[0]}"

    def index_genome(self):
        """Execute terminal command to index genome."""
        try:
            subprocess.run(f"{self.samtools_path} index {self.query_genome}", shell=True, check=True)
            subprocess.run(f"{self.bwa_mem2_path} index {self.query_genome}", shell=True, check=True)
            return {'status': f'{self.query_genome} indexed'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to index genome: {e}'}

    def map_reads(self):
        """Execute terminal command to map reads."""
        try:
            cmd = f"{self.bwa_mem2_path} mem -t {self.cpus} {self.query_genome} {self.fread} {self.rread} | {self.samtools_path} view - -Sb -@{self.cpus} | {self.samtools_path} view -b -@{self.cpus} -F 4 | {self.samtools_path} sort - -@{self.cpus} -o {self.output_dir}/{self.out}.allMapped.sorted.bam"
            subprocess.run(cmd, shell=True, check=True)
            return {'status' : f'Mapped reads: {self.output_dir}/{self.out}.allMapped.sorted.bam'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to map reads: {e}'}

    def remove_duplicates(self, memory_size, picard_jar):
        """Execute terminal command to remove duplicate reads."""
        bam_file = f"{self.output_dir}/{self.out}.allMapped.sorted.bam"
        try:
            if picard_jar and os.path.exists(picard_jar):
                cmd = f"java -Xmx{memory_size} -jar {picard_jar} MarkDuplicate INPUT={bam_file} O={self.output_dir}/{self.out}.allMapped.sorted.markdup.bam M={self.output_dir}/{self.out}.allMapped.sorted.markdup.bam.metrics.txt REMOVE_DUPLICATES=True"
                subprocess.run(cmd, shell=True, check=True)
                return {'status' : 'PCR duplicates removed with Picard Duplicates'}
            else:
                print('Picard jar is not provided or missing, using SAMtools for duplicates removal')
                cmd0 = f"{self.samtools_path} collate -o {self.output_dir}/{self.out}.allMapped.sorted.collate.bam {bam_file} -@{self.cpus}"
                subprocess.run(cmd0, shell=True, check=True)
                os.remove(bam_file)
                cmd1 = f"{self.samtools_path} fixmate -m {self.output_dir}/{self.out}.allMapped.sorted.collate.bam {self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam -@{self.cpus}"
                subprocess.run(cmd1, shell=True, check=True)
                os.remove(f'{self.output_dir}/{self.out}.allMapped.sorted.collate.bam')
                cmd2 = f"{self.samtools_path} sort -o {self.output_dir}/{self.out}.allMapped.sorted.02.bam {self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam -@{self.cpus}"
                subprocess.run(cmd2, shell=True, check=True)
                os.remove(f'{self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam')
                cmd3 = f"{self.samtools_path} markdup {self.output_dir}/{self.out}.allMapped.sorted.02.bam {self.output_dir}/{self.out}.allMapped.sorted.markdup.bam -@{self.cpus}"
                subprocess.run(cmd3, shell=True, check=True)
                os.remove(f"{self.output_dir}/{self.out}.allMapped.sorted.02.bam")
                return {'status' : 'Duplicates removed with SAMtools'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to remove duplicates: {e}'}

    def index_bam(self):
        """Execute terminal command to index bam file."""
        bam_file = f"{self.output_dir}/{self.out}.allMapped.sorted.markdup.bam"
        try:
            cmd = f"{self.samtools_path} index {bam_file} -@{self.cpus}"
            subprocess.run(cmd, shell=True, check=True)
            return {'status' : f'Indexed bam file: {bam_file}'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to index bam file: {e}'}

class CoverageAnalysis:
    def __init__(self, bam_file, output_dir, cds_gene,sra_ids,query_genome_prefix,cpus=28,bedtools_path=None, qualimap_path=None,bamcoverage_path=None):
        self.bam_file = bam_file
        self.output_dir = output_dir
        self.sra_ids = sra_ids
        self.cds_gene = cds_gene
        self.query_genome_prefix = query_genome_prefix
        # self.genome_file = genome_file
        self.cpus = cpus
        self.bamcoverage_path = shutil.which(bamcoverage_path) or bamcoverage_path or 'bamcoverage'
        self.qualimap_path = shutil.which(qualimap_path) or qualimap_path or 'qualimap'
        self.bedtools_path = shutil.which(bedtools_path) or bedtools_path or 'bedtools'

    def run_qualimap(self, memory_size):
        """Execute qualimap command."""
        qualimap_dir = os.path.join(self.output_dir, 'qualimap')
        qualimap_by_region_dir = os.path.join(self.output_dir, 'qualimap_by_regions')
        try:
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_bed:
                for chr, start, end, strand, gene_id in self.cds_gene:
                    temp_bed.write(f"{chr}\t{start}\t{end}\t{strand}\t{gene_id}\n")
                temp_bed_path = temp_bed.name
            cmd0 = f"{self.qualimap_path} bamqc -bam {self.bam_file} -output {qualimap_dir} -outfile {self.sra_ids}_{self.query_genome_prefix}.qualimap -sd -c -nt {self.cpus} -outformat PDF:HTML ip --jave-mem-size={memory_size}"
            subprocess.run(cmd0, shell=True, check=True)
            os.rename(os.path.join(qualimap_dir, 'genome_results.txt'), os.path.join(qualimap_dir, f"{self.sra_ids}_{self.query_genome_prefix}.genome_results.txt"))
            cmd1 = f"{self.qualimap_path} bamqc -bam {self.bam_file} -output {qualimap_by_region_dir}.qualimap_genes -outfile {self.sra_ids}_{self.query_genome_prefix}.qualimap -sd -c -nt {self.cpus} -gff {temp_bed_path} -oc {self.reads}_{self.query_genome_prefix}.qualimap_genes_cov.txt -os {self.reads}_{self.query_genome_prefix}.qualimap_repeats -outformat PDF:HTML ip --jave-mem-size={memory_size}"
            subprocess.run(cmd1, shell=True, check=True)
            os.rename(os.path.join(qualimap_by_region_dir, 'genome_results.txt'), os.path.join(qualimap_by_region_dir, f"{self.sra_ids}_{self.query_genome_prefix}.genome_results.txt"))
            os.remove(temp_bed_path)
            return {'status' : 'Qualimap anaysis completed'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to run qualimap command: {e}'}

    def generate_bigWigCoverage(self):
        """Execute commands to generate bigWig coverage file."""
        bw_file = os.path.join(self.output_dir, f"{os.path.basename(self.bam_file).replace('.bam', '')}.coverage.bw")
        try:
            cmd = f"{self.bamcoverage_path} -b {self.bam_file} --numberOfProcessors {self.cpus} -o {bw_file}"
            subprocess.run(cmd, shell=True, check=True)
            return {'status' : 'BigWig coverage generated'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to generate bigWig coverage: {e}'}

    def generate_gene_coverage(self, memory_size):
        """Execute commands to generate gene coverage file."""
        cov_file = os.path.join(self.output_dir, f"{self.sra_ids}_{self.query_genome_prefix}.allMapped.reads_gene.cov.bed")
        try:
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_bed:
                for chrom, start, end, strand, gene_id in self.cds_gene:
                    temp_bed.write(f"{chrom}\t{start}\t{end}\t{strand}\t{gene_id}")
                temp_bed_path = temp_bed.name
            cmd = f"{self.bedtools_path} bamtobed -i {self.bam_file} | {self.bedtools_path} coverage -a {temp_bed_path} -iobuf {memory_size} -b - > {cov_file}"
            subprocess.run(cmd, shell=True, check=True)
            return {'status' : 'Gene coverage generated'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to generate gene coverage: {e}'}

class CDSPresenceAbsence:
    def __init__(self, output_dir, query_genome_prefix, query_genome, query_genome_gff, sra_ids, cds_gene, cpus=28, bedtools_path=None):
        self.output_dir = output_dir
        self.query_genome_prefix = query_genome_prefix
        self.query_genome = query_genome
        self.gff_file = query_genome_gff
        self.sra_ids = sra_ids
        self.cpus = cpus
        self.cds_gene = cds_gene
        self.bedtools_path = shutil.which("bedtools") or bedtools_path or "bedtools"

    def extract_readsPAVcds_sequences(self):
        """Extract CDS coordinate as a bed file and a final CDS sequences of those genes whose reads mapping coverage were less than 0.245 or manually given reads coverage value."""
        reads_mapping_CDS_frac_file = os.path.join(self.output_dir, f"{self.sra_ids}_{self.query_genome_prefix}.l25p_reads_cov_gene.bed")
        query_cds_fasta = os.path.join(self.output_dir, f"{self.query_genome_prefix}.l25p_reads_cov_gene.fasta")
        try:
            low_cov_gene = set()
            for sra_id in self.sra_ids:
                cov_file = os.path.join(self.output_dir, f"{self.sra_ids}_{self.query_genome_prefix}.allMapped.reads_gene.cov.bed")
                if os.path.exists(cov_file):
                    cov_data = pd.read_csv(cov_file, sep="\t", header=None)
                    low_cov_gene.update(cov_data[cov_data.iloc[:,-1] < 0.245].iloc[:,3].unique())
                if not low_cov_gene:
                    return {'status': 'No genes with coverage < 0.245 was found'}
            gff = pybedtools.BedTool(self.gff_file).filter(lambda x: x[2]=="gene") and parse_attributes(x.attrs).get("ID",'') in low_cov_gene
            cds_bed = gff.to_dataframe(names=['seqid', 'source', 'type', 'start', 'end' 'score', 'strand', 'phase' 'attributes'])[('seqid', 'start', 'end', 'attributes')]
            cds_bed['gene_id'] = cds_bed['attributes'].apply(lambda x: parse_attributes(x).get("ID",''))
            cds_bed = cds_bed[cds_bed['gene_id'] != ''][['seqid', 'start', 'end', 'gene_id']]
            cds_bed.to_csv(reads_mapping_CDS_frac_file, sep='\t', index=False)
            cmd = f"{self.bedtools_path} getfasta -fi {self.query_genome} -bed {reads_mapping_CDS_frac_file} -fo {query_cds_fasta}"
            subprocess.run(cmd, shell=True, check=True)
            return {'status': 'CDS sequences extracted'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to extract CDS sequences from {e}'}

class CDSAlignment:
    def __init__(self, output_dir, sra_ids,query_genome_prefix, query_genome, ref_genome_prefix, ref_genome,ref_genome_gff,reads_mapping_cds_dir, cpus=28, nucmer_path=None, lastz_path=None, bedtools_path=None, samtools_path=None):
        """Initialize with output directory, reads prefix, ..."""
        self.output_dir = output_dir
        self.sra_ids = sra_ids
        self.query_genome_prefix = query_genome_prefix
        self.query_genome = query_genome
        self.ref_genome_prefix = ref_genome_prefix
        self.ref_genome = ref_genome
        self.ref_genome_gff = ref_genome_gff
        self.reads_mapping_cds_dir = reads_mapping_cds_dir
        self.cpus = cpus
        self.nucmer_path = shutil.which("nucmer") or nucmer_path or "nucmer"
        self.lastz_path = shutil.which("lastz") or lastz_path or "lastz"
        self.bedtools_path = shutil.which("bedtools") or bedtools_path or "bedtools"
        self.samtools_path = shutil.which("samtools") or samtools_path or samtools_path

    def mask_repeats(self,cds_bed,genome):
        self.cds_bed = cds_bed
        self.genome = genome
        try:
            cds_complement = os.path.join(self.output_dir, f"{self.genome}.cds_complement.bed")
            repeats_masked_genome = os.path.join(self.output_dir, f"{self.genome}.repeats_masked.fasta")
            genome_file = get_gene(self.genome)
            if cds_bed and genome_file:
                cmd0 = f"{self.bedtools_path} complement -i {self.cds_bed} -g {genome_file} -o {genome_file} > {cds_complement}"
                subprocess.run(cmd0, shell=True, check=True)
                cmd1 = f"{self.bedtools_path} makefasta -if {self.genome} -bed  {cds_complement} -fo {repeats_masked_genome}"
                subprocess.run(cmd1, shell=True, check=True)
                return {'status': 'Repeats masked genome generated'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to generate genome mask: {e}'}

    def generate_ref_cds_fasta(self):
        try:
            ref_cds_bed = os.path.join(self.output_dir, self.ref_genome_prefix + "_ref_cds.bed")
            ref_cds_fasta = os.path.join(self.output_dir, self.ref_genome_prefix + "_ref_cds.fasta")
            gff_result = get_gene(self.ref_genome_gff)
            if not gff_result['cds_gene']:
                return {'status': f'No genes with CDS gene was found in : {self.ref_genome_prefix}'}
            with open(ref_cds.bed, 'w') as ref_cds_bed_file:
                for chrom, start, end, strand in gff_result['cds_gene']:
                    ref_cds_bed_file.write(f'{chrom}\t{start}\t{end}\t{strand}\t{gene_id}\n')
            cmd = f"{self.bedtools_path} getfatsa -fi {ref_cds_bed} -bed {ref_cds_bed_file} -fo {ref_cds_fasta}"
            subprocess.run(cmd, shell=True, check=True)
            return {'status': f'CDS sequences extracted for: {self.ref_genome_prefix}'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to generate reference CDS fasta: {e}'}

    def align_nucmer(self, ref_genome_prefix):
        query_cds_fasta = os.path.join(self.output_dir, f"{self.query_genome_prefix}.l25p_reads_cov_gene.fasta")
        ref_cds_fasta = os.path.join(self.output_dir, f"{self.ref_genome_prefix}_ref_cds.fasta")
        prefix = f"{self.query_genome_prefix}_{self.ref_genome_prefix}"
        delta_file = os.path.join(self.output_dir, f"{prefix}.delta")
        coords_file = os.path.join(self.output_dir, f"{prefix}.coords")
        ref_nucmer_bed = os.path.join(self.output_dir, f"{prefix}.nucmer.bed")
        nucmer_sorted_bed = os.path.join(self.output_dir, f"{prefix}.nucmer.sorted.bed")
        ref_specific_bed = os.path.join(self.output_dir, f"{prefix}.ref_specific.cds.bed")
        query_nucmer_bed = os.path.join(self.output_dir, f"{prefix}.query_nucmer.bed")
        query_nucmer_sorted_bed = os.path.join(self.output_dir, f"{prefix}.query_nucmer.sorted.bed")
        query_specific_bed = os.path.join(self.output_dir, f"{prefix}.query_specific.cds.bed")
        try:
            cmd0 = f"{self.nucmer_path} --prefix {prefix} {ref_cds_fasta} {query_cds_fasta} --batch 1 --threads {self.cpus}"
            subprocess.run(cmd0, shell=True, check=True)
            cmd1 = f"show-coords -c -d -l -r -o -T {delta_file} > {coords_file}"
            subprocess.run(cmd1, shell=True, check=True)
            # Convert coords to BED (reference)
            cmd2 = f"cat {coords_file} | grep -v '=\\|/\\|NUCMER' | sed 's/|//g' | grep -v \"\\[\" | awk '{{print $14,$1,$2,$5,$7,$10,$8}}' | grep \"[a-zA-Z]\" | awk '{{ if ($2>$3) print $1,$3,$2,\"-\", $4,$5,$6,$7; else print $1,$2,$3,\"+\", $4,$5,$6,$7}}' OFS='\\t' > {ref_nucmer_bed}"
            subprocess.run(cmd2, shell=True, check=True)
            # Sort and merge reference BED
            cmd3 = f"{self.bedtools_path} sort -i {ref_nucmer_bed} | {self.bedtools_path} merge > {nucmer_sorted_bed}"
            subprocess.run(cmd3, shell=True, check=True)
            # Intersect for reference-specific CDS
            cmd4 = f"{self.bedtools_path} intersect -a {self.reads_mapping_cds_dir}/EV_mazia_panma_cov.sorted.merged.bed -b {nucmer_sorted_bed} -wao | awk '{{if ($7/($3-$2+1)*100 < 24.4) print $1,$2,$3,$4,$5,$6,$7,$7/($3-$2)*100}}' OFS='\\t' > {ref_specific_bed}"
            subprocess.run(cmd4, shell=True, check=True)
            # Convert coords to BED (query)
            cmd5 = f"cat {coords_file} | grep -v '=\\|/\\|NUCMER' | sed 's/|//g' | grep -v \"\\[\" | awk '{{print $15,$3,$4,$6,$7,$11,$9}}' | grep \"[a-zA-Z]\" | awk '{{ if ($2>$3) print $1,$3,$2,\"-\", $4,$5,$6,$7; else print $1,$2,$3,\"+\", $4,$5,$6,$7}}' OFS='\\t' > {query_nucmer_bed}"
            subprocess.run(cmd5, shell=True, check=True)
            # Sort and merge query BED
            cmd6 = f"{self.bedtools_path} sort -i {query_nucmer_bed} | {self.bedtools_path} merge | sed 's/:/\\t/g' | sed 's/-/\t/g' | awk '{{print $1, $2+$4-1,$2+$5-1}}' OFS='\\t' | {self.bedtools_path} sort | {self.bedtools_path} merge > {query_nucmer_sorted_bed}"
            subprocess.run(cmd6, shell=True, check=True)
            # Intersect for query-specific CDS
            cmd7 = f"{self.bedtools_path} intersect -a {self.reads_mapping_cds_dir}/MA_panev_cov.sorted.merged.bed -b {query_nucmer_sorted_bed} -wao | awk '{{if ($7/($3-$2+1)*100 < 25) print $1,$2,$3,$4,$5,$6,$7,$7/($3-$2)*100}}' OFS='\\t' > {query_specific_bed}"
            subprocess.run(cmd7, shell=True, check=True)
            return {'status': f'Nucmer alignment completed: {ref_specific_bed}, {query_specific_bed}'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed Nucmer alignment: {e}'}

    def align_lastz(self):


























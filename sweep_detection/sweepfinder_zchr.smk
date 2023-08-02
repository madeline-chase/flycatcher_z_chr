## Snakemake workflow to run recombination rate estimation workflow using LDhelmet
configfile: "config.yaml"
localrules: all

import pandas as pd

scaffolds = pd.read_table(config['scaff_file'], names = ['chr_name','scaff_name', 'category'])

rule all:
    input:
        sf_out = expand('results/zscaffs/sf_out/{spec}/{spec}.{scaff}.fem_adjusted.sf2.ug.anc_fix_rmvd.out.txt', spec = config['spec'], scaff = scaffolds.scaff_name[scaffolds['category']=='z_scaff'])

## Get frequency data for polarized sites for all species
rule get_counts:
    input:
        vcf = '../allsites_vc/data/vcf/z_scaffs/missing_filt/all_sites.{scaff}.indel_rmvd.max_allele2.hard_filt.pass_only.RM.CM.minGQ.minDP.maxDP.max_miss_10_perc.vcf.gz',
        samples = '../allsites_vc/data/samples/{spec}.{sex}.txt',
        polar_sites = '../allsites_vc/data/counts/polar_sites/no_cpg_filt/polar_sites.{scaff}.all_filt.missing_rmvd.z_scaffs.var_only_counts.groups_combined.txt'
    output:
        counts = 'intermediate/zscaffs/counts/{spec}.{scaff}.{sex}.polar_sites.counts.frq.count'
    params:
        prefix = 'intermediate/zscaffs/counts/{spec}.{scaff}.{sex}.polar_sites.counts'
    shell:
        """
        vcftools --gzvcf {input.vcf} --counts --keep {input.samples} --positions {input.polar_sites} --out {params.prefix}
        """

## Split count output data
rule split_counts:
    input:
        counts = 'intermediate/zscaffs/counts/{spec}.{scaff}.{sex}.polar_sites.counts.frq.count'
    output:
        counts_split = 'intermediate/zscaffs/counts/{spec}.{scaff}.{sex}.polar_sites.counts.split.txt'
    shell:
        """
        awk -v OFS='\t' '{{if(NR>1) {{split($5,a,":") split($6,b,":"); print $1,$2,$3,$4,a[1],a[2],b[1],b[2]}}}}' {input.counts} > {output.counts_split}
        """

rule combine_male_fem_counts:
    input:
        fem_counts = 'intermediate/zscaffs/counts/{spec}.{scaff}.fem.polar_sites.counts.split.txt',
        male_counts = 'intermediate/zscaffs/counts/{spec}.{scaff}.male.polar_sites.counts.split.txt'
    output:
        counts_combined = 'intermediate/zscaffs/counts/{spec}.{scaff}.fem_adjusted.polar_sites.counts.split.txt'
    shell:
        """
        paste {input.fem_counts} {input.male_counts} | awk -v OFS='\t' '{{fem_p=$6/2; fem_q=$8/2; fem_total=$4/2; print $1,$2,$3,fem_total + $12,$5,fem_p+$14,$7,fem_q+$16}}' > {output.counts_combined}
        """


rule get_sf_data:
    input:
        polar_sites = '../allsites_vc/data/counts/polar_sites/no_cpg_filt/polar_sites.{scaff}.all_filt.missing_rmvd.z_scaffs.var_only_counts.groups_combined.txt',
        split_counts = 'intermediate/zscaffs/counts/{spec}.{scaff}.fem_adjusted.polar_sites.counts.split.txt',
        sf_header = 'data/sf_header.txt'
    output:
        sf_count_input = 'intermediate/zscaffs/counts/sf_in/{spec}.{scaff}.fem_adjusted.polar_sites.sf2_format.anc_fix_rmvd.txt'
    shell:
        """
        paste {input.polar_sites} {input.split_counts} | awk -v OFS='\t' '{{if($3==$8) print $5,$11,$7,"0"; else if($3==$10) print $5,$9,$7,"0"}}' | awk '{{if($2>0) print}}' | cat {input.sf_header} -  > {output.sf_count_input}
        """

rule get_grid_file:
    input:
        sf_count_input = 'intermediate/zscaffs/counts/sf_in/{spec}.{scaff}.fem_adjusted.polar_sites.sf2_format.anc_fix_rmvd.txt'
    output:
        sf_grid_file = 'intermediate/zscaffs/counts/sf_in/{spec}.{scaff}.fem_adjusted.polar_sites.sf2_format.anc_fix_rmvd.grid_file.txt'
    run:
        if len(pd.read_csv(input.sf_count_input)) > 0:
            shell("grep -v 'folded' {input.sf_count_input} | awk '{{print $1}}' > {output.sf_grid_file}")
        else:
            shell("touch {output.sf_grid_file}")

rule get_whole_genome_input:
    input:
        counts = expand('intermediate/zscaffs/counts/sf_in/{{spec}}.{scaff}.fem_adjusted.polar_sites.sf2_format.anc_fix_rmvd.txt', scaff = scaffolds.scaff_name[scaffolds['category']=='z_scaff']),
        sf_header = 'data/sf_header.txt'
    output:
        whole_gen_counts = 'intermediate/zscaffs/counts/sf_in/{spec}.whole_genome.fem_adjusted.polar_sites.sf2_format.anc_fix_rmvd.txt'
    shell:
        """
        cat {input.counts} | grep -v 'folded' | cat {input.sf_header} - > {output.whole_gen_counts}
        """

rule get_whole_genome_sf:
    input:
        whole_gen_counts = 'intermediate/zscaffs/counts/sf_in/{spec}.whole_genome.fem_adjusted.polar_sites.sf2_format.anc_fix_rmvd.txt'
    output:
        whole_gen_spec = 'intermediate/zscaffs/sf_spec/{spec}.whole_genome.fem_adjusted.polar_sites.sf2_format.anc_fix_rmvd.sf_spec'
    shell:
        """
        scripts/SF2/SweepFinder2 -f {input.whole_gen_counts} {output.whole_gen_spec}
        """

rule run_sf:
    input:
        sf_counts = 'intermediate/zscaffs/counts/sf_in/{spec}.{scaff}.fem_adjusted.polar_sites.sf2_format.anc_fix_rmvd.txt',
        whole_gen_spec = 'intermediate/zscaffs/sf_spec/{spec}.whole_genome.fem_adjusted.polar_sites.sf2_format.anc_fix_rmvd.sf_spec',
        sf_grid_file = 'intermediate/zscaffs/counts/sf_in/{spec}.{scaff}.fem_adjusted.polar_sites.sf2_format.anc_fix_rmvd.grid_file.txt'
    output:
        sf_out = 'results/zscaffs/sf_out/{spec}/{spec}.{scaff}.fem_adjusted.sf2.ug.anc_fix_rmvd.out.txt'
    run:
        if len(pd.read_csv(input.sf_counts)) > 15:
            shell("scripts/SF2/SweepFinder2 -lu {input.sf_grid_file} {input.sf_counts} {input.whole_gen_spec} {output.sf_out}")
        else:
            shell("touch {output.sf_out}")


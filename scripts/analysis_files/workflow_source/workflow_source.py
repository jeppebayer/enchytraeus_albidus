#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def poolsnp_vcf_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Create :format:`mpileup`, :format:`sync` and :format:`VCF` files using :script:`samtools`, :script:`popoolation2` and :script:`PoolSNP`.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    config = yaml.safe_load(open(config_file))
    ACCOUNT: str = config['account']
    SPECIES_NAME: str = config['species_name']
    SAMPLE_LIST: list = config['sample_list']
    sample_names: list = [os.path.basename(os.path.dirname(path)) for path in SAMPLE_LIST]
    REFERENCE_GENOME: str = config['reference_genome_path']
    WORKING_DIR: str = config['working_directory_path']
    OUTPUT_DIR: str = config['output_directory_path']
    ALLSITES: bool = config['all_sites']
    MAXCOV: float | list = config['max_cov']
    MINCOV: int | list = config['min_cov']
    MINCOUNT: int | list = config['min_count']
    MINFREQ: float | list = config['min_freq']
    MISSFRAC: float | list = config['miss_frac']
    BASEQUAL: int | list = config['bq']
    
    partition_size = 200000

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    if os.path.exists(f'reference_partitions.{partition_size}bp.txt') and os.path.exists('reference_sequences.txt'):
        # If files exists reads data directly from files
        # Loads reference genome partitioning
        with open(f'reference_partitions.{partition_size}bp.txt', 'r') as infile:
            partitions = [{'num': entry.split(sep='\t')[0].strip(), 'region': entry.split(sep='\t')[1].strip(), 'start': entry.split(sep='\t')[2].strip(), 'end': entry.split(sep='\t')[3].strip()} for entry in infile]
            npadding = len(str(sum(1 for line in partitions)))
        # Loads list of contigs in reference genome
        with open('reference_sequences.txt', 'r') as infile:
            contigs = [{'contig': entry.split(sep='\t')[0].strip()} for entry in infile]
    else:
        # If files don't exist, generate data and write files
        sequences = parse_fasta(REFERENCE_GENOME)
        with open(f'reference_sequences.txt', 'w') as outfile:
            outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))
        # Partitions reference genome
        npadding = padding_calculator(parse_fasta=sequences, size=partition_size)
        partitions = partition_chrom(parse_fasta=sequences, size=partition_size, npad=npadding)
        with open(f'reference_partitions.{partition_size}bp.txt', 'w') as outfile:
            outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in partitions))
        # Creates list of contigs in reference genome
        contigs = [{'contig': contig['sequence_name']} for contig in sequences]

    top_dir = f'{WORKING_DIR}/{SPECIES_NAME.replace(" ", "_")}/analysis_files'
    output_dir = f'{OUTPUT_DIR}/{species_abbreviation(SPECIES_NAME)}'
    
    # Creates list of dictionaries with contigs and corresponding mpileup files
    mpileup_filelist = mpileup_partitions_filelist(partitions=partitions, top_dir=top_dir, species_name=SPECIES_NAME, npad=npadding)

    mpileup = gwf.map(
        name=name_mpileup,
        template_func=mpileup_parts,
        inputs=partitions,
        extra={'bam_files': SAMPLE_LIST,
            'reference_genome_file': REFERENCE_GENOME,
            'species_name': SPECIES_NAME,
            'output_directory': top_dir}
    )

    sync = gwf.map(
        name=name_sync,
        template_func=mpileup2sync,
        inputs=collect(mpileup.outputs, ['mpileup'])['mpileups'],
        extra={'output_directory': top_dir}
    )

    concat_mpileup = gwf.target_from_template(
        name='concatenate_mpileup',
        template=concat(
            files=collect(mpileup.outputs, ['mpileup'])['mpileups'],
            output_name=f'{species_abbreviation(SPECIES_NAME)}',
            output_directory=output_dir
        )
    )

    concat_sync = gwf.target_from_template(
        name='concatenate_sync',
        template=concat(
            files=collect(sync.outputs, ['sync'])['syncs'],
            output_name=f'{species_abbreviation(SPECIES_NAME)}',
            output_directory=output_dir
        )
    )

    if ALLSITES == 1:
        sitestate = 'allsites'
    elif ALLSITES == 0:
        sitestate = 'variants'
    if type(MAXCOV) is list:
        prev = None
        for i in range(0, len(MAXCOV)):
            if MAXCOV[i] != prev:

                coverage_threshold = gwf.map(
                    name=name_cov,
                    template_func=max_cov_threshold,
                    inputs=mpileup_filelist,
                    extra={'cutoff': MAXCOV[i],
                        'output_directory': top_dir}
                )

                concat_coverage = gwf.target_from_template(
                    name=f'concatenate_coverage_{MAXCOV[i]}',
                    template=concat(
                        files=collect(coverage_threshold.outputs, ['cutoff'])['cutoffs'],
                        output_name=f'{species_abbreviation(SPECIES_NAME)}-cov-{MAXCOV[i]}',
                        output_directory=top_dir
                    )
                )

            run_poolsnp = gwf.target_from_template(
                name=f'poolsnp_{sitestate}_maxcov{MAXCOV[i]}_mincov{MINCOV[i]}_mincnt{MINCOUNT[i]}_minfrq{MINFREQ[i]}_missfrc{MISSFRAC[i]}_bq{BASEQUAL[i]}',
                template=poolsnp(
                    mpileup_file=concat_mpileup.outputs['concat_file'],
                    max_cov_file=concat_coverage.outputs['concat_file'],
                    sample_list=sample_names,
                    reference_genome_file=REFERENCE_GENOME,
                    working_directory=top_dir,
                    species_name=SPECIES_NAME,
                    output_directory=f'{output_dir}/vcf_{sitestate}_maxcov{MAXCOV[i]}_mincov{MINCOV[i]}_mincnt{MINCOUNT[i]}_minfrq{MINFREQ[i]}_missfrc{MISSFRAC[i]}_bq{BASEQUAL[i]}',
                    min_cov=MINCOV[i],
                    min_count=MINCOUNT[i],
                    min_freq=MINFREQ[i],
                    miss_frac=MISSFRAC[i],
                    bq=BASEQUAL[i],
                    sites=ALLSITES
                )
            )

            biallelic_vcf = gwf.target_from_template(
                name=f'biallelic_filter_{sitestate}_maxcov{MAXCOV[i]}_mincov{MINCOV[i]}_mincnt{MINCOUNT[i]}_minfrq{MINFREQ[i]}_missfrc{MISSFRAC[i]}_bq{BASEQUAL[i]}',
                template=vcf_filter(
                    vcf_file=run_poolsnp.outputs['vcf'],
                    output_directory=f'{output_dir}/vcf_{sitestate}_maxcov{MAXCOV[i]}_mincov{MINCOV[i]}_mincnt{MINCOUNT[i]}_minfrq{MINFREQ[i]}_missfrc{MISSFRAC[i]}_bq{BASEQUAL[i]}',
                    species_name=SPECIES_NAME
                )
            )

            prev = MAXCOV[i]
        
    elif type(MAXCOV) is float:
        coverage_threshold = gwf.map(
            name=name_cov,
            template_func=max_cov_threshold,
            inputs=mpileup_filelist,
            extra={'cutoff': MAXCOV,
                'output_directory': top_dir}
        )

        concat_coverage = gwf.target_from_template(
            name=f'concatenate_coverage_{MAXCOV}',
            template=concat(
                files=collect(coverage_threshold.outputs, ['cutoff'])['cutoffs'],
                output_name=f'{species_abbreviation(SPECIES_NAME)}-cov-{MAXCOV}',
                output_directory=top_dir
            )
        )

        run_poolsnp = gwf.target_from_template(
            name=f'poolsnp_{sitestate}_maxcov{MAXCOV}_mincov{MINCOV}_mincnt{MINCOUNT}_minfrq{MINFREQ}_missfrc{MISSFRAC}_bq{BASEQUAL}',
            template=poolsnp(
                mpileup_file=concat_mpileup.outputs['concat_file'],
                max_cov_file=concat_coverage.outputs['concat_file'],
                sample_list=sample_names,
                reference_genome_file=REFERENCE_GENOME,
                working_directory=top_dir,
                species_name=SPECIES_NAME,
                output_directory=output_dir,
                min_cov=MINCOV,
                min_count=MINCOUNT,
                min_freq=MINFREQ,
                miss_frac=MISSFRAC,
                bq=BASEQUAL,
                sites=ALLSITES
            )
        )

        biallelic_vcf = gwf.target_from_template(
            name=f'biallelic_filter_{sitestate}_maxcov{MAXCOV}_mincov{MINCOV}_mincnt{MINCOUNT}_minfrq{MINFREQ}_missfrc{MISSFRAC}_bq{BASEQUAL}',
            template=vcf_filter(
                vcf_file=run_poolsnp.outputs['vcf'],
                output_directory=output_dir,
                species_name=SPECIES_NAME
            )
        )

    return gwf
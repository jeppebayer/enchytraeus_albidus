#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def draft_assembly_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Creates draft assembly from PacBio HiFi data.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    CONFIG = yaml.safe_load(open(config_file))
    ACCOUNT: str = CONFIG['account']
    SPECIES_NAME: str = CONFIG['species_name']
    OUTPUT_DIR: str = CONFIG['output_directory_path']
    SETUPS: list = CONFIG['setup_list']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    for setup in SETUPS:
        NAME: str = setup['name']
        HIFI_SEQ: str = setup['pacbio_hifi_sequence_file']
        USE_HIC: int = setup['use_hic_data']
        HIC_READS: list = setup['hic_sequence_files']
        BUSCO_DATASET: str = setup['busco_dataset_name']
        PURGE_ROUNDS: int = setup['purge_dups_rounds']

        top_dir = f'{OUTPUT_DIR}/{SPECIES_NAME.replace(" ", "_")}/{NAME}/draft_assembly'

        remove_adapters = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_adapterremoval',
            template=hifiadapterfilt(
                pacbio_hifi_file=HIFI_SEQ,
                output_directory=top_dir,
            )
        )
        
        # Create draft assembly using only PacBio HiFi data
        if USE_HIC == 2 or USE_HIC == 3:

            hifiasm_primary_assembly = gwf.target_from_template(
                name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_hifiasm_primary',
                template=hifiasm_primary(
                    hifi_sequence_file=remove_adapters.outputs['filt'],
                    output_directory=top_dir,
                    species_name=SPECIES_NAME
                )
            )

            busco_primary = gwf.target_from_template(
            name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_busco_primary',
            template=busco_genome(
                genome_assembly_file=hifiasm_primary_assembly.outputs['fasta'],
                busco_dataset=BUSCO_DATASET
                )
            )

            for i in range(1, PURGE_ROUNDS + 1):
                if i == 1:
                    in_file = hifiasm_primary_assembly.outputs['fasta']
                else:
                    in_file = step3_primary.outputs['purged']
                step1_primary = gwf.target_from_template(
                    name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_primary_pdups_s1_r{i:02}',
                    template=purge_dups_1_map_hifi_to_genome(
                        gemone_assembly_file=in_file,
                        hifi_sequence_file=HIFI_SEQ,
                        output_directory=top_dir,
                        species_name=SPECIES_NAME,
                        directory_addition='primary',
                        round_number=i
                    )
                )
                step2_primary = gwf.target_from_template(
                    name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_primary_pdups_s2_r{i:02}',
                    template=purge_dups_2_map_to_self(
                        genome_assembly_file=in_file,
                        output_directory=top_dir,
                        species_name=SPECIES_NAME,
                        directory_addition='primary',
                        round_number=i
                    )
                )
                step3_primary = gwf.target_from_template(
                    name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_primary_pdups_s3_r{i:02}',
                    template=purge_dups_3_purge_duplicates(
                        pb_stat_file=step1_primary.outputs['stat'],
                        pb_base_cov_file=step1_primary.outputs['cov'],
                        self_alignment_paf=step2_primary.outputs['paf'],
                        genome_assembly_file=in_file,
                        output_directory=top_dir,
                        species_name=SPECIES_NAME,
                        directory_addition='primary',
                        round_number=i
                    )
                )
                purge_busco_primary = gwf.target_from_template(
                    name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_primary_pdups_busco_r{i:02}',
                    template=busco_genome(
                        genome_assembly_file=step3_primary.outputs['purged'],
                        busco_dataset=BUSCO_DATASET
                    )
                )

        # Create draft assembly with PacBio HiFi data and Hi-C data
        if USE_HIC == 1 or USE_HIC == 3:

            hifiasm_hic_assembly = gwf.target_from_template(
                name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_hifiasm_hic',
                template=hifiasm_hic(
                    hifi_sequence_file=HIFI_SEQ,
                    hic_sequence_files=HIC_READS,
                    output_directory=top_dir,
                    species_name=SPECIES_NAME
                )
            )

            busco_hic = gwf.target_from_template(
                name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_busco_hic',
                template=busco_genome(
                    genome_assembly_file=hifiasm_hic_assembly.outputs['fasta'][0],
                    busco_dataset=BUSCO_DATASET
                )
            )

            for i in range(1, PURGE_ROUNDS + 1):
                if i == 1:
                    in_file = hifiasm_hic_assembly.outputs['fasta'][0]
                else:
                    in_file = step3_hic.outputs['purged']
                step1_hic = gwf.target_from_template(
                    name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_hic_pdups_s1_r{i:02}',
                    template=purge_dups_1_map_hifi_to_genome(
                        gemone_assembly_file=in_file,
                        hifi_sequence_file=HIFI_SEQ,
                        output_directory=top_dir,
                        species_name=SPECIES_NAME,
                        directory_addition='hic',
                        round_number=i
                    )
                )
                step2_hic = gwf.target_from_template(
                    name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_hic_pdups_s2_r{i:02}',
                    template=purge_dups_2_map_to_self(
                        genome_assembly_file=in_file,
                        output_directory=top_dir,
                        species_name=SPECIES_NAME,
                        directory_addition='hic',
                        round_number=i
                    )
                )
                step3_hic = gwf.target_from_template(
                    name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_hic_pdups_s3_r{i:02}',
                    template=purge_dups_3_purge_duplicates(
                        pb_stat_file=step1_hic.outputs['stat'],
                        pb_base_cov_file=step1_hic.outputs['cov'],
                        self_alignment_paf=step2_hic.outputs['paf'],
                        genome_assembly_file=in_file,
                        output_directory=top_dir,
                        species_name=SPECIES_NAME,
                        directory_addition='hic',
                        round_number=i
                    )
                )
                purge_busco = gwf.target_from_template(
                    name=f'{species_abbreviation(SPECIES_NAME)}_{NAME}_hic_pdups_busco_r{i:02}',
                    template=busco_genome(
                        genome_assembly_file=step3_hic.outputs['purged'],
                        busco_dataset=BUSCO_DATASET
                    )
                )

    return gwf
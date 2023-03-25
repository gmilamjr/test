#Common functions that are used on different locations

import os
import logging
import pandas as pd
logger = logging.getLogger(__name__)

from os.path import join
from filesplit.split import Split
import shutil

def make_directory(directory):
    try:
        os.mkdir(directory)
        logger.info(f'Created {directory} directory')
    except OSError:
        pass

def get_dfs_from_disk(data_dir):
    all_dfs = []
    for fn in os.listdir(data_dir):
        if fn.endswith('.tsv'):
            all_dfs.append(pd.read_table(os.path.join(data_dir, fn)))

    df = pd.concat(all_dfs, ignore_index=True)
    return df

def save_sequences(fn,seq_list,singleline=False,seqid=None):
    with open(fn,'wt') as f:
        for i in range(len(seq_list)):
            us = seq_list[i]
            if singleline:
                f.write(f"{us}\n")
            elif seqid is None:
                f.write(f"> {us}\n{us}\n")
            else:
                f.write(f"> {seqid}\n{us}\n")
    
#Creation of final master list. Will read data from disk
def create_master_list(taxonomy_id, data_dir):
    logging.info("Getting Peptides")
    peptides_df = pd.read_table(join(data_dir,f'mhci_peptides_{taxonomy_id}.tsv'))
    assert not peptides_df.isnull().values.any(), 'Found NaN values on mhci peptides this should not have occurred'
    peptides_df.head()

    logging.info("Setting number of binding alleles per peptide")
    allele_peptide_agg = peptides_df.groupby('peptide')['allele'].agg('count').reset_index(name="Number of Binding Alleles")
    assert not allele_peptide_agg.isnull().values.any(), 'Found NaN values when grouping peptides for allele counting. This should not have occurred'

    master_list = pd.merge(peptides_df, allele_peptide_agg, how='left', on='peptide')
    assert not master_list.isnull().values.any(), 'Found NaN values after merging count of alleles with MHCI peptides'

    logging.info("Adding protein sequences")
    uniprot_df = pd.read_table(join(data_dir,f'uniprot_protein_sequences_{taxonomy_id}.tsv'), usecols=['protein_id','protein_sequence'])

    master_list = pd.merge(master_list, uniprot_df, how='inner', on=['protein_id'])
    assert not master_list.isnull().values.any(), 'Found NaN values after merging with protein sequences'

    logging.info("Adding immunogenicity")
    immunogenicity_df = pd.read_table(join(data_dir,f'immunogenicity_peptides_{taxonomy_id}.tsv'))
    master_list = pd.merge(master_list, immunogenicity_df[['immunogenicity_score','peptide']], how='inner', on=['peptide'])

    logging.info("Adding antigenicity")
    antigenicity_df = pd.read_table(join(data_dir,f'antigenicity_peptides_{taxonomy_id}.tsv'))
    antigenicity_df.rename(columns={'sequence':'peptide'}, inplace=True)
    master_list = pd.merge(master_list, antigenicity_df, how='inner', on=['peptide'])

    return master_list

def get_rank(df, aff_w=0.3,pro_w=0.1,imm_w=0.3,aff_flurry=1,aff_netmhc=1,aff_netctl=0.2,pro_flurry=1,pro_netmhc=1,vaxijen_w=1,iedb_w=1,prime_w=1,num_all_w=0.3):
    
    #Standardize within 0 - 500 range?
    df['aff_score'] = ((aff_flurry * df['mhcflurry_aff']) + (aff_netmhc * df['iedb_aff'])) / 2#Waiting for netctl confirmation
    
    #Processing score 0 - 1?
    df['pro_score'] = (pro_flurry * df['mhcflurry_processing_score'])#Need to add netmhc processing score waiting for confirmation
    
    #Standardize within -1 to 1 for immunogenicity. ? Antigenicity has a larger range, might need to be dynamic
    df['imm_score'] = ((vaxijen_w * df['antigen_score']) + (iedb_w * df['immunogenicity_score'])) / 2 #Prime waiting for confirmation
    
    df['rank_score'] = (aff_w * df['aff_score']) + (pro_w * df['pro_score']) + (imm_w * df['imm_score']) + df['Number of Binding Alleles']

    df.sort_values('rank_score', inplace=True, ascending=False, ignore_index=False)
    
    return df

def check_split_alleles(allele_filepath, allele_dir):
    
    make_directory(allele_dir)
    
    with open(allele_filepath) as f:
        x = len(f.readlines())
        if x > 8:
            
            #break file down
            split = Split(inputfile = allele_filepath, outputdir = allele_dir)
            split.bylinecount(linecount = 4)
            os.remove(join(allele_dir,'manifest')) #Manifest file not needed by us
        else:
            shutil.copy2(allele_filepath, allele_dir)

    return
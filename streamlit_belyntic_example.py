
# import libraries
import streamlit as st
import pandas as pd
from os.path import join
from utils.utils import create_master_list, get_rank

st.set_page_config(layout="wide")
# add title and instructions
st.title("Ranking of epitope candidates")

#Below to remove the plus minus button on number input
st.markdown("""
<style>
    button.step-up {display: none;}
    button.step-down {display: none;}
    div[data-baseweb] {border-radius: 4px;}
</style>""",
unsafe_allow_html=True)

state = st.session_state

def set_session_states():
    
    if 'aff_w' not in state:
        state['aff_w'] = 0.3
    if 'aff_flurry' not in state:
        state['aff_flurry'] = 1
    if 'aff_netmhc' not in state:
        state['aff_netmhc'] = 1
    if 'pro_w' not in state:
        state['pro_w'] = 0.1
    if 'pro_flurry' not in state:
        state['pro_flurry'] = 1
    if 'pro_netmhc' not in state:
        state['pro_netmhc'] = 1
    if 'imm_w' not in state:
        state['imm_w'] = 0.3
    if 'vaxijen_w' not in state:
        state['vaxijen_w'] = 1 
    if 'iedbimm_w' not in state:
        state['iedbimm_w'] = 1
    if 'prime_w' not in state:
        state['prime_w'] = 1
    if 'num_all_w' not in state:  
        state['num_all_w'] = 0.1
    if 'actions_submit' not in state:
        state['actions_submit'] = False

def _change_callback(data_dir, taxonomy_id, include_assay_results):
    state.actions_submit = True
    
    if 'df' not in state:
        state['df'] = get_df(data_dir, taxonomy_id, include_assay_results)
    

def main(data_dir = 'data'):
    
    set_session_states()
    
    #Use select box for organism that have been already searched for.
    with st.sidebar:
        
        taxonomy_id = st.number_input("Enter taxonomy ID to provide epitope candidates for it: ", 
                                      value = 0,
                                      format = "%d",
                                      help= "https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi")
        
        allele_file = st.file_uploader("Upload file containing alleles and respective lengths that should be predicted for",
                                       help = "E.g. HLA-A*01:01,8")
        
        if taxonomy_id > 0:
            data_dir = join(data_dir, f'{taxonomy_id}') 
            
            pipeline = st.selectbox("Choose which pipeline to run: ", ['ALL', 'MHCI', 'MHCII', None])
            
            if pipeline is None:
                
                mhc_class = st.radio("What type of (MHC) class wish to focus on?", ('MHCI', 'MCHII'))
                if mhc_class == 'MHCI':
                    tools = st.multiselect("Which data source need updating", 
                                         ["UniProt","BV-BRC", "IEDB MHCI", "MHCFLurry", "IEDB Immunogenecity", "Vaxijen", "Prime", "Epitope Ranking"])
                else:
                    tools = st.multiselect("Which data source need updating", 
                                         ["UniProt","BV-BRC", "IEDB MHCII", "Vaxijen", "Epitope Ranking"])
            
            update_rank_dict = {
                "data_dir" : data_dir,
                "taxonomy_id": taxonomy_id,
                "include_assay_results": True                
                }
            
            st.button("Click to perform chosen actions", on_click=_change_callback, kwargs=update_rank_dict)
                
    if state.actions_submit:
        with st.form(key="weights"):
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.number_input("Overall Binding Affinity (BA) Weight", value=state.aff_w, key="aff_w")
                #st.text("Individual Binding Affinity Weights")
                st.number_input("MHC Flurry BA Weight", value=state.aff_w, key="aff_flurry")
                st.number_input("NetMHC BA Weight", value=state.aff_w, key="aff_netmhc")
            with col2: 
                st.number_input("Overall Processing Score Weight", value = state.pro_w, key="pro_w")
                st.number_input("MHC Flurry PS Weight", value=state.aff_w, key="pro_flurry")
                st.number_input("NetMHC PS Weight", value=state.aff_w, key="pro_netmhc")
            with col3:                
                st.number_input("Overall Immunogenicity/Antigenicity Weight", value = state.imm_w, key="imm_w")
                st.number_input("Vaxijen Weight", value=state.aff_w, key="vaxijen_w")
                st.number_input("IEDB Immunogenicity Weight", value=state.aff_w, key="iedbimm_w")
                st.number_input("Prime Weight", value=state.aff_w, key="prime_w")
            with col4:
                st.number_input("Number of Binding Allele Weight", value = state.num_all_w, key="num_all_w")
            
            update_rank_dict = {
                "data_dir" : data_dir,
                "taxonomy_id": taxonomy_id,
                "include_assay_results": True                
                }
            st.form_submit_button('Calculate Ranking',
                                  on_click = update_rank,
                                  kwargs = update_rank_dict)
            
        st.dataframe(state.df)
 
def update_rank(data_dir, taxonomy_id, include_assay_results):
    state.df = get_df(data_dir, taxonomy_id, include_assay_results)

def get_df(data_dir, taxonomy_id, include_assay_results):
    master_list_df = create_master_list(taxonomy_id, data_dir)
    ranked_ml_df = get_rank(master_list_df,
                    aff_w=state.aff_w,
                    pro_w=state.pro_w,
                    imm_w=state.imm_w,
                    num_all_w=state.num_all_w,
                    aff_flurry=state.aff_flurry,
                    aff_netmhc=state.aff_netmhc,
                    aff_netctl=0.2,
                    pro_flurry=state.pro_flurry,
                    pro_netmhc=state.pro_netmhc,
                    vaxijen_w=state.vaxijen_w,
                    iedb_w=state.iedbimm_w,
                    prime_w=state.prime_w)
    
    if include_assay_results:
        #Following in case not starting from the beginning. SHould put a check
        bv_brc_df = pd.read_table(join(data_dir,f'bvbrc_assay_results_{taxonomy_id}.tsv'), usecols =['peptide','protein_id','bcell_assays', 'tcell_assays', 'mch_assays'])
        ranked_ml_df = pd.merge(ranked_ml_df,bv_brc_df, how="left", on=["peptide", "protein_id"])
    
    return ranked_ml_df

if __name__ == '__main__':
    main()

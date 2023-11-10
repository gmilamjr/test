import streamlit as st
from os.path import join
from PIL import Image
from styles import streamlit_style

#Get images
logo = Image.open(join("imgs","microorganisms.jpg"))
icon = Image.open(join("imgs","brachionus.jpg"))

streamlit_style('Microorganism Detection in Water', layout = 'wide', page_icon=icon)

st.markdown('''
<style>
[data-testid="stMarkdownContainer"] ul{
    padding-left:40px;
}
</style>
''', unsafe_allow_html=True)    

def main():
    
    st.title("Microorganism Detection in Water")
    col1, col2 = st.columns(2)
    
    with col1:
        st.image(logo)
        
    with col2:
        st.subheader("Project Description")
        st.write("In the aqueous tapestry that sustains life, the role of environmental microorganisms (EMs) in water \
                 takes center stage, especially when it comes to the paramount concern of public health. Every drop of \
                 water tells a story, and within it, microscopic actors play a crucial part. Detecting these environmental \
                 microorganisms is not just a scientific endeavor but a vital step in safeguarding communities against \
                 waterborne threats. From identifying potential pathogens that could jeopardize public health to ensuring \
                 the purity of our drinking water sources, the microscopic world within water holds the answers. By \
                 delving into the unseen realms, we equip ourselves with the knowledge needed to protect the most \
                 fundamental resource we all share—clean and safe water for the well-being of everyone.")
        
        st.subheader("Dataset")
        st.write("The overall dataset is made up of the following two datasets:") 
        st.markdown("- Environmental Microorganism Image Dataset Seventh Version is a microscopic image data set \
                    containing 41 types of EMs, with 2,65 images and 13,216 labeled objects in XML Format.")
    
        st.markdown("- Environmental Microorganism Image Dataset Sixth Version is a microscopic image data set containing \
                21 types of EMs. Each type of EM contains 40 original and 40 GT images, in total 1680 EM ages. \
                Box annotations were generated using the Roboflow tool, generating XML and YoloV8 files.")
        
    st.write("© 2023 Omdena. All rights reserved.")        
    
if __name__ == '__main__':
    
    main()
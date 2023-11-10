import streamlit as st

def streamlit_style(page_title, page_icon, layout):

    st.set_page_config(initial_sidebar_state = "auto", layout = layout, page_title = page_title, page_icon=page_icon)

    st.markdown('<style> div.block-container{padding-top:0rem;} </style>', unsafe_allow_html=True)
    
    st.markdown("""
    <style>
        button.step-up {display: none;}
        button.step-down {display: none;}
        div[data-baseweb] {border-radius: 4px;}
    </style>
    """,
    unsafe_allow_html=True)
    
    st.markdown(
        """
       <style>
       [data-testid="stSidebar"][aria-expanded="true"]{
           min-width: 200px;
           max-width: 200px;
       }
       """,
        unsafe_allow_html=True,
    )
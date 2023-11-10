###Import libs
import streamlit as st
from utils import get_binary_model, get_image
from styles import streamlit_style

streamlit_style('Water Contaminate Detection', layout = 'centered', page_icon=None)

st.title('Water Contaminate Detection')
model = get_binary_model("models/binary/yolo.pt")

# Gettting image
file = st.file_uploader("Upload image for prediction")

if file is not None:
    img = get_image(file)
    st.image(img, use_column_width=True, clamp = True)       
    
    results = model.predict(source=img)  # Predict contamination
    results_numpy = results[0].probs.numpy()
    probability_scores = results_numpy.data
           
    if probability_scores[0] > probability_scores[1]: 
       st.error(f"This image ***contains*** a harmful microorganisms. It predicted it with a probability of {probability_scores[0]:.4f}")
    elif probability_scores[1] > probability_scores[0]: 
       st.info(f"This image ***does not contain*** a harmful microorganism. It predicted it with a probability of {probability_scores[1]:.4f}")
    else:
       st.warning("Cannot be decided")





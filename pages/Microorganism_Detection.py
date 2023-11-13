###Import libs
import streamlit as st
from utils import get_detection_model, get_image, detection_img
from styles import streamlit_style

streamlit_style('Microorganism Detection', layout = 'wide', page_icon=None)

st.title('Microorganism Detection')
model = get_detection_model("models/object_detection/detectron2_weights.pth")

# Gettting image
file = st.file_uploader("Upload image for detections")

c1,c2=st.columns(2)
with c1:
    conf_threshold=float(st.slider("Confidence Threshold", min_value=0.0,max_value=1.0,value=0.3,step=0.02))
with c2:
    iou_threshold=float(st.slider("IOU Threshold", min_value=0.0,max_value=1.0,value=0.7,step=0.02))
    
if file is not None:
    
    img = get_image(file)

    col1,col2=st.columns(2)
    with col1:
        st.write("Uploaded Image")
        st.image(img)
    with col2:
        col3,col4=st.columns(2)
        with col3:
            st.write("Microorganims Detected")
        with col4:
            viz_name = st.toggle('Display Microorganism Name')
        with st.spinner("Predicting"):
            result_img = detection_img(model, img, conf_threshold, iou_threshold, viz_name)
            st.image(result_img)
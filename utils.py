import bysp
import cv2
import glob
import os
import streamlit as st
import numpy as np
import torch
import torchvision
import pandas as pd

from io import BytesIO
from skimage import io
from ultralytics import YOLO

from detectron2 import model_zoo
from detectron2.config import get_cfg
from detectron2.engine import DefaultPredictor

class_df = pd.read_csv("resources/class_mapping_alok.csv")

detectron2_path = "models/object_detection/detectron2_weights.pth"

if not os.path.exists(detectron2_path):
    bysp.combine_file(filename=detectron2_path, parts=glob.glob("models/object_detection/detectron2_weights.pth.*"), save=True)

##Binary classification model. Model Author::Alok
@st.cache_resource()
def get_binary_model(model_path):
    # load model
    model = YOLO(model_path)
    return model

### load Image
@st.cache_data(max_entries=1, show_spinner=False, ttl = 2*60)
def load_image(image):
     img = io.imread(image)
     return img
 
def get_image(file):
    file_bytes = BytesIO(file.read())    
    img = load_image(file_bytes)

    if (len(img) < 3):
       img = np.stack((img,) * 3, axis=-1)
       
    img = img[:,:,:3]

    return img

@st.cache_resource()
def get_detection_model(model_path=detectron2_path):
    # load model
    cfg = get_cfg()
    cfg.merge_from_file(model_zoo.get_config_file('COCO-Detection/retinanet_R_101_FPN_3x.yaml'))
    cfg.MODEL.WEIGHTS = model_path
    cfg.MODEL.DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    predictor = DefaultPredictor(cfg)
    
    return predictor

@st.cache_data(max_entries=1, show_spinner=False, ttl = 2*60)
def get_img_results(_model, img, conf_threshold, iou_threshold):
    results = _model(img)
    processed_results = preprocess_bbox(results['instances'], conf_threshold, iou_threshold)
    return processed_results

def detection_img(model, img, conf_threshold, iou_threshold, viz_name=False):
    
    processed_results = get_img_results(model, img, conf_threshold, iou_threshold)
    img = show_bbox(img, processed_results, class_df['Class Name' if viz_name else 'Class ID'])
    
    return img
    
#Modified from https://github.com/sudhanshu2198/Environmental-Microorganism-Detection/blob/master/utils.py
def preprocess_bbox(predictions,conf_threshold,iou_threshold):
    indeces = predictions.scores>=conf_threshold
    
    processed_bbox={}
    boxes=predictions.pred_boxes.tensor[indeces]
    scores=predictions.scores[indeces]
    labels=predictions.pred_classes[indeces]

    nms=torchvision.ops.nms(boxes,scores,iou_threshold=iou_threshold)

    processed_bbox["boxes"]=boxes[nms]
    processed_bbox["scores"]=scores[nms]
    processed_bbox["labels"]=labels[nms]
    
    return processed_bbox

#Modified from https://github.com/sudhanshu2198/Environmental-Microorganism-Detection/blob/master/utils.py
#Add color mapping, different colors for different classes
def show_bbox(img,target,classes,color=(0,0,255)):
    boxes=target["boxes"].numpy().astype("int")
    labels=target["labels"].numpy()
    scores=target["scores"].numpy()
    img=img.copy()
    
    for i,box in enumerate(boxes):
        text=f"{classes[labels[i]]}-{scores[i]:.2f}"
        cv2.rectangle(img,(box[0],box[1]),(box[2],box[3]),color,4)
        y=box[1]-10 if box[1]-40>40 else box[1]+40
        cv2.putText(img,text,(box[0],y),cv2.FONT_HERSHEY_SIMPLEX,1,color,2)
    return img
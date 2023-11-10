import cv2
import streamlit as st
import numpy as np
import torch
import torchvision

from io import BytesIO
from skimage import io
from ultralytics import YOLO

from detectron2.config import get_cfg
from detectron2.engine import DefaultPredictor
from detectron2 import model_zoo

classes= [ "G" + str(i).zfill(3) for i in range(1,63)]  # Best to take from the class mapping file

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
def get_detection_model(model_path):
    # load model
    cfg = get_cfg()
    cfg.merge_from_file(model_zoo.get_config_file('COCO-Detection/retinanet_R_101_FPN_3x.yaml'))
    cfg.MODEL.WEIGHTS = model_path
    cfg.MODEL.DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    predictor = DefaultPredictor(cfg)
    
    return predictor


def detection_img(model, img, conf_threshold, iou_threshold):
    
    results = model(img)
    processed_results = preprocess_bbox(results['instances'], conf_threshold, iou_threshold)
    img = show_bbox(img, processed_results, classes)
    
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
    
    
    '''
    look at unique labels -> objects
    N = len(objects)
    colors = sns.color_palette(palette = 'plasma', n_colors=N)
    c_dict = dict(map(lambda i,j : (i,tuple([int(c * 255) for c in j])) , objects, colors))'''
    
    for i,box in enumerate(boxes):
        text=f"{classes[labels[i]]}-{scores[i]:.2f}"
        cv2.rectangle(img,(box[0],box[1]),(box[2],box[3]),color,4)
        y=box[1]-10 if box[1]-40>40 else box[1]+40
        cv2.putText(img,text,(box[0],y),cv2.FONT_HERSHEY_SIMPLEX,1,color,2)
    return img
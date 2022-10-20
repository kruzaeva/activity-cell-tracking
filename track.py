import numpy as np


from tifffile import imread


from scipy import optimize

from utils import activity_map, connect, generatedic
import glob
from tifffile import imsave

folder="data/12/"

imgslist=sorted(glob.glob(folder+"im/*.tif"))
maskslist=sorted(glob.glob(folder+"seg/*.tif"))

n=len(imgslist)
masks=[]
imgs=[]

print("loading data")
for i in range(0,n):
    masks.append(imread(maskslist[i]).astype('float32'))
    imgs.append(imread(imgslist[i]).astype('float32'))
print(np.shape(masks))
print("data is loaded")
print("calculate activity map")
masks=np.array(masks)
imgs=np.array(imgs)
masksnew=activity_map(masks,imgs)
print("activity map is calculated")

celldict=generatedic(masks,masksnew)

for key in celldict[0].keys():
    celldict[0][key]["index"]=key

for frame in list(celldict.keys())[:-1]:

    celldict=connect(celldict,imgs, frame,frame+1)

txt=[]
for frame in list(celldict.keys()):
    for key in list(celldict[frame].keys()):
        if ("dod" in celldict[frame][key]) or frame==list(celldict.keys())[-1]:
            txt.append([celldict[frame][key]["index"], celldict[frame][key]["dob"],frame,celldict[frame][key]["initial_mother"]])

txt=np.array(txt)
with open(folder+"res/res_track.txt" , 'wb') as f:
    np.savetxt(f, txt.astype(int), newline='\n', header='', footer='', comments='# ',fmt="%i")

masks_gt=1*masks
imsave(folder+"res/mask"+str(frame).zfill(3)+".tif", masks_gt[0].astype(np.uint16))
for frame in list(celldict.keys())[:-1]:
    print(frame)
    masks2=1*masks[frame+1]
    k=frame
    j=frame+1
    masks_connected=1*masks2*1000
    for key2 in list(celldict[j].keys()):
        if "index" in celldict[j][key2]:
            masks_gt[frame+1][masks_connected==key2*1000]=int(celldict[j][key2]["index"])
            imsave(folder+"res/mask"+str(frame).zfill(3)+".tif", masks_gt[frame+1].astype(np.uint16))
import numpy as np


import scipy
import operator

from numba import njit


def gaussian_heatmap(center = (2, 2), image_size = (10, 10), sig = 1):
    """
    It produces single gaussian at expected center
    :param center:  the mean position (X, Y) - where high value expected
    :param image_size: The total image size (width, height)
    :param sig: The sigma value
    :return:
    """
    x_axis = np.linspace(0, image_size[0]-1, image_size[0]) - center[0]
    y_axis = np.linspace(0, image_size[1]-1, image_size[1]) - center[1]
    xx, yy = np.meshgrid(x_axis, y_axis)
    kernel = np.exp(-0.5 * (np.square(xx) + np.square(yy)) / np.square(sig))
    return kernel.T

@njit
def activity_map(masks,imgs):
    masksnew=np.zeros(masks.shape)
    for i in range(0,masks.shape[0]-1):
        print(i)
        activemap=0.5*np.abs(imgs[i]-imgs[i+1])
        for j in range(1,int(masks[i].max()+1)):
            ma=(masks[i]==j)
            area=np.sum(ma)
            active_inten=np.sum(activemap*ma)/(area+0.0001)
            masksnew[i]=masksnew[i]+active_inten*ma

    return masksnew
def find_maxind(dic):
    maxind=[]
    for c in dic.keys():
        maxind.append(dic[c]["index"])
    return np.max(maxind)

def generatedic(masks,masksnew):
    all_dic={}
    for m in range(0,np.shape(masks)[0]):
        celldict1={}
        #print(m)
        for i in range(1,int(masks[m].max()+1)):
            if np.sum(masks[m]==i)>0:
                celllowerlevel={}
                celllowerlevel["cm"]=np.mean(np.where(masks[m]==i),axis=1)
                celllowerlevel["activity"]=masksnew[m,celllowerlevel["cm"][0].astype(int),celllowerlevel["cm"][1].astype(int)]
                celllowerlevel["area"]=np.sum(masks[m]==i)
                celllowerlevel["mother"]=None
                celllowerlevel["index"]=None
                celllowerlevel["dob"]=0
                celllowerlevel["initial_mother"]=0
                if celllowerlevel["area"]>0:
                    celldict1[i]=celllowerlevel
        all_dic[m]=celldict1
        print(m)
    return all_dic

def connect(celldict, imgs, frame1,frame2, k=2.5, areaind=0.9):
    adjacdic={}
    newIndex = sorted(celldict[frame1], key=lambda x: celldict[frame1][x]['activity'])
    leftover1=list(newIndex)
    leftover2=list(celldict[frame2].keys())
    maxmax=[find_maxind(celldict[frame1])+1]
    nonekeys=[]
    for key1 in newIndex:
        subdic={}
        gm=gaussian_heatmap(center = (int(celldict[frame1][key1]["cm"][0]), int(celldict[frame1][key1]["cm"][1])), image_size = imgs[0].shape, sig = celldict[frame1][key1]["activity"]/k)
        for key2 in leftover2:
            pr=gm[int(celldict[frame2][key2]['cm'][0]),int(celldict[frame2][key2]['cm'][1])]
            if pr>0.01:
                
                subdic[key2]=pr

        if len(subdic) == 0:
            leftover1.remove(int(key1))
        else:
            key2=max(subdic, key=subdic.get)
            if areaind*celldict[frame1][key1]['area']<celldict[frame2][key2]['area']:

                if  celldict[frame1][key1]["index"] is not None:
                    celldict[frame2][key2]['mother']=1*celldict[frame1][key1]['index']
                    celldict[frame2][key2]["index"]=celldict[frame2][key2]['mother']
                    celldict[frame2][key2]["dob"]=celldict[frame1][key1]["dob"]
                    celldict[frame2][key2]["initial_mother"]=celldict[frame1][key1]['initial_mother']
                    maxmax.append(celldict[frame1][key1]['index'])

                else:
                    nonekeys.append(key2)
                leftover2.remove(max(subdic, key=subdic.get))
                leftover1.remove(int(key1))
            subdic=sorted(subdic.items(),key=operator.itemgetter(1),reverse=True)
            adjacdic[key1]=dict(subdic)

    for key1 in adjacdic.copy().keys():
        if key1 not in leftover1: 
            del adjacdic[key1]
        else:
            for key2 in leftover2:
                if key2 not in leftover2: 
                    del adjacdic[key1][key2]

    d=leftover2
    m=leftover1
    dm=m+m
    forhung2=np.zeros([len(dm),len(d)])
    i=0
    for mi in dm:
        j=0
        for dj in d:
            if dj in set(adjacdic[mi].keys()):
                forhung2[i,j]=1-adjacdic[mi][dj]
            else:
                forhung2[i,j]=1
            j=j+1
        i=i+1
    pairsm,pairsd=scipy.optimize.linear_sum_assignment(forhung2)
    maxmax=np.max(maxmax)+1
    #print(leftover2)
    foundd=[]
    for i in range(len(pairsm)):

        celldict[frame2][d[pairsd[i]]]['mother']=celldict[frame1][dm[pairsm[i]]]["index"]
        celldict[frame2][d[pairsd[i]]]['dob']=frame2
        celldict[frame2][d[pairsd[i]]]["index"]=maxmax
        celldict[frame2][d[pairsd[i]]]['initial_mother']=celldict[frame1][dm[pairsm[i]]]["index"]
        celldict[frame1][dm[pairsm[i]]]["dod"]=True
        foundd.append(d[pairsd[i]])
        maxmax=maxmax+1
    dif=set(leftover2).difference(set(foundd))
    #print(dif)
    nonekeys=nonekeys+list(dif)
    #print(nonekeys)
    for k in nonekeys:
        celldict[frame2][k]["index"]=maxmax
        maxmax=maxmax+1
    return celldict

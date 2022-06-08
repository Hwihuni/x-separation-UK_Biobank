from os.path import splitext
from os import listdir
import numpy as np
import scipy.io
import torch
from torch.utils.data import Dataset
import logging
from PIL import Image
import hdf5storage
from torchvision import transforms
import random

class BasicDataset(Dataset):
    def __init__(self,path,istrain,transform = None):
        load = hdf5storage.loadmat(path)
        self.transform = transform
        self.t1= load['t1_norm'].astype('float32')
        self.flair= load['flair_norm'].astype('float32')
        self.r2 = load['r2_norm'].astype('float32')
        
        logging.info(f'Creating dataset with {self.r2.shape[2]} examples')\
        
    @classmethod
    def preprocess(cls, x,y):
        img_trans_x = np.transpose(x,(2,0,1))
        img_trans_y = np.transpose(y,(2,0,1))
        return img_trans_x,img_trans_y
    
    def __len__(self):
        return self.r2.shape[2]
    
    def __getitem__(self, i):


        if self.transform:
            seed = np.random.randint(2147483647)
            random.seed(seed)
            torch.manual_seed(seed)
            t1 = self.transform(self.t1[:,:,i])
            random.seed(seed)
            torch.manual_seed(seed)
            flair = self.transform(self.flair[:,:,i])
            random.seed(seed)
            torch.manual_seed(seed)
            r2 = self.transform(self.r2[:,:,i])
            img = torch.cat((t1,flair),0)
            target = r2

            return {
                'image': img,
                'target': target
            }
        else:
            t1 = self.t1[:,:,i]
            flair = self.flair[:,:,i]
            r2 = self.r2[:,:,i]
            img = np.concatenate((np.expand_dims(t1,0),np.expand_dims(flair,0)),0)
            target = np.expand_dims(r2,0)

            assert img.size == 2*target.size, \
                f'Image and mask {i} should be the same size, but are {img.size} and {target.size}'
            return {
                'image': torch.from_numpy(img).type(torch.FloatTensor),
                'target': torch.from_numpy(target).type(torch.FloatTensor)
            }
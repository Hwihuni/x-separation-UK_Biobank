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
        self.img= load['img'].astype('float32')
        self.mask= load['Mask'].astype('float32')
        self.r2star = load['r2star_norm'].astype('float32')
        
        logging.info(f'Creating dataset with {self.r2star.shape[2]} examples')\
       
    #@classmethod  
       
    def __len__(self):
        return self.r2star.shape[2]
    
    def __getitem__(self, i):


        if self.transform:
            seed = np.random.randint(2147483647)
            random.seed(seed)
            torch.manual_seed(seed)
            #img0 = torch.unsqueeze(self.transform(np.squeeze(self.img[:,:,i,0])), 1)
            #img1 = torch.unsqueeze(self.transform(np.squeeze(self.img[:,:,i,0])), 1)
            img0 = self.transform(np.squeeze(self.img[:,:,i,0]))
            img1 = self.transform(np.squeeze(self.img[:,:,i,1]))
            img = torch.cat((img0,img1), dim=0)

            random.seed(seed)
            torch.manual_seed(seed)
            r2star = self.transform(self.r2star[:,:,i])
            target = r2star
            return {
                'image': img,
                'target': target
            }
        else:
            img = np.transpose(np.squeeze(self.img[:,:,i,:]),(2,0,1))
            r2star = self.r2star[:,:,i]
            target = np.expand_dims(r2star,0)

            assert img.size == 2*target.size, \
                f'Image and mask {i} should be the same size, but are {img.size} and {target.size}'
            return {
                'image': torch.from_numpy(img).type(torch.FloatTensor),
                'target': torch.from_numpy(target).type(torch.FloatTensor)
            }
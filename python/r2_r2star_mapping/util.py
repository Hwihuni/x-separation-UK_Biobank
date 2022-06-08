import torch
import torch.nn as nn
import numpy as np

def init_weights(m):
    if type(m) == nn.Conv2d or type(m) ==nn.ConvTranspose2d or type(m) ==nn.Linear:
        nn.init.xavier_uniform_(m.weight)

def psnr(img1, img2, max_value=255):
    mse = np.mean((np.array(img1, dtype=np.float32) - np.array(img2, dtype=np.float32)) ** 2)
    if mse == 0:
        return 100
    return 20 * np.log10(max_value / (np.sqrt(mse)))

    
def rmse(predictions, targets,max_value):
    return np.sqrt(((predictions - targets) ** 2).mean())/max_value
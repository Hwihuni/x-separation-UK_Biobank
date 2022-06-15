import argparse
import logging
import os
import scipy.io
import glob
import numpy as np
import torch

from unet import *
from utils.dataset import BasicDataset
from torch.utils.data import DataLoader
dir_checkpoint = 'checkpoints/'


inf_dir = './inf/'

def predict_net(net_u,net_fc,device):

    net_u.eval()
    net_fc.eval()
    val = BasicDataset(path_val,istrain = False)
    n_val = len(val)
    batch_size=6
    val_loader = DataLoader(val, batch_size=batch_size, shuffle=False, num_workers=8, pin_memory=True, drop_last=True)
    out = np.zeros((n_val,1,198,256))
    i = 0
    for batch in val_loader:
            imgs, true_masks = batch['image'], batch['target']
            
            imgs = imgs.to(device=device, dtype=torch.float32)
            true_masks = true_masks.to(device=device, dtype=torch.float32)

            with torch.no_grad():
                b1_pre = net_u(imgs)
                x_mid = torch.cat([b1_pre, imgs], dim=1)
                mask_pred = net_fc(x_mid)

            im_pred = mask_pred.cpu().detach().numpy()
            out[i*batch_size:(i+1)*batch_size,:,:,:] = im_pred
            i = i+1
    logging.info('Inference done')
    return out


def get_args():
    parser = argparse.ArgumentParser(description='Train the UNet on images and target masks',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-b', '--batch-size', metavar='B', type=int, nargs='?', default=30,
                        help='Batch size', dest='batchsize')
    parser.add_argument('-l', '--learning-rate', metavar='LR', type=float, nargs='?', default=1e-4,
                        help='Learning rate', dest='lr')
    parser.add_argument('-f', '--load', dest='load', type=str, default=True,
                        help='Load model from a .pth file')
    parser.add_argument('-s', '--scale', dest='scale', type=float, default=1,
                        help='Downscaling factor of the images')
    parser.add_argument('-v', '--validation', dest='val', type=float, default=10.0,
                        help='Percent of the data that is used as validation (0-100)')

    return parser.parse_args()


if __name__ == '__main__':
    os.environ["CUDA_VISIBLE_DEVICES"] = '5'
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    args = get_args()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logging.info(f'Using device {device}')

    net_u = UNet(n_channels=2, n_classes=1, bilinear=True)
    net_fc = Fc(n_channels=3, n_classes=1)
    
    net_u.load_state_dict(torch.load('./checkpoints/DS_Unet_prc_normal_10_23_21_11.pth', map_location=device))
    net_fc.load_state_dict(torch.load('./checkpoints/DS_Fc_prc_normal_10_23_21_11.pth', map_location=device))

    net_u.to(device=device)
    net_fc.to(device=device)


    try:
        path = "./data/Data_to_gpu_t2map*.mat"
        file_list = glob.glob(path)
        print ("file_list_py: {}".format(file_list))
        for step,path_val in enumerate(file_list):
            print(inf_dir+f'inference_r2_{step}.mat')
            img = predict_net(net_u=net_u,net_fc = net_fc,device=device)
            savedict = {'T2':np.squeeze(img[:,0,:,:])}
            load = scipy.io.savemat(inf_dir+f'inference_r2_{step}.mat',savedict)
            logging.info('Inference saved')

    except KeyboardInterrupt:
        logging.info('Saved interrupt')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
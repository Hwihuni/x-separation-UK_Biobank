""" Full assembly of the parts to form the complete network """

import torch.nn.functional as F
import math
from .unet_parts import *


class UNet(nn.Module):
    def __init__(self, n_channels, n_classes, bilinear=True):
        super(UNet, self).__init__()
        self.n_channels = n_channels
        self.n_classes = n_classes

        self.bilinear = bilinear

        self.inc = DoubleConv(self.n_channels, 64)
        self.down1 = Down(64, 128)
        self.down2 = Down(128, 256)
        self.down3 = Down(256, 512)
        factor = 2 if bilinear else 1
        self.down4 = Down(512, 1024 // factor, 1024)
        self.up1 = Up(1024, 512 // factor, bilinear)
        self.up2 = Up(512, 256 // factor, bilinear)
        self.up3 = Up(256, 128 // factor, bilinear)
        self.up4 = Up(128, 64, bilinear)
        self.outc = OutConv(64, self.n_classes)


    def forward(self, x):

        x1 = self.inc(x)
        x2 = self.down1(x1)
        x3 = self.down2(x2)
        x4 = self.down3(x3)
        x5 = self.down4(x4)
        x6 = self.up1(x5, x4)
        x7 = self.up2(x6, x3)
        x8 = self.up3(x7, x2)
        #
        x9 = self.up4(x8, x1)
        logits = self.outc(x9)

        return logits
    
class Fc(nn.Module):
    def __init__(self,n_channels, n_classes):
        super(Fc, self).__init__()
        self.n_channels2 = n_channels
        self.n_classes2 = n_classes
        self.fc1 = Conv_1D(self.n_channels2, 160)
        self.fc2 = Conv_1D(160, 240)
        self.fc3 = Conv_1D(240+self.n_channels2, 320)
        self.fc4 = Conv_1D(320, 360)
        self.fc5 = Conv_1D(360+self.n_channels2, 480)
        self.fc6 = Conv_1D(480, 520)
        self.fc7 = Conv_1D(520+self.n_channels2, 600)
        self.outfc = Conv_1D(600, self.n_classes2)
    def forward(self, x_mid):
        x21 = self.fc1(x_mid)
        x22 =  torch.cat([x_mid,self.fc2(x21)], dim=1)
        x23 = self.fc3(x22)
        x24 = torch.cat([x_mid,self.fc4(x23)], dim=1)
        x25 = self.fc5(x24)
        x26 = torch.cat([x_mid,self.fc6(x25)], dim=1)
        x27 = self.fc7(x26)
        logits = self.outfc(x27)
        return logits

class Model_int(nn.Module):
    def __init__(self, n_channels, n_classes, bilinear=True):
        super(Model_int, self).__init__()
        self.n_channels1 = 1
        self.n_classes1 = 1

        self.n_channels2 = 3
        self.n_classes2 = 1

        self.bilinear = bilinear

        self.inc = DoubleConv(self.n_channels1, 64)
        self.down1 = Down(64, 128)
        self.down2 = Down(128, 256)
        self.down3 = Down(256, 512)
        factor = 2 if bilinear else 1
        self.down4 = Down(512, 1024 // factor, 1024)
        self.up1 = Up(1024, 512 // factor, bilinear)
        self.up2 = Up(512, 256 // factor, bilinear)
        self.up3 = Up(256, 128 // factor, bilinear)
        self.up4 = Up(128, 64, bilinear)
        self.outc = OutConv(64, self.n_classes1)
        self.pool = nn.AvgPool2d(3)
        self.intep = nn.Upsample(scale_factor=3)
        self.fc1 = Conv_1D(self.n_channels2, 160)
        self.fc2 = Conv_1D(160, 240)
        self.fc3 = Conv_1D(240+self.n_channels2, 320)
        self.fc4 = Conv_1D(320, 360)
        self.fc5 = Conv_1D(360+self.n_channels2, 480)
        self.fc6 = Conv_1D(480, 520)
        self.fc7 = Conv_1D(520+self.n_channels2, 600)
        self.outfc = Conv_1D(600, self.n_classes2)

    def forward(self, xin):

        x = xin[:,1:2,:,:]/(1e-12+xin[:,0:1,:,:])
        x1 = self.inc(x)
        x2 = self.down1(x1)
        x3 = self.down2(x2)
        x4 = self.down3(x3)
        x5 = self.down4(x4)
        x6 = self.up1(x5, x4)
        x7 = self.up2(x6, x3)
        x8 = self.up3(x7, x2)
        #
        x9 = self.up4(x8, x1)
        #mid = self.intep(self.pool(self.outc(x9)))
        mid = self.outc(x9)
        x_mid = torch.cat([mid, xin], dim=1)
        x21 = self.fc1(x_mid)
        x22 =  torch.cat([x_mid,self.fc2(x21)], dim=1)
        x23 = self.fc3(x22)
        x24 = torch.cat([x_mid,self.fc4(x23)], dim=1)
        x25 = self.fc5(x24)
        x26 = torch.cat([x_mid,self.fc6(x25)], dim=1)
        x27 = self.fc7(x26)
        logits = torch.cat([mid, self.outfc(x27)], dim=1)
        return logits

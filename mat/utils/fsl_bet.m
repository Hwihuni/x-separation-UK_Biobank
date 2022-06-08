function [ image_bet, mask_bet ] = fsl_bet( image,f,res)
%fsl bet2 function
%   filename : string filename WITHOUT nii, f: f value for bet
switch nargin
    case 1
        f = 0.5;
        res = [1,1,2];
    case 2
        res = [1,1,2];
end
[~, ~, ~] = mkdir('D:\forFSL');
save_nii(make_nii(abs(rot90(image,-1)),res),'D:\forFSL\temp.nii')
command = ['bash -c "source ~/.profile && bet2 /mnt/d/forFSL/temp /mnt/d/forFSL/temp_bet_f' num2str(f) ' -m -f ' num2str(f),'"'];
[status, cmdout] = system(command);
gunzip(['D:\forFSL\temp_bet_f' num2str(f) '.nii.gz']);
gunzip(['D:\forFSL\temp_bet_f' num2str(f) '_mask.nii.gz']);
imstr = load_nii(['D:\forFSL\temp_bet_f' num2str(f) '.nii']);
image_bet = rot90(imstr.img);
maskstr = load_nii(['D:\forFSL\temp_bet_f' num2str(f) '_mask.nii']);
mask_bet = double(rot90(maskstr.img));
delete('D:\forFSL\*.nii')
delete('D:\forFSL\*.nii.gz')
end

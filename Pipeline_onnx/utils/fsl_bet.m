function [ image_bet, mask_bet ] = fsl_bet( image,f,res)
%fsl bet2 function
%   filename : string filename WITHOUT nii, C: f value for bet
switch nargin
    case 1
        f = 0.5;
        res = [1,1,2];
    case 2
        res = [1,1,2];
end
[~, ~, ~] = mkdir('C:\forFSL');
save_nii(make_nii(abs(rot90(image,-1)),res),'C:\forFSL\temp.nii')
command = ['bash -c "source ~/.profile && bet2 /mnt/c/forFSL/temp /mnt/c/forFSL/temp_bet_f' num2str(f) ' -m -f ' num2str(f),'"'];
[status, cmdout] = system(command);
gunzip(['C:\forFSL\temp_bet_f' num2str(f) '.nii.gz']);
gunzip(['C:\forFSL\temp_bet_f' num2str(f) '_mask.nii.gz']);
imstr = load_nii(['C:\forFSL\temp_bet_f' num2str(f) '.nii']);
image_bet = rot90(imstr.img);
maskstr = load_nii(['C:\forFSL\temp_bet_f' num2str(f) '_mask.nii']);
mask_bet = double(rot90(maskstr.img));
delete('C:\forFSL\*.nii')
delete('C:\forFSL\*.nii.gz')
end

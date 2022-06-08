function [ out_image, R_mat] = fsl_flirt2( ref_image,in_inmage,DOF,res,init,apply)
%fsl fast function
%   filename : string filename WITHOUT nii, f: f value for fast
switch nargin
    case 2
        DOF = 12;
        res = [1,1,2];
        init = false;
        apply = false;
    case 3
        res = [1, 1, 2];
        init = false;
        apply = false;
    case 4
        init = false;
        apply = false;
    case 5
        init = false;
        apply = false;
        
end
[~, ~, ~] = mkdir('D:\forFSL');
save_nii(make_nii((rot90(in_inmage,-1)),res),'D:\forFSL\temp1.nii')
save_nii(make_nii((rot90(ref_image,-1)),res),'D:\forFSL\temp2.nii')
if init
    if apply
        command = ['bash -c "source ~/.profile && flirt -in /mnt/d/forFSL/temp1 -ref /mnt/d/forFSL/temp2 -dof ' num2str(DOF) ' -applyxfm -out /mnt/d/forFSL/temp3  -init /mnt/d/forFSL/R_mat -interp sinc"'];
    else
        command = ['bash -c "source ~/.profile && flirt -in /mnt/d/forFSL/temp1 -ref /mnt/d/forFSL/temp2 -dof ' num2str(DOF) ' -out /mnt/d/forFSL/temp3  -init /mnt/d/forFSL/R_mat -interp sinc"'];
    end
else
    command = ['bash -c "source ~/.profile && flirt -in /mnt/d/forFSL/temp1 -ref /mnt/d/forFSL/temp2 -out /mnt/d/forFSL/temp3 -omat /mnt/d/forFSL/R_mat -dof ' num2str(DOF) '"'];
end
system(command);
load('D:\forFSL\R_mat','-ascii','R_mat');

gunzip('D:\forFSL\*.nii.gz');
imstr = load_nii('D:\forFSL\temp3.nii');
out_image = rot90(imstr.img);

delete('D:\forFSL\*.nii')
delete('D:\forFSL\*.nii.gz')
end
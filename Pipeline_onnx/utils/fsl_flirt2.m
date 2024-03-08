function [ out_image, R_mat] = fsl_flirt2( ref_image,in_inmage,DOF,res,init,apply)
%fsl fast function
%   filename : string filename WITHOUT nii, C: f value for fast
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
[~, ~, ~] = mkdir('C:\forFSL');
save_nii(make_nii((rot90(in_inmage,-1)),res),'C:\forFSL\temp1.nii')
save_nii(make_nii((rot90(ref_image,-1)),res),'C:\forFSL\temp2.nii')
if init
    if apply
        command = ['bash -c "source ~/.profile && flirt -in /mnt/c/forFSL/temp1 -ref /mnt/c/forFSL/temp2 -dof ' num2str(DOF) ' -applyxfm -out /mnt/c/forFSL/temp3  -init /mnt/c/forFSL/R_mat -interp sinc"'];
    else
        command = ['bash -c "source ~/.profile && flirt -in /mnt/c/forFSL/temp1 -ref /mnt/c/forFSL/temp2 -dof ' num2str(DOF) ' -out /mnt/c/forFSL/temp3  -init /mnt/c/forFSL/R_mat -interp sinc"'];
    end
else
    command = ['bash -c "source ~/.profile && flirt -in /mnt/c/forFSL/temp1 -ref /mnt/c/forFSL/temp2 -out /mnt/c/forFSL/temp3 -omat /mnt/c/forFSL/R_mat -dof ' num2str(DOF) '"'];
end
system(command);
load('C:\forFSL\R_mat','-ascii','R_mat');

gunzip('C:\forFSL\*.nii.gz');
imstr = load_nii('C:\forFSL\temp3.nii');
out_image = rot90(imstr.img);

delete('C:\forFSL\*.nii')
delete('C:\forFSL\*.nii.gz')
end
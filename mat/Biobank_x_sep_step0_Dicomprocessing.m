path = 'D:\UKB\Real_data\1004671_3\';

di =[path 'QSM\' dir([path 'QSM\15*']).name '\'];
flist=dir([di '*.dcm']); 

echoes_num = 2;
coils_num = size(flist,1)/96;
slices_num = size(flist,1)/coils_num/echoes_num;
img_data = zeros([ size(dicomread([di flist(1).name])), slices_num*echoes_num,coils_num]);
coil = [];

for i = 1:size(flist,1)
    fname=flist(i).name;
    info1 = dicominfo([di fname]);
    if info1.InstanceNumber==1
    coil = cat(1,coil, string(info1.Private_0051_100f));
    end
    if size(coil,1) == coils_num
        break
    end
end
for i = 1:size(flist,1)
    fname=flist(i).name;
    info1 = dicominfo([di fname]);
    img_data(:,:,info1.InstanceNumber,coil == string(info1.Private_0051_100f))=dicomread(info1);
end
img = make_nii(rot90(img_data(:,:,1:slices_num,1),-1),[0.8 0.8 3]);
save_nii(img,'tmp.nii');
nii_info = niftiinfo('tmp.nii');
nii_info.Datatype = 'double';
nii_info.ImageSize = size(rot90(img_data(:,:,1:slices_num,1),-1));
nii_info.PixelDimensions = [str2num(info1.Private_0051_100c(9:11))/size(img_data,1) str2num(info1.Private_0051_100c(5:7))/size(img_data,2) 3];
mkdir([di(1:end-12) 'NII']);
mkdir([di(1:end-12) 'NII\MAG_TE1']);
mkdir([di(1:end-12) 'NII\MAG_TE2']);
for i = 1:coils_num
    niftiwrite(rot90(img_data(:,:,1:slices_num,i),-1), [di(1:end-12) 'NII\MAG_TE1\' num2str(i) '.nii'], nii_info, 'Compressed', true);
    niftiwrite(rot90(img_data(:,:,slices_num+1:end,i),-1), [di(1:end-12) 'NII\MAG_TE2\' num2str(i) '.nii'], nii_info, 'Compressed', true);
end
%% 

di =[path 'QSM\' dir([path 'QSM\18*']).name '\'];
flist=dir([di '*.dcm']); 

echoes_num = 2;
coils_num = size(flist,1)/96;
slices_num = size(flist,1)/coils_num/echoes_num;
img_data = zeros([ size(dicomread([di flist(1).name])), slices_num*echoes_num,coils_num]);
coil = [];

for i = 1:size(flist,1)
    fname=flist(i).name;
    info1 = dicominfo([di fname]);
    if info1.InstanceNumber==1
    coil = cat(1,coil, string(info1.Private_0051_100f));
    end
    if size(coil) == [coils_num, 1]
        break
    end
end
for i = 1:size(flist,1)
    fname=flist(i).name;
    info1 = dicominfo([di fname]);
    img_data(:,:,info1.InstanceNumber,coil == string(info1.Private_0051_100f))=dicomread(info1);
end
img = make_nii(rot90(img_data(:,:,1:slices_num,1),-1),[0.8 0.8 3]);
save_nii(img,'tmp.nii');
nii_info = niftiinfo('tmp.nii');
nii_info.Datatype = 'double';
nii_info.ImageSize = size(rot90(img_data(:,:,1:slices_num,1),-1));
nii_info.PixelDimensions = [str2num(info1.Private_0051_100c(9:11))/size(img_data,1) str2num(info1.Private_0051_100c(5:7))/size(img_data,2) 3];
mkdir([di(1:end-12) 'NII']);
mkdir([di(1:end-12) 'NII\PHA_TE1']);
mkdir([di(1:end-12) 'NII\PHA_TE2']);
for i = 1:coils_num
    niftiwrite(rot90(img_data(:,:,1:slices_num,i),-1), [di(1:end-12) 'NII\PHA_TE1\' num2str(i) '.nii'], nii_info, 'Compressed', true);
    niftiwrite(rot90(img_data(:,:,slices_num+1:end,i),-1), [di(1:end-12) 'NII\PHA_TE2\' num2str(i) '.nii'], nii_info, 'Compressed', true);
end
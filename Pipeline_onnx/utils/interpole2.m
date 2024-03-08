function Vq = interpole2(input,scale,option)
[X,Y,C,B]=ndgrid(0:1/(size(input,1)-1):1,0:1/(size(input,2)-1):1,0:1/(size(input,3)-1):1,0:1/(size(input,4)-1):1);
[Xq,Yq,Cq,Bq]=ndgrid(0:1/(scale*size(input,1)-1):1,0:1/(scale*size(input,2)-1):1,0:1/(size(input,3)-1):1,0:1/(size(input,4)-1):1);
Vq = interpn(X,Y,C,B,input,Xq,Yq,Cq,Bq,option);
end
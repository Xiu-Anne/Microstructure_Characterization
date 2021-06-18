%Pure pore + its interfaces => pore voxels satisfying total pore volume
Pore_ad=0;
Add_Pore=nnz(and(P_Pore,Fiber));
se=strel('sphere',2);P_Pore_Dilate=imdilate(P_Pore,se);P_Pore_Dilate=P_Pore_Dilate-and(P_Pore_Dilate,(Fiber));P_Pore_Dilate=P_Pore_Dilate>0;
%
Pore_Matrix=(OI==50);
M_Pore_Matrix=MeanGray_3(Pore_Matrix,OriginalImage);
N1_Inter=(M_Matrix-M_Pore_Matrix)/(M_Matrix-M_Pore)*nnz(Pore_Matrix);
%N1_Inter = 2*10^5;
Pore_Max=and(P_Pore_Dilate,Pore_Matrix);
in1=find(Pore_Max);
if round(N1_Inter)<=size(in1,1)
    random1=randperm(size(in1,1),round(N1_Inter));
else
    random1=randperm(size(in1,1),size(in1,1));
    Pore_ad=round(N1_Inter)-size(in1,1);
end
ind_p1=in1(random1);
%
Pore_Fiber = Fiber_Pore;
M_Pore_Fiber=MeanGray_3(Pore_Fiber,OriginalImage);
N2_Inter=(M_Fiber-M_Pore_Fiber)/(M_Fiber-M_Pore)*nnz(Pore_Fiber);
Pore_Fib=and(P_Pore_Dilate,Pore_Fiber);
in2=find(Pore_Fib);
if round(N2_Inter)<size(in2,1)
    random2=randperm(size(in2,1),round(N2_Inter));
else
    random2=randperm(size(in2,1),size(in2,1));
    Pore_ad=Pore_ad+round(N2_Inter)-size(in2,1); 
end
ind_p2=in2(random2);
%
Pore_add=zeros(size(OriginalImage));
Pore_add(ind_p1)=1;Pore_add(ind_p2)=1;
%
Pore_add_2=zeros(size(OriginalImage));
se=strel('sphere',1);Remaining_Pore=imdilate(Pore_add,se)-Pore_add;Remaining_Pore=Remaining_Pore-and(Remaining_Pore,P_Pore);Remaining_Pore=Remaining_Pore-and(Remaining_Pore,Fiber);
ind=[];random=[];ind_add=[];
ind=find(or(and(Remaining_Pore,Pore_Matrix),and(Remaining_Pore,Pore_Fiber)));
if (Add_Pore+Pore_ad)<size(ind,1)
    random=randperm(size(ind,1),(Add_Pore+Pore_ad));
else
    random=randperm(size(ind,1),size(ind,1));
    rest=(Add_Pore+Pore_ad)-size(ind,1)
end
ind_add=ind(random);
Pore_add_2(ind_add)=1;
%
Pore=[];
Pore=P_Pore-and(P_Pore,Fiber)+Pore_add+Pore_add_2;
Verify_Pore=OriginalImage-OriginalImage.*uint16(Pore);
%%
Matrix = OriginalImage>=0; Matrix = Matrix - Fiber_final; Matrix = Matrix - Pore;
Syn_Image = Pore*0 + Matrix*100 + Fiber_final*250;
Directory='C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\Synthetic_Images\Syn_Images\Syn_';
for i=1:size(Touch_Final,3)
    g{i}=uint8(Syn_Image(:,:,i));
    imwrite(g{i},sprintf('%s%05d.tif', Directory, i))
end 
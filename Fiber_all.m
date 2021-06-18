Directory = 'C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\Synthetic_Images\All_WS\WS_';
i1=1;i2=200;
for i=i1:i2
fname{i-i1+1}= sprintf('%s%05d.tif', Directory, i-1);
I{i-i1+1}=imread(fname{i-i1+1});
end
%
for i=1:(i2-i1)+1% Number of 2Dimages
    WS_Fib(:,:,i)=I{1,i};
end
%%
Gray_all= Gray_Fiber_Matrix + Gray_Fiber_Pore;%+Gray_Fiber_BP
P_BWPF=P_Pore;%+P_Matrix;BP+WP++Fiber
%
Sub_Fib=[];Rest_Fiber=[];
Fiber_add=zeros(size(OriginalImage));
%%
for i=1:max(max(max(WS_Fib)))
    WS1=(WS_Fib==i);
    se=strel('sphere',1);Eliminate=WS1-imerode(WS1,se);
    %
    Fib_WS=Fiber.*WS1;
    Fib_Fiber=nnz(Fib_WS);
    Fib_PFiber=nnz(P_Fiber.*WS1);
    %
    ind=[];Mass=[];
    ind=find(Gray_all.*double(WS1));
    Mass=Gray_all(ind);
    Fib_all=sum(Mass);
    
    Add_Fib=Fib_PFiber+Fib_all-Fib_Fiber;
    %
    if Add_Fib>0
        Sub_Fib=[Sub_Fib 0];
        j=1;ind_add=[];Fib_Before=zeros(size(OriginalImage));
        while and(j<5,Add_Fib>0)
            se=strel('sphere',j);Fib_Di=imdilate(Fib_WS,se);Fib_Dil=and(Fib_Di,WS1);%Fib_DilD=double(Fib_Dil);
            Fib_Dila=Fib_Dil-Fib_Before;Fib_Dilat=Fib_Dila-and(Fib_Dila,Eliminate);Fib_Dilate=Fib_Dilat-and(Fib_Dilat,P_BWPF);
            Fib_Before=Fib_Dil;
            ind=[];random=[];
            ind=find(Fib_Dilate);
            if round(Add_Fib)<=size(ind,1)
                random=randperm(size(ind,1),round(Add_Fib));
                Ad_B=round(Add_Fib);
            else
                random=randperm(size(ind,1),size(ind,1));
                Ad_B=size(ind,1);
            end
            Add_Fib=Add_Fib-Ad_B;
            ind_add=[ind_add;ind(random)];
            j=j+1;
        end
        Fiber_add(ind_add)=1; %Rest_Fiber=[Rest_Fiber Rest_F];
        clc; i 
    else Sub_Fib=[Sub_Fib -Add_Fib];
    end
end
%%
Eliminate=zeros(size(OriginalImage));
for i=1:max(max(max(WS_Fib)))
    WS1=(WS_Fib==i);
    se=strel('sphere',1);Eliminate=Eliminate+(WS1-imerode(WS1,se));
end
Eliminate = Eliminate >0;
%%
Fiber_final=(Fiber_add+Fiber)>0;
Fiber_final=Fiber_final-and(Fiber_final,Eliminate);
N_final=(sum(sum(sum(Gray_all+P_Fiber-Fiber_final))))%Take into account real fiber pure phases!!!!
%%
se=strel('sphere',1);
FF_dil=imdilate(Fiber_final,se)-Fiber_final;
Eliminate_all=(Eliminate+P_Pore+Fiber_final)>0;
FF_dilate=FF_dil-and(FF_dil,Eliminate_all);
%
Fiber_add_f=zeros(size(OriginalImage));
ind=[];random=[];ind_add=[];
ind=find(FF_dilate);
size(ind,1)
%%
random=randperm(1785441,1.5508e+06);%randperm(size(ind,1),N_final);
ind_add=ind(random);
Fiber_add_f(ind_add)=1;
%%
Fiber_final=(Fiber_add+Fiber+Fiber_add_f)>0;
(sum(sum(sum(Fiber_final))))/(500*500*200)*100
%%
Directory='C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\Synthetic_Images\Fiber_Final\F_';
for i=1:size(Fiber_final,3)
    g{i}=uint8(255*Fiber_final(:,:,i));
    imwrite(g{i},sprintf('%s%05d.tif', Directory, i))
end
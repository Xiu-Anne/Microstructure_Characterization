% NonTouch_All + Touch_Final=> WS sub-zones
Directory = 'C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\Synthetic_Images\Interfaces_WS\Inter_WS_';
i1=1;i2=size(Touch_Final,3);
for i=i1:i2
fname{i-i1+1}= sprintf('%s%05d.tif', Directory, i-1);
I{i-i1+1}=imread(fname{i-i1+1});
end
%
for i=1:(i2-i1)+1% Number of 2Dimages
    WS_Fib(:,:,i)=I{1,i};
end
%% Verify if each particle is separated into only 1 or multiple zones
All_Fibers=NonTouch_All + Touch_Final;
Fiberlabel=bwlabeln(All_Fibers,26);%label all particles
Siz_1=[];Siz_2=[];
for i=1:max(max(max(Fiberlabel)))
    Test=(Fiberlabel==i);
    Test_Verify=WS_Fib.*single(Test);
    Siz_1=[Siz_1 size(unique(nonzeros(Test_Verify)),1)];
end
%%
L_Inter=find(Siz_1>1);%particle separated into multiple zones
L_fib=1:max(max(max(Fiberlabel)));
L_fib_WS=nonzeros(L_fib-L_fib.*ismember(L_fib,L_Inter));
%% Check 1 examples: particle separated into multiple zones
Test=(Fiberlabel==L_Inter(2));
ind=[];x=[];y=[];z=[];
ind=find(Test);
[x, y, z] = ind2sub(size(OriginalImage),ind);
xmin=min(x);xmax=max(x);
ymin=min(y);ymax=max(y);
zmin=min(z);zmax=max(z);
%% Pure pore/matrix/fiber & their mean gray value
P_Pore=(OI==0);M_Pore=MeanGray_3(P_Pore,OriginalImage)
P_Matrix=(OI==150);M_Matrix=MeanGray_3(P_Matrix,OriginalImage)
P_Fiber=(OI==250);M_Fiber=MeanGray_3(P_Fiber,OriginalImage) %42966
%%
Fiber_Pore=(OI==100);%fiber/pore interfaces
Fiber_Matrix=(OI==200);%fiber/matrix interfaces
%
Gray_Fiber_Matrix=Fiber_Matrix - and(All_Fibers,Fiber_Matrix);
Gray_Fiber_Matrix=(double(Gray_Fiber_Matrix).*double(OriginalImage)-M_Matrix*double(Gray_Fiber_Matrix))./(double(Gray_Fiber_Matrix)*(M_Fiber-M_Matrix));
Gray_Fiber_Matrix(isnan(Gray_Fiber_Matrix))=0;
%
Gray_Fib_Pore=(double(Fiber_Pore).*double(OriginalImage)-M_Pore*double(Fiber_Pore))./(double(Fiber_Pore)*(M_Fiber-M_Pore));
Gray_Fib_Pore(isnan(Gray_Fib_Pore))=0;
Gray_Fiber_Pore=Gray_Fib_Pore-Gray_Fib_Pore.*double(and(All_Fibers,Fiber_Pore));
%
Gray_Fiber=Gray_Fiber_Matrix+Gray_Fiber_Pore+All_Fibers;%Matrix containing all fiber porportions: 1 for pure fiber voxels & fiber proportion in corresponding interfaces for the rest
%% particle separated into only 1 zone
Volume_Fib=[];
V1_fib=[]; V2_fib=[]; V3_fib=[];X_fib=[];Y_fib=[];Z_fib=[];
Vecteur1_fib=[];Vecteur3_fib=[];Vecteur2_fib=[];
%
for k =1:size(L_fib_WS,1)
    
    Img=[];ind=[];Mass=[];Ma=[];m=[];n=[];p=[];X=[];Y=[];Z=[];
    Img=(Fiberlabel==L_fib_WS(k));WS_c=WS_Fib.*single(Img);WS_corres=(WS_Fib==unique(nonzeros(WS_c)));
    se=strel('sphere',4);Img_Dil=imdilate(Img,se)-Img;Img_Dilate=and(Img_Dil,WS_corres);
    P_Im=and(Img_Dilate,Fiber_Pore);M_Im=and(Img_Dilate,Fiber_Matrix);
    D_Img=P_Im+M_Im+Img;
    Ma=Gray_Fiber.*D_Img;
    
    ind=find(D_Img);
    Mass=Ma(ind);
    [m, n, p] = ind2sub(size(D_Img),ind);

    x_c=(m'*Mass)/sum(Mass); y_c=n'*Mass/sum(Mass); z_c=p'*Mass/sum(Mass); 
    X_fib=[X_fib x_c]; Y_fib=[Y_fib y_c]; Z_fib=[Z_fib z_c]; Siz= size(m,1);
    Volume_Fib=[Volume_Fib sum(sum(sum(Mass)))];

    X=m'-ones(1,Siz).*x_c;Y=n'-ones(1,Siz).*y_c;Z=p'-ones(1,Siz).*z_c;

    XX=X.^2; XY=X.*Y; XZ=X.*Z;
    YY=Y.^2; YZ=Y.*Z; ZZ=Z.^2;

    M_2=[sum(Mass'.*(XX)) sum(Mass'.*XY) sum(Mass'.*XZ);sum(Mass'.*XY) sum(Mass'.*(YY)) sum(Mass'.*YZ);sum(Mass'.*XZ) sum(Mass'.*YZ) sum(Mass'.*(ZZ))];

    [V,Values_2] = eig(M_2);

    Dig=diag(Values_2);
    [D_a, idx]=sort(Dig,'ascend');
    Vec=V(:,idx);

    Values=sqrt(D_a);
  
    V1_fib=[V1_fib Values(1)]; V2_fib=[V2_fib Values(2)]; V3_fib=[V3_fib Values(3)]; 
    
    Vecteur1_fib=[Vecteur1_fib Vec(:,1)]; Vecteur2_fib=[Vecteur2_fib Vec(:,2)]; Vecteur3_fib=[Vecteur3_fib Vec(:,3)]; 
    
    clc; k
end
V11_fib=(V1_fib).^2;V22_fib=(V2_fib).^2;V33_fib=(V3_fib).^2;
%
a_fib=(12*V11_fib./(V3_fib./V1_fib)./(V2_fib./V1_fib)).^0.2;
b_fib=a_fib.*(V2_fib./V1_fib);
c_fib=a_fib.*(V3_fib./V1_fib);
ca_fib=V3_fib./V1_fib;
ba_fib=V2_fib./V1_fib;
% particle separated into multiple zones
Volume_Fib_2=[];
V1_fib_2=[]; V2_fib_2=[]; V3_fib_2=[];X_fib_2=[];Y_fib_2=[];Z_fib_2=[];
Vecteur1_fib_2=[];Vecteur3_fib_2=[];Vecteur2_fib_2=[];
for k =1:size(L_Inter,2)
    
    Img=[];ind=[];Mass=[];Ma=[];m=[];n=[];p=[];X=[];Y=[];Z=[];
    Img=(Fiberlabel==L_Inter(k));WS_c=WS_Fib.*single(Img);WS_corres=ismember(WS_Fib,unique(nonzeros(WS_c)));
    se=strel('sphere',4);Img_Dil=imdilate(Img,se)-Img;Img_Dilate=and(Img_Dil,WS_corres);
    P_Im=and(Img_Dilate,Fiber_Pore);M_Im=and(Img_Dilate,Fiber_Matrix);
    D_Img=P_Im+M_Im+Img;
    Ma=Gray_Fiber.*D_Img;
    
    ind=find(D_Img);
    Mass=Ma(ind);
    [m, n, p] = ind2sub(size(D_Img),ind);

    x_c=(m'*Mass)/sum(Mass); y_c=n'*Mass/sum(Mass); z_c=p'*Mass/sum(Mass); 
    X_fib_2=[X_fib_2 x_c]; Y_fib_2=[Y_fib_2 y_c]; Z_fib_2=[Z_fib_2 z_c]; Siz= size(m,1);
    Volume_Fib_2=[Volume_Fib_2 sum(sum(sum(Mass)))];

    X=m'-ones(1,Siz).*x_c;Y=n'-ones(1,Siz).*y_c;Z=p'-ones(1,Siz).*z_c;

    XX=X.^2; XY=X.*Y; XZ=X.*Z;
    YY=Y.^2; YZ=Y.*Z; ZZ=Z.^2;

    M_2=[sum(Mass'.*(XX)) sum(Mass'.*XY) sum(Mass'.*XZ);sum(Mass'.*XY) sum(Mass'.*(YY)) sum(Mass'.*YZ);sum(Mass'.*XZ) sum(Mass'.*YZ) sum(Mass'.*(ZZ))];

    [V,Values_2] = eig(M_2);

    Dig=diag(Values_2);
    [D_a, idx]=sort(Dig,'ascend');
    Vec=V(:,idx);

    Values=sqrt(D_a);
  
    V1_fib_2=[V1_fib_2 Values(1)]; V2_fib_2=[V2_fib_2 Values(2)]; V3_fib_2=[V3_fib_2 Values(3)]; 
    
    Vecteur1_fib_2=[Vecteur1_fib_2 Vec(:,1)]; Vecteur2_fib_2=[Vecteur2_fib_2 Vec(:,2)]; Vecteur3_fib_2=[Vecteur3_fib_2 Vec(:,3)]; 
    
end
V11_fib_2=(V1_fib_2).^2;V22_fib_2=(V2_fib_2).^2;V33_fib_2=(V3_fib_2).^2;
%
a_fib_2=(12*V11_fib_2./(V3_fib_2./V1_fib_2)./(V2_fib_2./V1_fib_2)).^0.2;
b_fib_2=a_fib_2.*(V2_fib_2./V1_fib_2);
c_fib_2=a_fib_2.*(V3_fib_2./V1_fib_2);
ca_fib_2=V3_fib_2./V1_fib_2;
ba_fib_2=V2_fib_2./V1_fib_2;
%% all fibers
X_fib_all=[X_fib X_fib_2];Y_fib_all=[Y_fib Y_fib_2];Z_fib_all=[Z_fib Z_fib_2];
a_fib_all=[a_fib a_fib_2];b_fib_all=[b_fib b_fib_2];c_fib_all=[c_fib c_fib_2];
Vecteur1_fib_all=[Vecteur1_fib Vecteur1_fib_2];Vecteur2_fib_all=[Vecteur2_fib Vecteur2_fib_2];Vecteur3_fib_all=[Vecteur3_fib Vecteur3_fib_2];
Volume_Fiber=[Volume_Fib Volume_Fib_2];
%% if exist nan in b_fib (or a_fib, c_fib)
X_fib_all = X_fib_all(find(~isnan(b_fib_all))); Y_fib_all = Y_fib_all(find(~isnan(b_fib_all))); Z_fib_all = Z_fib_all(find(~isnan(b_fib_all)));
a_fib_all = a_fib_all(find(~isnan(b_fib_all))); b_fib_all = b_fib_all(find(~isnan(b_fib_all))); c_fib_all = c_fib_all(find(~isnan(b_fib_all)));
Vecteur1_fib_all = Vecteur1_fib_all(:,find(~isnan(b_fib_all))); Vecteur2_fib_all = Vecteur2_fib_all(:,find(~isnan(b_fib_all))); Vecteur3_fib_all = Vecteur3_fib_all(:,find(~isnan(b_fib_all)));
Volume_Fiber = Volume_Fiber(find(~isnan(b_fib_all)));
%% plot aspect ratio
abc=c_fib_all;%./a_fib_all;
r0=0;
N=[];N_R=[];
for i = 1:199
    r=r0+i*1;
    R=and(abc>=r,abc<r+1);
    V_R=R.*Volume_Fiber;
    N=[N nnz(R)];
    N_R=[N_R sum(V_R)];
end
L_E=1.5:1:199.5;
figure, 
yyaxis left
plot(L_E, N)
yyaxis right
plot(L_E, N_R)
%% Orientation
B_Z=acos(Vecteur3_fib_all(3,:));
B_Y=acos(Vecteur3_fib_all(2,:));
B_XZ=atan(Vecteur3_fib_all(1,:)./Vecteur3_fib_all(3,:));
B_YZ=atan(Vecteur3_fib_all(2,:)./Vecteur3_fib_all(3,:));
%% plot fiber orientations
B_X=acos(Vecteur3_fib_all(1,:));
r0=-pi/45;
N=[];N_R=[];
for i = 1:92
    r=r0-pi/2+i*pi/90;
    R1=and(B_X>r,B_X<=r+pi/90);
    V_R1=R1.*Volume_Fiber;
    N=[N nnz(R1)];
    N_R=[N_R sum(sum(V_R1))];
end
L_E=-pi/2-pi/180:pi/90:pi/2+pi/180;
figure, 
yyaxis left
plot(L_E, N)
yyaxis right
plot(L_E, N_R)
ax=gca;
set(gca,'XTick',-pi/2:pi/4:pi/2)
ax.XTickLabel = {'-pi/2','-pi/4','0','pi/4','pi/2'};
%% save 
save('C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\Volume_Fiber','Volume_Fiber');
save('C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\X_cen_1','X_fib_all');
save('C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\Y_cen_1','Y_fib_all');
save('C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\Z_cen_1','Z_fib_all');
save('C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\R_a_1','a_fib_all');
save('C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\R_b_1','b_fib_all');
save('C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\R_c_1','c_fib_all');
save('C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\Vecteur1_1','Vecteur1_fib_all');
save('C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\Vecteur2_1','Vecteur2_fib_all');
save('C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\Vecteur3_1','Vecteur3_fib_all');
%% synthesize all parameters to create equivalent fibers
XX_min=[];XX_max=[];YY_min=[];YY_max=[];ZZ_min=[];ZZ_max=[];
for p=1:size(a_fib_all,2)
    X=[];Y=[];Z=[];E_OR=[];E_Cx=[];E_Cy=[];E_Cz=[];
    x1=-a_fib_all(p)/2+1;x2=a_fib_all(p)/2;
    y1=-b_fib_all(p)/2+1;y2=b_fib_all(p)/2;
    z1=-c_fib_all(p)/2+1;z2=c_fib_all(p)/2;
    X=[x1 x1 x1 x1 x2 x2 x2 x2];
    Y=[y1 y1 y2 y2 y1 y1 y2 y2];
    Z=[z1 z2 z1 z2 z1 z2 z1 z2];
    Rotate=[Vecteur1_fib_all(:,p) Vecteur2_fib_all(:,p) Vecteur3_fib_all(:,p)];
    E_OR=Rotate*[X' Y' Z']';
    E_Cx=E_OR(1,:)+ones(size(E_OR(1,:))).*X_fib_all(p);E_Cy=E_OR(2,:)+ones(size(E_OR(1,:))).*Y_fib_all(p);E_Cz=E_OR(3,:)+ones(size(E_OR(1,:))).*Z_fib_all(p);
    XX_min=[XX_min min(E_Cx)];YY_min=[YY_min min(E_Cy)];ZZ_min=[ZZ_min min(E_Cz)];
    XX_max=[XX_max max(E_Cx)];YY_max=[YY_max max(E_Cy)];ZZ_max=[ZZ_max max(E_Cz)];
end
%
XX_min(find(XX_min<=0.5))=1;YY_min(find(YY_min<=0.5))=1;ZZ_min(find(ZZ_min<=0.5))=1;
XX_max(find(XX_max>size(Touch_Final,1)))=size(Touch_Final,1);YY_max(find(YY_max>size(Touch_Final,2)))=size(Touch_Final,2);ZZ_max(find(ZZ_max>size(Touch_Final,3)))=size(Touch_Final,3);
%%
Fiber_R=zeros(size(P_Fiber));
for p=1:size(a_fib_all,2)
    clc; p
    for i=round(XX_min(p)):round(XX_max(p))
        for j=round(YY_min(p)):round(YY_max(p))
            for k=round(ZZ_min(p)):round(ZZ_max(p))
                L=[X_fib_all(1,p)-i Y_fib_all(1,p)-j Z_fib_all(1,p)-k];
                L1=dot(L,Vecteur1_fib_all(:,p));L2=dot(L,Vecteur2_fib_all(:,p));L3=dot(L,Vecteur3_fib_all(:,p));
                if and(and(abs(L1)<=a_fib_all(p)/2,abs(L2)<=b_fib_all(p)/2),abs(L3)<=c_fib_all(p)/2)                    
                     Fiber_R(i,j,k)=1;
                else Fiber_R(i,j,k)=0;
                end
            end
        end
    end
end
%% Save equivalent fibers
Fiber=Fiber_R>0;
Directory='C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\Synthetic_Images\Fiber_Inter_Final\';
for i=1:size(Touch_Final,3)
    g{i}=uint8(255*Fiber(:,:,i));
    imwrite(g{i},sprintf('%s%05d.tif', Directory, i))
end 
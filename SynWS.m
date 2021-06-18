% Center cordinates + aspect ratios + orientations of non-touching particles separated via watershed
XX_min=[];XX_max=[];YY_min=[];YY_max=[];ZZ_min=[];ZZ_max=[];
for p=1:size(a_WS,2)
    X=[];Y=[];Z=[];E_OR=[];E_Cx=[];E_Cy=[];E_Cz=[];
    x1=-a_WS(p)/2+1;x2=a_WS(p)/2;
    y1=-b_WS(p)/2+1;y2=b_WS(p)/2;
    z1=-c_WS(p)/2+1;z2=c_WS(p)/2;
    X=[x1 x1 x1 x1 x2 x2 x2 x2];
    Y=[y1 y1 y2 y2 y1 y1 y2 y2];
    Z=[z1 z2 z1 z2 z1 z2 z1 z2];
    Rotate=[Vecteur1_WS(:,p) Vecteur2_WS(:,p) Vecteur3_WS(:,p)];
    E_OR=Rotate*[X' Y' Z']';
    E_Cx=E_OR(1,:)+ones(size(E_OR(1,:))).*X_WS(p);E_Cy=E_OR(2,:)+ones(size(E_OR(1,:))).*Y_WS(p);E_Cz=E_OR(3,:)+ones(size(E_OR(1,:))).*Z_WS(p);    
    XX_min=[XX_min min(E_Cx)];YY_min=[YY_min min(E_Cy)];ZZ_min=[ZZ_min min(E_Cz)];
    XX_max=[XX_max max(E_Cx)];YY_max=[YY_max max(E_Cy)];ZZ_max=[ZZ_max max(E_Cz)];
end
%
XX_min(find(XX_min<=0.5))=1;YY_min(find(YY_min<=0.5))=1;ZZ_min(find(ZZ_min<=0.5))=1;
XX_max(find(XX_max>size(Touch_Final,1)))=size(Touch_Final,1);YY_max(find(YY_max>size(Touch_Final,2)))=size(Touch_Final,2);ZZ_max(find(ZZ_max>size(Touch_Final,3)))=size(Touch_Final,3);
%% Synthesize equivalent fibers
Fiber_RR=zeros(size(P_Fiber));%k=100;
for p=1:size(a_WS,2)
    for i=round(XX_min(p)):round(XX_max(p))
        for j=round(YY_min(p)):round(YY_max(p))
            for k=round(ZZ_min(p)):round(ZZ_max(p))
                L=[X_WS(1,p)-i Y_WS(1,p)-j Z_WS(1,p)-k];
                L1=dot(L,Vecteur1_WS(:,p));L2=dot(L,Vecteur2_WS(:,p));L3=dot(L,Vecteur3_WS(:,p));
                if and(and(abs(L1)<=a_WS(p)/2,abs(L2)<=b_WS(p)/2),abs(L3)<=c_WS(p)/2)                    
                     Fiber_RR(i,j,k)=1;
                else Fiber_RR(i,j,k)=0;
                end
            end
        end
    end
end
%% Saving synthetic fibers
FiberRR=Fiber_RR>0;
Directory='C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\Synthetic_Images\WS_Without_Interfaces\WS_WI_';
for i=1:size(Touch,3)
    g{i}=uint8(255*FiberRR(:,:,i));
    imwrite(g{i},sprintf('%s%05d.tif', Directory, i))
end 
%% Save synthetic WS fibers
NonTouch_WS = FiberRR;
NonTouch_All = NonTouch_PF + NonTouch_WS;
Directory='C:\Users\thi-xiu.le\Documents\Postdoc\Code Segmentation\Test3_1751_1950\NonTouch_All\NT_';
for i=1:size(Touch,3)
    g{i}=uint8(255*NonTouch_All(:,:,i));
    imwrite(g{i},sprintf('%s%05d.tif', Directory, i))
end 
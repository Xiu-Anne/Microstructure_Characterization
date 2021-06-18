%% Touching particles, separated via watershed with markers
Directory = 'WS_';
i1=1;i2=size(Touch,3);
for i=i1:i2
fname{i-i1+1}= sprintf('%s%05d.tif', Directory, i-1);
I{i-i1+1}=imread(fname{i-i1+1});
end
% 3D images
for i=1:(i2-i1)+1% Number of 2Dimages
    WS(:,:,i)=I{1,i};
end
%% Finding bouding box & volume of each particles
x=[];y=[];z=[];V=[];
xmin=[];ymin=[];zmin=[];xmax=[];ymax=[];zmax=[];
for i = 1:size(L_TA,1)
    ind=find(Image_label==L_TA(i));
    [xi, yi, zi] = ind2sub(size(Image_label),ind);
    [x]=[x; xi];y=[y; yi];z=[z; zi];
    xmin=[xmin min(xi)];ymin=[ymin min(yi)];zmin=[zmin min(zi)];
    xmax=[xmax max(xi)];ymax=[ymax max(yi)];zmax=[zmax max(zi)];
    V=[V size(ind,1)];
end
%% Verify whether each previous touching particle is separated into at least 2?
N_TA=[];
for i=1:size(L_TA,1)
    Veri=[];Verif=[];WS_T=[];WS_Test=[];
    %Test=Touch(xmin(1,i):xmax(1,i),ymin(1,i):ymax(1,i),zmin(1,i):zmax(1,i));
    Veri=Image_label(xmin(1,i):xmax(1,i),ymin(1,i):ymax(1,i),zmin(1,i):zmax(1,i));
    Verif=(Veri==L_TA(i));
    WS_T=WS(xmin(1,i):xmax(1,i),ymin(1,i):ymax(1,i),zmin(1,i):zmax(1,i));
    WS_Test=WS_T.*Verif;
    List=unique(nonzeros(WS_Test))';
    N_TA=[N_TA size(List,2)];
end
%% Verify if still exist touching particles or eliminate too small particles
NO=[];
L_TA_0=[];L_TA_01=[];Vol=[];V_in=[];%V_eq=[];
%X_WS=[];Y_WS=[];Z_WS=[];
%a_WS=[];b_WS=[];c_WS=[];
%Vecteur1_WS=[];Vecteur2_WS=[];Vecteur3_WS=[];
L_TA_1=L_TA(find(N_TA>1));
xmin_1=xmin(find(N_TA>1));ymin_1=ymin(find(N_TA>1));zmin_1=zmin(find(N_TA>1));
xmax_1=xmax(find(N_TA>1));ymax_1=ymax(find(N_TA>1));zmax_1=zmax(find(N_TA>1));
for i =1:size(L_TA_1,1)
    Veri=[];Verif=[];WS_Te=[];Tes=[];WS_Test=[];Test=[];
    Veri=Image_label(xmin_1(1,i):xmax_1(1,i),ymin_1(1,i):ymax_1(1,i),zmin_1(1,i):zmax_1(1,i));
    Verif=(Veri==L_TA_1(i));
    WS_Te=WS(xmin_1(1,i):xmax_1(1,i),ymin_1(1,i):ymax_1(1,i),zmin_1(1,i):zmax_1(1,i));
    Tes=Touch(xmin_1(1,i):xmax_1(1,i),ymin_1(1,i):ymax_1(1,i),zmin_1(1,i):zmax_1(1,i));
    WS_Test=WS_Te.*Verif;Test=Tes.*Verif;
    No=[];List=unique(nonzeros(WS_Test))';
    for j=1:size(List,2)
        No=[No nnz(WS_Test==List(1,j))];
    end
    No=No./sum(No);
    [N_a, idx]=sort(No,'descend');
    List_Na{i}=List(:,idx);
    %NO=[NO max(find(N_a>0.1))]
    if max(find(N_a>0.1))>1
        %
        Check=[];
        for l=1:max(find(N_a>0.1))
                WS_T0{l}=[];WS_TD{l}=[];
                WS_T0{l}=(WS_Test==List_Na{i}(l));
                se=strel('sphere',3);
                WS_TD{l}=imdilate(WS_T0{l},se);
                WS_TD{l}=and(WS_TD{l},Test);
        end
        for m=2:max(find(N_a>0.1))
                Check=WS_TD{m-1}+WS_TD{m};
        end
        %
        Ch=Check>1;
        for n =1:max(find(N_a>0.1))
            WS_TF =[];
            WS_TF{n}=WS_TD{n}-and(Ch,WS_TD{n});
        end
        %
        for k=1:max(find(N_a>0.1))
            ind=[];x=[];y=[];z=[];X=[];Y=[];Z=[];
            ind=find(WS_TF{k});
            [x, y, z] = ind2sub(size(WS_Test), ind);
            Siz= size(x,1);

            x_c=mean(x);y_c=mean(y);z_c=mean(z);
            %X_WS=[X_WS x_c+xmin_1(1,i)-1];Y_WS=[Y_WS y_c+ymin_1(1,i)-1];Z_WS=[Z_WS z_c+zmin_1(1,i)-1];

            Z=z'-ones(1,Siz).*z_c;X=x'-ones(1,Siz).*x_c;Y=y'-ones(1,Siz).*y_c;

            XX=sum(sum(X.^2)); XY=sum(sum(X.*Y)); XZ=sum(sum(X.*Z));
            YY=sum(sum(Y.^2)); YZ=sum(sum(Y.*Z)); ZZ=sum(sum(Z.^2));

            Moment_2=[XX XY XZ ; XY YY YZ ; XZ YZ ZZ];

            [V,Values_2] = eig(Moment_2);
            Dig=diag(Values_2);
            [D_a, idx]=sort(Dig,'ascend');
            Vec=V(:,idx);

            %
            Values=sqrt(D_a);

            %Vecteur1_WS=[Vecteur1_WS Vec(:,1)];Vecteur2_WS=[Vecteur2_WS Vec(:,2)];Vecteur3_WS=[Vecteur3_WS Vec(:,3)];

            a=(12*D_a(1,1)/(Values(3)/Values(1))/(Values(2)/Values(1)))^0.2;
            b=a*Values(2)/Values(1);c=a*Values(3)/Values(1);
            
            %c_WS=[c_WS c];b_WS=[b_WS b];c_WS=[c_WS c];

            Verify=[];
            for p = 1:size(X,2) 
                L=[X(1,p) Y(1,p) Z(1,p)];
                L1=dot(L,Vec(:,1));
                L2=dot(L,Vec(:,2));
                L3=dot(L,Vec(:,3));
                if or(or(abs(L1)>a/2,abs(L2)>b/2),abs(L3)>c/2)
                    Verify=[Verify 0];
                else Verify=[Verify 1];
                end
            end
            V_in(i,k)=nnz(Verify);
            Vol(i,k)=Siz;
        end
    else L_TA_0=[L_TA_0 L_TA_1(i)];
    end    
end
Ver=(Vol-V_in)./Vol;
[ik kk] = ind2sub(size(Ver), find(Ver>0.32));%threhold to distinguish between touching & non-touching particles
L_TA_01=[L_TA(find(N_TA==1)); L_TA_0'; L_TA_1(ik)];
L_TA_01=unique(L_TA_01);
%% Find aspect ratios & orientation of only non-touching particles
L_NT_01=nonzeros(L_TA-L_TA.*ismember(L_TA,L_TA_01));
id=1:size(L_TA,1);idx=id'.*ismember(L_TA,L_NT_01);idxx=nonzeros(idx);
xxmin=xmin(idxx);yymin=ymin(idxx);zzmin=zmin(idxx);
xxmax=xmax(idxx);yymax=ymax(idxx);zzmax=zmax(idxx);
X_WS=[];Y_WS=[];Z_WS=[];
a_WS=[];b_WS=[];c_WS=[];
Vecteur1_WS=[];Vecteur2_WS=[];Vecteur3_WS=[];
for i =1:size(L_NT_01,1)
    WS_Test=[];Test=[];Verif=[];Veri=[];WS_Te=[];Tes=[];
    Veri=Image_label(xxmin(1,i):xxmax(1,i),yymin(1,i):yymax(1,i),zzmin(1,i):zzmax(1,i));
    Verif=(Veri==L_NT_01(i));
    WS_Te=WS(xxmin(1,i):xxmax(1,i),yymin(1,i):yymax(1,i),zzmin(1,i):zzmax(1,i));
    Tes=Touch(xxmin(1,i):xxmax(1,i),yymin(1,i):yymax(1,i),zzmin(1,i):zzmax(1,i));
    WS_Test=WS_Te.*Verif;Test=Tes.*Verif;
    No=[];List=unique(nonzeros(WS_Test))';
    for j=1:size(List,2)
        No=[No nnz(WS_Test==List(1,j))];
    end
    No=No./sum(No);
    [N_a, idx]=sort(No,'descend');
    List_N{i}=List(:,idx);
    
    %
    Check=[];
    for l=1:max(find(N_a>0.1))
            WS_T0{l}=[];WS_TD{l}=[];
            WS_T0{l}=(WS_Test==List_N{i}(l));
            se=strel('sphere',3);
            WS_TD{l}=imdilate(WS_T0{l},se);
            WS_TD{l}=and(WS_TD{l},Test);
    end
    for m=2:max(find(N_a>0.1))
            Check=WS_TD{m-1}+WS_TD{m};
    end
    %
    Ch=Check>1;
    for n =1:max(find(N_a>0.1))
        WS_TF{n}=WS_TD{n}-and(Ch,WS_TD{n});
    end
    %
    for k=1:max(find(N_a>0.1))
        ind=[];x=[];y=[];z=[];X=[];Y=[];Z=[];
        ind=find(WS_TF{k});
        [x, y, z] = ind2sub(size(WS_Test), ind);
        Siz= size(x,1);

        x_c=mean(x);y_c=mean(y);z_c=mean(z);
        X_WS=[X_WS x_c+xxmin(1,i)-1];Y_WS=[Y_WS y_c+yymin(1,i)-1];Z_WS=[Z_WS z_c+zzmin(1,i)-1];

        Z=z'-ones(1,Siz).*z_c;X=x'-ones(1,Siz).*x_c;Y=y'-ones(1,Siz).*y_c;

        XX=sum(sum(X.^2)); XY=sum(sum(X.*Y)); XZ=sum(sum(X.*Z));
        YY=sum(sum(Y.^2)); YZ=sum(sum(Y.*Z)); ZZ=sum(sum(Z.^2));

        Moment_2=[XX XY XZ ; XY YY YZ ; XZ YZ ZZ];

        [V,Values_2] = eig(Moment_2);
        Dig=diag(Values_2);
        [D_a, idx]=sort(Dig,'ascend');
        Vec=V(:,idx);

        %
        Values=sqrt(D_a);

        Vecteur1_WS=[Vecteur1_WS Vec(:,1)];Vecteur2_WS=[Vecteur2_WS Vec(:,2)];Vecteur3_WS=[Vecteur3_WS Vec(:,3)];

        a=(12*D_a(1,1)/(Values(3)/Values(1))/(Values(2)/Values(1)))^0.2;
        b=a*Values(2)/Values(1);c=a*Values(3)/Values(1);

        a_WS=[a_WS a];b_WS=[b_WS b];c_WS=[c_WS c];
    end
end
%% Recover aspect ratios & orientation of separated particles
ba_WS=b_WS./a_WS;ca_WS=c_WS./a_WS;
L_ws=1:size(a_WS,2);
L_WS=find(L_ws-or(isnan(a_WS),or(a_WS<1.2,and(ba_WS<=1.2,ca_WS<=1.2))).*L_ws);
%
a_WS=a_WS(L_WS);b_WS=b_WS(L_WS);c_WS=c_WS(L_WS);
X_WS=X_WS(L_WS);Y_WS=Y_WS(L_WS);Z_WS=Z_WS(L_WS);
Vecteur1_WS=Vecteur1_WS(:,L_WS);Vecteur2_WS=Vecteur2_WS(:,L_WS);Vecteur3_WS=Vecteur3_WS(:,L_WS);
%% Save remaing touching particles
Touch_2=ismember(Image_label,L_TA_01);
Touch_Final=Touch_2;%Touch_add+
Directory='TR_';
for i=1:size(Touch_Final,3)
    g{i}=uint8(255*Touch_Final(:,:,i));
    imwrite(g{i},sprintf('%s%05d.tif', Directory, i))
end


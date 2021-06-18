clear;clc;close all
%% Original_Image reading
Directory = 'Images_L_500x500x200_';
%for i=n-1:-1:1
i1=1;
i2=200;
for i=i1:i2
fname{i-i1+1}= sprintf('%s%05d.tif', Directory, i-1);
I{i-i1+1}=imread(fname{i-i1+1});
end
%
for i=1:(i2-i1)+1% Number of 2Dimages
    OriginalImage(:,:,i)=I{1,i};
end
%
% Read segmented images
Directory = '';
i1=1;
i2=200;
for i=i1:i2
fname{i-i1+1}= sprintf('%s%05d.tif', Directory, i-1);
I{i-i1+1}=imread(fname{i-i1+1});
end
% 3D images
for i=1:(i2-i1)+1% Number of 2Dimages
    OI(:,:,i)=I{1,i};
end
%%
P_Fiber=(OI==250);%Pure fibers
Image_label=bwlabeln(P_Fiber,26);%Label particles
n=max(max(max(Image_label)))%number of particles
%% Calcualte volume, eigenvetors & eigenvalues
%Calculate nb of voxels of each particle belonging to corresponding equivalent parallelepiped: (V_in)
Volume=[];V_in_R=[]; V_eq_R=[];
V1=[]; V2=[]; V3=[];X_c=[];Y_c=[];Z_c=[];
Vecteur1=[];Vecteur3=[];Vecteur2=[];

for i = 1:n%number of particles
    
    ind=[];x=[];y=[];z=[];X=[];Y=[];Z=[];
    ind=find(Image_label==i);
    [x, y, z] = ind2sub(size(Image_label), ind);
    Siz= size(x,1);
    Volume=[Volume Siz];
    
    z_c=mean(z);x_c=mean(x);y_c=mean(y);
    X_c=[X_c x_c]; Y_c=[Y_c y_c]; Z_c=[Z_c z_c]; 

    Z=z'-ones(1,Siz).*z_c;X=x'-ones(1,Siz).*x_c;Y=y'-ones(1,Siz).*y_c;

    XX=sum(sum(X.^2)); XY=sum(sum(X.*Y)); XZ=sum(sum(X.*Z));
    YY=sum(sum(Y.^2)); YZ=sum(sum(Y.*Z)); ZZ=sum(sum(Z.^2));

    Moment_2=[XX XY XZ ; XY YY YZ ; XZ YZ ZZ];

    [V,Values_2] = eig(Moment_2);

    Dig=diag(Values_2);
    [D_a, idx]=sort(Dig,'ascend');
    Vec=V(:,idx);

    Values=sqrt(D_a);
  
    V1=[V1 Values(1)]; V2=[V2 Values(2)]; V3=[V3 Values(3)]; 
    
    Vecteur1=[Vecteur1 Vec(:,1)]; Vecteur2=[Vecteur2 Vec(:,2)]; Vecteur3=[Vecteur3 Vec(:,3)]; 
    
    a=(12*D_a(1,1)/(Values(3)/Values(1))/(Values(2)/Values(1)))^0.2;
    b=a*Values(2)/Values(1);c=a*Values(3)/Values(1);
    
    Verify_R=[];
    for ii = 1:size(X,2) 
        L=[X(1,ii) Y(1,ii) Z(1,ii)];
        L1=dot(L,Vec(:,1));
        L2=dot(L,Vec(:,2));
        L3=dot(L,Vec(:,3));
        if or(or(abs(L1)>a/2,abs(L2)>b/2),abs(L3)>c/2)
            Verify_R=[Verify_R 0];
        else Verify_R=[Verify_R 1];
        end
    end
    V_in_R=[V_in_R nnz(Verify_R)];
    V_eq_R=[V_eq_R a*b*c];
end
V11=(V1).^2;V22=(V2).^2;V33=(V3).^2;
%% Curves: 1-V_in/V_total
Verify_Volume_R=Volume-V_in_R;
ca=Verify_Volume_R./Volume;
r0=-0.05;
N=[];N_R=[];
for i = 1:20
    r=r0+i*0.05;
    R=and(ca>r,ca<r+0.05);
    V_R=R.*Volume;
    N=[N nnz(R)];
    N_R=[N_R sum(sum(V_R))];
end
L_E=0.025:0.05:0.975;
figure, 
yyaxis left
plot(L_E, N)
yyaxis right
plot(L_E, N_R)
%% rectangle parra
a_R=(12*V11./(V3./V1)./(V2./V1)).^0.2;
b_R=a_R.*(V2./V1);
c_R=a_R.*(V3./V1);
ca_R=V3./V1;
ba_R=V2./V1;
%% plot aspect ratio
abc=c_R;%./a_fib_all;
r0=0;
N=[];N_R=[];
for i = 1:199
    r=r0+i*1;
    R=and(abc>=r,abc<r+1);
    V_R=R.*Volume;
    N=[N nnz(R)];
    N_R=[N_R sum(V_R)];
end
L_E=1.5:1:199.5;
figure, 
yyaxis left
plot(L_E, N)
yyaxis right
plot(L_E, N_R)
%%
seuil=0.32;
L_ABC=find(or(isnan(a_R),or(a_R<1.5,and(ba_R<=1.2,ca_R<=1.2))));%find very small fibers (~some voxles) =>eliminate
L_All=1:n;
L_Fiber=L_All(find(L_All-L_All.*ismember(L_All,L_ABC)));%nonzeros(L_All-L_All.*ismember(L_All,L_ABC));
L_T=find(Verify_Volume_R./Volume>seuil);%threshold to distinguish betweem touching particles and non
L_TA=nonzeros(L_T-L_T.*ismember(L_T,L_ABC));
L_NT=nonzeros(L_Fiber-L_Fiber.*ismember(L_Fiber,L_TA));
%% Verify whether touching particles?
Touch=ismember(Image_label,L_TA);
T_label=Image_label.*double(Touch);
imtool(T_label(:,:,100))
imtool(permute(T_label(:,250,:),[1 3 2]))
%% Center cordonates & aspect ratios + orientations of non-touching particles
X_cen=X_c(L_NT);Y_cen=Y_c(L_NT);Z_cen=Z_c(L_NT);
R_a=a_R(L_NT);R_b=b_R(L_NT);R_c=c_R(L_NT);
Vecteur1_a=Vecteur1(:,L_NT);Vecteur2_b=Vecteur2(:,L_NT);Vecteur3_c=Vecteur3(:,L_NT);
%% Saving touching particles
Directory='';
for i=1:size(Touch,3)
    g{i}=uint8(255*Touch(:,:,i));
    imwrite(g{i},sprintf('%s%05d.tif', Directory, i))
end
%%
NonTouch_PF=ismember(Image_label,L_NT);%non-touching particles
imtool(NonTouch_PF(:,:,100))
imtool(permute(NonTouch_PF(:,250,:),[1 3 2]))

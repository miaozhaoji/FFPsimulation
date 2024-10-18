close all;clc;clear all;
tic;                                                                        %计算tic~toc之间的执行时间
%% 系统参数
kB = 1.3806488e-23;     TT = 26.5+273.15;     u0 = 4*pi*1e-7;                	%外部常量：玻尔兹曼常数kB，绝对温度（室温)TT，真空磁导率u0
D = 30e-9;      V=D^3*pi/6;     Ms = 0.6/u0;    m=Ms*V;                     %样品信息：个数浓度c，粒径D，粒子体积V，粒子饱和磁化强度(Ms=0.6T/u0)，粒子磁矩m

Fs=1e6;     Fn=1e4;     t=0.1/Fs:1/Fs:Fn/Fs-0.1/Fs;      Len_t=Fn;                  	%采样信息：采样率Fs，采样数Fn,离散时间点t（单位s）
Gx=3/u0;fx=24.51e3;   Tx=1/fx;  Ax=12e-3/u0;  xmax=Ax/Gx;                                
Gy=3/u0;fy=26.04e3;   Ty=1/fy;  Ay=12e-3/u0;  ymax=Ay/Gy;            

G=[Gx;Gy];
K=u0*m/(kB*TT);                                                             %K：郎之万方程的比例系数   
Hb=1e-3*1e-4/u0;
Hsat=100/K;                                                                 %单位A/m
toc;
%% 成像空间参数
step=xmax/16;                                                    %样品分布空间离散点的最小间距为0.1mm
x=-xmax+step/2:step:xmax;	n_c=length(x);	x=x';	Len_x=length(x);   	%空间信息：2D成像空间，长度单位m
y=-ymax+step/2:step:ymax;	n_r=length(y);  y=y';	Len_y=length(y);   	c0=1;
c=zeros(n_r,n_c); n=0;
aa=double(imresize(imread('10narrowblood_1.png'),[32,32]));
c=aa/max(max(aa));

% c(3:15,7)=ones();c(3:15,15)=ones();

figure(1);
imagesc(x*1e3,y*1e3,c(end:-1:1,:));axis xy;axis image;set(gca,'FontSize',15,'tickdir','out');
axis([-xmax xmax -ymax ymax]*1e3);
xlabel('X/mm','Fontsize',20);ylabel('Y/mm','Fontsize',20);
colormap(gray);

toc;
%% 信号生成
%% 系统函数标定以及信号生成
Hx_t=Ax*cos(2*pi*fx*t);             Hy_t=Ay*cos(2*pi*fy*t);                                         
H_t=[Hx_t;Hy_t]; 
H_Amp_t=sqrt(sum(H_t.^2));  H_Dir_t=bsxfun(@rdivide,H_t,H_Amp_t);          	%H_t_Amp和H_t_Dir分别是H_t的幅值和方向
dHx_t=-Ax*2*pi*fx*sin(2*pi*fx*t);   dHy_t=-Ay*2*pi*fy*sin(2*pi*fy*t);
dH_t=[dHx_t;dHy_t];                                     %H_t和dH_t表示2维向量，，每一列代表一个离散时刻的激励合磁场向量H(t)和dH(t)/dt
dH_Amp_t=sqrt(sum(dH_t.^2));  dH_Dir_t=bsxfun(@rdivide,dH_t,dH_Amp_t);          	%H_t_Amp和H_t_Dir分别是H_t的幅值和方向

%求101*101个理想点的准确坐标

[X, Y] = meshgrid(x, y);
x_cord = X(:);  % 将 X 转换为列向量
y_cord = Y(:);  % 将 Y 转换为列向量

X_cord=[x_cord';y_cord'];     
Hs=bsxfun(@times,G,X_cord);
c_reshape=reshape(c,length(c(:)),1);
SYM_t=zeros(Len_t,length(c(:)),2);
dr2_t=K*dH_t;
SIG_t=zeros(Len_t,2);
for i=1:length(c(:))
            disp(i/length(c(:)))
            r2=K*(bsxfun(@minus,H_t,Hs(:,i)));       r2_Mag=sqrt(sum(r2.^2));                      %合磁场=H(t)-GX
%             L_r2Mag=coth(r2_Mag)-1./r2_Mag; dL_r2Mag=1./(r2_Mag.^2)-(coth(r2_Mag)).^2+1;
            [dL_r2Mag, L_r2Mag] = Langevin(r2_Mag);
            dr2_t_parl=bsxfun(@times,sum(dr2_t.*r2)./(r2_Mag.^2),r2);
            dr2_t_perp=dr2_t-dr2_t_parl; 
            temp=c0*m*1e10*(bsxfun(@times,dL_r2Mag,dr2_t_parl)+bsxfun(@times,L_r2Mag./r2_Mag,dr2_t_perp));
%             temp1 = awgn(temp,10,'measured');
            SYM_t(:,i,:)=temp';
%             figure(1)
%             plot(temp(1,:).');
            if c_reshape(i)~=0 
                SIG_t=SIG_t+c0*c_reshape(i)*m*1e10*(bsxfun(@times,dL_r2Mag,dr2_t_parl)+bsxfun(@times,L_r2Mag./r2_Mag,dr2_t_perp))'; 
            end
end
toc;
%% 信号处理

temp1=fft(awgn(SYM_t,40,'measured'));


%sym1=awgn(SYM_t,10,'measured');
%sym2=awgn(SYM_t,20,'measured');
%sym3=awgn(SYM_t,30,'measured');
%temp_2 = fft(sym1);
%ttt = squeeze(temp_2(:,100,1));
%temp_1 = fft(sym1);
%temp_2 = fft(sym2);
%temp_3 = fft(sym3);
% data=importdata('D:\大孔径实测\0609_1mg\20220609\202206091.txt');
% temp_big = SYM_t(:,:,1)+data;
% clear SYM_t;




SYM = zeros(0.5*Len_t,2,32*32);
SYM(:,1,:) = real(temp1(2:end/2+1,:,1));
SYM(:,2,:) = imag(temp1(2:end/2+1,:,1));
SYM = reshape(SYM,[0.5*Len_t,2,32,32]);

writeNPY(SYM,'7.npy')


SYM_f=zeros(Len_t,length(c(:)));            
SYM_f(1:0.5*Len_t,:)=real(temp1(2:end/2+1,:,1));
SYM_f(0.5*Len_t+1:Len_t,:)=imag(temp1(2:end/2+1,:,1));

%for n = 1:20
%figure
%colormap(gray);
%imagesc(squeeze(SYM(n*fix(fx/100),:,:)));
%end

a = 1;
%save('sm_select_3t_size_30nm.h5','SYM','-v7.3');

%SYM_1f=zeros(Len_t,length(c(:)));            %带噪系统矩阵
%SYM_1f(1:0.5*Len_t,:)=real(temp_1(2:end/2+1,:,1));
%SYM_1f(0.5*Len_t+1:Len_t,:)=imag(temp_1(2:end/2+1,:,1));

%SYM_2f=zeros(Len_t,length(c(:)));            %带噪系统矩阵
%SYM_2f(1:0.5*Len_t,:)=real(temp_2(2:end/2+1,:,1));
%SYM_2f(0.5*Len_t+1:Len_t,:)=imag(temp_2(2:end/2+1,:,1));

%SYM_3f=zeros(Len_t,length(c(:)));            %带噪系统矩阵
%SYM_3f(1:0.5*Len_t,:)=real(temp_3(2:end/2+1,:,1));
%SYM_3f(0.5*Len_t+1:Len_t,:)=imag(temp_3(2:end/2+1,:,1));

%snr_1f = SYM_1f ./ (SYM_1f-SYM_f);
%SNR_1P =abs(snr_1f(:,10));

%figure
%plot(SNR_1P(1:50000))

%snr_2f = SYM_2f ./ (SYM_2f-SYM_f);
%SNR_2P =abs(snr_2f(:,10));
%figure
%plot(SNR_2P(1:50000))

%snr_3f = SYM_3f ./ (SYM_3f-SYM_f);
%SNR_3P =abs(snr_1f(:,10));
%figure
%plot(SNR_3P(1:50000))


% rSNR=10*log10(mean(abs(SYM_f').^2)./mean(abs(SYM_2f'-SYM_f').^2)); 
% rSNR_1f=abs(mean(SYM_1f./(SYM_1f-SYM_f),2));
% rSNR_2f=abs(mean(SYM_2f./(SYM_2f-SYM_f),2));
% rSNR_3f=abs(mean(SYM_3f./(SYM_3f-SYM_f),2));
% SNR_1f = rSNR_1f>3;
% SNR_2f = rSNR_2f>3;
% SNR_3f = rSNR_3f>3;
% S_SNR_1f = SYM_1f(SNR_1f(:),:);
% SYM_1fc = SYM_f(SNR_1f(:),:);
% S_SNR_2f = SYM_2f(SNR_2f(:),:);
% SYM_2fc = SYM_f(SNR_1f(:),:);
% S_SNR_3f = SYM_3f(SNR_3f(:),:);
% SYM_3fc = SYM_f(SNR_1f(:),:);
% SYM_f = SYM_f(SNR(:),:);
% SYM_1f = SYM_1f(SNR(:),:);
% SYM_2f = SYM_2f(SNR(:),:);
% SYM_3f = SYM_3f(SNR(:),:);

% SIG_f = SYM_f * reshape(c,19*19,1);
%SIG_f = SYM_f * reshape(c,32*32,1);

%SYM_f_ = reshape(SYM_f,[4000,1,40,40]);
%SYM_f_ = repmat(SYM_f_, [1, 1, 1, 1, 40]);

%save('sm_select_3t_size_30nm.h5','SYM_f_','-v7.3');

%SYM_f_ = SYM_f_(:,:,:,:,1);
%SYM_f = reshape(SYM_f_,[4000,40*40]);
%save('sm_select_1t_size_60nm_tem_37.5.h5','SYM_f','-v7.3');

temp=fft(SIG_t);
SIG_f=zeros(Len_t,1);                     %单维度线圈
SIG_f(1:0.5*Len_t,1)=real(temp(2:end/2+1,1));
SIG_f(0.5*Len_t+1:Len_t,1)=imag(temp(2:end/2+1,1));
%SIG_f(Len_t+1:1.5*Len_t,1)=real(temp(2:end/2+1,2));
%SIG_f(1.5*Len_t+1:2*Len_t,1)=imag(temp(2:end/2+1,2));


% % 
% SIG_f=zeros(Len_t,1);                     %单维度线圈
% SIG_f(1:0.5*Len_t,1)=temp2(2:end/2+1,1);
% SIG_f(0.5*Len_t+1:Len_t,1)=temp2(2:end/2+1,2);
% % SIG_fr=SIG_f(SNR(:),:);
% for i = 1:5
%     img=SYM_1f(i*1000,:);
%     img=reshape(img,32,32);
%     img=abs(img);
%     figure
%     colormap gray;
%     imagesc(img)
%     
%     img=SYM_f(i*1000,:);
%     img=reshape(img,32,32);
%     img=abs(img);
%     figure
%     colormap gray;
%     imagesc(img)
% end
%% 图像重建:反问题求解
% theta=1e-8;
% % ma =size(S_SNR);
% SYM=SYM_f'*SYM_f+theta*eye(size(c(:)));
% SIG=SYM_f'*SIG_f;
% 
% tol=1e-10;
% maxit=500;
% [c_solution,flag,relres,iter,resvec]=cgs(SYM,SIG,tol,maxit);
% 
% SYM2=SYM_1f'*SYM_1f+theta*eye(size(c(:)));
% SIG2=SYM_1f'*SIG_f;
% [c_solution2,flag2,relres2,iter2,resvec2]=cgs(SYM2,SIG2,tol,maxit);
% 
% SYM3=SYM_2f'*SYM_2f+theta*eye(size(c(:)));
% SIG3=SYM_2f'*SIG_f;
% [c_solution3,flag3,relres3,iter3,resvec3]=cgs(SYM3,SIG3,tol,maxit);
% 
% SYM4=SYM_3f'*SYM_3f+theta*eye(size(c(:)));
% SIG4=SYM_3f'*SIG_f;
% [c_solution4,flag4,relres4,iter4,resvec4]=cgs(SYM4,SIG4,tol,maxit);
iterations=10;
lambd=1e-5;


[c_solution1]= kaczmarzReg((SYM_f).',(SIG_f),iterations,lambd,0,1,1);

%[c_solution2]= kz_select((SYM_1f).',(SIG_f),iterations,lambd,0,1,1, 1e-5,32);

%[c_solution3]= kz_select((SYM_2f).',(SIG_f),iterations,lambd,0,1,1, 1e-5,32);

%[c_solution4]= kz_select((SYM_3f).',(SIG_f),iterations,lambd,0,1,1, 1e-5,32);

figure(1)
imagesc(c_solution1)
%figure(2)
%imagesc(c_solution2)
%figure(3)
%imagesc(c_solution3)
%figure(4)
%imagesc(c_solution4)

%% 重建图像
figure(2);
subplot(1,2,1);
plot(c_reshape);
subplot(1,2,2);
plot(c_solution1);
c_reconstruct1=reshape(c_solution1,size(c));
%c_reconstruct2=reshape(c_solution2,size(c));
%c_reconstruct3=reshape(c_solution3,size(c));
%c_reconstruct4=reshape(c_solution4,size(c));
figure(3);
subplot(1,2,1);
imagesc(x*1e3,y*1e3,c(end:-1:1,:));axis xy;axis image;set(gca,'FontSize',15,'tickdir','out');
axis([-xmax xmax -ymax ymax]*1e3);
xlabel('x/mm','Fontsize',20);ylabel('z/mm','Fontsize',20);title('原始图像');%set(gca, 'XDir','reverse');%axis square;
colormap(gray);
subplot(1,2,2);

imagesc(x*1e3,y*1e3,c_reconstruct1(end:-1:1,:));axis xy;axis image;set(gca, 'FontSize', 15,'tickdir','out'); 
axis([-xmax xmax -ymax ymax]*1e3);
xlabel('x/mm','Fontsize',20);ylabel('z/mm','Fontsize',20);
title('幅值重建图像(CGS)');%set(gca, 'XDir','reverse');%axis square;
colormap(gray);


% 计算PSNR
psnr_value = psnr(c, c_reconstruct1);

% 计算SSIM
[ssim_value, ssim_map] = ssim(c, c_reconstruct1);

fprintf('PSNR: %f dB\n', psnr_value);
fprintf('SSIM: %f\n', ssim_value);
%set(gca,'xcolor','black','xtick',-inf:inf:inf);
%set(gca,'ycolor','black','ytick',-inf:inf:inf);
%set(gcf,'Position',get(0,'screensize')); 
%figure(6);
%imagesc(x*1e3,y*1e3,c_reconstruct2(end:-1:1,:));axis xy;axis image;set(gca, 'FontSize', 15,'tickdir','out'); 
%axis([-xmax xmax -ymax ymax]*1e3);
% xlabel('x/mm','Fontsize',20);ylabel('z/mm','Fontsize',20);
%title('幅值重建图像(CGS)');%set(gca, 'XDir','reverse');%axis square;
%colormap(gray);
%set(gca,'xcolor','black','xtick',-inf:inf:inf);
%set(gca,'ycolor','black','ytick',-inf:inf:inf);
% set(gcf,'Position',get(0,'screensize')); 
%figure(7);
%imagesc(x*1e3,y*1e3,c_reconstruct3(end:-1:1,:));axis xy;axis image;set(gca, 'FontSize', 15,'tickdir','out'); 
%axis([-xmax xmax -ymax ymax]*1e3);
% xlabel('x/mm','Fontsize',20);ylabel('z/mm','Fontsize',20);
%title('幅值重建图像(CGS)');%set(gca, 'XDir','reverse');%axis square;
%colormap(gray);
%set(gca,'xcolor','black','xtick',-inf:inf:inf);
%set(gca,'ycolor','black','ytick',-inf:inf:inf);
% set(gcf,'Position',get(0,'screensize')); 
%figure(8);
%imagesc(x*1e3,y*1e3,c_reconstruct4(end:-1:1,:));axis xy;axis image;set(gca, 'FontSize', 15,'tickdir','out'); 
%axis([-xmax xmax -ymax ymax]*1e3);
% xlabel('x/mm','Fontsize',20);ylabel('z/mm','Fontsize',20);
%title('幅值重建图像(CGS)');%set(gca, 'XDir','reverse');%axis square;
%colormap(gray);
%set(gca,'xcolor','black','xtick',-inf:inf:inf);
%set(gca,'ycolor','black','ytick',-inf:inf:inf);
% set(gcf,'Position',get(0,'screensize')); 
%figure(4);
%surf(x*1e3,y*1e3,c_reconstruct1./max(c_reconstruct1(:)));view(3);hold on;
%xlabel('x/mm','Fontsize',20);ylabel('z/mm','Fontsize',20);title('幅值重建图像(CGS)');%set(gca, 'XDir','reverse');%axis square;
%colormap(hot);
%x_MAT=repmat(x',length(y),1)*1e3;
%y_MAT=repmat(y,1,length(x))*1e3;
%hold on;
%plot3(x_MAT,y_MAT,c./max(c(:)),'b*');%axis image;
%hold off;
toc;

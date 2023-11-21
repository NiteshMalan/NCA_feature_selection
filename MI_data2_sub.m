clc; clear; close all;
Fs=250;
[s, h] = sload('B0403T.gdf'); %load the data, need biosig toolbox to load gdf data
cl_left=find(h.EVENT.TYP==769);
cl_right=find(h.EVENT.TYP==770);
N_trails=length(cl_left);
posCL_left=h.EVENT.POS(cl_left);
for k=1:length(cl_left)
X{k}=s((posCL_left(k)+305:posCL_left(k)+1072),:);
X{k}(isnan(X{k}))=0;
end

posCL_right=h.EVENT.POS(cl_right);
for k1=1:length(cl_left)
Y{k1}=s((posCL_right(k1)+305:posCL_right(k1)+1072),:);
Y{k1}(isnan(Y{k1}))=0;
end

%% Temporal Filtering

for dq=1:N_trails;
    L=4;
X_uc3{dq}= dtcwtu(X{dq}(:,1),L);
X_ucz{dq}= dtcwtu(X{dq}(:,2),L);
X_uc4{dq}= dtcwtu(X{dq}(:,3),L);

Y_uc3{dq}= dtcwtu(Y{dq}(:,1),L);
Y_ucz{dq}= dtcwtu(Y{dq}(:,2),L);
Y_uc4{dq}= dtcwtu(Y{dq}(:,3),L);

X_u2c3{dq}= dtcwtu2(X{dq}(:,1),L);
X_u2cz{dq}= dtcwtu2(X{dq}(:,2),L);
X_u2c4{dq}= dtcwtu2(X{dq}(:,3),L);

Y_u2c3{dq}= dtcwtu2(Y{dq}(:,1),L);
Y_u2cz{dq}= dtcwtu2(Y{dq}(:,2),L);
Y_u2c4{dq}= dtcwtu2(Y{dq}(:,3),L);

X_bc3{dq}= dtcwtb(X{dq}(:,1),L);
X_bcz{dq}= dtcwtb(X{dq}(:,2),L);
X_bc4{dq}= dtcwtb(X{dq}(:,3),L);

Y_bc3{dq}= dtcwtb(Y{dq}(:,1),L);
Y_bcz{dq}= dtcwtb(Y{dq}(:,2),L);
Y_bc4{dq}= dtcwtb(Y{dq}(:,3),L);

meanuxc3{dq}=mean(X_uc3{dq});
meanuxcz{dq}=mean(X_ucz{dq});
meanuxc4{dq}=mean(X_uc4{dq});

meanuyc3{dq}=mean(Y_uc3{dq});
meanuycz{dq}=mean(Y_ucz{dq});
meanuyc4{dq}=mean(Y_uc4{dq});

meanu2xc3{dq}=mean(X_u2c3{dq});
meanu2xcz{dq}=mean(X_u2cz{dq});
meanu2xc4{dq}=mean(X_u2c4{dq});

meanu2yc3{dq}=mean(Y_u2c3{dq});
meanu2ycz{dq}=mean(Y_u2cz{dq});
meanu2yc4{dq}=mean(Y_u2c4{dq});

meanbxc3{dq}=mean(X_bc3{dq});
meanbxcz{dq}=mean(X_bcz{dq});
meanbxc4{dq}=mean(X_bc4{dq});

meanbyc3{dq}=mean(Y_bc3{dq});
meanbycz{dq}=mean(Y_bcz{dq});
meanbyc4{dq}=mean(Y_bc4{dq});

stduxc3{dq}=std(X_uc3{dq});
stduxcz{dq}=std(X_ucz{dq});
stduxc4{dq}=std(X_uc4{dq});

stduyc3{dq}=std(Y_uc3{dq});
stduycz{dq}=std(Y_ucz{dq});
stduyc4{dq}=std(Y_uc4{dq});

stdu2xc3{dq}=std(X_u2c3{dq});
stdu2xcz{dq}=std(X_u2cz{dq});
stdu2xc4{dq}=std(X_u2c4{dq});

stdu2yc3{dq}=std(Y_u2c3{dq});
stdu2ycz{dq}=std(Y_u2cz{dq});
stdu2yc4{dq}=std(Y_u2c4{dq});

stdbxc3{dq}=std(X_bc3{dq});
stdbxcz{dq}=std(X_bcz{dq});
stdbxc4{dq}=std(X_bc4{dq});

stdbyc3{dq}=std(Y_bc3{dq});
stdbycz{dq}=std(Y_bcz{dq});
stdbyc4{dq}=std(Y_bc4{dq});

varuxc3{dq}=var(X_uc3{dq});
varuxcz{dq}=var(X_ucz{dq});
varuxc4{dq}=var(X_uc4{dq});

varuyc3{dq}=var(Y_uc3{dq});
varuycz{dq}=var(Y_ucz{dq});
varuyc4{dq}=var(Y_uc4{dq});

varu2xc3{dq}=var(X_u2c3{dq});
varu2xcz{dq}=var(X_u2cz{dq});
varu2xc4{dq}=var(X_u2c4{dq});

varu2yc3{dq}=var(Y_u2c3{dq});
varu2ycz{dq}=var(Y_u2cz{dq});
varu2yc4{dq}=var(Y_u2c4{dq});

varbxc3{dq}=var(X_bc3{dq});
varbxcz{dq}=var(X_bcz{dq});
varbxc4{dq}=var(X_bc4{dq});

varbyc3{dq}=var(Y_bc3{dq});
varbycz{dq}=var(Y_bcz{dq});
varbyc4{dq}=var(Y_bc4{dq});

N = length(X_uc3{dq});
dftuxc3{dq} = fft(X_uc3{dq});
dftuxc3{dq} = dftuxc3{dq}(1:N/2+1);
psduxc3{dq} = (1/(Fs*N)) * abs(dftuxc3{dq}).^2;
psduxc3{dq}(2:end-1) = 2*psduxc3{dq}(2:end-1);
psduxc3{dq}=sum(psduxc3{dq});

dftuxcz{dq} = fft(X_ucz{dq});
dftuxcz{dq} = dftuxcz{dq}(1:N/2+1);
psduxcz{dq} = (1/(Fs*N)) * abs(dftuxcz{dq}).^2;
psduxcz{dq}(2:end-1) = 2*psduxcz{dq}(2:end-1);
psduxcz{dq}=sum(psduxcz{dq});

dftuxc4{dq} = fft(X_uc4{dq});
dftuxc4{dq} = dftuxc4{dq}(1:N/2+1);
psduxc4{dq} = (1/(Fs*N)) * abs(dftuxc4{dq}).^2;
psduxc4{dq}(2:end-1) = 2*psduxc4{dq}(2:end-1);
psduxc4{dq}=sum(psduxc4{dq});

dftu2xc3{dq} = fft(X_u2c3{dq});
dftu2xc3{dq} = dftu2xc3{dq}(1:N/2+1);
psdu2xc3{dq} = (1/(Fs*N)) * abs(dftu2xc3{dq}).^2;
psdu2xc3{dq}(2:end-1) = 2*psdu2xc3{dq}(2:end-1);
psdu2xc3{dq}=sum(psdu2xc3{dq});

dftu2xcz{dq} = fft(X_u2cz{dq});
dftu2xcz{dq} = dftu2xcz{dq}(1:N/2+1);
psdu2xcz{dq} = (1/(Fs*N)) * abs(dftu2xcz{dq}).^2;
psdu2xcz{dq}(2:end-1) = 2*psdu2xcz{dq}(2:end-1);
psdu2xcz{dq}=sum(psdu2xcz{dq});

dftu2xc4{dq} = fft(X_u2c4{dq});
dftu2xc4{dq} = dftu2xc4{dq}(1:N/2+1);
psdu2xc4{dq} = (1/(Fs*N)) * abs(dftu2xc4{dq}).^2;
psdu2xc4{dq}(2:end-1) = 2*psdu2xc4{dq}(2:end-1);
psdu2xc4{dq}=sum(psdu2xc4{dq});

dftbxc3{dq} = fft(X_bc3{dq});
dftbxc3{dq} = dftbxc3{dq}(1:N/2+1);
psdbxc3{dq} = (1/(Fs*N)) * abs(dftbxc3{dq}).^2;
psdbxc3{dq}(2:end-1) = 2*psdbxc3{dq}(2:end-1);
psdbxc3{dq}=sum(psdbxc3{dq});

dftbxcz{dq} = fft(X_bcz{dq});
dftbxcz{dq} = dftbxcz{dq}(1:N/2+1);
psdbxcz{dq} = (1/(Fs*N)) * abs(dftbxcz{dq}).^2;
psdbxcz{dq}(2:end-1) = 2*psdbxcz{dq}(2:end-1);
psdbxcz{dq}=sum(psdbxcz{dq});

dftbxc4{dq} = fft(X_bc4{dq});
dftbxc4{dq} = dftbxc4{dq}(1:N/2+1);
psdbxc4{dq} = (1/(Fs*N)) * abs(dftbxc4{dq}).^2;
psdbxc4{dq}(2:end-1) = 2*psdbxc4{dq}(2:end-1);
psdbxc4{dq}=sum(psdbxc4{dq});

dftuyc3{dq} = fft(Y_uc3{dq});
dftuyc3{dq} = dftuyc3{dq}(1:N/2+1);
psduyc3{dq} = (1/(Fs*N)) * abs(dftuyc3{dq}).^2;
psduyc3{dq}(2:end-1) = 2*psduyc3{dq}(2:end-1);
psduyc3{dq}=sum(psduyc3{dq});

dftuycz{dq} = fft(Y_ucz{dq});
dftuycz{dq} = dftuycz{dq}(1:N/2+1);
psduycz{dq} = (1/(Fs*N)) * abs(dftuycz{dq}).^2;
psduycz{dq}(2:end-1) = 2*psduycz{dq}(2:end-1);
psduycz{dq}=sum(psduycz{dq});

dftuyc4{dq} = fft(Y_uc4{dq});
dftuyc4{dq} = dftuyc4{dq}(1:N/2+1);
psduyc4{dq} = (1/(Fs*N)) * abs(dftuyc4{dq}).^2;
psduyc4{dq}(2:end-1) = 2*psduyc4{dq}(2:end-1);
psduyc4{dq}=sum(psduyc4{dq});

dftu2yc3{dq} = fft(Y_u2c3{dq});
dftu2yc3{dq} = dftu2yc3{dq}(1:N/2+1);
psdu2yc3{dq} = (1/(Fs*N)) * abs(dftu2yc3{dq}).^2;
psdu2yc3{dq}(2:end-1) = 2*psdu2yc3{dq}(2:end-1);
psdu2yc3{dq}=sum(psdu2yc3{dq});

dftu2ycz{dq} = fft(Y_u2cz{dq});
dftu2ycz{dq} = dftu2ycz{dq}(1:N/2+1);
psdu2ycz{dq} = (1/(Fs*N)) * abs(dftu2ycz{dq}).^2;
psdu2ycz{dq}(2:end-1) = 2*psdu2ycz{dq}(2:end-1);
psdu2ycz{dq}=sum(psdu2ycz{dq});

dftu2yc4{dq} = fft(Y_u2c4{dq});
dftu2yc4{dq} = dftu2yc4{dq}(1:N/2+1);
psdu2yc4{dq} = (1/(Fs*N)) * abs(dftu2yc4{dq}).^2;
psdu2yc4{dq}(2:end-1) = 2*psdu2yc4{dq}(2:end-1);
psdu2yc4{dq}=sum(psdu2yc4{dq});

dftbyc3{dq} = fft(Y_bc3{dq});
dftbyc3{dq} = dftbyc3{dq}(1:N/2+1);
psdbyc3{dq} = (1/(Fs*N)) * abs(dftbyc3{dq}).^2;
psdbyc3{dq}(2:end-1) = 2*psdbyc3{dq}(2:end-1);
psdbyc3{dq}=sum(psdbyc3{dq});

dftbycz{dq} = fft(Y_bcz{dq});
dftbycz{dq} = dftbycz{dq}(1:N/2+1);
psdbycz{dq} = (1/(Fs*N)) * abs(dftbycz{dq}).^2;
psdbycz{dq}(2:end-1) = 2*psdbycz{dq}(2:end-1);
psdbycz{dq}=sum(psdbycz{dq});

dftbyc4{dq} = fft(Y_bc4{dq});
dftbyc4{dq} = dftbyc4{dq}(1:N/2+1);
psdbyc4{dq} = (1/(Fs*N)) * abs(dftbyc4{dq}).^2;
psdbyc4{dq}(2:end-1) = 2*psdbyc4{dq}(2:end-1);
psdbyc4{dq}=sum(psdbyc4{dq});

sampuxc3{dq}=SampEn(2,0.2*std(X_uc3{dq}),X_uc3{dq});
sampuxcz{dq}=SampEn(2,0.2*std(X_ucz{dq}),X_uc4{dq});
sampuxc4{dq}=SampEn(2,0.2*std(X_uc4{dq}),X_uc4{dq});

sampuyc3{dq}=SampEn(2,0.2*std(Y_uc3{dq}),Y_uc3{dq});
sampuycz{dq}=SampEn(2,0.2*std(Y_ucz{dq}),Y_ucz{dq});
sampuyc4{dq}=SampEn(2,0.2*std(Y_uc4{dq}),Y_uc4{dq});

sampu2xc3{dq}=SampEn(2,0.2*std(X_u2c3{dq}),X_u2c3{dq});
sampu2xcz{dq}=SampEn(2,0.2*std(X_u2cz{dq}),X_u2c4{dq});
sampu2xc4{dq}=SampEn(2,0.2*std(X_u2c4{dq}),X_u2c4{dq});

sampu2yc3{dq}=SampEn(2,0.2*std(Y_u2c3{dq}),Y_u2c3{dq});
sampu2ycz{dq}=SampEn(2,0.2*std(Y_u2cz{dq}),Y_u2cz{dq});
sampu2yc4{dq}=SampEn(2,0.2*std(Y_u2c4{dq}),Y_u2c4{dq});

sampbxc3{dq}=SampEn(2,0.2*std(X_bc3{dq}),X_bc3{dq});
sampbxcz{dq}=SampEn(2,0.2*std(X_bcz{dq}),X_bc4{dq});
sampbxc4{dq}=SampEn(2,0.2*std(X_bc4{dq}),X_bc4{dq});

sampbyc3{dq}=SampEn(2,0.2*std(Y_bc3{dq}),Y_bc3{dq});
sampbycz{dq}=SampEn(2,0.2*std(Y_bcz{dq}),Y_bcz{dq});
sampbyc4{dq}=SampEn(2,0.2*std(Y_bc4{dq}),Y_bc4{dq});

uxhc3{dq} = hilbert(X_uc3{dq});
uxhcz{dq} = hilbert(X_ucz{dq});
uxhc4{dq} = hilbert(X_uc4{dq});
uxc3phase{dq} = (unwrap(angle(uxhc3{dq})))';
uxczphase{dq} = (unwrap(angle(uxhcz{dq})))';
uxc4phase{dq} = (unwrap(angle(uxhc4{dq})))';

uyhc3{dq} = hilbert(Y_uc3{dq});
uyhcz{dq} = hilbert(Y_ucz{dq});
uyhc4{dq} = hilbert(Y_uc4{dq});
uyc3phase{dq} = (unwrap(angle(uyhc3{dq})))';
uyczphase{dq} = (unwrap(angle(uyhcz{dq})))';
uyc4phase{dq} = (unwrap(angle(uyhc4{dq})))';

u2xhc3{dq} = hilbert(X_u2c3{dq});
u2xhcz{dq} = hilbert(X_u2cz{dq});
u2xhc4{dq} = hilbert(X_u2c4{dq});
u2xc3phase{dq} = (unwrap(angle(u2xhc3{dq})))';
u2xczphase{dq} = (unwrap(angle(u2xhcz{dq})))';
u2xc4phase{dq} = (unwrap(angle(u2xhc4{dq})))';

u2yhc3{dq} = hilbert(Y_u2c3{dq});
u2yhcz{dq} = hilbert(Y_u2cz{dq});
u2yhc4{dq} = hilbert(Y_u2c4{dq});
u2yc3phase{dq} = (unwrap(angle(u2yhc3{dq})))';
u2yczphase{dq} = (unwrap(angle(u2yhcz{dq})))';
u2yc4phase{dq} = (unwrap(angle(u2yhc4{dq})))';

plvux1{dq}=(abs(sum(exp(-i*(uxc3phase{dq}-uxczphase{dq})))))/758;
plvux2{dq}=(abs(sum(exp(-i*(uxc4phase{dq}-uxczphase{dq})))))/758;

plvuy1{dq}=(abs(sum(exp(-i*(uyc3phase{dq}-uyczphase{dq})))))/758;
plvuy2{dq}=(abs(sum(exp(-i*(uyc4phase{dq}-uyczphase{dq})))))/758;

plvu2x1{dq}=(abs(sum(exp(-i*(u2xc3phase{dq}-u2xczphase{dq})))))/758;
plvu2x2{dq}=(abs(sum(exp(-i*(u2xc4phase{dq}-u2xczphase{dq})))))/758;

plvu2y1{dq}=(abs(sum(exp(-i*(u2yc3phase{dq}-u2yczphase{dq})))))/758;
plvu2y2{dq}=(abs(sum(exp(-i*(u2yc4phase{dq}-u2yczphase{dq})))))/758;

bxhc3{dq} = hilbert(X_bc3{dq});
bxhcz{dq} = hilbert(X_bcz{dq});
bxhc4{dq} = hilbert(X_bc4{dq});
bxc3phase{dq} = (unwrap(angle(bxhc3{dq})))';
bxczphase{dq} = (unwrap(angle(bxhcz{dq})))';
bxc4phase{dq} = (unwrap(angle(bxhc4{dq})))';

byhc3{dq} = hilbert(Y_bc3{dq});
byhcz{dq} = hilbert(Y_bcz{dq});
byhc4{dq} = hilbert(Y_bc4{dq});
byc3phase{dq} = (unwrap(angle(byhc3{dq})))';
byczphase{dq} = (unwrap(angle(byhcz{dq})))';
byc4phase{dq} = (unwrap(angle(byhc4{dq})))';

plvbx1{dq}=(abs(sum(exp(-i*(bxc3phase{dq}-bxczphase{dq})))))/758;
plvbx2{dq}=(abs(sum(exp(-i*(bxc4phase{dq}-bxczphase{dq})))))/758;

plvby1{dq}=(abs(sum(exp(-i*(byc3phase{dq}-byczphase{dq})))))/758;
plvby2{dq}=(abs(sum(exp(-i*(byc4phase{dq}-byczphase{dq})))))/758;

rmsuxc3{dq}=rms(X_uc3{dq});
rmsuxcz{dq}=rms(X_ucz{dq});
rmsuxc4{dq}=rms(X_uc4{dq});

rmsuyc3{dq}=rms(Y_uc3{dq});
rmsuycz{dq}=rms(Y_ucz{dq});
rmsuyc4{dq}=rms(Y_uc4{dq});

rmsu2xc3{dq}=rms(X_u2c3{dq});
rmsu2xcz{dq}=rms(X_u2cz{dq});
rmsu2xc4{dq}=rms(X_u2c4{dq});

rmsu2yc3{dq}=rms(Y_u2c3{dq});
rmsu2ycz{dq}=rms(Y_u2cz{dq});
rmsu2yc4{dq}=rms(Y_u2c4{dq});


rmsbxc3{dq}=rms(X_bc3{dq});
rmsbxcz{dq}=rms(X_bcz{dq});
rmsbxc4{dq}=rms(X_bc4{dq});

rmsbyc3{dq}=rms(Y_bc3{dq});
rmsbycz{dq}=rms(Y_bcz{dq});
rmsbyc4{dq}=rms(Y_bc4{dq});

% sparsex=sparsefilt(X{ds},10,'IterationLimit',20);
% sparsey=sparsefilt(Y{ds},10,'IterationLimit',20);
end
%% feature vector
% first Mean

meanuxxc3=cell2mat(meanuxc3);
meanuxxcz=cell2mat(meanuxcz);
meanuxxc4=cell2mat(meanuxc4);

meanuyyc3=cell2mat(meanuyc3);
meanuyycz=cell2mat(meanuycz);
meanuyyc4=cell2mat(meanuyc4);

meanu2xxc3=cell2mat(meanu2xc3);
meanu2xxcz=cell2mat(meanu2xcz);
meanu2xxc4=cell2mat(meanu2xc4);

meanu2yyc3=cell2mat(meanu2yc3);
meanu2yycz=cell2mat(meanu2ycz);
meanu2yyc4=cell2mat(meanu2yc4);

meanbxxc3=cell2mat(meanbxc3);
meanbxxcz=cell2mat(meanbxcz);
meanbxxc4=cell2mat(meanbxc4);

meanbyyc3=cell2mat(meanbyc3);
meanbyycz=cell2mat(meanbycz);
meanbyyc4=cell2mat(meanbyc4);

stduxxc3=cell2mat(stduxc3);
stduxxcz=cell2mat(stduxcz);
stduxxc4=cell2mat(stduxc4);

stduyyc3=cell2mat(stduyc3);
stduyycz=cell2mat(stduycz);
stduyyc4=cell2mat(stduyc4);

stdu2xxc3=cell2mat(stdu2xc3);
stdu2xxcz=cell2mat(stdu2xcz);
stdu2xxc4=cell2mat(stdu2xc4);

stdu2yyc3=cell2mat(stdu2yc3);
stdu2yycz=cell2mat(stdu2ycz);
stdu2yyc4=cell2mat(stdu2yc4);

stdbxxc3=cell2mat(stdbxc3);
stdbxxcz=cell2mat(stdbxcz);
stdbxxc4=cell2mat(stdbxc4);

stdbyyc3=cell2mat(stdbyc3);
stdbyycz=cell2mat(stdbycz);
stdbyyc4=cell2mat(stdbyc4);

varuxxc3=cell2mat(varuxc3);
varuxxcz=cell2mat(varuxcz);
varuxxc4=cell2mat(varuxc4);

varuyyc3=cell2mat(varuyc3);
varuyycz=cell2mat(varuycz);
varuyyc4=cell2mat(varuyc4);

varu2xxc3=cell2mat(varu2xc3);
varu2xxcz=cell2mat(varu2xcz);
varu2xxc4=cell2mat(varu2xc4);

varu2yyc3=cell2mat(varu2yc3);
varu2yycz=cell2mat(varu2ycz);
varu2yyc4=cell2mat(varu2yc4);

varbxxc3=cell2mat(varbxc3);
varbxxcz=cell2mat(varbxcz);
varbxxc4=cell2mat(varbxc4);

varbyyc3=cell2mat(varbyc3);
varbyycz=cell2mat(varbycz);
varbyyc4=cell2mat(varbyc4);

psduxxc3=cell2mat(psduxc3);
psduxxcz=cell2mat(psduxcz);
psduxxc4=cell2mat(psduxc4);

psduyyc3=cell2mat(psduyc3);
psduyycz=cell2mat(psduycz);
psduyyc4=cell2mat(psduyc4);

psdu2xxc3=cell2mat(psdu2xc3);
psdu2xxcz=cell2mat(psdu2xcz);
psdu2xxc4=cell2mat(psdu2xc4);

psdu2yyc3=cell2mat(psdu2yc3);
psdu2yycz=cell2mat(psdu2ycz);
psdu2yyc4=cell2mat(psdu2yc4);

psdbxxc3=cell2mat(psdbxc3);
psdbxxcz=cell2mat(psdbxcz);
psdbxxc4=cell2mat(psdbxc4);

psdbyyc3=cell2mat(psdbyc3);
psdbyycz=cell2mat(psdbycz);
psdbyyc4=cell2mat(psdbyc4);

sampuxxc3=cell2mat(sampuxc3);
sampuxxcz=cell2mat(sampuxcz);
sampuxxc4=cell2mat(sampuxc4);

sampuyyc3=cell2mat(sampuyc3);
sampuyycz=cell2mat(sampuycz);
sampuyyc4=cell2mat(sampuyc4);

sampu2xxc3=cell2mat(sampu2xc3);
sampu2xxcz=cell2mat(sampu2xcz);
sampu2xxc4=cell2mat(sampu2xc4);

sampu2yyc3=cell2mat(sampu2yc3);
sampu2yycz=cell2mat(sampu2ycz);
sampu2yyc4=cell2mat(sampu2yc4);

sampbxxc3=cell2mat(sampbxc3);
sampbxxcz=cell2mat(sampbxcz);
sampbxxc4=cell2mat(sampbxc4);

sampbyyc3=cell2mat(sampbyc3);
sampbyycz=cell2mat(sampbycz);
sampbyyc4=cell2mat(sampbyc4);

plvuxx1=cell2mat(plvux1);
plvuxx2=cell2mat(plvux2);

plvuyy1=cell2mat(plvuy1);
plvuyy2=cell2mat(plvuy2);

plvu2xx1=cell2mat(plvu2x1);
plvu2xx2=cell2mat(plvu2x2);

plvu2yy1=cell2mat(plvu2y1);
plvu2yy2=cell2mat(plvu2y2);

plvbxx1=cell2mat(plvbx1);
plvbxx2=cell2mat(plvbx2);

plvbyy1=cell2mat(plvby1);
plvbyy2=cell2mat(plvby2);

rmsuxxc3=cell2mat(rmsuxc3);
rmsuxxcz=cell2mat(rmsuxcz);
rmsuxxc4=cell2mat(rmsuxc4);

rmsuyyc3=cell2mat(rmsuyc3);
rmsuyycz=cell2mat(rmsuycz);
rmsuyyc4=cell2mat(rmsuyc4);

rmsu2xxc3=cell2mat(rmsu2xc3);
rmsu2xxcz=cell2mat(rmsu2xcz);
rmsu2xxc4=cell2mat(rmsu2xc4);

rmsu2yyc3=cell2mat(rmsu2yc3);
rmsu2yycz=cell2mat(rmsu2ycz);
rmsu2yyc4=cell2mat(rmsu2yc4);

rmsbxxc3=cell2mat(rmsbxc3);
rmsbxxcz=cell2mat(rmsbxcz);
rmsbxxc4=cell2mat(rmsbxc4);

rmsbyyc3=cell2mat(rmsbyc3);
rmsbyycz=cell2mat(rmsbycz);
rmsbyyc4=cell2mat(rmsbyc4);

FVX=[meanuxxc3' meanuxxc4' meanu2xxc3' meanu2xxc4' meanbxxc3' meanbxxc4' stduxxc3' stduxxc4' stdu2xxc3' stdu2xxc4' stdbxxc3' stdbxxc4' varuxxc3' varuxxc4' varu2xxc3' varu2xxc4' varbxxc3' varbxxc4' psduxxc3' psduxxc4' psdu2xxc3' psdu2xxc4' psdbxxc3' psdbxxc4' sampuxxc3' sampuxxc4' sampu2xxc3' sampu2xxc4' sampbxxc3' sampbxxc4' plvuxx1' plvuxx2' plvu2xx1' plvu2xx2' plvbxx1' plvbxx2' rmsuxxc3' rmsuxxc4' rmsu2xxc3' rmsu2xxc4' rmsbxxc3' rmsbxxc4'];
FVY=[meanuyyc3' meanuyyc4' meanu2yyc3' meanu2yyc4' meanbyyc3' meanbyyc4' stduyyc3' stduyyc4' stdu2yyc3' stdu2yyc4' stdbyyc3' stdbyyc4' varuyyc3' varuyyc4' varu2yyc3' varu2yyc4' varbyyc3' varbyyc4' psduyyc3' psduyyc4' psdu2yyc3' psdu2yyc3' psdbyyc3' psdbyyc4' sampuyyc3' sampuyyc4' sampu2yyc3' sampu2yyc4' sampbyyc3' sampbyyc4' plvuyy1' plvuyy2' plvu2yy1' plvu2yy2' plvbyy1' plvbyy2' rmsuyyc3' rmsuyyc4' rmsu2yyc3' rmsu2yyc4' rmsbyyc3' rmsbyyc4'];
 
label = [repmat('0',N_trails,1) ; repmat('1',N_trails,1)];
FV=[FVX;FVY];FV1=table(FV,label);
Label_fs=[ones(N_trails,1); -ones(N_trails,1)];

%%fisher 
%[index,featureScore] = feature_rank(FV,label)

% First using NCA
% lbfgs   fsrnca  sgd
tic
mdl = fscnca(FV,Label_fs,'FitMethod','exact','Verbose',1,...
              'Solver','sgd');
figure(7)
plot(mdl.FeatureWeights,'ro')
grid on
xlabel('Feature index')
ylabel('Feature weight')
title('NCA')

cvp           = cvpartition(Label_fs,'kfold',5);
numtestsets   = cvp.NumTestSets;
lambdavalues  = linspace(0,2,100)/length(Label_fs);
lossvalues    = zeros(length(lambdavalues),numtestsets);

for i = 1:length(lambdavalues)
    for k = 1:numtestsets

        % Extract the training set from the partition object
        Xtrain = FV(cvp.training(k),:);
        ytrain = Label_fs(cvp.training(k),:);

        % Extract the test set from the partition object
        Xtest  = FV(cvp.test(k),:);
        ytest  = Label_fs(cvp.test(k),:);

        % Train an nca model for classification using the training set
%          ncaMdl = fitcsvm(Xtrain,ytrain);
        
        ncaMdl = fscnca(Xtrain,ytrain,'FitMethod','exact',...
            'Solver','sgd','Lambda',lambdavalues(i));

        % Compute the classification loss for the test set using the nca
        % model
        lossvalues(i,k) = loss(ncaMdl,Xtest,ytest,...
            'LossFunction','quadratic');

    end
end

figure(100)
plot(lambdavalues,mean(lossvalues,2),'ro-');
xlabel('Lambda values');
ylabel('Loss values');
grid on;


[~,idx] = min(mean(lossvalues,2)); % Find the index
bestlambda = lambdavalues(idx) % Find the best lambda value

% A2 = unique(A(:));
% out = A2(2); 

ncaMdl = fscnca(FV,Label_fs,'FitMethod','exact','Verbose',1,...
     'Solver','sgd','Lambda',bestlambda);
figure(102)
plot(ncaMdl.FeatureWeights,'ro');
xlabel('Feature index');
ylabel('Feature weight');
grid on;
 toc
% Second using RNCA

mdl_rnca = fscnca(FV,Label_fs,'Solver','lbfgs','Verbose',1);
figure(8)
plot(mdl_rnca.FeatureWeights,'ro')
grid on
xlabel('Feature index')
ylabel('Feature weight')
title('RNCA')

% Third Relieff
tic 
[ranks,mdl_relief] = relieff(FV,Label_fs,10);
relieff_rank= ranks(1:42)
figure(9)
bar(mdl_relief(ranks))
xlabel('Predictor rank')
ylabel('predictor weight')
title('Relieff')
toc
% fourth PCA 
tic
PCAxy= pca(FV);
pcavar= var(PCAxy);
pcmean= mean(PCAxy);
% pca_feat=pcmean-pcavar;
% pca_feat = max(pca_feat,0);
figure(10)
hold on
plot(1:42,pcavar);
plot(1:42,pcmean);
hold off
title('PCA ')
toc

% save data for ga algorithm
Data.X=FV;
Data.Y=cell(N_trails,1)
for ga=1:N_trails
Data.Y{ga}='p';
end
for ga1=N_trails+1:2*N_trails
    Data.Y{ga1}='n';
end
save gaB0102T Data
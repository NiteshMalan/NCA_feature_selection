function out = dtcwtb(X,L)

[Faf, Fsf] = FSfarras; % 1st stage anal. & synth. filters
[af, sf] = dualfilt1;
dtcwtX = dddtree('cplxdt',X,L,Faf,af);

dtcwtX.cfs{1,1}(:,:,1)=0;
dtcwtX.cfs{1,1}(:,:,2)=0;
% dtcwtX.cfs{1,2}(:,:,1)=0;
% dtcwtX.cfs{1,2}(:,:,2)=0;
dtcwtX.cfs{1,3}(:,:,1)=0;
dtcwtX.cfs{1,3}(:,:,2)=0;
dtcwtX.cfs{1,4}(:,:,1)=0;
dtcwtX.cfs{1,4}(:,:,2)=0;
dtcwtX.cfs{1,5}(:,:,1)=0;
dtcwtX.cfs{1,5}(:,:,2)=0;

out=(dddtreecfs('r',dtcwtX));
end
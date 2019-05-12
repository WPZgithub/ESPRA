% For GreeneData
clc

% Input data
% Select from£ºbirthdeath, expandcontract, hide, mergesplit
nets = importdata('./Dataset/GreeneData/hide/networklist.txt');
T = length(nets);
N = 1000;
adjMatrix=cell(1,T);
for i = 1:T
    % Select from£ºbirthdeath, expandcontract, hide, mergesplit
    adjMatrix{i}=edgelistFile2adjMatrix(['./Dataset/GreeneData/hide/',nets{i}], N);
end

% evolutionary clustering, alpha=0.8
[result] = ESPRA( adjMatrix, 0.8, 0.5);

% save to files
% for t=1:T
%     fp = fopen(['./GreeneData/hide/result_', num2str(t), '.txt'],'wt');
%     r = result{t};
%     for i=1:size(r, 1)
%         fprintf(fp, '%d\t%d\n', r(i, 1), r(i, 2));
%     end
% end
% fclose(fp); 

% performance
nmi=zeros(T, 1);
errate=zeros(T, 1);
for i = 1:T
    % Ground truth. Select from: birthdeath, expandcontract, hide, mergesplit
    benchmark_cluster=load(['./Dataset/GreeneData/hide/t',num2str(i),'.txt']);
    benchmark_cluster=benchmark_cluster(benchmark_cluster(:,1)~=0,:);

    nmi(i) = NMI( result{i}, benchmark_cluster );
    errate(i) = ErrorRate( result{i}, benchmark_cluster );
end
disp('---------------------------')

% figure of performance
figure();
%subplot(1,2,1);
plot(nmi, 'sr-',  'LineWidth',2, 'MarkerFaceColor','r');
grid on
%axis([0 T+1 0 1.05]);
xlabel('Time steps');
ylabel('NMI');
title('NMI over Time steps');
% subplot(1,2,2);
% plot(errate, '*b-');
% grid on
% axis([0 T+1 0 max(errate)+0.05]);
% xlabel('Time steps');
% ylabel('Error rate');
% title('Error rate over Time steps');


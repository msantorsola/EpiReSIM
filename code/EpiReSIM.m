function EpiReSIM(Case_Num,Control_Num,SNP_num,MAF,PD,H2,n,repeat,SNP_name,data_mat,data_txt,realdata_name)
% SNP_num:the number of  SNPs
% Case_Num:the the number of case samples
% Control_Num:the the number of control samples
% MAF:the 
% PD:the value of prevalence of model
% H2:the value of heritability of model.If user calculates model with heritability,
%then user should give the value of heritability,otherwise the values
%shuould be 0.
% n: the order of the model
% repeat:the number of simulated dataset;
% EpiReSIM(100,100,50,[0.2,0.3,0.1],0.2,0.3,3,10,'SNP1',1,1,'T1D.mat')
a = Data(realdata_name,SNP_num);
pt.order= n;
%% Select the control dataset and calculate the MAF of them
idx = find(all(a.class(1,:)==2,1));
% idx = find(all(a.class(1,:)==1,1));
d = a.pts(idx,:);
maf =zeros(1,size(d,2));
for j = 1:size(d,2)
    Aa = 0;
    aa = 0;
    for i = 1:size(d,1)
        if d(i,j)==2
            Aa= Aa+1;
        end
        if d(i,j)==3
            aa=aa+1;
        end
        i = i+1;
    end
    maf(j)= (Aa+aa*2)/(size(d,1)*2);
    j = j+1;
end
%% Select the site of the model
for j = 1:n
    s = 0.01;%Search in steps of 0.01
    candidates = [];
    while 1
        for i = 1:size(d,2)
            if maf(i) >= (MAF(j)-s) && maf(i) <= (MAF(j)+s)
               candidates= [candidates,i];
            end
        end
        candidates_size = size(candidates);
        if candidates_size ~= [0 0] 
            break;
        else s = s + 0.01;
        end
    end
    site(j) = candidates(randperm(numel(candidates),1));
end
for i = 1:n
    pt.MAF(i) = maf(site(i));
end
pt.loci =site(1,:);
%% Calculate model using prevalence or heritability
% p is prevalence
% h is heritability
% x is baseline penetrance
% n is the order of the model
[pt.penetrance,pt.pd,pt.h] = calculation(pt.MAF,PD,H2,n);
%% Generate simulated data and calculate sample labels
for k = 1:repeat
    SNP = simulation(d,Control_Num,Case_Num,pt,SNP_num);
    if data_mat
        filename1 = strcat(SNP_name,'_',num2str(k),'.mat');
        save(filename1,'SNP');
    end
    if data_txt
        filename2 = strcat(SNP_name,'_',num2str(k),'.txt');
        [Row,Col]=size(SNP);
        fid=fopen(filename2,'a');
        for i=1:Row
            for j=1:Col
                fprintf(fid,'%s\t',num2str(SNP(i,j)));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
    end
end
%% 保存模型信息
fid=fopen('log.txt','w');
fprintf(fid,'%s', 'order:');
fprintf(fid,'%d\n', pt.order);
fprintf(fid,'\n');
fprintf(fid,'%s', 'MAF:');
for i=1:n
    fprintf(fid,'%f\t', pt.MAF(i));
end
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%s', 'loci:');
for i=1:n
    fprintf(fid,'%d\t', pt.loci(i));
end
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%s\n', 'Penetrance:');
for i=1:3^n
    fprintf(fid,'%f\t', pt.penetrance(i));
    if mod(i,3)==0
        fprintf(fid,'\n');
    end
    if mod(i,9)==0
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');
fprintf(fid,'%s', 'Prevalence:');
fprintf(fid,'%f\n', pt.pd);
fprintf(fid,'\n');
fprintf(fid,'%s', 'Heritability:');
fprintf(fid,'%f\n', pt.h);
fclose(fid);
end


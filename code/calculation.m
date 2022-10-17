function [x_result,pd,h]=calculation(MAF,PD,H2,n)
% n=n; the order of the model
k= n-1;
syms x [3^n 1]
X= x;
% PD=p;
% H2=h;
D=double(genotype_probabilities(MAF)');
%% the calculation of 2-oder model
if n==2
    %  the marginal penetrance values of A
    PD_AA=0;
    PD_Aa=0;
    PD_aa=0;
    PD_BB=0;
    PD_Bb=0;
    PD_bb=0;
    for i=1:3^k
        PD_AA=PD_AA+X(i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_Aa=PD_Aa+X(3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_aa=PD_aa+X(2*3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
    end
    PD_AA=PD_AA - PD;
    PD_Aa=PD_Aa - PD;
    PD_aa=PD_aa - PD;
    %  the marginal penetrance values of A
    for j=1:3^k
        i=3*(j-1)+1;
        PD_BB=PD_BB + X(i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        PD_Bb=PD_Bb + X(1+i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        PD_bb=PD_bb + X(2+i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
    end
    PD_BB=PD_BB - PD;
    PD_Bb=PD_Bb - PD;
    PD_bb=PD_bb - PD;
    f10= D*X-PD;
    if H2==0 % using the prevalence to calculate the eNME model
        [E,f] = equationsToMatrix([PD_AA;PD_Aa;PD_aa;PD_BB;PD_Bb;PD_bb;f10],[X]);
        E = double(E);
        f = double(f);
        x_result = lsqminnorm(E,f);
    else % When using both prevalence and heritability to calculate the eNME model
        f11=0;
        for i =1:3^n
            f11 = f11+D(i)*(( X(i)- PD)^2);
        end
        f11=f11-H2*PD*(1-PD);
        f = [PD_AA;PD_Aa;PD_aa;PD_BB;PD_Bb;PD_bb;f10;f11];
        % zero vector
        x0 = zeros(3^n,1);
        error_dxk = 0.01;
        error_fkk = 0.05;
        num = 10;
        % jacobi1 = [diff(f1,x1) diff(f1,x2);diff(f2,x1) diff(f2,x2)]
        % 直接用自带函数求雅克比矩阵:
        %jacobi = jacobian([f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11],[x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16; x17; x18; x19; x20; x21; x22; x23; x24; x25; x26; x27]);
        jacobi = jacobian(f,X);
        for k = 1:num
            Ak = double( subs(jacobi, X, x0) );
            bk = double( subs(f, X, x0) );
            dxk = pre_seidel(Ak,-bk,k);  % 
            x0 = x0 + dxk;
            fkk = double( subs(f, X, x0) );  % fk+1
            if norm(dxk) < error_dxk | norm(fkk) < error_fkk
                break;
            end
        end
        
        if k < num
            x_result = x0;
            x_result;
        else
            x_result = x0;
            x_result;
        end
    end
    for i=1:3^n
        if x_result(i)<0
            x_result(i)=0;
        end
    end
end
%% the calculation of 3-oder model
if n==3
    % 计算边际外显%% A位点的等位基因的边际外显率
    PD_AA=0;
    PD_Aa=0;
    PD_aa=0;
    i =0;
    for i=1:3^k
        PD_AA=PD_AA+X(i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        %     PD_Aa=PD_Aa+RTable(3^k+i,1)*(Combin_Freq(i,1)+Combin_Freq(3^k+i,1)+Combin_Freq(2*3^k+i,1));
        %     PD_Aa=PD_Aa+RTable(3^k+i,1)*(Combin_Freq(1,i)+Combin_Freq(1,3^k+i)+Combin_Freq(1,2*3^k+i));
        %     PD_aa=PD_aa+RTable(2*3^k+i,1)*(Combin_Freq(1,i)+Combin_Freq(1,2*3^k+i)+Combin_Freq(1,2*3^k+i));
        PD_Aa=PD_Aa+X(3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_aa=PD_aa+X(2*3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
    end
    PD_AA=PD_AA - PD;
    PD_Aa=PD_Aa - PD;
    PD_aa=PD_aa - PD;
    % B位点的等位基因的边际外显率
    % 
    i=0;
    j=0;
    PD_BB=0;
    PD_Bb=0;
    PD_bb=0;
    for o=0:2
        % for o=0:k
        for j = 1:3^(k-1)
            %     for j = 1:3
            %         i =j+ 9*o;
            i =j+ o*(3^k);
            %         PD_BB=PD_BB + RTable(i,1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            %         PD_Bb=PD_Bb + RTable(i+3^(k-1),1)*(Combin_Freq(1,i)+Combin_Freq(1,3^(k-1)+i)+Combin_Freq(1,2*3^(k-1)+i));
            %         PD_bb=PD_bb + RTable(i+2*3^(k-1),1)*(Combin_Freq(1,i)+Combin_Freq(1,3^(k-1)+i)+Combin_Freq(1,2*3^(k-1)+i));
            PD_BB=PD_BB + X(i,1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            PD_Bb=PD_Bb + X(i+3^(k-1),1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            PD_bb=PD_bb + X(i+2*3^(k-1),1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            
        end
    end
    PD_BB=PD_BB - PD;
    PD_Bb=PD_Bb - PD;
    PD_bb=PD_bb - PD;
    % C位点的等位基因频率
    i=0;
    j=0;
    PD_CC=0;
    PD_Cc=0;
    PD_cc=0;
    for o=0:2
        for j=0:k
            i = 1+o*3^k+3*j;
            %         PD_CC=PD_CC + RTable(i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
            %         PD_Cc=PD_Cc + RTable(i+1,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
            %         PD_cc=PD_cc + RTable(i+2,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
            PD_CC=PD_CC + X(i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
            PD_Cc=PD_Cc + X(i+1,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
            PD_cc=PD_cc + X(i+2,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        end
    end
    PD_CC=PD_CC - PD;
    PD_Cc=PD_Cc - PD;
    PD_cc=PD_cc - PD;
    f10= D*X-PD;
    if H2==0
        [E,f] = equationsToMatrix([PD_AA;PD_Aa;PD_aa;PD_BB;PD_Bb;PD_bb;PD_CC;PD_Cc;PD_cc;f10;],[X]);
        E = double(E);
        f = double(f);
        x_result = lsqminnorm(E,f);
    else
        f11=0;
        for i =1:3^n
            f11 = f11+( X(i)- PD)^2* D(i);
        end
        f11=f11-H2*PD*(1-PD);
        f = [PD_AA;PD_Aa;PD_aa;PD_BB;PD_Bb;PD_bb;PD_CC;PD_Cc;PD_cc;f10;f11];
        % 初始值: 统一用列向量
        x0 = zeros(3^n,1);
        error_dxk = 0.05;
        error_fkk = 0.05;
        num = 10;
        % jacobi1 = [diff(f1,x1) diff(f1,x2);diff(f2,x1) diff(f2,x2)]
        % 
        %jacobi = jacobian([f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11],[x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16; x17; x18; x19; x20; x21; x22; x23; x24; x25; x26; x27]);
        jacobi = jacobian(f,X);
        for k = 1:num
            Ak = double( subs(jacobi, X, x0) );
            bk = double( subs(f, X, x0) );
            dxk = pre_seidel(Ak,-bk,k);  %
            x0 = x0 + dxk;
            fkk = double( subs(f, X, x0) );  % fk+1
            if norm(dxk) < error_dxk | norm(fkk) < error_fkk
                break;
            end
        end
        
        if k < num
            x_result = x0;
            x_result;
        else
            x_result = x0;
            x_result;
        end
    end
    for i=1:3^n
        if x_result(i)<0
            x_result(i)=0;
        end
    end
end
%% the calculation of 4-oder model
if n==4
    % A的边际外显率的计算公式
    i =0;
    j=0;
    PD_AA=0;
    PD_Aa=0;
    PD_aa=0;
    for i=1:3^k
        PD_AA=PD_AA+X(i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_Aa=PD_Aa+X(3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_aa=PD_aa+X(2*3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
    end
    PD_AA=PD_AA - PD;
    PD_Aa=PD_Aa - PD;
    PD_aa=PD_aa - PD;
    % B位点的等位基因的边际外显率
    i=0;
    j=0;
    PD_BB=0;
    PD_Bb=0;
    PD_bb=0;
    for o=0:2
        for j = 1:3^(k-1)
            %         i =j+ 9*o;
            i =j+ o*(3^k);
            PD_BB=PD_BB + X(i,1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            PD_Bb=PD_Bb + X(i+3^(k-1),1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            PD_bb=PD_bb + X(i+2*3^(k-1),1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            
        end
    end
    PD_BB=PD_BB - PD;
    PD_Bb=PD_Bb - PD;
    PD_bb=PD_bb - PD;
    % C位点的等位基因的边际外显率
    i=0;
    j=0;
    PD_CC=0;
    PD_Cc=0;
    PD_cc=0;
    for o=0:3^(k-1)-1
        for j=1:3
            i = j+(3^(n-2))*o;
            PD_CC=PD_CC + X(i,1)*(D(1,i)+D(1,3^1+i)+D(1,2*3^1+i));
            PD_Cc=PD_Cc + X(i+3*1,1)*(D(1,i)+D(1,3^1+i)+D(1,2*3^1+i));
            PD_cc=PD_cc + X(i+3*2,1)*(D(1,i)+D(1,3^1+i)+D(1,2*3^1+i));
        end
    end
    PD_CC=PD_CC - PD;
    PD_Cc=PD_Cc - PD;
    PD_cc=PD_cc - PD;
    % D位点等位基因的边际外显率
    i=0;
    j=0;
    PD_DD=0;
    PD_Dd=0;
    PD_dd=0;
    for o=0:2
        for j=0:k
            i = 1+o*3^k+3*j;
            PD_DD=PD_DD + X(i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
            PD_Dd=PD_Dd + X(i+1,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
            PD_dd=PD_dd + X(i+2,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        end
    end
    PD_DD=PD_DD - PD;
    PD_Dd=PD_Dd - PD;
    PD_dd=PD_dd - PD;
    f10= D*X-PD;
    if H2==0
        % 不使用遗传力计算 转换成系数矩阵
        [E,f] = equationsToMatrix([PD_AA;PD_Aa;PD_aa;PD_BB;PD_Bb;PD_bb;PD_CC;PD_Cc;PD_cc;PD_DD;PD_Dd;PD_dd;f10;],X);
        E = double(E);
        f = double(f);
        x_result = lsqminnorm(E,f);
    else
        % 使用患病率计算 牛顿法
        f11=0;
        for i =1:3^n
            f11 = f11+( X(i)- PD)^2* D(i);
        end
        f11=f11-H2*PD*(1-PD);
        f = [PD_AA;PD_Aa;PD_aa;PD_BB;PD_Bb;PD_bb;PD_CC;PD_Cc;PD_cc;PD_DD;PD_Dd;PD_dd;f10;f11];
        % 初始值: 统一用列向量
        x0 = zeros(3^n,1);
        error_dxk = 0.03;
        error_fkk = 0.03;
        num = 10;
        % jacobi1 = [diff(f1,x1) diff(f1,x2);diff(f2,x1) diff(f2,x2)]
        % 直接用自带函数求雅克比矩阵:
        %jacobi = jacobian([f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11],[x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16; x17; x18; x19; x20; x21; x22; x23; x24; x25; x26; x27]);
        jacobi = jacobian(f,X);
        for k = 1:num
            Ak = double( subs(jacobi, X, x0) );
            bk = double( subs(f, X, x0) );
            dxk = pre_seidel(Ak,-bk,k);  % 步长
            x0 = x0 + dxk;
            fkk = double( subs(f, X, x0) );  % fk+1单纯用来判断
            if norm(dxk) < error_dxk | norm(fkk) < error_fkk
                break;
            end
        end
        if k < num
            x_result = x0;
            x_result;
        else
            x_result = x0;
            x_result;
        end
    end
    for i=1:3^n
        if x_result(i)<0
            x_result(i)=0;
        end
    end
end
%% the calculation of 5-oder model
if n==5
    % A的边际外显率的计算公式
    i =0;
    j=0;
    PD_AA=0;
    PD_Aa=0;
    PD_aa=0;
    for i=1:3^k
        PD_AA=PD_AA+X(i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_Aa=PD_Aa+X(3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_aa=PD_aa+X(2*3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
    end
    PD_AA=PD_AA - PD;
    PD_Aa=PD_Aa - PD;
    PD_aa=PD_aa - PD;
    % B位点的等位基因的边际外显率
    i=0;
    j=0;
    PD_BB=0;
    PD_Bb=0;
    PD_bb=0;
    for o=0:2
        for j = 1:3^(k-1)
            %         i =j+ 9*o;
            i =j+ o*(3^k);
            PD_BB=PD_BB + X(i,1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            PD_Bb=PD_Bb + X(i+3^(k-1),1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            PD_bb=PD_bb + X(i+2*3^(k-1),1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            
        end
    end
    PD_BB=PD_BB - PD;
    PD_Bb=PD_Bb - PD;
    PD_bb=PD_bb - PD;
    % C位点的等位基因的边际外显率
    i=0;
    j=0;
    PD_CC=0;
    PD_Cc=0;
    PD_cc=0;
    for o=0:3^(n-3)-1%3^(k-2)-1
        for j=1:3^2
            i = j+(3^(n-2))*o;
            PD_CC=PD_CC + X(i,1)*(D(1,i)+D(1,3^(k-2)+i)+D(1,2*3^(k-2)+i));
            PD_Cc=PD_Cc + X(i+3^(k-2),1)*(D(1,i)+D(1,3^(k-2)+i)+D(1,2*3^(k-2)+i));
            PD_cc=PD_cc + X(i+2*3^(k-2),1)*(D(1,i)+D(1,3^(k-2)+i)+D(1,2*3^(k-2)+i));
        end
    end
    PD_CC=PD_CC - PD;
    PD_Cc=PD_Cc - PD;
    PD_cc=PD_cc - PD;
    % D位点等位基因的边际外显率
    i=0;
    j=0;
    PD_DD=0;
    PD_Dd=0;
    PD_dd=0;
    for o=0:3^(k-1)-1
        for j=1:3
            i = j+(3^(n-3))*o;
            PD_DD=PD_DD +X(i,1)*(D(1,i)+D(1,3^(k-3)+i)+D(1,2*3^(k-3)+i));
            PD_Dd=PD_Dd +X(i+3^(k-3),1)*(D(1,i)+D(1,3^(k-3)+i)+D(1,2*3^(k-3)+i));
            PD_dd=PD_dd +X(i+2*3^(k-3),1)*(D(1,i)+D(1,3^(k-3)+i)+D(1,2*3^(k-3)+i));
        end
    end
    PD_DD=PD_DD - PD;
    PD_Dd=PD_Dd - PD;
    PD_dd=PD_dd - PD;
    % E位点等位基因的边际外显率
    i=0;
    j=0;
    PD_EE=0;
    PD_Ee=0;
    PD_ee=0;
    for o=0:3^k-1
        %     for j=1:1
        j=1;
        i = j+(3^(n-4))*o;
        PD_EE=PD_EE + X(i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        PD_Ee=PD_Ee + X(i+1,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        PD_ee=PD_ee + X(i+2,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        %     end
    end
    PD_EE=PD_EE - PD;
    PD_Ee=PD_Ee - PD;
    PD_ee=PD_ee - PD;
    f10= D*X-PD;
    if H2==0
        % 不使用遗传力计算 转换成系数矩阵
        [E,f] = equationsToMatrix([PD_AA;PD_Aa;PD_aa;PD_BB;PD_Bb;PD_bb;PD_CC;PD_Cc;PD_cc;PD_DD;PD_Dd;PD_dd;PD_EE;PD_Ee;PD_ee;f10;],X);
        E = double(E);
        f = double(f);
        x_result = lsqminnorm(E,f);
    else
        % 使用患病率计算 牛顿法
        f11=0;
        for i =1:3^n
            f11 = f11+( X(i)- PD)^2* D(i);
        end
        f11=f11-H2*PD*(1-PD);
        f = [PD_AA;PD_Aa;PD_aa;PD_BB;PD_Bb;PD_bb;PD_CC;PD_Cc;PD_cc;PD_DD;PD_Dd;PD_dd;PD_EE;PD_Ee;PD_ee;f10;f11];
        % 初始值: 统一用列向量
        x0 = zeros(3^n,1);
        error_dxk = 0.05;
        error_fkk = 0.05;
        num = 10;
        % jacobi1 = [diff(f1,x1) diff(f1,x2);diff(f2,x1) diff(f2,x2)]
        % 直接用自带函数求雅克比矩阵:
        %jacobi = jacobian([f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11],[x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16; x17; x18; x19; x20; x21; x22; x23; x24; x25; x26; x27]);
        jacobi = jacobian(f,X);
        for k = 1:num
            Ak = double( subs(jacobi, X, x0) );
            bk = double( subs(f, X, x0) );
            dxk = pre_seidel(Ak,-bk,k);  % 步长
            x0 = x0 + dxk;
            fkk = double( subs(f, X, x0) );  % fk+1单纯用来判断
            if norm(dxk) < error_dxk | norm(fkk) < error_fkk
                break;
            end
        end
        if k < num
            x_result = x0;
            x_result;
        else
            x_result = x0;
            x_result;
        end
    end
    for i=1:3^n
        if x_result(i)<0
            x_result(i)=0;
        end
    end
end

[v,pd,h] = verification(x_result,n,D);
if v==0
    x_result =zeros(3^n,1);
    ME = MException("MYFUN:Badcalculation", "There is no solution to the problem defined.");
    throwAsCaller(ME);
%     fprintf( "There is no solution to the problem defined.\n");
%     x_result =zeros(3^n,1);
end
end

function [x] = pre_seidel(A,b,n)
%  this fuction is cited from:https://github.com/GaoBoYu599/Num_Func/tree/master/Nonlinear_Equations
% 预处理: 就是这么简单
b = A'*b;
A = A'*A;
% 下面是正常的高斯-赛德尔操作:
D = diag(diag(A));
L = tril(A,-1);    
U = triu(A,1);      
B2 = -inv(D+L)*U;   
g2 = inv(D+L)*b;

radius = max(abs(eig(B2)));  
x = zeros(length(b),1);  
error = 0.0001;
count = 0;      
while 1
    tmp = B2*x + g2;
    if max(abs(tmp - x)) < error
        break;
    end
    x = tmp;
    count = count + 1;
end
end

function [v,pd,h] = verification(x_result,n,D)
PD_=zeros(n,3);
if n==2
    %A的边际外显率的计算公式
    k= n-1;;%order-1
    i =0;
    j=0;
    for i=1:3^k
        PD_(1,1)=PD_(1,1)+x_result(i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_(1,2)=PD_(1,2)+x_result(3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_(1,3)=PD_(1,3)+x_result(2*3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
    end
    %B位点的等位基因的边际外显率
    i=0;
    j=0;
    for j=1:3^k
        i=3*(j-1)+1;
        PD_(2,1)=PD_(2,1) + x_result(i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        PD_(2,2)=PD_(2,2) + x_result(1+i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        PD_(2,3)=PD_(2,3) + x_result(2+i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
    end
end
if n==3
    %A的边际外显率的计算公式
    k= n-1;;%order-1
    i =0;
    j=0;
    for i=1:3^k
        PD_(1,1)=PD_(1,1)+x_result(i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_(1,2)=PD_(1,2)+x_result(3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_(1,3)=PD_(1,3)+x_result(2*3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
    end
    %B位点的等位基因的边际外显率
    i=0;
    j=0;
    for o=0:2
        % for o=0:k
        for j = 1:3^(k-1)
            %     for j = 1:3
            %         i =j+ 9*o;
            i =j+ o*(3^k);
            PD_(2,1)=PD_(2,1) + x_result(i,1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            PD_(2,2)=PD_(2,2) + x_result(i+3^(k-1),1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            PD_(2,3)=PD_(2,3) + x_result(i+2*3^(k-1),1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            
        end
    end
    %C位点的等位基因的边际外显率
    i=0;
    j=0;
    for o=0:2
        for j=0:k
            i = 1+o*3^k+3*j;
            PD_(3,1)=PD_(3,1) + x_result(i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
            PD_(3,2)=PD_(3,2) + x_result(i+1,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
            PD_(3,3)=PD_(3,3) + x_result(i+2,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        end
    end
end
if n==4
    %A的边际外显率的计算公式
    k=3;%order-1
    i =0;
    j=0;
    for i=1:3^k
        PD_(1,1)=PD_(1,1)+x_result(i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_(1,2)=PD_(1,2)+x_result(3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_(1,3)=PD_(1,3)+x_result(2*3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
    end
    %B位点的等位基因的边际外显率
    i=0;
    j=0;
    for o=0:2
        for j = 1:3^(k-1)
            %         i =j+ 9*o;
            i =j+ o*(3^k);
            PD_(2,1)=PD_(2,1) + x_result(i,1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            PD_(2,2)=PD_(2,2) + x_result(i+3^(k-1),1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            PD_(2,3)=PD_(2,3) + x_result(i+2*3^(k-1),1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            
        end
    end
    %C位点的等位基因的边际外显率
    i=0;
    j=0;
    for o=0:3^(k-1)-1
        for j=1:3
            i = j+(3^(n-2))*o;
            PD_(3,1)=PD_(3,1) + x_result(i,1)*(D(1,i)+D(1,3^1+i)+D(1,2*3^1+i));
            PD_(3,2)=PD_(3,2) + x_result(i+3*1,1)*(D(1,i)+D(1,3^1+i)+D(1,2*3^1+i));
            PD_(3,3)=PD_(3,3) + x_result(i+3*2,1)*(D(1,i)+D(1,3^1+i)+D(1,2*3^1+i));
        end
    end
    %D位点等位基因的边际外显率
    i=0;
    j=0;
    for o=0:2
        for j=0:k
            i = 1+o*3^k+3*j;
            PD_(4,1)=PD_(4,1) + x_result(i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
            PD_(4,2)=PD_(4,2) + x_result(i+1,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
            PD_(4,3)=PD_(4,3) + x_result(i+2,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        end
    end
end
if n==5
    % A的边际外显率的计算公式
    k=4;%order-1
    i =0;
    j=0;
    for i=1:3^k
        PD_(1,1)=PD_(1,1)+x_result(i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_(1,2)=PD_(1,2)+x_result(3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
        PD_(1,3)=PD_(1,3)+x_result(2*3^k+i,1)*(D(1,i)+D(1,3^k+i)+D(1,2*3^k+i));
    end
    % B位点的等位基因的边际外显率
    i=0;
    j=0;
    for o=0:2
        for j = 1:3^(k-1)
            %         i =j+ 9*o;
            i =j+ o*(3^k);
            PD_(2,1)=PD_(2,1) + x_result(i,1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            PD_(2,2)=PD_(2,2) + x_result(i+3^(k-1),1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            PD_(2,3)=PD_(2,3) + x_result(i+2*3^(k-1),1)*(D(1,i)+D(1,3^(k-1)+i)+D(1,2*3^(k-1)+i));
            
        end
    end
    %% C位点的等位基因的边际外显率
    i=0;
    j=0;
    for o=0:3^(n-3)-1%3^(k-2)-1
        for j=1:3^2
            i = j+(3^(n-2))*o;
            PD_(3,1)=PD_(3,1) + x_result(i,1)*(D(1,i)+D(1,3^(k-2)+i)+D(1,2*3^(k-2)+i));
            PD_(3,2)=PD_(3,2) + x_result(i+3^(k-2),1)*(D(1,i)+D(1,3^(k-2)+i)+D(1,2*3^(k-2)+i));
            PD_(3,3)=PD_(3,3) + x_result(i+2*3^(k-2),1)*(D(1,i)+D(1,3^(k-2)+i)+D(1,2*3^(k-2)+i));
        end
    end
    % D位点等位基因的边际外显率
    i=0;
    j=0;
    for o=0:3^(k-1)-1
        for j=1:3
            i = j+(3^(n-3))*o;
            PD_(4,1)=PD_(4,1) +x_result(i,1)*(D(1,i)+D(1,3^(k-3)+i)+D(1,2*3^(k-3)+i));
            PD_(4,2)=PD_(4,2) +x_result(i+3^(k-3),1)*(D(1,i)+D(1,3^(k-3)+i)+D(1,2*3^(k-3)+i));
            PD_(4,3)=PD_(4,3) +x_result(i+2*3^(k-3),1)*(D(1,i)+D(1,3^(k-3)+i)+D(1,2*3^(k-3)+i));
        end
    end
    % E位点等位基因的边际外显率
    i=0;
    j=0;
    for o=0:3^k-1
        %     for j=1:1
        j=1;
        i = j+(3^(n-4))*o;
        PD_(5,1)=PD_(5,1) + x_result(i,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        PD_(5,2)=PD_(5,2) + x_result(i+1,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        PD_(5,3)=PD_(5,3) + x_result(i+2,1)*(D(1,i)+D(1,1+i)+D(1,2+i));
        %     end
    end

end
pd=D*x_result;
h=0;
for i =1:3^n
    h = h +D(i)*(( x_result(i)- pd)^2);
end
h = h/(pd*(1-pd));
if max(abs(pd-min(PD_)))<0.05 & max(abs(pd-max(PD_)))<0.05
    v=1;
else
    v=0;
end
end
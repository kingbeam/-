function [x_opt, fx_opt, iter] = Simplex_eye(A, b, c)
    % 初始化
    [m, n] = size(A);
    table = [0,c;b,A];
    index=zeros(m+1,1);
    theta=zeros(m,1);
    x_opt =zeros(n,1);
    iter = 0;
    c_b=zeros(m+1,1);
    for i=1:m
        index(i+1,1)=n-m+i;
    end
    table=[c_b,index,table];
    sigma=zeros(1,n);
    for i=1:n
        sigma(1,i)=table(1,i+3)-c_b'*table(:,i+3);
        if (sigma(1,i)>=0) && all(table(2:m+1,i+3)<0)
            error('无界解');
        end
    end
    while(max(sigma)>0)
        iter = iter +1;
        [~,max_index]=max(sigma);
        for i=1:m
           theta(i,1)=table(i+1,3)/table(i+1,max_index+3);
        end
        for i=1:m
            if(theta(i,1)<0)
                theta(i,1)=100000000;
            end
        end
        [~,min_index]=min(theta);
        table(min_index+1,3:n+3)=table(min_index+1,3:n+3)/table(min_index+1,max_index+3);
        table(min_index+1,2)=max_index;
        table(min_index+1,1)=c(1,max_index);
        for i=1:m
          if i~=min_index
            table(i+1,3:n+3)=table(i+1,3:n+3)-table(i+1,max_index+3)/table(min_index+1,max_index+3)*table(min_index+1,3:n+3);
        end
        end
   
        c_b=table(2:1+m,1);
       for i=1:n
        sigma(1,i)=table(1,i+3)-c_b'*table(2:m+1,i+3);
        if (sigma(1,i)>=0) && all(table(2:m+1,i+3)<0)
            error('无界解');
        end
       end
    end
    tmp=sigma;
    for i=1:m
        x_opt(table(i+1,2),1)=table(i+1,3);
        tmp(1,table(i+1,2))=-1;
    end
    if any(tmp==0)
        disp('无穷多个最优解');
    end
    fx_opt = c*x_opt;
end

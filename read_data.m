function [x,y]=read_data(flag)
    if flag==1
        load YearPredictionMSD;
        x=full(x);
    end
    if flag==2
        load usps;
        x=full(x);
%         x = x(:,[1:80])
    end
    if flag==3
        load cardiac;
%         x= x(:,[1:400]);
        y=y(:,1);
    end
    if flag==4
        load cardiac;
%         x= x(:,[1:400]);
        y=y(:,2);
    end
    x=zscore(x);
    y=zscore(y);
end

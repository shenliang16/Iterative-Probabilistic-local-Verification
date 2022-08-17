


function  [X, Y, normal] =AffineCPD_normalize(x, y, flag, normal)

if exist('normal', 'var')
    D = size(x, 2);  D_normal = size(normal.xd, 2);
    x = x - repmat(normal.xd, 1, D/D_normal);
    y = y - repmat(normal.yd, 1, D/D_normal);

    X=x./normal.xscale;
    Y=y./normal.yscale;

elseif exist('flag', 'var')
    %% reshape
    D = size(x,2);
    
    %% NaN
    x_ = x';             x_(isnan(x_)) = [];
    y_ = y';             y_(isnan(y_)) = [];
    
    x_reshape = reshape(x_, 2, [])';
    y_reshape = reshape(y_, 2, [])';
    
    n=size(x_reshape,1);
    m=size(y_reshape,1);

    %% 
    normal.xd = mean(x_reshape);   
    normal.yd = mean(y_reshape); %
    
    x = x - repmat(normal.xd, 1, D/2);
    y = y - repmat(normal.yd, 1, D/2);

    %% 
    normal.xscale=sqrt(sum(sum(x.^2,2))/n); 
    normal.yscale=sqrt(sum(sum(y.^2,2))/m);

    X=x./normal.xscale;
    Y=y./normal.yscale;
    
elseif size(x,2)==2
    
    n=size(x,1);
    m=size(y,1);

    %% 
    normal.xd=mean(x);   normal.yd=mean(y); %
    
    x=x-repmat(normal.xd,n,1);
    y=y-repmat(normal.yd,m,1);

    %% 
    normal.xscale=sqrt(sum(sum(x.^2,2))/n); 
    normal.yscale=sqrt(sum(sum(y.^2,2))/m);

    X=x/normal.xscale;
    Y=y/normal.yscale;

elseif size(x,2)==4
    n=size(x,1);
    m=size(y,1);

    %% 
    normal.xd=mean([x(:,1:2); x(:,3:4)]);   normal.yd=mean([y(:,1:2); y(:,3:4)]); %
    
    x=x-repmat(normal.xd,n,2);
    y=y-repmat(normal.yd,m,2);

    %% 
    normal.xscale=sqrt(sum(sum(x.^2,2))/2/n); 
    normal.yscale=sqrt(sum(sum(y.^2,2))/2/m);

    X=x/normal.xscale;
    Y=y/normal.yscale;
    
elseif size(x,2)==6

    n=size(x,1);
    m=size(y,1);

    %% 
    normal.xd=mean([x(:,1:2); x(:,3:4); x(:,5:6)]);   normal.yd=mean([y(:,1:2); y(:,3:4); y(:,5:6)]); %
    
    x=x-repmat(normal.xd,n,3);
    y=y-repmat(normal.yd,m,3);

    %% 
    normal.xscale=sqrt(sum(sum(x.^2,2))/3/n); 
    normal.yscale=sqrt(sum(sum(y.^2,2))/3/m);

    X=x/normal.xscale;
    Y=y/normal.yscale;
    
    
end


function [L mL] = waveletLocalizers(this, lengths, mLengths, model)
normBase = 1.2;

oL = model.distanceMatrix();
h = model.measureOp();
H = h(eye(model.deterministicDimension));
omL = H*oL;

mL = zeros(sum(mLengths(1:end-1)));
pos = 1;
for i = 1:length(mLengths)-1
    M = mLengths(i);
    
    [X Y] = meshgrid(0:1/(length(omL)-1):1);
    [X1 Y1] = meshgrid(0:1/(M-1):1);
    A = interp2(X,Y,omL,X1,Y1,'spline') + (length(mLengths) - i);
    
    
    
    factor = (length(mLengths)-i);
    A = this.localizationFunction(A, normBase^factor*model.decorrelationLength()*2);
    %                 A = A .* normBase^(-length(mLengths)+i);
    
    mL(pos:pos+M-1,pos:pos+M-1) = A;
    
    % compute number of upstream levels
    N = length(mLengths)-i-1;
    
    for r = 1:N
        for n = 1:M
            foo = interp1(0:1/(mLengths(i)-1):1,(A(:,n)),0:1/(mLengths(i+r)-1):1);
            mL(sum(mLengths(1:i+r-1))+1:sum(mLengths(1:i+r-1))+1+mLengths(i+r)-1,pos+n-1) = foo;
            mL(pos+n-1,sum(mLengths(1:i+r-1))+1:sum(mLengths(1:i+r-1))+1+mLengths(i+r)-1) = foo;
        end
    end
    
    pos = pos + M;
end



L = zeros(sum(lengths(1:end-1)), sum(mLengths(1:end-1)));

xpos = 1;
ypos = 1;
for i = 1:length(mLengths)-1
    xM = mLengths(i);
    yM = lengths(i);
    
    [X Y] = meshgrid(0:1/(size(oL,2)-1):1, 0:1/(size(oL,1)-1):1);
    [X1 Y1] = meshgrid(0:1/(xM-1):1, 0:1/(yM-1):1);
    A = interp2(X,Y,oL,X1,Y1,'spline') + (length(mLengths) - i);
    
    
    % FIXME normalization is probably incorrect - change??
    factor = (length(mLengths)-i);
    A = this.localizationFunction(A, normBase^factor*model.decorrelationLength()*2);
    %                 A = A .* normBase^(-length(mLengths)+i);
    
    L(ypos:ypos+yM-1,xpos:xpos+xM-1) = A;
    
    % compute number of upstream levels
    N = length(mLengths)-i-1;
    
    for r = 1:N
        for n = 1:yM
            foo = interp1(0:1/(mLengths(i)-1):1,(A(n,:)),0:1/(mLengths(i+r)-1):1);
            L(ypos+n-1,sum(mLengths(1:i+r-1))+1:sum(mLengths(1:i+r-1))+1+mLengths(i+r)-1) = foo;
        end
        for m = 1:xM
            foo = interp1(0:1/(lengths(i)-1):1,(A(:,m)),0:1/(lengths(i+r)-1):1);
            L(sum(lengths(1:i+r-1))+1:sum(lengths(1:i+r-1))+1+lengths(i+r)-1, xpos+m-1) = foo;
        end
    end
    
    xpos = xpos + xM;
    ypos = ypos + yM;
end
end


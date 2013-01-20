function [ X ] = iwtrans(this, wX, lengths)
s = waverec(wX(:,1), lengths, this.Lo_R, this.Hi_R);
X = zeros(length(s), size(wX,2));
for i = 1:size(wX,2)
    X(:,i) = waverec(wX(:,i), lengths, this.Lo_R, this.Hi_R);
end
end


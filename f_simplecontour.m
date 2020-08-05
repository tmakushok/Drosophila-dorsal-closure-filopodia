%% Adapted from the code by Jérôme Solon
function [inter1D] = f_simplecontour(I)
BW=bwlabel(I,4);
[bx by]=find(BW==0);
for i=1:length(bx)
   if bx(i)~=1 & bx(i)~=size(BW,1) & by(i)~=1 & by(i)~=size(BW,2)
       box=BW(bx(i)-1:bx(i)+1,by(i)-1:by(i)+1);
       m=1;
       for k=1:3
           for l=1:3
           connexions(i,m)=box(k,l);
           m=m+1;
           end
       end
   end   
    
end

for i=1:max(max(connexions))
    k=1;
    l=1;
    for j=1:length(connexions)
        [a,b,c]=find(connexions(j,:)==i);
        c=unique(c);
        if c==1
            boundary(k,1)=bx(j);
            boundary(k,2)=by(j);
            k=k+1;
            [a,b,c]=find(connexions(j,:));
            c=unique(c);
            [xi]=find(c==i);
            c(xi)=[];
            neighbour(l:l-1+length(c))=c(1:length(c));
            l=l+length(c);
        end
    end
    Cell{i}=boundary;
    boundary=0;
    ngbr{i}=unique(neighbour);
    neighbour=0;
end
indInter1D = 1;
for i=1:length(ngbr)
    neighbour=ngbr{i};
    for j=1:length(neighbour)
        inter1D{indInter1D, 1} = intersect(Cell{i},Cell{neighbour(j)},'rows');
        indInter1D = indInter1D + 1;
    end    

end 

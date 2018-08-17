function [auxcoord]=distortedramd
global coord bedge

auxcoord=zeros(size(coord,1),4);

a=-1;
b=1;

ex = a + (b-a)*rand(size(coord,1),1);
ey = a + (b-a)*rand(size(coord,1),1);
% qquando mas proximo a 1 eh mas ditorcido
alpha=0.004; 
h=1/54;

for icoord=1:size(coord,1)
  
     if icoord>size(bedge,1) %&& abs(coord(icoord,1)-0.5)>1e-10 
        x=coord(icoord,1)+alpha*ex(icoord,1);
        
        y=coord(icoord,2)+alpha*ey(icoord,1);
        auxcoord(icoord,1)=icoord;
        auxcoord(icoord,2)=x;
        auxcoord(icoord,3)=y;
    else
        x=coord(icoord,1);
        
        y=coord(icoord,2);
        auxcoord(icoord,1)=icoord;
        auxcoord(icoord,2)=x;
        auxcoord(icoord,3)=y;
    end
end

end

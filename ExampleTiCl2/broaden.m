function yout = broaden(xin,yin,xout,w,btype)
% yout = BROADEN(XIN,YIN,XOUT,W,BTYPE)
%
% Broadens data by a Lorenzian or Gaussian and returns result.
% The input data are column vectors (xin,yin) specifying amplitudes
% yin at position xin.
%
% xout is a column vector of for the x-points
% desired in the broadened function, which is returned in yout and
% has the same size as xout.
%
% w is the broadening width.  
%
% btype = 1 is for Lorenzian and 2 for Gaussian broadening

nin = size(xin,1);
nout = size(xout',1);

yout = zeros(size(xout));
  if (btype == 1)
    if nin < nout       
        for j=1:nin
          yout = yout + yin(j)./((xin(j)-xout).^2+w^2);
        end
    else
        for j=1:nout
          yout(j) = 1./((xin(j)-xout).^2+w^2)'*yin;
        end          
    end
    yout = yout*w/pi;   
    
  else if (btype == 2)
    if nin < nout       
        for j=1:nin
          yout = yout + exp(-(xin(j)-xout).^2/(2*w^2))*yin(j);
        end
    else
        for j=1:nout
          yout(j) = exp(-(xin-xout(j)).^2/(2*w^2))'*yin;
        end          
    end
    yout = yout/sqrt(2*pi*w^2);
  else
    error('Unknown broadening type')
  end
end

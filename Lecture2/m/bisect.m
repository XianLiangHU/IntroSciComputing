function x = bisect(f,a,b,tol)
  xlow = a;  plow  = f(xlow);
  xhigh = b; phigh = f(xhigh);  
  while xhigh - xlow > 2*tol,
    xmid = (xlow + xhigh)/2; pmid = f(xmid);    
    if pmid*plow < 0,
       xhigh = xmid;   phigh = pmid;
    elseif pmid*phigh < 0,
       xlow = xmid;     plow = pmid;
    else
       xlow = xmid;    xhigh = xmid;
    end
  end
  x = [xlow, xhigh];
end

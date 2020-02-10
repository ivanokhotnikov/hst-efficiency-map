function out=fit
load('Data.mat','cat');
y=@(c,x) c(1).*x+c(2);
y1=@(c,x) c(1).*exp(x.*(-c(2)))+c(3);
mp=[cat.brevini.mp;cat.rexroth.mp;cat.eaton.mp;cat.parker.mp;...
    cat.danfoss.mp;cat.kawasaki.mp;cat.hydac.mp;cat.linde.mp;cat.poclain.mp];
Dp=[cat.brevini.Dp;cat.rexroth.Dp;cat.eaton.Dp;cat.parker.Dp;...
    cat.danfoss.Dp;cat.kawasaki.Dp;cat.hydac.Dp;cat.linde.Dp;cat.poclain.Dp];
np=[cat.brevini.np;cat.rexroth.np;cat.eaton.np;cat.parker.np;...
    cat.danfoss.np;cat.kawasaki.np;cat.hydac.np;cat.linde.np;cat.poclain.np];
mm=[cat.rexroth.mm;cat.eaton.mm;cat.danfoss.mm;cat.hydac.mm;...
    cat.parker.mm;cat.kawasaki.mm;cat.linde.mm;cat.poclain.mm];
Dm=[cat.rexroth.Dm;cat.eaton.Dm;cat.danfoss.Dm;cat.hydac.Dm;...
    cat.parker.Dm;cat.kawasaki.Dm;cat.linde.Dm;cat.poclain.Dm];
nm=[cat.rexroth.nm;cat.eaton.nm;cat.danfoss.nm;cat.hydac.nm;...
    cat.parker.nm;cat.kawasaki.nm;cat.linde.nm;cat.poclain.nm];
olsmp = @(c) sum((y(c,Dp) - mp).^2);
olsmm = @(c) sum((y(c,Dm) - mm).^2);
olsnp = @(c) sum((y1(c,Dp) - np).^2);
olsnm = @(c) sum((y1(c,Dm) - nm).^2);
opts = optimset('MaxFunEvals',1e6,'MaxIter',1e6);
Cp = fminsearch(olsmp, rand(2,1), opts);
Cm = fminsearch(olsmm, rand(2,1), opts);
Cnp = fminsearch(olsnp, rand(3,1), opts);
Cnm = fminsearch(olsnm, rand(3,1), opts);
out.C=[Cp,Cm];
out.Cn=[Cnp,Cnm];
rmsmp=sqrt(olsmp(Cp)/(length(mp)-numel(Cp)));
rmsmm=sqrt(olsmp(Cm)/(length(mm)-numel(Cm)));
rmsnp=sqrt(olsnp(Cnp)/(length(np)-numel(Cnp)));
rmsnm=sqrt(olsnm(Cnm)/(length(nm))-numel(Cnm));
resmp=mp-y(Cp,Dp);
resnp=np-y1(Cnp,Dp);
resmm=mm-y(Cm,Dm);
resnm=nm-y1(Cnm,Dm);
r=[corrcoef(Dp,mp),corrcoef(Dm,mm),corrcoef(Dp,np),corrcoef(Dm,nm)];
out.r=r(1,2:2:end);
rres=[corrcoef(Dp,resmp),corrcoef(Dm,resmm),corrcoef(Dp,resnp),...
    corrcoef(Dm,resnm)];
out.rres=rres(1,2:2:end);
out.rmse=[rmsmp,rmsmm,rmsnp,rmsnm];
out.resp=[resmp,resnp];
out.resm=[resmm,resnm];
fit=out;
save('Data.mat','fit','-append')
end

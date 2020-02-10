function printEffMap(s,h,oil,const)
%%  SIZING
% Piston/ cylinder block
s.d=(4*s.Vdmax*1e-6*const.k1/(s.z^2*tan(deg2rad(s.sa))))^(1/3);
s.Ap=pi*s.d^2/4;
s.pcd=s.z*s.d/(pi*const.k1);
s.h=s.pcd*tan(deg2rad(s.sa));
s.eng=1.4*s.d;
s.Ak=const.k3*s.Ap;
s.w=2*(sqrt(s.d^2+(pi-4)*s.Ak)-s.d)/(pi-4);
s.t=const.k2*s.z*s.Ap/(pi*s.pcd)-s.w;
s.rbo=(s.pcd+s.w)/2;
s.Rbo=s.rbo+s.t;
s.Rbi=(s.pcd-s.w)/2;
s.rbi=s.Rbi-s.t;
% Slipper
s.As=const.k4*s.Ap/cos(deg2rad(s.sa));
s.Rs=pi*s.pcd*const.k5/(2*s.z);
s.rs=sqrt(s.Rs^2-s.As/pi);
%% EFFICIENCY CALCULATION
h.n=linspace(h.nmin,h.nmax);
h.pmin=h.pmin*1e5;
h.pmax=h.pmax*1e5;
h.p2=linspace(h.pmin,h.pmax);
h.p1=h.p1*1e5;
for ii=1:max(size(h.n))
    for jj=1:max(size(h.p2))
        eff.ql1(jj)=pi*s.h1^3*0.5*(h.p2(jj)+h.p1)*...
            (1/log(s.Rbo/s.rbo)+1/log(s.Rbi/s.rbi))/(6*oil.mu(4));
        eff.ql2(jj)=s.z*pi*s.h2^3*0.5*(h.p2(jj)+h.p1)/(6*oil.mu(4)*...
            log(s.Rs/s.rs));
        for kk=1:s.z
            eff.ql3i(jj,kk)=s.z*pi*s.d*s.h3^3*0.5*(h.p2(jj)+h.p1)*...
                (1+1.5*s.e^3)*(1/(s.eng+s.h*sin(pi*(kk-1)/s.z)))/...
                (12*oil.mu(4));
        end
        eff.ql3(jj)=sum(eff.ql3i(jj,:));
        eff.ql(jj)=eff.ql1(jj)+eff.ql2(jj)+eff.ql3(jj);
        eff.qth(ii)=h.n(ii)*s.Vdmax/6e7;
        eff.vol_p(ii,jj)=(1-(h.p2(jj)-h.p1)/oil.b-...
            eff.ql(jj)/eff.qth(ii))*100;
        eff.vol_m(ii,jj)=(1-eff.ql(jj)/(eff.qth(ii)))*100;
        eff.vol_hsu(ii,jj)=eff.vol_p(ii,jj)*eff.vol_m(ii,jj)*1e-2;
        eff.mech_p(ii,jj)=(1-const.A*exp(-const.Bp*oil.mu(4)*1e3*h.n(ii)/...
            (s.sa*(h.p2(jj)-h.p1)*1e-5))-...
            const.Cp*sqrt(oil.mu(4)*1e3*h.n(ii)/...
            (s.sa*(h.p2(jj)-h.p1)*1e-5))-...
            const.D/(s.sa*(h.p2(jj)-h.p1)*1e-5))*100;
        eff.mech_m(ii,jj)=(1-const.A*exp(-const.Bm*oil.mu(4)*1e3*h.n(ii)*...
            eff.vol_hsu(ii,jj)*1e-2/(s.sa*(h.p2(jj)-h.p1)*1e-5))-...
            const.Cm*sqrt(oil.mu(4)*1e3*h.n(ii)*eff.vol_hsu(ii,jj)*1e-2/...
            (s.sa*(h.p2(jj)-h.p1)*1e-5))-const.D/...
            (s.sa*(h.p2(jj)-h.p1)*1e-5))*100;
        eff.mech_hsu(ii,jj)=eff.mech_p(ii,jj)*eff.mech_m(ii,jj)*1e-2;
        eff.p(ii,jj)=eff.vol_p(ii,jj)*eff.mech_p(ii,jj)*1e-2;
        eff.m(ii,jj)=eff.vol_m(ii,jj)*eff.mech_m(ii,jj)*1e-2;
        eff.hsu(ii,jj)=eff.p(ii,jj)*eff.m(ii,jj)*1e-2;
        eff.pow_in(ii,jj)=(h.p2(jj)-h.p1)*h.n(ii)*s.Vdmax*(1/6e7)/...
            (eff.mech_p(ii,jj)*1e-2)*1e-3;
        eff.pow_out(ii,jj)=eff.pow_in(ii,jj)*eff.hsu(ii,jj)*1e-2;
        
        eff.n_out(ii,jj)=h.n(ii).*eff.vol_hsu(ii,jj)*1e-2;
        eff.t_in(ii,jj)=(h.p2(jj)-h.p1).*s.Vdmax*1e-6./...
            (2*pi*eff.mech_p(ii,jj)*1e-2);
        eff.t_out(ii,jj)=(h.p2(jj)-h.p1).*s.Vdmax*1e-6./...
            (2*pi*eff.mech_p(ii,jj)*1e-2).*...
            (eff.mech_hsu(ii,jj)*1e-2);
        eff.t_out_check(ii,jj)=eff.pow_out(ii,jj)*1e3./...
            (eff.n_out(ii,jj)*pi/30);
    end
end
%% EFFICIENCY PLOT
figure('Name','HST efficiency','NumberTitle','off');
set(groot,'defaultTextInterpreter','latex');
[C,j]=contour((h.p2-h.p1)*1e-5,h.n,eff.hsu,40:1:100,'ShowText','on',...
    'LineWidth',0.75);
clabel(C,j,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
xlabel('Pressure difference, bar');
ylabel('Pump speed, rpm');
axis([(h.p2(1)-h.p1)*1e-5 (h.p2(end)-h.p1)*1e-5 h.n(1) h.n(end)]);
grid on
box off
%% POWER PLOT
figure('Name','HST output power','NumberTitle','off');
set(groot,'defaultTextInterpreter','latex');
[C,j]=contour((h.p2-h.p1)*1e-5,h.n,eff.pow_out,100:100:1200,...
    'ShowText','on','LineWidth',0.75);
clabel(C,j,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
axis([(h.p2(1)-h.p1)*1e-5 (h.p2(end)-h.p1)*1e-5 h.n(1) h.n(end)]);
xlabel('Pressure difference, bar');
ylabel('Pump speed, rpm');
grid on
box off
end

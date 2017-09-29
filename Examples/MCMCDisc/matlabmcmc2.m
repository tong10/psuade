% This file shows posteriors plots
% ns  - set to 1 for 1-step smoothing of 2D contours
% ns1 - set to 1 for 1-step smoothing of 1D histgrams
ns  = 0;
ns1 = 0;
if exist('noCLF') 
   hold off
else
   clf
end;
active = [
0
1
];
L = [
6.000000e-01 0.000000e+00 ];
U = [
9.000000e-01 1.000000e+00 ];
iStr = {
'X1','X2'};
X = zeros(2,20);
D = zeros(2,20);
NC = zeros(2,2,20,20);
X(1,:) = [
6.075000e-01 6.225000e-01 6.375000e-01 6.525000e-01 6.675000e-01 6.825000e-01 6.975000e-01 7.125000e-01 7.275000e-01 7.425000e-01 7.575000e-01 7.725000e-01 7.875000e-01 8.025000e-01 8.175000e-01 8.325000e-01 8.475000e-01 8.625000e-01 8.775000e-01 8.925000e-01 ];
D(1,:) = [
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 ];
NC(1,2,:,:) = [
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 15000 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
]';
NC(2,1,:,:) = [
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 15000 0 0 0 0 0 0 0 0 0 
]';
X(2,:) = [
2.500000e-02 7.500000e-02 1.250000e-01 1.750000e-01 2.250000e-01 2.750000e-01 3.250000e-01 3.750000e-01 4.250000e-01 4.750000e-01 5.250000e-01 5.750000e-01 6.250000e-01 6.750000e-01 7.250000e-01 7.750000e-01 8.250000e-01 8.750000e-01 9.250000e-01 9.750000e-01 ];
D(2,:) = [
0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 ];
nInps  = length(active);
nPlots = 0;
for ii = 1 : nInps
   if (active(ii) == 1)
      nPlots = nPlots + 1;
      active(ii) = nPlots;
   end;
end;
for ii = 1 : nInps
  for jj = ii : nInps
    if (active(ii) ~= 0 & active(jj) ~= 0)
      index = (active(ii)-1) * nPlots + active(jj);
      subplot(nPlots,nPlots,index)
      if (ii == jj)
        n = length(D(ii,:));
        DN = D(ii,:);
        for kk = 1 : ns1
          DN1 = DN;
          for ll = 2 : n-1
            DN(ll) = DN(ll) + DN1(ll+1);
            DN(ll) = DN(ll) + DN1(ll-1);
            DN(ll) = DN(ll) / 3;
          end;
        end;
        bar(X(ii,:), DN, 1.0);
        xmin = min(X(ii,:));
        xmax = max(X(ii,:));
        xwid = xmax - xmin;
        xmin = xmin - 0.5 * xwid / 20;
        xmax = xmax + 0.5 * xwid / 20;
        ymax = max(DN);
        axis([xmin xmax 0 ymax])
        set(gca,'linewidth',2)
        set(gca,'fontweight','bold')
        set(gca,'fontsize',12)
        xlabel(iStr(ii),'FontWeight','bold','FontSize',12)
        ylabel('Probabilities','FontWeight','bold','FontSize',12)
        grid on
        box on
      else
        n = length(X(jj,:));
        XT = X(jj,:);
        YT = X(ii,:);
        HX = (XT(n) - XT(1)) / (n-1);
        HY = (YT(n) - YT(1)) / (n-1);
        ZZ = squeeze(NC(ii,jj,:,:));
        for kk = 1 : ns
          ZZ1 = ZZ;
          for ll = 2 : n-1
            for mm = 2 : n-1
              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm);
              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm);
              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll,mm+1);
              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll,mm-1);
              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm+1);
              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm-1);
              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll-1,mm+1);
              ZZ(ll,mm) = ZZ(ll,mm) + ZZ1(ll+1,mm-1);
              ZZ(ll,mm) = ZZ(ll,mm) / 9;
            end;
          end;
        end;
        ZZ = ZZ / (sum(sum(ZZ)));
%      [YY,XX]=meshgrid(XT(1):HX:XT(n),YT(1):HY:YT(n));
%      HX = 0.01 * (XT(n) - XT(1));
%      HY = 0.01 * (YT(n) - YT(1));
%      [YI,XI]=meshgrid(XT(1):HX:XT(n),YT(1):HY:YT(n));
%      ZI=interp2(YY, XX, ZZ, YI, XI, 'spline');
%      pcolor(XI,YI,ZI)
%      shading interp
%      hold on
%      contour(XI,YI,ZI,5,'k')
        imagesc(ZZ')
        xtick = L(ii):(U(ii)-L(ii))/4:U(ii);
        set(gca,'XTick',0:n/4:n);
        set(gca,'XTickLabel', xtick);
        ytick = L(jj):(U(jj)-L(jj))/4:U(jj);
        set(gca,'YTick',0:n/4:n);
        set(gca,'YTickLabel', ytick);
        set(gca,'YDir', 'normal');
        xlabel(iStr(jj),'FontWeight','bold','FontSize',12)
        ylabel(iStr(ii),'FontWeight','bold','FontSize',12)
set(gca,'linewidth',2)
set(gca,'fontweight','bold')
set(gca,'fontsize',12)
box on
      end;
    end;
  end;
end;
      subplot(nPlots,nPlots,1)
set(gcf,'NextPlot','add');
axes;
h=title('MCMC Posterior Distributions, neg. log likelihood=5.554075e+03','fontSize',12,'fontWeight','bold');
set(gca,'Visible','off');
set(h,'Visible','on');
negll = 5.554075e+03;

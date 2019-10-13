% Plot collagen distribution cdf & pdf
function [pdf1,pdf2,cdf1,cdf2] = collagen_disp(LR,DR,a,b)
% Set up the lamina propria distribution
pd1 = makedist('Triangular','a',LR(1),'b',LR(3),'c',LR(2));
% Set up the destrusor layer distribution
pd2 = makedist('Triangular','a',DR(1),'b',DR(3),'c',DR(2));
Xplot = linspace(1, 2);
pdf1 = pdf(pd1,Xplot);
pdf2 = pdf(pd2,Xplot);
cdf1 = cdf(pd1,Xplot);
cdf2 = cdf(pd2,Xplot);
%% Plot the pdf
figure(a)
hold on
plot(Xplot,pdf1,'--b','LineWidth',4)
plot(Xplot,pdf2,'--r','LineWidth',4)
xlim([1,2])
xlabel('Collagen recruitment stretch')
ylabel('Probability density')
set(gca,'fontsize',15)
%% Plot the cdf
figure(b)
hold on
plot(Xplot,cdf1,'--b','LineWidth',4)
plot(Xplot,cdf2,'--r','LineWidth',4)
xlim([1,2])
xlabel('Stretch')
ylabel('Collagen Fiber Recruitment Fraction')
set(gca,'fontsize',15)
end

clear all
load('borisColl.mat')
GP=load('group_impGitHub.mat'); % generelle informationen ueber das hefe genom

for i=1:numel(GP.TFs)
    aidOM{GP.TFs(i),i}=NaN;
    aidOP{GP.TFs(i),i}=NaN;
    
    delOM{GP.TFs(i),i}=NaN;
    delOP{GP.TFs(i),i}=NaN;

    ideaO{GP.TFs(i),i}=NaN;
end

%%
ideaZ=(ideaO{:,:}-mean(ideaO{:,:},2,'omitnan'))./std(ideaO{:,:},[],2,'omitnan');
delZ=(delOM{:,:}-mean(delOM{:,:},2,'omitnan'))./std(delOM{:,:},[],2,'omitnan');
TFstats.ZRegIdea=diag(corr(ideaZ,combBind,'rows','pairwise'))
TFstats.ZRegDel=diag(corr(delZ,combBind,'rows','pairwise'))


%% figure 1
foldTH=log2(1.5)
pTHaid=0.05;
pTHdel=0.05;
%aidOM{:,:}=-aidOM{:,:};
%delOM{:,:}=-delOM{:,:};


TFstats=table();

TFstats.isDel=(double(sum(~isnan(delOM{:,:}),'omitnan')>2000))';
TFstats.isAid=(double(sum(~isnan(aidOM{:,:}),'omitnan')>2000))';
TFstats.isIdea=(double(sum(~isnan(ideaO{:,:}),'omitnan')>2000))';


TFstats.DELup=sum((delOM{:,:}>foldTH)&(delOP{:,:}<pTHdel))';
TFstats.DELdown=sum((delOM{:,:}<-foldTH)&(delOP{:,:}<pTHdel))';
TFstats{TFstats.isDel==0,{'DELdown','DELup'}}=nan;
TFstats.Properties.RowNames=delOM.Properties.VariableNames;
TFstats.AIDup=sum((aidOM{:,:}>foldTH)&(aidOP{:,:}<pTHaid))';
TFstats.AIDdown=sum((aidOM{:,:}<-foldTH)&(aidOP{:,:}<pTHaid))';
TFstats{TFstats.isAid==0,{'AIDdown','AIDup'}}=NaN;
TFstats.IDEAup=sum(ideaO{:,:}>foldTH)';
TFstats.IDEAdown=sum(ideaO{:,:}<-foldTH)';
TFstats{TFstats.isIdea==0,{'IDEAup','IDEAdown'}}=NaN;

TFstats.AIDfrac=(TFstats.AIDup-TFstats.AIDdown)./(TFstats.AIDup+TFstats.AIDdown+1)
TFstats.IDEAfrac=(TFstats.IDEAup-TFstats.IDEAdown)./(TFstats.IDEAup+TFstats.IDEAdown+1)
TFstats.DELfrac=(TFstats.DELup-TFstats.DELdown)./(TFstats.DELup+TFstats.DELdown+1)

TFdesc=readtable("TFdesc.xlsx")
[~,idx]=ismember(TFdesc.Var1,GP.gene_infoR64.orf(GP.TFs))
TFstats.sgd(idx(idx>0))=TFdesc.Act(idx>0)
TFstats.sgd(isnan(TFstats.sgd))=0;
commonTFs=find(sum(TFstats{:,1:3},2)==3);

ideaOS=ideaO{:,:};
ideaOS(ideaOS==0)=nan;

delOMS=delOM{:,:};
delOMS(delOP{:,:}>pTHdel)=nan;

aidOMS=aidOM{:,:};
aidOMS(aidOP{:,:}>pTHaid)=nan;


figure
subplot(3,5,1)
hold off
randGenes=rand(6701,1)>0.75;
scatter(delOM.MSN2(randGenes),min(-log10(delOP.MSN2(randGenes)),20),'.')
xlim([-2 2])
hold on
plot(xlim(),-log10([pTHdel]).*[1,1],'k-')
selGenes=(delOM.MSN2>foldTH)&(delOP.MSN2<pTHdel)
scatter(delOM.MSN2(selGenes),min(-log10(delOP.MSN2(selGenes)),20),'.')
selGenes=(delOM.MSN2<-foldTH)&(delOP.MSN2<pTHdel)
scatter(delOM.MSN2(selGenes),min(-log10(delOP.MSN2(selGenes)),20),'.')
hold on
plot(foldTH.*[-1 1;-1 1],ylim()','k-')
title('MSN2')

subplot(3,5,2)
hold off
scatter(aidOM.GCN4(randGenes),min(-log10(aidOP.GCN4(randGenes)),20),'.')
xlim([-2 2])
hold on
plot(xlim(),-log10([pTHaid]).*[1,1],'k-')
selGenes=(aidOM.GCN4>foldTH)&(aidOP.GCN4<pTHaid);
scatter(aidOM.GCN4(selGenes),min(-log10(aidOP.GCN4(selGenes)),20),'.')

selGenes=(aidOM.GCN4<-foldTH)&(aidOP.GCN4<pTHaid);
scatter(aidOM.GCN4(selGenes),min(-log10(aidOP.GCN4(selGenes)),20),'.')
hold on
plot(foldTH.*[-1 1;-1 1],ylim()','k-')
title('GCN4')

subplot(3,5,3)
hold off
[val,idx]=sort(ideaO.NRG1,'MissingPlacement','last')
idx(isnan(val))=[];
val(isnan(val))=[];
scatter(find(randGenes(1:numel(idx))),val(randGenes(1:numel(idx))),'.')
hold on
scatter(find(val>foldTH),val(val>foldTH),'.')
scatter(find(val<-foldTH),val(val<-foldTH),'.')
plot(xlim()',foldTH.*[-1 1;-1 1],'k-')
title('NRG1')

% subplot(3,5,5+[1:4])
% 
% [~,idx]=sort(median(TFstats{commonTFs,4:9},2))
% 
% selCols=4:9;
% hold off
% imagesc(log10(TFstats{commonTFs(idx),selCols}+0.3)')
% set(gca,'XTick',1:numel(commonTFs),'XTickLabel',TFstats.Properties.RowNames(commonTFs(idx)),'YTick',1:numel(selCols),'YTickLabel',TFstats.Properties.VariableNames(selCols))
% colorbar()
% colormap(gca,brighten(brewermap(128,'greens'),0.3))
% hold on
% plot(xlim()',[1:5]+[.5;.5],'Color',[1 1 1].*0.75)
% plot([1:numel(idx)]+[.5;.5],ylim()','Color',[1 1 1].*0.75)
% 
% 
% subplot(3,5,5+5)
% hold off
% imagesc(corr(log2(TFstats{commonTFs(idx),selCols}+0.1)),[0 1]);
% colorbar()
% colormap(gca,brighten(brewermap(128,'blues'),0.3))
% hold on
% plot(xlim()',[1:5]+[.5;.5],'Color',[1 1 1].*0.75)
% plot([1:numel(idx)]+[.5;.5],ylim()','Color',[1 1 1].*0.75)


%% 


[~,selCols]=ismember({'AIDfrac','DELfrac','IDEAfrac'},TFstats.Properties.VariableNames);
[~,idx]=sort(median(TFstats{commonTFs,selCols},2))
subplot(3,5,10+[1:4])
hold off
imagesc(TFstats{commonTFs(idx),selCols}')
hold on
set(gca,'XTick',1:numel(commonTFs),'XTickLabel',TFstats.Properties.RowNames(commonTFs(idx)),'YTick',1:numel(selCols),'YTickLabel',TFstats.Properties.VariableNames(selCols))
colormap(gca,brighten([flipud(brewermap(128,'Reds'));brewermap(128,'Greens')],0.3))
plot(xlim()',[1:5]+[.5;.5],'Color',[1 1 1].*0.75)
plot([1:numel(idx)]+[.5;.5],ylim()','Color',[1 1 1].*0.75)
colorbar()
subplot(3,5,15)
imagesc(corr(TFstats{commonTFs,selCols}),[0 1])
colormap(gca,brighten(brewermap(128,'blues'),0.3))
colorbar()
hold on
plot(xlim()',[1:numel(selCols)]+[.5;.5],'Color',[1 1 1].*0.75)
plot([1:numel(selCols)]+[.5;.5],ylim()','Color',[1 1 1].*0.75)
saveas(gcf,'Figure1partAimp.svg')
saveas(gcf,'Figure1partA.svg')


%% second part

TFstats.AidDel=diag(corr(aidOM{:,:},delOM{:,:},'rows','pairwise'));
TFstats.DelIdea=diag(corr(delOM{:,:},ideaO{:,:},'rows','pairwise'));
TFstats.AidIdea=diag(corr(aidOM{:,:},ideaO{:,:},'rows','pairwise'));

sigTarget=(aidOMS>foldTH)+(delOMS>foldTH)+(ideaOS>foldTH)-(aidOMS<-foldTH)-(delOMS<-foldTH)-(ideaOS<-foldTH);

figure

%% MIG1
TFs={'MIG1','MIG2'}
[~,tfPos]=ismember(TFs,GP.gene_infoR64.nameNew(GP.TFs));

randGenes=rand(6701,1)>0.9;
subplot(3,5,1)
hold off
scatter(delOM.(TFs{1})(randGenes),ideaO.(TFs{1})(randGenes),'.')
hold on
highGenes=abs(ideaO.(TFs{1}))>0.05 | abs(delOM.(TFs{1}))>foldTH/2;
scatter(delOM.(TFs{1})(highGenes),ideaO.(TFs{1})(highGenes),'.')

scatter(delOM.(TFs{1})(abs(sigTarget(:,tfPos(1)))>1),ideaO.(TFs{1})(abs(sigTarget(:,tfPos(1)))>1),'.')

[sum(sigTarget(:,tfPos(1))>1) sum(sigTarget(:,tfPos(1))<-1)]

xlim([-1.5 1])
ylim([-3 1.5])
plot([0 0],ylim(),'k-')
plot(xlim(),[0 0],'k-')
xlabel('DEL')
ylabel('IDEA')
title(sprintf('%s %.2f',TFs{1},corr(delOM.(TFs{1}),ideaO.(TFs{1}),'rows','pairwise')))
%% MGA1
randGenes=rand(6701,1)>0.9;
subplot(3,5,6)
hold off
scatter(delOM.(TFs{2})(randGenes),aidOM.(TFs{2})(randGenes),'.')
hold on
highGenes=abs(aidOM.(TFs{2}))>0.05 | abs(delOM.(TFs{2}))>foldTH/2;
scatter(delOM.(TFs{2})(highGenes),aidOM.(TFs{2})(highGenes),'.')

hold on
xlim([-.8 .8])
ylim([-4 2.2])
plot([0 0],ylim(),'k-')
plot(xlim(),[0 0],'k-')
xlabel('DEL')
ylabel('AID')
title(sprintf('%s %.2f',TFs{2},corr(delOM.(TFs{2}),aidOM.(TFs{2}),'rows','pairwise')))




% subplot(3,5,[1:4]+1)
% [~,idx]=sort(mean(TFstats{commonTFs,10:12},2))
% selCols=[11,10,12];
% imagesc(TFstats{commonTFs(idx),selCols}')
% colorbar()
% set(gca,'XTick',1:numel(commonTFs),'XTickLabel',TFstats.Properties.RowNames(commonTFs(idx)),'YTick',1:numel(selCols),'YTickLabel',TFstats.Properties.VariableNames(selCols))
% colormap(gca,brighten(brewermap(128,'Rdbu'),0.3))
% caxis([-.5 .5])
% hold on
% plot(xlim()',[1:numel(selCols)-1]+[.5;.5],'Color',[1 1 1].*0.75)
% plot([1:numel(idx)]+[.5;.5],ylim()','Color',[1 1 1].*0.75)


subplot(3,5,10)
hold off
[~,selCols]=ismember({'AidDel','DelIdea','AidIdea'},TFstats.Properties.VariableNames);
plot([-1:0.1:1],movmean(histcounts(TFstats{commonTFs,selCols(1)}',[-1.05:0.1:1.05],'Normalization','count'),2),'LineWidth',2,'DisplayName',TFstats.Properties.VariableNames{10})
hold on
plot([-1:0.1:1],movmean(histcounts(TFstats{commonTFs,selCols(2)}',[-1.05:0.1:1.05],'Normalization','count'),2),'LineWidth',2,'DisplayName',TFstats.Properties.VariableNames{11})
plot([-1:0.1:1],movmean(histcounts(TFstats{commonTFs,selCols(3)}',[-1.05:0.1:1.05],'Normalization','count'),2),'LineWidth',2,'DisplayName',TFstats.Properties.VariableNames{12})
xlim([-.2 1])
saveas(gcf,'Figure1partB.svg')






TFstats.dblUp=sum(sigTarget>1)';
TFstats.dblDown=sum(sigTarget<-1)';
TFstats.dblScore=(TFstats.dblUp-TFstats.dblDown)./(TFstats.dblUp+TFstats.dblDown+1);
subplot(3,5,[6 12])
scatter(log2(TFstats.dblUp(commonTFs)+TFstats.dblDown(commonTFs)+.5),TFstats.dblScore(commonTFs),200,TFstats.sgd(commonTFs),'filled')
text(log2(TFstats.dblUp(commonTFs)+TFstats.dblDown(commonTFs)+.5),TFstats.dblScore(commonTFs),TFstats.Properties.RowNames(commonTFs))
saveas(gcf,'Figure1partB.svg')

%TFstats.AidDelS=diag(corr(aidOMS,delOMS,'rows','pairwise'));
%TFstats.AidIdeaS=diag(corr(aidOMS,ideaOS,'rows','pairwise'));
%TFstats.DelIdeaS=diag(corr(delOMS,ideaOS,'rows','pairwise'));

%% Figure S1

subplot(3,5,[1:4]+10)
selCols=[13:15];
[~,idx]=sort(mean(TFstats{commonTFs,selCols},2,'omitmissing'),'missingplacement','first')
imagesc(TFstats{commonTFs(idx),selCols}')
caxis([-1.1 1])
colormap(gca,brighten(brewermap(128,'Rdbu'),0.3))
colorbar()


set(gca,'XTick',1:numel(commonTFs),'XTickLabel',TFstats.Properties.VariableNames(commonTFs(idx)))
subplot(3,5,15)
hold off
plot([-1:0.1:1],histcounts(TFstats{commonTFs,13}',[-1.05:0.1:1.05],'Normalization','cumcount'))
hold on
plot([-1:0.1:1],histcounts(TFstats{commonTFs,14}',[-1.05:0.1:1.05],'Normalization','cumcount'))
plot([-1:0.1:1],histcounts(TFstats{commonTFs,15}',[-1.05:0.1:1.05],'Normalization','cumcount'))

%% figure2

%combine expression datasets
clear combExp
combExp(:,:,1)=ideaO{:,:}./max(abs(ideaO{:,:}));
combExp(:,:,2)=aidOM{:,:}./max(abs(aidOM{:,:}));
%combExp(:,:,3)=delOM{:,:}./max(abs(delOM{:,:}));
combExp=mean(combExp,3,'omitmissing');

%combined binding DS
combBind(:,:,1)=barkaiSPorder{:,:};
combBind(:,:,2)=hahnSPorder{:,:};
combBind=mean(combBind,3,'omitmissing');

TFstats.isBarkai=sum(~isnan(barkaiSPorder{:,:}))'>1000;
TFstats.isHahn=sum(~isnan(hahnSPorder{:,:}))'>1000;
TFstats.isBH=sum(~isnan(combBind))'>1000;

forTarget=combBind;
ideaOT=table2array(ideaO);
ideaOT((forTarget<1e4)|isnan(forTarget))=nan;
delOMT=table2array(delOM);
delOMT((forTarget<1e4)|isnan(forTarget))=nan;
aidOMT=table2array(aidOM);
aidOMT((forTarget<1e4)|isnan(forTarget))=nan;
combExpT=combExp;
combExpT((forTarget<1e4)|isnan(forTarget))=nan;

TFstats.tgtMeanIdea=mean(ideaOT,'omitnan')';
TFstats.tgtMeanAid=mean(aidOMT,'omitnan')';
TFstats.tgtMeanDel=mean(delOMT,'omitnan')';
TFstats.tgtMeanExp=mean(combExpT,'omitnan')';

allCmp=find(sum(TFstats{:,{'isAid','isDel','isIdea','isBH'}},2)==4);

ideaOST=ideaOS;
ideaOST((forTarget<1e4)|isnan(forTarget))=nan;
delOMST=delOMS;
delOMST((forTarget<1e4)|isnan(forTarget))=nan;
aidOMST=aidOMS;
aidOMST((forTarget<1e4)|isnan(forTarget))=nan;

TFstats.tgtUpfracIdea=(sum(ideaOST>foldTH)./sum(~isnan(ideaOT)))'; %singificantly  upregulated targets
TFstats.tgtDownfracIdea=(sum(ideaOST<-foldTH)./sum(~isnan(ideaOT)))'; %singificantly  upregulated targets
TFstats.tgtFracIdea=TFstats.tgtUpfracIdea-TFstats.tgtDownfracIdea;

TFstats.tgtUpfracDel=(sum(delOMST>foldTH)./sum(~isnan(delOMT)))'; %singificantly  upregulated targets
TFstats.tgtDownfracDel=(sum(delOMST<-foldTH)./sum(~isnan(delOMT)))'; %singificantly  upregulated targets
TFstats.tgtFracDel=TFstats.tgtUpfracDel-TFstats.tgtDownfracDel;

TFstats.tgtUpfracAid=(sum(aidOMST>foldTH)./sum(~isnan(aidOMT)))'; %singificantly  upregulated targets
TFstats.tgtDownfracAid=(sum(aidOMST<-foldTH)./sum(~isnan(aidOMT)))'; %singificantly  upregulated targets
TFstats.tgtFracAid=TFstats.tgtUpfracAid-TFstats.tgtDownfracAid;

figure;
subplot(1,2,1)
histogram((sum(combBind(:,allCmp)>1e4)),8)
xlabel('# binding targets')
ylabel('# TFs')

subplot(1,2,2)
hold off
[~,selTF]=ismember('SOK2',TFstats.Properties.RowNames)
scatter(combBind(:,selTF)/1e4,delOM{:,selTF},'.');
ylim([-2 2])
hold on
plot([1 1],ylim,'k--')
plot(xlim,[1 1].*foldTH,'k--')
plot(xlim,[1 1].*-foldTH,'k--')

selGenes=delOMS(:,selTF)>foldTH&combBind(:,selTF)>1e4;
scatter(combBind(selGenes,selTF)/1e4,delOM{selGenes,selTF},'filled');
selGenes=delOMS(:,selTF)<-foldTH&combBind(:,selTF)>1e4;
scatter(combBind(selGenes,selTF)/1e4,delOM{selGenes,selTF},'filled');
title(sprintf('%s %d %d %d',TFstats.Properties.RowNames{selTF},sum(delOMS(:,selTF)>foldTH&combBind(:,selTF)>1e4),sum(delOMS(:,selTF)<-foldTH&combBind(:,selTF)>1e4),sum(combBind(:,selTF)>1e4)))
saveas(gcf,'Figure2pA.svg')


figure;
subplot(3,5,1)
scatter(TFstats.tgtUpfracIdea(allCmp),TFstats.tgtDownfracIdea(allCmp),[],TFstats.tgtMeanIdea(allCmp),'filled')
xlabel('fraction up')
ylabel('fraction down')
ylabel(colorbar(),'meanExp')
title('Idea')
caxis([-1 1].*0.9.*max(abs(caxis)))
colormap(gca,brewermap(128,'RdYlGn'))

subplot(3,5,2)
scatter(TFstats.tgtUpfracDel(allCmp),TFstats.tgtDownfracDel(allCmp),[],TFstats.tgtMeanDel(allCmp),'filled')
xlabel('fraction up')
ylabel('fraction down')
ylabel(colorbar(),'meanExp')
title('Del')
caxis([-1 1].*0.9.*max(abs(caxis)))
colormap(gca,brewermap(128,'RdYlGn'))

subplot(3,5,3)
scatter(TFstats.tgtUpfracAid(allCmp),TFstats.tgtDownfracAid(allCmp),[],TFstats.tgtMeanAid(allCmp),'filled')
xlabel('fraction up')
ylabel('fraction down')
ylabel(colorbar(),'meanExp')
title('AID')
caxis([-1 1].*0.9.*max(abs(caxis)))
colormap(gca,brewermap(128,'RdYlGn'))


[~,selCols]=ismember({'tgtFracIdea','tgtFracDel','tgtFracAid';'tgtMeanIdea','tgtMeanDel','tgtMeanAid'},TFstats.Properties.VariableNames)


subplot(3,5,4)
hold off
imagesc(tril(corr(TFstats{allCmp,selCols(1,:)}),-1)+triu(corr(TFstats{allCmp,selCols(2,:)}),1)+triu(tril(corr(TFstats{allCmp,selCols(1,:)},TFstats{allCmp,selCols(2,:)}))),[0 1])
colorbar()
colormap(gca,brighten(brewermap(128,'blues'),0.3))
hold on
plot(xlim()',[1:numel(selCols)-1]+[.5;.5],'Color',[1 1 1].*0.75)
plot([1:numel(selCols)-1]+[.5;.5],ylim()','Color',[1 1 1].*0.75)


% subplot(3,5,5)
% hold off
% imagesc(corr(TFstats{allCmp,selCols}),[0 1])
% colorbar()
% colormap(gca,brighten(brewermap(128,'blues'),0.3))
% hold on
% plot(xlim()',[1:numel(selCols)-1]+[.5;.5],'Color',[1 1 1].*0.75)
% plot([1:numel(selCols)-1]+[.5;.5],ylim()','Color',[1 1 1].*0.75)

saveas(gcf,'Figure2pB.svg')

% [~,selCols]=ismember({'tgtMeanIdea','tgtMeanAid','tgtMeanDel'},TFstats.Properties.VariableNames);
% figure
% subplot(3,5,1:4)
% imagesc((TFstats{allCmp(idx),selCols}./max(abs(TFstats{allCmp(idx),selCols})))',[-1 1])
% caxis([-1 1])
% subplot(3,5,[1:4]+5)
% bar(log10(sum(combBind(:,allCmp(idx))>=1e4)))
% ylim([.5 2.3])
% colorbar()

%% Figure 2 part 2
% TFstats.CrBarkaiAid=diag(corr(barkaiSPorder{:,:},aidOM{:,:},'rows','pairwise'));
% TFstats.CrBarkaiIdea=diag(corr(barkaiSPorder{:,:},ideaO{:,:},'rows','pairwise'));
% TFstats.CrBarkaiDel=diag(corr(barkaiSPorder{:,:},delOM{:,:},'rows','pairwise'));
% 
% TFstats.CrBarkaiAidS=diag(corr(barkaiSPorder{:,:},aidOMS,'rows','pairwise'));
% TFstats.CrBarkaiIdeaS=diag(corr(barkaiSPorder{:,:},ideaOS,'rows','pairwise'));
% TFstats.CrBarkaiDelS=diag(corr(barkaiSPorder{:,:},delOMS,'rows','pairwise'));
% 
% TFstats.CrHahnAid=diag(corr(hahnSPorder{:,:},aidOM{:,:},'rows','pairwise'));
% TFstats.CrHahnIdea=diag(corr(hahnSPorder{:,:},ideaO{:,:},'rows','pairwise'));
% TFstats.CrHahnDel=diag(corr(hahnSPorder{:,:},delOM{:,:},'rows','pairwise'));
% 
% TFstats.CrHahnAidS=diag(corr(hahnSPorder{:,:},aidOMS,'rows','pairwise'));
% TFstats.CrHahnIdeaS=diag(corr(hahnSPorder{:,:},ideaOS,'rows','pairwise'));
% TFstats.CrHahnDelS=diag(corr(hahnSPorder{:,:},delOMS,'rows','pairwise'));
% 
% TFstats.CrHahnComb=diag(corr(hahnSPorder{:,:},combExp,'rows','pairwise'));
% TFstats.CrBarkaiComb=diag(corr(barkaiSPorder{:,:},combExp,'rows','pairwise'));
% 
% TFstats.CrBHComb=diag(corr(combBind,combExp,'rows','pairwise'));



TFstats.CrBHIdea=diag(corr(combBind,ideaO{:,:},'rows','pairwise'));
TFstats.CrBHDel=diag(corr(combBind,delOM{:,:},'rows','pairwise'));
TFstats.CrBHAid=diag(corr(combBind,aidOM{:,:},'rows','pairwise'));

TFstats.CrBHIdeaS=diag(corr(combBind,ideaOS,'rows','pairwise'));
TFstats.CrBHDelS=diag(corr(combBind,delOMS,'rows','pairwise'));
TFstats.CrBHAidS=diag(corr(combBind,aidOMS,'rows','pairwise'));

%% add known data
temp=load('tfActivationDomains.mat')
[~,bidx]=ismember(regexprep(temp.tf.names,{'CAD1','CIN5','PUL4','WTM2'},{'YAP2','YAP4','YNR063W','YOR229W'}),TFstats.Properties.RowNames)
temp.tf.actScore=cellfun(@(x)sum(x(x>10),'omitnan'),temp.tf.actProfile)
TFstats.actScore(:)=nan;
TFstats.actScore(bidx(bidx>0))=temp.tf.actScore(bidx>0);

biogrid=readtable('delAnalysis.xlsx','Sheet','Paper','ReadVariableNames',false)
[~,biogrid.bid]=ismember(regexprep(biogrid.Var1,{'CAD1','CIN5','PUL4','WTM2'},{'YAP2','YAP4','YNR063W','YOR229W'}),GP.gene_infoR64.nameNew)
TFstats.cyc8=ismember(GP.TFs,biogrid.bid)

delAnalysis=readtable('delAnalysis.xlsx')
delAnalysis(cellfun('isempty',delAnalysis.Gene),:)=[];
delAnalysis.citeScore=contains(delAnalysis.Literature_combinedFromSGD_Http___www_yeastgenome_org_CherryEtA,'activator')-contains(delAnalysis.Literature_combinedFromSGD_Http___www_yeastgenome_org_CherryEtA,'repressor')
[~,delAnalysis.bid]=ismember(regexprep(delAnalysis.Gene,{'CAD1','CIN5','PUL4','WTM2'},{'YAP2','YAP4','YNR063W','YOR229W'}),GP.gene_infoR64.nameNew)
delAnalysis(delAnalysis.bid==0,:)=[];
TFstats.litScore(:)=0;
[~,bidx]=ismember(GP.TFs,delAnalysis.bid);
TFstats.litScore(bidx>0)=delAnalysis.citeScore(bidx(bidx>0))

clearvars delAnalysis bidx biogrid temp

%% Figure 2 part 3
figure
[~,selCols]=ismember({'CrBHIdea','CrBHDel','CrBHAid'},TFstats.Properties.VariableNames)

subplot(3,5,[6 9])
[~,idx]=sort(mean(TFstats{allCmp,selCols},2))
hold off
imagesc(TFstats{allCmp(idx),selCols}')
colormap(gca,brighten([flipud(brewermap(128,'Reds'));brewermap(128,'Greens')],0.3))
caxis([-1 1].*max(abs(caxis)))
hold on
plot(xlim()',[1:numel(selCols)-1]+[.5;.5],'Color',[1 1 1].*0.75)
plot([1:numel(idx)-1]+[.5;.5],ylim()','Color',[1 1 1].*0.75)
set(gca,'XTick',1:numel(commonTFs),'XTickLabel',TFstats.Properties.RowNames(commonTFs(idx)),'YTick',1:numel(selCols),'YTickLabel',TFstats.Properties.VariableNames(selCols))
colorbar()

subplot(3,5,10)
hold off
imagesc(corr(TFstats{allCmp,selCols}),[0 1])
colorbar()
colormap(gca,brighten(brewermap(128,'blues'),0.3))
hold on
plot(xlim()',[1:numel(selCols)-1]+[.5;.5],'Color',[1 1 1].*0.75)
plot([1:numel(selCols)-1]+[.5;.5],ylim()','Color',[1 1 1].*0.75)


clear temp
subplot(9,5,6*5+[1:4])
hold off
imagesc(log2(TFstats.actScore(allCmp(idx))'+10),[3 12])
hold on
plot([1:numel(idx)-1]+[.5;.5],ylim()','Color',[1 1 1].*0.75)
colorbar()
colormap(gca,[.75 .75 .75;brighten(brewermap(128,'greens'),0.3)])

subplot(9,5,6*5+5)
imagesc(corr(log2(TFstats.actScore(allCmp)+0.1),TFstats{allCmp,selCols},'rows','pairwise'),[0 1])
colorbar()
colormap(gca,brighten(brewermap(128,'blues'),0.3))



subplot(9,5,7*5+[1:4])
imagesc(TFstats.sgd(allCmp(idx))',[-1.5 1.5])
colormap(gca,brighten([flipud(brewermap(128,'Reds'));brewermap(128,'Greens')],0.3))
hold on
plot([1:numel(idx)-1]+[.5;.5],ylim()','Color',[1 1 1].*0.75)
colorbar()
subplot(9,5,7*5+5)
imagesc(corr(TFstats.sgd(allCmp),TFstats{allCmp,selCols},'rows','pairwise'),[0 1])
colorbar()
colormap(gca,brighten(brewermap(128,'blues'),0.3))


subplot(9,5,8*5+[1:4])
imagesc(-TFstats.cyc8(allCmp(idx))',[-1.5 1.5])
colormap(gca,brighten([flipud(brewermap(128,'Reds'));brewermap(128,'Greens')],0.3))
hold on
plot([1:numel(idx)-1]+[.5;.5],ylim()','Color',[1 1 1].*0.75)
colorbar()

subplot(9,5,8*5+5)
imagesc(corr(-TFstats.cyc8(allCmp),TFstats{allCmp,selCols},'rows','pairwise'),[0 1])
colorbar()
colormap(gca,brighten(brewermap(128,'blues'),0.3))
saveas(gcf,'Figure2pC.svg')

%% compare binding against motifs

allBS=readtable("TFanalysis-BorisMherHuett\bestPBMonScer.txt")
GP=load('group_impGitHub.mat'); % generelle informationen ueber das hefe genom
[~,allBS.tfid]=ismember(allBS.motif_alt_id,GP.gene_infoR64.orf)
[~,allBS.tfPos]=ismember(allBS.tfid,GP.TFs)
load('promoterIDXvecFimp.mat')
[~,allBS.chrID]=ismember(allBS.sequence_name,GP.chrNames)
allBS.pos=GP.chrIdx(allBS.chrID)+round((allBS.start+allBS.stop)/2);
allBS.type=promoterIDXvecF(allBS.pos)
selBS=allBS.type==1;
motifMat=sparse(allBS.pos(selBS),allBS.tfPos(selBS),1,GP.chrIdx(end),numel(GP.TFs));
motifGeneMat=nan(6701,width(motifMat));

for i=find(full(sum(motifMat)>0)
    motifGeneMat(:,i)=sum(metaProfilePromLenDivya(full(motifMat(:,i)),'promEnd','position','promLen',promoterLengthsORF,'afterTSS',0),2,'omitnan');
end
clearvars motifMat allBS i promoterIDXvecF promoterLengthsORF

TFstats.hasMotif=~all(isnan(motifGeneMat))';
TFstats.hasBH=~all(isnan(combBind))';

for i=find(TFstats.hasMotif & TFstats.hasBH)'
    targets=combBind(:,i)>1e4;
    TFstats.motTarget(i)=mean(motifGeneMat(targets,i));
    TFstats.motNonTarget(i)=mean(motifGeneMat(~targets,i));
    [~,TFstats.pValMot(i)]=ttest2(motifGeneMat(targets,i),motifGeneMat(~targets,i));
    TFstats.motCorr(i)=corr(combBind(:,i),motifGeneMat(:,i),'rows','complete');
end
figure
selTFs=find(sum(TFstats{:,1:3},2)==3 & TFstats.hasMotif & TFstats.hasBH);
[~,idx]=sort(diff(log2(TFstats{selTFs,{'motNonTarget','motTarget'}}+0.01),1,2))
scatter(max(diff(log2(TFstats{selTFs(idx),{'motNonTarget','motTarget'}}+0.01),1,2),0),1:numel(idx),100,-log10(TFstats.pValMot(selTFs(idx))),'filled')
caxis([0 20])
colormap(gca,brighten(brewermap(128,'pubu'),0.3))
set(gca,'YTick',1:numel(idx),'YTickLabel',TFstats.Properties.RowNames(selTFs(idx)))
ylabel(colorbar(),'-log10 pValue')
saveas(gcf,'FigureS2.svg')
%% Figure 3
regTH=0.05;
TFstats.regType(:)=nan;
TFstats.regType(~isnan(TFstats.CrBHIdea))=(TFstats.CrBHIdea(~isnan(TFstats.CrBHIdea))>regTH) - (TFstats.CrBHIdea(~isnan(TFstats.CrBHIdea))<-regTH) 
TFstats.nBind=sum(combBind>1e4)'

temp=load('ChecProm.mat');
[~,idx]=ismember({'Gal11_TS','Cyc8_OL'},temp.barkaiSamples)
TFstats.CrCyc=corr(combBind,temp.sumProm(:,idx(2)),'rows','pairwise');
TFstats.CrMed=corr(combBind,temp.sumProm(:,idx(1)),'rows','pairwise')
clear temp idx

figure;

subplot(3,5,[1 2])
selTFs=find(~isnan(TFstats.CrBHIdea))
[~,idx]=sort(TFstats.CrBHIdea(selTFs))
hold off
scatter(1:numel(idx),TFstats.CrBHIdea(selTFs(idx)),[],TFstats.sgd(selTFs(idx)),'filled')
hold on
%colormap(gca,[brighten(brewermap(128,'greens'),0.3)])
%caxis([3 13])
xlim([0 numel(idx)]+0.5)
plot(xlim',[-.05 0 .05].*[1;1],'k-')
xlabel('TFs ordered')
ylabel('regScore')

subplot(3,5,6)
hold off
scatter(1+randn(sum(TFstats.regType==-1),1)/10,TFstats.nBind(TFstats.regType==-1),'filled')
hold on
scatter(1,median(TFstats.nBind(TFstats.regType==-1)),'filled')
scatter(2+randn(sum(TFstats.regType==1),1)/10,TFstats.nBind(TFstats.regType==1),'filled')
scatter(2,median(TFstats.nBind(TFstats.regType==1)),'filled')
ylabel('n Targets')

ohnologs=unique(sort([find(~isnan(GP.gene_infoR64.dubl)),GP.gene_infoR64.dubl(~isnan(GP.gene_infoR64.dubl))],2),'rows');
TFohno=ohnologs(all(ismember(ohnologs,GP.TFs),2),:)
[~,idxOhno]=ismember(TFohno,GP.TFs)
subplot(3,5,7)
hold off
scatter(TFstats.CrBHIdea(idxOhno(:,1)),TFstats.CrBHIdea(idxOhno(:,2)),100,TFstats.regType(idxOhno(:,1))~=TFstats.regType(idxOhno(:,2)),'filled')
text(TFstats.CrBHIdea(idxOhno(:,1)),TFstats.CrBHIdea(idxOhno(:,2)),TFstats.Properties.RowNames(idxOhno(:,1)))
hold on
plot([-0 0],ylim,'k-')
plot(xlim,[0 0],'k-')
xlabel('regScore-Paralog1')
ylabel('regScore-Paralog2')

subplot(3,4,9)
hold off
scatter(TFstats.CrCyc,TFstats.CrBHIdea,[],TFstats.cyc8,'filled')
selTFs=~isnan(TFstats.CrCyc)&~isnan(TFstats.CrBHIdea)
p=linortfit2(TFstats.CrCyc(selTFs),TFstats.CrBHIdea(selTFs))
hold on
plot(xlim,xlim.*p(1)+p(2))
caxis([-1 1].*.8)
xlabel('binding similarity Cyc8')
ylabel('regScore')
title(sprintf('%.2f',corr(TFstats.CrCyc(selTFs),TFstats.CrBHIdea(selTFs))))
colorbar()


subplot(3,4,10)
scatter(TFstats.CrMed,TFstats.CrBHIdea,[],log2(TFstats.actScore+10),'filled')
caxis([3 13])
colormap(gca,[brighten(brewermap(128,'greens'),0.3)])
selTFs=~isnan(TFstats.CrMed)&~isnan(TFstats.CrBHIdea)
p=linortfit2(TFstats.CrMed(selTFs),TFstats.CrBHIdea(selTFs))
hold on
title(sprintf('%.2f',corr(TFstats.CrMed(selTFs),TFstats.CrBHIdea(selTFs))))

colorbar()
xlabel('binding similarity Med15')
ylabel('regScore')
saveas(gcf,'figure3a.svg')
%% Figure S3
load('TFstats.mat')
load('borisColl.mat')
GP=load('group_impGitHub.mat'); % generelle informationen ueber das hefe genom

figure

subplot(1,3,1)
hold off
scatter(TFstats.CrCyc,TFstats.CrBHDel,[],TFstats.cyc8,'filled')
selTFs=~isnan(TFstats.CrCyc)&~isnan(TFstats.CrBHDel)
% p=linortfit2(TFstats.CrCyc(selTFs),TFstats.CrBHIdea(selTFs))
% hold on
% plot(xlim,xlim.*p(1)+p(2))
caxis([-1 1].*.8)
xlabel('binding similarity Cyc8')
ylabel('regScore (DEL)')
title(sprintf('%.2f',corr(TFstats.CrCyc(selTFs),TFstats.CrBHDel(selTFs))))
colorbar()


subplot(1,3,2)
scatter(TFstats.CrMed,TFstats.CrBHDel,[],log2(TFstats.actScore+10),'filled')
caxis([3 13])
colormap(gca,[brighten(brewermap(128,'greens'),0.3)])
selTFs=~isnan(TFstats.CrMed)&~isnan(TFstats.CrBHDel)
p=linortfit2(TFstats.CrMed(selTFs),TFstats.CrBHDel(selTFs))
hold on
xlabel('binding similarity Gal11')
ylabel('regScore (DEL)')
title(sprintf('%.2f',corr(TFstats.CrMed(selTFs),TFstats.CrBHDel(selTFs))))
colorbar()

subplot(1,3,3)
hold off
imagesc(corr(TFstats{:,{'CrBHAid','CrBHDel','CrBHIdea'}},TFstats{:,{'CrCyc','CrMed'}},'rows','pairwise'),[-.45 .45])
set(gca,'xtick',1:2,'xticklabel',{'CrCyc','CrMed'},'YTick',1:3,'YTickLabel',{'CrBHAid','CrBHDel','CrBHIdea'})
colormap(gca,brighten([flipud(brewermap(128,'Reds'));brewermap(128,'Blues')],0.3))
hold on
plot(xlim',[1 2]+[.5;.5],'-','Color',[ 1 1 1].*.75)
plot([1]+[.5;.5],ylim','-','Color',[ 1 1 1].*.75)
colorbar()
saveas(gcf,'FigureS3pA.svg')


%%
load('pughSP.mat')
load('coFactList.mat')
repressors=unique(coFact.Gene_Complex(coFact.act==-1&~ismember(coFact.bid,GP.TFs)));
pughRep=unique(pughNames(ismember(lower(pughNames),lower(repressors))));
clear meanRepPro
for i=1:numel(pughRep)
    meanRepPro(:,i)=mean(spPugh(:,ismember(pughNames,pughRep(i))),2,'omitnan');
end
repCorr=corr(combBind,meanRepPro,'rows','pairwise');

activators=unique(coFact.Gene_Complex(coFact.act==1&~ismember(coFact.bid,GP.TFs)&~ismember(coFact.Gene_Complex,repressors)));
pughAct=unique(pughNames(ismember(lower(pughNames),lower(activators))));
clear meanActPro
for i=1:numel(pughAct)
    meanActPro(:,i)=mean(spPugh(:,ismember(pughNames,pughAct(i))),2,'omitnan');
end
actCorr=corr(combBind,meanActPro,'rows','pairwise');

figure;
subplot(2,1,1)
hold off
imagesc(corr(TFstats{:,{'CrBHAid','CrBHDel','CrBHIdea'}},repCorr,'rows','pairwise'),[-.45 .45])
set(gca,'xtick',1:numel(pughRep),'xticklabel',upper(pughRep),'YTick',1:3,'YTickLabel',{'CrBHAid','CrBHDel','CrBHIdea'},'XTickLabelRotation',90)
colormap(gca,brighten([flipud(brewermap(128,'Reds'));brewermap(128,'Blues')],0.3))
hold on
plot(xlim',[1 2]+[.5;.5],'-','Color',[ 1 1 1].*.75)
plot([1:numel(pughRep)]+[.5;.5],ylim','-','Color',[ 1 1 1].*.75)

subplot(2,1,2)
hold off
imagesc(corr(TFstats{:,{'CrBHAid','CrBHDel','CrBHIdea'}},actCorr,'rows','pairwise'),[-.45 .45])
set(gca,'xtick',1:numel(pughAct),'xticklabel',upper(pughAct),'YTick',1:3,'YTickLabel',{'CrBHAid','CrBHDel','CrBHIdea'})
colormap(gca,brighten([flipud(brewermap(128,'Reds'));brewermap(128,'Blues')],0.3))
hold on
plot(xlim',[1 2]+[.5;.5],'-','Color',[ 1 1 1].*.75)
plot([1:numel(pughAct)]+[.5;.5],ylim','-','Color',[ 1 1 1].*.75)
saveas(gcf,'FigureS3pB.svg')

%% figure 4
clear all
load('TFstats.mat')
GP=load('group_impGitHub.mat'); % generelle informationen ueber das hefe genom
clearvars -except GP TFstats
load('seqOrd.mat','seqOrd')
temp=load('tfActivationDomains.mat')

[~,cyc8]=ismember('RPD3',GP.gene_infoR64.nameNew)
%load('AFactVsrep.mat')
load('AFoutDatasmall.mat','AFoutData')
AFoutData(cellfun('isempty',{AFoutData.ProtA}))=[];
AFoutData=struct2table(AFoutData)

[~,AFoutData.idA]=ismember(AFoutData.ProtA,GP.gene_infoR64.orf);
[~,AFoutData.idB]=ismember(AFoutData.ProtB,GP.gene_infoR64.orf);
AFoutData=AFoutData(AFoutData.idB==cyc8,:);

%% get DBDs
load('seqOrd.mat','seqOrd')
cisbpInfo=readtable('./generalYeast/motifs/TF_Information_all_motifs.txt','ReadVariableNames',false,'Delimiter','\t');
cisbpInfo=unique(cisbpInfo(:,[1,6]))
[~,cisbpInfo.id]=ismember(cisbpInfo.Var6,GP.gene_infoR64.orf)
[~,idx]=ismember(GP.TFs,cisbpInfo.id)
TFstats.cisbpID(idx>0)=cisbpInfo.Var1(idx(idx>0))
for i=find(~cellfun('isempty',TFstats.cisbpID))'
    cisbp=webread(sprintf('https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF=%s',TFstats.cisbpID{i}));
    cisbp=regexp(cisbp,'(?<=.h4.DNA Binding Domains).*?(?=\/table)','match','once');
    cisbp=regexp(cisbp,'(?<=\<textarea ).*?(?=\/textarea)','match');
    cisbp=regexp(cisbp,'[A-Z]{5,}','match','once');
    [dbdStart,dbdEnd]=regexp(seqOrd.newSeq{GP.TFs(i)},join(cisbp,'|'),'start','end');
    TFstats.dbd{i}=[dbdStart{1};dbdEnd{1}];
end
clearvars dbdStart dbdEnd cisbp cisbpInfo

repDomains=readtable('./TFanalysis-BorisMherHuett/RepDomains.xlsx','ReadVariableNames',false)
[~,repDomains.id]=ismember(upper(repDomains.Var1),GP.gene_infoR64.nameNew)
repDomains.num=str2double(split(repDomains.Var2,'–'));
[~,idx]=ismember(GP.TFs,repDomains.id);
TFstats.repDomain(idx>0,:)=repDomains.num(idx(idx>0),:);

AFoutData.nintA=cellfun('prodofsize',AFoutData.intA);
AFoutData.nintB=cellfun('prodofsize',AFoutData.intB);
[~,idx]=ismember(GP.TFs,AFoutData.idA);
TFstats{:,{'nIntA','nIntB'}}=NaN;
TFstats{idx>0,{'nIntA','nIntB'}}=AFoutData{idx(idx>0),{'nintA','nintB'}}
TFstats.intA(idx>0)=AFoutData.intA(idx(idx>0))
TFstats.disOrder=seqOrd.Var2(GP.TFs);
TFstats.seq=seqOrd.newSeq(GP.TFs);

TFstats.cyc8DO(:)=nan;
selTFs=cellfun('prodofsize',TFstats.intA)>0;
TFstats.cyc8DO(selTFs)=arrayfun(@(x)mean(TFstats.disOrder{x}(movmean(sparse(TFstats.intA{x},1,1,numel(TFstats.disOrder{x}),1),3)>0)),find(selTFs))
TFstats.cyc8Region(selTFs)=arrayfun(@(x)TFstats.seq{x}(movmean(sparse(TFstats.intA{x},1,1,numel(TFstats.disOrder{x}),1),3)>0),find(selTFs),'UniformOutput',false)
TFstats.cyc8comp(selTFs)=cellfun(@(x)accumarray(aa2int(x)',1,[24 1])',TFstats.cyc8Region(selTFs),'UniformOutput',false)

temp=load('tfActivationDomains.mat')
[~,idx]=ismember(GP.TFs,temp.tf.bid);
TFstats.actProfile(idx>0)=temp.tf.actProfile(idx(idx>0))
TFstats.actDO(:)=nan;
selTFs=cellfun('prodofsize',TFstats.actProfile)>0
TFstats.actDO(selTFs)=arrayfun(@(x)mean(TFstats.disOrder{x}(TFstats.actProfile{x}>10)),find(selTFs));
TFstats.actRegion(selTFs)=arrayfun(@(x)TFstats.seq{x}(TFstats.actProfile{x}>10),find(selTFs),'UniformOutput',false)
selTFs=cellfun('prodofsize',TFstats.actRegion)>0
TFstats.actcomp(selTFs)=cellfun(@(x)accumarray(aa2int(x)',1,[24 1])',TFstats.actRegion(selTFs),'UniformOutput',false)



figure
subplot(3,5,1)
hist=histogram(log2(TFstats.nIntB(TFstats.regType<0)+0.5),10,'normalization','probability')
hold on
histogram(log2(TFstats.nIntB(TFstats.regType>0)+0.5),hist.BinEdges,'normalization','probability')
xlabel('number of interacting AAs')
ylabel('fraction of activators/repressorsw')
xlim([-1.1 6])

subplot(3,5,[2 4])
selInts=ismember(AFoutData.idA,GP.TFs(TFstats.regType==-1))&AFoutData.nintA>0;
intMat=cell2mat(arrayfun(@(x)full(movmean(sparse(1,AFoutData.intB{x},1,1,numel(seqOrd.newSeq{cyc8})),5)>0),find(selInts),'UniformOutput',false))

tdrs=readtable("delAnalysis.xlsx",'Sheet','cyc8domains','ReadVariableNames',false)
tdrs.start=str2double(extractBefore(tdrs.Var2,'-'));
tdrs.stop=str2double(extractAfter(tdrs.Var2,'-'));

imagesc(sum(intMat))
ylabel(colorbar(),'npartners')
colormap(gca,brewermap(128,'reds'))
hold on
for i=1:height(tdrs)
    plot([tdrs{i,{'start','stop'}}],[1 1],'k-')
end

subplot(3,5,[6 11])
selReps=find(TFstats.regType==-1 & TFstats.nIntB>0)
[~,idx]=sort(cellfun('prodofsize',TFstats.disOrder(selReps)))
selReps=selReps(idx);
hold off
c=0
for i=selReps'
    c=c+1;
    disRescaled=rescale(TFstats.disOrder{i},c-0.4,c+0.4,'InputMax',1,'InputMin',0);
    plot(1:numel(TFstats.disOrder{i}),disRescaled,'Color',[1 1 1].*0.3)
    hold on
    scatter(TFstats.intA{i},disRescaled(TFstats.intA{i}),10,[1 0 0],'filled')
    plot(TFstats.dbd{i},repmat(c+0.45,height(TFstats.dbd{i}),width(TFstats.dbd{i})),'b-')
    plot(TFstats.repDomain(i,:),c+[0.3 .3],'r-')
end
ylim([0 numel(selReps)]+0.5)
set(gca,'Ytick',1:numel(selReps),'YTickLabel',TFstats.Properties.RowNames(selReps))


subplot(3,5,7)
hold off
scatter(1+randn(sum(TFstats.regType==-1),1)/10,TFstats.cyc8DO(TFstats.regType==-1),'filled')
hold on
scatter(1,median(TFstats.cyc8DO(TFstats.regType==-1),'omitnan'),'filled')

scatter(2+randn(sum(TFstats.regType==1),1)/10,TFstats.actDO(TFstats.regType==1),'filled')
scatter(2,median(TFstats.actDO(TFstats.regType==1),'omitnan'),'filled')

selTFs=cellfun('prodofsize',TFstats.cyc8Region)>0;
xlim([.5 2.5])

repComp=(cell2mat(TFstats.cyc8comp(~cellfun('isempty',TFstats.cyc8comp)&TFstats.regType==-1)));
actComp=(cell2mat(TFstats.actcomp(~cellfun('isempty',TFstats.actcomp)&TFstats.regType==1)));

idrSelComp=sum(cell2mat(arrayfun(@(x)accumarray([aa2int(seqOrd.newSeq{x}(seqOrd.Var2{x}>0.35 & seqOrd.Var2{x}<0.55))';24],1,[24,1])',find(~cellfun('isempty',seqOrd.Var2)&~cellfun('isempty',seqOrd.newSeq)),'UniformOutput',false)))

subplot(3,5,12)
hold off
scatter(idrSelComp./sum(idrSelComp),mean(repComp./sum(repComp,2)),100,'filled')
text(idrSelComp./sum(idrSelComp),mean(repComp./sum(repComp,2)),int2aa([1:24]'))
hold on
plot(xlim,xlim,'k--')
saveas(gcf,'figure4pA.svg')

figure
molviewer('AF-Cyc8-F1-model_v4.pdb')
evalrasmolscript(gcf,'set background white')
evalrasmolscript(gcf,'select 1-966;color cartoons [200,200,250]')
evalrasmolscript(gcf,'select 300-400;color cartoons [255,125,125]')
evalrasmolscript(gcf,'select 130-170;color cartoons [255,175,175]')

saveas(gcf,'figure4pB.svg')


%% Figure S4
clear all
load('TFstats.mat')
% t1=load('AFrepressors.mat');
% t2=load('AFactVsrep.mat')
% AFoutData=[t1.AFoutData,t2.AFoutData];
load('AFoutDatabig.mat','AFoutData')
GP=load('group_impGitHub.mat'); % generelle informationen ueber das hefe genom
[~,intReps]=ismember({'CYC8','TUP1','RPD3','SIN3','SAP30'},GP.gene_infoR64.nameNew)
AFoutData(cellfun('isempty',{AFoutData.ProtA}))=[];
AFoutData=struct2table(AFoutData)
[~,AFoutData.idA]=ismember(AFoutData.ProtA,GP.gene_infoR64.orf);
[~,AFoutData.idB]=ismember(AFoutData.ProtB,GP.gene_infoR64.orf);
AFoutData=AFoutData(ismember(AFoutData.idB,intReps)&ismember(AFoutData.idA,GP.TFs(TFstats.regType==-1)),:)
[protB,~,AFoutData.row]=unique(AFoutData.idB);
[protA,~,AFoutData.col]=unique(AFoutData.idA);
AFoutData=sortrows(AFoutData,{'idA','idB'});
intMat=full(sparse(AFoutData.row,AFoutData.col,cellfun('prodofsize',AFoutData.intB),numel(protB),numel(protA)))./full(sparse(AFoutData.row,AFoutData.col,1,numel(protB),numel(protA)));
protBOrder=[1,2,4,3,5];
[~,tfOrder]=sort(sum((intMat>0).*[4;2;0;1;0],1),'descend')
figure;
subplot(9,1,[1 6])
hold off
imagesc(log2(intMat(protBOrder,tfOrder)+0.5),[-1.1 6]);
colormap(gca,[.75 .75 .75;brighten(brewermap(128,'Purples'),0.3)])
hold on
plot(xlim,[1:numel(protB)]+[.5;.5],'k-')
plot([1:numel(protA)]+[.5;.5],ylim','k-')
set(gca,'ytick',1:numel(protB),'ytickLabel',GP.gene_infoR64.nameNew(protB(protBOrder)),'xtick',1:numel(protA),'xtickLabel',GP.gene_infoR64.nameNew(protA(tfOrder)),'XTickLabelRotation',90)
colorbar()
repList={'CTI6'	'DAL80'	'FKH1'	'GAL80'	'MIG1'	'MOT3'	'NRG1'	'OPI1'	'RDR1'	'ROX1'	'SKO1'	'UME6'	'URE2'	'XBP1'	'YHP1'	'YOX1'	'WHI5'};
[~,repID]=ismember(repList,GP.gene_infoR64.nameNew)
subplot(9,1,9)
hold off
imagesc(ismember(protA(tfOrder),repID)')
colorbar()
colormap(gca,brighten(brewermap(128,'reds'),0.5))
hold on
plot([1:numel(tfOrder)]+[.5;.5],ylim,'k-')

subplot(9,1,8)
hold off
imagesc(ismember(protA(tfOrder),GP.gene_infoR64.dubl)')
colorbar()
colormap(gca,brighten(brewermap(128,'reds'),0.5))

saveas(gcf,'FigureS4pA.svg')
function y = SaveFigPDF(fig,filename)

% この下の行のコメントアウトを解除するとより余白を削ったPDFを作成します。
% set(gca, 'LooseInset', get(gca, 'TightInset'));

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(fig,filename,'-dpdf','-r0')
end
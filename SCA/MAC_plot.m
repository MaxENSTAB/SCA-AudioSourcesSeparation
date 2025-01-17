function MAC = MAC_plot(x,y)

% xºÍyÎªÒª¼ÆËã²¢»æÖÆMAC¾ØÕóµÄÁ½¸öÕñÐÍ¾ØÕó,xÓëyÓ¦ÎªÍ¬ÐÍ¾ØÕó

flag  = size(x,2);
for i = 1:flag
    for j = 1:flag
        MAC(i,j) = (x(:,i)'*y(:,j))^2/((x(:,i)'*x(:,i))*(y(:,j)'*y(:,j)));
    end
end

set(gcf,'color','w')
h=bar3(MAC);
for n=1:numel(h)
    cdata=get(h(n),'zdata');
    cdata=repmat(max(cdata,[],2),1,4);
    set(h(n),'cdata',cdata,'facecolor','flat');
end
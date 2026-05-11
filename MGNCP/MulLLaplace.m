

function [Lapk,dma,ma]=MulLLaplace(ngmar)
Lapk={};
num=length(size(ngmar));
%     for i=1:3
%         Lapk{i}=zeros(size(ngmar,num),size(ngmar,num));
%     end
        

        ngtemp=tenmat(ngmar,num);
        ngtemp=double(ngtemp');
        T=1:size(ngtemp,2);
        sizea=size(ngtemp,2);
        T=sort(repmat(T',sizea,1));
        for j=1:sizea
            Ta=T((j-1)*sizea+1:j*sizea);

            Tatemp=ngtemp(:,Ta);
            ngsum1=Tatemp-ngtemp;
            ma(j,:)=double(exp(-vecnorm(ngsum1,2).^2));
            ma2(j,:)=ngtemp(:,j)'*ngtemp;
            ma3(j,:)=sum(min(ngtemp,Tatemp));
%             ma1(:,j)=double(exp(-vecnorm(ngsum1,2).^2));
        end

%         Lapk{1}=ones(size(Lapk{1}));    


        ma=ma-eye(size(ma));
        dma=diag(sum(ma));
        Lapk{2}=dma-ma;


%         ma2=ma2-eye(size(ma2));
        dma2=diag(sum(ma2));
        Lapk{3}=dma2-ma2;



%         ma3=ma3-eye(size(ma3));
        dma3=diag(sum(ma3));
        Lapk{4}=dma3-ma3;

        
%         clear ngs ngsum1 ngtemp 
        Lapk(cellfun(@isempty,Lapk))=[];
end







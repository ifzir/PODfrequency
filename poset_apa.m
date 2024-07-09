
global idx_trans
idx_trans=@(i,j,d)((i-1)*d+j);
function idr=trans_dr(d,dr)
    t=zeros(d);cnt=0;
    cnt=length(find(dr<80));
    for i=1:d-1
        for j=(i+1):d
            if dr(i)>80||dr(j)>80||(dr(i)==dr(j))
            t(i,j)=0;
            else if dr(i)>dr(j)
                t(i,j)=1;
                else
                t(i,j)=-1;
                end
            end
        end
    end
    idr=triu(t)-triu(t)';
    idr(1,1)=cnt;
end
function Y=randmtx(d,min_num)
    X=randperm(d);
    global M
    M=zeros(d);
    global visited
    item_num=min_num-1+randi(d-min_num+1);
    a_sort=randsrc(1,item_num-1,[1 -1]);
    for i=1:item_num-1
        M(X(i),X(i+1))=a_sort(i);
        M(X(i+1),X(i))=-a_sort(i);
    end
    M_O=M;
    for i=1:d
        for j=1:d
              if M(i,j)==1
                visited=logical(zeros(1,d));
                graphDFS(i,j,d,1);
              else if M(i,j)==-1
                   visited=logical(zeros(1,d));
                    graphDFS(i,j,d,-1);
                  end
            end
        end
    end
    Y={M,M_O};
end
function graphDFS(s,v,num,flag)
    global visited;
    visited(v)=true;
    global M;
    M(s,v)=flag;
    M(v,s)=-flag;
    for i=1:num
        if (M(v,i)==flag)&&(~visited(i))
            graphDFS(s,i,num,flag);
        end
    end
end
function B=prj(A,prj_list)
    n=size(A,2);
    A=squeeze(A);
    B=zeros(n);
    for i=1:n
        for j=1:n
            B(prj_list(i),prj_list(j))=A(i,j);
        end
    end
end
function i_prj=list_trans(prj_list)
    n=size(prj_list,2);
    i_prj=zeros(1,n);
    for i=1:n
        i_prj(prj_list(i))=i;
    end
end

function O2=step2ptb(O1_o,n,d,eps2,prj_list)
    O2=zeros(n,d,d);
    p=exp(eps2)./(exp(eps2)+1);
    for i=1:n
        r=rand;
        if r<p
            O2(i,:,:)=O1_o(i,:,:);
        else
            O2(i,:,:)=prj(O1_o(i,:,:),prj_list);
        end
    end
    O2_o=O2;
    O2={O2,O2_o};
end
function O2hat=step2clb(O2,n,d,eps2,i_prj)
    O2(find(O2<0))=0;
    agg=squeeze(sum(O2));
    tnij=zeros((d^2-d)./2,1);
    nij=zeros(d^2,1);
    p=exp(eps2)./(exp(eps2)+1);
    q=1-p;
    A=p.*eye(d^2);
    global idx_trans
    for i=1:d
        for j=1:d
            tnij(idx_trans(i,j))=agg(i,j)+agg(j,i);
            nij(idx_trans(i,j,d))=agg(i,j);
            A(idx_trans(i,j,d),idx_trans(i_prj(i),i_prj(j),d))=q;
        end
    end
    nij_hat=A\nij;
    O2hat=reshape(nij_hat,[d d])';
    O2hat(find(O2hat<0))=0;
end
function O1=step1ptb(O,n,d,eps1)
    O1=O;
    p=exp(eps1)./(exp(eps1)+1);
    for i=1:n
        r=rand;
        if r>p
            t=squeeze(O(i,:,:));
            k=find(triu(t,1));
            if length(k)==0
                O1(i,:,:)=O1(i-1,:,:);
            else
            r=randi(length(k));
            O1(i,k(r))=-O1(i,k(r));
            O1(i,:,:)=triu(squeeze(O1(i,:,:)))-triu(squeeze(O1(i,:,:)))';
            end
        end
    end
    O_op=O1;
    O1={O1,O_op};
end
function O1hat=step1clb1(O1_o,O1_all,n,d,eps1)
    O1_all(find(O1_all<0))=0;
    O1_o(find(O1_o<0))=0;
    counter=zeros(d,d);
    p=exp(eps1)./(exp(eps1)+1);
    O1hat=zeros(d,d);
    for i=1:n
         counter=counter+O1_all(i,:,:);
    end
    k=d-1;
    q=(1-p)./k;
        for i=1:d-1
            for j=(i+1):d
                answer=[p q;q p]\[counter(i,j);counter(j,i)];
                O1hat(i,j)=answer(1);O1hat(j,i)=answer(2);
            end
        end
     O1hat(find(O1hat<0))=0;
end

function O1hat=step1clb(O1_o,O1_all,n,d,eps1,cnt)
    O1_all(find(O1_all<0))=0;
    O1_o(find(O1_o<0))=0;
    counter=zeros(d,d,d);
    p=exp(eps1)./(exp(eps1)+1);
    O1hat=zeros(d,d,d);
    for i=1:n

        counter(cnt(i),:,:)=counter(cnt(i),:,:)+O1_all(cnt(i),:,:);
    end
    for k=2:d
        for i=1:d-1
            for j=(i+1):d
                q=(1-p)./(k-1);
                answer=[p q;q p]\[counter(k,i,j);counter(k,j,i)];
             O1hat(k,i,j)=answer(1);O1hat(k,j,i)=answer(2);
            end
        end
    end
    O1hat=squeeze(sum(O1hat,1));
end





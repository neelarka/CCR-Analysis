clear
reads = fastaread('C:\Users\eajoseph\Desktop\UNCG\Projects\Cas9\emulsion\bigcut.fasta')
ot = fastaread('C:\Users\eajoseph\Desktop\UNCG\Projects\Cas9\emulsion\ot.fasta')
scores=zeros(4,length(reads));
assign=0;
c1=0;
c2=0;
c3=0;
c4=0;

for x = 1:length(reads)
    s=reads(x).Sequence;
    scores(1,x)=swalign(s(44:100),ot(1).Sequence(100:170));
    scores(2,x)=swalign(s(44:100),ot(2).Sequence(100:170));
    scores(3,x)=swalign(s(44:100),ot(3).Sequence(100:170));
    scores(4,x)=swalign(s(44:100),ot(4).Sequence(100:170));
    sc=scores(:,x);
    if (sc(1)>110)&&(max([sc(2) sc(3) sc(4)]<110))
        assign(x)=1;
    elseif (sc(2)>110)&&(max([sc(1) sc(3) sc(4)]<110))
        assign(x)=2;
    elseif (sc(3)>110)&&(max([sc(1) sc(2) sc(4)]<110))
        assign(x)=3;
    elseif (sc(4)>110)&&(max([sc(1) sc(3) sc(2)]<110))
        assign(x)=4;
    else
        assign(x)=0;
    end
    if assign(x)~=0
   
        s1 = strfind(s,'GACTCACTATAG');
        s2 = strfind(s,'GAGTCCGAGCAGAAGAAGAA');
        if (numel(s1)>0)&&(numel(s2)>0)
        if (s2-s1)==20;
            if assign(x)==1
                c1=c1+1;
                hp1(c1,:)=s((s1+12):(s2-1));
            elseif assign(x)==2
                c2=c2+1;
                hp2(c2,:)=s((s1+12):(s2-1));
            elseif assign(x)==3
                c3=c3+1;
                hp3(c3,:)=s((s1+12):(s2-1));
                elseif assign(x)==4
                c4=c4+1;
                hp4(c4,:)=s((s1+12):(s2-1));
            end
        end
        end
    end
end
u1 = unique(hp1,'rows');
u2 = unique(hp2,'rows');
u3 = unique(hp3,'rows');
u4 = unique(hp4,'rows');
fid = fopen('C:\Users\eajoseph\Desktop\UNCG\Projects\Cas9\emulsion\hp1.txt','w');
for x = 1:size(u1,1)
n = num2str(sum(sum(hp1==u1(x,:),2)==8)); 
fprintf(fid,'%s\t%s\n',u1(x,:),n);
end
fclose(fid)
fid = fopen('C:\Users\eajoseph\Desktop\UNCG\Projects\Cas9\emulsion\hp2.txt','w');
for x = 1:size(u2,1)
n = num2str(sum(sum(hp2==u2(x,:),2)==8)); 
fprintf(fid,'%s\t%s\n',u2(x,:),n);
end
fclose(fid)
fid = fopen('C:\Users\eajoseph\Desktop\UNCG\Projects\Cas9\emulsion\hp3.txt','w');
for x = 1:size(u3,1)
n = num2str(sum(sum(hp3==u3(x,:),2)==8)); 
fprintf(fid,'%s\t%s\n',u3(x,:),n);
end
fclose(fid)
fid = fopen('C:\Users\eajoseph\Desktop\UNCG\Projects\Cas9\emulsion\hp4.txt','w');
for x = 1:size(u4,1)
n = num2str(sum(sum(hp4==u4(x,:),2)==8)); 
fprintf(fid,'%s\t%s\n',u2(x,:),n);
end
fclose(fid)
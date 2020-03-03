function check = check_DLMI(P,dP,R,G,A)
l=length(dP)/4;
check=1;
P0=P;
a=max(eig(P0-R))<0;
a1=min(eig(P0-G(1:4,1:4)))>0;
count=0;
a2=issymmetric(P0);
for i=2:l
    dP1=dP((i-1)*4+1:(i-1)*4+4,1:4);
    P1=dP1*1/l+P0;
    a5=issymmetric(P1);
    a3=min(eig(P1-G((i-1)*4+1:(i-1)*4+4,1:4)))>0;
    a4=max(eig(dP1+A'*P1+P1*A))<0;
    P0=P1;
    if (a==0 | a1==0 | a3==0 | a4==0 | a2==0 | a5==0)
        count=count+1;
        check=0;
        [a,a1,a2,a3,a4,a5]
        i;
    end
end
if check==1
    disp('System is FTS')
else
    disp('System is unstable')
    count
end
    
end
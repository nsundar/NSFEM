clear all
clc
Nodes=load('Node_COORDS.txt');
Elements=load('Elements.txt');
NNODES=size(Nodes,1);
COORDS=Nodes(:,2:3);
NCA=Elements(:,2:4);
NELM=size(NCA,1);
Nod_AdjEles=get_nod_adjele(NNODES,NELM,NCA);

B=cell(NNODES,1);

for inn=1:NNODES   % loop for all nodes   
    len=size(Nod_AdjEles{inn},2);
    for iadj=1:len  % loop for adjacent elements
       if iadj==1 % Add the nodes of first adjacent element
           inc=0;
           for k=1:3
           inc = inc+1;
           B{inn}(inc) = NCA(Nod_AdjEles{inn}(iadj),k);
           end
       else   
           for i=1:3
                flag=0;
                for j=1:inc
                	if B{inn}(j) == NCA(Nod_AdjEles{inn}(iadj),i)
                    flag=1;
                    break;
                    end
                end
                if(flag==1)
                inc=inc;
                else
                inc = inc+1;
                B{inn}(inc) = NCA(Nod_AdjEles{inn}(iadj),i);
                end
           end
        end
    end
end

i3=0;
i4=0;
i5=0;
i5=0;
i6=0;
i7=0;
i8=0;
i9=0;
i10=0;
i11=0;
i12=0;
i13=0;
U3(1) = 0; 
U4(1) = 0; 
U5(1) = 0; 
U6(1) = 0; 
U7(1) = 0; 
U8(1) = 0;
U9(1) = 0; 
U10(1) = 0;
U11(1) = 0;
U12(1) = 0;
U13(1) = 0;
        
for inn=1:NNODES
    Bnods_len=length(B{inn});
    switch Bnods_len
           case 3
           i3 = i3+1;
           U3(i3) = inn; 
           case 4
           i4 = i4+1;
           U4(i4) = inn; 
           case 5
           i5 = i5+1;
           U5(i5) = inn; 
           case 6
           i6 = i6+1;
           U6(i6) = inn; 
           case 7
           i7 = i7+1;
           U7(i7) = inn; 
           case 8
           i8 = i8+1;
           U8(i8) = inn; 
           case 9
           i9 = i9+1;
           U9(i9) = inn; 
           case 10
           i10 = i10+1;
           U10(i10) = inn; 
           case 11
           i11 = i11+1;
           U11(i11) = inn;
           case 12
           i12 = i12+1;
           U12(i12) = inn;
           case 13
           i13 = i13+1;
           U13(i13) = inn;
    end
end
if(U3(1)>=1)
    for i=1:length(U3)
        for j=1:4
            if j==1
                m3(i,j) = U3(i);
            else
                m3(i,j) = B{U3(i)}(j-1);
            end
        end
    end
    csvwrite('U1_Nodes.dat',m3);
else
end

if(U4(1)>=1) 
    for i=1:length(U4)
        for j=1:5
            if j==1
                m4(i,j) = U4(i);
            else
                m4(i,j) = B{U4(i)}(j-1);
            end
        end
    end
    csvwrite('U2_Nodes.dat',m4);
else
end

if(U5(1)>=1) 
    for i=1:length(U5)
        for j=1:6
            if j==1
                m5(i,j) = U5(i);
            else
                m5(i,j) = B{U5(i)}(j-1);
            end
        end
    end
    csvwrite('U3_Nodes.dat',m5);
else
end

if(U6(1)>=1) 
 for i=1:length(U6)
        for j=1:7
            if j==1
                m6(i,j) = U6(i);
            else
                m6(i,j) = B{U6(i)}(j-1);
            end
        end
 end
csvwrite('U4_Nodes.dat',m6);
else
end

if(U7(1)>=1) 
 for i=1:length(U7)
        for j=1:8
            if j==1
                m7(i,j) = U7(i);
            else
                m7(i,j) = B{U7(i)}(j-1);
            end
        end
 end
csvwrite('U5_Nodes.dat',m7);
else
end

if(U8(1)>=1) 
 for i=1:length(U8)
        for j=1:9
            if j==1
                m8(i,j) = U8(i);
            else
                m8(i,j) = B{U8(i)}(j-1);
            end
        end
 end
csvwrite('U6_Nodes.dat',m8);
else
end

if(U9(1)>=1) 
 for i=1:length(U9)
        for j=1:10
            if j==1
                m9(i,j) = U9(i);
            else
                m9(i,j) = B{U9(i)}(j-1);
            end
        end
 end
csvwrite('U7_Nodes.dat',m9);
else
end

if(U10(1)>=1) 
 for i=1:length(U10)
        for j=1:11
            if j==1
                m10(i,j) = U10(i);
            else
                m10(i,j) = B{U10(i)}(j-1);
            end
        end
 end
csvwrite('U8_Nodes.dat',m10);
else
end

if(U11(1)>=1) 
 for i=1:length(U11)
        for j=1:12
            if j==1
                m11(i,j) = U11(i);
            else
                m11(i,j) = B{U11(i)}(j-1);
            end
        end
 end
csvwrite('U9_Nodes.dat',m11);
else
end

if(U12(1)>=1) 
 for i=1:length(U12)
        for j=1:13
            if j==1
                m12(i,j) = U12(i);
            else
                m12(i,j) = B{U12(i)}(j-1);
            end
        end
 end
csvwrite('U10_Nodes.dat',m12);
else
end

if(U13(1)>=1) 
 for i=1:length(U13)
        for j=1:14
            if j==1
                m13(i,j) = U13(i);
            else
                m13(i,j) = B{U13(i)}(j-1);
            end
        end
 end
csvwrite('U11_Nodes.dat',m13);
else
end

% Filename: iteration.m
% Author: Yifei Li
% Queensland University of Technology, Brisbane, Australia, Jan 2021
% Reference:  Y. Li, S.T. Johnston, P.R. Buenzli, P. van Heijster, M.J. Simpson (2021) 
% Dimensionality affects extinction of bistable populations.
% The script contains:
%   - function 'iteration': iterate one time step in discrete simulations
%   - function 'judge_M': judge whether a movement event is successful
%   - function 'judge_P': judge whether a growth event is successful
%   - function 'Func_P': Allee kinetics
%   - function 'isoverBC': apply periodic boundary condition for all boundaries

function [count,indexA,A]=iteration(A,N,M,indexA,count,Pi,r,neighbour,AlleeParameter)
    j=1;
    randomchoose=randi(count,count,1);
    while j<count+1
            row=indexA(randomchoose(j),1);
            col=indexA(randomchoose(j),2);
            [newrow,newcol]=judge_M(A,row,col,N,M);
            if newrow>0
                A(row,col)=0;
                A(newrow,newcol)=1;
                indexA(randomchoose(j),1)=newrow;
                indexA(randomchoose(j),2)=newcol;
            end
        j=j+1;
    end
    j=1;
    count_fix=count;
    while j<count_fix+1
        if rand(1)<Pi
            randomchooseindex=randi(count);
            row=indexA(randomchooseindex,1);
            col=indexA(randomchooseindex,2);
            [newrow,newcol,occu]=judge_P(A,row,col,N,M,r,neighbour);
            Pii=Func_P(occu,AlleeParameter);
            %with death events
            if Pii>0
                if rand(1)<Pii
                    A(newrow,newcol)=1;
                    count=count+1;
                    indexA(count,1)=newrow;
                    indexA(count,2)=newcol;
                end
            elseif Pii<0
                if rand(1)<-Pii
                    A(row,col)=0;
                    indexA(randomchooseindex,:)=[];
                    count=count-1;
                end
            end
        end
        j=j+1;
    end
end
%movement mechanism
function [newrow,newcol]=judge_M(A,row,col,N,M)
    a=randi(6);
    if a<2
        newrow=row-1;
        newcol=col-1;
        if newrow<1
            newrow=N;
        end
        if newcol<1
            newcol=M;
        end
    elseif a<3
        newrow=row-1;
        newcol=col+1;
        if newrow<1
            newrow=N;
        end
        if newcol>M
            newcol=1;
        end
    elseif a<4
        newrow=row;
        newcol=col-2;
        if newcol<1
            newcol=M+newcol;
        end
    elseif a<5
        newrow=row;
        newcol=col+2;
        if newcol>M
            newcol=newcol-M;
        end
    elseif a<6
        newrow=row+1;
        newcol=col-1;
        if newrow>N
            newrow=1;
        end
        if newcol<1
            newcol=M;
        end
    else
        newrow=row+1;
        newcol=col+1;
        if newrow>N
            newrow=1;
        end
        if newcol>M
            newcol=1;
        end
    end
    if A(newrow,newcol)>0
        newrow=0;
    end
end
%proliferation mechanism
function [newrow,newcol,occu]=judge_P(A,row,col,N,M,r,neighbour)
    coordinate=zeros(neighbour,2);
    count=0;
    for i=row-r:row-1
        for j=col+row-2*r-i:2:col-row+2*r+i
            [newi,newj]=isoverBC(i,j,N,M);
            if A(newi,newj)<1
                count=count+1;
                coordinate(count,1)=newi;
                coordinate(count,2)=newj;
            end
        end
    end
    for j=col-2*r:2:col+2*r
        i=row;
        if j~=col
          [newi,newj]=isoverBC(i,j,N,M);
          if A(newi,newj)<1
              count=count+1;
              coordinate(count,1)=newi;
              coordinate(count,2)=newj;
          end
        end
    end
    for i=row+1:row+r
        for j=col-row-2*r+i:2:col+row+2*r-i
            [newi,newj]=isoverBC(i,j,N,M);
            if A(newi,newj)<1
                count=count+1;
                coordinate(count,1)=newi;
                coordinate(count,2)=newj;
            end
        end
    end
    occu=(neighbour-count)/(neighbour);
    if count>0
        selectedsite=randi(count);
        newrow=coordinate(selectedsite,1);
        newcol=coordinate(selectedsite,2);
    else
        newrow=0;
        newcol=0;
    end
end
%growth crowding function
function a=Func_P(C,AlleeParameter)
    %Allee effect
    a=2.5*(1-C)*(C-AlleeParameter);
end
%periodic boundary condition
function [newrow,newcol]=isoverBC(newrow,newcol,N,M)
    if newrow<1
        newrow=N+newrow;
    end
    if newrow>N
        newrow=newrow-N;
    end
    if newcol<1
        newcol=M+newcol;
    end
    if newcol>M
        newcol=newcol-M;
    end
end%Reflected boundary condition

function yout = poolData(yin,nVars,polyorder)
n = size(yin,1);

ind = 1;
% poly order 0
yout(:,ind) = ones(n,1);
ind = ind+1;

% poly order 1
for i=1:nVars
    yout(:,ind) = yin(:,i);
    ind = ind+1;
end

if(polyorder>=2)    % poly order 2
    for i=1:nVars
         for j=i:nVars
            if i == j
                yout(:,ind) = yin(:,i).*yin(:,j);
                ind = ind + 1;
            end
         end
    end
end

if(polyorder>=3)    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                if i == j && j == k
                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
                ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=4)    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l);
                    ind = ind+1;
                end
            end
        end
    end
end
if(polyorder>=5)    % poly order 5
    for i=1:nVars
            yout(:,ind) = sin(yin(:,i));
            ind = ind+1;
    end
end
if(polyorder>=6)    % poly order 6
    for i=1:nVars
        for j=1:nVars
            yout(:,ind) =yin(:,j).*sin(yin(:,i));
            ind = ind+1;
        end
    end
end
if(polyorder>=7)    % poly order 7
    for i=1:nVars
        for j=1:nVars
            for k=1:nVars
            yout(:,ind) =yin(:,j).*yin(:,k).*sin(yin(:,i));
            ind = ind+1;
            end
        end
    end
end
if(polyorder>=8)    % poly order 8
    for i=1:nVars
        for j=1:nVars
            for k=1:nVars
                for l=1:nVars
                yout(:,ind) =yin(:,j).*yin(:,k).*yin(:,l).*sin(yin(:,i));
                ind = ind+1;
                end
            end
        end
    end
end
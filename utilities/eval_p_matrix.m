function p_matrix = eval_p_matrix(nReach,nSpecies,pD_vec,downNode,V,depth,width)
%EVAL_P_MATRIX Summary of this function goes here
%   Detailed explanation goes here
p_matrix=zeros(nReach,nReach,nSpecies-1); % 3D matrix
            for sp=1:length(pD_vec)
                pD=pD_vec(sp);
                for i=1:nReach
                    j=downNode(i);
                    if j~=0
                        p_matrix(j,i,sp)=pD*V(i)/V(j); % ratio of volumes justifies that we work with abundances
                    else
                        p_matrix(i,i,sp)=pD; % quota of fish that do not exit from the outlet
                    end
                    k = find(downNode==i);
                    if ~isempty(k)
                        den=sum(width(k).*depth(k));
                        for indk=1:length(k)
                            kk=k(indk);
                            p_matrix(kk,i,sp)=(1-pD)*width(kk)*depth(kk)/den*V(i)/V(kk);
                        end
                    else
                        p_matrix(i,i,sp)=1-pD; % quota of fish that do not exit from a headwater
                    end
                end
            end
end


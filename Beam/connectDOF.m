function Td = connectDOF(data,Tn)


for e=1:size(Tn,1)
    for i=1:data.nne
        for j=1:data.ni
            Td(e,nod2dof(data.ni,i,j))=nod2dof(data.ni,Tn(e,i),j);
        end
    end
end
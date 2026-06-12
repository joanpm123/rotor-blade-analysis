function f = assemblyLoads(data,Td,Fe,Fel)

f=accumarray(reshape(Td',[],1),Fel(:),[data.ndof 1]) + Fe;

end
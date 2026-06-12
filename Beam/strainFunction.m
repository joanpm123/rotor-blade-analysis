function [ea,eb,es,et] = strainFunction(data,ue,Ba,Bb,Bs,Bt,R)

for e=1:data.nel
    ea(1,e) = Ba{e}*R(:,:,e)*ue(:,e);
    eb(:,e) = Bb{e}*R(:,:,e)*ue(:,e);
    es(:,e) = Bs{e}*R(:,:,e)*ue(:,e);
    et(1,e) = Bt{e}*R(:,:,e)*ue(:,e);
end
end
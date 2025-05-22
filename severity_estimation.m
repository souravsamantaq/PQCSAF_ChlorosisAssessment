function [ severity_score ] = severity_estimation(spx_pre )

[L]=size(spx_pre);

type_1=0;
type_2=0;
type_3=0;
type_4=0;

w1=1;
w2=2;
w3=3;
w4=4;

for i=1:L
    if(spx_pre(i)==1)
        type_1=type_1+1;
    end
    
     if(spx_pre(i)==2)
        type_2=type_2+1;
     end
    
      if(spx_pre(i)==3)
        type_3=type_3+1;
      end
    
       if(spx_pre(i)==4)
        type_4=type_4+1;
    end
end

total_type=type_1+type_2+type_3+type_4;

severity_score=w1*(type_1/total_type)+w2*(type_2/total_type)+w3*(type_3/total_type)+w4*(type_4/total_type);

end


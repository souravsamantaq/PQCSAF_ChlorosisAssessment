function [class_response] = classification_superpixel(spx_featureset)
[h1]=size(spx_featureset,1);
class_response=zeros(h1,1);

for i=1:h1
    
    x1=spx_featureset(i,:);
    [y1]=function_neuralnetwork(x1);
    
    [max_val,max_ind]=max(y1);
    
    class_response(i)=max_ind;
    
end

end
function [leaf_spx] = leaf_superpixel(l,Sp,im)
%==============================================
[h w]=size(l);
super_pix_l=zeros(1,600);
k=1; super_pix_l(k)=l(1,1);
for i=2:h
    for j=2:w
              x=l(i,j);
           for p=1:k
               if (super_pix_l(p)==x)
                  sw=1;
               end
           end
            if(sw==0)
                k=k+1;
                super_pix_l(k)=x;
            else
                 sw=0;
            end
             
    end
end


super_pix_no(1:k)=super_pix_l(1:k);
super_pix_im=zeros(h,w,3);
super_stats=zeros(4,k);

spxsum=zeros(1,k);


%----------------------------------------------------------------
v=1;
for i=1:k                                          
     
x=super_pix_no(i);
[h w]=size(im(:,:,1));
super_pix_im=zeros(h,w,3);

for i=1:h
    for j=1:w
        if (l(i,j)==x)
         super_pix_im(i,j,:)=im(i,j,:) ;
        end
    end
end
  

    super_pix_im=uint8(super_pix_im);
    %imshow((super_pix_im))
    spi=super_pix_im(:,:,2);
    max_pix=max(spi(:));
    %imshow((super_pix_im))
     if(max_pix > 0)
         sim=im2bw((rgb2gray(super_pix_im)));
        st=regionprops(sim, 'BoundingBox' );
        [x1]=size(st,1);
        if (x1==1) 
           csim=imcrop(super_pix_im,[st.BoundingBox(1) st.BoundingBox(2) st.BoundingBox(3) st.BoundingBox(4)]);
          %imshow(csim)
                [h1,w1]=size(csim(:,:,1));
                 if((h1*w1)>50)
                   %imshow(csim)
                    plant(v).superpixel=csim;                       %need to open
                    v=v+1;
                 end
        end
     end
     

% Sp(i).L
%    if( Sp(i).L~=0 )
%     %imshow(super_pix_im)
%     sim=im2bw((rgb2gray(super_pix_im)));
%    
%     st=regionprops(sim, 'BoundingBox' );
%    
%     [x1]=size(st,1)
%    
%     
%         if (x1==1) 
%           csim=imcrop(super_pix_im,[st.BoundingBox(1) st.BoundingBox(2) st.BoundingBox(3) st.BoundingBox(4)]);
%         end
%     end
%      
%      mp=csim;
%      %----------------------------------------------------------------
%      %if(Sp(i).L~=0 && Sp(i).L >20)
%      %mp=superpixel_selection(im,l,super_pix_no,i); %need to open
%      plant(v).superpixel=mp;                       %need to open
%      v=v+1;
%      cp=mp(:,:,2);                                 %need to open
%      %imshow(mp)
%      %end
%       % rr=(C(i).r);
%      %  cc=(C(i).c);
%      %  text(cc,rr,num2str(i),'horizontalAlignment','center','color',[1 1 1])
%      %im1=insertText(im,[rr cc],num2str(i),'FontSize',20,'textColor',[255 255 255],'BoxOpacity', 0.01);
%      %text(rr,cc,'k')
%      %entropyinfo(i)=sum(cp(:));
%      %plantsuperpixelfeature(i,:)=im_glcm(mp);
%      %dim(uint8(C(i).r),uint8(C(i).c))=mp;
%    %super_stats(3,i)=mp;%superpixel_selection(im,l,super_pix_no,i);
 end                                                %need to open




leaf_spx=plant;
v;
end





function [csim]=superpixel_selection(im,l,super_pix_no,xind) %spxv
x=super_pix_no(xind);
[h w]=size(im(:,:,1));
super_pix_im=zeros(h,w,3);

for i=1:h
    for j=1:w
        if (l(i,j)==x)
         super_pix_im(i,j,:)=im(i,j,:) ;
        end
    end
end

sim=im2bw((rgb2gray(super_pix_im)));
% figure(3)
imshow(uint8((sim)))
% %size(sim)
% %   if (max(sim(:))==0)
% %     csim=(sim);
% %   else    
    st=regionprops(sim, 'BoundingBox' );
%     st.BoundingBox(1)
%     st.BoundingBox(2)
%     st.BoundingBox(3)
%     st.BoundingBox(4)
%    
   size(st)
    csim=imcrop(super_pix_im,[st(1).BoundingBox(1) st(1).BoundingBox(2) st(1).BoundingBox(3) st(1).BoundingBox(4)]);
%  end
%
% csim=ones(45,45);
csim=uint8(csim);
%spxv=entropy(csim);
%figure(4)
%imshow(uint8(super_pix_im))
% st = regionprops(sim, 'BoundingBox' );
%  %rectangle('Position',[st.BoundingBox(1),st.BoundingBox(2),st.BoundingBox(3),st.BoundingBox(4)])
% csim=imcrop(super_pix_im,[st.BoundingBox(1) st.BoundingBox(2) st.BoundingBox(3) st.BoundingBox(4)]);
% size(csim)
% imshow(uint8((csim)))
%[out] = GLCM_Features1(csim,pairs)
end



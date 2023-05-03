function max_like_model = sort_like(sam, like_temp)

for i = 1:3
    
   temp = permute(sam(:,i,:), [3,1,2]);
   
   temp = [like_temp', temp];
   
   temp = sortrows(temp,'descend');
   
   sam_sort(:,i,:) = permute( temp(:,2:end),[2,3,1] );
    
    
end

max_like_model = sam_sort(:,:,1);

end

function model_final = model_final(cand,likelihood)

[mm,ind] = max(likelihood);

temp = cand(:,:,ind);

temp = temp(1:nnz(temp(:,1)),1:6);

model_final(:,1) = temp(:,1);

model_final(:,2:3) = temp(:,3:4);

end
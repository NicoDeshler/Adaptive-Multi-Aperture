function tp_out = zero_padding(tp, n_max)

    if size(tp,1) == n_max
       
        tp_out = tp;
        
    else
        
        tp_out = [tp; zeros(n_max-size(tp,1), size(tp,2), size(tp,3))];
        
    end
    

end
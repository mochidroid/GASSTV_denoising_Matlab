function result = Dlt_circulant(z)
    result = z(:,:,[end,1:end-1],:) - z(:,:,:,:);
end
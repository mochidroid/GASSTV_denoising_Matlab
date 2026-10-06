function result = Dv_Neumann(z)
    result = z([2:end, end], :, :) - z;
end
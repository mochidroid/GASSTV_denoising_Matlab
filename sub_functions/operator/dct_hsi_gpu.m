function[z] = dct_hsi_gpu(z,sizeof)

zr = reshape(z, sizeof(1)*sizeof(2), sizeof(3));
zr = dct(transpose(zr));
z = reshape(transpose(zr),sizeof);
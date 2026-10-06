from torchvision import transforms
from torch.utils.data import Dataset
from os import listdir, path
from PIL import Image
import torch
import torchvision.transforms.functional as TF
import scipy.io as scio
import numpy as np
import h5py
class MyToTensor(object):
    """Convert a ``PIL Image`` or ``numpy.ndarray`` to tensor.

    Converts a PIL Image or numpy.ndarray (H x W x C) in the range
    [0, 255] to a torch.FloatTensor of shape (C x H x W) in the range [0.0, 1.0]
    if the PIL Image belongs to one of the modes (L, LA, P, I, F, RGB, YCbCr, RGBA, CMYK, 1)
    or if the numpy.ndarray has dtype = np.uint8

    In the other cases, tensors are returned without scaling.
    """

    def __call__(self, pic):
        """
        Args:
            pic (PIL Image or numpy.ndarray): Image to be converted to tensor.

        Returns:
            Tensor: Converted image.
        """
        return TF.to_tensor(pic.copy())

    def __repr__(self):
        return self.__class__.__name__ + '()'



def load_mat_v73(file_name):
    """MATLAB v7.3形式のファイルを読み込む関数"""
    with h5py.File(file_name, 'r') as f:
        data = f['HSI_noisy'][:]  # 'DataCube'はファイル内の変数名に合わせてください
        data = np.array(data).astype(np.float32)
    return data

class Dataset(Dataset):
    def __init__(self, root_dirs, transform=None, verbose=False):
        self.root_dirs = root_dirs
        self.transform = transform
        self.images_path = []
        for cur_path in root_dirs:
            self.images_path += [path.join(cur_path, file) for file in listdir(cur_path) if file.endswith(('tif','png','jpg','jpeg','bmp','mat'))]
        self.verbose = verbose

    def __len__(self):
        return len(self.images_path)

    def __getitem__(self, idx):
        img_name = self.images_path[idx]
        # image = scio.loadmat(img_name)['HSI_noisy'].astype(np.float32)
        try:
            # MATLAB v7.3以外のファイル形式をscipyで読み込む
            image = scio.loadmat(img_name)['HSI_noisy'].astype(np.float32)
        except NotImplementedError:
            # MATLAB v7.3形式のファイルをh5pyで読み込む
            image = load_mat_v73(img_name)
        if self.transform:
            image = self.transform(image)
        if self.verbose:
            return image, img_name.split('/')[-1]
        return image
def get_gt(gt_path, img_name,verbose=False):
    tfs = []
    tfs += [
    MyToTensor()
    ]
    gt_transforms = transforms.Compose(tfs)
    image = scio.loadmat(gt_path+img_name)['HSI_clean'].astype(np.float32)
    image = gt_transforms(image)
    image=image/image.max()
    return image
def get_dataloaders(test_path_list,
                    drop_last=True, n_worker=0, verbose=False):
    batch_sizes = {'test':1}
    test_transforms = transforms.Compose([MyToTensor()])
    data_transforms = {'test': test_transforms}
    image_datasets = {'test': Dataset(test_path_list, data_transforms['test'], verbose=verbose)}
    dataloaders = {x: torch.utils.data.DataLoader(image_datasets[x], batch_size=batch_sizes[x],
                                                  num_workers=n_worker, drop_last=drop_last, shuffle=False)
                   for x in ['test']}
    return dataloaders

import torch
import torch.nn as nn
from ops.im2col import Cube2Col,Col2Cube
from collections import namedtuple
from ops.utils import soft_threshold,sparsity,kronecker,Init_DCT
from tqdm import tqdm
from pysptools.material_count import HySime
SMDSNetParams = namedtuple('SMDSNetParams', ['kernel_size', 'num_filters', 'stride', 'unfoldings','threshold', 'multi_lmbda','verbose'])
import  numpy as np
from skimage.restoration import  denoise_nl_means,estimate_sigma

def estimate_sigma_multichannel(img):
    """マルチチャネル（例: RGB画像）の場合、各チャネルごとに sigma を推定し、その平均を返す"""
    if img.ndim == 3:  # 画像が複数チャネルの場合 (H, W, C)
        sigmas = [estimate_sigma(img[..., i], average_sigmas=True) for i in range(img.shape[-1])]
        return np.mean(sigmas)
    else:  # シングルチャネルの場合
        return estimate_sigma(img, average_sigmas=True)

def denoise_nl_means_multichannel(img, patch_size, patch_distance, h, fast_mode, sigma):
    """マルチチャネル（例: RGB画像）の場合、各チャネルごとに denoise_nl_means を適用"""
    if img.ndim == 3:  # 画像が複数チャネルの場合 (H, W, C)
        denoised_channels = [
            denoise_nl_means(img[..., i], patch_size=patch_size, patch_distance=patch_distance, h=h, fast_mode=fast_mode, sigma=sigma)
            for i in range(img.shape[-1])
        ]
        return np.stack(denoised_channels, axis=-1)  # 各チャネルの結果を結合して戻す
    else:  # シングルチャネルの場合
        return denoise_nl_means(img, patch_size=patch_size, patch_distance=patch_distance, h=h, fast_mode=fast_mode, sigma=sigma)

class SMDSNet(nn.Module):
    def __init__(self, params: SMDSNetParams):
        super(SMDSNet, self).__init__()
        D=[]

        # # ここで num_filters の値と型を確認する
        # print(f'params.num_filters: {params.num_filters}, type: {type(params.num_filters)}')
        # print(f'params.num_filters[0]: {params.num_filters[0]}, type: {type(params.num_filters[0])}')

        # # 各フィルターごとに Init_DCT に整数を渡す
        # for i in range(3):  # num_filters の各要素に対応
        #     D.append(Init_DCT(params.kernel_size, params.num_filters[i]))

        # params.num_filters のコピーを取得
        num_filters = params.num_filters

        # もし num_filters がリストのリストになっている場合、最初のリストを使う
        if isinstance(num_filters[0], list):
            num_filters = num_filters[0]

        # デバッグ出力で確認
        print(f'num_filters: {num_filters}, type: {type(num_filters)}')

        for i in range(3):
            print(f'params.num_filters[{i}]: {num_filters[i]}, type: {type(num_filters[i])}')  # デバッグ出力
            # D.append(Init_DCT(params.kernel_size, 9))  # ここで整数を渡す
            D.append(Init_DCT(params.kernel_size, num_filters[i]))  # ここで整数を渡す


        # D.append(Init_DCT(params.kernel_size, params.num_filters[0]))
        # D.append(Init_DCT(params.kernel_size, params.num_filters[1]))
        # D.append(Init_DCT(params.kernel_size, params.num_filters[2]))
        A=[]
        B=[]
        W=[]
        self.apply_A=nn.ParameterList()
        self.apply_D=nn.ParameterList()
        self.apply_W=nn.ParameterList()
        Dic = kronecker(kronecker(D[0], D[1]), D[2])
        self.dom = torch.pinverse(Dic)
        for i in range(3):
            dtd = D[i].t() @ D[i]
            _, s, _ = dtd.svd()
            l = torch.max(s)
            D[i] /= torch.sqrt(l)
            A.append(D[i].transpose(0, 1))
            B.append(torch.clone(A[i].transpose(0, 1)))
            W.append(torch.clone(A[i].transpose(0, 1)))
            self.apply_A.append(nn.Parameter(A[-1]))
            self.apply_D.append(nn.Parameter(B[-1]))
            self.apply_W.append(nn.Parameter(W[-1]))
        self.params = params
        # total_filters=9*9*9
        total_filters=num_filters[0]*num_filters[1]*num_filters[2]
        # total_filters=params.num_filters[0]*params.num_filters[1]*params.num_filters[2]
        if params.multi_lmbda:
            self.lmbda = nn.ParameterList(
            [nn.Parameter(torch.zeros(1,1, 1,1,total_filters)) for _ in range(params.unfoldings)])
            for x in self.lmbda:
                nn.init.constant_(x, params.threshold)
            # [nn.init.constant_(x, params.threshold) for x in self.lmbda]
        else:
            self.lmbda = nn.Parameter(torch.zeros(1,1, 1,1,total_filters))
            nn.init.constant_(self.lmbda, params.threshold)

        self.soft_threshold = soft_threshold
    def forward(self, I):
        R, Ek, I_sub = self.pro_sub(I)
        output = self.denoise_sub(I_sub)
        bs=len(R)
        im=torch.Tensor([]).to(device=I.device)
        for i in range(bs):
            im_sub=output[i].permute([0, 2, 3, 1])
            _im=torch.matmul(im_sub, Ek[i].T)
            im=torch.cat((im,_im),0)
        output=im.permute([0,3,1,2])
        return output
    

    def pro_sub(self,I):
        hs = HySime()
        bands=I.shape[1]
        R=[]
        Ek=[]
        I_sub=[]
        sigma_est=0
        for _I in I:
            _I=_I.permute([1, 2, 0])
            sigma_est = estimate_sigma_multichannel(_I.cpu().numpy())  # ここを修正
            # sigma_est = estimate_sigma(_I.cpu().numpy(), multichannel=True, average_sigmas=True)
            # 各チャネルに対して denoise_nl_means を適用
            I_nlm = denoise_nl_means_multichannel(_I.cpu().numpy(), patch_size=7, patch_distance=9, h=0.08, fast_mode=True, sigma=sigma_est)
            _R, _Ek = hs.count(I_nlm)
            # I_nlm=denoise_nl_means(_I.cpu().numpy(), patch_size=7, patch_distance=9, h=0.08, multichannel=True,
            #                  fast_mode=True,sigma=sigma_est)
            _R, _Ek = hs.count(I_nlm)
            _Ek = torch.FloatTensor(_Ek).to(device=_I.device)
            if _R < self.params.kernel_size:
                _Ek = torch.cat((_Ek, torch.zeros([bands, self.params.kernel_size - _R],dtype=_Ek.dtype).to(device=_I.device)),1)
            _Ek=_Ek.to(device=_I.device)
            I_sub.append(torch.matmul(_I, _Ek).permute([2,0,1]))
            R.append(_R)
            Ek.append(_Ek)
        return R, Ek, I_sub
    def denoise_sub(self,I):
        params = self.params
        thresh_fn = self.soft_threshold
        bs = len(I)
        I_col=torch.Tensor([])
        batch_ind=[]
        I_size=[]
        for _I in I:
            padding = (params.stride - (_I.shape[2] - params.kernel_size) % params.stride) % params.stride
            im=Cube2Col(_I.unsqueeze(0), kernel_size=params.kernel_size, stride=params.stride, padding=padding, tensorized=True)
            I_size.append([im.shape[2]-1+params.kernel_size,_I.shape[1],_I.shape[2]])
            batch_ind.append(im.shape[2])
            I_col= torch.cat((I_col,im),2)
        I_col=I_col.to(_I.device)
        mean_patch = I_col.mean(dim=1, keepdim=True)
        I_col = I_col - mean_patch
        if I_col.is_cuda:
            self.apply_A = self.apply_A.cuda()
            self.apply_D = self.apply_D.cuda()
            self.apply_W = self.apply_W.cuda()
            self.dom=self.dom.cuda()
        kr_A = kronecker(kronecker(self.apply_A[0], self.apply_A[1]), self.apply_A[2])
        kr_D = kronecker(kronecker(self.apply_D[0], self.apply_D[1]), self.apply_D[2])
        kr_W = kronecker(kronecker(self.apply_W[0], self.apply_W[1]), self.apply_W[2])
        I_col = I_col.permute([0, 2, 3, 4, 1])
        gamma_k=thresh_fn(torch.matmul(I_col,self.dom.t()),0.001)
        num_unfoldings = params.unfoldings
        N = I_col.shape[1] * I_col.shape[2] * I_col.shape[0]
        for k in range(num_unfoldings):
            x_k = torch.matmul(gamma_k, kr_D.t())
            res = x_k - I_col
            r_k = torch.matmul(res, kr_A.t())
            lmbda_ = self.lmbda[k] if params.multi_lmbda else self.lmbda
            gamma_k = thresh_fn(gamma_k - r_k, lmbda_)
            if params.verbose:
                residual = 0.5 * (x_k - I_col).pow(2).sum() / N
                loss = residual + (lmbda_ * gamma_k).abs().sum() / N
                tqdm.write(
                    f'patch res {residual.item():.2e} | sparsity {sparsity(gamma_k):.2f} | loss {loss.item():0.4e} ')

        output_all = torch.matmul(gamma_k, kr_W.t())
        output_all = output_all.permute([0, 4, 1, 2, 3])
        output_all = output_all + mean_patch
        inds=np.cumsum(batch_ind)
        inds=np.hstack(([0],inds))
        output=[]
        for i in range(bs):
            output.append(Col2Cube(output_all[:,:,inds[i]:inds[i+1]], I_size[i], kernel_size=params.kernel_size, stride=params.stride, padding=0,
                              avg=True, input_tensorized=True))
        if params.verbose:
            tqdm.write('')
        return output


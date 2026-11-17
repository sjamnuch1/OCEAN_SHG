#!/usr/bin/env python3
"""Synthetic staging and optional compiled complex contraction tests."""
import argparse
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import numpy as np

HERE=Path(__file__).resolve().parents[1]


def run(exe=None):
    rng=np.random.default_rng(8293)
    def rand(shape):
        return rng.normal(size=shape)+1j*rng.normal(size=shape)
    for nb,ns in [(1,1),(4,1),(4,2)]:
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            nc,nv,nk,nvs,np_,l=3,2,2,4,3,1
            nm=2*l+1
            proj=rand((np_,nvs,nk,ns))
            mel=rand((np_,nm))
            cks=root/'parcksv.Ti0001'
            cks.write_bytes(struct.pack('=3i',np_,nvs*nk,ns)+proj.tobytes(order='F'))
            mf=root/'mels.test'
            mf.write_text(''.join(f'{z.real:.17g} {z.imag:.17g}\n' for z in mel.flatten(order='F')))
            xs=[rand((nc,nv,nk,nb)) for _ in range(2)]
            lines=['status XC_MINUS_COMPLETE','sector VALENCE_TDA','right_axis z',
                   'binary_layout raw_complex_DP_stream_Fortran_ncon_nval_nk_nbeta_nbw_no_padding',
                   'solution_units bohr_per_hartree',f'dimensions ncon nval nk nbeta nbw {nc} {nv} {nk} {nb} 1',
                   'brange 2 3 5 7','kmesh 1 1 2','q_reduced 0 0 .0001','eta_hartree .01','energy_count 2']
            for i,x in enumerate(xs,1):
                name=f'X{i}.bin'
                (root/name).write_bytes(x.tobytes(order='F'))
                lines.append(f'Gminus_omega z {i} {i}.0 {i*27}.0 1e-8 2e-8 1.0 4 {name}')
            manifest=root/'xc.info'
            manifest.write_text('\n'.join(lines)+'\n')
            out=root/'out'
            args=[sys.executable,str(HERE/'prepare_rbxc_cv.py'),'--manifest',str(manifest),
                  '--parcksv',str(cks),'--mels',str(mf),'--projector-first-valence','1',
                  '--znl','22','2','1','--element','Ti','--site','1','--axis','z','--out',str(out)]
            subprocess.run(args,check=True,capture_output=True,text=True)
            assert subprocess.run(args,capture_output=True).returncode!=0
            bad=args.copy(); bad[-1]=str(root/'bad'); bad[bad.index('--projector-first-valence')+1]='3'
            assert subprocess.run(bad,capture_output=True).returncode!=0
            # Build the full dense operator independently of Fortran's MATMUL loop.
            d=np.conj(np.einsum('pm,pvks->vkms',mel,proj[:,1:3,:,:]))
            A=np.zeros((nc*nk*4*nm,nc*nv*nk*nb),complex)
            for h in range(2):
                for m in range(nm):
                    for e in range(2):
                        if nb==1 and h!=e:
                            continue
                        beta=0 if nb==1 else 2*h+e
                        alpha=2*(h*nm+m)+e
                        for k in range(nk):
                            for c in range(nc):
                                for v in range(nv):
                                    oi=np.ravel_multi_index((c,k,alpha),(nc,nk,4*nm),order='F')
                                    ii=np.ravel_multi_index((c,v,k,beta),(nc,nv,nk,nb),order='F')
                                    A[oi,ii]=-d[v,k,m,min(h,ns-1)]
            probe=rand((A.shape[0],))
            flat=xs[0].flatten(order='F')
            assert np.allclose(np.vdot(probe,A@flat),np.vdot(A.conj().T@probe,flat))
            if exe:
                subprocess.run([str(exe),'cv.in'],cwd=out,check=True,capture_output=True,text=True)
                assert 'completion RBXC_CV_COMPLETE' in (out/'shg_rbxc_cv.info').read_text()
                for i,x in enumerate(xs,1):
                    actual=np.fromfile(out/f'shg_cache_RbXc_cv_zz_{i:06d}.bin',dtype=np.complex128)
                    np.testing.assert_allclose(actual,A@x.flatten(order='F'),rtol=2e-13,atol=2e-13)
                assert subprocess.run([str(exe),'cv.in'],cwd=out,capture_output=True).returncode!=0
    print('PASS: preparation, band slicing, overwrite rejection, dense-reference adjoint identity')
    print('PASS: compiled outputs for three spin layouts and two energies' if exe else
          'SKIP: compiled Fortran tests; supply --exe /path/shg_rbxc_cv.x')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--exe',type=Path)
    a=p.parse_args()
    run(a.exe.resolve(strict=True) if a.exe else None)

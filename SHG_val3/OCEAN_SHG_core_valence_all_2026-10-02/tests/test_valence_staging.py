#!/usr/bin/env python3
"""Synthetic staging checks; optional actual independent-particle solver test.

python3 tests/test_staging.py
python3 tests/test_staging.py --exe /path/to/shg_xc_minus_valence.x
python3 tests/test_staging.py --exe /path/to/shg_xc_minus_valence.x --mpi 2

Fixture bands: valence 1:2, conduction 3:6; two k points; one spin.
The optional solve has bubble/ladder/ALDA disabled, allowing exact comparison
with d/(omega+i*eta-(epsilon_n-epsilon_v)) on every retained component.
"""
import argparse
import json
import math
import shutil
import struct
import subprocess
import tempfile
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]
PREPARE = PACKAGE / 'prepare_valence.pl'
HA = 27.211386245988


def put_json(path, obj):
    path.write_text(json.dumps(obj))


def fixture(root):
    for p in ('Common', 'DFT', 'PREP/BSE', 'support/Common', 'OPF/zpawinfo'):
        (root / p).mkdir(parents=True)
    controls = dict(bande=True, bflag=False, lflag=False, aldaf=False,
                    bwflg=False, backf=False, semitda=False, qpflg=False,
                    broaden=0.3, strength=1, decut=10, disable_intraband=False,
                    epsilon_threshold=0, solver='gmres', write_rhs=False,
                    gmres=dict(nloop=12, gprc=10, ffff=1e-16, estyle='list',
                               elist=[20], erange=[], echamp=False),
                    gw=dict(control='none'),
                    haydock=dict(niter=100, converge=dict(thresh=-0.001,spacing=5)),
                    plot=dict(points=101, range=[0,30]))
    q = [0, 0, 1e-4]
    common = dict(calc=dict(mode='val', photon_q=q), screen=dict(mode='hll'),
                  computer=dict(para_prefix=''), structure=dict(metal=False),
                  bse=dict(val=controls, nbands=6,
                           occupation=dict(type='none',value=0)))
    b = dict(complete=True, xmesh=[2,2,2], kmesh=[1,1,2], kshift=[0,0,0],
             nspin=1, photon_q=q, brange=[1,2,3,6], nelec=4, fermi=0,
             cbm=HA*0.2, max_occ_band=2, min_unocc_band=3,
             cks=dict(enabled=True), edges=['1 1 0'],
             xred=[[0,0,0]], typat=[1], znucl=[4], elname=['Be'], epsilon=2,
             avecs=[[10,0,0],[0,10,0],[0,0,10]],
             bvecs=[[2*math.pi/10,0,0],[0,2*math.pi/10,0],[0,0,2*math.pi/10]])
    put_json(root/'Common/postDefaultsOceanDatafile', common)
    put_json(root/'DFT/dft.json',dict(bse=dict(complete=True)))
    put_json(root/'PREP/prep.json',dict(bse=b,dft_time=0,time=0))
    (root/'DFT/rhoofr').write_text('1\n1\n2 2 2\n')
    for name in ('Pquadrature','sphpts'):
        (root/'support/Common'/name).write_text('synthetic unused table\n')
    meljtv=root/'support/Common/meljtv.x'
    meljtv.write_text('#!/bin/sh\nprintf "1.0 0.0\\n" > mels\n')
    meljtv.chmod(0o755)
    for name in ('coreorbz004n01l00','prjfilez004','corezetaz004','radfilez004'):
        (root/'OPF/zpawinfo'/name).write_text('synthetic OPF input\n')
    pb=root/'PREP/BSE'
    # Native complex(real64), valence index fastest as in PREP M_vn.
    vals=[]
    for k in range(2):
        for n in range(4):
            for v in range(2):
                vals.extend([1e-4*(1+n+2*k+v),1e-4*(0.4+n-v)])
    (pb/'ptmels.dat').write_bytes(struct.pack('='+'d'*len(vals),*vals))
    (pb/'tmels.info').write_text('2 3 6 2\n')
    # Read by enk_selector=0, energies in Rydberg.
    (pb/'enkfile').write_text('-0.8 -0.4\n0.4 0.6 0.8 1.0\n'*2)
    # These wavefunction files are staged, but not read with bflag=lflag=aldaf=0.
    for name in ('con.u2.dat','val.u2.dat','wvfcninfo','wvfvainfo'):
        (pb/name).write_bytes(b'fixture-unused')
    for name in ('parcksv.Be0001','parcksc.Be0001'):
        (pb/name).write_bytes(b'synthetic-projector')
    return common,b


def run(root,target,extra=(),ok=True):
    cmd=['perl',str(PREPARE),'--root',str(root),'--target',str(root/target),
         '--axis','z','--ocean-source',str(root/'support'),*extra]
    r=subprocess.run(cmd,text=True,capture_output=True)
    assert (r.returncode==0)==ok, (cmd,r.stdout,r.stderr)
    return r


def complex_file(path):
    data=path.read_bytes()
    return [complex(r,i) for r,i in struct.iter_unpack('=dd',data)]


def solver_check(root,exe,mpi):
    target=root/'n2'
    cmd=([ 'mpirun','-n',str(mpi)] if mpi else [])+[str(Path(exe).resolve())]
    r=subprocess.run(cmd,cwd=target,text=True,capture_output=True,timeout=180)
    if r.returncode:
        raise AssertionError(r.stdout+'\n'+r.stderr)
    d=complex_file(target/'shg_cache_d_valence_z.bin')
    x=complex_file(target/'shg_cache_Gm_valence_w_z_000001.bin')
    expected_d=[]
    expected_x=[]
    q=2*math.pi/10*1e-4
    for k in range(2):
        for v,ev in enumerate([-0.4,-0.2]):
            for n,en in enumerate([0.2,0.3]):
                m=complex(1e-4*(1+n+2*k+v),1e-4*(0.4+n-v))
                dip=-1j*m.conjugate()/q
                expected_d.append(dip)
                expected_x.append(dip/(20/HA+0.3j/HA-(en-ev)))
    assert len(d)==len(x)==8
    for actual,expected in ((d,expected_d),(x,expected_x)):
        relative=math.sqrt(sum(abs(a-e)**2 for a,e in zip(actual,expected)) /
                           sum(abs(e)**2 for e in expected))
        assert relative < 2e-6, relative
    print('PASS: compiled solver source phase, normalization, layout and analytic resolvent')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--exe')
    ap.add_argument('--mpi',type=int,default=0)
    args=ap.parse_args()
    with tempfile.TemporaryDirectory(prefix='shg-val-test-') as temp:
        root=Path(temp)
        common,b=fixture(root)
        original=(root/'PREP/BSE/ptmels.dat').read_bytes()
        run(root,'dry',['--con-bands','3:4','--dry-run'])
        assert not (root/'dry').exists()
        run(root,'n2',['--con-bands','3:4'])
        run(root,'n4',['--con-bands','3:6','--energy-range','20 21 0.3'])
        # Side-by-side valence data must work without renaming or replacing
        # an existing core Common/DFT/PREP tree.
        shutil.copytree(root/'Common', root/'Common_VAL')
        shutil.copytree(root/'DFT', root/'DFT_VAL')
        shutil.copytree(root/'PREP', root/'PREP_VAL')
        run(root,'side_by_side',[
            '--common',str(root/'Common_VAL/postDefaultsOceanDatafile'),
            '--dft-dir',str(root/'DFT_VAL'),
            '--prep-dir',str(root/'PREP_VAL'),
            '--con-bands','3:4'])
        assert (root/'n2/brange.ipt').read_text()=='1 2\n3 4\n'
        assert (root/'n4/brange.ipt').read_text()=='1 2\n3 6\n'
        for name in ('n2','n4'):
            assert (root/name/'file_brange.ipt').read_text()=='1 2\n3 6\n'
            assert (root/name/'ptmels.dat').read_bytes()==original
            assert (root/name/'tmel_selector').read_text().strip()=='1'
            assert (root/name/'bwflg').read_text().strip()=='0'
            assert (root/name/'con.u2.dat').resolve()==root/'PREP/BSE/con.u2.dat'
            assert (root/name/'parcksv.Be0001').read_bytes()==b'synthetic-projector'
            assert all((root/name/f'mels.z004n01l00p{i:02d}').is_file() for i in (1,2,3))
        assert (root/'side_by_side/con.u2.dat').resolve()==root/'PREP_VAL/BSE/con.u2.dat'
        side_manifest=json.loads((root/'side_by_side/shg_valence_stage.json').read_text())
        assert Path(side_manifest['common_file'])==root/'Common_VAL/postDefaultsOceanDatafile'
        assert Path(side_manifest['dft_dir'])==root/'DFT_VAL'
        assert Path(side_manifest['prep_dir'])==root/'PREP_VAL'
        grid=(root/'n4/bse.in').read_text().splitlines()
        assert grid[4:]==['list','4','20','20.3','20.6','20.9'],grid
        assert (root/'PREP/BSE/ptmels.dat').read_bytes()==original
        run(root,'n2',ok=False)
        run(root,'too_many',['--con-bands','3:7'],ok=False)
        run(root,'missing_cbm',['--con-bands','4:6'],ok=False)
        run(root,'bad_step',['--energy-range','1 2 0'],ok=False)
        # q=0 remains fatal. An oblique q now stages with a prominent warning
        # and records that the source is directional rather than pure-z.
        b['photon_q']=[0,0,0]
        put_json(root/'PREP/prep.json',dict(bse=b,dft_time=0,time=0))
        run(root,'zero_q',ok=False)
        b['photon_q']=[1e-4,0,0]
        put_json(root/'PREP/prep.json',dict(bse=b,dft_time=0,time=0))
        common['calc']['photon_q']=[1e-4,0,0]
        put_json(root/'Common/postDefaultsOceanDatafile',common)
        result=run(root,'wrong_q')
        assert 'WARNING: stored PREP q is not parallel' in result.stderr
        oblique=json.loads((root/'wrong_q/shg_valence_stage.json').read_text())
        assert oblique['q_parallel_to_requested_axis'] is False
        assert abs(oblique['q_perpendicular_fraction']-1.0) < 1e-12
        print('PASS: staging two band windows, preserved file strides, energy grid, validation, no upstream writes')
        if args.exe:
            solver_check(root,args.exe,args.mpi)
        else:
            print('SKIP: compiled GMRES test (supply --exe and optionally --mpi)')


if __name__=='__main__':
    main()

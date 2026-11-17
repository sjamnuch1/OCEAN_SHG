#!/usr/bin/env python3
"""Prepare a serial R_CV X_V calculation from an existing valence manifest."""
import argparse
import json
import math
from pathlib import Path
import struct


def parse_manifest(path):
    fields, entries = {}, []
    for line in path.read_text().splitlines():
        w = line.split()
        if not w:
            continue
        if w[0] == 'Gminus_omega':
            if len(w) != 10:
                raise ValueError('Malformed Gminus_omega record')
            entries.append(w)
        else:
            if w[0] in fields:
                raise ValueError('Duplicate manifest field: '+w[0])
            fields[w[0]] = w[1:]
    if fields.get('status') != ['XC_MINUS_COMPLETE'] or fields.get('sector') != ['VALENCE_TDA']:
        raise ValueError('Requires a completed shg_xc_minus_valence manifest')
    layout = ['raw_complex_DP_stream_Fortran_ncon_nval_nk_nbeta_nbw_no_padding']
    if fields.get('binary_layout') != layout:
        raise ValueError('Unsupported input binary layout')
    if fields.get('solution_units') != ['bohr_per_hartree']:
        raise ValueError('Input must retain physical normalization')
    nc,nv,nk,nb,bw = map(int,fields['dimensions'][-5:])
    if min(nc,nv,nk) < 1 or nb not in (1,4) or bw != 1:
        raise ValueError('Unsupported dimensions/spin/TDA sector')
    if int(fields['energy_count'][0]) != len(entries) or not entries:
        raise ValueError('Incomplete energy list')
    if [int(w[2]) for w in entries] != list(range(1,len(entries)+1)):
        raise ValueError('Expected sequential energy indices')
    for w in entries:
        if w[1] != fields['right_axis'][0] or not all(math.isfinite(float(x)) for x in w[3:8]):
            raise ValueError('Invalid source axis or numerical diagnostics')
    return fields,entries,(nc,nv,nk,nb)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',required=True,type=Path)
    p.add_argument('--parcksv',type=Path,help='override auto-detected same-site valence OPF projector')
    p.add_argument('--mels',type=Path,help='override auto-detected length-dipole mels file')
    p.add_argument('--projector-first-valence',type=int,
                   help='override first global valence band stored in parcksv')
    p.add_argument('--znl',required=True,type=int,nargs=3,metavar=('Z','N','L'))
    p.add_argument('--element',required=True)
    p.add_argument('--site',required=True,type=int)
    p.add_argument('--axis',required=True,choices=['x','y','z'])
    p.add_argument('--out',required=True,type=Path,
                   help='calculation directory; may already exist')
    a=p.parse_args()
    manifest=a.manifest.resolve(strict=True)
    out=a.out.absolute()
    if out.exists():
        if not out.is_dir():
            raise ValueError('Output path exists but is not a directory')
    elif not out.parent.is_dir():
        raise ValueError('Output parent must exist')
    if not a.element.isalpha() or len(a.element)>2 or a.site<1:
        raise ValueError('Invalid element/site')
    z,n,l=a.znl
    if not (1<=z<=118 and n>=1 and 0<=l<n and l<=8):
        raise ValueError('Invalid ZNL')
    f,entries,(nc,nv,nk,nb)=parse_manifest(manifest)
    staged=manifest.parent
    cks=(a.parcksv if a.parcksv is not None
         else staged/f'parcksv.{a.element}{a.site:04d}').resolve(strict=True)
    photon_index={'x':1,'y':2,'z':3}[a.axis]
    auto_mels=staged/f'mels.z{z:03d}n{n:02d}l{l:02d}p{photon_index:02d}'
    mels=(a.mels if a.mels is not None else auto_mels).resolve(strict=True)
    if a.projector_first_valence is None:
        if 'file_brange' not in f or len(f['file_brange']) != 4:
            raise ValueError('Manifest lacks file_brange; supply --projector-first-valence')
        projector_first_valence=int(f['file_brange'][0])
    else:
        projector_first_valence=a.projector_first_valence
    if any('\n' in str(x) or len(str(x))>4000 for x in (manifest,cks,mels,out)):
        raise ValueError('Unsupported path')
    bands=list(map(int,f['brange']))
    mesh=list(map(int,f['kmesh']))
    if len(bands)!=4 or len(mesh)!=3 or min(mesh)<1 or math.prod(mesh)!=nk:
        raise ValueError('Invalid band/k metadata')
    if bands[1]-bands[0]+1!=nv or bands[3]-bands[2]+1!=nc:
        raise ValueError('Band bounds disagree with dimensions')
    with cks.open('rb') as inp:
        hdr=inp.read(12)
    if len(hdr)!=12:
        raise ValueError('Short parcksv header')
    np,nt,ns=struct.unpack('=3i',hdr)
    if min(np,nt)<1 or ns not in (1,2) or nt%nk:
        raise ValueError('Invalid parcksv header (requires native int32/complex128)')
    if nb==1 and ns!=1:
        raise ValueError('Scalar valence source requires one DFT spin channel')
    nvs=nt//nk
    skip=bands[0]-projector_first_valence
    if projector_first_valence<1 or skip<0 or skip+nv>nvs:
        raise ValueError('Projector file does not cover active valence bands')
    if cks.stat().st_size!=12+16*np*nt*ns:
        raise ValueError('Projector byte count mismatch')
    rows=[line.split() for line in mels.read_text().splitlines() if line.strip()]
    if len(rows)!=np*(2*l+1) or any(len(r)!=7 for r in rows):
        raise ValueError(
            'mels must have exactly nproj*(2*lcore+1) rows and 7 columns per row')
    if not all(math.isfinite(float(v.replace('D','E').replace('d','e')))
               for r in rows for v in r[:2]):
        raise ValueError('Nonfinite mels matrix element')
    try:
        for r in rows:
            tuple(int(v) for v in r[2:])
    except ValueError as exc:
        raise ValueError('mels columns 3-7 must be integer metadata') from exc
    q=list(map(float,f['q_reduced']))
    eta=float(f['eta_hartree'][0])
    if len(q)!=3 or not all(math.isfinite(v) for v in q+[eta]) or eta<=0:
        raise ValueError('Invalid eta/q')
    lines=[' '.join(map(str,[nc,nv,nk,nb,l,ns,nvs,skip,len(entries)])),
           a.axis+' '+f['right_axis'][0],
           ' '.join(map(str,bands+mesh+a.znl+[a.site,a.element])),
           ' '.join(format(v,'.17g') for v in [eta]+q),str(cks),str(mels)]
    for row in entries:
        src=Path(row[-1])
        if not src.is_absolute():
            src=manifest.parent/src
        src=src.resolve(strict=True)
        if src.stat().st_size!=16*nc*nv*nk*nb:
            raise ValueError('Incorrect X file size: '+str(src))
        if '\n' in str(src) or len(str(src))>4000:
            raise ValueError('Unsupported source path')
        filename=f'shg_cache_RbXc_cv_{a.axis}{row[1]}_{int(row[2]):06d}.bin'
        lines.extend([' '.join(row[3:5]),str(src),filename])
    out.mkdir(exist_ok=True)
    (out/'cv.in').write_text('\n'.join(lines)+'\n')
    (out/'provenance.json').write_text(json.dumps({
        'input_manifest':str(manifest),'input_metadata':f,
        'projectors':str(cks),'mels':str(mels),'projector_header':[np,nt,ns],
        'projector_first_valence':projector_first_valence,
        'auto_detected_projectors':a.parcksv is None,
        'auto_detected_mels':a.mels is None,
        'auto_detected_projector_first_valence':a.projector_first_valence is None,
        'znl':a.znl,'element':a.element,'site':a.site,'middle_axis':a.axis,
        'frame_requirement':'Projectors must share the input valence Bloch gauge; output electron frame is unchanged.',
        'not_verified_from_binary':'site identity, polarization, radial units, Bloch gauge and actual k coordinates'
    },indent=2)+'\n')
    print(f'Prepared {len(entries)} energies in {out}')
    print('Run there: /path/to/shg_rbxc_cv.x cv.in')


if __name__=='__main__':
    try:
        main()
    except (ValueError,KeyError,OSError,struct.error) as exc:
        raise SystemExit(str(exc))


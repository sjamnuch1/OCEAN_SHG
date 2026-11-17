#!/usr/bin/env python3
"""Mock executable integration test; does not test numerical PREP."""
import json
from pathlib import Path
import subprocess
import tempfile

SCRIPT=Path(__file__).resolve().parents[1]/'prep_val.pl'
with tempfile.TemporaryDirectory() as t:
    root=Path(t)
    common=root/'Common_VAL'
    dft=root/'DFT_VAL'
    binpath=root/'bin'
    for p in [common,dft,binpath,root/'PREP',root/'Common',root/'OPF/zpawinfo']:
        p.mkdir(parents=True)
    sentinel=root/'PREP'/'core_result'
    sentinel.write_text('KEEP CORE DATA\n')
    config={'calc':{'mode':'val','photon_q':[0,0,.001],'nonzero_q':True},
            'dft':{'program':'qe'},'computer':{'para_prefix':''},
            'bse':{'kmesh':[1,1,1],'kshift':[0,0,0],'xmesh':[2,2,2],'val_start':1}}
    data={'bse':{'complete':True,'hash':'fixture','split':True,'brange':[1,2,3,5],
                 'max_occ_band':2,'min_unocc_band':3,'kmesh':[1,1,1],
                 'kshift':[0,0,0],'photon_q':[0,0,.001],'nonzero_q':True},
          'scf':{'fermi':0,'nelec':4},'general':{'nspin':1},'time':0,
          'structure':{'avecs':[[10,0,0],[0,10,0],[0,0,10]],
                       'bvecs':[[.628,0,0],[0,.628,0],[0,0,.628]],
                       'typat':[1],'znucl':[14],'elname':['Si'],
                       'xred':[[0,0,0]],'epsilon':10}}
    cf=common/'postDefaultsOceanDatafile'
    cf.write_text(json.dumps(config))
    core_config={'calc':{'mode':'xas','edges':['1 1 0']}}
    (root/'Common/postDefaultsOceanDatafile').write_text(json.dumps(core_config))
    (root/'OPF/opf.json').write_text(json.dumps(
        {'completed':{'14':{'psp_hash':'fixture-psp','input_hash':'fixture-input'}}}))
    (dft/'dft.json').write_text(json.dumps(data))
    for shift in ['0.000000','0.999000']:
        folder=dft/f'k1_1_1q0.000000_0.000000_{shift}'
        save=folder/'Out'/'system.save'
        save.mkdir(parents=True)
        (save/'data-file-schema.xml').write_text('<fixture/>')
        (folder/'QE_EIGS.txt').write_text('fixture')
    for name in ['specpnt.5','Pquadrature','sphpts']:
        (binpath/name).write_text('fixture')
    executable=binpath/'ocean_prep.x'
    executable.write_text('''#!/usr/bin/env perl
use strict;
open my $t,'<','prep.tmels' or die $!;
die 'tmels disabled' unless <$t> =~ /true/;
open my $c,'<','prep.cks' or die $!;
die 'cks disabled' unless <$c> =~ /^1/;
open my $e,'<','edges' or die $!;
die 'wrong edge' unless <$e> =~ /^1 1 0/;
for my $f (qw(val.u2.dat con.u2.dat ptmels.dat tmels.info enkfile wvfcninfo wvfvainfo)) {
  open my $o,'>',$f or die $!; print $o "mock numerical output\\n";
}
for my $f (qw(parcksv.Si0001 parcksc.Si0001)) {
  open my $o,'>',$f or die $!; print $o "mock projector\\n";
}
open my $o,'>','cbm.out' or die $!; print $o "0.1\\n";
''')
    executable.chmod(0o755)
    args=['perl',str(SCRIPT),'--root',str(root),'--common',str(cf),
          '--dft-dir',str(dft),'--ocean-bin',str(binpath)]
    def call(*extra,ok=True):
        p=subprocess.run(args+list(extra),capture_output=True,text=True)
        if ok and p.returncode:
            raise RuntimeError(p.stdout+p.stderr)
        if not ok:
            assert p.returncode!=0,p.stdout
        return p
    call('--dry-run')
    assert not (root/'PREP_VAL').exists()
    original=cf.read_bytes()
    config['calc']['mode']='xas'
    cf.write_text(json.dumps(config))
    call('--dry-run',ok=False)
    config['calc']['mode']='val'
    config['calc']['photon_q']=[0,0,.002]
    cf.write_text(json.dumps(config))
    call('--dry-run',ok=False)
    cf.write_bytes(original)
    call()
    prep=root/'PREP_VAL'
    result=json.loads((prep/'prep.json').read_text())
    assert result['bse']['complete'] and result['bse']['tmels']
    assert result['bse']['cks']['enabled']
    assert result['bse']['edges']==['1 1 0']
    assert result['calc']['mode']=='val'
    assert (prep/'BSE'/'Out').resolve()==dft/'k1_1_1q0.000000_0.000000_0.999000'/'Out'
    assert (prep/'BSE'/'Out_shift').resolve()==dft/'k1_1_1q0.000000_0.000000_0.000000'/'Out'
    assert sentinel.read_text()=='KEEP CORE DATA\n'
    assert cf.read_bytes()==original
    assert json.loads((dft/'dft.json').read_text())==data
    call(ok=False)
print('PASS: VAL guard, edge import, simultaneous tmels/CKS, q validation, mock launch, core preservation')

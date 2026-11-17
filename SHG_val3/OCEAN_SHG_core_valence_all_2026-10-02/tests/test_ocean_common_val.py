#!/usr/bin/env python3
"""Mock integration test for the Common-only wrapper."""
import json
from pathlib import Path
import subprocess
import tempfile

SCRIPT=Path(__file__).resolve().parents[1]/'ocean_common_val.pl'

def write_helper(path,body):
    path.write_text('#!/usr/bin/env perl\nuse strict; use warnings;\n'+body)
    path.chmod(0o755)

with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp)
    ocean_bin=root/'bin'
    ocean_bin.mkdir()
    (ocean_bin/'oparse.json').write_text('{}')
    (ocean_bin/'oparse.type.json').write_text('{}')
    (ocean_bin/'Version').write_text('test\nhash\n')
    write_helper(ocean_bin/'parseJSON.pl',r'''
my $input=$ARGV[0]; open my $i,'<',$input or die $!;
my $mode=''; my @q=(0,0,0); my @k=(1,1,1);
while(<$i>) { $mode=$1 if /^\s*calc\.mode\s*\{\s*(\w+)/;
  @q=split ' ',$1 if /^\s*photon\.q\s*\{([^}]*)/;
  @k=split ' ',$1 if /^\s*bse\.kmesh\s*\{([^}]*)/; }
open my $o,'>','oceanDatafile' or die $!;
require JSON::PP;
print $o JSON::PP->new->encode({calc=>{mode=>$mode,photon_q=>[map {0+$_} @q]},
                                  bse=>{kmesh=>[map {0+$_} @k]}});
open my $c,'>','calc' or die $!; print $c "$mode\n";
open my $d,'>','dft' or die $!; print $d "qe\n";
''')
    write_helper(ocean_bin/'extractPsp.pl',r"""open my $o,'>','psp8'; print $o "fixture\n";
""")
    write_helper(ocean_bin/'psp_parser.pl',"exit 0;\n")
    write_helper(ocean_bin/'defaults.pl',r'''
open my $i,'<','oceanDatafile' or die $!; local $/; my $x=<$i>;
open my $o,'>','postDefaultsOceanDatafile' or die $!; print $o $x;
''')
    for helper in ['structure.pl','edges.pl','photon.pl']:
        write_helper(ocean_bin/helper,f"""open my $o,'>','{helper}.ran'; print $o "yes\\n";
""")
    val=root/'val.in'
    val.write_text('calc.mode { val }\nphoton.q { 0 0 .001 }\nbse.kmesh { 2 2 1 }\n')
    cmd=['perl',str(SCRIPT),'--root',str(root),'--input',str(val),'--ocean-bin',str(ocean_bin)]
    p=subprocess.run(cmd,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+p.stderr
    common=root/'Common_VAL'
    data=json.loads((common/'postDefaultsOceanDatafile').read_text())
    assert data['calc']['mode']=='val' and data['bse']['kmesh']==[2,2,1]
    assert (common/'common_val.complete').read_text().startswith('COMMON_VAL_COMPLETE')
    assert all((common/(x+'.ran')).is_file() for x in ['structure.pl','edges.pl','photon.pl'])
    assert not any((root/x).exists() for x in ['DFT','OPF','PREP','SCREEN','CNBSE','Common'])
    assert subprocess.run(cmd,capture_output=True).returncode!=0

    core=root/'core.in'
    core.write_text('calc.mode { xas }\nphoton.q { 0 0 0 }\nbse.kmesh { 1 1 1 }\n')
    bad=root/'BadCommon'
    p=subprocess.run(cmd[:4]+[str(core),'--ocean-bin',str(ocean_bin),'--output',str(bad)],capture_output=True,text=True)
    assert p.returncode!=0
    assert not (bad/'common_val.complete').exists()

print('PASS: Common-only execution, VAL guard, completion marker, overwrite protection, no later stages')

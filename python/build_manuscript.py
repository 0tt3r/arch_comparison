"""Build the main manuscript, or rebuild an archived note without replacing it."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

from project_paths import ROOT


def tectonic_environment():
    env=os.environ.copy()
    executable=shutil.which('tectonic')
    if executable:return executable,env
    # Local workshop environment: the conda package cache contains Tectonic,
    # although it is not installed into qb_gsee. Other machines use PATH.
    cache=Path.home()/'soft/mambaforge/pkgs'
    candidates=sorted(cache.glob('tectonic-*/bin/tectonic'))
    if not candidates:
        raise SystemExit('Install tectonic on PATH, or run latexmk -pdf from tex/.')
    from packaging.version import Version
    libraries=[]
    for name in ('icu','openssl','libstdcxx','libgcc','harfbuzz','libfreetype6','libglib',
                 'fontconfig','libpng','libzlib','libexpat','pcre2','graphite2','bzip2',
                 'brotli','libbrotlidec','libbrotlicommon','libiconv','libffi','libuuid'):
        options=[]
        for info in cache.glob(name+'-*/info/index.json'):
            row=json.loads(info.read_text())
            if row['name']==name:
                options.append((Version(row['version']),row.get('build_number',0),info.parent.parent/'lib'))
        if options:libraries.append(str(max(options)[2]))
    env['LD_LIBRARY_PATH']=':'.join(libraries+([env['LD_LIBRARY_PATH']] if env.get('LD_LIBRARY_PATH') else []))
    return str(candidates[-1]),env


def check_log(path):
    log=path.read_text()
    bad=[line for line in log.splitlines() if any(term in line for term in
         ('Overfull','undefined references','undefined citations','Citation `','Reference `','multiply defined'))]
    if bad:raise RuntimeError('TeX needs attention:\n'+'\n'.join(bad))


def build(archive=None):
    executable,env=tectonic_environment()
    if archive:
        if archive not in ('analytical_forms','nisq_refinement','refinement_v3'):
            raise ValueError('Unknown archived note')
        directory=ROOT/'build/tex/archive'/archive;directory.mkdir(parents=True,exist_ok=True)
        for path in (ROOT/'tex/archive').glob('*.bib'):shutil.copy2(path,directory/path.name)
        shutil.copy2(ROOT/f'tex/archive/{archive}.tex',directory/f'{archive}.tex')
        link=directory/'results'
        if not link.exists():link.symlink_to(ROOT/'results',target_is_directory=True)
        name=archive;out=directory
    else:
        directory=ROOT/'tex';name='main';out=ROOT/'build/tex/main';out.mkdir(parents=True,exist_ok=True)
    subprocess.run([executable,'--keep-logs','--keep-intermediates','--outdir',str(out),name+'.tex'],
                   cwd=directory,env=env,check=True)
    if not archive:
        # The supplement resolves main-text references through xr and build/tex/main/main.aux.
        sup=ROOT/'build/tex/supplement';sup.mkdir(parents=True,exist_ok=True)
        subprocess.run([executable,'--keep-logs','--keep-intermediates','--outdir',str(sup),'supplement.tex'],
                       cwd=directory,env=env,check=True)
        check_log(sup/'supplement.log')
        shutil.copy2(sup/'supplement.pdf',ROOT/'tex/supplement.pdf')
    check_log(out/(name+'.log'))
    result=out/(name+'.pdf')
    if not archive:
        shutil.copy2(result,ROOT/'tex/main.pdf');result=ROOT/'tex/main.pdf'
    print('Built',result.relative_to(ROOT))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive',choices=('analytical_forms','nisq_refinement','refinement_v3'))
    args=parser.parse_args();build(args.archive)

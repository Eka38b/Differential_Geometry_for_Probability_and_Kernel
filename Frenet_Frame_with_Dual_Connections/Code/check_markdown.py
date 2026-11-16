"""Check GitHub math conventions, local links, and LaTeX syntax in Markdown."""
import argparse,json,re,shutil,subprocess,tempfile
from pathlib import Path


def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=root/'Data/markdown_validation.json')
    args=ap.parse_args();errors=[];formulas=[];files=list(root.rglob('*.md'))
    for path in files:
        source=path.read_text()
        if '\\operatorname' in source:errors.append(f'{path}: rejected macro')
        body=re.sub(r'```[\s\S]*?```','',source)
        for href in re.findall(r'\]\(([^)]+)\)',body):
            if href.startswith(('https://','http://','#')):continue
            target=path.parent/href.split('#',1)[0]
            if target.resolve()!=args.out.resolve() and not target.exists():errors.append(f'{path}: missing link {href}')
        body=re.sub(r'`[^`\n]*`','',body)
        local=re.findall(r'^```math[ \t]*\n([\s\S]*?)^```[ \t]*$',source,re.MULTILINE)
        def display(match):
            local.append(match.group(1));return ''
        body=re.sub(r'\$\$([\s\S]*?)\$\$',display,body)
        local.extend(re.findall(r'(?<!\\)\$([^$\n]+?)(?<!\\)\$',body))
        for formula in local:
            if '\\begin{pmatrix}' in formula:errors.append(f'{path}: use the archive array matrix convention')
            if chr(92)*2 in formula:errors.append(f'{path}: use the archive explicit row separator')
            if re.search(r'(?<!\\)_(?!\{)',formula):errors.append(f'{path}: unbraced subscript {formula}')
            depth=0
            for k,char in enumerate(formula):
                if k and formula[k-1]=='\\':continue
                if char=='{':depth+=1
                if char=='}':depth-=1
                if depth<0:errors.append(f'{path}: unbalanced braces');break
            if depth:errors.append(f'{path}: unbalanced braces')
        formulas.extend(local)
    if errors:raise AssertionError('\n'.join(errors))
    compiler=shutil.which('pdflatex');compiled=False
    finder=shutil.which('kpsewhich')
    if compiler and finder:
        available=subprocess.run([finder,'pdflatex.fmt'],capture_output=True,text=True)
        if available.returncode or not available.stdout.strip():compiler=None
    elif compiler:
        compiler=None
    if compiler:
        with tempfile.TemporaryDirectory()as tmp:
            tmp=Path(tmp)
            tex='\\documentclass{article}\n\\usepackage{amsmath,amssymb}\n\\begin{document}\n'
            tex+='\n\n'.join('\\[\n'+formula+'\n\\]' for formula in formulas)
            tex+='\n\\end{document}\n'
            (tmp/'math.tex').write_text(tex)
            run=subprocess.run([compiler,'-halt-on-error','-interaction=nonstopmode','math.tex'],cwd=tmp,capture_output=True,text=True)
            if run.returncode:raise AssertionError(run.stdout[-5000:])
            compiled=True
    result=dict(date='2026-11-16',passed=True,markdown_files=len(files),math_expressions=len(formulas),
                forbidden_macro_found=False,unbraced_subscripts=0,broken_local_links=0,
                validation_kind='structural source checks including fenced math',local_latex_compilation=compiled,
                live_github_rendering_checked=False)
    args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()

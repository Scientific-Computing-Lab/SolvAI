"""Editable SVG primitives for the locally built SolvAI journal illustrations.

Molecular connectivity and conformers come from retained structure files. Solvent
shells and cavity outlines are explanatory schematics, not simulation results.
"""
from __future__ import annotations
from html import escape
from pathlib import Path
import io
import math
import re
import xml.etree.ElementTree as ET

import cairosvg
import cairocffi as cairo
import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / 'figures/source/fig1_assets'
INK = '#172D3A'
MID = '#657782'
GRID = '#DCE4E8'
BLUE = '#3F7FAD'
TEAL = '#078774'
AMBER = '#CB8C39'
PURPLE = '#9776AE'
ROSE = '#B44C7D'
COLORS = [BLUE, AMBER, '#599FB0', '#8EAAB7', '#557FBE', PURPLE]
COUNTS = [1, 5, 1, 1, 1, 6]


class Art:
    """720 drawing units correspond to a 180-mm figure width."""
    def __init__(self, height, title):
        self.height = height
        self.parts = []
        self.defs = []
        self.title_text = title
        self.nma = Chem.MolFromPDBFile(str(ASSETS / 'nma_openbabel.pdb'), removeHs=False)
        assert Chem.MolToSmiles(Chem.RemoveHs(self.nma)) == 'CNC(C)=O'
        self.conformers = [m for m in Chem.SDMolSupplier(
            str(ASSETS / 'dimethoxyethane_selected_conformers.sdf'), removeHs=False) if m]
        for symbol, pale, color, dark in [
            ('C','#CFD7DD','#5E6F7B','#243D4B'),
            ('N','#B7DAF6','#3C83B9','#19496F'),
            ('O','#F5C4B9','#D76D58','#923E31'),
            ('H','#FFFFFF','#EFF2F4','#B8C3CB'),
        ]:
            self.defs.append(f'<radialGradient id="atom{symbol}" cx="30%" cy="25%" r="75%">'
                f'<stop offset="0" stop-color="{pale}"/><stop offset=".45" stop-color="{color}"/>'
                f'<stop offset="1" stop-color="{dark}"/></radialGradient>')
        self.rect(0, 0, 720, height, 'white')

    def add(self, value):
        self.parts.append(value)

    def text(self, x, y, text, size=12, color=INK, weight=400, anchor='start', italic=False):
        # Baseline coordinates make optical alignment explicit.
        text = re.sub(r'(?<![\w-])-(?=\d)', '−', text)
        for i,line in enumerate(text.split('\n')):
            style = ' font-style="italic"' if italic else ''
            tx, line_anchor = x, anchor
            markup = escape(line)
            if 'mol⁻¹' in markup:
                before, after = markup.split('mol⁻¹', 1)
                # Relative baseline shifts create separate SVG text chunks.
                # Anchor the complete measured line once, then use start anchors
                # for all chunks so a superscript cannot escape the right edge.
                raw_before, raw_after = line.split('mol⁻¹', 1)
                context = cairo.Context(cairo.ImageSurface(cairo.FORMAT_ARGB32, 1, 1))
                context.select_font_face('Arial', cairo.FONT_SLANT_NORMAL,
                    cairo.FONT_WEIGHT_BOLD if weight >= 600 else cairo.FONT_WEIGHT_NORMAL)
                context.set_font_size(size)
                width = context.text_extents(raw_before+'mol')[4] + context.text_extents(raw_after)[4]
                context.set_font_size(size*.7)
                width += context.text_extents('−1')[4]
                tx -= width if anchor == 'end' else width/2 if anchor == 'middle' else 0
                line_anchor = 'start'
                markup = (before + 'mol' +
                    f'<tspan dy="{-size*.35}" font-size="{size*.7}">−1</tspan>' +
                    f'<tspan dy="{size*.35}">{after}</tspan>')
            self.add(f'<text x="{tx}" y="{y+i*size*1.28:.3f}" font-size="{size}" '
                f'fill="{color}" font-weight="{weight}" text-anchor="{line_anchor}"{style}>{markup}</text>')

    def title(self, x, y, letter, title):
        self.text(x,y,letter,size=20,weight=700)
        self.text(x+24,y-1,title,size=14,weight=700)

    def line(self,x1,y1,x2,y2,color=GRID,width=1,dash=None):
        extra = f' stroke-dasharray="{dash}"' if dash else ''
        self.add(f'<path d="M{x1},{y1} L{x2},{y2}" fill="none" stroke="{color}" stroke-width="{width}"{extra}/>')

    def path(self,d,color=INK,width=1,fill='none',opacity=1):
        self.add(f'<path d="{d}" fill="{fill}" stroke="{color}" stroke-width="{width}" opacity="{opacity}"/>')

    def rect(self,x,y,w,h,fill,stroke='none',radius=0,opacity=1):
        self.add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" '
                 f'fill="{fill}" stroke="{stroke}" opacity="{opacity}"/>')

    def circle(self,x,y,r,fill,stroke='none',width=1,opacity=1):
        self.add(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{r:.3f}" fill="{fill}" '
                 f'stroke="{stroke}" stroke-width="{width}" opacity="{opacity}"/>')

    def ellipse(self,x,y,rx,ry,fill,stroke='none',opacity=1):
        self.add(f'<ellipse cx="{x}" cy="{y}" rx="{rx}" ry="{ry}" fill="{fill}" stroke="{stroke}" opacity="{opacity}"/>')

    def arrow(self,x1,y1,x2,y2,color=INK,width=1.6):
        self.line(x1,y1,x2,y2,color,width)
        angle=math.atan2(y2-y1,x2-x1)
        tip=np.array([x2,y2]);back=tip-6*np.array([math.cos(angle),math.sin(angle)])
        normal=2.5*np.array([-math.sin(angle),math.cos(angle)])
        p,q=back+normal,back-normal
        self.path(f'M{p[0]},{p[1]} L{x2},{y2} L{q[0]},{q[1]} Z',color,0,fill=color)

    def lock(self,x,y,color=TEAL,s=1):
        self.path(f'M{x-3*s},{y} v{-3*s} a{3*s},{3*s} 0 0 1 {6*s},0 v{3*s}',color,1.2)
        self.rect(x-4*s,y,8*s,6*s,'white',color,1)
        self.circle(x,y+2.5*s,.7*s,color)

    def molecule(self,x,y,w,h,mol=None,angle=-25,hydrogens=True):
        m=mol if mol is not None else self.nma
        xyz=m.GetConformer().GetPositions().copy()
        rad=math.radians(angle)
        # Orthographic projection with fixed, modest tilt preserves connectivity.
        xyz=xyz @ np.array([[math.cos(rad),-math.sin(rad),0],
                           [math.sin(rad), math.cos(rad),0],[0,0,1]])
        xyz[:,1]=.92*xyz[:,1]+.38*xyz[:,2]
        keep=[a.GetIdx() for a in m.GetAtoms() if hydrogens or a.GetAtomicNum()>1]
        xy=xyz[:,:2];span=np.ptp(xy[keep],axis=0)
        scale=min(w/(span[0]+.8),h/(span[1]+.8))
        xy=(xy-(xy[keep].max(0)+xy[keep].min(0))/2)*scale
        xy[:,1]*=-1;xy+=np.array([x+w/2,y+h/2])
        items=[]
        for b in m.GetBonds():
            i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
            if i not in keep or j not in keep:continue
            p,q=xy[i],xy[j];v=q-p;n=np.array([-v[1],v[0]])/max(np.linalg.norm(v),1)
            offsets=[-.9,.9] if b.GetBondTypeAsDouble()==2 else [0]
            for off in offsets:
                p1,p2=p+off*n,q+off*n
                items.append(((xyz[i,2]+xyz[j,2])/2-.1,
                    f'<path d="M{p1[0]},{p1[1]} L{p2[0]},{p2[1]}" stroke="#80919B" stroke-width="{max(1.1,scale*.17)}" stroke-linecap="round"/>'))
        for i in keep:
            sym=m.GetAtomWithIdx(i).GetSymbol();r=scale*(.30 if sym=='H' else .39)
            if sym=='O':r*=1.07
            p=xy[i]
            items.append((xyz[i,2],f'<circle cx="{p[0]}" cy="{p[1]}" r="{r}" '
                f'fill="url(#atom{sym})" stroke="white" stroke-width=".35"/>'))
        for _,svg in sorted(items,key=lambda v:v[0]):self.add(svg)

    def water(self,x,y,scale=1,angle=0,opacity=1):
        self.add(f'<g transform="translate({x},{y}) rotate({angle}) scale({scale})" opacity="{opacity}">')
        self.line(0,0,-5,4,'#9AA8B1',2);self.line(0,0,5,4,'#9AA8B1',2)
        self.circle(-5,4,2.5,'url(#atomH)','#C8D0D5',.3)
        self.circle(5,4,2.5,'url(#atomH)','#C8D0D5',.3)
        self.circle(0,0,3.3,'url(#atomO)','white',.3)
        self.add('</g>')

    def shell(self,x,y,w,h,explicit=True):
        cx,cy=x+w/2,y+h/2
        if not explicit:
            for i,opacity in [(2,.07),(1,.10),(0,.16)]:
                ww=w/2+i*3;hh=h/2+i*3
                d=f'M{cx-ww*.87},{cy-hh*.18} C{cx-ww},{cy-hh*.9} {cx-ww*.35},{cy-hh} {cx},{cy-hh*.7} C{cx+ww*.6},{cy-hh*1.1} {cx+ww},{cy-hh*.35} {cx+ww*.82},{cy+hh*.2} C{cx+ww*.85},{cy+hh*.85} {cx+ww*.2},{cy+hh*.83} {cx-ww*.08},{cy+hh*.66} C{cx-ww*.7},{cy+hh} {cx-ww},{cy+hh*.5} {cx-ww*.87},{cy-hh*.18}Z'
                self.path(d,BLUE,.6,BLUE,opacity)
        self.molecule(x+w*.22,y+h*.18,w*.56,h*.63)
        if explicit:
            for j,angle in enumerate([-165,-117,-63,-14,38,90,143]):
                a=math.radians(angle)
                self.water(cx+w*.45*math.cos(a),cy+h*.45*math.sin(a),scale=w/95,angle=angle-80,opacity=.92)

    def responses(self,x,y,width=250,height=12,counts=COUNTS):
        total=sum(counts);gap=3;groupgap=7
        cell=(width-gap*(total-len(counts))-groupgap*(len(counts)-1))/total
        cur=x
        for n,color in zip(counts,COLORS):
            if not n:continue
            for j in range(n):
                self.rect(cur,y,cell,height,color,radius=.6);cur+=cell+gap
            cur+=groupgap-gap

    def feature_stripes(self,x,y,width=80,height=13):
        rng=np.random.default_rng(42)
        for i in range(40):
            self.rect(x+i*width/40,y,width/44,height,INK if rng.uniform()>.45 else GRID)

    def ai_model(self,x,y,w=90,h=36,color=TEAL):
        """Architecture-neutral learned mapping; no tree or neural topology implied."""
        self.rect(x+w*.18,y,w*.64,h,'#FFFFFF',color,4)
        for fraction in (.25,.5,.75):
            self.line(x,y+h*fraction,x+w*.18,y+h*fraction,color,1.4)
            self.circle(x,y+h*fraction,1.5,color)
        self.line(x+w*.82,y+h*.5,x+w,y+h*.5,color,1.4)
        self.circle(x+w,y+h*.5,1.5,color)
        self.text(x+w*.5,y+h*.5+5,'AI',16,color,600,'middle')

    def trees(self,x,y,w=90,h=36,color=TEAL):
        # Shared node coordinates keep each edge attached to the actual branch
        # nodes, with four distinct leaves rather than near-duplicate middle tips.
        stride=w/3
        for k in range(3):
            cx=x+(k+.5)*stride
            nodes=[(cx,y),(cx-.18*stride,y+.37*h),(cx+.18*stride,y+.37*h),
                   (cx-.30*stride,y+.74*h),(cx-.07*stride,y+.74*h),
                   (cx+.07*stride,y+.74*h),(cx+.30*stride,y+.74*h)]
            for parent,child in [(0,1),(0,2),(1,3),(1,4),(2,5),(2,6)]:
                self.line(*nodes[parent],*nodes[child],color,.9)
            for xx,yy in nodes:self.circle(xx,yy,1.5,color)

    def embed_plot(self,fig,x,y,w,h,prefix):
        from matplotlib.text import Text
        from journal_style import polish_figure
        polish_figure(fig)
        for item in fig.findobj(match=Text):
            text = item.get_text().replace('mol⁻¹', 'mol$^{-1}$')
            item.set_text(re.sub(r'(?<![\w{}$-])-(?=\d)', '−', text))
        buff=io.StringIO();fig.savefig(buff,format='svg',transparent=True)
        source=buff.getvalue();root=ET.fromstring(source)
        view=root.attrib['viewBox'].split();sw,sh=float(view[2]),float(view[3])
        inner=''.join(ET.tostring(child,encoding='unicode') for child in root)
        ids=re.findall(r'\bid="([^"]+)"',inner)
        for ident in ids:
            inner=inner.replace(f'id="{ident}"',f'id="{prefix}_{ident}"')
            inner=inner.replace(f'#{ident})',f'#{prefix}_{ident})').replace(f'="#{ident}"',f'="#{prefix}_{ident}"')
        self.add(f'<g transform="translate({x},{y}) scale({w/sw},{h/sh})">{inner}</g>')

    def svg(self):
        return (f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="180mm" height="{self.height/4}mm" viewBox="0 0 720 {self.height}">'
            f'<title>{escape(self.title_text)}</title><defs>'+''.join(self.defs)+'</defs>'
            '<g font-family="Arial, Helvetica, sans-serif" stroke-linecap="round" stroke-linejoin="round">'
            +''.join(self.parts)+'</g></svg>')

    def save(self,path,preview=None):
        path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
        source='\n'.join(line.rstrip() for line in self.svg().splitlines());path.with_suffix('.svg').write_text(source+'\n')
        cairosvg.svg2pdf(bytestring=source.encode(),write_to=str(path.with_suffix('.pdf')))
        if preview:
            preview=Path(preview);preview.parent.mkdir(parents=True,exist_ok=True)
            cairosvg.svg2png(bytestring=source.encode(),write_to=str(preview),output_width=1800)
        return path.with_suffix('.pdf')

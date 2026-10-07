"""Compose checked Azure motifs with native scientific labels and connectors."""
from pathlib import Path
import base64
import html
import argparse
import cairosvg

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'paper/supplementary/figures'
ASSET=OUT/'endpoint_components'
ASSET.mkdir(exist_ok=True)
parser=argparse.ArgumentParser()
parser.add_argument('--preview-dir',type=Path)
QA=parser.parse_args().preview_dir
if QA:QA.mkdir(parents=True,exist_ok=True)
INK='#18303C';BLUE='#417FA8';TEAL='#008C7A';AMBER='#C78935';GRAY='#647581'

# Selected label-free raster motifs are supplied in endpoint_components/.
# Every scientific label and connector is authored below as native SVG.

class SVG:
    def __init__(self,h=840):
        self.markers=set()
        self.items=[f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="720" height="{h}" viewBox="0 0 720 {h}">',
          '<defs><marker id="arrow" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0,0 L6,3 L0,6 Z" fill="context-stroke"/></marker></defs>',
          f'<rect width="720" height="{h}" fill="#FFFFFF"/>']
    def text(self,x,y,t,size=12.5,color=INK,bold=False,anchor='start'):
        self.items.append(f'<text x="{x}" y="{y}" fill="{color}" font-family="Arial" font-size="{size}" font-weight="{700 if bold else 400}" text-anchor="{anchor}">{html.escape(t)}</text>')
    def rect(self,x,y,w,h,stroke=BLUE,fill='white',r=5):
        self.items.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="1.4"/>')
    def path(self,d,color=GRAY,arrow=False,dash=False):
        marker='arrow_'+color.lstrip('#')
        if arrow and marker not in self.markers:
            self.items.append(f'<defs><marker id="{marker}" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0,0 L6,3 L0,6 Z" fill="{color}"/></marker></defs>')
            self.markers.add(marker)
        self.items.append(f'<path d="{d}" stroke="{color}" stroke-width="1.6" fill="none"'+(f' marker-end="url(#{marker})"' if arrow else '')+(' stroke-dasharray="4 3"' if dash else '')+'/>')
    def arrow(self,x,y,xx,yy,color=BLUE,dash=False):self.path(f'M{x},{y} L{xx},{yy}',color,True,dash)
    def img(self,key,x,y,w,h):
        b=base64.b64encode((ASSET/(key+'.png')).read_bytes()).decode()
        self.items.append(f'<image x="{x}" y="{y}" width="{w}" height="{h}" preserveAspectRatio="xMidYMid meet" xlink:href="data:image/png;base64,{b}"/>')
    def box(self,x,y,w,h,title,color=BLUE):
        self.rect(x,y,w,h,color,r=3);self.text(x+w/2,y+h/2+4,title,12.5,color,True,'middle')
    def module(self,x,y,w,h,title,color=TEAL):
        path=ROOT/'paper/figures/journal/teaser_components/endpoint.png'
        b=base64.b64encode(path.read_bytes()).decode()
        self.items.append(f'<image x="{x}" y="{y}" width="{w}" height="{h}" preserveAspectRatio="xMidYMid meet" xlink:href="data:image/png;base64,{b}"/>')
        self.text(x+w/2,y+h+18,title,12,color,True,'middle')
    def lock(self,x,y,color=BLUE):
        self.rect(x-4,y,8,8,color,r=1)
        self.path(f'M{x-3},{y} v-4 a3,3 0 0 1 6,0 v4',color)

    def plus(self,x,y):
        self.items.append(f'<circle cx="{x}" cy="{y}" r="13" fill="white" stroke="{TEAL}" stroke-width="1.4"/>');self.text(x,y+5,'+',20,TEAL,False,'middle')
    def row(self,y,letter,title,question,evidence,limit):
        if y:self.path(f'M10,{y} L710,{y}','#D6E0E5')
        self.text(12,y+26,letter,21,INK,True);self.text(39,y+25,title,17,INK,True)
        self.text(39,y+48,question,12,GRAY)
        self.text(39,y+246,evidence,12,TEAL,True);self.text(39,y+266,limit,11.5,GRAY)
    def finish(self):return '\n'.join(self.items+['</svg>'])

def diagram(c,key,y):
    Y=lambda v:y+v
    if key=='tree':
        c.text(28,Y(138),'Features',13,BLUE,True)
        c.path(f'M91,{Y(134)} L125,{Y(134)} L125,{Y(77)} L462,{Y(77)}',BLUE)
        for j in range(3):
            x=162+j*119;c.arrow(x+43,Y(77),x+43,Y(88),BLUE)
            c.img('tree',x,Y(86),108,77)
            c.path(f'M{x+42},{Y(157)} L{x+42},{Y(192)} L573,{Y(192)}',TEAL)
        c.text(350,Y(217),'One selected leaf value per tree',12,GRAY,False,'middle')
        c.box(574,Y(170),74,44,'Mean',TEAL);c.arrow(651,Y(192),681,Y(192),TEAL)
        c.text(690,Y(198),'ŷ',21,TEAL,True)
    elif key=='resnet':
        c.text(21,Y(153),'Features',13,BLUE,True);c.arrow(85,Y(149),126,Y(149))
        for x in [138,358]:
            c.module(x,Y(128),118,45,'Learned correction',BLUE)
            c.arrow(x+120,Y(149),x+151,Y(149),BLUE);c.plus(x+168,Y(149))
            c.path(f'M{x-11},{Y(149)} L{x-11},{Y(90)} L{x+168},{Y(90)} L{x+168},{Y(133)}',TEAL,True)
            c.text(x+72,Y(80),'Identity path',12,TEAL,False,'middle')
            c.arrow(x+183,Y(149),x+211,Y(149),TEAL)
        c.box(578,Y(126),105,45,'Readout',TEAL)
    elif key=='dual':
        for yy,name,col in [(104,'Structure',BLUE),(174,'Responses',AMBER)]:
            c.text(22,Y(yy+4),name,13,col,True);c.arrow(107,Y(yy),154,Y(yy),col)
            c.module(162,Y(yy-24),125,45,'Branch representation',col)
            c.path(f'M294,{Y(yy)} L392,{Y(yy)} L392,{Y(140)}',col)
        c.text(391,Y(83),'Concatenate',12,TEAL,False,'middle')
        c.arrow(394,Y(140),439,Y(140),TEAL)
        c.module(449,Y(114),135,51,'Joint head',TEAL)
        c.arrow(591,Y(140),662,Y(140),TEAL);c.text(677,Y(146),'ŷ',22,TEAL)
    elif key=='ridge_residual':
        c.text(23,Y(102),'Continuous',12,BLUE,True);c.text(23,Y(119),'descriptors',12,BLUE,True)
        c.arrow(115,Y(111),172,Y(111),BLUE);c.box(180,Y(90),163,43,'Ridge baseline',BLUE)
        c.text(23,Y(183),'All features',12,TEAL,True);c.arrow(114,Y(179),172,Y(179),TEAL)
        c.module(180,Y(154),163,49,'Neural residual',TEAL)
        c.path(f'M349,{Y(111)} L532,{Y(111)} L532,{Y(131)}',BLUE,True)
        c.path(f'M349,{Y(179)} L532,{Y(179)} L532,{Y(161)}',TEAL,True)
        c.plus(532,Y(146));c.arrow(550,Y(146),666,Y(146),TEAL);c.text(681,Y(152),'ŷ',22,TEAL)
        c.text(386,Y(75),'Training target − baseline',12,AMBER,False,'middle')
        c.path(f'M387,{Y(82)} L387,{Y(187)} L350,{Y(187)}',AMBER,True,True)
    elif key=='graph':
        c.box(22,Y(118),78,46,'Graph',BLUE);c.arrow(106,Y(141),147,Y(141),BLUE)
        c.img('graph',156,Y(74),63,125);c.text(187,Y(219),'Atom states',12,BLUE,False,'middle')
        c.arrow(227,Y(140),278,Y(140),TEAL)
        c.box(287,Y(117),161,46,'Sum | mean | count',TEAL)
        c.text(366,Y(96),'Concatenate',12,TEAL,False,'middle')
        c.arrow(454,Y(140),506,Y(140),TEAL)
        c.module(515,Y(116),146,49,'Neural head',TEAL)
        c.text(464,Y(221),'RDKit + response branch',12,AMBER,False,'middle')
        c.path(f'M463,{Y(206)} L501,{Y(206)} L501,{Y(161)} L516,{Y(161)}',AMBER,True)
        c.arrow(665,Y(141),693,Y(141),TEAL)
    elif key in ['molformer_deterministic','molformer_fixed_tree']:
        c.text(22,Y(127),'SMILES',13,BLUE,True);c.arrow(85,Y(123),119,Y(123),BLUE)
        c.module(129,Y(100),170,51,'Frozen MoLFormer',BLUE);c.lock(280,Y(98))
        c.arrow(306,Y(123),356,Y(123),BLUE)
        c.path(f'M363,{Y(111)} L363,{Y(193)} M296,{Y(188)} L363,{Y(188)}',TEAL)
        c.text(257,Y(213),'Structure + 15 responses',12,AMBER,False,'middle')
        c.text(370,Y(85),'Concatenate',12,TEAL,False,'middle')
        if key=='molformer_deterministic':
            c.arrow(365,Y(150),420,Y(150),TEAL)
            c.module(432,Y(124),148,53,'Neural head',TEAL)
            c.arrow(588,Y(150),674,Y(150),TEAL)
        else:
            c.path(f'M365,{Y(150)} L412,{Y(150)} L412,{Y(91)} L486,{Y(91)} L486,{Y(103)}',TEAL,True)
            c.img('tree',431,Y(92),115,85)
            c.path(f'M476,{Y(172)} L476,{Y(198)} L568,{Y(198)} L568,{Y(151)}',TEAL,True)
            c.box(558,Y(112),110,36,'Leaf mean',TEAL)
            c.arrow(670,Y(132),696,Y(132),TEAL)
            c.text(532,Y(80),'ExtraTrees',12,TEAL,True,'middle')
    elif key=='tabm':
        c.text(22,Y(146),'Features',13,BLUE,True);c.arrow(87,Y(142),129,Y(142),BLUE)
        c.rect(128,Y(70),302,150,BLUE,r=3)
        c.img('tabm',140,Y(80),87,123);c.text(184,Y(235),'Shared weights',12,BLUE,True,'middle')
        for j in range(3):
            yy=98+j*43;c.arrow(234,Y(yy),279,Y(yy),BLUE)
            c.box(286,Y(yy-14),127,29,'Modulate rᵢ, sᵢ',AMBER)
            c.path(f'M418,{Y(yy)} L494,{Y(yy)} L494,{Y(142)}',TEAL)
        c.text(278,Y(65),'Within-layer modulation',12,AMBER,False,'middle')
        c.text(499,Y(216),'32 predictions',12,TEAL,True,'middle')
        c.arrow(497,Y(142),545,Y(142),TEAL);c.box(552,Y(120),74,43,'Mean',TEAL)
        c.arrow(632,Y(142),678,Y(142),TEAL);c.text(689,Y(148),'ŷ',21,TEAL)
    elif key=='tabpfn':
        c.text(22,Y(80),'Labelled context',13,BLUE,True)
        for row in range(3):
            for col in range(5):c.rect(24+col*29,Y(91+row*22),26,19,AMBER if col==4 else BLUE,'white',1)
        c.text(153,Y(82),'y',12,AMBER,True,'middle')
        c.text(22,Y(179),'Unlabelled query',12,TEAL,True)
        for col in range(5):c.rect(24+col*29,Y(190),26,20,AMBER if col==4 else BLUE,'white',1)
        c.text(154,Y(205),'?',14,AMBER)
        c.path(f'M175,{Y(123)} L225,{Y(123)} L225,{Y(140)} L282,{Y(140)}',BLUE,True)
        c.path(f'M175,{Y(201)} L238,{Y(201)} L238,{Y(154)} L282,{Y(154)}',TEAL,True)
        c.module(292,Y(117),169,62,'Pretrained TabPFN',BLUE);c.lock(443,Y(115))
        c.arrow(468,Y(148),510,Y(148),TEAL)
        c.box(518,Y(121),181,51,'Distribution → mean',TEAL)
        c.text(602,Y(199),'No gradient fine-tuning',12,GRAY,False,'middle')


ROWS=[
('tree','ExtraTrees','Local partitions: can averaging nearby labels predict hydration?',
 'ARROW MAE 0.202; strong original matched-response controls.', 'Finite leaf-value range limits extrapolation of the mean.'),
('resnet','Residual MLP','Learn nonlinear corrections while preserving an identity path.',
 'Tests a direct neural replacement; ARROW MAE 0.340.', 'Neural flexibility alone does not recover the tree accuracy.'),
('dual','Dual structure / response branches','Let each input group learn its own representation before combining them.',
 'Size-held-out MAE 1.560 versus 1.772 for refitted trees.', 'Separate branches do not impose physical additivity.'),
('ridge_residual','Ridge + neural residual','Fit a continuous baseline; learn the remaining target with a network.',
 'Size-held-out MAE 1.391; no tree-leaf ceiling in the neural correction.', 'A linear trend is not a validated asymptotic free-energy law.'),
('graph','Size-sensitive message passing','Aggregate learned atom states rather than only recording fragment presence.',
 'Responses lower external MAE from 1.311 to 1.168.', 'Shown: sum / mean readout; atomic-plus-global variant also tested.'),
('molformer_deterministic','Frozen MoLFormer + neural head','Reuse a pretrained sequence representation alongside response information.',
 'Responses help in all three cohorts; external MAE 1.111.', 'The encoder is frozen; the endpoint is trained on hydration labels.'),
('tabm','TabM: parameter-sharing ensemble','Obtain several neural predictors while sharing most of their weights.',
 'External MAE 1.132; size-held-out MAE 1.389.', '32 members share parameters; they are not 32 independent fits.'),
('tabpfn','TabPFN-3.5: prediction from context','Condition a pretrained transformer on eligible labelled training rows.',
 'ARROW MAE 0.199; size-held-out MAE 1.318.', 'Native intervals under-cover; input recipes differ across partitions.'),
('molformer_fixed_tree','Frozen MoLFormer + ExtraTrees','Change the representation while retaining the original endpoint family.',
 'ARROW MAE 0.234 with responses; 0.348 without.', 'This control retains the bounded tree endpoint.')]

for group in range(3):
    c=SVG()
    for j,row in enumerate(ROWS[group*3:group*3+3]):
        key,title,question,evidence,limit=row;y=j*280
        c.row(y,'abc'[j],title,question,evidence,limit);diagram(c,key,y)
        solo=SVG(280);solo.row(0,'',title,question,evidence,limit);diagram(solo,key,0)
        (ASSET/(key+'.svg')).write_text(solo.finish())
    stem=f'Supp_Fig{8+group}_endpoint_mechanisms'
    svg=c.finish();(OUT/(stem+'.svg')).write_text(svg)
    cairosvg.svg2pdf(bytestring=svg.encode(),write_to=str(OUT/(stem+'.pdf')))
    if QA:cairosvg.svg2png(bytestring=svg.encode(),write_to=str(QA/(stem+'.png')),scale=1.5)
print('Nine checked conceptual composites; native text/connectors, Azure raster motifs, white backgrounds.')

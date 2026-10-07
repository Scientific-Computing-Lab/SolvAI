"""Native, port-anchored mechanism plates: what each endpoint changes."""
from pathlib import Path
import argparse
from journal_vector_art import Art, ROOT, INK, MID, GRID, BLUE, TEAL, AMBER
from journal_diagram import wire, block, concat, add, vector, tree, forest

OUT=ROOT/'paper/supplementary/figures'
ASSET=OUT/'endpoint_components'

ROWS=[
('tree','ExtraTrees','Partition feature space; average one stored leaf value from every tree.',
 'ARROW MAE 0.202 · Matched and shuffled-response controls',
 'A finite range of leaf values bounds the prediction.'),
('resnet','Residual MLP','Learn nonlinear corrections along an identity-preserving pathway.',
 'ARROW MAE 0.340 · Tests a direct neural replacement',
 'A flexible function alone does not recover the tree accuracy.'),
('dual','Dual structure / response branches','Learn a representation for each input group before combining them.',
 'Size-held-out MAE 1.560 · Refitted ExtraTrees: 1.772',
 'Separate branches do not impose physical additivity.'),
('ridge_residual','Ridge + neural residual','Combine a continuous baseline with a learned correction.',
 'Size-held-out MAE 1.391 · No tree-leaf ceiling in the correction',
 'A linear trend is not a validated asymptotic free-energy law.'),
('graph','Size-sensitive message passing','Count and aggregate atom states, not only fragment presence.',
 'External MAE 1.168 · Without responses: 1.311',
 'Shown: sum / mean / count readout; an atomic-plus-global variant was also tested.'),
('molformer_deterministic','Frozen MoLFormer + neural head','Reuse a pretrained sequence representation alongside response information.',
 'External MAE 1.111 · Responses help in all three cohorts',
 'The encoder stays fixed; the endpoint learns from hydration labels.'),
('tabm','TabM: a parameter-sharing ensemble','Use member-specific factors within layers that share their main weights.',
 'External MAE 1.132 · Size-held-out MAE 1.389',
 '32 members share parameters; these are not 32 independent training runs.'),
('tabpfn','TabPFN-3.5: prediction from context','Condition a pretrained transformer on eligible labelled training rows.',
 'ARROW MAE 0.199 · Size-held-out MAE 1.318',
 'Native intervals under-cover; input recipes differ across partitions.'),
('molformer_fixed_tree','Frozen MoLFormer + ExtraTrees','Change the representation while keeping the tree endpoint family.',
 'ARROW MAE 0.234 · Without responses: 0.348',
 'The endpoint remains bounded by its stored leaf values.')]

def frame(c,y,letter,title,reason,evidence,limit):
    if y:c.line(10,y,710,y,GRID,1)
    c.text(12,y+27,letter,20,INK,700)
    c.text(39,y+26,title,15,INK,600)
    c.text(39,y+48,reason,11.5,MID)
    c.line(39,y+210,85,y+210,TEAL,2)
    c.text(39,y+230,evidence,11.8,TEAL,600)
    c.text(39,y+249,limit,10.7,MID)

def result(c,start,end_y,x=683):
    wire(c,start,(x-16,end_y),TEAL)
    c.text(x,end_y+6,'ŷ',24,TEAL,500,'middle')

def joined_head(c,join,head,y):
    wire(c,join.port('r'),head.port('l'),TEAL)
    result(c,head.port('r'),y)

def diagram(c,key,y):
    Y=lambda v:y+v
    if key=='tree':
        c.text(23,Y(126),'Features',12,BLUE,600)
        wire(c,(82,Y(122)),(132,Y(82)),BLUE,via=[(105,Y(122)),(105,Y(82))],arrow=False)
        roots,ends=forest(c,138,Y(90),311,69)
        c.line(132,Y(82),roots[-1][0],Y(82),BLUE,1.4)
        for root in roots:wire(c,(root[0],Y(82)),(root[0],root[1]-3.5),BLUE)
        for p in ends:wire(c,p,(p[0],Y(184)),TEAL,arrow=False)
        c.line(ends[0][0],Y(184),519,Y(184),TEAL,1.4)
        mean=block(c,537,Y(162),87,44,'Mean')
        wire(c,(519,Y(184)),mean.port('l'),TEAL)
        c.text(302,Y(203),'Selected leaf values',10.8,MID,400,'middle')
        result(c,mean.port('r'),Y(184))
    elif key=='resnet':
        c.text(24,Y(139),'Features',12,BLUE,600)
        wire(c,(85,Y(135)),(125,Y(135)),BLUE)
        for x in [125,337]:
            f=block(c,x,Y(114),113,42,'F1' if x==125 else 'F2',color=BLUE)
            p=add(c,x+147,Y(135),10)
            wire(c,f.port('r'),p.port('l'),BLUE)
            wire(c,(x-10,Y(135)),p.port('t'),TEAL,via=[(x-10,Y(84)),(x+147,Y(84))])
            c.text(x+56,Y(74),'Identity',10.8,TEAL,500,'middle')
            c.text(x+57,Y(178),'Learned correction',10.8,MID,400,'middle')
            wire(c,p.port('r'),(337 if x==125 else 553,Y(135)),TEAL)
        h=block(c,553,Y(114),85,42,'Readout')
        result(c,h.port('r'),Y(135))
    elif key=='dual':
        join=concat(c,408,Y(106),35,65)
        for yy,label,col in [(99,'Structure',BLUE),(174,'Responses',AMBER)]:
            c.text(25,Y(yy+4),label,12,col,600)
            b=block(c,155,Y(yy-20),166,40,'Encoder',color=col)
            wire(c,(109,Y(yy)),b.port('l'),col)
            target=(408,Y(yy+19 if yy==99 else yy-17))
            wire(c,b.port('r'),target,col,via=[(378,Y(yy)),(378,target[1])])
        h=block(c,503,Y(117),128,42,'Joint head')
        joined_head(c,join,h,Y(138))
    elif key=='ridge_residual':
        baseline=block(c,169,Y(82),159,41,'Ridge baseline',color=BLUE)
        residual=block(c,169,Y(156),159,41,'Neural residual')
        c.text(23,Y(99),'Continuous',11.5,BLUE,600)
        c.text(23,Y(114),'descriptors',11.5,BLUE,600)
        c.text(23,Y(181),'All features',11.5,TEAL,600)
        wire(c,(116,Y(102.5)),baseline.port('l'),BLUE)
        wire(c,(116,Y(176.5)),residual.port('l'),TEAL)
        p=add(c,557,Y(139),11)
        wire(c,baseline.port('r'),p.port('t'),BLUE,via=[(557,Y(102.5))])
        wire(c,residual.port('r'),p.port('b'),TEAL,via=[(557,Y(176.5))])
        c.text(407,Y(91),'b(x)',11.5,BLUE,500,'middle')
        c.text(422,Y(195),'r(x)',11.5,TEAL,500,'middle')
        c.text(423,Y(142),'Train r on y − b(x)',10.7,AMBER,500,'middle')
        result(c,p.port('r'),Y(139))
    elif key=='graph':
        # Schematic graph is explicitly an abstract topology, not a molecule.
        pts=[(45,110),(75,91),(101,110),(95,149),(56,155)]
        for a,b in [(0,1),(1,2),(2,3),(3,4),(4,0),(1,4)]:
            c.line(pts[a][0],Y(pts[a][1]),pts[b][0],Y(pts[b][1]),BLUE,1.3)
        for x,yy in pts:c.circle(x,Y(yy),5,'white',BLUE,1.5)
        c.text(75,Y(186),'Graph',11.2,BLUE,500,'middle')
        m=block(c,144,Y(104),125,51,'Message passing',color=BLUE)
        wire(c,(108,Y(129.5)),m.port('l'),BLUE)
        read=block(c,308,Y(104),146,51,'Σhᵢ | mean hᵢ | N',color=TEAL)
        wire(c,m.port('r'),read.port('l'),TEAL)
        c.text(383,Y(91),'Readout',10.7,TEAL,500,'middle')
        join=concat(c,489,Y(110),29,58)
        wire(c,read.port('r'),(489,Y(127)),TEAL)
        c.text(323,Y(191),'RDKit + responses',10.8,AMBER,500)
        wire(c,(441,Y(184)),(489,Y(155)),AMBER,via=[(471,Y(184)),(471,Y(155))])
        h=block(c,557,Y(118),85,42,'Head')
        joined_head(c,join,h,Y(139))
    elif key in ['molformer_deterministic','molformer_fixed_tree']:
        c.text(25,Y(111),'SMILES',12,BLUE,600)
        enc=block(c,129,Y(85),159,47,'MoLFormer',color=BLUE,frozen=True)
        wire(c,(85,Y(108.5)),enc.port('l'),BLUE)
        join=concat(c,354,Y(96),34,75)
        wire(c,enc.port('r'),(354,Y(111)),BLUE)
        c.text(128,Y(167),'Structure + responses',11.2,AMBER,500)
        wire(c,(284,Y(162)),(354,Y(157)),AMBER,via=[(331,Y(162)),(331,Y(157))])
        if key=='molformer_deterministic':
            h=block(c,462,Y(112),153,43,'Neural head')
            joined_head(c,join,h,Y(133.5))
        else:
            roots,ends=forest(c,431,Y(96),130,55)
            wire(c,join.port('r'),(418,Y(83)),TEAL,via=[(408,Y(133.5)),(408,Y(83))],arrow=False)
            c.line(418,Y(83),roots[-1][0],Y(83),TEAL,1.3)
            for root in roots:wire(c,(root[0],Y(83)),(root[0],root[1]-3.5),TEAL)
            for p in ends:wire(c,p,(p[0],Y(168)),TEAL,arrow=False)
            c.line(ends[0][0],Y(168),570,Y(168),TEAL,1.3)
            c.text(496,Y(70),'ExtraTrees',11.2,TEAL,600,'middle')
            h=block(c,582,Y(146),65,44,'Mean')
            wire(c,(570,Y(168)),h.port('l'),TEAL)
            result(c,h.port('r'),Y(168))
    elif key=='tabm':
        c.text(23,Y(142),'Features',12,BLUE,600)
        wire(c,(84,Y(138)),(234,Y(138)),BLUE,arrow=False)
        c.text(153,Y(122),'Within a layer',10.7,MID,400,'middle')
        c.line(234,Y(92),234,Y(184),BLUE,1.3)
        c.text(366,Y(69),'Same W · different rᵢ and sᵢ',11.2,AMBER,500,'middle')
        for j in range(3):
            yy=92+j*46
            b=block(c,263,Y(yy-15),207,30,'diag(sᵢ) W diag(rᵢ) x',color=AMBER)
            wire(c,(234,Y(yy)),b.port('l'),BLUE)
            wire(c,b.port('r'),(505,Y(yy)),TEAL,arrow=False)
        c.line(505,Y(92),505,Y(184),TEAL,1.3)
        c.text(585,Y(184),'32 member outputs',10.7,MID,400,'middle')
        h=block(c,548,Y(116),86,44,'Mean')
        wire(c,(505,Y(138)),h.port('l'),TEAL)
        result(c,h.port('r'),Y(138))
    elif key=='tabpfn':
        c.text(25,Y(77),'Training context',11.5,BLUE,600)
        for i in range(3):
            for j in range(4):c.rect(25+j*23,Y(88+i*19),19,14,'white',BLUE,1)
            c.rect(127,Y(88+i*19),22,14,'white',AMBER,1)
        c.text(137,Y(77),'y',11.5,AMBER,500,'middle')
        c.text(25,Y(165),'Query',11.5,TEAL,600)
        for j in range(4):c.rect(25+j*23,Y(177),19,14,'white',TEAL,1)
        c.text(138,Y(190),'?',16,AMBER,500,'middle')
        enc=block(c,251,Y(113),164,52,'TabPFN',color=BLUE,frozen=True)
        wire(c,(158,Y(114)),(251,Y(130)),BLUE,via=[(207,Y(114)),(207,Y(130))])
        wire(c,(158,Y(184)),(251,Y(150)),TEAL,via=[(222,Y(184)),(222,Y(150))])
        c.text(333,Y(94),'Pretrained transformer',11.2,BLUE,500,'middle')
        wire(c,enc.port('r'),(453,Y(139)),TEAL)
        c.text(499,Y(146),'p(y | x, D)',17,TEAL,500,'middle')
        c.text(499,Y(176),'Predictive distribution',10.3,MID,400,'middle')
        h=block(c,573,Y(117),65,44,'Mean')
        wire(c,(551,Y(139)),h.port('l'),TEAL)
        result(c,h.port('r'),Y(139))

def main():
    p=argparse.ArgumentParser();p.add_argument('--preview-dir',type=Path);a=p.parse_args()
    ASSET.mkdir(exist_ok=True)
    for group in range(3):
        c=Art(788,'Endpoint mechanisms: '+', '.join(r[1] for r in ROWS[group*3:group*3+3]))
        for j,row in enumerate(ROWS[group*3:group*3+3]):
            key,title,reason,evidence,limit=row;y=j*263
            frame(c,y,'abc'[j],title,reason,evidence,limit);diagram(c,key,y)
            solo=Art(263,title);frame(solo,0,'',title,reason,evidence,limit);diagram(solo,key,0)
            (ASSET/(key+'.svg')).write_text(solo.svg()+'\n')
        name=f'Supp_Fig{8+group}_endpoint_mechanisms'
        c.save(OUT/name,a.preview_dir/(name+'.png') if a.preview_dir else None)
    print('Nine native vector mechanisms; connected ports; no generated scientific topology.')

if __name__=='__main__':main()

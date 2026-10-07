"""Exact vector topology with an Azure-generated paper-texture accent.

Texture is non-evidential: native paths define all six source streams. Molecular
connectivity, coordinate counts, labels and values come from retained sources.
"""
import base64
from journal_vector_art import Art, ROOT, INK, MID, GRID, BLUE, TEAL, AMBER, COLORS, COUNTS
from journal_diagram import wire, block, concat, bank, heading

OUT=ROOT/'paper/figures/journal'
ASSET=OUT/'teaser_components'

def output(c,x,y,size=25):
    c.text(x,y,'ΔG',size,TEAL,500)
    c.text(x+size*1.58,y+size*.22,'hyd',size*.43,TEAL)

def texture(c):
    for j in range(6):
        b=base64.b64encode((ASSET/f'paper_texture_{j+1}.png').read_bytes()).decode()
        c.defs.append(f'<pattern id="paper{j}" width="100" height="30" patternUnits="userSpaceOnUse"><image width="100" height="30" preserveAspectRatio="none" xlink:href="data:image/png;base64,{b}"/></pattern>')

def overview(data,preview):
    c=Art(640,'SolvAI: turn repeated solvent calculations into reusable molecular information')
    texture(c)
    heading(c,10,27,'a','Repeated calculation')
    c.text(36,49,'Gas-to-water free energy',11.7,MID)
    c.molecule(28,82,147,105)
    c.text(103,207,'New molecule',12,INK,500,'middle')
    wire(c,(104,218),(104,252),AMBER,width=1.8)
    c.shell(19,266,167,111,explicit=True)
    c.text(103,401,'Sample solvent configurations',11.3,AMBER,500,'middle')
    c.text(103,418,'Repeat the free-energy calculation',11.3,MID,400,'middle')
    c.line(213,61,213,421,GRID,1)

    heading(c,235,27,'b','Learn once, reuse across molecules')
    c.text(260,49,'Six sources → six frozen structure-to-response maps',11.7,MID)
    names=['COSMOtherm','Abraham','OpenFF + δ','GBn2 + δ','SMD(water)','ConfSolv']
    desc=['continuum water','empirical axes','explicit water','implicit water','continuum water','conformer response']
    for j,(name,sub,col,n) in enumerate(zip(names,desc,COLORS,COUNTS)):
        yy=87+j*47
        c.text(246,yy,name,12.3,col,600)
        c.text(246,yy+15,sub,10.2,MID)
        start=381;stop=568;end=104+j*39;th=12
        path=(f'M{start},{yy-7} C447,{yy-7} 453,{end-15} {stop},{end-15} '
              f'L{stop},{end-15+th} C453,{end-15+th} 447,{yy+5} {start},{yy+5} Z')
        c.path(path,col,.6,f'url(#paper{j})')
        c.circle(577,end-9,12,'white',col,1.25)
        c.text(577,end-5,'φ'+str(j+1),10.3,col,500,'middle')
        wire(c,(589,end-9),(611,end-9),col,width=1.25)
        for k in range(n):c.rect(619+k*13,end-17,9,16,col,radius=1)
    c.text(560,367,'Learned maps',11.4,INK,500,'end')
    c.text(659,367,'15 coordinates',12.3,TEAL,600,'middle')
    c.text(249,413,'δ: learned correction to a physical calculation',10.7,MID)
    c.line(10,443,710,443,GRID,1)

    heading(c,10,472,'c','Predict hydration from structure alone')
    c.text(35,495,'No new solvent calculation at inference',11.5,TEAL)
    c.text(45,532,'SMILES',11.5,BLUE,600,'middle')
    b=bank(c,104,505,82,53)
    wire(c,(75,530),b.port('l'),BLUE)
    c.text(145,577,'Six fixed maps',10.9,TEAL,500,'middle')
    wire(c,b.port('r'),(208,530),TEAL)
    c.responses(217,522,183,16)
    c.text(309,558,'15 response descriptors',10.9,TEAL,500,'middle')
    wire(c,(83,530),(419,566),MID,via=[(83,595),(409,595),(409,566)])
    c.rect(217,583,171,17,'white')
    c.feature_stripes(228,585,151,10)
    c.text(304,615,'Structural descriptors',10.7,MID,400,'middle')
    join=concat(c,419,515,32,52)
    wire(c,(402,530),(419,530),TEAL)
    h=block(c,482,514,113,54,'Hydration model',symbol='h',frozen=True)
    wire(c,join.port('r'),h.port('l'),TEAL)
    wire(c,h.port('r'),(630,541),TEAL)
    output(c,647,549,25)
    return c.save(OUT/'F1_overview',preview/'F1_overview.png' if preview else None)

def learning(data,preview):
    c=Art(706,'SolvAI: source supervision, endpoint fitting and evaluation exposure')
    heading(c,10,27,'a','Distil each source into a structure-to-response map')
    c.text(35,49,'Source targets are computed quantities or empirical measurements',11.6,MID)
    c.text(190,77,'Frozen surrogate',10.4,MID,500,'middle')
    c.text(359,77,'Predicted coordinates',10.4,MID,500,'middle')
    c.text(488,77,'Source training targets',10.4,AMBER,500)
    labels=['COSMOtherm water','Abraham E, S, A, B, L','OpenFF + correction',
            'GBn2 + correction','SMD(water)','ConfSolv summaries']
    types=['D-MPNN','ExtraTrees','ExtraTrees','ExtraTrees','D-MPNN','LightGBM']
    c.molecule(15,145,83,64)
    c.text(57,235,'Molecular',11.5,INK,400,'middle')
    c.text(57,251,'structure',11.5,INK,400,'middle')
    wire(c,(97,176),(115,176),BLUE,arrow=False)
    c.line(115,110,115,287,BLUE,1.2)
    for j,(name,n,col,typ) in enumerate(zip(labels,COUNTS,COLORS,types)):
        y=94+j*35.4
        b=block(c,145,y,119,29,'φ'+str(j+1)+'  ·  '+typ,color=col)
        wire(c,(115,y+14.5),b.port('l'),col,width=1.2)
        wire(c,b.port('r'),(311,y+14.5),col,width=1.2)
        for k in range(n):c.rect(322+k*22,y+6,16,17,col,radius=1)
        c.text(488,y+20,name,11.8,INK)
    c.text(359,328,'1 + 5 + 1 + 1 + 1 + 6 = 15',12.2,TEAL,600,'middle')
    c.line(10,348,710,348,GRID,1)

    heading(c,10,377,'b','Learn the endpoint separately')
    c.text(35,400,'Experimental hydration labels supervise the endpoint, not its source targets',11.6,MID)
    c.responses(35,431,205,15)
    c.text(137,467,'Response descriptors',11.4,TEAL,500,'middle')
    c.feature_stripes(35,490,205,12)
    c.text(137,522,'Structural descriptors',11.4,MID,400,'middle')
    join=concat(c,280,437,34,62)
    wire(c,(248,438),(280,450),TEAL,via=[(264,438),(264,450)])
    wire(c,(248,496),(280,486),MID,via=[(264,496),(264,486)])
    h=block(c,378,440,154,56,'Hydration model',symbol='h')
    wire(c,join.port('r'),h.port('l'),TEAL)
    c.text(455,425,'Experimental ΔG labels · training folds',10.3,AMBER,500,'middle')
    wire(c,(455,430),h.port('t'),AMBER)
    wire(c,h.port('r'),(586,468),TEAL)
    output(c,609,477,28)
    c.text(455,526,'Swap the endpoint; keep the response layer fixed',11.3,MID,400,'middle')
    c.line(10,546,710,546,GRID,1)

    heading(c,10,575,'c','Distinguish endpoint exclusion from source exclusion')
    c.text(35,608,'85 ARROW molecules',14,BLUE,600)
    c.text(35,632,'Absent from all six supervised sources',11.4,INK)
    c.text(35,657,'Each hydration prediction is out of fold:',11.3,MID)
    c.text(35,676,'its label is excluded from that endpoint fit.',11.3,MID)
    c.line(336,595,336,688,GRID,1)
    c.text(355,608,'220 external molecules',14,BLUE,600)
    c.text(355,629,'All absent from endpoint training',11.4,INK)
    for i in range(220):
        c.circle(360+(i%22)*6.4,645+(i//22)*4.7,1.45,MID if i<123 else TEAL)
    c.text(519,656,'123 source-exposed',11.2,MID,600)
    c.text(519,677,'97 absent from all six',11.2,TEAL,600)
    c.text(710,701,'97 is a subset of 220, not a separate fit.',10.5,MID,400,'end')
    return c.save(OUT/'F2_learning_and_evaluation',preview/'F2_learning_and_evaluation.png' if preview else None)

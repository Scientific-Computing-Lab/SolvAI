"""Compose the two opening figures from verified chemistry, data and concept art.

Selected raster objects depict information reuse, not molecules, fields or a
particular architecture. Every scientific label, connector and value is vector.
"""
from pathlib import Path
import base64
from journal_vector_art import Art, ROOT, INK, MID, GRID, BLUE, TEAL, AMBER, PURPLE, COLORS, COUNTS

OUT=ROOT/'paper/figures/journal'
ASSET=OUT/'teaser_components'

def image(c,name,x,y,w,h):
    data=base64.b64encode((ASSET/(name+'.png')).read_bytes()).decode()
    c.add(f'<image x="{x}" y="{y}" width="{w}" height="{h}" preserveAspectRatio="xMidYMid meet" xlink:href="data:image/png;base64,{data}"/>')

def heading(c,x,y,letter,title,subtitle=None):
    c.text(x,y,letter,20,INK,700)
    c.text(x+24,y-1,title,14.2,INK,700)
    if subtitle:c.text(x+24,y+18,subtitle,11.6,MID)

def output(c,x,y,size=25):
    c.text(x,y,'ΔG',size,TEAL,600)
    c.text(x+size*1.65,y+size*.23,'hyd',size*.46,TEAL)

def response_bar(c,x,y,w,h=13,labels=False):
    c.responses(x,y,w,h)
    if labels:
        gap=3;gg=7;cell=(w-gap*(sum(COUNTS)-len(COUNTS))-gg*(len(COUNTS)-1))/sum(COUNTS)
        cur=x
        for count,color in zip(COUNTS,COLORS):
            span=count*cell+(count-1)*gap
            c.text(cur+span/2,y+h+16,str(count),11.5,color,600,'middle')
            cur+=span+gg

def ai(c,x,y,w,h,locked=False):
    image(c,'endpoint',x,y,w,h)
    c.text(x+w*.50,y+h*.66,'AI',19,TEAL,700,'middle')
    if locked:c.lock(x+w*.88,y+8,TEAL,1.25)

def overview(data,preview):
    c=Art(690,'SolvAI: learn solvent responses once and reuse them from structure')
    heading(c,10,26,'a','The recurring physical cost')
    c.text(110,65,'Gas → water',14,INK,600,'middle')
    c.molecule(53,79,117,81)
    c.text(110,180,'One new molecule',12,MID,400,'middle')
    c.arrow(109,192,109,215,AMBER,2)
    c.shell(26,222,164,109,explicit=True)
    c.text(110,350,'Sample and calculate again',12,AMBER,600,'middle')
    c.path('M181,249 C199,217 192,167 174,148',AMBER,1.3)
    c.arrow(174,148,170,144,AMBER,1.3)
    c.line(213,45,213,357,GRID,1)

    heading(c,235,26,'b','Learn reusable solvent responses')
    c.text(260,60,'Six complementary source families',12,MID)
    labels=[('COSMOtherm','Continuum water'),('Abraham','Five empirical axes'),
            ('OpenFF + δ','Explicit water'),('GBn2 + δ','Implicit water'),
            ('SMD(water)','Continuum water'),('ConfSolv','Conformer summaries')]
    ys=[95,140,185,230,275,320]
    for (name,desc),y,color in zip(labels,ys,COLORS):
        c.circle(244,y-4,3.3,color)
        c.text(255,y,name,12.4,color,600)
        c.text(255,y+15,desc,10.7,MID)
    image(c,'response_bank',365,75,334,265)
    # Each colored ribbon is conceptual; exact scientific counts are native.
    c.responses(413,339,279,13)
    c.text(553,371,'15 structure-predicted descriptors',12.7,TEAL,600,'middle')
    c.text(270,371,'δ: learned residual',10.5,MID)
    c.line(10,393,710,393,GRID,1)

    heading(c,10,422,'c','Reuse the learned mappings for a new molecule')
    c.text(52,451,'SMILES',11.6,INK,600,'middle')
    c.molecule(12,456,79,63)
    c.text(53,537,'CNC(C)=O',11.1,MID,400,'middle')
    c.arrow(100,492,143,492,BLUE,1.7)
    image(c,'response_core',150,446,62,95)
    c.lock(201,446,TEAL,1.2)
    c.text(181,559,'Six frozen mappings',11.5,TEAL,600,'middle')
    c.arrow(217,492,261,492,TEAL,1.7)
    c.responses(273,482,128,20)
    c.text(337,524,'15 responses',11.5,TEAL,600,'middle')
    c.path('M408,492 L433,492 L433,512',TEAL,1.6)
    c.path('M112,492 L112,567 L434,567 L434,525',MID,1.4)
    c.circle(112,492,2.5,BLUE)
    c.feature_stripes(275,555,103,10)
    c.text(327,590,'Structural descriptors',11.2,MID,400,'middle')
    c.circle(434,518,6,'white',TEAL,1.1)
    c.text(434,522,'+',12,TEAL,600,'middle')
    c.text(434,548,'Concat.',10.5,TEAL,400,'middle')
    c.arrow(442,518,470,518,TEAL,1.7)
    ai(c,478,485,130,60,locked=True)
    c.text(543,567,'Learned endpoint',11.5,TEAL,600,'middle')
    c.arrow(613,518,642,518,TEAL,1.7);output(c,650,527,22)
    c.text(700,590,'No new solvent calculation',12,TEAL,600,'end')
    c.line(10,607,710,607,GRID,1)

    heading(c,10,635,'d','Original matched evidence')
    c.text(34,657,'ExtraTrees endpoint',11.4,MID)
    c.text(34,677,'MAE (kcal mol⁻¹)',11.4,MID)
    # Exact data-derived values, no generated numeric content.
    for cohort,cx in [('ARROW-85',346),('External-220',591)]:
        f=data['modern'];g=f.loc[f.cohort.eq(cohort)]
        a=float(g.loc[g.model.eq('structure'),'mae'].item())
        b=float(g.loc[g.model.eq('solvai'),'mae'].item())
        c.text(cx,634,cohort,12,INK,600,'middle')
        c.text(cx-43,660,f'{a:.3f}',19,MID,500,'end')
        c.arrow(cx-34,654,cx-2,654,TEAL,1.5)
        c.text(cx+8,660,f'{b:.3f}',19,TEAL,700)
        c.text(cx-68,680,'Structure',10.5,MID,400,'middle')
        c.text(cx+35,680,'+ responses',10.5,TEAL,600,'middle')
    return c.save(OUT/'F1_overview',preview/'F1_overview.png' if preview else None)

def learning(data,preview):
    c=Art(690,'SolvAI: separate source learning, endpoint fitting and evaluation exposure')
    heading(c,10,26,'a','Learn the source mappings, then freeze them')
    c.text(34,50,'Source supervision: calculated quantities or empirical measurements',11.8,MID)
    c.molecule(20,119,83,72)
    c.text(62,215,'Molecular\nstructure',12,INK,400,'middle')
    c.arrow(104,151,133,151,BLUE)
    labels=['COSMOtherm water','Abraham E, S, A, B, L','OpenFF + correction',
            'GBn2 + correction','SMD(water)','ConfSolv summaries']
    for j,(name,count,col) in enumerate(zip(labels,COUNTS,COLORS)):
        y=83+j*29
        c.path(f'M133,151 C153,151 140,{y+8} 172,{y+8}',col,1.2)
        c.rect(180,y-3,81,23,'#FFFFFF',col,4)
        c.text(220,y+13,'φ'+str(j+1),12,col,600,'middle')
        c.lock(279,y+4,col,.95)
        c.arrow(295,y+8,325,y+8,col,1.2)
        for k in range(count):c.rect(335+k*19,y,14,15,col,radius=1)
        c.text(483,y+11,name,12,INK)
    c.text(221,271,'Six distinct surrogate families',11.4,TEAL,600,'middle')
    c.text(485,271,'1 + 5 + 1 + 1 + 1 + 6 = 15 coordinates',11.4,TEAL,600,'middle')
    c.line(10,291,710,291,GRID,1)

    heading(c,10,321,'b','Learn hydration from responses and structure')
    c.text(34,344,'Endpoint supervision: experimental hydration free energy',11.8,MID)
    response_bar(c,35,373,194,16)
    c.text(132,413,'15 frozen response predictions',11.7,TEAL,600,'middle')
    c.feature_stripes(37,432,191,13)
    c.text(132,466,'Structural descriptors',11.7,MID,400,'middle')
    c.path('M239,381 L267,381 L267,412',TEAL,1.4)
    c.path('M239,439 L267,439 L267,424',MID,1.4)
    c.circle(267,418,6,'white',TEAL,1.1);c.text(267,422,'+',12,TEAL,600,'middle')
    c.text(292,372,'Concatenate',10.7,MID,400,'middle')
    c.arrow(275,418,330,418,TEAL,1.8)
    ai(c,343,378,193,87)
    c.text(439,477,'Endpoint model h',12.4,TEAL,600,'middle')
    c.arrow(543,418,585,418,TEAL,1.8);output(c,601,426,27)
    c.text(612,465,'Hydration prediction',11.7,TEAL,600,'middle')
    c.text(695,354,'Training-fold labels only',11.2,AMBER,600,'end')
    c.path('M590,363 L440,363',AMBER,1.3)
    c.arrow(440,363,440,380,AMBER,1.3)
    c.text(34,498,'The response layer stays fixed when the endpoint architecture changes.',11.7,MID)
    c.line(10,519,710,519,GRID,1)

    heading(c,10,548,'c','Separate endpoint exclusion from source exclusion')
    c.text(34,581,'85 ARROW molecules',14,BLUE,600)
    c.text(34,605,'Absent from all six supervised sources',11.6,INK)
    c.text(34,630,'Each prediction is out of fold:',11.6,MID)
    c.text(34,650,'its hydration label is not used in that fit.',11.6,MID)
    c.line(335,567,335,676,GRID,1)

    c.text(354,581,'220 external molecules',14,BLUE,600)
    c.text(354,601,'All absent from endpoint training',11.6,INK)
    # Each native dot is one molecule: exact123/97partition, not two cohorts.
    for i in range(220):
        col=MID if i<123 else TEAL
        c.circle(361+(i%22)*7.1,621+(i//22)*5.7,1.75,col)
    c.text(544,626,'123 source-exposed',11.5,MID,600)
    c.text(544,648,'97 absent from all six',11.5,TEAL,600)
    c.text(544,671,'Nested subsets, same fitted model',10.4,MID)
    return c.save(OUT/'F2_learning_and_evaluation',preview/'F2_learning_and_evaluation.png' if preview else None)

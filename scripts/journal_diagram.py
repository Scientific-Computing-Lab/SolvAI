"""Native diagram grammar with explicitly anchored ports and connected geometry."""
from dataclasses import dataclass
import math
from journal_vector_art import Art, INK, MID, GRID, BLUE, TEAL, AMBER, COLORS, COUNTS

@dataclass(frozen=True)
class Box:
    x: float
    y: float
    w: float
    h: float
    def port(self,side):
        return {'l':(self.x,self.y+self.h/2),'r':(self.x+self.w,self.y+self.h/2),
                't':(self.x+self.w/2,self.y),'b':(self.x+self.w/2,self.y+self.h)}[side]

def wire(c,start,end,color=TEAL,via=(),arrow=True,dashed=False,width=1.5):
    points=[start,*via,end]
    for a,b in zip(points,points[1:]):c.line(*a,*b,color,width,'4 3' if dashed else None)
    if arrow:
        a,b=points[-2:];dx,dy=b[0]-a[0],b[1]-a[1];d=math.hypot(dx,dy)
        assert d>0
        u=(dx/d,dy/d);v=(-u[1],u[0]);length=5.5;half=2.3
        p=(b[0]-length*u[0]+half*v[0],b[1]-length*u[1]+half*v[1])
        q=(b[0]-length*u[0]-half*v[0],b[1]-length*u[1]-half*v[1])
        c.path(f'M{p[0]},{p[1]} L{b[0]},{b[1]} L{q[0]},{q[1]} Z',color,0,color)

def block(c,x,y,w,h,title,subtitle=None,color=TEAL,symbol=None,frozen=False):
    """Function module, not a claim about network layer count or physical form."""
    b=Box(x,y,w,h)
    c.rect(x,y,w,h,'#FFFFFF',color,5)
    c.line(x+1,y+6,x+1,y+h-6,color,2.5)
    if symbol:
        c.text(x+w/2,y+26,symbol,23,color,500,'middle',True)
        c.text(x+w/2,y+h-13,title,11.7,color,600,'middle')
    else:
        c.text(x+w/2,y+h/2+(0 if subtitle else 4),title,12,color,600,'middle')
        if subtitle:c.text(x+w/2,y+h/2+17,subtitle,10.3,MID,400,'middle')
    if frozen:
        c.rect(x+w-48,y-8,42,16,'white')
        c.text(x+w-27,y+3,'fixed',9.5,color,500,'middle')
    return b

def concat(c,x,y,w=44,h=45):
    b=Box(x,y,w,h)
    c.path(f'M{x+6},{y} H{x} V{y+h} H{x+6}',INK,1.4)
    c.path(f'M{x+w-6},{y} H{x+w} V{y+h} H{x+w-6}',INK,1.4)
    for row in range(2):
        for j in range(3):
            c.rect(x+6+j*(w-12)/3,y+h*(.27+.36*row),(w-15)/3,5,INK if row==0 else MID,radius=.5)
    c.text(x+w/2,y+h+18,'join',10.5,MID,400,'middle')
    return b

def add(c,x,y,r=11):
    c.circle(x,y,r,'white',TEAL,1.4);c.text(x,y+5,'+',20,TEAL,400,'middle')
    return Box(x-r,y-r,2*r,2*r)

def vector(c,x,y,w,n=6,color=BLUE,h=14):
    gap=2.5;cell=(w-(n-1)*gap)/n
    for j in range(n):c.rect(x+j*(cell+gap),y,cell,h,color,radius=.6)
    return Box(x,y,w,h)

def response_matrix(c,x,y,w=138,h=68):
    row=h/6
    for j,(n,col) in enumerate(zip(COUNTS,COLORS)):
        vector(c,x,y+j*row,w,n,col,row-4)
    return Box(x,y,w,h)

def bank(c,x,y,w=126,h=72):
    """Six separate source functions with the established color identities."""
    for j,col in enumerate(COLORS):
        yy=y+j*h/6
        c.rect(x,yy,w,h/6-3,'white',col,2)
        if h>=72:c.text(x+w/2,yy+h/6-4,'φ'+str(j+1),9.5,col,500,'middle')
    return Box(x,y,w,h-3)

def tree(c,x,y,w=100,h=60,color=TEAL,selected=1):
    """Deterministic tree geometry; edge endpoints meet the actual native nodes."""
    root=(x+w/2,y+5)
    inner=[(x+w*.25,y+h*.45),(x+w*.75,y+h*.45)]
    leaves=[(x+w*f,y+h-5) for f in [.08,.37,.63,.92]]
    nodes=[root,*inner,*leaves]
    edges=[(0,1),(0,2),(1,3),(1,4),(2,5),(2,6)]
    chosen=3+selected;parent=1 if selected<2 else 2
    for a,b in edges:
        col=color if (a,b) in [(0,parent),(parent,chosen)] else '#AEBCC3'
        c.line(*nodes[a],*nodes[b],col,1.35 if col==color else 1)
    for j,p in enumerate(nodes):
        active=j in [0,parent,chosen]
        if j<3:c.circle(*p,3.5,'white',color if active else '#AEBCC3',1.1)
        else:
            lw=min(10,w*.17)
            c.rect(p[0]-lw/2,p[1]-3.5,lw,7,color if active else 'white',color if active else '#AEBCC3',1)
    return root,(leaves[selected][0],leaves[selected][1]+3.5)

def forest(c,x,y,w=126,h=54,color=TEAL):
    roots=[];outputs=[]
    for j in range(3):
        a,b=tree(c,x+j*w/3,y,w/3-6,h,color,selected=j%3)
        roots.append(a);outputs.append(b)
    return roots,outputs

def heading(c,x,y,letter,title,subtitle=None):
    c.text(x,y,letter,19,INK,700)
    c.text(x+25,y-1,title,14.4,INK,600)
    if subtitle:c.text(x+25,y+18,subtitle,11.3,MID)

def chip(c,x,y,w=92,h=44,color=TEAL,label='Learned model'):
    return block(c,x,y,w,h,label,color=color)

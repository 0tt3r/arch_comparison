"""Editable vector cartoons of three ways to allocate the same qubit array.

These illustrations explain resource roles. Patch sizes and allocation ratios
do not encode the numerical floorplans used by the cost models.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Rectangle,FancyArrowPatch,Polygon
import numpy as np

from project_paths import ROOT
from plot_capacity import COLORS

INK='#273240';MUTED='#627080';GRID='#dce2e8'


def block(ax,x,y,w,h,color,label=None,face=None,fontsize=9,dashed=False):
    patch=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.16,rounding_size=0.5',
                        linewidth=1.35,edgecolor=color,facecolor=face or color,
                        alpha=1 if face else .10,linestyle='--' if dashed else '-',zorder=2)
    ax.add_patch(patch)
    if face is None:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.16,rounding_size=0.5',
                     linewidth=1.35,edgecolor=color,facecolor='none',linestyle='--' if dashed else '-',zorder=4))
    if label:ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=fontsize,
                     color=INK,zorder=6,bbox={'facecolor':'white','edgecolor':'none','alpha':.87,'pad':1.6})


def arrow(ax,a,b,color=INK,style='-|>'):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=11,color=color,lw=1.5,zorder=7))


def hardware(ax):
    for x in range(44):ax.plot([x,x],[0,31],color=GRID,lw=.28,zorder=0)
    for y in range(32):ax.plot([0,43],[y,y],color=GRID,lw=.28,zorder=0)
    xx,yy=np.meshgrid(np.arange(44),np.arange(32))
    ax.scatter(xx.ravel(),yy.ravel(),s=3.0,color='#9da9b6',alpha=.66,zorder=3,linewidths=0)
    ax.add_patch(Rectangle((-.7,-.7),44.4,32.4,facecolor='none',edgecolor='#7b8997',lw=1,zorder=4))
    ax.set_xlim(-1.5,44.5);ax.set_ylim(-10,38);ax.set_aspect('equal');ax.axis('off')


def data_and_surgery(ax,color):
    # Operation space connects code patches to rotation or magic-state supply.
    ax.add_patch(Rectangle((1,17.1),40,3.8,facecolor='#e3e7ed',edgecolor='none',zorder=1))
    for x in (9.4,19.4):ax.add_patch(Rectangle((x,9),2.4,21,facecolor='#e3e7ed',edgecolor='none',zorder=1))
    for i,(x,y) in enumerate((x,y) for y in (22,9) for x in (1.5,11.5,21.5)):
        block(ax,x,y,7.2,7.2,color,rf'$\bar q_{i}$',fontsize=12)
    ax.text(15,19,'surgery / routing',ha='center',va='center',fontsize=8,color=MUTED,zorder=6,
            bbox={'facecolor':'#e3e7ed','edgecolor':'none','pad':1})
    block(ax,1,1,40.5,4.8,'#98a5b3','Spare area → additional protected samplers',fontsize=8,dashed=True)


def overview():
    fig,axes=plt.subplots(1,3,figsize=(15.8,6.9))
    titles=['NISQ','STAR','FTQC']
    subtitles=['Physical circuits + repeated sampling','Protected Clifford gates + analog rotations','Protected gates + high-fidelity magic states']
    for ax,title,subtitle,arm in zip(axes,titles,subtitles,('nisq','star','ft')):
        hardware(ax)
        ax.text(0,36,title,fontsize=19,fontweight='bold',color=COLORS[arm],ha='left')
        ax.text(0,33.3,subtitle,fontsize=9.0,color=MUTED,ha='left')
    ax=axes[0];color=COLORS['nisq']
    for i,(x,y) in enumerate((x,y) for y in (22,12,2) for x in (1,15,29)):
        block(ax,x,y,12,8,color,f'Sampler {i+1}',fontsize=9)
    ax.text(21.5,-3.1,'Prepare → evolve → measure',ha='center',fontsize=11,color=INK)
    ax.text(21.5,-6.0,'Pool shots; select and mitigate',ha='center',fontsize=10,color=MUTED)
    ax=axes[1];color=COLORS['star'];data_and_surgery(ax,color)
    block(ax,32,21.5,9.4,8.0,color,'Analog\ninjection',fontsize=8.5)
    block(ax,32,9.2,9.4,7.0,color,'Grow / retry',fontsize=8)
    arrow(ax,(36.5,21.2),(36.5,16.6),color)
    arrow(ax,(34.0,19.0),(29.4,19.0),color)
    ax.text(21.5,-3.1,'Logical patches + rotation delivery',ha='center',fontsize=11,color=INK)
    ax.text(21.5,-6.0,'Mitigate residual analog errors',ha='center',fontsize=10,color=MUTED)
    ax=axes[2];color=COLORS['ft'];data_and_surgery(ax,color)
    block(ax,32,21.5,9.4,8.0,color,'T-state\nbuffers',fontsize=8.5)
    block(ax,32,9.2,9.4,7.0,color,'Factory\ncultivate / distill',fontsize=7.6)
    arrow(ax,(39.7,16.5),(39.7,21.2),color)
    block(ax,31,17.6,3,3,color,'T',fontsize=9,face='#f0ebfa')
    arrow(ax,(36.4,22),(33.1,20.6),color)
    arrow(ax,(30.7,19.0),(28.6,19.0),color)
    ax.text(21.5,-3.1,'Logical patches + magic-state supply',ha='center',fontsize=11,color=INK)
    ax.text(21.5,-6.0,'Synthesize rotations; budget logical faults',ha='center',fontsize=10,color=MUTED)
    fig.suptitle('The same sea of physical qubits, allocated in three ways',fontsize=18,x=.04,ha='left',y=.985)
    fig.text(.5,.028,'Dots: physical qubits · Gray corridors: operation space · Colored outlines: resource groups · Schematic allocations, not optimized floorplans',
             ha='center',fontsize=9,color=MUTED)
    fig.subplots_adjust(left=.035,right=.985,top=.92,bottom=.065,wspace=.11)
    return fig


def patch_detail():
    fig,axes=plt.subplots(1,2,figsize=(11.6,5.2),gridspec_kw={'width_ratios':[1,1.35]})
    ax=axes[0];ax.set_aspect('equal');ax.axis('off')
    # A distance-three rotated-code illustration. Rotating these coordinates
    # makes each data-check connection an edge of a square nearest-neighbor grid.
    def xy(x,y):return ((x+y)/2,(y-x)/2)
    checks=[(1,1,'X'),(3,3,'X'),(3,1,'Z'),(1,3,'Z'),
            (1,5,'X'),(3,-1,'X'),(-1,1,'Z'),(5,3,'Z')]
    data=[(x,y) for x in (0,2,4) for y in (0,2,4)]
    for x,y,kind in checks:
        neighbors=[q for q in data if abs(x-q[0])==1 and abs(y-q[1])==1]
        for q in neighbors:ax.plot(*zip(xy(x,y),xy(*q)),color='#c5cdd6',lw=1.25,zorder=1)
    for x,y in data:ax.scatter(*xy(x,y),s=115,color=COLORS['ft'],zorder=3,edgecolor='white',linewidth=.8)
    for x,y,kind in checks:
        color=COLORS['nisq'] if kind=='X' else COLORS['star']
        ax.scatter(*xy(x,y),s=150,facecolor='white',edgecolor=color,linewidth=2,zorder=3)
        ax.text(*xy(x,y),kind,ha='center',va='center',fontsize=8,color=color,zorder=4)
    ax.text(2,3.0,'Inside one logical patch',ha='center',fontsize=14,color=INK)
    ax.text(2,-3.2,'Data qubits + repeated syndrome measurements',ha='center',fontsize=9.4,color=MUTED)
    ax.text(2,-3.7,'Distance-three illustration: 9 data + 8 check qubits',ha='center',fontsize=8.5,color=MUTED)
    ax.set_xlim(-1.7,5.7);ax.set_ylim(-4.0,3.6)
    ax=axes[1];ax.axis('off');ax.set_xlim(0,12);ax.set_ylim(0,8)
    ax.text(0,7.6,'At the architecture scale',ha='left',fontsize=14,color=INK)
    block(ax,.2,3.7,3.1,2.6,COLORS['ft'],'Logical data\npatches',fontsize=10)
    block(ax,4.25,3.7,3.0,2.6,'#8493a3','Surgery space\n+ buffers',fontsize=10,face='#edf0f4')
    block(ax,8.25,5,3.25,1.5,COLORS['star'],'Analog injection',fontsize=9)
    block(ax,8.25,2.5,3.25,1.5,COLORS['ft'],'T-state supply',fontsize=9)
    arrow(ax,(8.1,5.7),(7.4,5.4),COLORS['star'])
    arrow(ax,(8.1,3.3),(7.4,4.5),COLORS['ft'])
    arrow(ax,(4.05,5.0),(3.45,5.0),INK,style='<->')
    ax.text(.2,2.55,'STAR uses the analog path; FTQC uses the T-state path.',fontsize=9.5,color=INK)
    ax.text(.2,1.65,'The hybrid uses both: fixed T states + analog small angles.',fontsize=9.5,color=INK)
    ax.text(.2,.55,'More area can buy greater distance, faster supply,\nmore operation space, or more complete sampler replicas.',fontsize=10,color=MUTED,linespacing=1.5)
    fig.subplots_adjust(left=.025,right=.985,bottom=.03,top=.96,wspace=.1)
    return fig


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
    out=ROOT/'figures';out.mkdir(exist_ok=True)
    for name,make in [('architecture_overview',overview),('architecture_detail',patch_detail)]:
        fig=make()
        for suffix in ('pdf','svg','png'):
            path=out/f'{name}.{suffix}'
            fig.savefig(path,dpi=220)
            if suffix=='svg':path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
        plt.close(fig)
    print('Wrote two architecture cartoons as PDF, editable SVG, and PNG.')


if __name__=='__main__':main()

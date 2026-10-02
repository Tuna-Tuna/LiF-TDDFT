"""Reproduce PRA figures 2–6 from packaged, unmodified author data.

Requires Python >=3.10, numpy, pandas, matplotlib and Pillow.
Run: python scripts/reproduce_pra_figures.py --data-dir /path/to/source_data
The optional --dependency-dir only supports an existing local matplotlib install.
No new simulation, fit, or probability rescaling is performed.
Display-only interpolation is used for the current-field arrows and flux markers.
"""
from pathlib import Path
import argparse, hashlib, json, sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data-dir', type=Path, required=True, help='Local author-supplied source_data directory (figure2 through figure6)')
parser.add_argument('--output', type=Path, default=Path('results/pra_figures'))
parser.add_argument('--font-family', default='Times New Roman', help='Installed font family; missing fonts fail instead of falling back')
parser.add_argument('--verify-eps', action='store_true', help='Also convert each checked EPS to PDF and render that PDF with Ghostscript')
parser.add_argument('--ghostscript', help='Ghostscript executable for --verify-eps')
parser.add_argument('--dependency-dir', type=Path)
parser.add_argument('--figures', type=int, nargs='+', choices=[2,3,4,5,6], default=[2,3,4,5,6])
parser.add_argument('--preview-dir', type=Path)
args = parser.parse_args()
if args.dependency_dir:
    sys.path.insert(0, str(args.dependency_dir))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties, findfont
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from lif_tddft.eps_export import find_ghostscript, render_eps_for_review, save_eps
DATA = args.data_dir.expanduser().resolve()
if not DATA.is_dir():
    parser.error(f'Source-data directory does not exist: {DATA}')
try:
    FONT_FILE = findfont(FontProperties(family=args.font_family), fallback_to_default=False)
except ValueError:
    parser.error(f'Font family is not installed: {args.font_family}. Install it or explicitly select --font-family.')
if args.verify_eps:
    args.ghostscript = find_ghostscript(args.ghostscript)
OUT = args.output.resolve()
OUT.mkdir(parents=True, exist_ok=True)
WIDTH = 17.2/2.54
B, O, G, C, M, Y, K = '#0072B2','#D55E00','#009E73','#56B4E9','#CC79A7','#B79B00','#171717'
SPEEDS = [.1,.15,.2,.3,.4,.5]
plt.rcParams.update({
    'font.family':args.font_family, 'font.size':11.5,
    'mathtext.fontset':'custom', 'mathtext.rm':args.font_family,
    'mathtext.it':args.font_family+':italic', 'mathtext.bf':args.font_family+':bold',
    'axes.labelsize':12,'axes.titlesize':11.5,'xtick.labelsize':11,'ytick.labelsize':11,
    'legend.fontsize':10.5,'legend.title_fontsize':11,'axes.linewidth':.9,
    'lines.linewidth':1.7,'lines.markersize':5.7,'lines.markeredgewidth':1.2,
    'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,
    'xtick.major.width':.9,'ytick.major.width':.9,'xtick.minor.width':.7,'ytick.minor.width':.7,
    'pdf.fonttype':42,'ps.fonttype':3,'svg.fonttype':'none','figure.dpi':150,
    'savefig.facecolor':'white','figure.facecolor':'white',
})
AUDIT = {'font_file':FONT_FILE,'figures':{},'numeric_checks':{}}
if (OUT/'figure_audit.json').exists() and set(args.figures)!={2,3,4,5,6}:
    AUDIT=json.loads((OUT/'figure_audit.json').read_text(encoding='utf-8'))
    AUDIT['font_file']=FONT_FILE

def read(n, f):
    return pd.read_csv(DATA/f'figure{n}'/f)

def tidy(ax, label=None, minor=True):
    if minor:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        if ax.get_yscale() != 'log':
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    if label:
        ax.text(-.19,1.05,f'({label})',transform=ax.transAxes,fontsize=13,fontweight='bold')

def four(height=15.2, left=.115, right=.975, wspace=.41, hspace=.43):
    fig, axs=plt.subplots(2,2,figsize=(WIDTH,height/2.54))
    fig.subplots_adjust(left=left,right=right,bottom=.10,top=.95,wspace=wspace,hspace=hspace)
    return fig,axs.flatten()

def checked_line(ax,x,y,**kwargs):
    xx,yy=np.asarray(x,float),np.asarray(y,float)
    valid=np.isfinite(xx)&np.isfinite(yy)
    line,=ax.plot(xx[valid],yy[valid],**kwargs)
    np.testing.assert_array_equal(line.get_xdata(),xx[valid])
    np.testing.assert_array_equal(line.get_ydata(),yy[valid])
    return line

def save(fig,n,notes):
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    outside=[]
    sizes=[]
    hidden_tick_labels=set()
    for ax in fig.axes:
        for axis in (ax.xaxis,ax.yaxis):
            drawn=set(axis._update_ticks())
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if tick not in drawn:
                    hidden_tick_labels.update([tick.label1,tick.label2])
    for txt in fig.findobj(matplotlib.text.Text):
        if not txt.get_visible() or not txt.get_text() or txt in hidden_tick_labels:continue
        box=txt.get_window_extent(renderer)
        # Tick labels outside the chosen axis limits are not actually drawn.
        if box.width and box.height and (box.x0<-.7 or box.y0<-.7 or box.x1>fig.bbox.x1+.7 or box.y1>fig.bbox.y1+.7):
            outside.append(txt.get_text())
        sizes.append(txt.get_fontsize())
    assert not outside, (n,outside)
    minsize=min(sizes)
    cap=min(TextPath((0,0),char,size=minsize,prop=FontProperties(family=args.font_family)).get_extents().height for char in 'H0123456789')
    cap_mm=cap/72*25.4*(16.51/17.2)
    assert cap_mm>=2.0,(n,minsize,cap_mm)
    for ext in ['pdf','png']:
        fig.savefig(OUT/f'figure{n}.{ext}',dpi=600)
    eps_audit = save_eps(fig, OUT/f'figure{n}.eps', dpi=600)
    if args.verify_eps:
        eps_audit = render_eps_for_review(OUT/f'figure{n}.eps', OUT/'eps_review', ghostscript=args.ghostscript)
    with Image.open(OUT/f'figure{n}.png') as im:
        im.convert('RGB').save(OUT/f'figure{n}.tif',compression='tiff_lzw',dpi=(600,600))
        dimensions=im.size
    if args.preview_dir:
        args.preview_dir.mkdir(parents=True,exist_ok=True)
        fig.savefig(args.preview_dir/f'figure{n}_preview.png',dpi=190)
    AUDIT['figures'][str(n)]={'size_cm':(fig.get_size_inches()*2.54).tolist(),
        'eps_validation':eps_audit,'pixel_dimensions':list(dimensions),'raster_dpi':600,'minimum_point_size':minsize,
        'minimum_cap_digit_height_mm_after_Word_scaling':cap_mm,
        'all_text_inside_canvas':True,'axes_bounds_figure_fraction':[list(ax.get_position().bounds) for ax in fig.axes],'notes':notes}
    plt.close(fig)

def figure2():
    fig,(a,b,c,d)=four(height=15.4,left=.13,wspace=.50)
    model=[('embedded',B,'-',r'Embedded Li$_9$F$_9$'),
           ('li25',O,'--',r'All-real Li$_{25}$F$_{25}$'),('slab',G,'-.','Periodic slab')]
    for key,color,ls,label in model:
        f=read(2,f'a_potential/panel_a_{key}.csv')
        checked_line(a,f.h_au,f.potential_eV,color=color,ls=ls,label=label)
        f=read(2,f'b_density/profile_{key}.csv')
        checked_line(b,f.h_au,f.rho_a0m3,color=color,ls=ls,label=label)
    a.set(xlabel=r'Height, $h$ ($a_0$)',ylabel='External + Hartree\npotential energy (eV)',xlim=(1.2,10),ylim=(-1.5,1.4))
    a.set_yticks([-1,-.5,0,.5,1.])
    a.legend(frameon=False,loc='upper right',fontsize=10.5,handlelength=2,borderpad=.2,labelspacing=.25)
    b.set(xlabel=r'Height, $h$ ($a_0$)',ylabel=r'Cylinder-averaged density ($a_0^{-3}$)',xlim=(1.2,10),yscale='log',ylim=(3e-10,.1))
    b.set_yticks([1e-9,1e-7,1e-5,1e-3,1e-1])
    b.legend(frameon=False,loc='upper right',fontsize=10.5,handlelength=2,borderpad=.2,labelspacing=.25)
    f=read(2,'c_radial/radial_differences.csv')
    for key,color,ls,label in model[1:]:
        checked_line(c,f.r_au,f[f'{key}_vacuum_minus_embedded_a0m3'],color=color,ls=ls,label=label)
    c.axhline(0,color=K,lw=.8)
    c.set(xlabel=r'Radial distance, $r$ ($a_0$)',ylabel=r'Density difference, $\Delta\bar{\rho}^{+}$ ($a_0^{-3}$)',xlim=(0,3),ylim=(-.25,.85))
    c.legend(frameon=False,loc='upper right',fontsize=10.5,handlelength=2)
    f=read(2,'d_integrated/panel_d.csv')
    for key,color,ls,label,marker in [('li25',O,'-',model[1][3],'o'),('slab',G,'--',model[2][3],'s')]:
        checked_line(d,f.R_au,f[f'{key}_minus_embedded_e'],color=color,ls=ls,marker=marker,mfc='white',ms=4.7,markevery=2,label=label)
    d.axvline(1.889726125,color='.35',ls=':',lw=1)
    d.text(1.925,.015,r'$R=1.89$ $a_0$',rotation=90,color='.4',fontsize=10.5)
    d.set(xlabel=r'Cylinder radius, $R$ ($a_0$)',ylabel='Integrated difference,\n'+r'$\Delta N$ (electrons)',xlim=(0,3),ylim=(-.006,.031))
    sf=ScalarFormatter(useMathText=True);sf.set_powerlimits((-2,-2));d.yaxis.set_major_formatter(sf)
    d.legend(frameon=False,loc='upper left',fontsize=10.5,handlelength=2)
    for ax,label in zip((a,b,c,d),'abcd'):tidy(ax,label)
    for ax in (a,b):ax.set_xticks([2,4,6,8,10])
    save(fig,2,['Original four-panel layout retained; no extra density inset added.',
        'Potential_eV, native density grids, vacuum-hemisphere differences and signed electron differences used without scaling.'])

def figure3():
    fig,(a,b,c,d)=four(height=15.4,left=.10,wspace=.35)
    p=read(3,'probability.csv');pc=read(3,'point_charge_distance.csv')
    colors=[G,'#B27C00',B,C,O,M];marks=['o','P','s','^','D','v'];styles=['-',(0,(3,1,1,1)),'--','-.',':',(0,(5,1,1,1))]
    for i,v in enumerate(SPEEDS):
        f=p[(p.model=='LiF(100)')&np.isclose(p.velocity_au,v)].sort_values('height_a0')
        checked_line(a,f.height_a0,f.P_det_1,color=colors[i],ls=styles[i],marker=marks[i],mfc='white',label=f'{v:.2f}')
    a.legend(title=r'$v$ (a.u.)',frameon=False,ncol=2,loc='upper right',columnspacing=.5,handletextpad=.55,handlelength=1.6,labelspacing=.3)
    a.set(xlim=(1.05,10.15),ylim=(-.02,.96),xlabel=r'$h$ ($a_0$)',ylabel=r'$P_{\rm det}$')
    a.set_yticks([0,.2,.4,.6,.8]);a.set_xticks([2,4,6,8,10])
    f=p[(p.model=='LiF(100)')&np.isclose(p.velocity_au,.3)&(p.height_a0>=2.5)].sort_values('height_a0')
    checked_line(b,f.height_a0,f.P_det_1,color=B,marker='o',mfc='white',label='Explicit electrons')
    checked_line(b,pc.height_a0,pc.P_det_1,color=O,marker='s',mfc='white',ls='--',label='Point charge')
    interp=pc[pc.node_type=='interpolated']
    b.legend(frameon=False,loc='upper right',handlelength=1.7)
    b.set(xlim=(2,10.15),ylim=(-.008,.185),xlabel=r'$h$ ($a_0$)',ylabel=r'$P_{\rm det}$')
    b.set_yticks([0,.05,.10,.15]);b.set_xticks([2,4,6,8,10])
    for h,i in zip([2.5,3,3.5,4,5],[0,2,3,4,5]):
        f=p[(p.model=='LiF(100)')&np.isclose(p.height_a0,h)].sort_values('velocity_au')
        checked_line(c,f.velocity_au,f.P_det_1,color=colors[i],ls=styles[i],marker=marks[i],mfc='white',label=f'{h:g}')
    c.legend(title=r'$h$ ($a_0$)',frameon=False,ncol=2,loc='upper left',columnspacing=.6,handlelength=1.7,labelspacing=.3)
    for model,color,mark,ls in [('LiF(100)',B,'o','-'),('Point charge',O,'s','--')]:
        f=p[(p.model==model)&np.isclose(p.height_a0,3.5)].sort_values('velocity_au')
        checked_line(d,f.velocity_au,f.P_det_1,color=color,marker=mark,ls=ls,mfc='white')
    for ax in (c,d):
        ax.set(xlabel=r'$v$ (a.u.)',ylabel=r'$P_{\rm det}$',xlim=(.085,.515));ax.set_xticks([.1,.2,.3,.4,.5])
    c.set_ylim(-.02,.425);d.set_ylim(-.008,.147)
    for ax,label in zip((a,b,c,d),'abcd'):tidy(ax,label)
    AUDIT['numeric_checks']['figure3']={'author_confirmation':'Source provenance annotations declared outdated by the author; only numeric plotting columns used.',
        'near_surface_nodes':p[p.height_a0<2][['height_a0','velocity_au','P_det_1']].to_dict('records'),
        'point_charge_h3_linear_check':bool(np.isclose(interp.P_det_1.iloc[0],np.mean(pc.loc[pc.height_a0.isin([2.5,3.5]),'P_det_1']),atol=1e-16))}
    save(fig,3,['Author confirmed outdated provenance annotations should be ignored; numeric columns and original open-marker styles preserved.'])

def figure4():
    fig,(a,b,c,d)=four(height=15.2,left=.12,wspace=.43)
    p=read(4,'merged.csv')
    for key,color,mark,ls,label in [('sigma_max',B,'o','-',r'$\sigma_{\max}$'),('SAB_frobenius',O,'s','--',r'$\Vert S_{AB}\Vert_F$')]:
        f=p[p.series_key==key].sort_values('RFF_a0');checked_line(a,f.RFF_a0,f.value,color=color,marker=mark,ls=ls,mfc='white',label=label)
    a.legend(frameon=False,loc='upper right',fontsize=12)
    for ax,key,color,mark in [(b,'lowdin_half_L1_e',G,'o'),(c,'vW_delta_cutoff_1e10_eV',O,'D')]:
        f=p[p.series_key==key].sort_values('RFF_a0');checked_line(ax,f.RFF_a0,f.value,color=color,marker=mark,mfc='white')
    a.set_ylabel('Occupied-space overlap');b.set_ylabel(r'$Q_{\rm orth}$ (electrons)');c.set_ylabel(r'$\Delta T_{\rm vW}$ (eV)')
    for ax in (a,b,c):
        ax.set(xlabel=r'$R_{FF}$ ($a_0$)',xlim=(3.2,19.2));ax.set_xticks([5,10,15])
    a.set_ylim(-.015,.33);b.set_ylim(-.008,.175);c.set_ylim(-.45,9.25)
    for i,(key,color,label) in enumerate([('orth_half_L1_e',G,r'$\rho_L-\rho_0$'),('relax_half_L1_e',B,r'$\rho_{\rm SCF}-\rho_L$'),('total_half_L1_e',O,r'$\rho_{\rm SCF}-\rho_0$')]):
        f=p[p.series_key==key].sort_values('plot_order')
        d.bar(np.arange(3)+(i-1)*.24,f.value,width=.235,color=color,label=label)
    d.set_xticks(range(3),['6.143','4.347','3.500']);d.set(xlabel=r'$R_{FF}$ ($a_0$)',ylabel=r'Half-$L^1$ difference (electrons)',ylim=(0,.57))
    d.legend(frameon=False,loc='upper left',fontsize=11.5,labelspacing=.1,handlelength=1.8)
    for ax,label in zip((a,b,c,d),'abcd'):tidy(ax,label)
    d.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    save(fig,4,['All 53 supplied rows used; actual RFF_a0 values used, with rounded categorical labels only in panel d.'])

def interp2(x,z,f,xq,zq):
    assert xq.min()>=x[0] and xq.max()<=x[-1] and zq.min()>=z[0] and zq.max()<=z[-1]
    ix=np.clip(np.searchsorted(x,xq)-1,0,len(x)-2);iz=np.clip(np.searchsorted(z,zq)-1,0,len(z)-2)
    tx=(xq-x[ix])/(x[ix+1]-x[ix]);tz=(zq-z[iz])/(z[iz+1]-z[iz])
    return (1-tz)*((1-tx)*f[iz,ix]+tx*f[iz,ix+1])+tz*((1-tx)*f[iz+1,ix]+tx*f[iz+1,ix+1])

def figure5():
    meta=json.loads((DATA/'figure5/radial_current_manifest.json').read_text(encoding='utf-8'))
    figure_height_cm=13.0
    frame_width=.335
    frame_width_cm=frame_width*17.2
    reference_ratio=(17.2*(.975-.12)/(2+.43))/(15.2*(.95-.10)/(2+.43))
    top_height_cm=frame_width_cm*14.5/21
    bottom_height_cm=frame_width_cm/reference_ratio
    bottom_y=1.35/figure_height_cm
    top_y=(1.35+bottom_height_cm+1.75)/figure_height_cm
    top_height=top_height_cm/figure_height_cm
    bottom_height=bottom_height_cm/figure_height_cm
    fig=plt.figure(figsize=(WIDTH,figure_height_cm/2.54))
    a=fig.add_axes([.09,top_y,frame_width,top_height]);b=fig.add_axes([.55,top_y,frame_width,top_height])
    c=fig.add_axes([.09,bottom_y,frame_width,bottom_height]);d=fig.add_axes([.55,bottom_y,frame_width,bottom_height])
    cmap=LinearSegmentedColormap.from_list('inward_white_outward',[B,'#FBFBF9',O])
    scale=meta['color_scale'];norm=SymLogNorm(linthresh=scale['linear_threshold'],linscale=scale['linear_scale'],vmin=-scale['symmetric_limit'],vmax=scale['symmetric_limit'],base=10)
    spatial_checks=[]
    for ax,panel,rec in zip([a,b],'ab',meta['panels']):
        with np.load(DATA/f'figure5/panel_{panel}_radial_field.npz') as src:f={k:src[k] for k in src.files}
        x,z=f['x_active_a0'],f['z_a0'];xx,zz=np.meshgrid(x,z)
        xp,zp=rec['projectile_x_active_a0'],rec['projectile_z_a0'];radius=np.hypot(xx-xp,zz-zp)
        jx,jz=f['jx_e_a0m2_fsm1'],f['jz_e_a0m2_fsm1']
        radial_check=np.divide(jx*(xx-xp)+jz*(zz-zp),radius,out=np.zeros_like(radius),where=radius>1e-12)
        np.testing.assert_array_equal(radial_check,f['jr_e_a0m2_fsm1'])
        mesh=ax.pcolormesh(xx,zz,np.ma.masked_where(radius<.12,f['jr_e_a0m2_fsm1']),cmap=cmap,norm=norm,shading='nearest',rasterized=True)
        # Deliberately no contour layer: the author requested its removal.
        xs,zs=np.meshgrid(np.arange(-13.5,7,.95),np.arange(-3.1,10.8,.95))
        sx,sz=interp2(x,z,jx,xs,zs),interp2(x,z,jz,xs,zs);mag=np.hypot(sx,sz)
        mask=(mag>=.0015)&(np.hypot(xs-xp,zs-zp)>=.65)
        length=.25*np.arcsinh(mag/.002)
        u=np.divide(sx,mag,out=np.zeros_like(sx),where=mag>0)*length;w=np.divide(sz,mag,out=np.zeros_like(sz),where=mag>0)*length
        assert int(mask.sum())==rec['arrows_displayed']
        ax.quiver(xs[mask],zs[mask],u[mask],w[mask],color='.12',angles='xy',scale_units='xy',scale=1,width=.0044,pivot='mid',headwidth=3.8,headlength=4.5,zorder=4)
        ax.add_patch(Circle((xp,zp),meta['sphere_radius_a0'],fill=False,color='.42',lw=1.05,ls=(0,(4,3)),zorder=5))
        ax.scatter(0,0,s=46,marker='D',c='#E69F00',edgecolors='black',linewidths=1.1,zorder=8)
        ax.scatter(xp,zp,s=46,c='#00BFC4',edgecolors='black',linewidths=1.1,zorder=8)
        ax.axhline(0,color='.55',lw=.75,ls='--');ax.axvline(0,color='.65',lw=.75,ls=':')
        # Restore the original visible range; the frame is rectangular, preserving x/z unit scale.
        ax.set(xlim=(-14,7),ylim=(-3.5,11),xlabel=r'$x$ ($a_0$)',ylabel=r'$z$ ($a_0$)')
        ax.set_aspect('equal',adjustable='box');ax.set_xticks([-12,-6,0,6]);ax.set_yticks([0,4,8])
        ax.text(.025,.97,r'$t-t_F=$'+f"{rec['time_minus_force_fs']:.2f} fs",va='top',transform=ax.transAxes,fontsize=10.5,zorder=10,bbox={'facecolor':'white','edgecolor':'none','pad':1.2})
        spatial_checks.append({'panel':panel,'radial_projection_exact':True,'arrow_count':int(mask.sum()),'density_contours_removed':True})
    cbax=fig.add_axes([.91,top_y,.016,top_height]);cb=fig.colorbar(mesh,cax=cbax)
    cb.set_ticks([-1,-.1,-.01,0,.01,.1,1]);cb.set_ticklabels(['−1','−0.1','−0.01','0','0.01','0.1','1'])
    cb.ax.tick_params(labelsize=10.5,length=2.8);cb.ax.set_title(r'$\Delta j_r$',fontsize=12,pad=7)
    fig.text(.928,top_y-.78/figure_height_cm,r'$e/(a_0^2\,\mathrm{fs})$',ha='center',fontsize=10.5)
    flux=read(5,'panel_c_signed_flux_v030.csv').sort_values('t_minus_force_fs');phases=read(5,'selected_phases.csv');balance=read(5,'panel_d_balance.csv')
    t,phi=flux.t_minus_force_fs.to_numpy(),flux.phi_direct_e_fs.to_numpy()
    checked_line(c,t,phi,color='.15',lw=1.7)
    c.fill_between(t,0,phi,where=phi>=0,interpolate=True,color='#F3C5A8')
    c.fill_between(t,0,phi,where=phi<=0,interpolate=True,color='#A6CEE1')
    c.axhline(0,color='.3',lw=.9)
    for label,(_,rec),color in zip('ab',phases.iterrows(),[O,B]):
        tx=rec.t_minus_force_fs;value=float(np.interp(tx,t,phi))
        c.axvline(tx,color='.55',lw=.9,ls='--');c.plot(tx,value,'o',color=color,mfc='white',ms=6,zorder=5)
        c.annotate(f'({label})',(tx,value),xytext=(5,3 if value>0 else -12),textcoords='offset points',color=color,fontsize=11.5)
    c.set(xlabel=r'$t-t_F$ (fs)',ylabel=r'$\Delta\Phi_R$ (electrons/fs)',xlim=(t.min(),1.2),ylim=(-.35,.17))
    c.set_xticks([-.4,0,.4,.8,1.2]);c.set_yticks([-.3,-.2,-.1,0,.1])
    c.legend(handles=[Patch(facecolor='#F3C5A8',label='Outward (+)'),Patch(facecolor='#A6CEE1',label='Return (−)')],frameon=False,loc='lower left',handlelength=1.0,fontsize=10.5,labelspacing=.25)
    for col,label,color,marker,ls in [('Outward_e','Outward',O,'o','-'),('Return_e','Return',B,'s','--'),('Net_e','Net',M,'D','-.')]:
        checked_line(d,balance.velocity_au,balance[col],label=label,color=color,marker=marker,ls=ls,mfc='white')
    np.testing.assert_allclose(balance.Outward_e-balance.Return_e,balance.Net_e,atol=1e-15,rtol=0)
    d.set(xlabel=r'$v$ (a.u.)',ylabel='Cumulative transfer\n(electrons)',xlim=(.075,.525),ylim=(0,.235))
    d.set_xticks([.1,.2,.3,.4,.5]);d.set_yticks([0,.05,.10,.15,.20]);d.set_yticklabels(['0','0.05','0.10','0.15','0.20'])
    d.legend(frameon=False,ncol=2,loc='upper left',columnspacing=.6,handlelength=1.6,handletextpad=.4,fontsize=10.5,labelspacing=.3)
    for ax,label in zip([a,b,c,d],'abcd'):tidy(ax,label)
    AUDIT['numeric_checks']['figure5']=spatial_checks
    fig.canvas.draw()
    bounds=np.array([ax.get_position().bounds for ax in [a,b,c,d]])
    np.testing.assert_allclose(bounds[:,2],bounds[0,2],rtol=0,atol=1e-12)
    np.testing.assert_allclose(bounds[[0,1],3],bounds[0,3],rtol=0,atol=1e-12)
    np.testing.assert_allclose(bounds[[2,3],3],bounds[2,3],rtol=0,atol=1e-12)
    np.testing.assert_allclose(bounds[[0,2],0],bounds[0,0],rtol=0,atol=1e-12)
    np.testing.assert_allclose(bounds[[1,3],0],bounds[1,0],rtol=0,atol=1e-12)
    np.testing.assert_allclose(bounds[[0,1],1],bounds[0,1],rtol=0,atol=1e-12)
    np.testing.assert_allclose(bounds[[2,3],1],bounds[2,1],rtol=0,atol=1e-12)
    assert a.get_aspect()==b.get_aspect()==1.0
    physical_bounds=np.array([ax.get_window_extent().bounds[2:] for ax in [a,b,c,d]])/fig.dpi
    np.testing.assert_allclose(physical_bounds[:2,0]/physical_bounds[:2,1],21/14.5,rtol=0,atol=1e-12)
    np.testing.assert_allclose(physical_bounds[2:,0]/physical_bounds[2:,1],reference_ratio,rtol=0,atol=1e-12)
    save(fig,5,['Radial-current NPZ arrays from the adopted figure chain, not the older density-colormap candidate.',
        'Symmetric-log color scale, arrow locations/directions/lengths, sphere and times preserved; density contours removed.',
        'Local flux not rescaled to full-window totals; balance Outward − Return = Net verified.',
        'Panels a,b restore x=[-14,7], z=[-3.5,11] a0, with equal physical x/z scale and no blank band outside the field grid. Panels c,d match the width/height ratio of Figure4; all column edges align. Data arrays unchanged; no contours or c,d titles.'])

def figure6():
    fig,(a,b,c,d)=four(height=15.2,left=.12,wspace=.43)
    reference_ratio=(17.2*(.975-.12)/(2+.43))/(15.2*(.95-.10)/(2+.43))
    force=read(6,'Fig6a_force.csv');track=read(6,'Fig6b_trajectories.csv');capture=read(6,'Fig6c_capture.csv')
    colors=[K,O,B,G,M,Y];styles=['-','--','-.',':',(0,(4,2)),(0,(5,1,1,1))]
    a.plot([.375,.51],[.375,.51],color='.5',ls='--',lw=.9)
    geometries=['Overlap onset','Intermediate overlap','Near force maximum']
    for i,v in enumerate(SPEEDS):
        f=force[np.isclose(force.v_au,v)]
        for geometry,mark in zip(geometries,['o','s','D']):
            s=f[f.geometry==geometry]
            assert len(s)==1,(geometry,s)
            a.plot(s.TD_Fz_au_X,s.Static_Fz_au_Y,ls='none',marker=mark,mfc='none',color=colors[i],ms=5.8,mew=1.2)
    leg=[Line2D([],[],color=colors[i],marker='o',mfc='white',ls='none',label=f'{v:.2f}') for i,v in enumerate(SPEEDS)]
    a.legend(handles=leg,title=r'$v$ (a.u.)',frameon=False,ncol=3,loc='upper left',columnspacing=.6,handlelength=.6,handletextpad=.35,fontsize=10.5,borderpad=.2,labelspacing=.3)
    label_box={'facecolor':'white','edgecolor':'none','pad':.7}
    a.text(.378,.397,'Overlap onset',fontsize=10.5,bbox=label_box)
    a.text(.390,.434,'Intermediate overlap',fontsize=10.5,bbox=label_box)
    a.text(.482,.456,'Near force\nmaximum',fontsize=10.5,ha='center',bbox=label_box)
    a.set(xlabel=r'TD $F_z$ (a.u.)',ylabel=r'Static $F_z$ (a.u.)',xlim=(.375,.51),ylim=(.375,.51))
    a.set_xticks([.38,.42,.46,.50]);a.set_yticks([.38,.42,.46,.50])
    for i,v in enumerate(SPEEDS):
        vs=f'{v:.2f}'
        checked_line(b,track[f'x_v{vs}_a0_X'],track[f'h_v{vs}_a0_Y'],color=colors[i],ls=styles[i],label=vs)
        checked_line(c,capture[f'h_v{vs}_a0_X'],capture[f'Pcap_v{vs}_Y'],color=colors[i],ls=styles[i],label=vs)
    b.set(xlabel=r'Distance along surface ($a_0$)',ylabel=r'$h$ ($a_0$)',xlim=(0,1140),ylim=(1,10.5));b.set_xticks([0,300,600,900])
    b.legend(title=r'$v$ (a.u.)',frameon=False,ncol=2,loc='upper left',bbox_to_anchor=(.12,1.015),columnspacing=.35,handlelength=.85,handletextpad=.35,labelspacing=.3)
    c.set(xlabel=r'$h$ ($a_0$)',ylabel=r'Capture probability, $P_{\rm cap}$',xlim=(1.2,10),ylim=(0,.4));c.set_xticks([2,4,6,8,10]);c.set_yticks([0,.1,.2,.3,.4])
    c.legend(title=r'$v$ (a.u.)',frameon=False,ncol=2,loc='upper right',columnspacing=.7,handlelength=1.7,labelspacing=.3)
    p=read(6,'Fig6d_force_sensitivity.csv');no=read(6,'Fig6d_no_detachment.csv');exp=read(6,'Fig6d_experiment.csv')
    spec=[('baseline_fraction','Theory (1.0, 1.0)',B,'s','-'),('low_070_high_100_fraction','Low 0.7 (0.7, 1.0)',G,'D','--'),
          ('low_100_high_070_fraction','High 0.7 (1.0, 0.7)',M,'^',':'),('low_100_high_130_fraction','High 1.3 (1.0, 1.3)',Y,'v','-.'),
          ('low_070_high_130_fraction','Both (0.7, 1.3)',K,'x','--')]
    for j,(key,label,color,mark,ls) in enumerate(spec):
        checked_line(d,p.v_au,p[key],color=color,marker=mark,ls=ls,mfc=color if j==0 else 'white',label=label,ms=5.2,lw=1.7 if j==0 else 1.25)
    checked_line(d,no.v_au_X,no.No_detachment_Fminus_fraction_Y,color='.35',marker='^',mfc='white',ls='--',label='No detachment',ms=5.3)
    d.errorbar(exp.v_au_X,exp.Experiment_Fminus_fraction_Y,yerr=np.vstack([exp.Y_error_minus,exp.Y_error_plus]),color=O,marker='o',mfc='white',ls='-',lw=1.25,elinewidth=.8,capsize=2.5,capthick=.8,label='Experiment [21]',ms=5.2)
    d.set(xlabel=r'$v$ (a.u.)',ylabel=r'F$^{-}$ fraction',xlim=(.05,.525),ylim=(0,1.05));d.set_xticks([.1,.2,.3,.4,.5])
    legend_handles,legend_labels=d.get_legend_handles_labels()
    legend_order=[5,0,1,2,3,4,6]
    d.legend([legend_handles[i] for i in legend_order],[legend_labels[i] for i in legend_order],frameon=False,loc='lower left',bbox_to_anchor=(.155,.01),fontsize=10.5,handlelength=1.0,handletextpad=.3,labelspacing=0,borderpad=0,borderaxespad=.1)
    for ax,label in zip([a,b,c,d],'abcd'):tidy(ax,label)
    for ax in [a,b,c,d]:ax.set_box_aspect(1/reference_ratio)
    fig.canvas.draw()
    physical_bounds=np.array([ax.get_window_extent().bounds[2:] for ax in [a,b,c,d]])/fig.dpi
    np.testing.assert_allclose(physical_bounds[:,0]/physical_bounds[:,1],reference_ratio,rtol=0,atol=1e-12)
    AUDIT['numeric_checks']['figure6']={'force_max_absolute_relative_deviation_percent':float((abs(force.TD_Fz_au_X-force.Static_Fz_au_Y)/abs(force.Static_Fz_au_Y)*100).max()),
        'experiment_points':len(exp),'asymmetric_error_lengths_preserved':True,'force_scenarios':5,
        'd_source':'Author-supplied final fractions; source README excluded by explicit author correction; no new calculation performed.'}
    save(fig,6,['Force/trajectory/capture arrays preserved; all 5 force-sensitivity series, no-detachment reference, 13 experimental points and asymmetric error bars included.',
        'Only author-supplied numeric plotting columns used; no new physical calculation performed.',
        'All four axis frames match Figure4 width/height ratio. Panel d legend remains inside its frame below the curves; panel b velocity heading remains an in-frame legend title.'])

for n in args.figures:
    fn=globals()[f'figure{n}']
    fn()
    print(fn.__name__+' complete',flush=True)
AUDIT['source_hashes']={str(p.relative_to(DATA)):hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.rglob('*') if p.is_file()}
AUDIT['output_hashes']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.name!='figure_audit.json'}
(OUT/'figure_audit.json').write_text(json.dumps(AUDIT,indent=2,ensure_ascii=False),encoding='utf-8')
print('All figures and numerical checks complete.',flush=True)

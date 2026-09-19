import numpy as np, pickle
from stripping import *

def run_case(name, dt=10., nl=41, nn=240, **kw):
    cfg = Config(**kw)
    res = classify(run(cfg, n_lat=nl, n_lon=nn, dt=dt, r_stop=8.0, verbose=False))
    res['name'] = name
    return res

def stats(res):
    w, cat = res['wmass'], res['cat']
    lifted = cat != 'non_leve'
    win = cat == 'fenetre_lunaire'
    out = dict(name=res['name'],
               leve=w[lifted].sum(),
               fen=w[win].sum(),
               fen_rel=w[win].sum()/max(w[lifted].sum(), 1e-30),
               embB=w[cat=='emportee_B'].sum(),
               embB_rel=w[cat=='emportee_B'].sum()/max(w[lifted].sum(),1e-30))
    # statistique robuste sur les particules bornees issues des basses latitudes
    b = np.isin(cat, ['fenetre_lunaire','disque_interne']) & (np.abs(res['lat0'])<np.radians(5))
    out['rcirc_med'] = float(np.nanmedian(res['rcirc'][b])) if b.any() else np.nan
    out['rcirc_p90'] = float(np.nanpercentile(res['rcirc'][b], 90)) if b.any() else np.nan
    out['L_med'] = np.sqrt(out['rcirc_med']) if b.any() else np.nan
    out['L_p90'] = np.sqrt(out['rcirc_p90']) if b.any() else np.nan
    # part des particules "fenetre" ayant frole B (< 1.5 R_terre)
    if win.any():
        out['fen_close'] = float((res['dminB'][win] < 1.5*RE).mean())
        out['fen_qperi'] = float(np.nanmedian(res['qperi'][win]))
    else:
        out['fen_close'] = np.nan; out['fen_qperi'] = np.nan
    return out

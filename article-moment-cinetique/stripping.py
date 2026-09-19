"""
Arrachement de matiere depuis la bande intertropicale d'une proto-Terre en
rotation rapide, lors du passage rasant non collisionnel d'un corps B.

Modele :
  - proto-Terre = modele de Roche (masse centralement condensee, J2 = 0),
    surface = equipotentielle Phi = -GM/r - 0.5 Omega^2 r^2 cos^2(lat)
  - B suit la vraie conique hyperbolique (mu_tot = G(M + MB))
  - champ differentiel exact (terme direct + terme indirect), sans
    developpement multipolaire
  - particules-test lachees depuis la surface, en corotation, au moment ou
    la composante verticale de (gravite + centrifuge + maree) devient positive

M. Debailleul -- verification numerique, 2026
"""

import numpy as np
from scipy.optimize import brentq

# ---------------------------------------------------------------- constantes
G     = 6.67430e-11
ME    = 5.9722e24
RE    = 6.371e6
ML    = 7.342e22
LEM   = 3.500e34          # moment cinetique du systeme Terre-Lune actuel
GM    = G * ME
AROCHE = 2.903            # en R_terre


class Config:
    def __init__(self, wtilde=0.97, a=1.70, muB=0.95, b=3.0,
                 iB_deg=25.0, omegaB_deg=90.0, vinf=2.0e3,
                 retrograde=False):
        self.wtilde = wtilde
        self.a      = a                      # Req / RE
        self.muB    = muB
        self.b      = b                      # qB / RE
        self.vinf   = vinf
        self.Req    = a * RE
        self.qB     = b * RE
        self.MB     = muB * ME
        self.mu_tot = G * (ME + self.MB)
        self.Omega  = wtilde * np.sqrt(GM / self.Req**3)
        self.iB     = np.radians(180.0 - iB_deg if retrograde else iB_deg)
        self.omegaB = np.radians(omegaB_deg)
        self.retro  = retrograde
        self.eB     = 1.0 + self.qB * vinf**2 / self.mu_tot
        # constante de l'equipotentielle de surface
        self.C      = -GM / self.Req - 0.5 * self.Omega**2 * self.Req**2

    def __repr__(self):
        return (f"<Config w={self.wtilde} a={self.a} muB={self.muB} b={self.b} "
                f"iB={np.degrees(self.iB):.0f} retro={self.retro} eB={self.eB:.3f}>")


# ------------------------------------------------------- surface de Roche
def surface_radius(cfg, lat):
    """rayon de l'equipotentielle de surface a la latitude lat (rad)."""
    lat = np.atleast_1d(lat)
    out = np.empty_like(lat, dtype=float)
    for k, la in enumerate(lat):
        c2 = np.cos(la)**2
        f = lambda r: -GM / r - 0.5 * cfg.Omega**2 * r**2 * c2 - cfg.C
        lo = 0.5 * GM / (-cfg.C)          # borne basse sure
        out[k] = brentq(f, lo, cfg.Req * (1 + 1e-12))
    return out


def eff_gravity(cfg, r):
    """acceleration effective (gravite + centrifuge) dans le repere tournant.
    r : (N,3). retourne (N,3)."""
    rn = np.linalg.norm(r, axis=1, keepdims=True)
    g = -GM * r / rn**3
    cent = np.zeros_like(r)
    cent[:, 0] = cfg.Omega**2 * r[:, 0]
    cent[:, 1] = cfg.Omega**2 * r[:, 1]
    return g + cent


# ---------------------------------------------------- trajectoire de B
def rot_matrix(cfg):
    """plan orbital de B -> repere equatorial (noeud ascendant sur Ox)."""
    i, w = cfg.iB, cfg.omegaB
    Rw = np.array([[np.cos(w), -np.sin(w), 0],
                   [np.sin(w),  np.cos(w), 0],
                   [0, 0, 1]])
    Ri = np.array([[1, 0, 0],
                   [0, np.cos(i), -np.sin(i)],
                   [0, np.sin(i),  np.cos(i)]])
    return Ri @ Rw


def B_state(cfg, r_start_factor=10.0):
    """etat initial de B, entrant, a r = r_start_factor * qB."""
    e, q, mu = cfg.eB, cfg.qB, cfg.mu_tot
    p = q * (1 + e)
    r0 = r_start_factor * q
    cosnu = (p / r0 - 1.0) / e
    cosnu = np.clip(cosnu, -1.0, 1.0)
    nu = -np.arccos(cosnu)                      # branche entrante
    h = np.sqrt(mu * p)
    # position / vitesse dans le plan orbital (periapse sur Ox)
    rv = np.array([r0 * np.cos(nu), r0 * np.sin(nu), 0.0])
    vv = (mu / h) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])
    M = rot_matrix(cfg)
    return M @ rv, M @ vv


def accel_B(cfg, rB):
    return -cfg.mu_tot * rB / np.linalg.norm(rB)**3


def tidal_field(cfg, r, rB):
    """champ differentiel exact de B sur les particules (terme direct + indirect)."""
    d = r - rB
    dn = np.linalg.norm(d, axis=1, keepdims=True)
    rBn = np.linalg.norm(rB)
    return -G * cfg.MB * (d / dn**3 + rB / rBn**3)


# ---------------------------------------------------------- particules
def seed_particles(cfg, n_lat=41, n_lon=240, lat_max_deg=26.0):
    lats = np.radians(np.linspace(-lat_max_deg, lat_max_deg, n_lat))
    lons = np.linspace(0.0, 2 * np.pi, n_lon, endpoint=False)
    rs = surface_radius(cfg, lats)
    LA, LO = np.meshgrid(lats, lons, indexing='ij')
    R = np.repeat(rs[:, None], n_lon, axis=1)
    x = R * np.cos(LA) * np.cos(LO)
    y = R * np.cos(LA) * np.sin(LO)
    z = R * np.sin(LA)
    pos = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
    # poids de masse ~ surface elementaire cos(lat)
    w = (np.cos(LA)).ravel()
    return pos, LA.ravel(), LO.ravel(), w / w.sum()


def corotation_velocity(cfg, r):
    v = np.zeros_like(r)
    v[:, 0] = -cfg.Omega * r[:, 1]
    v[:, 1] = +cfg.Omega * r[:, 0]
    return v


def spin_rotate(r0, Omega, t):
    c, s = np.cos(Omega * t), np.sin(Omega * t)
    out = np.empty_like(r0)
    out[:, 0] = c * r0[:, 0] - s * r0[:, 1]
    out[:, 1] = s * r0[:, 0] + c * r0[:, 1]
    out[:, 2] = r0[:, 2]
    return out


# ------------------------------------------------------------ integration
def run(cfg, n_lat=41, n_lon=240, dt=10.0, r_stop=10.0, verbose=True,
        lat_max_deg=26.0, seed_pos=None, checkpoints=()):
    """integre le passage. retourne un dict de resultats par particule."""
    if seed_pos is None:
        pos0, lat0, lon0, wmass = seed_particles(cfg, n_lat, n_lon, lat_max_deg)
    else:
        pos0, lat0, lon0, wmass = seed_pos
    N = len(pos0)

    rB, vB = B_state(cfg, r_stop)
    t = 0.0

    released = np.zeros(N, dtype=bool)
    dead = np.zeros(N, dtype=bool)          # retombee sur la proto-Terre
    t_rel = np.full(N, np.nan)
    lat_rel = np.full(N, np.nan)
    r = np.zeros((N, 3))
    v = np.zeros((N, 3))
    dminB = np.full(N, np.inf)

    # rayon de surface par latitude, pour le test de retombee
    lat_grid = np.radians(np.linspace(-90, 90, 181))
    rsurf_grid = surface_radius(cfg, lat_grid)

    def rsurf_of(rvec):
        la = np.arcsin(np.clip(rvec[:, 2] / np.linalg.norm(rvec, axis=1), -1, 1))
        return np.interp(la, lat_grid, rsurf_grid)

    def deriv(rr, vv, rBc):
        a = -GM * rr / np.linalg.norm(rr, axis=1, keepdims=True)**3
        a += tidal_field(cfg, rr, rBc)
        return vv, a

    snaps = {}
    cp = sorted(checkpoints)
    icp = 0
    past_peri = False

    nstep = 0
    while True:
        # --- position des particules encore attachees (corotation) ---
        att = ~released & ~dead
        if att.any():
            rat = spin_rotate(pos0[att], cfg.Omega, t)
            gvec = eff_gravity(cfg, rat)                 # gravite + centrifuge
            nhat = -gvec / np.linalg.norm(gvec, axis=1, keepdims=True)
            atide = tidal_field(cfg, rat, rB)
            lift = np.einsum('ij,ij->i', gvec + atide, nhat) > 0.0
            if lift.any():
                idx = np.where(att)[0][lift]
                r[idx] = rat[lift]
                v[idx] = corotation_velocity(cfg, rat[lift])
                released[idx] = True
                t_rel[idx] = t
                lat_rel[idx] = np.arcsin(rat[lift][:, 2] /
                                         np.linalg.norm(rat[lift], axis=1))

        # --- RK4 sur les particules libres ---
        live = released & ~dead
        if live.any():
            rr, vv = r[live], v[live]
            k1r, k1v = deriv(rr, vv, rB)
            rB2 = rB + 0.5 * dt * vB
            vB2 = vB + 0.5 * dt * accel_B(cfg, rB)
            k2r, k2v = deriv(rr + 0.5 * dt * k1r, vv + 0.5 * dt * k1v, rB2)
            k3r, k3v = deriv(rr + 0.5 * dt * k2r, vv + 0.5 * dt * k2v, rB2)
            rB3 = rB + dt * vB2
            k4r, k4v = deriv(rr + dt * k3r, vv + dt * k3v, rB3)
            r[live] = rr + dt / 6 * (k1r + 2 * k2r + 2 * k3r + k4r)
            v[live] = vv + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)

        # --- B (RK4 two-body) ---
        k1r, k1v = vB, accel_B(cfg, rB)
        k2r, k2v = vB + .5*dt*k1v, accel_B(cfg, rB + .5*dt*k1r)
        k3r, k3v = vB + .5*dt*k2v, accel_B(cfg, rB + .5*dt*k2r)
        k4r, k4v = vB + dt*k3v,    accel_B(cfg, rB + dt*k3r)
        rB = rB + dt/6*(k1r + 2*k2r + 2*k3r + k4r)
        vB = vB + dt/6*(k1v + 2*k2v + 2*k3v + k4v)
        t += dt
        nstep += 1

        # --- distance minimale a B ---
        live = released & ~dead
        if live.any():
            dd = np.linalg.norm(r[live] - rB, axis=1)
            idx = np.where(live)[0]
            dminB[idx] = np.minimum(dminB[idx], dd)

        # --- retombee ---
        live = released & ~dead
        if live.any():
            rn = np.linalg.norm(r[live], axis=1)
            hit = rn < rsurf_of(r[live])
            if hit.any():
                dead[np.where(live)[0][hit]] = True

        dB = np.linalg.norm(rB)
        if not past_peri and np.dot(rB, vB) > 0:
            past_peri = True
        if past_peri and icp < len(cp) and dB > cp[icp] * cfg.qB:
            live_now = released & ~dead
            rn = np.linalg.norm(r, axis=1)
            eps = 0.5*np.sum(v**2, axis=1) - GM/rn
            jv = np.cross(r, v)
            jj = np.linalg.norm(jv, axis=1)
            snaps[cp[icp]] = dict(eps=eps.copy(), rcirc=(jj**2/GM/RE),
                                  live=live_now.copy(), dead=dead.copy())
            icp += 1
        if past_peri and dB > r_stop * cfg.qB and t > 0:
            break
        if nstep > 200000:
            break

    if verbose:
        print(f"  {nstep} pas, t = {t/3600:.1f} h, laches = {released.sum()}/{N}, "
              f"retombees en vol = {dead.sum()}")

    return dict(snaps=snaps, dminB=dminB, cfg=cfg, r=r, v=v, rB=rB, vB=vB, released=released, dead=dead,
                t_rel=t_rel, lat0=lat0, lon0=lon0, lat_rel=lat_rel,
                wmass=wmass, pos0=pos0, t_end=t)


# ------------------------------------------------------------ classement
def classify(res):
    cfg = res['cfg']
    N = len(res['lat0'])
    cat = np.array(['non_leve'] * N, dtype=object)
    rcirc = np.full(N, np.nan)
    jz = np.full(N, np.nan)
    qperi = np.full(N, np.nan)

    rel = res['released']
    cat[rel & res['dead']] = 'retombee'

    live = rel & ~res['dead']
    if live.any():
        r, v = res['r'][live], res['v'][live]
        rn = np.linalg.norm(r, axis=1)
        eps = 0.5 * np.sum(v**2, axis=1) - GM / rn
        jvec = np.cross(r, v)
        j = np.linalg.norm(jvec, axis=1)
        rc = j**2 / GM
        # lie a B ?
        dr = r - res['rB']
        dv = v - res['vB']
        epsB = 0.5 * np.sum(dv**2, axis=1) - G * cfg.MB / np.linalg.norm(dr, axis=1)
        withB = epsB < 0

        e = np.sqrt(np.maximum(0.0, 1 + 2 * eps * j**2 / GM**2))
        qp = np.where(eps < 0, rc / (1 + e), np.nan)

        lab = np.empty(live.sum(), dtype=object)
        lab[:] = 'disque_interne'
        lab[(eps < 0) & (rc >= AROCHE * RE)] = 'fenetre_lunaire'
        lab[(eps >= 0) & ~withB] = 'echappee'
        lab[withB] = 'emportee_B'
        cat[np.where(live)[0]] = lab
        rcirc[np.where(live)[0]] = rc / RE
        jz[np.where(live)[0]] = jvec[:, 2]
        qperi[np.where(live)[0]] = qp / RE

    res.update(cat=cat, rcirc=rcirc, jz=jz, qperi=qperi)
    return res


def summary(res):
    w = res['wmass']
    cats = ['fenetre_lunaire', 'disque_interne', 'retombee', 'echappee',
            'emportee_B', 'non_leve']
    return {c: float(w[res['cat'] == c].sum()) for c in cats}
